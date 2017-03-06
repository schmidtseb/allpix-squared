/**
 * @author Koen Wolters <koen.wolters@cern.ch>
 * @author Mathieu Benoit <benoit@lal.in2p3.fr>
 */

#include "GeometryConstructionG4.hpp"

//#include "G4SDManager.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4VisAttributes.hh"
#include "G4PVDivision.hh"
#include "G4VSolid.hh"
#include "G4ThreeVector.hh"
#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4UserLimits.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4PVParameterised.hh"

#include "CLHEP/Units/SystemOfUnits.h"

#include "BumpsParameterizationG4.hpp"
#include "DetectorModelG4.hpp"

#include "../../core/geometry/PixelDetectorModel.hpp"
#include "../../core/utils/log.h"

using namespace CLHEP;
using namespace std;

using namespace allpix;

G4VPhysicalVolume *GeometryConstructionG4::Construct(){
    // Materials
    // vacuum (FIXME: useless code, only here for later) 
    double z,a,density;
    G4Material *vacuum = new G4Material("Vacuum", z=1 , a=1.01*g/mole, density= 0.0001*g/cm3);
    
    // air
    G4NistManager * nistman = G4NistManager::Instance();
    G4Material *air = nistman->FindOrBuildMaterial("G4_AIR");
    
    // FIXME: stick to vacuum as world material now
    world_material_ = vacuum;
    LOG(DEBUG) << "Material of world: " << world_material_->GetName();
    
    //build the world
    G4Box* world_box = new G4Box("World",
                                   world_size_[0],
                                   world_size_[1],
                                   world_size_[2]);
    world_log_= new G4LogicalVolume(world_box,
                          world_material_,
                          "World",
                          0,0,0);
    
    // set the world to invisible in the viewer
    // FIXME: this should strictly not be defined here, but simplifies things a lot
    G4VisAttributes * invisibleVisAtt = new G4VisAttributes(G4Color(1.0, 0.65, 0.0, 0.1));
    invisibleVisAtt->SetVisibility(false);
    invisibleVisAtt->SetForceSolid(false);
    world_log_->SetVisAttributes ( invisibleVisAtt );
    
    // place the world at the center
    world_phys_ = new G4PVPlacement(0,
                        G4ThreeVector(0.,0.,0.),
                        world_log_,
                        "World",
                        0x0,
                        false,
                        0);
    
    // build the pixel devices
    BuildPixelDevices();
    
    // WARNING: some debug printing here about the placement of the detectors
    
    return world_phys_;
}

// WARNING: A DEFINE HERE THAT SHOULD PROBABLY BE A PARAMETER 
// NOTE: MULTIPLE ALERT HERE TO IDENTIFY PARTS THAT ARE NOT COPIED FROM ALLPIX 1
void GeometryConstructionG4::BuildPixelDevices() {
    
    /* MATERIALS
     * fetch and build the required materials
     */
    
    G4NistManager * nistman = G4NistManager::Instance();
    
    // fetch standard materials
    G4Material *Air = nistman->FindOrBuildMaterial("G4_AIR");
    G4Material * Silicon = nistman->FindOrBuildMaterial("G4_Si");
    G4Material * Epoxy = nistman->FindOrBuildMaterial("G4_PLEXIGLASS"); // FIXME: this is supposed to be epoxy
    
    // fetch elements needed
    G4Element* Sn = new G4Element("Tin", "Sn", 50., 118.710*g/mole);
    G4Element* Pb = new G4Element("Lead","Pb", 82., 207.2*g/mole);

    // create custom materials
    G4Material* Solder = new G4Material("Solder", 8.4*g/cm3, 2);
    Solder->AddElement(Sn,63);
    Solder->AddElement(Pb,37);

    /* VISIBILITY
     * set the visibility of several of the volumes 
     * FIXME: should strictly not be here but simplifies visualization
     */
    
    G4VisAttributes * BoxVisAtt= new G4VisAttributes(G4Color(0,1,1,1));
    BoxVisAtt->SetLineWidth(2);
    BoxVisAtt->SetForceSolid(true);

    G4VisAttributes * ChipVisAtt= new G4VisAttributes(G4Color::Gray());
    ChipVisAtt->SetLineWidth(2);
    ChipVisAtt->SetForceSolid(true);

    G4VisAttributes * BumpBoxVisAtt = new G4VisAttributes(G4Color(0,1,0,1.0));
    BumpBoxVisAtt->SetLineWidth(1);
    BumpBoxVisAtt->SetForceSolid(false);
    BumpBoxVisAtt->SetVisibility(true);

    G4VisAttributes * BumpVisAtt = new G4VisAttributes(G4Color::Yellow());
    BumpVisAtt->SetLineWidth(2);
    BumpVisAtt->SetForceSolid(true);

    G4VisAttributes * pcbVisAtt = new G4VisAttributes(G4Color::Green());
    pcbVisAtt->SetLineWidth(1);
    pcbVisAtt->SetForceSolid(true);

    G4VisAttributes * guardRingsVisAtt = new G4VisAttributes(G4Color(0.5,0.5,0.5,1));
    guardRingsVisAtt->SetLineWidth(1);
    guardRingsVisAtt->SetForceSolid(true);

    G4VisAttributes * wrapperVisAtt = new G4VisAttributes(G4Color(1,0,0,0.9));
    wrapperVisAtt->SetLineWidth(1);
    wrapperVisAtt->SetForceSolid(false);
    wrapperVisAtt->SetVisibility(false);

    /* NAMES
     * define the global names for all the elements in the setup
     */
    pair<G4String, G4String> wrapperName = make_pair("wrapper", "");
    pair<G4String, G4String> PCBName = make_pair("PCB", "");
    pair<G4String, G4String> BoxName = make_pair("Box", "");
    //pair<G4String, G4String> CoverlayerName = make_pair("Coverlayer", "");
    pair<G4String, G4String> SliceName = make_pair("Slice", "");
    pair<G4String, G4String> GuardRingsExtName = make_pair("GuardRingsExt", "");
    pair<G4String, G4String> GuardRingsName = make_pair("GuardRings", "");
    pair<G4String, G4String> PixelName = make_pair("Pixel", "");
    pair<G4String, G4String> ChipName = make_pair("Chip", "");
    pair<G4String, G4String> SDName = make_pair("BoxSD", "");
    pair<G4String, G4String> BumpName = make_pair("Bump", "");
    pair<G4String, G4String> BumpBoxName = make_pair("BumpBox", "");

    // WARNING: THIS SHOULD PROBABLY BE REMOVED
    // User limits applied only to Si wafers.  Setting step.
#define PARAM_MAX_STEP_LENGTH DBL_MAX
    G4UserLimits *user_limits = new G4UserLimits(PARAM_MAX_STEP_LENGTH);

    /* CONSTRUCTION
     * construct the detectors part by part
     */
    
    std::vector<std::shared_ptr<Detector>> detectors = geo_manager_->getDetectors();
    LOG(DEBUG) << "Building " << detectors.size() << " device(s) ..." ;
    
    auto detItr = detectors.begin();
    for( ; detItr != detectors.end() ; detItr++){
        // get pointers for the model of the detector
        std::shared_ptr<DetectorModelG4> model_g4 = std::make_shared<DetectorModelG4>();
        std::shared_ptr<PixelDetectorModel> model = std::dynamic_pointer_cast<PixelDetectorModel>((*detItr)->getModel());
        
        // ignore all non-pixel detectors
        if(model == 0){
            LOG(WARNING) << "cannot build a G4 model for any non-pixel detectors yet... ignoring detector named " << (*detItr)->getName();
            continue;
        }
        
        // NOTE: create temporary internal integer
        // FIXME: this is not necessary unique! if this is necessary generate it internally!
        std::hash<std::string> hash_func;
        size_t temp_g4_id = hash_func((*detItr)->getName());
        
        LOG(DEBUG) << "start creating G4 detector " << (*detItr)->getName() << " (" << temp_g4_id << ")";
        
        /* POSITIONS
         * calculate the positions of all the elements
         */
        G4double bump_height = model->GetBumpHeight();

        // positions of all the elements
        G4ThreeVector posCoverlayer(0,0,0);
        G4ThreeVector posDevice(0,0,0);
        G4ThreeVector posBumps(0,0,0);
        G4ThreeVector posChip(0,0,0);
        G4ThreeVector posPCB(0,0,0);
        
        // ALERT: NO WRAPPER ENHANCEMENTS
        // apply position offset for the detector due to the enhancement
        /*if ( m_vectorWrapperEnhancement.find(temp_g4_id) != m_vectorWrapperEnhancement.end() ) {
            posDevice.setX(posDevice.x() - m_vectorWrapperEnhancement[*(*detItr)].x()/2.);
            posDevice.setY(posDevice.y() - m_vectorWrapperEnhancement[*(*detItr)].y()/2.);
            posDevice.setZ(posDevice.z() - m_vectorWrapperEnhancement[*(*detItr)].z()/2.);
        } else {*/
            posDevice.setX(posDevice.x() );
            posDevice.setY(posDevice.y() );
            posDevice.setZ(posDevice.z() );
        //}
        
        // calculation of position of the different physical volumes
        if ( model->GetHalfChipZ() != 0 ) {
            posBumps.setX( posDevice.x() );
            posBumps.setY( posDevice.y() );
            posBumps.setZ( posDevice.z()
            //wrapperHZ
            - model->GetHalfSensorZ()
            - 2.*model->GetHalfCoverlayerZ()
            - (bump_height/2.)
            );
            
            posChip.setX( posDevice.x() + model->GetChipXOffset() );
            posChip.setY( posDevice.y() + model->GetChipYOffset() );
            posChip.setZ( posDevice.z()
            //wrapperHZ
            - model->GetHalfSensorZ()
            - 2.*model->GetHalfCoverlayerZ()
            - bump_height
            - model->GetHalfChipZ()
            + model->GetChipZOffset()
            );
            
        } else {
            // make sure no offset because of bumps for PCB is calculated if chip is not included
            bump_height = 0;
        }
        posPCB.setX(posDevice.x()-1*model->GetSensorXOffset());
        posPCB.setY(posDevice.y()-1*model->GetSensorYOffset());
        posPCB.setZ(posDevice.z()
        //wrapperHZ
        - model->GetHalfSensorZ()
        - 2.*model->GetHalfCoverlayerZ()
        - bump_height
        - 2.*model->GetHalfChipZ()
        - model->GetHalfPCBZ()
        );
        
        LOG(DEBUG) << "local relative positions of the elements";
        LOG(DEBUG) << "- Coverlayer position  : " << posCoverlayer ;
        LOG(DEBUG) << "- Sensor position      : " << posDevice ;
        LOG(DEBUG) << "- Bumps position       : " << posBumps ;
        LOG(DEBUG) << "- Chip position        : " << posChip ;
        LOG(DEBUG) << "- PCB position         : " << posPCB ;
        
        /* NAMES
         * define the local names of the specific detectors
         */
        
        std::string name = (*detItr)->getName();
        wrapperName.second = wrapperName.first + "_" + name;
        PCBName.second = PCBName.first + "_" + name;
        BoxName.second = BoxName.first + "_" + name;
        //CoverlayerName.second = CoverlayerName.first + "_" + name;
        GuardRingsExtName.second = GuardRingsExtName.first  + "_" + name;
        GuardRingsName.second = GuardRingsName.first  + "_" + name;
        SliceName.second = BoxName.first + "_" + name;
        PixelName.second = BoxName.first + "_" + name;
        ChipName.second = ChipName.first + "_" + name;
        SDName.second = SDName.first + "_" + name; 
        BumpName.second = BumpName.first + "_" + name;
        BumpBoxName.second = BumpBoxName.first + "_" + name;
        
        /* WRAPPER
         * the wrapper is the box around all of the detector 
         * 
         * FIXME: this should be centered at the correct place
         */
        
        // The wrapper might be enhanced when the user set up
        //  Appliances to the detector (extra layers, etc).
        G4double wrapperHX = model->GetHalfWrapperDX();
        G4double wrapperHY = model->GetHalfWrapperDY();
        G4double wrapperHZ = model->GetHalfWrapperDZ();
        
        LOG(DEBUG) << "Wrapper Dimensions [mm] : " << wrapperHX/mm << " " << wrapperHY/mm << " " << wrapperHZ/mm;
        
        G4Box* wrapper_box = new G4Box(wrapperName.second,
                                       2.*wrapperHX,
                                       2.*wrapperHY,
                                       2.*wrapperHZ);
        model_g4->wrapper_log = new G4LogicalVolume(wrapper_box,
                                                       world_material_,
                                                       wrapperName.second+"_log");
        model_g4->wrapper_log->SetVisAttributes(wrapperVisAtt);
        
        // WARNING: get a proper geometry lib
        std::tuple<double, double, double> pos_tup = (*detItr)->getPosition();
        G4ThreeVector posWrapper;
        posWrapper[0] = std::get<0>(pos_tup);
        posWrapper[1] = std::get<1>(pos_tup);
        posWrapper[2] = std::get<2>(pos_tup);
        
        // ALERT: NO WRAPPER ENHANCEMENTS
        // Apply the enhancement to the medipixes
        // We can have N medipixes and K enhancements, where K<=N.
        // For instance, for 2 medipixes.  We can have.
        // medipix 1 --> with enhancement
        // medipix 2 --> no enhancement
        /* if ( m_vectorWrapperEnhancement.findtemp_g4_id != m_vectorWrapperEnhancement.end() ) {
            wrapperHX += m_vectorWrapperEnhancement[*(*detItr)].x()/2.; // half
            wrapperHY += m_vectorWrapperEnhancement[*(*detItr)].y()/2.;
            wrapperHZ += m_vectorWrapperEnhancement[*(*detItr)].z()/2.;
        } */
    
        // Apply position Offset for the wrapper due to the enhancement
        /*if ( m_vectorWrapperEnhancement.findtemp_g4_id != m_vectorWrapperEnhancement.end() ) {
            posWrapper.setX(posWrapper.x() + m_vectorWrapperEnhancement[*(*detItr)].x()/2.);
            posWrapper.setY(posWrapper.y() + m_vectorWrapperEnhancement[*(*detItr)].y()/2.);
            posWrapper.setZ(posWrapper.z() + m_vectorWrapperEnhancement[*(*detItr)].z()/2.);
        } else {*/
            //posWrapper.setX(posWrapper.x() - sensorOffsetX );
            //posWrapper.setY(posWrapper.y() - sensorOffsetY );
            posWrapper.setX(posWrapper.x() );
            posWrapper.setY(posWrapper.y() );
            posWrapper.setZ(posWrapper.z() );
        //}
                
        // WARNING: get a proper geometry lib
        std::tuple<double, double, double> rot_tup = (*detItr)->getOrientation();
        G4RotationMatrix *rotWrapper = new G4RotationMatrix(std::get<0>(rot_tup), std::get<1>(rot_tup), std::get<2>(rot_tup));
            
        // starting at user position --> vector pos
        model_g4->wrapper_phys = new G4PVPlacement(rotWrapper,
                                                posWrapper,
                                                model_g4->wrapper_log,
                                                wrapperName.second+"_phys",
                                                world_log_,
                                                false,
                                                temp_g4_id, // copy number
                                                true);
        
        /* DEVICE
         * the sensitive detector 
         */
        
        // create the general box containing the sensor
        G4Box * Box_box = new G4Box(BoxName.second,
                                    model->GetHalfSensorX(),
                                    model->GetHalfSensorY(),
                                    model->GetHalfSensorZ());
        
        // create the box containing the slices and pixels
        // The Si wafer is placed respect to the wrapper.
        // Needs to be pushed -half Si wafer in z direction        
        model_g4->box_log = new G4LogicalVolume(Box_box,
                                                   Silicon,
                                                   BoxName.second+"_log");
        
        model_g4->box_log->SetVisAttributes(BoxVisAtt);
        
        model_g4->box_phys = new G4PVPlacement( 0,
                                                posDevice,
                                                model_g4->box_log,
                                                BoxName.second+"_phys",
                                                model_g4->wrapper_log, // mother log
                                                false,
                                                temp_g4_id, // copy number
                                                true); // check overlap
        
        
        // create the slices and pixels (replicas)
        G4Box * Box_slice = new G4Box(SliceName.first,
                                      model->GetHalfPixelX(),
                                      model->GetHalfSensorY(),
                                      model->GetHalfSensorZ());
        model_g4->slice_log = new G4LogicalVolume(Box_slice,
                                                     Silicon,
                                                     SliceName.second); // 0,0,0);
                                                     
        G4Box * Box_pixel = new G4Box(PixelName.first,
        model->GetHalfPixelX(),
        model->GetHalfPixelY(),
        model->GetHalfPixelZ());                
        model_g4->pixel_log = new G4LogicalVolume(Box_pixel,
                                                    Silicon,
                                                    PixelName.second); // 0,0,0);
                       
        // set the user limit (FIXME: is this needed / this is currently fixed)
        if (user_limits) model_g4->pixel_log->SetUserLimits(user_limits);
        
        // place the slices
        new G4PVDivision(SliceName.second,
                        model_g4->slice_log,
                        model_g4->box_log,
                        kXAxis,
                        //model->GetPixelX(),
                        model->GetNPixelsX(),
                        0); // offset
        
        // place the pixels
        new G4PVDivision(PixelName.second,
                        model_g4->pixel_log,
                        model_g4->slice_log,
                        kYAxis,
                        //model->GetPixelY(),
                        model->GetNPixelsY(),
                        0); // offset
        
        /* BUMPS
         * the construction of the bump bonds
         */
        
        // construct the bumps only if necessary
        if ( model->GetBumpHeight() != 0.0 and model->GetHalfChipZ() != 0 ) {
            // define types from parameters
            G4double bump_radius = model->GetBumpRadius();
            G4double bump_dr =  model->GetBumpDr();
            G4Sphere * aBump_Sphere = new G4Sphere(BumpName.first+"sphere",0,bump_radius,0,360*deg,0,360*deg);
            G4Tubs * aBump_Tube= new G4Tubs(BumpName.first+"Tube", 0., bump_radius-bump_dr, bump_height/2., 0., 360 *deg);
            G4UnionSolid * aBump = new G4UnionSolid(BumpName.first,aBump_Sphere,aBump_Tube);
            
            //create the volume containing the bumps
            G4Box * Bump_Box = new G4Box(BumpBoxName.first,
                                        model->GetHalfSensorX(),
                                        model->GetHalfSensorY(),
                                        bump_height/2.);
            
            //create the logical volume
            model_g4->bumps_log = new G4LogicalVolume(Bump_Box,
                                                            Air,
                                                            BumpBoxName.second+"_log");
            model_g4->bumps_log->SetVisAttributes(BumpBoxVisAtt);
            
            // place the general bumps volume
            model_g4->bumps_phys = new G4PVPlacement(0,
                                                    posBumps,
                                                    model_g4->bumps_log,
                                                    BumpBoxName.second+"_phys",
                                                    model_g4->wrapper_log, // mother log
                                                    false,
                                                    temp_g4_id, // copy number
                                                    true); // check overlap
            
             
            // create the individual bumps through a parameterization
            model_g4->bumps_cell_log = new G4LogicalVolume(aBump,Solder,BumpBoxName.second+"_log" );
            model_g4->bumps_cell_log->SetVisAttributes(BumpVisAtt);
            
            model_g4->parameterization_ = new BumpsParameterizationG4( model );
            G4int NPixTot = model->GetNPixelsX()*model->GetNPixelsY();
            new G4PVParameterised(BumpName.second+"phys",
                                  model_g4->bumps_cell_log,        // logical volume
                                  model_g4->bumps_log,             // mother volume
                                  kUndefined,                         // axis
                                  NPixTot,                            // replicas
                                  model_g4->parameterization_);       // G4VPVParameterisation
        }
        

        /* COVER LAYER 
         * ALERT: REMOVED SUPPORT COVER LAYERS
         */
        
        // If the coverlayer is requested.  It is forced to fit the sensor in X,Y.
        // The user can pick the thickness.
        
        /* G4Box * Coverlayer_box = nullptr;
        if ( model->IsCoverlayerON() ) {
            Coverlayer_box = new G4Box(
                CoverlayerName.second,
                model->GetHalfSensorX(),
                                       model->GetHalfSensorY(),
                                       model->GetHalfCoverlayerZ()
            );
        }*/
        
        /*if ( model->IsCoverlayerON() ) {
            
            posCoverlayer.setX( posDevice.x() );
            posCoverlayer.setY( posDevice.y() );
            posCoverlayer.setZ(
                posDevice.z()
                -model->GetHalfSensorZ()
                -model->GetHalfCoverlayerZ()
            );
        }*/
        
        // coverlayer material
        /*if ( model->IsCoverlayerON() ) {
            
            // Find out about the material that the user requested
            G4Material * covermat = nistman->FindOrBuildMaterial( model->GetCoverlayerMat() );
            
            // If this is an inexistent choice, them force Aluminum
            if ( covermat == nullptr ) {
                covermat = nistman->FindOrBuildMaterial( "G4_Al" );
            }
            
            m_Coverlayer_log[temp_g4_id] = new G4LogicalVolume( Coverlayer_box,
                                                               covermat,
                                                               CoverlayerName.second+"_log");
            
            m_Coverlayer_log[temp_g4_id]->SetVisAttributes(CoverlayerVisAtt);
        }*/
        
        // coverlayer placement
        /*if ( model->IsCoverlayerON() ) {
            
            m_Coverlayer_phys[temp_g4_id] = new G4PVPlacement(
                0,
                posCoverlayer,
                m_Coverlayer_log[temp_g4_id],
                                                             CoverlayerName.second+"_phys",
                                                             model_g4->wrapper_log, // mother log
                                                             false,
                                                             temp_g4_id, // copy number
                                                             true); // check overlap
            
        }*/
        
        /* GUARD RINGS
         * rings around the sensitive device 
         */
        
        // create the box volumes for the guard rings
        G4Box * Box_GuardRings_Ext = new G4Box(
            GuardRingsExtName.second,
            model->GetHalfSensorX() + (model->GetSensorExcessHRight() + model->GetSensorExcessHLeft()),
                                               model->GetHalfSensorY() + (model->GetSensorExcessHTop() + model->GetSensorExcessHBottom()),
                                               model->GetHalfSensorZ()); // same depth as the sensor
        
        G4VSolid * Solid_GuardRings =  new G4SubtractionSolid(GuardRingsName.second,
                                                              Box_GuardRings_Ext,
                                                              Box_box);
        
        // create the logical volume for the guard rings
        model_g4->guard_rings_log = new G4LogicalVolume(Solid_GuardRings,
                                                        Silicon,
                                                        GuardRingsName.second+"_log");
        model_g4->guard_rings_log->SetVisAttributes(guardRingsVisAtt);
        
        // place the guard rings
        model_g4->guard_rings_phys = new G4PVPlacement(0,
                                                        posDevice,
                                                        model_g4->guard_rings_log,
                                                        GuardRingsName.second+"_phys",
                                                        model_g4->wrapper_log, // mother log
                                                        false,
                                                        0,//temp_g4_id, // copy number
                                                        true); // check overlap
        
        /* CHIPS
         * the chips connected to the PCBs
         */
        
        model_g4->chip_log = 0;
        model_g4->chip_phys = 0;
        if ( model->GetHalfChipZ() != 0 ) {
            // create the chip box volume
            G4Box * Chip_box = new G4Box(ChipName.second,
                                 model->GetHalfChipX(),
                                 model->GetHalfChipY(),
                                 model->GetHalfChipZ());
            
            // create the logical volume
            // The Si wafer is placed respect to the wrapper.
            // Needs to be pushed -half Si wafer in z direction
            model_g4->chip_log = new G4LogicalVolume(Chip_box,
                                                        Silicon,
                                                        ChipName.second+"_log");
            model_g4->chip_log->SetVisAttributes(ChipVisAtt);
            
            // place the chip
            model_g4->chip_phys = new G4PVPlacement(0,
                                                    posChip,
                                                    model_g4->chip_log,
                                                    ChipName.second+"_phys",
                                                    model_g4->wrapper_log, // mother log
                                                    false,
                                                    temp_g4_id, // copy number
                                                    true); // check overlap
        }
        
        /* PCB
         * global pcb chip 
         */
        
        model_g4->PCB_log = 0;
        model_g4->PCB_phys = 0;
        if ( model->GetHalfPCBZ() != 0 ) {
            // create the box containing the PCB
            G4Box * PCB_box = new G4Box(PCBName.second,
                                        model->GetHalfPCBX(),
                                        model->GetHalfPCBY(),
                                        model->GetHalfPCBZ());
        
            // create the logical volume for the PCB
            // The PCB is placed respect to the wrapper.
            // Needs to be pushed -half Si wafer in z direction
            model_g4->PCB_log = new G4LogicalVolume(PCB_box,
                                                       Epoxy,
                                                       PCBName.second+"_log");
            model_g4->PCB_log->SetVisAttributes(pcbVisAtt);
            
            // place the PCB
            model_g4->PCB_phys = new G4PVPlacement( 0,
                                                    posPCB,
                                                    model_g4->PCB_log,
                                                    PCBName.second+"_phys",
                                                    model_g4->wrapper_log, // mother log
                                                    false,
                                                    temp_g4_id,
                                                    true); // copy number
        }
        
        // WARNING: the construction of the sensitive devices is currently done by the deposition module
        
        // WARNING: temperature, flux, magnetic and electric field are at the moment not part of the geometry
        
        // add this geant4 model to the detector
        (*detItr)->setExternalModel(model_g4);
        
        LOG(DEBUG) << "detector " << (*detItr)->getName() << " ... done" ;
                                                                                                
    }
}