/**
 * @file
 * @brief Implementation of Geant4 deposition module
 * @remarks Based on code from Mathieu Benoit
 * @copyright Copyright (c) 2017 CERN and the Allpix Squared authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#include "DepositionGeant4Module.hpp"

#include <limits>
#include <string>
#include <utility>

// #include <G4TrackingAction.hh>

#include <G4EmParameters.hh>
#include <G4HadronicProcessStore.hh>
#include <G4LogicalVolume.hh>
#include <G4PhysListFactory.hh>
#include <G4RunManager.hh>
#include <G4StepLimiterPhysics.hh>
#include <G4UImanager.hh>
#include <G4UserLimits.hh>

#include "core/config/exceptions.h"
#include "core/geometry/GeometryManager.hpp"
#include "core/module/exceptions.h"
#include "core/utils/log.h"
#include "objects/DepositedCharge.hpp"
#include "tools/ROOT.h"
#include "tools/geant4.h"

#include "GeneratorActionG4.hpp"
#include "GeneratorActionCustomG4.hpp"
#include "SensitiveDetectorActionG4.hpp"

#include "PhysicsList.hpp"

#define G4_NUM_SEEDS 10

using namespace allpix;

/**
 * Includes the particle source point to the geometry using \ref GeometryManager::addPoint.
 */
DepositionGeant4Module::DepositionGeant4Module(Configuration config, Messenger* messenger, GeometryManager* geo_manager)
    : Module(std::move(config)), messenger_(messenger), geo_manager_(geo_manager), last_event_num_(1),
      run_manager_g4_(nullptr) {
    // Create user limits for maximum step length in the sensor
    user_limits_ = std::make_unique<G4UserLimits>(config_.get<double>("max_step_length", Units::get(1.0, "um")));

    // Set default physics list
    config_.setDefault("physics_list", "FTFP_BERT_LIV");

    config_.setDefault<bool>("output_plots", false);
    config_.setDefault<int>("output_plots_scale", Units::get(100, "ke"));

    // Add the particle source position to the geometry
    geo_manager_->addPoint(config_.get<ROOT::Math::XYZPoint>("beam_position"));
}

/**
 * Module depends on \ref GeometryBuilderGeant4Module loaded first, because it owns the pointer to the Geant4 run manager.
 */
void DepositionGeant4Module::init() {
    // Load the G4 run manager (which is owned by the geometry builder)
    run_manager_g4_ = G4RunManager::GetRunManager();
    if(run_manager_g4_ == nullptr) {
        throw ModuleError("Cannot deposit charges using Geant4 without a Geant4 geometry builder");
    }

    // Suppress all output from G4
    SUPPRESS_STREAM(G4cout);

    // Get UI manager for sending commands
    G4UImanager* ui_g4 = G4UImanager::GetUIpointer();

    // Apply optional PAI model
    if(config_.get<bool>("enable_pai", false)) {
        LOG(TRACE) << "Enabling PAI model on all detectors";
        G4EmParameters::Instance();

        for(auto& detector : geo_manager_->getDetectors()) {
            // Get logical volume
            auto logical_volume = detector->getExternalObject<G4LogicalVolume>("sensor_log");
            if(logical_volume == nullptr) {
                throw ModuleError("Detector " + detector->getName() + " has no sensitive device (broken Geant4 geometry)");
            }
            // Create region
            G4Region* region = new G4Region(detector->getName() + "_sensor_region");
            region->AddRootLogicalVolume(logical_volume.get());

            auto pai_model = config_.get<std::string>("pai_model", "pai");
            auto lcase_model = pai_model;
            std::transform(lcase_model.begin(), lcase_model.end(), lcase_model.begin(), ::tolower);
            if(lcase_model == "pai") {
                pai_model = "PAI";
            } else if(lcase_model == "paiphoton") {
                pai_model = "PAIphoton";
            } else {
                throw InvalidValueError(config_, "pai_model", "model has to be either 'pai' or 'paiphoton'");
            }

            ui_g4->ApplyCommand("/process/em/AddPAIRegion all " + region->GetName() + " " + pai_model);
        }
    }

    // Find the physics list    
    /*
    bool use_decay_physics;
    if(config_.has("gps")) {
        use_decay_physics = false;
    } else {
        use_decay_physics = true;
    }
    G4VModularPhysicsList* physicsList = new PhysicsList(use_decay_physics); 
    */

    G4PhysListFactory physListFactory;
    G4VModularPhysicsList* physicsList = physListFactory.GetReferencePhysList(config_.get<std::string>("physics_list"));

    // Tracking Messenger
    // G4TrackingAction* trackingAction;
    // trackingAction->SetFullChain(false);

    if(physicsList == nullptr) {
        std::string message = "specified physics list does not exists";
        std::vector<G4String> base_lists = physListFactory.AvailablePhysLists();
        message += " (available base lists are ";
        for(auto& base_list : base_lists) {
            message += base_list;
            message += ", ";
        }
        message = message.substr(0, message.size() - 2);
        message += " with optional suffixes for electromagnetic lists ";
        std::vector<G4String> em_lists = physListFactory.AvailablePhysListsEM();
        for(auto& em_list : em_lists) {
            if(em_list.empty()) {
                continue;
            }
            message += em_list;
            message += ", ";
        }
        message = message.substr(0, message.size() - 2);
        message += ")";

        throw InvalidValueError(config_, "physics_list", message);
    } else {
        LOG(INFO) << "Using G4 physics list \"" << config_.get<std::string>("physics_list") << "\"";
    }
    // Register a step limiter (uses the user limits defined earlier)
    physicsList->RegisterPhysics(new G4StepLimiterPhysics());

    // Set the range-cut off threshold for secondary production:
    double production_cut;
    if(config_.has("range_cut")) {
        production_cut = config_.get<double>("range_cut");
        LOG(INFO) << "Setting configured G4 production cut to " << Units::display(production_cut, {"mm", "um"});
    } else {
        // Define the production cut as one fifth of the minimum size (thickness, pitch) among the detectors
        double min_size = std::numeric_limits<double>::max();
        std::string min_detector;
        for(auto& detector : geo_manager_->getDetectors()) {
            auto model = detector->getModel();
            double prev_min_size = min_size;
            min_size =
                std::min({min_size, model->getPixelSize().x(), model->getPixelSize().y(), model->getSensorSize().z()});
            if(min_size != prev_min_size) {
                min_detector = detector->getName();
            }
        }
        production_cut = min_size / 5;
        LOG(INFO) << "Setting G4 production cut to " << Units::display(production_cut, {"mm", "um"})
                  << ", derived from properties of detector \"" << min_detector << "\"";
    }
    ui_g4->ApplyCommand("/run/setCut " + std::to_string(production_cut));

    // Initialize the physics list
    LOG(TRACE) << "Initializing physics processes";
    run_manager_g4_->SetUserInitialization(physicsList);
    // run_manager_g4_->SetUserInitialization(new PhysicsList);
    run_manager_g4_->InitializePhysics();

    // Initialize the full run manager to ensure correct state flags
    run_manager_g4_->Initialize();

    // Build particle generator
    LOG(TRACE) << "Constructing particle source";
    LOG(INFO) << "Using custom? " << config_.get<bool>("custom");
    LOG(INFO) << "Has gps? " << config_.has("gps");

    if(config_.get<bool>("custom")) {
        run_manager_g4_->SetUserAction(new GeneratorActionCustomG4(config_));
    } else {
        run_manager_g4_->SetUserAction(new GeneratorActionG4(config_));
    }

    if(config_.has("gps")) {
        LOG(INFO) << "Using gps source";
        auto macro_infile = config_.get<std::string>("gps");
        if(config_.has("macro_path")) {
            std::string macro_path = config_.get<std::string>("macro_path");
            ui_g4->ApplyCommand("/control/macroPath " + macro_path);
            LOG(INFO) << "Set macro path to: " << macro_path << "\"";
        }
        ui_g4->ApplyCommand("/control/execute " + macro_infile);
    }

    // Get the creation energy for charge (default is silicon electron hole pair energy)
    auto charge_creation_energy = config_.get<double>("charge_creation_energy", Units::get(3.64, "eV"));

    // Loop through all detectors and set the sensitive detector action that handles the particle passage
    bool useful_deposition = false;
    for(auto& detector : geo_manager_->getDetectors()) {
        // Do not add sensitive detector for detectors that have no listeners for the deposited charges
        // FIXME Probably the MCParticle has to be checked as well
        if(!messenger_->hasReceiver(this,
                                    std::make_shared<DepositedChargeMessage>(std::vector<DepositedCharge>(), detector))) {
            LOG(INFO) << "Not depositing charges in " << detector->getName()
                      << " because there is no listener for its output";
            continue;
        }
        useful_deposition = true;

        // Get model of the sensitive device
        auto sensitive_detector_action = new SensitiveDetectorActionG4(this, detector, messenger_, charge_creation_energy);
        auto logical_volume = detector->getExternalObject<G4LogicalVolume>("sensor_log");
        if(logical_volume == nullptr) {
            throw ModuleError("Detector " + detector->getName() + " has no sensitive device (broken Geant4 geometry)");
        }

        // Apply the user limits to this element
        logical_volume->SetUserLimits(user_limits_.get());

        // Add the sensitive detector action
        logical_volume->SetSensitiveDetector(sensitive_detector_action);
        sensors_.push_back(sensitive_detector_action);

        // If requested, prepare output plots
        if(config_.get<bool>("output_plots")) {
            LOG(TRACE) << "Creating output plots";

            // Plot axis are in kilo electrons - convert from framework units!
            int maximum = static_cast<int>(Units::convert(config_.get<int>("output_plots_scale"), "ke"));
            int nbins = 5 * maximum;

            // Create histograms if needed
            std::string plot_name = "deposited_charge_" + sensitive_detector_action->getName();
            charge_per_event_[sensitive_detector_action->getName()] =
                new TH1D(plot_name.c_str(), "deposited charge per event;deposited charge [ke];events", nbins, 0, maximum);
        }
    }

    if(!useful_deposition) {
        LOG(ERROR) << "Not a single listener for deposited charges, module is useless!";
    }

    // Disable verbose messages from processes
    ui_g4->ApplyCommand("/process/verbose 0");
    ui_g4->ApplyCommand("/process/em/verbose 0");
    ui_g4->ApplyCommand("/process/eLoss/verbose 0");
    G4HadronicProcessStore::Instance()->SetVerbose(0);

    // Set the random seed for Geant4 generation
    // NOTE Assumes this is the only Geant4 module using random numbers
    std::string seed_command = "/random/setSeeds ";
    for(int i = 0; i < G4_NUM_SEEDS; ++i) {
        seed_command += std::to_string(static_cast<uint32_t>(getRandomSeed() % INT_MAX));
        if(i != G4_NUM_SEEDS - 1) {
            seed_command += " ";
        }
    }
    ui_g4->ApplyCommand(seed_command);

    // Release the output stream
    RELEASE_STREAM(G4cout);
}

void DepositionGeant4Module::run(unsigned int event_num) {
    // Suppress output stream if not in debugging mode
    IFLOG(DEBUG);
    else {
        SUPPRESS_STREAM(G4cout);
    }

    // Start a single event from the beam
    LOG(TRACE) << "Enabling beam";
    run_manager_g4_->BeamOn(static_cast<int>(config_.get<unsigned int>("number_of_particles", 1)));
    last_event_num_ = event_num;

    // Release the stream (if it was suspended)
    RELEASE_STREAM(G4cout);

    // Dispatch the necessary messages
    for(auto& sensor : sensors_) {
        sensor->dispatchMessages();

        // Fill output plots if requested:
        if(config_.get<bool>("output_plots")) {
            double charge = static_cast<double>(Units::convert(sensor->getDepositedCharge(), "ke"));
            charge_per_event_[sensor->getName()]->Fill(charge);
        }
    }
}

void DepositionGeant4Module::finalize() {
    size_t total_charges = 0;
    for(auto& sensor : sensors_) {
        total_charges += sensor->getTotalDepositedCharge();
    }

    if(config_.get<bool>("output_plots")) {
        // Write histograms
        LOG(TRACE) << "Writing output plots to file";
        for(auto& plot : charge_per_event_) {
            plot.second->Write();
        }
    }

    // Print summary or warns if module did not output any charges
    if(!sensors_.empty() && total_charges > 0 && last_event_num_ > 0) {
        size_t average_charge = total_charges / sensors_.size() / last_event_num_;
        LOG(INFO) << "Deposited total of " << total_charges << " charges in " << sensors_.size() << " sensor(s) (average of "
                  << average_charge << " per sensor for every event)";
    } else {
        LOG(WARNING) << "No charges deposited";
    }
}