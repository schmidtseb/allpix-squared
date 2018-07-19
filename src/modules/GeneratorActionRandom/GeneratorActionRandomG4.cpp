/**
 * @file
 * @brief Implements the particle generator
 * @remark Based on code from John Idarraga
 * @copyright Copyright (c) 2017 CERN and the Allpix Squared authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#include "GeneratorActionRandomG4.hpp"
#include "Math/Point3D.h"

#include <limits>
#include <memory>
#include <Eigen/Core>

#include <G4Event.hh>
#include <G4GeneralParticleSource.hh>
#include <G4SingleParticleSource.hh>
#include <G4ParticleDefinition.hh>
#include <G4ParticleTable.hh>
#include <G4IonTable.hh>

#include "core/config/exceptions.h"
#include "core/utils/log.h"
#include "core/utils/unit.h"
#include "tools/geant4.h"

#include <TBranchElement.h>
#include <TClass.h>
#include <TDirectory.h>

using namespace allpix;

GeneratorActionRandomG4::GeneratorActionRandomG4(const Configuration& config)
    : particle_source_(std::make_unique<G4GeneralParticleSource>()) {
    // === Load from config ===
    // Set verbosity of source to off
    particle_source_->SetVerbosity(0);

    // Find Geant4 particle
    auto particle_type = config.get<std::string>("particle_type", "");
    std::transform(particle_type.begin(), particle_type.end(), particle_type.begin(), ::tolower);
    auto particle_code = config.get<int>("particle_code", 0);
    G4ParticleDefinition* particle = nullptr;

    auto pdg_table = G4ParticleTable::GetParticleTable();
        
    if(!particle_type.empty() && particle_code != 0) {
        if(pdg_table->FindParticle(particle_type) == pdg_table->FindParticle(particle_code)) {
            LOG(WARNING) << "particle_type and particle_code given. Continuing because they match.";
            particle = pdg_table->FindParticle(particle_code);
            if(particle == nullptr) {
                throw InvalidValueError(config, "particle_code", "particle code does not exist.");
            }
        } else {
            throw InvalidValueError(
                config, "particle_type", "Given particle_type does not match particle_code. Please remove one of them.");
        }
    } else if(particle_type.empty() && particle_code == 0) {
        throw InvalidValueError(config, "particle_code", "Please set particle_code or particle_type.");
    } else if(particle_code != 0) {
        particle = pdg_table->FindParticle(particle_code);
        if(particle == nullptr) {
            throw InvalidValueError(config, "particle_code", "particle code does not exist.");
        }
    } else {
        particle = pdg_table->FindParticle(particle_type);
        if(particle == nullptr) {
            throw InvalidValueError(config, "particle_type", "particle type does not exist.");
        }
    }  

    LOG(DEBUG) << "Using particle " << particle->GetParticleName() << " (ID " << particle->GetPDGEncoding() << ").";

    // Get angular distribution
    ang_dist_type = config.get<std::string>("ang_dist_type", "None");
    if(ang_dist_type == "pyramid") {
        G4ThreeVector phant_pos_conf = config.get<G4ThreeVector>("phantom_position");
        phantom_position = Eigen::Vector3d(phant_pos_conf.x(), phant_pos_conf.y(), phant_pos_conf.z());

        // Half length of phantom
        phantom_size = config.get<double>("phantom_size", Units::get(15.0, "mm"));

        // Source position
        G4ThreeVector beam_pos_conf = config.get<G4ThreeVector>("beam_position");

        source_position = Eigen::Vector3d(beam_pos_conf.x(), beam_pos_conf.y(), beam_pos_conf.z());
    }

    // Get number of events to simulate
    number_of_events = config.get<int>("number_of_events");

    // Use energy only?
    energy_only = config.get<bool>("energy_only", false);

    // Invert z?
    if (config.get<bool>("invert_z", false)) {
        invert_z = true;
    } else {
        invert_z = false;
    } 

    // Initialize random stuff
    InitRandom(config);

    // Initialize source
    InitSource(config, particle);
}

void GeneratorActionRandomG4::InitSource(const Configuration& config, G4ParticleDefinition* particle) {
    // Get source
    single_source = particle_source_->GetCurrentSource();
    particle_source_->SetCurrentSourceIntensity(1.);

    // Set global parameters of the source
    single_source->SetNumberOfParticles(1);
    single_source->SetParticleDefinition(particle);
    // Set the primary track's start time in for the current event to zero:
    single_source->SetParticleTime(0.0);

    // Set energy parameters
    single_source->GetEneDist()->SetEnergyDisType(config.get<std::string>("energy_dis_type", "Mono"));

    // Get source offset
    G4ThreeVector source_offset = config.get<G4ThreeVector>("offset", G4ThreeVector(0, 0, 0));
    x_offset = source_offset.x();
    y_offset = source_offset.y();
    z_offset = source_offset.z();

    // Parallel source energy
    E_parallel = config.get<double>("parallel_source_energy");

    // Invert the z-axis? 
    // No need to invert x- and y-axes due to symmetries
    if (invert_z) {
        z_invert = -1;
    } else {
        z_invert = 1;
    }
}

void GeneratorActionRandomG4::InitRandom(const Configuration& config) {
    // Load ROOT-File
    root_file = TFile::Open(config.get<std::string>("file").c_str(), "READ");
    LOG(INFO) << "Loaded ROOT-file: " << config.get<std::string>("file");

    // Get tree
    root_tree = static_cast<TTree*>(root_file->Get(config.get<std::string>("tree").c_str()));
    LOG(INFO) << "Loaded ROOT-tree: " << config.get<std::string>("tree");

    num_entries = root_tree->GetEntries();
    LOG(INFO) << "Number of entries: " << num_entries;

    // Get bounding box dimensions
    // Position of top right corner
    bounding_box_x = root_tree->GetMaximum("x");
    bounding_box_y = root_tree->GetMaximum("y");

    // Get attributes
    root_tree->SetBranchAddress(config.get<std::string>("E_attribute").c_str(), &E_branch, nullptr);
    if (!energy_only) {
        root_tree->SetBranchAddress(config.get<std::string>("x_attribute").c_str(), &x_branch, nullptr);
        root_tree->SetBranchAddress(config.get<std::string>("y_attribute").c_str(), &y_branch, nullptr);
        root_tree->SetBranchAddress(config.get<std::string>("z_attribute").c_str(), &z_branch, nullptr);
        root_tree->SetBranchAddress(config.get<std::string>("phi_attribute").c_str(), &phi_branch, nullptr);
        root_tree->SetBranchAddress(config.get<std::string>("theta_attribute").c_str(), &theta_branch, nullptr);
        LOG(INFO) << "Set branches";
    } 

    if (ang_dist_type == "pyramid") {
        // Vector of source-phantom
        Eigen::Vector3d v = phantom_position - source_position;
        LOG(INFO) << phantom_position;
        // LOG(INFO) << source_position;

        // Convert to spherical coordinates
        r_p = v.norm();
        phi_p = std::atan2(v.x(), v.z());
        theta_p = std::acos(v.y() / r_p);

        // Get angular range
        theta_rand = std::atan(phantom_size/r_p);
        phi_rand = theta_rand;

        // Init phi and theta ranges
        minPhi = phi_p - phi_rand;
        maxPhi = phi_p + phi_rand;

        minTheta = theta_p - theta_rand;
        maxTheta = theta_p + theta_rand;

        // Create uniform random generator
        uniform_distribution_pyramid = std::uniform_real_distribution<double>(0, 1);
    }

    // Create binomial random generator
    execute_parallel = false;

    double probability = 1.;
    if(config.has("poly_parameters")) {
        poly_parameters = config.getArray<double>("poly_parameters");
        double bounding_box_area = 2 * bounding_box_x * 2 * bounding_box_y / (std::pow(Units::get(1., "mm"), 2));
        double phantom_area = config.get<double>("phantom_area", Units::get(30., "cm") * Units::get(30, "cm"));
        double energy = config.get<double>("parallel_source_energy") / Units::get(1., "keV");

        // Probability to hit the Dosemeter Bounding Box when irradiating the phantom
        double pDBB = GetPolynomValue(energy, poly_parameters) / 100;
        LOG(INFO) << "Probability to hit the Dosemeter Bounding Box: " << pDBB;
        probability = bounding_box_area / (bounding_box_area * (1 - pDBB) + phantom_area * pDBB);
        LOG(INFO) << "Fraction of parallel and backscattered events: " << probability;

        double psim = pDBB + bounding_box_area/phantom_area * (1 - pDBB);
        LOG(INFO) << "Total number of events: " << number_of_events / psim;
    } else if(config.has("parallel_source_intensity")) {
        probability = config.get<double>("parallel_source_intensity", 1.);
    } else if(config.has("phantom_events")) {
        probability = num_entries / config.get<double>("phantom_events");
    }

    std::binomial_distribution<int> binom_dist(number_of_events, probability);
    num_parallel = binom_dist(random_generator);
    LOG(INFO) << "Number of parallel events: " << num_parallel;

    // Create uniform random generator
    uniform_distribution = std::uniform_int_distribution<long long int>(1, number_of_events);
    LOG(INFO) << "Created random generator";
    idx = 1;
    main_idx = 1;
}

/* Parallal / Backscattering - Fraction Calculation */
double GeneratorActionRandomG4::GetPolynomValue(double x, std::vector<double> parameters) {
    LOG(INFO) << "Polynomials: ";
    double result = 0;
    int n = 0;
    for (auto p : parameters) {
        LOG(INFO) << p;
        result += p * std::pow(x, n);
        n++;
    }

    return result;
}

void GeneratorActionRandomG4::GeneratePrimariesRandom() {
    // Randomly select event from tree
    root_tree->GetEntry( idx ); // uniform_distribution(random_generator) );
    if (++idx > num_entries) {
        idx = 1;
    }

    // Set energy, position, and direction of particle source to sampled values
    // Energy
    single_source->GetEneDist()->SetMonoEnergy(E_branch / 1000.);
    if (energy_only) {
        return;
    }

    // Position
    single_source->GetPosDist()->SetCentreCoords(G4ThreeVector(x_branch + x_offset, y_branch + y_offset, z_branch*z_invert + z_offset));

    // Direction
    // Transform spherical to cartesian coordinates
    double x_dir = std::cos(phi_branch) * std::sin(theta_branch);
    double y_dir = std::sin(phi_branch) * std::sin(theta_branch);
    double z_dir = std::cos(theta_branch) * z_invert;

    single_source->GetAngDist()->SetParticleMomentumDirection(G4ThreeVector(x_dir, y_dir, z_dir));
}

void GeneratorActionRandomG4::GeneratePrimariesParallel() {
    execute_parallel = true;

    // Set source parameters
    // Energy
    single_source->GetEneDist()->SetMonoEnergy(E_parallel);

    // Position & Direction
    single_source->GetAngDist()->SetAngDistType("planar");

    single_source->GetPosDist()->SetPosDisType("Plane");
    single_source->GetPosDist()->SetPosDisShape("Square");
    if (invert_z) {
        single_source->GetPosDist()->SetCentreCoords(G4ThreeVector(0, 0, -10));
        single_source->GetAngDist()->SetParticleMomentumDirection(G4ThreeVector(0, 0, 1));
    } else {
        single_source->GetPosDist()->SetCentreCoords(G4ThreeVector(0, 0, 10));
        single_source->GetAngDist()->SetParticleMomentumDirection(G4ThreeVector(0, 0, -1));
    }
    single_source->GetPosDist()->SetHalfX(bounding_box_x);
    single_source->GetPosDist()->SetHalfY(bounding_box_y);
}

void GeneratorActionRandomG4::GeneratePrimariesPyramid() {
    // Generate random
    // Phi
    double randPhi = uniform_distribution_pyramid(random_generator_pyramid);
    double phi = minPhi + (maxPhi - minPhi) * randPhi;

    // Theta
    double randTheta = uniform_distribution_pyramid(random_generator_pyramid);
    double cosTheta = -randTheta * (std::cos(minTheta) - std::cos(maxTheta)) + std::cos(minTheta);
    double sinTheta = std::sqrt(1 - cosTheta*cosTheta);

    // Transform spherical to cartesian coordinates
    double z = sinTheta*std::cos(phi);
    double x = sinTheta*std::sin(phi);
    double y = cosTheta;

    // Set beam direction
    // auto single_source = particle_source_->GetCurrentSource();
    single_source->GetPosDist()->SetCentreCoords(G4ThreeVector(source_position.x(), source_position.y(), source_position.z()));
    single_source->GetAngDist()->SetAngDistType("planar");

    single_source->GetAngDist()->SetParticleMomentumDirection(G4ThreeVector(x, y, z));

    LOG(DEBUG) << "Particle data (x, y, z) = (" << x << ", " << y << ", " << z << ")";
}

/**
 * Called automatically for every event
 */
void GeneratorActionRandomG4::GeneratePrimaries(G4Event* event) {
    // Set event attributes if ROOT-source is used
    if (energy_only && ang_dist_type == "pyramid") {
        GeneratePrimariesRandom();
        GeneratePrimariesPyramid();
    } else {
        if (main_idx <= (number_of_events - num_parallel)) {
            GeneratePrimariesRandom();
        } else if (!execute_parallel) {
            GeneratePrimariesParallel();
        }
    }

    // Generate event
    particle_source_->GeneratePrimaryVertex(event);

    main_idx++;
}
