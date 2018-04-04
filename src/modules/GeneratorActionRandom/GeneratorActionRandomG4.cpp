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

    // Get bounding box dimensions
    // Position of top right corner
    bounding_box_x = root_tree->GetMaximum("x");
    bounding_box_y = root_tree->GetMaximum("y");

    // Get source offset
    G4ThreeVector source_offset = config.get<G4ThreeVector>("offset", G4ThreeVector(0, 0, 0));
    x_offset = source_offset.x();
    y_offset = source_offset.y();
    z_offset = source_offset.z();

    // Parallel source energy
    E_parallel = config.get<double>("parallel_source_energy");

    // Invert the z-axis? 
    // No need to invert x- and y-axes due to symmetries
    if (config.get<bool>("invert_z", false)) {
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

    // Get attributes
    root_tree->SetBranchAddress(config.get<std::string>("E_attribute").c_str(), &E_branch, nullptr);

    root_tree->SetBranchAddress(config.get<std::string>("x_attribute").c_str(), &x_branch, nullptr);
    root_tree->SetBranchAddress(config.get<std::string>("y_attribute").c_str(), &y_branch, nullptr);
    root_tree->SetBranchAddress(config.get<std::string>("z_attribute").c_str(), &z_branch, nullptr);

    root_tree->SetBranchAddress(config.get<std::string>("phi_attribute").c_str(), &phi_branch, nullptr);
    root_tree->SetBranchAddress(config.get<std::string>("theta_attribute").c_str(), &theta_branch, nullptr);
    LOG(INFO) << "Set branches";

    // Create binomial random generator
    execute_parallel = false;
    double probability = config.get<double>("parallel_source_intensity");
    std::binomial_distribution<int> binom_dist(static_cast<int>(num_entries), probability);
    num_parallel = binom_dist(random_generator);
    LOG(INFO) << "Number of parallel events: " << num_parallel;

    // Create uniform random generator
    uniform_distribution = std::uniform_int_distribution<long long int>(1, num_entries);
    LOG(INFO) << "Created random generator";
    idx = 1;
    main_idx = 1;
}

void GeneratorActionRandomG4::GeneratePrimariesRandom() {
    // Randomly select event from tree
    root_tree->GetEntry( idx ); // uniform_distribution(random_generator) );
    if (++idx > num_entries) {
        idx = 1;
    }

    // Set energy, position, and direction of particle source to sampled values
    // Energy
    single_source->GetEneDist()->SetMonoEnergy(E_branch / 1000);

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
    single_source->GetAngDist()->SetParticleMomentumDirection(G4ThreeVector(0, 0, 1));

    single_source->GetPosDist()->SetPosDisType("Plane");
    single_source->GetPosDist()->SetPosDisShape("Square");
    single_source->GetPosDist()->SetCentreCoords(G4ThreeVector(0, 0, -10));
    single_source->GetPosDist()->SetHalfX(bounding_box_x);
    single_source->GetPosDist()->SetHalfY(bounding_box_y);
}

/**
 * Called automatically for every event
 */
void GeneratorActionRandomG4::GeneratePrimaries(G4Event* event) {
    // Set event attributes if ROOT-source is used
    if (main_idx <= num_entries - num_parallel)
        GeneratePrimariesRandom();
    else if (!execute_parallel) {
        GeneratePrimariesParallel();
    }

    // Generate event
    particle_source_->GeneratePrimaryVertex(event);

    main_idx++;
}
