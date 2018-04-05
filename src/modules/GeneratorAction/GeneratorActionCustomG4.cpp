/**
 * @file
 * @brief Implements the particle generator
 * @remark Based on code from John Idarraga
 * @copyright Copyright (c) 2017 CERN and the Allpix Squared authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#include "GeneratorActionCustomG4.hpp"
#include "Math/Point3D.h"

#include <limits>
#include <memory>
#include <Eigen/Core>

#include <G4Event.hh>
#include <G4GeneralParticleSource.hh>
#include <G4ParticleDefinition.hh>
#include <G4ParticleTable.hh>
#include <G4IonTable.hh>
#include <G4UImanager.hh>

#include "core/config/exceptions.h"
#include "core/utils/log.h"
#include "tools/geant4.h"

using namespace allpix;

GeneratorActionCustomG4::GeneratorActionCustomG4(const Configuration& config)
    : particle_source_(std::make_unique<G4GeneralParticleSource>()) {
    // === Load from config ===
    // Source and phantom positions
    G4ThreeVector beam_pos_conf = config.get<G4ThreeVector>("beam_position");
    G4ThreeVector phant_pos_conf = config.get<G4ThreeVector>("phantom_position");

    source_position = Eigen::Vector3d(beam_pos_conf.x(), beam_pos_conf.y(), beam_pos_conf.z());
    phantom_position = Eigen::Vector3d(phant_pos_conf.x(), phant_pos_conf.y(), phant_pos_conf.z());

    // Half length of phantom
    phantom_size = config.get<double>("phantom_size", Units::get(15.0, "mm"));
    // === Load End ===

    // Initialize random stuff
    InitRandom();

    // Set verbosity of source to off
    particle_source_->SetVerbosity(0);

    // Get source specific parameters
    auto single_source = particle_source_->GetCurrentSource();
    // Set global parameters of the source
    single_source->SetNumberOfParticles(1);
    // Set the primary track's start time in for the current event to zero:
    single_source->SetParticleTime(0.0);

    // Get angle distribution parameter
    ang_dist_type = config.get<std::string>("ang_dist_type", "beam2d");


    if(config.has("gps")) {
        LOG(INFO) << "Using gps source";

        // Get UI manager for sending commands
        G4UImanager* ui_g4 = G4UImanager::GetUIpointer();

        auto macro_infile = config.get<std::string>("gps");
        if(config.has("macro_path")) {
            std::string macro_path = config.get<std::string>("macro_path");
            ui_g4->ApplyCommand("/control/macroPath " + macro_path);
            LOG(INFO) << "Set macro path to: " << macro_path << "\"";
        }
        ui_g4->ApplyCommand("/control/execute " + macro_infile);
    } else {
        // Find Geant4 particle
        auto particle_type = config.get<std::string>("particle_type", "");
        std::transform(particle_type.begin(), particle_type.end(), particle_type.begin(), ::tolower);
        auto particle_code = config.get<int>("particle_code", 0);
        G4ParticleDefinition* particle = nullptr;

        if(particle_type == "ion") {
            auto pdg_table = G4IonTable::GetIonTable();

            auto ion_type = config.get<G4TwoVector>("ion_type", G4TwoVector(0., 0.));

            int fAtomicNumber = static_cast<int>(ion_type.x());
            int fAtomicMass = static_cast<int>(ion_type.y());
            double fIonExciteEnergy = 0.0;
            // double fIonCharge = 0.0;

            if(fAtomicNumber == 0 && fAtomicMass == 0) {
                throw InvalidValueError(config, "ion_type", "Particle ion selected but no ion type specified.");
            }

            particle = pdg_table->GetIon(fAtomicNumber, fAtomicMass, fIonExciteEnergy);
            if(particle == nullptr) {
                throw InvalidValueError(config, "ion_type", "Specified ion is not defined.");
            }

            // single_source->SetParticleCharge(fIonCharge*eplus);
        } else {
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
        }  

        LOG(DEBUG) << "Using particle " << particle->GetParticleName() << " (ID " << particle->GetPDGEncoding() << ").";

        single_source->SetParticleDefinition(particle);

        // Set position parameters
        auto pos_dist_type = config.get<std::string>("pos_dist_type", "Beam");
        single_source->GetPosDist()->SetPosDisType(pos_dist_type);
        single_source->GetPosDist()->SetBeamSigmaInR(config.get<double>("beam_size", 0));
        single_source->GetPosDist()->SetCentreCoords(config.get<G4ThreeVector>("beam_position"));

        // Set angle distribution parameters
        single_source->GetAngDist()->SetAngDistType(ang_dist_type);
        single_source->GetAngDist()->DefineAngRefAxes("angref1", G4ThreeVector(-1., 0, 0));
        G4TwoVector divergence = config.get<G4TwoVector>("beam_divergence", G4TwoVector(0., 0.));
        single_source->GetAngDist()->SetBeamSigmaInAngX(divergence.x());
        single_source->GetAngDist()->SetBeamSigmaInAngY(divergence.y());
        G4ThreeVector direction = config.get<G4ThreeVector>("beam_direction");
        if(fabs(direction.mag() - 1.0) > std::numeric_limits<double>::epsilon()) {
            LOG(WARNING) << "Momentum direction is not a unit vector: magnitude is ignored";
        }
        single_source->GetAngDist()->SetParticleMomentumDirection(direction);

        // Set energy parameters
        if(particle_type == "ion") {
            single_source->GetEneDist()->SetEnergyDisType("Mono");
            single_source->GetEneDist()->SetMonoEnergy(0.);
        } else {
            single_source->GetEneDist()->SetEnergyDisType("Gauss");
            single_source->GetEneDist()->SetMonoEnergy(config.get<double>("beam_energy"));
            single_source->GetEneDist()->SetBeamSigmaInE(config.get<double>("beam_energy_spread", 0.));
        }
    }
}

void GeneratorActionCustomG4::InitRandom() {
    // Vector of source-phantom
    Eigen::Vector3d v = phantom_position - source_position;
    LOG(INFO) << phantom_position;
    LOG(INFO) << source_position;

    // Convert to spherical coordinates
    r_p = v.norm();
    phi_p = std::atan2(v.x(), v.z());
    theta_p = std::acos(v.y() / r_p);

    // Get angular range
    theta_rand = std::atan(phantom_size/r_p);
    phi_rand = theta_rand;

    // Create uniform random generator
    uniform_distribution = std::uniform_real_distribution<double>(0, 1);

    LOG(INFO) << "Spherical coordinates = (" << r_p << ", " << phi_p << ", " << theta_p << ")";
    LOG(INFO) << "Angular range = " << theta_rand;
}

void GeneratorActionCustomG4::GeneratePrimariesPyramid() {
    // Generate random
    // Phi
    double randPhi = uniform_distribution(random_generator);
    double minPhi = phi_p - phi_rand;
    double maxPhi = phi_p + phi_rand;
    double phi = minPhi + (maxPhi - minPhi) * randPhi;

    // Theta
    double randTheta = uniform_distribution(random_generator);
    double minTheta = theta_p - theta_rand;
    double maxTheta = theta_p + theta_rand;
    double cosTheta = -randTheta * (std::cos(minTheta) - std::cos(maxTheta)) + std::cos(minTheta);
    double sinTheta = std::sqrt(1 - cosTheta*cosTheta);

    // Transform spherical to cartesian coordinates
    double z = sinTheta*std::cos(phi);
    double x = -sinTheta*std::sin(phi);
    double y = -cosTheta;

    // Set beam direction
    auto single_source = particle_source_->GetCurrentSource();
    // single_source->GetPosDist()->SetCentreCoords(G4ThreeVector(source_position.x(), source_position.y(), source_position.z()));
    single_source->GetAngDist()->SetAngDistType("planar");
    single_source->GetAngDist()->SetParticleMomentumDirection(G4ThreeVector(x, y, z));
    // single_source->SetParticlePosition(G4ThreeVector(source_position.x(), source_position.y(), source_position.z()));

    LOG(DEBUG) << "Angle data: \"" << minPhi << ", " << maxPhi << ", " << minTheta << ", " << maxTheta;
    // LOG(INFO) << "" << ;
    LOG(DEBUG) << "Particle data (x, y, z) = (" << x << ", " << y << ", " << z << ")";
}

/**
 * Called automatically for every event
 */
void GeneratorActionCustomG4::GeneratePrimaries(G4Event* event) {
    if (ang_dist_type == "pyramid") {
        GeneratePrimariesPyramid();
    }

    // Generate event
    particle_source_->GeneratePrimaryVertex(event);
}
