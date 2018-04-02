/**
 * @file
 * @brief Defines the particle generator
 * @copyright Copyright (c) 2017 CERN and the Allpix Squared authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#ifndef ALLPIX_SIMPLE_DEPOSITION_MODULE_GENERATOR_ACTION_RANDOM_H
#define ALLPIX_SIMPLE_DEPOSITION_MODULE_GENERATOR_ACTION_RANDOM_H

#include <memory>
#include <random>
#include <Eigen/Core>

#include <G4GeneralParticleSource.hh>
#include <G4ParticleDefinition.hh>
#include <G4SDManager.hh>
#include <G4ThreeVector.hh>
#include <G4TwoVector.hh>
#include <G4VUserPrimaryGeneratorAction.hh>

#include <TFile.h>
#include <TTree.h>

#include "core/config/Configuration.hpp"

class G4SingleParticleSource;

namespace allpix {
    /**
     * @brief Generates the particles in every event
     */
    class GeneratorActionRandomG4 : public G4VUserPrimaryGeneratorAction {
    public:
        /**
         * @brief Constructs the generator action
         * @param config Configuration of the \ref DepositionGeant4Module module
         */
        explicit GeneratorActionRandomG4(const Configuration& config);

        /**
         * @brief Generate the particle for every event
         */
        void GeneratePrimaries(G4Event*) override;

    private:
        std::unique_ptr<G4GeneralParticleSource> particle_source_;

        // Defined in config

        // Root file
        TFile* root_file;
        TTree* root_tree;

        // Number of entries
        long long int num_entries;

        // Attributes
        double E_branch;
        double x_branch, y_branch, z_branch;
        double phi_branch, theta_branch;

        // Source offset
        double x_offset, y_offset, z_offset;

        // z-Inversion
        double z_invert;

        // Source
        G4SingleParticleSource* single_source;

        // Random generator
        void InitRandom(const Configuration& config);
        void GeneratePrimariesRandom();
        std::mt19937_64 random_generator;
        std::uniform_int_distribution<long long int> uniform_distribution;
    };
} // namespace allpix

#endif
