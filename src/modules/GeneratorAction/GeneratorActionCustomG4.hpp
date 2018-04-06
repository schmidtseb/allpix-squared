/**
 * @file
 * @brief Defines the particle generator
 * @copyright Copyright (c) 2017 CERN and the Allpix Squared authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#ifndef ALLPIX_SIMPLE_DEPOSITION_MODULE_GENERATOR_ACTION_CUSTOM_H
#define ALLPIX_SIMPLE_DEPOSITION_MODULE_GENERATOR_ACTION_CUSTOM_H

#include <memory>
#include <random>
#include <Eigen/Core>

#include <G4GeneralParticleSource.hh>
#include <G4ParticleDefinition.hh>
#include <G4SDManager.hh>
#include <G4ThreeVector.hh>
#include <G4TwoVector.hh>
#include <G4VUserPrimaryGeneratorAction.hh>

#include "core/config/Configuration.hpp"

class G4RadioactiveDecay;

namespace allpix {
    /**
     * @brief Generates the particles in every event
     */
    class GeneratorActionCustomG4 : public G4VUserPrimaryGeneratorAction {
    public:
        /**
         * @brief Constructs the generator action
         * @param config Configuration of the \ref DepositionGeant4Module module
         */
        explicit GeneratorActionCustomG4(const Configuration& config);

        /**
         * @brief Generate the particle for every event
         */
        void GeneratePrimaries(G4Event*) override;

    private:
        std::unique_ptr<G4GeneralParticleSource> particle_source_;
        // Decay physics
        bool use_decay_physics;
        G4RadioactiveDecay* radioactiveDecay;
        void ActivateDecay();

        // Defined in config
        std::string ang_dist_type;
        Eigen::Vector3d source_position;
        Eigen::Vector3d phantom_position;
        double phantom_size;

        // Random generator
		void InitRandom();
        void GeneratePrimariesPyramid();
        std::mt19937_64 random_generator;
        double r_p, phi_p, theta_p;
        double theta_rand, phi_rand;
        std::uniform_real_distribution<double> uniform_distribution;
    };
} // namespace allpix

#endif
