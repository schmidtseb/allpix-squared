/**
 * @file
 * @brief Definition of Monte-Carlo particle object
 * @copyright Copyright (c) 2017 CERN and the Allpix Squared authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#ifndef ALLPIX_MC_PARTICLE_H
#define ALLPIX_MC_PARTICLE_H

#include <Math/Point3D.h>
#include <TRef.h>

#include "Object.hpp"

namespace allpix {
    /**
     * @brief Monte-Carlo particle through the sensor
     */
    class MCParticle : public Object {
    public:
        /**
         * @brief Construct a Monte-Carlo particle
         * @param local_start_point Entry point of the particle in the sensor in local coordinates
         * @param global_start_point Entry point of the particle in the sensor in global coordinates
         * @param local_end_point Exit point of the particle in the sensor in local coordinates
         * @param global_end_point Exit point of the particle in the sensor in global coordinates
         * @param particle_id Identifier for the particle type
         */
        MCParticle(ROOT::Math::XYZPoint local_start_point,
                   ROOT::Math::XYZPoint global_start_point,
                   ROOT::Math::XYZPoint local_end_point,
                   ROOT::Math::XYZPoint global_end_point,
                   int particle_id);

        /**
         * @brief Get the entry point of the particle in local coordinates
         * @return Particle entry point
         */
        ROOT::Math::XYZPoint getLocalStartPoint() const;
        /**
         * @brief Get the entry point of the particle in global coordinates
         * @return Particle entry point
         */
        ROOT::Math::XYZPoint getGlobalStartPoint() const;

        /**
         * @brief Get the exit point of the particle in local coordinates
         * @return Particle exit point
         */
        ROOT::Math::XYZPoint getLocalEndPoint() const;
        /**
         * @brief Get the entry point of the particle in global coordinates
         * @return Particle entry point
         */
        ROOT::Math::XYZPoint getGlobalEndPoint() const;

        /**
         * @brief Get particle identifier
         * @return Particle identifier
         */
        int getParticleID() const;

        /**
         * @brief Set the Monte-Carlo particle
         * @param mc_particle The Monte-Carlo particle
         * @warning Special method because parent can only be set after creation, should not be replaced later.
         */
        void setParent(const MCParticle* mc_particle);
        /**
         * @brief Get the parent MCParticle if it has one
         * @return Parent MCParticle or null pointer if it has no parent
         * @warning No \ref MissingReferenceException is thrown, because a particle without parent should always be handled.
         */
        const MCParticle* getParent() const;

        /**
         * @brief ROOT class definition
         */
        ClassDef(MCParticle, 3);
        /**
         * @brief Default constructor for ROOT I/O
         */
        MCParticle() = default;

    private:
        ROOT::Math::XYZPoint local_start_point_{};
        ROOT::Math::XYZPoint global_start_point_{};
        ROOT::Math::XYZPoint local_end_point_{};
        ROOT::Math::XYZPoint global_end_point_{};

        int particle_id_{};

        TRef parent_;
    };

    /**
     * @brief Typedef for message carrying MC particles
     */
    using MCParticleMessage = Message<MCParticle>;
} // namespace allpix

#endif
