/**
 * @file
 * @brief Definition of [GeneratorAction] module
 * @copyright Copyright (c) 2017 CERN and the Allpix Squared authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 *
 * Contains minimal dummy module to use as a start for the development of your own module
 *
 * Refer to the User's Manual for more details.
 */

#include <string>

#include "core/config/Configuration.hpp"
#include "core/geometry/GeometryManager.hpp"
#include "core/messenger/Messenger.hpp"
#include "core/module/Module.hpp"

class G4RunManager;

namespace allpix {
    /**
     * @ingroup Modules
     * @brief Module to do function
     *
     * More detailed explanation of module
     */
    class GeneratorActionModule : public Module {
    public:
        /**
         * @brief Constructor for this unique module
         * @param config Configuration object for this module as retrieved from the steering file
         * @param messenger Pointer to the messenger object to allow binding to messages on the bus
         * @param geo_manager Pointer to the geometry manager, containing the detectors
         */
        GeneratorActionModule(Configuration config, Messenger* messenger, GeometryManager* geo_manager);

        /**
         * @brief [Initialise this module]
         */
        void init() override;

        /**
         * @brief [Run the function of this module]
         */
        void run(unsigned int) override;

    private:
        // General module members
        GeometryManager* geo_manager_;
        Messenger* messenger_;

        // Pointer to the Geant4 manager (owned by GeometryBuilderGeant4)
        G4RunManager* run_manager_g4_;
    };
} // namespace allpix
