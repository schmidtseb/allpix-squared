/**
 * @file
 * @brief Implementation of [GeneratorActionRandom] module
 * @copyright Copyright (c) 2017 CERN and the Allpix Squared authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#include "GeneratorActionRandomModule.hpp"

#include <string>
#include <utility>

#include "core/utils/log.h"

using namespace allpix;

GeneratorActionRandomModule::GeneratorActionRandomModule(Configuration config, Messenger* messenger, GeometryManager* geo_manager)
    : Module(std::move(config)), geo_manager_(geo_manager), messenger_(messenger) {

    // ... Implement ... (Typically bounds the required messages and optionally sets configuration defaults)
    // Input required by this module
    messenger_->bindMulti(this, &GeneratorActionRandomModule::messages_, MsgFlags::REQUIRED);
}

void GeneratorActionRandomModule::init() {
    // Loop over detectors and do something
    std::vector<std::shared_ptr<Detector>> detectors = geo_manager_->getDetectors();
    for(auto& detector : detectors) {
        // Get the detector name
        std::string detectorName = detector->getName();
        LOG(DEBUG) << "Detector with name " << detectorName;
    }
}

void GeneratorActionRandomModule::run(unsigned int) {
    // ... Implement ... (Typically uses the configuration to execute function and outputs an message)
    // Loop through all receieved messages and print some information
    for(auto& message : messages_) {
        std::string detectorName = message->getDetector()->getName();
        LOG(DEBUG) << "Picked up " << message->getData().size() << " objects from detector " << detectorName;
    }
}