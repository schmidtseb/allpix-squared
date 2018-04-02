/**
 * @file
 * @brief Implementation of [GeneratorAction] module
 * @copyright Copyright (c) 2017 CERN and the Allpix Squared authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#include "GeneratorActionModule.hpp"

#include <string>
#include <utility>

#include <G4RunManager.hh>

#include "core/utils/log.h"

#include "GeneratorActionG4.hpp"
#include "GeneratorActionCustomG4.hpp"

using namespace allpix;

GeneratorActionModule::GeneratorActionModule(Configuration config, Messenger* messenger, GeometryManager* geo_manager)
    : Module(std::move(config)), geo_manager_(geo_manager), messenger_(messenger),
      run_manager_g4_(nullptr) {

    config_.setDefault<bool>("random", false);
    config_.setDefault<bool>("custom", false);
}

void GeneratorActionModule::init() {
    run_manager_g4_ = G4RunManager::GetRunManager();
    if(run_manager_g4_ == nullptr) {
        throw ModuleError("Cannot deposit charges using Geant4 without a Geant4 geometry builder");
    }

    if(config_.get<bool>("custom")) {
        run_manager_g4_->SetUserAction(new GeneratorActionCustomG4(config_));
    } else {
        run_manager_g4_->SetUserAction(new GeneratorActionG4(config_));
    }
}

void GeneratorActionModule::run(unsigned int) {
}
