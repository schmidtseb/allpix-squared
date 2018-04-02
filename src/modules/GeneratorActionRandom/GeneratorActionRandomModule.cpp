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

#include <G4RunManager.hh>

#include "core/utils/log.h"

#include "GeneratorActionRandomG4.hpp"

using namespace allpix;

GeneratorActionRandomModule::GeneratorActionRandomModule(Configuration config, Messenger* messenger, GeometryManager* geo_manager)
    : Module(std::move(config)), geo_manager_(geo_manager), messenger_(messenger), 
      run_manager_g4_(nullptr) {
}

void GeneratorActionRandomModule::init() {
    run_manager_g4_ = G4RunManager::GetRunManager();
    if(run_manager_g4_ == nullptr) {
        throw ModuleError("Cannot deposit charges using Geant4 without a Geant4 geometry builder");
    }

    run_manager_g4_->SetUserAction(new GeneratorActionRandomG4(config_));
}

void GeneratorActionRandomModule::run(unsigned int) {
}
