/**
 *  @author Simon Spannagel <simon.spannagel@cern.ch>
 *  @author Koen Wolters <koen.wolters@cern.ch>
 */

#include "Configuration.hpp"

#include <fstream>
#include <iostream>
#include <cstdlib>

#include "../utils/exceptions.h"

using namespace allpix;

Configuration::Configuration(std::string name): name_(name) {}

std::string Configuration::getName() const {
    return name_;
}

void Configuration::print(std::ostream &out) const {
    for (auto iter = config_.begin(); iter != config_.end(); ++iter) {
        out << iter->first << " : " << iter->second << std::endl;
    }
}

void Configuration::print() const { 
    print(std::cout); 
}

std::string Configuration::get_string(const std::string &key) const {
    config_t::const_iterator iter = config_.find(key);
    if (iter != config_.end()) {
        return iter->second;
    }
    throw allpix::exception("Configuration: key not found");
}

void Configuration::set_string(const std::string &key,
                            const std::string &val) {
    config_[key] = val;
}
