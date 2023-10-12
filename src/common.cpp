/**
 * @file common.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements auxiliary classes and functions for executables
 * @version 0.1
 * @date 2023-10-12
 * 
 * @copyright Copyright (c) 2023
 * 
 * MIT License
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */


#include <iostream> 
#include <string>
#include <regex>

#include "common.hpp"

std::filesystem::path find_last_snapshot(const std::filesystem::path directory)
{
    namespace fs = std::filesystem;

    if (!fs::exists(directory)) {
        std::ostringstream oss;

        oss << "The path "<< directory<< " does not exist";
        throw std::runtime_error(oss.str());
    }

    if (!fs::is_directory(directory)) {
        std::ostringstream oss;

        oss << directory<< " is not a directory";
        throw std::runtime_error(oss.str());
    }

    std::regex re(std::string(directory / "snapshot_\\d+-\\d+.dat"));

    bool found{false};
    std::string last;
    for (const auto & entry : fs::directory_iterator(directory)) {
        if (entry.is_regular_file()) {
            const std::string e_string = entry.path();
            if (std::regex_match(e_string, re)) {
                if (!found || e_string.compare(last)>0) {
                    last = e_string;
                    found = true;
                }
            }
        }
    }

    if (!found) {
        std::ostringstream oss;

        oss << "No RACES simulation snapshot in "<< directory<< "";
        throw std::runtime_error(oss.str());
    }

    return last;
}

Races::Drivers::Simulation::Simulation 
load_drivers_simulation(const std::filesystem::path driver_directory, const bool quiet)
{
    auto last_snapshot_path = find_last_snapshot(driver_directory);

    Races::Archive::Binary::In archive(last_snapshot_path);

    Races::Drivers::Simulation::Simulation simulation;

    if (quiet) {
        archive & simulation;
    } else {
        archive.load(simulation, "simulation");
    }

    return simulation;
}
