/**
 * @file common.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements auxiliary classes and functions for executables
 * @version 0.3
 * @date 2023-10-16
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
#include <regex>

#include "common.hpp"

void BasicExecutable::print_help_and_exit(const std::string& message, const int& err_code) const
{ 
    std::ostream& os = (err_code==0 ? std::cout: std::cerr);
    
    os << message << std::endl << std::endl;
    
    print_help_and_exit(err_code);
}

void BasicExecutable::print_help_and_exit(const int& err_code) const
{ 
    std::ostream& os = (err_code==0 ? std::cout: std::cerr);
    
    os << "Syntax: " << program_name;
    
    for (unsigned int i=0;i<positional_options.max_total_count(); ++i) {
        os << " <" << positional_options.name_for_position(i) << ">"; 
    }

    os << std::endl << std::endl;

    for (const auto& [name,option]: visible_options) {
        os << option << std::endl;
    }

    exit(err_code);
}

BasicExecutable::BasicExecutable(const std::string& program_name,
                                 const std::vector<std::pair<std::string, std::string>>& option_sections): 
    program_name{program_name}, visible_options{}, hidden_options("Hidden options")
{
    for (const auto& [section_name, section_description]: option_sections) {
        visible_options.insert({section_name, {section_description}});
    }
}


std::filesystem::path 
BasicExecutable::get_last_snapshot_path(const std::string& simulation_dir,
                                        const std::string& parameter_name)
{
    if (!std::filesystem::is_directory(simulation_dir)) {
        print_help_and_exit("The \"" + parameter_name + "\" parameter must be a directory", 1);
    }

    try {
        return find_last_snapshot(simulation_dir);
    } catch (std::exception &except) {
        print_help_and_exit(except.what(), 1);

        exit(1);
    }
}

std::filesystem::path BasicExecutable::find_last_snapshot(const std::filesystem::path directory)
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
BasicExecutable::load_drivers_simulation(const std::filesystem::path snapshot_path, const bool quiet)
{
    Races::Archive::Binary::In archive(snapshot_path);

    Races::Drivers::Simulation::Simulation simulation;

    if (quiet) {
        archive & simulation;
    } else {
        archive.load(simulation, "simulation");
    }

    return simulation;
}
