/**
 * @file races_phylo.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Main file for the simulator
 * @version 0.7
 * @date 2023-10-14
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
#include <filesystem>
#include <regex>

#include <boost/program_options.hpp>

#include "common.hpp"

#include "simulation.hpp"

#include "phylogenetic_forest.hpp"

#include "phyloXML.hpp"

#include "json_config.hpp"


std::ostream& print_help(std::ostream& os, const std::string& program_name, 
                         const boost::program_options::options_description& options)
{
   os << "Syntax: "<< program_name <<" <directory> <passenger config>" << std::endl
      << options << std::endl;

   return os;   
}

Races::Drivers::PhylogeneticForest 
build_forest(const std::string& driver_directory, const nlohmann::json& simulation_cfg,
             const bool quiet=true)
{
    using namespace Races::Drivers;
    using ConfigReader = Races::Passengers::ConfigReader;

    if (!simulation_cfg.contains("sample region")) {
        throw std::runtime_error("The passengers simulation file must contain "
                                    "a \"sample region\" field");
    }

    const auto sampler_region = ConfigReader::get_sample_region(simulation_cfg["sample region"]);

    Simulation::BinaryLogger::CellReader reader(driver_directory);

    auto simulation = load_drivers_simulation(driver_directory, quiet);

    auto sample = simulation.sample_tissue(sampler_region);

    auto genotypes = simulation.tissue().get_genotypes();

    return grow_forest_from(sample, reader, genotypes);
}

int main(int argc, char* argv[])
{
    using namespace Races::Drivers;
    namespace po = boost::program_options;

    po::options_description visible("Options");
    visible.add_options()
        ("help,h", "get the help");
    
    std::string directory;
    std::string config;

    po::options_description hidden("Hidden options");
    hidden.add_options()
        ("directory", po::value<std::string>(&directory), 
         "the simulation log directory")
        ("passenger-config", po::value<std::string>(&config), 
         "the passenger simulation config")
    ;

    po::options_description generic;
    generic.add(visible).add(hidden);

    po::positional_options_description p;
    p.add("directory", 1);
    p.add("passenger-config", 1);

    po::variables_map vm;
    try {
        po::command_line_parser parser{argc, argv};
        po::store(parser.options(generic).positional(p).run(), vm);
        po::notify(vm);
    } catch (boost::wrapexcept<boost::program_options::unknown_option> &except) {
        std::cout << except.what() << std::endl;
        print_help(std::cout, argv[0], visible);
        exit(1);
    }

    if (vm.count("help")) {
        print_help(std::cout, argv[0], visible);
        return 1;
    }

    if (!vm.count("directory")) {
        std::cerr << "Missing simulation log directory!" << std::endl << std::endl;
        print_help(std::cout, argv[0], visible);
        return 1;
    }

    if (!vm.count("passenger-config")) {
        std::cerr << "Missing passenger configuration file!" << std::endl << std::endl;
        print_help(std::cout, argv[0], visible);
        return 1;
    }

    std::ifstream config_stream(config);
    nlohmann::json simulation_cfg = nlohmann::json::parse(config_stream);

    try {
        auto forest = build_forest(directory, simulation_cfg);

        IO::phyloXMLStream os;

        os << forest;
    } catch (std::runtime_error& ex) {
        std::cerr << ex.what() << std::endl << std::endl;
        print_help(std::cout, argv[0], visible); 

        return 1;
    }

    return 0;
}
