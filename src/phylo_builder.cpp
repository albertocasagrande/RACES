/**
 * @file phylo_builder.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Main file for the simulator
 * @version 0.2
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
#include <string>
#include <filesystem>

#include <boost/program_options.hpp>

#include "common.hpp"

#include "phylogenetic_forest.hpp"

#include "phyloXML.hpp"


class PhylogeneticBuilder : public BasicExecutable
{
    std::filesystem::path snapshot_path;
    std::filesystem::path drivers_directory;
    std::filesystem::path config;

    Races::Drivers::PhylogeneticForest 
    build_forest(const nlohmann::json& simulation_cfg, const bool quiet=true) const
    {
        using namespace Races::Drivers;
        using ConfigReader = Races::Passengers::ConfigReader;

        if (!simulation_cfg.contains("sample region")) {
            throw std::runtime_error("The passengers simulation file must contain "
                                        "a \"sample region\" field");
        }

        const auto sampler_region = ConfigReader::get_sample_region(simulation_cfg["sample region"]);

        Simulation::BinaryLogger::CellReader reader(drivers_directory);

        auto simulation = load_drivers_simulation(snapshot_path, quiet);

        auto sample = simulation.sample_tissue(sampler_region);

        auto genotypes = simulation.tissue().get_genotypes();

        return grow_forest_from(sample, reader, genotypes);
    }

public:

    PhylogeneticBuilder(int argc, char* argv[]):
        BasicExecutable(argv[0], {{"options", "Options"}})
    {
        namespace po = boost::program_options;

        (void)argc;

        visible_options.at("options").add_options()
            ("help,h", "get the help");

        hidden_options.add_options()
            ("simulation", po::value<std::filesystem::path>(&drivers_directory), 
            "the simulation directory")
            ("passenger-config", po::value<std::filesystem::path>(&config), 
            "the passenger simulation config")
        ;

        po::options_description program_options;
        for (const auto& [section_name, section]: visible_options) {
            program_options.add(section);
        }
        program_options.add(hidden_options);

        positional_options.add("simulation", 1);
        positional_options.add("passenger-config", 1);

        po::variables_map vm;
        try {
            po::command_line_parser parser{argc, argv};
            po::store(parser.options(program_options).positional(positional_options).run(), vm);
            po::notify(vm);
        } catch (boost::wrapexcept<boost::program_options::unknown_option> &except) {
            print_help_and_exit(except.what(), 1);
        }

        if (vm.count("help")) {
            print_help_and_exit(0);
        }

        if (!vm.count("simulation")) {
            print_help_and_exit("Missing simulation snapshot!", 1);
        }
    
        snapshot_path = get_last_snapshot_path(drivers_directory, "simulation");

        if (!vm.count("passenger-config")) {
            print_help_and_exit("Missing passenger configuration file!", 1);
        }
    }

    void run() const
    {
        try {
            std::ifstream config_stream(config);
            nlohmann::json simulation_cfg = nlohmann::json::parse(config_stream);

            auto forest = build_forest(simulation_cfg);

            Races::Drivers::IO::phyloXMLStream os;

            os << forest;
        } catch (std::runtime_error& ex) {
            print_help_and_exit(ex.what(), 1);
        }
    }
};

int main(int argc, char* argv[])
{
    PhylogeneticBuilder pbuilder(argc, argv);

    pbuilder.run();

    return 0;
}
