/**
 * @file phylo_builder.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Main file for the RACES phylogenetic tree builder
 * @version 0.5
 * @date 2023-10-20
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

    Races::Drivers::PhylogeneticForest 
    build_forest(const std::list<Races::Drivers::CellId>& sample, const bool quiet=true) const
    {
        using namespace Races::Drivers;

        Simulation::BinaryLogger::CellReader reader(drivers_directory);

        auto simulation = load_drivers_simulation(snapshot_path, quiet);
    
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
        ;

        po::options_description program_options;
        for (const auto& [section_name, section]: visible_options) {
            program_options.add(section);
        }
        program_options.add(hidden_options);

        positional_options.add("simulation", 1);

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
            print_help_and_exit("Missing driver simulation directory!", 1);
        }
    
        snapshot_path = get_last_snapshot_path(drivers_directory, "simulation");
    }

    void run() const
    {
        std::list<Races::Drivers::CellId> sample;
        try {
            sample = Races::Drivers::Simulation::BasicLogger::load_sampled_ids(drivers_directory);
        } catch (...) {
            print_help_and_exit("\""+std::string(drivers_directory)+"\" does not contain a list "
                                + "of sampled cell. Produce it by using \"tissue_sampler\".", 1);
        }

        try {
            auto forest = build_forest(sample);

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
