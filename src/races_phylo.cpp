/**
 * @file races_phylo.cpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Main file for the simulator
 * @version 0.4
 * @date 2023-09-17
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

#include "simulation.hpp"

#include "sampler.hpp"
#include "phylogenetic_forest.hpp"

#include "phyloXML.hpp"

std::ostream& print_help(std::ostream& os, const std::string& program_name, 
                         const boost::program_options::options_description& options)
{
   os << "Syntax: "<< program_name <<" <directory>" << std::endl
      << options << std::endl;

   return os;   
}

std::filesystem::path find_last_snapshot(const std::filesystem::path directory)
{
    namespace fs = std::filesystem;

    if (!fs::exists(directory)) {
        std::ostringstream oss;

        oss << "The path "<< directory<< " does not exists";
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
                if (!found || last.compare(e_string)>0) {
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

int main(int argc, char* argv[])
{
    using namespace Races::Drivers;
    namespace po = boost::program_options;

    po::options_description visible("Options");
    visible.add_options()
        ("help,h", "get the help");
    
    std::string directory;

    po::options_description hidden("Hidden options");
    hidden.add_options()
        ("directory", po::value<std::string>(&directory), 
         "the simulation log directory")
    ;

    po::options_description generic;
    generic.add(visible).add(hidden);

    po::positional_options_description p;
    p.add("directory", 1);

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

    std::filesystem::path last_snapshot_path;
    try {
        last_snapshot_path = find_last_snapshot(directory);
    } catch (std::runtime_error& ex) {
        std::cerr << ex.what() << std::endl << std::endl;
        print_help(std::cout, argv[0], visible); 

        return 1;
    }

    Simulation::BinaryLogger::CellReader reader(directory);
    
    Simulation::Simulation simulation;
    {
        Races::Archive::Binary::In archive(last_snapshot_path);

        archive & simulation;
    }

    const auto& tissue = simulation.tissue();

    RectangleSampler sampler(tissue,{496,495},{500,500});

    auto genotypes = tissue.get_genotypes();

    PhylogeneticForest forest = grow_forest_from(sampler, reader, genotypes);

    IO::phyloXMLStream os;

    os << forest;

    return 0;
}
