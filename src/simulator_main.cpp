/**
 * @file simulator_main.cpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Main file for the simulator
 * @version 0.2
 * @date 2023-06-28
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
#include <vector>

#include <boost/program_options.hpp>

#include "binary_logger.hpp"
#include "simulator.hpp"

#ifdef WITH_SDL2
#include "SDL_plot.hpp"
#endif

std::ostream& print_help(std::ostream& os, const std::string& program_name, 
                         const boost::program_options::options_description& options)
{
   os << "Syntax: "<< program_name <<" <time horizon>" << std::endl
      << options << std::endl;

   return os;   
}

int main(int argc, char* argv[])
{
    using namespace Races;
    namespace po = boost::program_options;

    po::options_description visible("Options");
    visible.add_options()
        ("help,h", "get the help")
        // ("bar,b", "show a bar for the elapsed simulation time")
#ifdef WITH_SDL2
        ("plot,p", 
         "plot a graphical representation of the simulation")
#endif
    ;
    
    po::options_description hidden("Hidden options");
    hidden.add_options()
        ("time-horizon", po::value<long double>(), 
         "the simulation time horizon")
    ;


    po::options_description generic;
    generic.add(visible).add(hidden);

    po::positional_options_description p;
    p.add("time-horizon", -1);

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

    if (!vm.count("time-horizon")) {
        std::cout << "Missing simulation time horizon!" << std::endl << std::endl;
        print_help(std::cout, argv[0], visible);
        return 1;
    }

    long double time_horizon = vm["time-horizon"].as<long double>();

    std::vector<DriverGenotype> genotypes;

    genotypes.push_back(DriverGenotype("A",{
            {CellEventType::DIE, 0.1},{CellEventType::DUPLICATE, 0.2}}));
    genotypes.push_back(DriverGenotype("A",{
            {CellEventType::DIE, 0.1},{CellEventType::DUPLICATE, 0.3}},true));

    genotypes.push_back(DriverGenotype("B",{
            {CellEventType::DIE, 0.1},{CellEventType::DUPLICATE, 0.2}}));
    genotypes.push_back(DriverGenotype("B",{
            {CellEventType::DIE, 0.01},{CellEventType::DUPLICATE, 0.02}},true));

    Tissue tissue("Liver", 1000,1000);

    for (const auto& gen: genotypes) {
        tissue.add_species(gen);
    }

    tissue.add_driver_epigenetic_mutation(genotypes[0].get_id(),genotypes[1].get_id(), 0.001);
    tissue.add_driver_epigenetic_mutation(genotypes[1].get_id(),genotypes[0].get_id(), 0.001);
    tissue.add_driver_epigenetic_mutation(genotypes[2].get_id(),genotypes[3].get_id(), 0.001);
    tissue.add_driver_epigenetic_mutation(genotypes[3].get_id(),genotypes[2].get_id(), 0.001);

    tissue.add_driver_somatic_mutation(genotypes[0].get_id(),genotypes[2].get_id(), 100);
    tissue.add_driver_somatic_mutation(genotypes[1].get_id(),genotypes[3].get_id(), 100);

    tissue.add(genotypes[0].get_id(), {250, 500});
    tissue.add(genotypes[2].get_id(), {750, 500});

#ifdef WITH_SDL2
    if (vm.count("plot")) {
        BinaryLogger logger(tissue.get_name());

        BasicSimulator<BinaryLogger, UI::SDLWindow> simulator(tissue, &logger);

        simulator.snapshot_interval = 20;

        simulator.run_up_to(time_horizon);
    } else {
#endif
        BinaryLogger logger(tissue.get_name());

        BasicSimulator<BinaryLogger> simulator(tissue, &logger);

        simulator.run_up_to(time_horizon);

#ifdef WITH_SDL2
    }
#endif

    return 0;
}
