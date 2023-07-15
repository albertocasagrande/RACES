/**
 * @file simulator_main.cpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Main file for the simulator
 * @version 0.11
 * @date 2023-07-15
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
#include <chrono>

#include <boost/program_options.hpp>

#include "simulation.hpp"

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
        ("frames-per-second,f", po::value<unsigned int>()->default_value(5), 
         "the number of frames per second")
#endif
        ("seed,s", po::value<int>()->default_value(0), "random generator seed")
        ("no-logging,n", "avoid logging");
    
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

    SomaticGenotype A("A",{{0.01,0.01}});
    A["-"].set_rates({{CellEventType::DIE, 0.1},
                      {CellEventType::DUPLICATE, 0.3},
                      {CellEventType::PASSENGER_MUTATION, 20}});
    A["+"].set_rates({{CellEventType::DIE, 0.1},
                      {CellEventType::DUPLICATE, 0.45},
                      {CellEventType::PASSENGER_MUTATION, 2}});

    SomaticGenotype B("B",{{0.01,0.01}});
    B["-"].set_rates({{CellEventType::DIE, 0.1},
                      {CellEventType::DUPLICATE, 0.2},
                      {CellEventType::PASSENGER_MUTATION, 1}});
    B["+"].set_rates({{CellEventType::DIE, 0.01},
                      {CellEventType::DUPLICATE, 0.02},
                      {CellEventType::PASSENGER_MUTATION, 20}});

    Simulation simulation(vm.count("seed"));

    simulation.set_tissue("Liver", {1000,1000})
              .add_species(A)
              .add_species(B)
              .add_cell(B["-"], {500, 500})
              .add_somatic_mutation(B,A,75);

    simulation.death_activation_level = 100;

    using namespace std::chrono_literals;
    simulation.set_interval_between_snapshots(5min);

    const bool logging = (vm.count("no-logging")==0);

#ifdef WITH_SDL2
    if (vm.count("plot")) {

        UI::TissuePlotter<UI::SDLWindow> plotter(simulation.tissue());

        plotter.set_frames_per_second(vm["frames-per-second"].as<unsigned int>());

        simulation.run_up_to(time_horizon, plotter);

        while (plotter.waiting_end()) {
            plotter.plot(simulation.get_statistics());
        }
    } else {
        simulation.run_up_to(time_horizon, logging);
    }
#else
    simulation.run_up_to(time_horizon, logging);
#endif

    return 0;
}
