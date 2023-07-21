/**
 * @file races_sim.cpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Main file for the simulator
 * @version 0.13
 * @date 2023-07-21
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
#include <csignal>

#include <boost/program_options.hpp>

#include "simulation.hpp"
#include "progress_bar.hpp"

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


Races::Drivers::Simulation::Simulation simulation;
Races::UI::ProgressBar *bar;

void termination_handling(int signal_num)
{
    if (bar != nullptr) {
        simulation.make_snapshot(bar);
    }

    if (bar != nullptr) {
        bar->set_message("Aborted");

        delete bar;
    }

    exit(signal_num);
}

int main(int argc, char* argv[])
{
    using namespace Races::Drivers;
    namespace po = boost::program_options;
    namespace UI = Races::UI;

    po::options_description visible("Options");
    visible.add_options()
        ("help,h", "get the help")
#ifdef WITH_INDICATORS
        ("quit,q", "disable progress bar")
#endif // WITH_INDICATORS
#ifdef WITH_SDL2
        ("plot,p", 
         "plot a graphical representation of the simulation")
        ("frames-per-second,f", po::value<unsigned int>()->default_value(5), 
         "the number of frames per second")
#endif
        ("seed,s", po::value<int>()->default_value(0), "random generator seed")
        ("no-logging,n", "disable logging");
    
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

    Genotype A("A",{{0.01,0.01}});
    A["-"].set_rates({{CellEventType::DIE, 0.1},
                      {CellEventType::DUPLICATE, 0.3}});
    A["+"].set_rates({{CellEventType::DIE, 0.1},
                      {CellEventType::DUPLICATE, 0.45}});

    Genotype B("B",{{0.01,0.01}});
    B["-"].set_rates({{CellEventType::DIE, 0.1},
                      {CellEventType::DUPLICATE, 0.2}});
    B["+"].set_rates({{CellEventType::DIE, 0.01},
                      {CellEventType::DUPLICATE, 0.02}});

    std::signal(SIGINT, termination_handling);

    if (vm.count("quiet")==0) {
        UI::ProgressBar::hide_console_cursor();

        bar = new UI::ProgressBar();
        bar->set_message("Creating tissue");
    }

    using namespace std::chrono_literals;

    simulation.set_tissue("Liver", {1000,1000})
              .add_species(A)
              .add_species(B)
              .add_cell(B["-"], {500, 500})
              .add_genomic_mutation(B,A,75)
              .random_generator_seed(vm.count("seed"))
              .set_interval_between_snapshots(5min);

    simulation.death_activation_level = 100;

    const bool logging = (vm.count("no-logging")==0);

#ifdef WITH_SDL2
    if (vm.count("plot")) {

        UI::TissuePlotter<UI::SDLWindow> plotter(simulation.tissue());

        plotter.set_frames_per_second(vm["frames-per-second"].as<unsigned int>());

        if (bar != nullptr) {
            simulation.run_up_to(time_horizon, plotter, *bar);
        } else {
            simulation.run_up_to(time_horizon, plotter);
        }

        while (plotter.waiting_end()) {
            plotter.plot(simulation.get_statistics());
        }
    } else {
#endif // WITH_SDL2

    if (bar != nullptr) {
        simulation.run_up_to(time_horizon, *bar, logging);
    } else {
        simulation.run_up_to(time_horizon, logging);
    }

#ifdef WITH_SDL2
    }
#endif // WITH_SDL2

    UI::ProgressBar::show_console_cursor();

    if (bar != nullptr) {
        delete bar;
    }
    
    return 0;
}
