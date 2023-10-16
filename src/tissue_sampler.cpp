/**
 * @file tissue_sampler.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Main file for a tool that sample the tissue
 * @version 0.1
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

#include <boost/program_options.hpp>

#include "common.hpp"

#include "simulation.hpp"

class TissueSampler : public BasicExecutable
{
    std::filesystem::path snapshot_path;
    std::filesystem::path drivers_directory;
    std::filesystem::path config;

    bool remove_sample_tissue;

    void perform_sampling(std::ostream& os, const nlohmann::json& sampling_cfg, const bool quiet=true) const
    {
        using namespace Races::Drivers;
        using ConfigReader = Races::Passengers::ConfigReader;

        if (!sampling_cfg.is_array()) {
            throw std::runtime_error("The sampling configuration does not consist in an array");
        }

        Simulation::BinaryLogger::CellReader reader(drivers_directory);

        auto simulation = load_drivers_simulation(snapshot_path, quiet);

        for (const auto& rectangle_json : sampling_cfg) {
            const auto sampler_region = ConfigReader::get_sample_region(rectangle_json);

            std::list<Races::Drivers::CellId> sample; 
            
            if (remove_sample_tissue) {
                sample = simulation.sample_and_remove_tissue(sampler_region);
            } else {
                sample = simulation.sample_tissue(sampler_region);
            }

            os << "# " << simulation.get_time()
               << " "<< sampler_region.lower_corner
               << " "<< sampler_region.upper_corner
               << " " << sample.size() 
               << " " << sampler_region.size() << std::endl;

            for (const auto& cell_id: sample) {
                os << cell_id << std::endl;
            }
        }

        if (remove_sample_tissue) {
            simulation.make_snapshot<Races::UI::ProgressBar>(nullptr);
        }
    }
public:

    TissueSampler(int argc, char* argv[]):
        BasicExecutable(argv[0], {{"options", "Options"}})
    {
        namespace po = boost::program_options;

        (void)argc;

        visible_options.at("options").add_options()
            ("remove-sample,r","remove the sample from the tissue")
            ("help,h", "get the help")
        ;

        hidden_options.add_options()
            ("simulation", po::value<std::filesystem::path>(&snapshot_path), 
            "the simulation directory/snapshot")
            ("sampling-config", po::value<std::filesystem::path>(&config), 
            "the sampling simulation config")
        ;

        po::options_description program_options;
        for (const auto& [section_name, section]: visible_options) {
            program_options.add(section);
        }
        program_options.add(hidden_options);

        positional_options.add("simulation", 1);
        positional_options.add("sampling-config", 1);

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

        if (std::filesystem::is_directory(snapshot_path)) {
            drivers_directory = snapshot_path;
            try {
                snapshot_path = find_last_snapshot(snapshot_path);
            } catch (...) {
                print_help_and_exit("\"" + std::string(snapshot_path)
                                    + "\" does not contain any snapshot", 1);
            }
        } else {
            drivers_directory = snapshot_path.parent_path();
        }

        if (!vm.count("sampling-config")) {
            print_help_and_exit("Missing passenger configuration file!", 1);
        }

        remove_sample_tissue = (vm.count("remove-sample")>0);
    }

    void run() const
    {
        try {
            std::ifstream config_stream(config);
            nlohmann::json simulation_cfg = nlohmann::json::parse(config_stream);

            perform_sampling(std::cout, simulation_cfg, true);
        } catch (std::runtime_error& ex) {
            print_help_and_exit(ex.what(), 1);
        }
    }
};

int main(int argc, char* argv[])
{
    TissueSampler sampler(argc, argv);

    sampler.run();

    return 0;
}
