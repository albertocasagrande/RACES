/**
 * @file drivers_sim.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Main file for the driver simulator
 * @version 0.7
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
#include <vector>
#include <chrono>
#include <csignal>
#include <fstream>
#include <regex>

#include <boost/program_options.hpp>
#include <nlohmann/json.hpp>

#include "common.hpp"
#include "simulation.hpp"
#include "progress_bar.hpp"

#ifdef WITH_SDL2
#include "SDL_plot.hpp"
#endif


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

class DriverSimulator : public BasicExecutable
{
    std::filesystem::path simulation_filename;
    std::string snapshot_path;
    long double time_horizon;
    unsigned int frames_per_second;
    bool plot;
    bool quiet;
    bool disable_storage;
    bool duplicate_internal_cells;

    static double get_rate(const nlohmann::json& rate_json)
    {
        double rate = rate_json.template get<double>();

        if (rate < 0) {
            throw std::domain_error("the rate must be non-negative");
        }

        return rate;
    }

    static double get_epigenetic_rate(const nlohmann::json& rate_json,
                                      const std::vector<std::string>& field_aliases)
    {
        for (const auto& field_name : field_aliases) {
            if (rate_json.contains(field_name)) {
                return get_rate(rate_json[field_name]);
            }
        }

        std::ostringstream oss;
        auto alias_it = std::begin(field_aliases);

        oss << "Every \"epigenetic rates\" element must contain a \""
            << *alias_it;
        while (++alias_it != std::end(field_aliases)) {
            oss << "/" << *alias_it;
        }
        oss << "\" field";

        throw std::domain_error(oss.str());
    }

    static Races::Drivers::Genotype create_genotype(const nlohmann::json& genotype_json)
    {
        using namespace Races::Drivers;
        std::vector<EpigeneticRates> epigenetic_rates;

        if (!genotype_json.is_object()) {
            throw std::domain_error("The genotype specification must be an object");
        }

        if (!genotype_json.contains("epigenetic rates")) {
            throw std::domain_error("The genotype specification must contain a \"epigenetic rates\" field");
        }

        if (!genotype_json["epigenetic rates"].is_array()) {
            throw std::domain_error("The \"epigenetic rates\" field must be an array of epigenetic rates");
        }

        if (!genotype_json.contains("epigenetic types")) {
            throw std::domain_error("The genotype specification must contain a \"epigenetic types\" field");
        }

        if (!genotype_json["epigenetic types"].is_array()) {
            throw std::domain_error("The \"epigenetic types\" field must be an array of epigenetic types");
        }

        if (static_cast<size_t>(1<<genotype_json["epigenetic rates"].size())!=genotype_json["epigenetic types"].size()) {
            throw std::domain_error("The epigenetic types must be exponential in number with respect "
                                    "to the number of promotor epigenetic rates");
        }

        for (const auto& rate_json : genotype_json["epigenetic rates"]) {

            if (!rate_json.is_object()) {
                throw std::domain_error("Every \"epigenetic rates\" must be an object");
            }

            const auto methylation_rate = get_epigenetic_rate(rate_json, 
                                                              {"methylation", "on"});
            const auto demethylation_rate = get_epigenetic_rate(rate_json, 
                                                                {"demethylation", "off"});
            epigenetic_rates.push_back({methylation_rate, demethylation_rate});
        }

        if (!genotype_json.contains("name")) {
            throw std::domain_error("The genotype specification must contain a \"name\" field");
        }

        Genotype genotype(genotype_json["name"].template get<std::string>(), epigenetic_rates);

        for (const auto& type_json : genotype_json["epigenetic types"]) {
            auto& type = genotype[type_json["status"].template get<std::string>()];
            type.set_rate(CellEventType::DIE, get_rate(type_json["death rate"]));
            type.set_rate(CellEventType::DUPLICATE, get_rate(type_json["duplication rate"]));
        }

        return genotype;
    }

    static void configure_tissue(const nlohmann::json& tissue_json)
    {
        using namespace Races::Drivers::Simulation;

        if (!tissue_json.contains("size")) {
            throw std::domain_error("The tissue specification must contain a \"size\" field");
        }

        if (!tissue_json["size"].is_array()) {
            throw std::domain_error("The tissue specification size must be an array");
        }

        std::vector<AxisSize> sizes;
        for (const auto& value : tissue_json["size"]) {
            sizes.push_back(value.template get<AxisSize>());
        }
        if (sizes.size()==2 || sizes.size()==3) {
            std::string name;
            if (tissue_json.contains("name")) {
                name = tissue_json["name"].template get<std::string>();
            }
            simulation.set_tissue(name, sizes);

            return;
        }

        throw std::domain_error("tissue is either a 2D or a 3D space");
    }

    static Races::Drivers::Simulation::PositionInTissue
    get_position(const nlohmann::json& position_json)
    {
        using namespace Races::Drivers::Simulation;

        if (!position_json.is_array()) {
            throw std::domain_error("The \"position\" field must be an array of Natural values");
        }
        std::vector<AxisPosition> position;
        for (const auto& value : position_json) {
            position.push_back(value.template get<AxisPosition>());
        }

        if (position.size()==2) {
            return {position[0], position[1]};
        }

        if (position.size()==3) {
            return {position[0], position[1], position[2]};
        }

        throw std::domain_error("tissue may be either a 2D or a 3D space");
    }

    static void configure_initial_cells(const nlohmann::json& initial_cells_json,
                                        const std::map<std::string, Races::Drivers::Genotype> genotypes)
    {
        if (!initial_cells_json.is_array()) {
            throw std::domain_error("The \"initial cells\" field must contain an array");
        }
        for (const auto& initial_cell_json : initial_cells_json) {

            if (!initial_cell_json.is_object()) {
                throw std::domain_error("Every element in the \"initial cells\" field must be an object");
            }

            if (!initial_cell_json.contains("genotype")) {
                throw std::domain_error("Every initial cell must contain a \"genotope\" field");
            }

            if (!initial_cell_json.contains("genotype")) {
                throw std::domain_error("Every initial cell must contain a \"genotope\" field");
            }
            auto& genotype_json = initial_cell_json["genotype"];

            if (!genotype_json.contains("name")) {
                throw std::domain_error("Every initial cell \"genotoype\" field must contain "
                                         "a \"name\" field");
            }
    
            std::string genotype_name = genotype_json["name"].template get<std::string>();

            if (!genotype_json.contains("epigenetic status")) {
                throw std::domain_error("Every initial cell \"genotoype\" field must contain "
                                         "a \"epigenetic status\" field");
            }

            std::string epigenetic_status = genotype_json["epigenetic status"].template get<std::string>();

            auto& epigenotype = genotypes.at(genotype_name)[epigenetic_status];

            if (!initial_cell_json.contains("genotype")) {
                throw std::domain_error("Every initial cell must contain a \"position\" field");
            }
            simulation.add_cell(epigenotype, get_position(initial_cell_json["position"]));
        }
    }

    static void configure_driver_mutations(const nlohmann::json& driver_mutations_json,
                                           const std::map<std::string, Races::Drivers::Genotype> genotypes)
    {
        using namespace Races;
        using namespace Races::Drivers;

        if (!driver_mutations_json.is_array()) {
            throw std::domain_error("The \"driver mutations\" field must contain an array");
        }

        for (const auto& mutation_json : driver_mutations_json) {
            if (!mutation_json.is_object()) {
                throw std::domain_error("Every element in the \"driver mutations\" field must be an object");
            }

            if (!mutation_json.contains("time")) {
                throw std::domain_error("Every driver mutation must contain a \"time\" field");
            }
            Time time = mutation_json["time"].template get<Time>();

            if (time <= 0) {
                throw std::domain_error("the mutation timing must be positive");
            }

            if (!mutation_json.contains("original genotype")) {
                throw std::domain_error("Every driver mutation must contain a \"original genotype\" field");
            }
            const auto& orig_genotype = genotypes.at(mutation_json["original genotype"].template get<std::string>());

            if (!mutation_json.contains("mutated genotype")) {
                throw std::domain_error("Every driver mutation must contain a \"mutated genotype\" field");
            }
            const auto& mutated_genotype = genotypes.at(mutation_json["mutated genotype"].template get<std::string>());

            simulation.add_genomic_mutation(orig_genotype,mutated_genotype,time);
        }
    }

    static void configure_simulation(const std::string& simulation_filename)
    {
        using namespace Races::Drivers;

        std::ifstream f(simulation_filename);
        nlohmann::json simulation_cfg = nlohmann::json::parse(f);

        std::map<std::string, Genotype> genotypes;

        if (!simulation_cfg.contains("tissue")) {
            throw std::domain_error("The simulation configuration file must contain a \"tissue\" field");
        }
        configure_tissue(simulation_cfg["tissue"]);

        if (!simulation_cfg.contains("genotypes")) {
            throw std::domain_error("The simulation configuration file must contain a \"genotypes\" field");
        }
        if (!simulation_cfg["genotypes"].is_array()) {
            throw std::domain_error("The \"genotypes\" field must be an array of genotype specifications");
        }
        for (const auto& genotype_json: simulation_cfg["genotypes"]) {
            auto genotype = create_genotype(genotype_json);
            simulation.add_species(genotype);

            std::string genotype_name = genotype.get_name();
            genotypes.emplace(genotype_name, std::move(genotype));
        }

        if (!simulation_cfg.contains("initial cells")) {
            throw std::domain_error("The simulation configuration file must contain a \"initial cells\" field");
        }
        configure_initial_cells(simulation_cfg["initial cells"], genotypes);

        if (!simulation_cfg.contains("driver mutations")) {
            throw std::domain_error("The simulation configuration file must contain a \"driver mutations\" field");
        }
        configure_driver_mutations(simulation_cfg["driver mutations"], genotypes);
        
        if (simulation_cfg.contains("death activation level")) {
            simulation.death_activation_level = simulation_cfg["death activation level"].template get<size_t>();
        }
    
        using namespace std::chrono_literals;

        simulation.set_interval_between_snapshots(5min);
    }

    void init_simulation()
    {
        using namespace Races;

        configure_simulation(simulation_filename);

        simulation.duplicate_internal_cells = duplicate_internal_cells;
    }

public:

    DriverSimulator(int argc, char* argv[]):
        BasicExecutable(argv[0], {{"options", "Options"}})
    {
        namespace po = boost::program_options;

        visible_options.at("options").add_options()
            ("recover-simulation,r", 
             po::value<std::string>(&snapshot_path)->default_value(""), 
             "recover a simulation")
            ("border-duplications-only,b", "admit duplications exclusively for cells on borders")
            ("disable-storage,d", "disable result storage")
            ("seed,s", po::value<int>()->default_value(0), "random generator seed")
    #ifdef WITH_INDICATORS
            ("quiet,q", "disable progress bar")
    #endif // WITH_INDICATORS
    #ifdef WITH_SDL2
            ("plot,p", 
            "plot a graphical representation of the simulation")
            ("frames-per-second,f", po::value<unsigned int>(&frames_per_second)->default_value(1), 
            "the number of frames per second")
    #endif
            ("help,h", "get the help");

        hidden_options.add_options()
            ("simulation file", po::value<std::filesystem::path>(&simulation_filename), 
            "the name of the file describing the simulation")
            ("time horizon", po::value<long double>(&time_horizon), 
            "the simulation time horizon");

        po::options_description program_options;
        for (const auto& [section_name, section]: visible_options) {
            program_options.add(section);
        }
        program_options.add(hidden_options);

        positional_options.add("simulation file", 1);
        positional_options.add("time horizon", 1);

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

        if (!vm.count("simulation file")) {
            print_help_and_exit("The simulation file is mandatory", 1);
        }

        if (!vm.count("time horizon")) {
            print_help_and_exit("The simulation time horizon is mandatory", 1);
        }

        quiet = vm.count("quiet")>0;
        plot = vm.count("plot")>0;
        disable_storage = vm.count("disable-storage")>0;
        duplicate_internal_cells = (vm.count("border-duplications-only")==0);
    }

    void run()
    {
        using namespace Races;

        if (!quiet) {
            UI::ProgressBar::hide_console_cursor();

            bar = new UI::ProgressBar();
            bar->set_message("Creating tissue");
        }

        if (snapshot_path != "") {
            snapshot_path = get_last_snapshot_path(snapshot_path, "simulation");

            simulation = load_drivers_simulation(snapshot_path, quiet);
        } else {
            init_simulation();
        }

#ifdef WITH_SDL2
        if (plot) {

            UI::TissuePlotter<UI::SDLWindow> plotter(simulation.tissue());

            plotter.set_frames_per_second(frames_per_second);

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
                simulation.run_up_to(time_horizon, *bar, disable_storage);
            } else {
                simulation.run_up_to(time_horizon, disable_storage);
            }

#ifdef WITH_SDL2
        }
#endif // WITH_SDL2

        UI::ProgressBar::show_console_cursor();

        if (bar != nullptr) {
            delete bar;
        }
    }

};

int main(int argc, char* argv[])
{
    std::signal(SIGINT, termination_handling);

    DriverSimulator simulator(argc, argv);

    simulator.run();

    return 0;
}
