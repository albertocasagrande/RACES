/**
 * @file passengers_sim.cpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Main file for the passenger mutations simulator
 * @version 0.2
 * @date 2023-08-14
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
#include <random>
#include <regex>

#include <boost/program_options.hpp>
#include <nlohmann/json.hpp>

#include "simulation.hpp"

#include "sampler.hpp"
#include "phylogenetic_forest.hpp"

#include "context_index.hpp"

#include "mutation_engine.hpp"

#include "fasta_utils.hpp"
#include "fasta_reader.hpp"

#include "sam_generator.hpp"

#include "progress_bar.hpp"

template<>
std::string boost::lexical_cast<std::string, Races::Passengers::IO::SAMGenerator<>::Mode>(const Races::Passengers::IO::SAMGenerator<>::Mode& mode)
{
    using namespace Races::Passengers::IO;
    switch (mode) {
        case SAMGenerator<>::Mode::CREATE:
            return "create";
        case SAMGenerator<>::Mode::OVERWRITE:
            return "overwrite";
        case SAMGenerator<>::Mode::APPEND:
            return "append";
        default:
            throw std::runtime_error("Unknown mode");
    }
}

template<>
Races::Passengers::IO::SAMGenerator<>::Mode boost::lexical_cast<Races::Passengers::IO::SAMGenerator<>::Mode, std::string>(const std::string& mode)
{
    std::string upper_token = mode;
    for (auto& t_char: upper_token) {
        t_char = toupper(t_char);
    }

    using namespace Races::Passengers::IO;

    if (upper_token == "CREATE") {
        return SAMGenerator<>::Mode::CREATE;
    }

    if (upper_token == "OVERWRITE") {
        return SAMGenerator<>::Mode::OVERWRITE;
    }

    if (upper_token == "APPEND") {
        return SAMGenerator<>::Mode::APPEND;
    }

    std::ostringstream oss;

    oss << "Invalid SAMGenerator mode \"" << mode << "\". "
        << "Supported \"create\", \"overwrite\", and \"append\".";

    namespace po = boost::program_options;

    throw po::validation_error(po::validation_error::kind_t::invalid_option_value, oss.str());
}

class PassengersSimulator
{
    const std::string program_name;
    boost::program_options::options_description visible_options;
    boost::program_options::positional_options_description positional_options;

    std::string simulation_filename;
    std::string drivers_directory;
    std::string context_index_filename;
    std::string ref_genome_filename;
    std::string SBS_filename;

    double coverage;
    std::string SAM_directory;
    Races::Passengers::IO::SAMGenerator<>::Mode SAM_mode;
    size_t read_size;
    std::string SVNs_csv_filename;
    std::string CNAs_csv_filename;

    int seed;
    size_t bytes_per_abs_position;
    bool quiet;

    static std::filesystem::path find_last_snapshot(const std::filesystem::path directory)
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

    Races::Drivers::Simulation::Simulation load_drivers_simulation(const std::string& driver_directory) const
    {
        auto last_snapshot_path = find_last_snapshot(driver_directory);

        Races::Archive::Binary::In archive(last_snapshot_path);

        Races::Drivers::Simulation::Simulation simulation;

        if (quiet) {
            archive & simulation;
        } else {
            archive.load(simulation, "simulation");
        }

        return simulation;
    }

    Races::Drivers::PhylogeneticForest 
    build_forest(const Races::Drivers::Simulation::Simulation& simulation,
                 const nlohmann::json& simulation_cfg) const
    {
        using namespace Races::Drivers;

        if (!simulation_cfg.contains("sample region")) {
            throw std::runtime_error("The passengers simulation file must contain "
                                        "a \"sample region\" field");
        }

        const auto sampler_region = get_sample_region(simulation_cfg["sample region"]);

        Simulation::BinaryLogger::CellReader reader(drivers_directory);

        RectangleSampler sampler(simulation.tissue(), sampler_region.first, sampler_region.second);

        return grow_forest_from(sampler, reader);
    }

    template<typename ABSOLUTE_GENOTYPE_POSITION = uint32_t>
    Races::Passengers::ContextIndex<ABSOLUTE_GENOTYPE_POSITION> load_context_index(const std::string& filename) const
    {
        Races::Archive::Binary::In archive(filename);

        Races::Passengers::ContextIndex<ABSOLUTE_GENOTYPE_POSITION> context_index;

        if (quiet) {
            archive & context_index;
        } else {
            archive.load(context_index, "context index");
        }

        return context_index;
    }

    static std::map<std::string, Races::Passengers::MutationalSignature> load_signatures(const std::string filename)
    {
        std::ifstream is(filename);

        return Races::Passengers::MutationalSignature::read_from_stream(is);
    }

    Races::Passengers::MutationStatistics 
    collect_statistics(const std::list<Races::Passengers::GenomeMutations>& cells_mutations) const
    {
        using namespace Races;
        using namespace Races::Passengers;

        MutationStatistics statistics;

        if (quiet) {
            for (const auto& mutations: cells_mutations) {
                statistics.record(mutations);
            }

            return statistics;
        }

        UI::ProgressBar progress_bar;
        progress_bar.set_message("Collecting statistics");

        size_t recorded{0};
        for (const auto& mutations: cells_mutations) {
            statistics.record(mutations);

            size_t percentage = 100*(++recorded)/cells_mutations.size();
            if (percentage>progress_bar.get_progress()) {
                progress_bar.set_progress(percentage);
            }
        }

        progress_bar.set_message("Statistics collected");

        return statistics;
    }

    template<typename ABSOLUTE_GENOME_POSITION>
    std::list<Races::Passengers::GenomeMutations>
    place_mutations(Races::Passengers::MutationEngine<ABSOLUTE_GENOME_POSITION>& engine,
                    const Races::Drivers::PhylogeneticForest& forest) const
    {
        if (quiet) {
            return engine.place_mutations(forest);
        }

        Races::UI::ProgressBar progress_bar;

        progress_bar.set_message("Placing mutations");

        auto leaves_mutations = engine.place_mutations(forest, progress_bar);

        progress_bar.set_message("Mutations placed");

        return leaves_mutations;
    }

    static std::vector<Races::Drivers::Simulation::AxisValue>
    get_corner(const nlohmann::json& sampler_region_json, const std::string& field_name)
    {
        if (!sampler_region_json.contains(field_name)) {
            throw std::runtime_error("The \"sampler region\" must contain a "
                                     "\"" + field_name + "\" field");
        }

        auto corner_json = sampler_region_json[field_name];

        if (!corner_json.is_array()) {
            throw std::runtime_error("The \"" + field_name + "\" must be an array of natural values");
        }

        std::vector<Races::Drivers::Simulation::AxisValue> corner;

        for (const auto& axis_value : corner_json) {
            corner.push_back(axis_value.template get<Races::Drivers::Simulation::AxisValue>());
        }

        if (corner.size()<2 || corner.size()>3) {
            throw std::runtime_error("\""+field_name+"\" must be either a 2D or 3D array");
        }

        return corner;
    }

    static std::pair<Races::Drivers::Simulation::PositionInTissue,
                     Races::Drivers::Simulation::PositionInTissue> 
    get_sample_region(const nlohmann::json& sampler_region_json)
    {
        if (!sampler_region_json.is_object()) {
            throw std::runtime_error("The \"sample region\" field must be objects");
        }

        auto lower_vector = get_corner(sampler_region_json, "lower corner");
        auto upper_vector = get_corner(sampler_region_json, "upper corner");

        if (lower_vector.size() != upper_vector.size()) {
            throw std::runtime_error("The \"lower corner\" and \"upper corner\" arrays must "
                                     "have the same size");
        }

        using namespace Races::Drivers::Simulation;

        PositionInTissue lower_corner, upper_corner;
        if (lower_vector.size()==2) {
            return {{lower_vector[0], lower_vector[1]},
                    {upper_vector[0], upper_vector[1]}};
        }

        return {{lower_vector[0], lower_vector[1], lower_vector[2]},
                {upper_vector[0], upper_vector[1], upper_vector[2]}};
    }

    static std::pair<std::string, double>
    get_mutation_rate(const nlohmann::json& mutation_rate_json)
    {  
        if (!mutation_rate_json.is_object()) {
            throw std::runtime_error("All the elements in \"mutation rates\" "
                                     "must be objects");
        }
    
        if (!mutation_rate_json.contains("epigenetic status")) {
            throw std::runtime_error("All the elements in \"mutation rates\" "
                                     "must contain a \"epigenetic status\" field");
        }

        if (!mutation_rate_json.contains("rate")) {
            throw std::runtime_error("All the elements in \"mutation rates\" "
                                    "must contain a \"rate\" field");
        }

        return {
            mutation_rate_json["epigenetic status"].template get<std::string>(),
            mutation_rate_json["rate"].template get<double>()
        };
    }

    static std::map<std::string, double>
    get_mutation_rates(const nlohmann::json& mutation_rates_json)
    {
        std::map<std::string, double> mutational_rates;

        if (!mutation_rates_json.is_array()) {
            throw std::runtime_error("The \"mutation rates\" field must be an array "
                                    "of driver mutation rates");
        }

        for (const auto& mutation_rate_json : mutation_rates_json) {
            mutational_rates.insert(get_mutation_rate(mutation_rate_json));
        }

        return mutational_rates;
    }
    
    static Races::Passengers::CopyNumberAlteration::Type
    get_CNA_type(const nlohmann::json& CNA_json)
    {
        if (!CNA_json.contains("subtype")) {
            throw std::runtime_error("All the CNA mutations must contain a \"subtype\" field");
        }

        auto subtype_str = CNA_json["subtype"].template get<std::string>();
        for (auto& subtype_chr:subtype_str) {
            subtype_chr = tolower(subtype_chr);
        }

        if (subtype_str == "amplification") {
            return Races::Passengers::CopyNumberAlteration::Type::AMPLIFICATION;
        }

        if (subtype_str == "deletion") {
            return Races::Passengers::CopyNumberAlteration::Type::DELETION;
        }

        throw std::runtime_error("Unknown CNA \""+
                                 CNA_json["subtype"].template get<std::string>()+"\"");
    }

    static void
    add_CNA(std::list<Races::Passengers::CopyNumberAlteration>& CNAs, 
            const nlohmann::json& CNA_json)
    {
        using namespace Races::Passengers;

        if (!CNA_json.contains("chromosome")) {
            throw std::runtime_error("All the CNA mutations "
                                    "must contain a \"chromosome\" field");
        }

        if (!CNA_json.contains("position")) {
            throw std::runtime_error("All the CNA mutations "
                                    "must contain a \"position\" field");
        }

        if (!CNA_json.contains("length")) {
            throw std::runtime_error("All the CNA mutations "
                                     "must contain a \"length\" field");
        }

        if (!CNA_json.contains("allele")) {
            throw std::runtime_error("All the CNA mutations "
                                     "must contain an \"allele\" field");
        }

        auto chr_str = CNA_json["chromosome"].template  get<std::string>();
        GenomicPosition genomic_position(GenomicPosition::stochr(chr_str),
                                         CNA_json["position"].template get<ChrPosition>());

        auto length = CNA_json["length"].template  get<ChrPosition>();
        auto allele = CNA_json["allele"].template  get<uint8_t>();

        switch (get_CNA_type(CNA_json)) {
            case CopyNumberAlteration::Type::AMPLIFICATION:
                CNAs.push_back(CopyNumberAlteration::new_amplification({genomic_position, length}, allele, 0));
                break;
            case CopyNumberAlteration::Type::DELETION:
                CNAs.push_back(CopyNumberAlteration::new_deletion({genomic_position, length}, allele));
                break;
            default:
                throw std::runtime_error("Unsupported CNA type");
        }
    }

    static void
    add_SNV(const std::string& name, std::list<Races::Passengers::SNV>& SNVs,
            const nlohmann::json& SNV_json)
    {
        using namespace Races::Passengers;

        if (!SNV_json.contains("chromosome")) {
            throw std::runtime_error("All the SNV mutations "
                                    "must contain a \"chromosome\" field");
        }

        if (!SNV_json.contains("position")) {
            throw std::runtime_error("All the SNV mutations "
                                    "must contain a \"position\" field");
        }

        if (!SNV_json.contains("context")) {
            throw std::runtime_error("All the SNV mutations "
                                    "must contain a \"context\" field");
        }

        if (!SNV_json.contains("mutated base")) {
            throw std::runtime_error("All the SNV mutations "
                                    "must contain a \"mutated base\" field");
        }

        auto chr_str = SNV_json["chromosome"].template  get<std::string>();

        GenomicPosition genomic_position(GenomicPosition::stochr(chr_str),
                                         SNV_json["position"].template get<ChrPosition>());

        MutationalContext context = SNV_json["context"].template  get<std::string>();
        std::string mutated_base = SNV_json["mutated base"].template  get<std::string>();

        SNVs.emplace_back(genomic_position, context, mutated_base[0], name);
    }

    static void
    add_driver_mutation(const std::string& name,
                        std::list<Races::Passengers::SNV>& SNVs,
                        std::list<Races::Passengers::CopyNumberAlteration>& CNAs,
                        const nlohmann::json& driver_mutation_json)
    {
        if (!driver_mutation_json.is_object()) {
            throw std::runtime_error("All the elements in \"mutations\" must be objects");
        }

        if (!driver_mutation_json.contains("type")) {
            throw std::runtime_error("All the elements in \"mutations\" "
                                     "must contain a \"type\" field");
        }

        auto type = driver_mutation_json["type"].template get<std::string>();
        if (type=="SNV") {
            add_SNV(name, SNVs, driver_mutation_json);

            return;
        }

        if (type=="CNA") {
            add_CNA(CNAs, driver_mutation_json);

            return;
        }

        throw std::runtime_error("Unsupported mutations type \""+type+"\"");
    }

    static void collect_driver_mutations(const std::string& name,
                                         std::list<Races::Passengers::SNV>& SNVs,
                                         std::list<Races::Passengers::CopyNumberAlteration>& CNAs,
                                         const nlohmann::json& driver_mutations_json)
    {
        if (!driver_mutations_json.is_array()) {
            throw std::runtime_error("The \"mutations\" field must be an array "
                                     "of mutations");
        }

        for (const auto& driver_mutation_json : driver_mutations_json) {
            add_driver_mutation(name, SNVs, CNAs, driver_mutation_json);
        }
    }

    static void add_driver_mutational_properties(Races::Passengers::SpeciesMutationalProperties& mutational_properties,
                                                 const Races::Drivers::Simulation::Simulation& drivers_simulation,
                                                 const nlohmann::json& driver_properties_json)
    {
        if (!driver_properties_json.is_object()) {
            throw std::runtime_error("All the elements in \"driver properties\" "
                                     "must be objects");
        }

        if (!driver_properties_json.contains("name")) {
            throw std::runtime_error("All the elements in \"driver properties\" "
                                     "must contain a \"name\" field");
        }

        auto name = driver_properties_json["name"].template get<std::string>();

        if (!driver_properties_json.contains("mutation rates")) {
            throw std::runtime_error("All the elements in \"driver properties\" "
                                    "must contain a \"mutation rates\" field");
        }
        auto mutation_rates = get_mutation_rates(driver_properties_json["mutation rates"]);

        using namespace Races::Passengers;
        std::list<SNV> SNVs;
        std::list<CopyNumberAlteration> CNAs;
        if (driver_properties_json.contains("mutations")) {
            collect_driver_mutations(name, SNVs, CNAs, driver_properties_json["mutations"]);
        }
        
        mutational_properties.add_species(drivers_simulation, name, mutation_rates, SNVs);
    }

    static Races::Passengers::SpeciesMutationalProperties
    get_mutational_properties(const Races::Drivers::Simulation::Simulation& drivers_simulation,
                              const nlohmann::json& simulation_cfg)
    {
        using namespace Races::Passengers;

        SpeciesMutationalProperties mutational_properties;

        if (!simulation_cfg.contains("driver properties")) {
            throw std::runtime_error("The passengers simulation configuration must contain "
                                     "a \"driver properties\" field");
        }

        auto& mutational_properties_json = simulation_cfg["driver properties"];

        if (!mutational_properties_json.is_array()) {
            throw std::runtime_error("The \"driver properties\" field must be an array "
                                     "of driver mutational properties");
        }

        for (const auto& driver_properties_json : mutational_properties_json) {
            add_driver_mutational_properties(mutational_properties, drivers_simulation,
                                             driver_properties_json);
        }

        return mutational_properties;
    }

    static double
    get_fraction(const nlohmann::json& fraction_json)
    {
        const double fraction = fraction_json.template get<double>();

        if (fraction<0) {
            throw std::runtime_error("The \"fraction\" field must be a non-negative value");
        }

        return fraction;
    }

    static uint8_t
    get_allele(const nlohmann::json& json)
    {
        if (!json.contains("allele")) {
            throw std::runtime_error("Expected a missing \"allele\" field");
        }

        return json.template get<uint8_t>();
    } 

    static double
    get_number_of_alleles(const nlohmann::json& simulation_cfg)
    {
        if (!simulation_cfg.contains("number of alleles")) {
            throw std::runtime_error("The passengers simulation configuration must contain "
                                     "a \"number of alleles\" field");
        }

        double value = simulation_cfg["number of alleles"].template get<double>();

        if (value <= 0) {
            throw std::runtime_error("The number of allele must be positive");
        }

        return value;
    }

    static std::map<std::string, double>
    get_mutational_coefficients(const nlohmann::json& mutational_coefficients_json)
    {
        if (!mutational_coefficients_json.is_array()) {
            throw std::runtime_error("Mutational coefficients must be an array of object");
        }

        double total=0;

        std::map<std::string, double> mutational_coefficients;
        for (const auto& mutational_coefficient_json : mutational_coefficients_json) {

            if (!mutational_coefficient_json.is_object()) {
                throw std::runtime_error("Mutational coefficient must be an object");
            }

            if (!mutational_coefficient_json.contains("SBS")) {
                throw std::runtime_error("Mutational coefficient must contain a \"SBS\" field");
            }

            const std::string SBS = mutational_coefficient_json["SBS"].template get<std::string>();
            if (mutational_coefficients.count(SBS)>0) {
                throw std::runtime_error("\""+SBS+"\" already among mutational coefficients");
            }

            if (!mutational_coefficient_json.contains("fraction")) {
                throw std::runtime_error("Mutational coefficient must contain a \"fraction\" field");
            }

            double fraction = get_fraction(mutational_coefficient_json["fraction"]);

            total += fraction;
            if (total>1) {
                std::ostringstream oss;

                oss << "The sum of the \"fraction\" fields must be 1. "
                    << "The difference is " << (1-total);

                throw std::runtime_error(oss.str());
            }

            mutational_coefficients[SBS] = fraction;
        }
        
        if (abs(1-total)>10*std::numeric_limits<double>::epsilon()) {
            std::ostringstream oss;

            oss << "The sum of the \"fraction\" fields must be 1. "
                << "The difference is " << (1-total);

            throw std::runtime_error(oss.str());
        }

        return mutational_coefficients;
    }

    static std::map<std::string, double>
    get_default_mutational_coefficients(const nlohmann::json& mutational_coeff_json)
    {
        if (!mutational_coeff_json.is_object()) {
            throw std::runtime_error("The \"mutational coefficients\" field must be an object");
        }

        if (!mutational_coeff_json.contains("default")) {
            throw std::runtime_error("The \"mutational coefficients\" field must contain "
                                     "a \"default\" field");
        }

        return get_mutational_coefficients(mutational_coeff_json["default"]);
    }

    template<typename ABSOLUTE_POSITION_TYPE>
    static void
    add_timed_mutational_coefficients(Races::Passengers::MutationEngine<ABSOLUTE_POSITION_TYPE>& engine, 
                                      const nlohmann::json& mutational_coeff_json)
    {
        std::map<double, std::map<std::string, double>> timed_coefficients;

        if (!mutational_coeff_json.is_object()) {
            throw std::runtime_error("The \"mutational coefficients\" field must be an object");
        }

        if (mutational_coeff_json.contains("timed")) {

            auto& timed_coeff_json = mutational_coeff_json["timed"];

            if (!timed_coeff_json.is_array()) {
                throw std::runtime_error("The optional \"timed\" field must be an array");
            }

            for (const auto& timed_json : timed_coeff_json) {
                if (!timed_json.is_object()) {
                    throw std::runtime_error("The elements of the \"timed\" field must be objects");
                }

                if (!timed_json.contains("time")) {
                    throw std::runtime_error("The elements of the \"timed\" field must contain "
                                             "a \"time\" field");
                }

                double time = timed_json["time"].template get<double>();

                if (timed_coefficients.count(time)>0) {
                    throw std::runtime_error("Two elements of the \"timed\" field have the "
                                             "same \"time\"");
                }

                if (!timed_json.contains("coefficients")) {
                    throw std::runtime_error("The elements of the \"timed\" field must contain "
                                             "a \"coefficients\" field");
                }

                engine.add(time, get_mutational_coefficients(mutational_coeff_json["default"]));
            }
        }
    }

    nlohmann::json get_simulation_json() const
    {
        std::ifstream simulation_stream(simulation_filename);

        return nlohmann::json::parse(simulation_stream);
    }

    void process_statistics(const std::list<Races::Passengers::GenomeMutations>& mutations) const
    {
        if ((SVNs_csv_filename != "") || (CNAs_csv_filename != "")) {
            auto statistics = collect_statistics(mutations);

            if (SVNs_csv_filename != "") {
                std::ofstream os(SVNs_csv_filename);

                statistics.write_SNVs_table(os);
            }

            if (CNAs_csv_filename != "") {
                std::ofstream os(CNAs_csv_filename);

                statistics.write_CNAs_table(os);
            }
        }
    }

    template<typename ABSOLUTE_POSITION_TYPE>
    void run_abs_position() const
    {
        using namespace Races;
        using namespace Races::Passengers;
        using namespace Races::Passengers::IO;

        nlohmann::json simulation_cfg = get_simulation_json();

        SAMGenerator<> sam_generator;

        if (coverage>0) {
            sam_generator = SAMGenerator<>(SAM_directory, ref_genome_filename, read_size, seed, 
                                           SAM_mode);
        }

        Drivers::PhylogeneticForest forest;
        SpeciesMutationalProperties mutational_properties;

        {
            auto drivers_simulation = load_drivers_simulation(drivers_directory);

            forest = build_forest(drivers_simulation, simulation_cfg);

            mutational_properties = get_mutational_properties(drivers_simulation, simulation_cfg);
        }

        auto signatures = load_signatures(SBS_filename);

        auto context_index = load_context_index<ABSOLUTE_POSITION_TYPE>(context_index_filename);

        if (!simulation_cfg.contains("mutational coefficients")) {
            throw std::runtime_error("The passengers simulation configuration must contain "
                                     "a \"mutational coefficients\" field");
        }

        const auto& mutational_coeff_json = simulation_cfg["mutational coefficients"];

        auto default_coefficients = get_default_mutational_coefficients(mutational_coeff_json);

        auto num_of_allele = get_number_of_alleles(simulation_cfg);

        MutationEngine<ABSOLUTE_POSITION_TYPE> engine(context_index, num_of_allele,
                                                      default_coefficients, signatures,
                                                      mutational_properties);

        add_timed_mutational_coefficients(engine, mutational_coeff_json);

        auto leaves_mutations = place_mutations(engine, forest);

        if (coverage) {
            sam_generator(leaves_mutations, coverage);
        }

        process_statistics(leaves_mutations);
    }

    void validate(boost::program_options::variables_map vm) const
    {
        namespace fs = std::filesystem;

        if (vm.count("help")) {
            print_help_and_exit(0);
        }

        if (!vm.count("simulation file")) {
            print_help_and_exit("The simulation file is mandatory", 1);
        }
    
        if (!vm.count("driver directory")) {
            print_help_and_exit("The driver directory is mandatory. "
                                "You can produce it by using `races_sim`", 1);
        }

        if (!fs::exists(drivers_directory)) {
            print_help_and_exit("\"" + drivers_directory + "\"  does not exist", 1);
        }

        if (!fs::is_directory(drivers_directory)) {
            std::cerr << "\"" << drivers_directory << "\"  is not a directory"
                      << std::endl << std::endl;
            print_help_and_exit(1);
        }

        try {
            find_last_snapshot(drivers_directory);
        } catch (std::runtime_error& ex) {
            std::cerr << ex.what() << std::endl << std::endl;
            print_help_and_exit(1);
        }

        if (!vm.count("context index")) {
            std::cerr << "The context index is mandatory. "
                      << "You can produce it by using `races_build_index`" << std::endl 
                      << std::endl;
            print_help_and_exit(1);
        }

        if (!fs::exists(context_index_filename)) {
            std::cerr << "\"" << context_index_filename << "\"  does not exist"
                      << std::endl << std::endl;
            print_help_and_exit(1);
        }

        if (!fs::is_regular_file(context_index_filename)) {
            std::cerr << "\"" << context_index_filename << "\"  is not a regular file"
                      << std::endl << std::endl;
            print_help_and_exit(1);
        }

        if (!vm.count("mutational signature")) {
            std::cerr << "The mutational signature is mandatory." << std::endl << std::endl;
            print_help_and_exit(1);
        }

        if (!fs::exists(SBS_filename)) {
            std::cerr << "\"" << SBS_filename << "\"  does not exist"
                      << std::endl << std::endl;
            print_help_and_exit(1);
        }

        if (!fs::is_regular_file(SBS_filename)) {
            std::cerr << "\"" << SBS_filename << "\"  is not a regular file"
                      << std::endl << std::endl;
            print_help_and_exit(1);
        }

        if (!vm.count("reference genome")) {
            std::cerr << "The reference genome FASTA file is mandatory." << std::endl << std::endl;
            print_help_and_exit(1);
        }

        if (!fs::exists(ref_genome_filename)) {
            std::cerr << "\"" << ref_genome_filename << "\"  does not exist"
                      << std::endl << std::endl;
            print_help_and_exit(1);
        }

        if (!fs::is_regular_file(ref_genome_filename)) {
            std::cerr << "\"" << ref_genome_filename << "\"  is not a regular file"
                      << std::endl << std::endl;
            print_help_and_exit(1);
        }

        if (vm.count("SNVs-csv")>0 && fs::exists(SVNs_csv_filename)) {
            std::cerr << "\"" << SVNs_csv_filename << "\" already exists" << std::endl << std::endl;
            print_help_and_exit(1);
        }

        if (vm.count("CNAs-csv")>0 && fs::exists(CNAs_csv_filename)) {
            std::cerr << "\"" << CNAs_csv_filename << "\" already exists" << std::endl << std::endl;
            print_help_and_exit(1);
        }

        if (coverage<0) {
            std::cerr << "coverage must be non-negative." << std::endl << std::endl;
            print_help_and_exit(1);
        }
        
        for (const auto SAM_option : {"SAM-directory"}) {
            if ((coverage>0) && (vm.count(SAM_option)==0)) {
                std::cerr << "\"--"<< SAM_option << "\" is mandatory when coverage differs from 0." 
                          << std::endl << std::endl;
                print_help_and_exit(1);
            }
        }
    }

public:
    void print_help_and_exit(const std::string& message, const int& err_code) const
    { 
        std::ostream& os = (err_code==0 ? std::cout: std::cerr);
        
        os << message << std::endl << std::endl;
        
        print_help_and_exit(err_code);
    }

    void print_help_and_exit(const int& err_code) const
    { 
        std::ostream& os = (err_code==0 ? std::cout: std::cerr);
        
        os << "Syntax: " << program_name;
        
        for (unsigned int i=0;i<positional_options.max_total_count(); ++i) {
            os << " <" << positional_options.name_for_position(i) << ">"; 
        }

        os << std::endl << visible_options << std::endl;

        exit(err_code);
    }

    PassengersSimulator(int argc, char* argv[]):
        program_name(argv[0]), visible_options("Options"),
        SVNs_csv_filename(""), CNAs_csv_filename("")
    {
        namespace po = boost::program_options;

        using namespace Races::Passengers::IO;

        visible_options.add_options()
            ("coverage,c", po::value<double>(&coverage)->default_value(0.0), 
             "coverage for the simulated reads (greater than 0 to produce SAM files)")
            ("SAM-directory,d", po::value<std::string>(&SAM_directory), 
             "the output directory for the produces SAM files")
            ("SAM-mode,m", po::value<SAMGenerator<>::Mode>(&SAM_mode)->default_value(SAMGenerator<>::Mode::CREATE), 
             "SAM generator mode (create, overwrite, append)")
            ("read-size,r", po::value<size_t>(&read_size)->default_value(150),
             "simulated read size")
            ("SNVs-CSV,S", po::value<std::string>(&SVNs_csv_filename),
             "the SNVs CSV output file")
            ("CNAs-CSV,C", po::value<std::string>(&CNAs_csv_filename),
             "the CNAs CSV output file")
            ("seed,s", po::value<int>(&seed)->default_value(0), 
             "random generator seed")
            ("quiet,q", "avoid output messages")
            ("help,h", "get the help")
        ;

        po::options_description hidden("Hidden options");
        hidden.add_options()
            ("simulation file", po::value<std::string>(&simulation_filename), 
             "the name of the file describing the passenger mutations simulation")
            ("driver directory", po::value<std::string>(&drivers_directory), 
             "the driver directory")
            ("context index", po::value<std::string>(&context_index_filename), 
             "the genome context index")
            ("mutational signature", po::value<std::string>(&SBS_filename), 
             "the mutational signature")
            ("reference genome", po::value<std::string>(&ref_genome_filename), 
             "the filename of the reference genome FASTA file")
        ;

        po::options_description generic;
        generic.add(visible_options).add(hidden);

        positional_options.add("simulation file", 1);
        positional_options.add("driver directory", 1);
        positional_options.add("context index", 1);
        positional_options.add("mutational signature", 1);
        positional_options.add("reference genome", 1);

        po::variables_map vm;

        try {
            po::command_line_parser parser{argc, argv};
            po::store(parser.options(generic).positional(positional_options).run(), vm);
            po::notify(vm);
        } catch (boost::wrapexcept<boost::program_options::unknown_option> &except) {
            print_help_and_exit(except.what(), 1);
        }

        quiet = vm.count("quiet")>0;

        try {
            using namespace Races::Passengers;

            Races::Archive::Binary::In archive(context_index_filename);

            bytes_per_abs_position = ContextIndex<>::read_bytes_per_absolute_position(archive);

        } catch (std::exception& except) {
            print_help_and_exit(except.what(), 1);
        }
    
        validate(vm);
    }

    void run() const
    {
        try {
            switch (bytes_per_abs_position) {
                case 2:
                    run_abs_position<uint16_t>();
                    break;
                case 4:
                    run_abs_position<uint32_t>();
                    break;
                case 8:
                    run_abs_position<uint64_t>();
                    break;
                default:
                    std::cerr << "Unsupported bits per absolute position."
                              << " Supported values are 2, 4, or 8." << std::endl;
                    exit(1);
            }
        } catch (std::exception& except) {
            print_help_and_exit(except.what(), 1);
        }
    }
};

int main(int argc, char* argv[])
{
    PassengersSimulator passenger_sim(argc, argv);

    passenger_sim.run();

    return 0;
}
