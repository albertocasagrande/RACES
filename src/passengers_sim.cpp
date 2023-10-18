/**
 * @file passengers_sim.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Main file for the RACES passenger mutations simulator
 * @version 0.16
 * @date 2023-10-18
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

#include "common.hpp"

#include "simulation.hpp"

#include "phylogenetic_forest.hpp"

#include "context_index.hpp"

#include "mutation_engine.hpp"

#include "fasta_utils.hpp"
#include "fasta_reader.hpp"

#include "read_simulator.hpp"

#include "filter.hpp"

#include "progress_bar.hpp"

#include "json_config.hpp"

template<>
std::string 
boost::lexical_cast<std::string, Races::Passengers::SequencingSimulations::ReadSimulator<>::Mode>(const Races::Passengers::SequencingSimulations::ReadSimulator<>::Mode& mode)
{
    using namespace Races::Passengers::SequencingSimulations;
    switch (mode) {
        case ReadSimulator<>::Mode::CREATE:
            return "create";
        case ReadSimulator<>::Mode::OVERWRITE:
            return "overwrite";
        case ReadSimulator<>::Mode::APPEND:
            return "append";
        default:
            throw std::runtime_error("Unknown mode");
    }
}

template<>
Races::Passengers::SequencingSimulations::ReadSimulator<>::Mode 
boost::lexical_cast<Races::Passengers::SequencingSimulations::ReadSimulator<>::Mode, std::string>(const std::string& mode)
{
    std::string upper_token = mode;
    for (auto& t_char: upper_token) {
        t_char = toupper(t_char);
    }

    using namespace Races::Passengers::SequencingSimulations;

    if (upper_token == "CREATE") {
        return ReadSimulator<>::Mode::CREATE;
    }

    if (upper_token == "OVERWRITE") {
        return ReadSimulator<>::Mode::OVERWRITE;
    }

    if (upper_token == "APPEND") {
        return ReadSimulator<>::Mode::APPEND;
    }

    std::ostringstream oss;

    oss << "Invalid ReadSimulator mode \"" << mode << "\". "
        << "Supported \"create\", \"overwrite\", and \"append\".";

    namespace po = boost::program_options;

    throw po::validation_error(po::validation_error::kind_t::invalid_option_value, oss.str());
}

class PassengersSimulator : public BasicExecutable
{
    std::filesystem::path simulation_filename;
    std::filesystem::path drivers_directory;
    std::filesystem::path snapshot_path;
    std::filesystem::path context_index_filename;
    std::filesystem::path ref_genome_filename;
    std::filesystem::path SBS_filename;

    double coverage;
    std::string seq_output_directory;
    Races::Passengers::SequencingSimulations::ReadSimulator<>::Mode read_simulator_output_mode;
    bool paired_read;
    size_t read_size;
    size_t insert_size;
    std::string SNVs_csv_filename;
    std::string CNAs_csv_filename;
    bool epigenetic_FACS;
    bool write_SAM;

    int seed;
    size_t bytes_per_abs_position;
    bool quiet;

    Races::Drivers::PhylogeneticForest 
    build_forest(const Races::Drivers::Simulation::Simulation& simulation) const
    {
        using namespace Races::Drivers;
        using namespace Races::Drivers::Simulation;

        std::list<Races::Drivers::CellId> sampled_ids;
        try {
            sampled_ids = BinaryLogger::load_sampled_ids(drivers_directory);
        } catch (...) {
            print_help_and_exit("\""+std::string(drivers_directory)+"\" does not contain a list "
                                + "of sampled cell. Produce it by using \"tissue_sampler\".", 1);
        }

        BinaryLogger::CellReader reader(drivers_directory);

        auto genotypes = simulation.tissue().get_genotypes();

        return grow_forest_from(sampled_ids, reader, genotypes);
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

    template<typename GENOME_MUTATION, typename FILTER, 
             std::enable_if_t<std::is_base_of_v<Races::Passengers::GenomeMutations, GENOME_MUTATION>
                                && std::is_base_of_v<Races::BaseFilter<Races::Drivers::EpigeneticGenotypeId>, FILTER>, bool> = true>
    Races::Passengers::MutationStatistics 
    collect_statistics(const std::list<GENOME_MUTATION>& cells_mutations, const FILTER& filter,
                       const std::string& data_name="") const
    {
        using namespace Races;
        using namespace Races::Passengers;

        MutationStatistics statistics;

        if (quiet) {
            for (const auto& cell_mutations: cells_mutations) {
                statistics.record(cell_mutations, filter);
            }

            return statistics;
        }

        UI::ProgressBar progress_bar;

        if (data_name.size()>0) {
            progress_bar.set_message("Collecting "+ data_name+ "'s data");
        } else {
            progress_bar.set_message("Collecting data");
        }
        size_t recorded{0};
        for (const auto& cell_mutations: cells_mutations) {
            statistics.record(cell_mutations, filter);

            size_t percentage = 100*(++recorded)/cells_mutations.size();
            if (percentage>progress_bar.get_progress()) {
                progress_bar.set_progress(percentage);
            }
        }

        if (data_name.size()>0) {
            progress_bar.set_message(data_name+ "'s data collected");
        } else {
            progress_bar.set_message("Data collected");
        }

        return statistics;
    }

    template<typename GENOME_MUTATION,
             std::enable_if_t<std::is_base_of_v<Races::Passengers::GenomeMutations, GENOME_MUTATION>, bool> = true>
    Races::Passengers::MutationStatistics 
    collect_statistics(const std::list<GENOME_MUTATION>& cells_mutations) const
    {
        Races::BaseFilter<Races::Drivers::EpigeneticGenotypeId> filter;

        return collect_statistics(cells_mutations, filter);
    }

    template<typename ABSOLUTE_GENOME_POSITION, typename RANDOM_GENERATOR>
    std::list<Races::Passengers::CellGenomeMutations>
    place_mutations(Races::Passengers::MutationEngine<ABSOLUTE_GENOME_POSITION,RANDOM_GENERATOR>& engine,
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

    nlohmann::json get_simulation_json() const
    {
        std::ifstream simulation_stream(simulation_filename);

        return nlohmann::json::parse(simulation_stream);
    }

    template<typename GENOME_MUTATION, std::enable_if_t<std::is_base_of_v<Races::Passengers::GenomeMutations, GENOME_MUTATION>, bool> = true>
    void process_statistics(const std::list<GENOME_MUTATION>& mutations) const
    {
        if ((SNVs_csv_filename != "") || (CNAs_csv_filename != "")) {
            auto statistics = collect_statistics(mutations);

            if (SNVs_csv_filename != "") {
                std::ofstream os(SNVs_csv_filename);

                statistics.write_SNVs_table(os);
            }

            if (CNAs_csv_filename != "") {
                std::ofstream os(CNAs_csv_filename);

                statistics.write_CNAs_table(os);
            }
        }
    }

    template<typename GENOME_MUTATION, std::enable_if_t<std::is_base_of_v<Races::Passengers::GenomeMutations, GENOME_MUTATION>, bool> = true>
    void process_statistics(const std::list<GENOME_MUTATION>& mutations,
                            const std::map<std::string, std::set<Races::Drivers::EpigeneticGenotypeId>>& epigenetic_classes) const
    {
        if ((SNVs_csv_filename != "") || (CNAs_csv_filename != "")) {
            std::string prefix_SNVs = SNVs_csv_filename.substr(0,SNVs_csv_filename.find_last_of('.'));
            std::string prefix_CNAs = CNAs_csv_filename.substr(0,CNAs_csv_filename.find_last_of('.'));

            for (const auto& [signature, id_set] : epigenetic_classes) {
                Races::FilterNotIn<Races::Drivers::EpigeneticGenotypeId> filter(id_set);

                auto statistics = collect_statistics(mutations, filter, signature);

                if (SNVs_csv_filename != "") {
                    std::ofstream os(prefix_SNVs+"_"+signature+".csv");

                    statistics.write_SNVs_table(os);
                }

                if (CNAs_csv_filename != "") {
                    std::ofstream os(prefix_CNAs+"_"+signature+".csv");

                    statistics.write_CNAs_table(os);
                }
            }
        }
    }

    static std::map<std::string, std::set<Races::Drivers::EpigeneticGenotypeId>>
    split_by_epigenetic_status(const std::vector<Races::Drivers::EpigeneticGenotype>& genotypes)
    {
        using namespace  Races::Drivers;
        std::map<std::string, std::set<EpigeneticGenotypeId>> epigenetic_status_classes;

        for (const auto& genotype : genotypes) {
            const auto signature = genotype.get_methylation_signature();

            const auto signature_string = Genotype::signature_to_string(signature);

            epigenetic_status_classes[signature_string].insert(genotype.get_id());
        }

        return epigenetic_status_classes;
    }

    template<typename GENOME_WIDE_POSITION>
    void run_abs_position() const
    {
        using namespace Races;
        using namespace Races::Passengers;
        using namespace Races::Passengers::SequencingSimulations;

        nlohmann::json simulation_cfg = get_simulation_json();

        ReadSimulator<> read_simulator;

        Drivers::PhylogeneticForest forest;
        SpeciesMutationalProperties mutational_properties;

        std::map<std::string, std::set<Drivers::EpigeneticGenotypeId>> epigenetic_status_classes;
        {
            auto drivers_simulation = load_drivers_simulation(snapshot_path, quiet);

            if (epigenetic_FACS) {
                auto tissue_genotypes = drivers_simulation.tissue().get_genotypes();

                epigenetic_status_classes = split_by_epigenetic_status(tissue_genotypes);
            }

            forest = build_forest(drivers_simulation);

            mutational_properties = ConfigReader::get_mutational_properties(drivers_simulation,
                                                                            simulation_cfg);
        }

        if (coverage>0) {
            if (paired_read) {
                read_simulator = ReadSimulator<>(seq_output_directory, ref_genome_filename, read_size, 
                                                 insert_size, read_simulator_output_mode, seed);
            } else {
                read_simulator = ReadSimulator<>(seq_output_directory, ref_genome_filename, read_size, 
                                                 read_simulator_output_mode, seed);
            }
        }

        auto signatures = load_signatures(SBS_filename);

        auto context_index = load_context_index<GENOME_WIDE_POSITION>(context_index_filename);

        if (!simulation_cfg.contains("mutational coefficients")) {
            throw std::runtime_error("The passengers simulation configuration must contain "
                                     "a \"mutational coefficients\" field");
        }

        const auto& mutational_coeff_json = simulation_cfg["mutational coefficients"];

        auto default_coefficients = ConfigReader::get_default_mutational_coefficients(mutational_coeff_json);

        auto num_of_allele = ConfigReader::get_number_of_alleles(simulation_cfg);

        MutationEngine<GENOME_WIDE_POSITION> engine(context_index, num_of_allele,
                                                    default_coefficients, signatures,
                                                    mutational_properties);

        ConfigReader::add_timed_mutational_coefficients(engine, mutational_coeff_json);

        auto leaves_mutations = place_mutations(engine, forest);

        if (coverage>0) {
            read_simulator.enable_SAM_writing(write_SAM);

            read_simulator(leaves_mutations, coverage, epigenetic_status_classes);
        }

        if (epigenetic_FACS) {
            process_statistics(leaves_mutations);
            process_statistics(leaves_mutations, epigenetic_status_classes);
        } else {
            process_statistics(leaves_mutations);
        }
    }

    void validate(boost::program_options::variables_map vm) const
    {
        namespace fs = std::filesystem;

        if (!vm.count("simulation file")) {
            print_help_and_exit("The simulation file is mandatory", 1);
        }
    
        if (!vm.count("driver simulation")) {
            print_help_and_exit("The driver simulation directory is mandatory. "
                                "You can produce it by using `drivers_sim`", 1);
        }
        
        if (!fs::exists(drivers_directory)) {
            print_help_and_exit("\"" + std::string(drivers_directory) + "\"  does not exist", 1);
        }

        if (!vm.count("context index")) {
            print_help_and_exit("The context index is mandatory. "
                                "You can produce it by using `races_build_index`", 1);
        }

        if (!fs::exists(context_index_filename)) {
            print_help_and_exit("\"" + std::string(context_index_filename)
                                + "\"  does not exist", 1);
        }

        if (!fs::is_regular_file(context_index_filename)) {
            print_help_and_exit("\"" + std::string(context_index_filename)
                                + "\"  is not a regular file", 1);
        }

        if (!vm.count("mutational signature")) {
            print_help_and_exit("The mutational signature is mandatory.", 1);
        }

        if (!fs::exists(SBS_filename)) {
            print_help_and_exit("\"" + std::string(SBS_filename)
                                + "\"  does not exist", 1);
        }

        if (!fs::is_regular_file(SBS_filename)) {
            print_help_and_exit("\"" + std::string(SBS_filename)
                                + "\"  is not a regular file", 1);
        }

        if (!vm.count("reference genome")) {
            print_help_and_exit("The reference genome FASTA file is mandatory.", 1);
        }

        if (!fs::exists(ref_genome_filename)) {
            print_help_and_exit("\"" + std::string(ref_genome_filename)
                                + "\"  does not exist", 1);
        }

        if (!fs::is_regular_file(ref_genome_filename)) {
            print_help_and_exit("\"" + std::string(ref_genome_filename)
                                + "\"  is not a regular file", 1);
        }

        if (vm.count("SNVs-csv")>0 && fs::exists(SNVs_csv_filename)) {
            print_help_and_exit("\"" + std::string(SNVs_csv_filename)
                                + "\" already exists", 1);
        }

        if (vm.count("CNAs-csv")>0 && fs::exists(CNAs_csv_filename)) {
            print_help_and_exit("\"" + std::string(CNAs_csv_filename)
                                + "\" already exists", 1);
        }

        if (coverage<0) {
            print_help_and_exit("coverage must be non-negative.", 1);
        }
        
        for (const auto sequencing_option : {"output-directory"}) {
            if ((coverage>0) && (vm.count(sequencing_option)==0)) {
                print_help_and_exit("\"--"+ std::string(sequencing_option) 
                                     + "\" is mandatory when coverage differs from 0.", 1);
            }
        }

        if ((coverage>0) && std::filesystem::exists(seq_output_directory)) {
            if (vm.count("overwrite")==0) {
                print_help_and_exit("The output directory \""+
                                    std::string(seq_output_directory)+"\" already exists", 1);
            }

            if (!fs::is_directory(seq_output_directory)) {
                print_help_and_exit("\""+seq_output_directory+"\" is not a directory", 1);
            }
        }
    }

public:
    PassengersSimulator(int argc, char* argv[]):
        BasicExecutable(argv[0], {{"passengers", "Passenger mutation simulation options"},
                                  {"sequencing", "Sequencing simulation options"},
                                  {"drivers", "Driver simulation related options"},
                                  {"generic", "Generic options"}}),
        paired_read(false), 
        insert_size(0), SNVs_csv_filename(""), CNAs_csv_filename(""), epigenetic_FACS(false)
    {
        namespace po = boost::program_options;

        using namespace Races::Passengers::SequencingSimulations;

        visible_options.at("passengers").add_options()
            ("SNVs-CSV,S", po::value<std::string>(&SNVs_csv_filename),
             "the SNVs CSV output file")
            ("CNAs-CSV,C", po::value<std::string>(&CNAs_csv_filename),
             "the CNAs CSV output file")
        ;

        visible_options.at("sequencing").add_options()
            ("coverage,c", po::value<double>(&coverage)->default_value(0.0), 
             "coverage for the simulated reads (>0 to simulate sequencing)")
            ("output-directory,d", po::value<std::string>(&seq_output_directory), 
             "the output directory for sequencing simulation")
            ("overwrite,o", "overwrite the output directory")
            ("read-size,r", po::value<size_t>(&read_size)->default_value(150),
             "simulated read size")
            ("insert-size,i", po::value<size_t>(&insert_size),
             "simulated insert size (enable paired-read)")
            ("write-SAM,w", "write the simulated read SAM files")
        ;

        visible_options.at("drivers").add_options()
            ("epigenetic-FACS,F", "distinguish between epigenetic states")
        ;

        visible_options.at("generic").add_options()
            ("seed,s", po::value<int>(&seed)->default_value(0), 
             "random generator seed")
            ("quiet,q", "avoid output messages")
            ("help,h", "get the help")
        ;

        hidden_options.add_options()
            ("simulation file", po::value<std::filesystem::path>(&simulation_filename), 
             "the name of the file describing the passenger mutations simulation")
            ("driver simulation", po::value<std::filesystem::path>(&drivers_directory), 
             "the driver simulation directory")
            ("context index", po::value<std::filesystem::path>(&context_index_filename), 
             "the genome context index")
            ("mutational signature", po::value<std::filesystem::path>(&SBS_filename), 
             "the mutational signature")
            ("reference genome", po::value<std::filesystem::path>(&ref_genome_filename), 
             "the filename of the reference genome FASTA file")
        ;

        po::options_description program_options;
        for (const auto& [section_name, section]: visible_options) {
            program_options.add(section);
        }
        program_options.add(hidden_options);

        positional_options.add("simulation file", 1);
        positional_options.add("driver simulation", 1);
        positional_options.add("context index", 1);
        positional_options.add("mutational signature", 1);
        positional_options.add("reference genome", 1);

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

        quiet = vm.count("quiet")>0;

        {
            using namespace Races::Passengers::SequencingSimulations;

            read_simulator_output_mode = ((vm.count("overwrite")>0)?
                                          ReadSimulator<>::Mode::OVERWRITE:
                                          ReadSimulator<>::Mode::CREATE);
        }

        paired_read = vm.count("insert-size")>0;

        epigenetic_FACS = vm.count("epigenetic-FACS")>0;
        write_SAM = vm.count("write-SAM")>0;

        try {
            using namespace Races::Passengers;

            Races::Archive::Binary::In archive(context_index_filename);

            bytes_per_abs_position = ContextIndex<>::read_bytes_per_absolute_position(archive);

        } catch (std::exception& except) {
            print_help_and_exit(except.what(), 1);
        }

        validate(vm);

        snapshot_path = get_last_snapshot_path(drivers_directory, "driver simulation");
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
