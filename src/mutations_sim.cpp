/**
 * @file mutations_sim.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Main file for the RACES mutations simulator
 * @version 0.3
 * @date 2023-12-12
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
boost::lexical_cast<std::string, Races::Mutations::SequencingSimulations::ReadSimulator<>::Mode>(const Races::Mutations::SequencingSimulations::ReadSimulator<>::Mode& mode)
{
    using namespace Races::Mutations::SequencingSimulations;
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
Races::Mutations::SequencingSimulations::ReadSimulator<>::Mode
boost::lexical_cast<Races::Mutations::SequencingSimulations::ReadSimulator<>::Mode, std::string>(const std::string& mode)
{
    std::string upper_token = mode;
    for (auto& t_char: upper_token) {
        t_char = toupper(t_char);
    }

    using namespace Races::Mutations::SequencingSimulations;

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

class MutationsSimulator : public BasicExecutable
{
    std::filesystem::path simulation_filename;
    std::filesystem::path species_directory;
    std::filesystem::path snapshot_path;
    std::filesystem::path context_index_filename;
    std::filesystem::path ref_genome_filename;
    std::filesystem::path SBS_filename;

    double coverage;
    std::string seq_output_directory;
    Races::Mutations::SequencingSimulations::ReadSimulator<>::Mode read_simulator_output_mode;
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

    std::list<Races::Mutants::Evolutions::TissueSample>
    get_samples(const Races::Mutants::Evolutions::Simulation& simulation, const nlohmann::json& simulation_cfg) const
    {
        using namespace Races::Mutants;
        using namespace Races::Mutants::Evolutions;

        std::list<Races::Mutants::Evolutions::TissueSample> samples;

        if (simulation_cfg.contains("sample regions")) {
            const auto& sample_regions_json = simulation_cfg["sample regions"];
            if (!sample_regions_json.is_array()) {
                throw std::domain_error("The \"sample regions\" field must be an array.");
            }

            for (const auto& sample_region_json: sample_regions_json) {
                auto sample_specification = Races::ConfigReader::get_sample_specification(sample_region_json);
                samples.push_back(simulation.simulate_sampling(sample_specification.get_name(),
                                                               sample_specification.get_region()));
            }

            return samples;
        }

        return simulation.get_tissue_samples();
    }

    template<typename ABSOLUTE_GENOTYPE_POSITION = uint32_t>
    Races::Mutations::ContextIndex<ABSOLUTE_GENOTYPE_POSITION> load_context_index(const std::string& filename) const
    {
        Races::Archive::Binary::In archive(filename);

        Races::Mutations::ContextIndex<ABSOLUTE_GENOTYPE_POSITION> context_index;

        if (quiet) {
            archive & context_index;
        } else {
            archive.load(context_index, "context index");
        }

        return context_index;
    }

    static std::map<std::string, Races::Mutations::MutationalSignature> load_signatures(const std::string filename)
    {
        std::ifstream is(filename);

        return Races::Mutations::MutationalSignature::read_from_stream(is);
    }

    template<typename ABSOLUTE_GENOME_POSITION, typename RANDOM_GENERATOR>
    std::list<Races::Mutations::SampleGenomeMutations>
    place_mutations(Races::Mutations::MutationEngine<ABSOLUTE_GENOME_POSITION,RANDOM_GENERATOR>& engine,
                    const Races::Mutants::DescendantsForest& forest) const
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

    void process_statistics(const std::list<Races::Mutations::SampleGenomeMutations>& mutations_list) const
    {
        using namespace Races::UI;
        using namespace Races::Mutations;

        if ((SNVs_csv_filename != "") || (CNAs_csv_filename != "")) {

            ProgressBar bar;

            MutationStatistics statistics;

            statistics.record(mutations_list, bar);

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

    static void
    split_epigenetic_status(std::list<Races::Mutations::SampleGenomeMutations>& FACS_samples,
                            const Races::Mutations::SampleGenomeMutations& sample_mutations,
                            std::map<Races::Mutants::SpeciesId, std::string> methylation_map)
    {
        using namespace Races::Mutants;
        using namespace Races::Mutants::Evolutions;
        using namespace Races::Mutations;

        std::map<SpeciesId, SampleGenomeMutations*> meth_samples;

        for (const auto& cell_mutations : sample_mutations.mutations) {
            auto found = meth_samples.find(cell_mutations.get_species_id());

            if (found == meth_samples.end()) {
                auto new_name = sample_mutations.get_name()+"_"+
                                    methylation_map.at(cell_mutations.get_species_id());

                FACS_samples.push_back(SampleGenomeMutations(sample_mutations));

                FACS_samples.back().mutations.push_back(cell_mutations);

                meth_samples.insert({cell_mutations.get_species_id(), &(FACS_samples.back())});
            } else {
                (found->second)->mutations.push_back(cell_mutations);
            }
        }
    }

    static std::list<Races::Mutations::SampleGenomeMutations>
    split_epigenetic_status(const std::list<Races::Mutations::SampleGenomeMutations>& sample_mutations_list,
                            const std::map<Races::Mutants::SpeciesId, std::string>& methylation_map)
    {
        using namespace Races::Mutants;
        using namespace Races::Mutants::Evolutions;

        std::list<Races::Mutations::SampleGenomeMutations> FACS_samples;

        for (const auto& sample_mutations : sample_mutations_list) {
            split_epigenetic_status(FACS_samples, sample_mutations, methylation_map);
        }

        return FACS_samples;
    }

    template<typename TISSUE_SAMPLE,
             std::enable_if_t<std::is_base_of_v<Races::Mutants::Evolutions::TissueSample, TISSUE_SAMPLE>, bool> = true>
    static std::map<Races::Time, std::list<TISSUE_SAMPLE>>
    split_list_by_time(const std::list<TISSUE_SAMPLE>& sample_list)
    {
        std::map<Races::Time, std::list<TISSUE_SAMPLE>> time_map;

        for (const auto& sample : sample_list) {
            time_map[sample.get_time()].push_back(sample);
        }

        return time_map;
    }

    template<typename GENOME_WIDE_POSITION>
    static std::map<Races::Mutations::ChromosomeId, size_t>
    get_number_of_alleles(const Races::Mutations::ContextIndex<GENOME_WIDE_POSITION>&context_index,
                          const nlohmann::json& simulation_cfg)
    {
        using namespace Races::Mutations;

        auto chromosome_ids = context_index.get_chromosome_ids();
        return Races::ConfigReader::get_number_of_alleles(simulation_cfg, chromosome_ids);
    }

    template<typename GENOME_WIDE_POSITION>
    void run_abs_position() const
    {
        using namespace Races;
        using namespace Races::Mutants;
        using namespace Races::Mutations;
        using namespace Races::Mutations::SequencingSimulations;

        nlohmann::json simulation_cfg = get_simulation_json();

        if (!quiet) {
            if (simulation_cfg.contains("sample region")) {
                std::cout << "Using \"sample region\" configuration field as sample region." << std::endl;
            } else {
                std::cout << "Using pre-collected samples." << std::endl;
            }
        }

        ReadSimulator<> read_simulator;

        Mutants::DescendantsForest forest;
        MutationalProperties mutational_properties;

        std::map<SpeciesId, std::string> methylation_map;
        {
            auto species_simulation = load_species_simulation(snapshot_path, quiet);

            for (const auto& species : species_simulation.tissue()) {
                const auto& signature = species.get_methylation_signature();
                methylation_map[species.get_id()] = MutantProperties::signature_to_string(signature);
            }

            auto samples = get_samples(species_simulation, simulation_cfg);

            forest = Mutants::DescendantsForest(species_simulation, samples);

            mutational_properties = ConfigReader::get_mutational_properties(species_simulation,
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

        auto num_of_alleles = get_number_of_alleles(context_index, simulation_cfg);

        MutationEngine<GENOME_WIDE_POSITION> engine(context_index, num_of_alleles,
                                                    default_coefficients, signatures,
                                                    mutational_properties);

        ConfigReader::add_timed_mutational_coefficients(engine, mutational_coeff_json);

        auto mutations_list = place_mutations(engine, forest);

        if (epigenetic_FACS) {
            mutations_list = split_epigenetic_status(mutations_list, methylation_map);
        }

        if (coverage>0) {
            read_simulator.enable_SAM_writing(write_SAM);

            read_simulator(mutations_list, coverage);
        }

        process_statistics(mutations_list);
    }

    void validate(boost::program_options::variables_map vm) const
    {
        namespace fs = std::filesystem;

        if (!vm.count("mutation file")) {
            print_help_and_exit("The mutation file is mandatory", 1);
        }

        if (!vm.count("species simulation")) {
            print_help_and_exit("The species simulation directory is mandatory. "
                                "You can produce it by using `mutants_sim`", 1);
        }

        if (!fs::exists(species_directory)) {
            print_help_and_exit("\"" + std::string(species_directory) + "\"  does not exist", 1);
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
    MutationsSimulator(int argc, char* argv[]):
        BasicExecutable(argv[0], {{"mutations", "Mutation simulation options"},
                                  {"sequencing", "Sequencing simulation options"},
                                  {"mutants", "Mutants evolution related options"},
                                  {"generic", "Generic options"}}),
        paired_read(false),
        insert_size(0), SNVs_csv_filename(""), CNAs_csv_filename(""), epigenetic_FACS(false)
    {
        namespace po = boost::program_options;

        using namespace Races::Mutations::SequencingSimulations;

        visible_options.at("mutations").add_options()
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

        visible_options.at("mutants").add_options()
            ("epigenetic-FACS,F", "distinguish between epigenetic states")
        ;

        visible_options.at("generic").add_options()
            ("seed,s", po::value<int>(&seed)->default_value(0),
             "random generator seed")
            ("quiet,q", "avoid output messages")
            ("help,h", "get the help")
        ;

        hidden_options.add_options()
            ("mutation file", po::value<std::filesystem::path>(&simulation_filename),
             "the name of the file describing the mutations simulation")
            ("species simulation", po::value<std::filesystem::path>(&species_directory),
             "the species simulation directory")
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

        positional_options.add("mutation file", 1);
        positional_options.add("species simulation", 1);
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
            using namespace Races::Mutations::SequencingSimulations;

            read_simulator_output_mode = ((vm.count("overwrite")>0)?
                                          ReadSimulator<>::Mode::OVERWRITE:
                                          ReadSimulator<>::Mode::CREATE);
        }

        paired_read = vm.count("insert-size")>0;

        epigenetic_FACS = vm.count("epigenetic-FACS")>0;
        write_SAM = vm.count("write-SAM")>0;

        try {
            using namespace Races::Mutations;

            Races::Archive::Binary::In archive(context_index_filename);

            bytes_per_abs_position = ContextIndex<>::read_bytes_per_absolute_position(archive);

        } catch (std::exception& except) {
            print_help_and_exit(except.what(), 1);
        }

        validate(vm);

        snapshot_path = get_last_snapshot_path(species_directory, "species simulation");
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
    MutationsSimulator passenger_sim(argc, argv);

    passenger_sim.run();

    return 0;
}
