/**
 * @file passengers_sim.cpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Main file for the passenger mutations simulator
 * @version 0.1
 * @date 2023-08-13
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

    Races::Drivers::Simulation::Simulation load_simulation(const std::string& driver_directory) const
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

    static Races::Drivers::PhylogeneticForest 
    build_forest(const Races::Drivers::Simulation::Simulation& simulation,
                 const std::string& driver_directory,
                 const Races::Drivers::Simulation::PositionInTissue& lower_corner, 
                 const Races::Drivers::Simulation::PositionInTissue& upper_corner)
    {
        using namespace Races::Drivers;

        Simulation::BinaryLogger::CellReader reader(driver_directory);

        RectangleSampler sampler(simulation.tissue(),lower_corner,upper_corner);

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

    std::list<Races::Passengers::GenomeMutations>
    place_mutations(Races::Passengers::MutationEngine<uint32_t>& engine,
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

    template<typename ABSOLUTE_POSITION_TYPE>
    void run_abs_position() const
    {
        using namespace Races;
        using namespace Races::Passengers;
        using namespace Races::Passengers::IO;

        SAMGenerator<> sam_generator;

        if (coverage>0) {
            sam_generator = SAMGenerator<>(SAM_directory, ref_genome_filename, read_size, seed, 
                                           SAM_mode);
        }

        Drivers::PhylogeneticForest forest;
        SpeciesMutationalProperties mutational_properties;

        {
            auto simulation = load_simulation(drivers_directory);

            forest = build_forest(simulation, drivers_directory, {450,450},{550,550});

            {
                mutational_properties.add_species(simulation, "A", {{"-", 3e-9},{"+", 6e-9}}, 
                                                  {{1, 11266, "TGC", 'A'}, {2, 10693, "CAG", 'G'}});
                mutational_properties.add_species(simulation, "B", {{"-", 3e-8},{"+", 6e-8}});
            }
        }

        auto signatures = load_signatures(SBS_filename);

        auto context_index = load_context_index<uint32_t>(context_index_filename);

        std::map<std::string, double> default_mutational_coefficients{{"SBS3", 0.6},
                                                                      {"SBS13", 0.3},
                                                                      {"SBS40", 0.1}};

        MutationEngine<uint32_t> engine(context_index, 2, default_mutational_coefficients,
                                        signatures, mutational_properties);

        engine.add(70, {{"SBS3", 0.5},{"SBS13", 0.3}, {"SBS40", 0.2}});

        auto leaves_mutations = place_mutations(engine, forest);

        if (coverage) {
            sam_generator(leaves_mutations, coverage);
        }

        if ((SVNs_csv_filename != "") || (CNAs_csv_filename != "")) {
            auto statistics = collect_statistics(leaves_mutations);

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

    void validate(boost::program_options::variables_map vm) const
    {
        namespace fs = std::filesystem;

        if (vm.count("help")) {
            print_help(0);
        }

        if (!vm.count("driver directory")) {
            std::cerr << "The driver directory is mandatory. " 
                    << "You can produce it by using `races_sim`" << std::endl << std::endl;
            print_help(1);
        }

        if (!fs::exists(drivers_directory)) {
            std::cerr << "\"" << drivers_directory << "\"  does not exist"
                      << std::endl << std::endl;
            print_help(1);
        }

        if (!fs::is_directory(drivers_directory)) {
            std::cerr << "\"" << drivers_directory << "\"  is not a directory"
                      << std::endl << std::endl;
            print_help(1);
        }

        try {
            find_last_snapshot(drivers_directory);
        } catch (std::runtime_error& ex) {
            std::cerr << ex.what() << std::endl << std::endl;
            print_help(1);
        }

        if (!vm.count("context index")) {
            std::cerr << "The context index is mandatory. "
                    << "You can produce it by using `races_build_index`" << std::endl << std::endl;
            print_help(1);
        }

        if (!fs::exists(context_index_filename)) {
            std::cerr << "\"" << context_index_filename << "\"  does not exist"
                    << std::endl << std::endl;
            print_help(1);
        }

        if (!fs::is_regular_file(context_index_filename)) {
            std::cerr << "\"" << context_index_filename << "\"  is not a regular file"
                    << std::endl << std::endl;
            print_help(1);
        }

        if (!vm.count("mutational signature")) {
            std::cerr << "The mutational signature is mandatory." << std::endl << std::endl;
            print_help(1);
        }

        if (!fs::exists(SBS_filename)) {
            std::cerr << "\"" << SBS_filename << "\"  does not exist"
                      << std::endl << std::endl;
            print_help(1);
        }

        if (!fs::is_regular_file(SBS_filename)) {
            std::cerr << "\"" << SBS_filename << "\"  is not a regular file"
                    << std::endl << std::endl;
            print_help(1);
        }

        if (!vm.count("reference genome")) {
            std::cerr << "The reference genome FASTA file is mandatory." << std::endl << std::endl;
            print_help(1);
        }

        if (!fs::exists(ref_genome_filename)) {
            std::cerr << "\"" << ref_genome_filename << "\"  does not exist"
                      << std::endl << std::endl;
            print_help(1);
        }

        if (!fs::is_regular_file(ref_genome_filename)) {
            std::cerr << "\"" << ref_genome_filename << "\"  is not a regular file"
                      << std::endl << std::endl;
            print_help(1);
        }

        if (vm.count("SNVs-csv")>0 && fs::exists(SVNs_csv_filename)) {
            std::cerr << "\"" << SVNs_csv_filename << "\" already exists" << std::endl << std::endl;
            print_help(1);
        }

        if (vm.count("CNAs-csv")>0 && fs::exists(CNAs_csv_filename)) {
            std::cerr << "\"" << CNAs_csv_filename << "\" already exists" << std::endl << std::endl;
            print_help(1);
        }

        if (coverage<0) {
            std::cerr << "coverage must be non-negative." << std::endl << std::endl;
            print_help(1);
        }
        
        for (const auto SAM_option : {"SAM-directory"}) {
            if ((coverage>0) && (vm.count(SAM_option)==0)) {
                std::cerr << "\"--"<< SAM_option << "\" is mandatory when coverage differs from 0." 
                        << std::endl << std::endl;
                print_help(1);
            }
        }
    }

public:

    void print_help(const int& err_code) const
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
            std::cout << except.what() << std::endl;
            print_help(1);
        }

        quiet = vm.count("quiet")>0;

        try {
            using namespace Races::Passengers;

            Races::Archive::Binary::In archive(context_index_filename);

            bytes_per_abs_position = ContextIndex<>::read_bytes_per_absolute_position(archive);

        } catch (std::exception& except) {
            std::cerr << except.what() << std::endl << std::endl;
            print_help(1);
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
            std::cerr << except.what() << std::endl << std::endl;
            print_help(1);
        }
    }
};

int main(int argc, char* argv[])
{
    PassengersSimulator passenger_sim(argc, argv);

    passenger_sim.run();

    return 0;
}
