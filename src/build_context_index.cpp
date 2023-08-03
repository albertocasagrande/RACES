/**
 * @file build_context_index.cpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Builds the context index
 * @version 0.1
 * @date 2023-08-03
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

#include "genome_mutations.hpp"
#include "context_index.hpp"
#include "progress_bar.hpp"

enum class SeqNameDecCode
{
    ENSEMBL_DECODER,
    NCBI_DECODER
};

template<>
std::string boost::lexical_cast<std::string>(const SeqNameDecCode& code)
{
    switch(code) {
        case SeqNameDecCode::ENSEMBL_DECODER:
            return "Ensembl";
        case SeqNameDecCode::NCBI_DECODER:
            return "NCBI";
        default:
            throw std::domain_error("Unknown sequence name decoder. "
                                    "Supported decoders are \"Ensembl\" and \"NCBI\".");
    }
}

void validate(boost::any& v, std::vector<std::string> const& values, 
                SeqNameDecCode*, int)
{
    using namespace boost::program_options;

    // Make sure no previous assignment to 'v' was made.
    validators::check_first_occurrence(v);

    // Extract the first string from 'values'. If there is more than
    // one string, it's an error, and exception will be thrown.
    std::string const& s = validators::get_single_string(values);

    if (s == "Ensembl") {
        v = boost::any(SeqNameDecCode::ENSEMBL_DECODER);

        return;
    }

    if (s == "NCBI") {
        v = boost::any(SeqNameDecCode::NCBI_DECODER);

        return;
    }

    throw validation_error(validation_error::invalid_option_value);
}

class IndexBuilder
{   
    const std::string program_name;
    boost::program_options::options_description visible_options;
    boost::program_options::positional_options_description positional_options;

    std::string output_filename;
    std::string genome_fasta_filename;
    unsigned int bits_per_abs_position;
    SeqNameDecCode seq_name_decoder;
    bool quiet;

    template<typename ABSOLUTE_GENOMIC_POSITION, typename SEQUENCE_NAME_DECODER = Races::IO::FASTA::EnsemblSeqNameDecoder>
    std::vector<Races::Passengers::GenomicRegion> build_and_save_context_index() const
    {
        using namespace Races;
        using namespace Races::Passengers;

        std::vector<GenomicRegion> chr_regions;
        {
            ContextIndex<ABSOLUTE_GENOMIC_POSITION> context_index;

            if (quiet) {
                context_index = ContextIndex<ABSOLUTE_GENOMIC_POSITION>::template build_index<SEQUENCE_NAME_DECODER>(genome_fasta_filename);
            } else {
                UI::ProgressBar::hide_console_cursor();

                UI::ProgressBar progress_bar;

                context_index = ContextIndex<ABSOLUTE_GENOMIC_POSITION>::template build_index<SEQUENCE_NAME_DECODER>(genome_fasta_filename, &progress_bar);
            }
            
            chr_regions = context_index.get_chromosome_regions();

            if (chr_regions.size() > 0) {
                Archive::Binary::Out archive(output_filename);

                archive.save(context_index, quiet);
            }

            if (!quiet) {
                std::cout << " Cleaning memory..." << std::flush;
            }
        }

        if (!quiet) {
            UI::ProgressBar::show_console_cursor();

            std::cout << "done" << std::endl;

            if (chr_regions.size()==0) {
                std::cout << " No chromosome processed. Try to select the proper "
                          << "sequence name decoder by using \"-s\"" << std::endl;
            } else {
                std::cout << " Processed:" << std::endl;
                for (const auto& chr_region: chr_regions) {
                    std::cout << "  - Chromosome " 
                                << GenomicPosition::chrtos(chr_region.get_chromosome_id()) 
                                << " (size: " << chr_region.size() << " bps)" << std::endl;
                }
            }
        }

        return chr_regions;
    }

    template<typename SEQUENCE_NAME_DECODER>
    std::vector<Races::Passengers::GenomicRegion> select_bit_per_abs_position() const
    {
        switch (bits_per_abs_position) {
            case 2:
                return build_and_save_context_index<uint16_t, SEQUENCE_NAME_DECODER>();
            case 4:
                return build_and_save_context_index<uint32_t, SEQUENCE_NAME_DECODER>();
            case 8:
                return build_and_save_context_index<uint64_t, SEQUENCE_NAME_DECODER>();
            default:
                std::cerr << "Unsupported bits per absolute position."
                        << " Supported values are 2, 4, or 8." << std::endl<< std::endl;
                print_help(std::cerr);

                exit(1);
        }
    }

public:

    std::ostream& print_help(std::ostream& os) const
    {
        os << "Syntax: " << program_name;
        
        for (unsigned int i=0;i<positional_options.max_total_count(); ++i) {
            os << " <" << positional_options.name_for_position(i) << ">"; 
        }

        os << std::endl << visible_options << std::endl;

        return os;   
    }

    IndexBuilder(const int argc, const char* argv[]):
        program_name(argv[0]), visible_options("Options")
    {
        namespace po = boost::program_options;

        visible_options.add_options()
            ("output-filename,o", po::value<std::string>(&output_filename), 
            "output filename")
            ("force-overwrite,f", "force overwriting output file")
            ("bytes-per-pos,b", po::value<unsigned int>(&bits_per_abs_position)->default_value(4),
            "bytes per absolute position (2, 4, or 8)")
            ("seq-name,s", po::value<SeqNameDecCode>(&seq_name_decoder)->default_value(SeqNameDecCode::ENSEMBL_DECODER),
            "sequence name decoder (i.e., \"Ensembl\" or \"NCBI\")")
#ifdef WITH_INDICATORS
            ("quiet,q", "disable output messages")
#endif // WITH_INDICATORS
            ("help,h", "get the help");

        po::options_description hidden("Hidden options");
        hidden.add_options()
            ("genome filename", po::value<std::string>(&genome_fasta_filename), 
            "the path to the genome in FASTA file format")
        ;

        po::options_description generic;
        generic.add(visible_options).add(hidden);

        positional_options.add("genome filename", 1);

        po::variables_map vm;
        try {
            po::command_line_parser parser{argc, argv};
            po::store(parser.options(generic).positional(positional_options).run(), vm);
            po::notify(vm);
        } catch (boost::wrapexcept<boost::program_options::unknown_option> &except) {
            std::cerr << except.what() << std::endl<< std::endl;
            print_help(std::cerr);
            exit(1);
        }

        if (vm.count("help")) {
            print_help(std::cout);
            exit(0);
        }

        quiet = vm.count("quit")>0;

        if (!vm.count("genome filename")) {
            std::cerr << "Missing genome FASTA filename" << std::endl << std::endl;
            print_help(std::cerr);
            exit(1);
        }

        if (!vm.count("output-filename")) {
            output_filename = std::filesystem::path(genome_fasta_filename + ".cif").filename();
        }

        if (std::ifstream(output_filename).good() && !vm.count("force-overwrite")) {
            std::cerr << "The output file \"" << output_filename << "\" already exists. " 
                    << "Use \"--force-overwrite\" to overwrite it." << std::endl<< std::endl;
            print_help(std::cerr);
            exit(1);
        }

    }

    std::vector<Races::Passengers::GenomicRegion> run() const
    {
        using namespace Races::IO::FASTA;

        switch(seq_name_decoder) {
            case SeqNameDecCode::ENSEMBL_DECODER:
                return select_bit_per_abs_position<EnsemblSeqNameDecoder>();
            case SeqNameDecCode::NCBI_DECODER:
                return select_bit_per_abs_position<NCBISeqNameDecoder>();
            default:
                std::cerr << "Unsupported FASTA sequence name decoder"
                          << " Supported decoder are \"Ensembl\" and \"NCBI\"." << std::endl
                          << std::endl;
                print_help(std::cerr);

                exit(1);
        };
    }
};

int main(const int argc, const char* argv[])
{
    IndexBuilder builder(argc, argv);

    builder.run();

    return 0;
}
