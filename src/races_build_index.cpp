/**
 * @file races_build_index.cpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Builds the context index
 * @version 0.1
 * @date 2023-07-29
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


std::ostream& print_help(std::ostream& os, const std::string& program_name, 
                         const boost::program_options::options_description& options)
{
   os << "Syntax: "<< program_name <<" <genome filename>" << std::endl
      << options << std::endl;

   return os;   
}

template<typename ABSOLUTE_GENOMIC_POSITION>
void build_and_save_context_index(const std::string& genome_fasta_filename, const std::string& output_filename, const bool quiet)
{
    using namespace Races;
    using namespace Races::Passengers;
    
    UI::ProgressBar* progress_bar{nullptr};
    if (!quiet) {
        progress_bar = new UI::ProgressBar();
    }
    auto context_index = ContextIndex<ABSOLUTE_GENOMIC_POSITION>::build_index(genome_fasta_filename, progress_bar);

    Races::Archive::Binary::Out archive(output_filename);

    if (!quiet) {
        progress_bar->set_message("Saving context index");
    }

    archive & context_index;
    
    if (!quiet) {
        progress_bar->set_message("Done");

        delete progress_bar;
    }
}

int main(int argc, char* argv[])
{
    namespace po = boost::program_options;
    
    std::string output_filename;
    std::string genome_fasta_filename;
    unsigned int bits_per_abs_position;

    po::options_description visible("Options");
    visible.add_options()
        ("output-filename,o", po::value<std::string>(&output_filename), 
         "output filename")
        ("force-overwrite,f", "force overwriting output file")
        ("bytes-per-pos,b", po::value<unsigned int>(&bits_per_abs_position)->default_value(4),
        "bytes per absolute position (2, 4, or 8)")
        ("quiet,q", "don't show the progress bar")
        ("help,h", "get the help");

    po::options_description hidden("Hidden options");
    hidden.add_options()
        ("genome filename", po::value<std::string>(&genome_fasta_filename), 
         "the path to the genome in FASTA file format")
    ;

    po::options_description generic;
    generic.add(visible).add(hidden);

    po::positional_options_description p;
    p.add("genome filename", 1);

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

    if (!vm.count("genome filename")) {
        std::cerr << "Missing genome FASTA filename" << std::endl << std::endl;
        print_help(std::cout, argv[0], visible);
        exit(1);
    }

    if (!vm.count("output-filename")) {
        output_filename = std::filesystem::path(genome_fasta_filename + ".cif").filename();
    }

    using namespace Races::Passengers;

    if (std::ifstream(output_filename).good() && !vm.count("force-overwrite")) {
        std::cerr << "The output file \"" << output_filename << "\" already exists. " 
                  << "Use \"--force-overwrite\" to force overwriting." << std::endl;
        exit(1);
    }

    switch (bits_per_abs_position) {
        case 2:
            build_and_save_context_index<uint16_t>(genome_fasta_filename, output_filename, vm.count("quiet"));
            break;
        case 4:
            build_and_save_context_index<uint32_t>(genome_fasta_filename, output_filename, vm.count("quiet"));
            break;
        case 8:
            build_and_save_context_index<uint64_t>(genome_fasta_filename, output_filename, vm.count("quiet"));
            break;
        default:
            std::cerr << "Not supported bits per absolute position."
                      << " Supported values are 2, 4, or 8." << std::endl;
            exit(1);
    }

    return 0;
}
