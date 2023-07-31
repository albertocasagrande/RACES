/**
 * @file sample_context_index.cpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Main file for the simulator
 * @version 0.1
 * @date 2023-07-31
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

#include <boost/program_options.hpp>

#include "context_index.hpp"


std::ostream& print_help(std::ostream& os, const std::string& program_name, 
                         const boost::program_options::options_description& options)
{
   os << "Syntax: "<< program_name 
      << " <context index filename> <output filename> <sampling ratio>" << std::endl
      << options << std::endl;

   return os;   
}

Races::Passengers::ContextIndex<> load_context_index(const std::string filename)
{
    Races::Archive::Binary::In archive(filename);

    return Races::Passengers::ContextIndex<>::load(archive);
}

boost::program_options::variables_map handle_options(int argc, char* argv[])
{
    namespace po = boost::program_options;

    po::options_description visible("Options");
    visible.add_options()
#ifdef WITH_INDICATORS
        ("quit,q", "disable progress bar")
#endif // WITH_INDICATORS
        ("help,h", "get the help");
    
    std::string context_index_filename;
    std::string output_filename;
    unsigned int sampling_ratio;

    po::options_description hidden("Hidden options");
    hidden.add_options()
        ("context index filename", po::value<std::string>(&context_index_filename), 
         "the genome context index filename")
        ("output filename", po::value<std::string>(&output_filename), 
         "the resulting genome context index filename")
        ("sampling ratio", po::value<unsigned int>(&sampling_ratio), 
         "the genome context index filename")
    ;

    po::options_description generic;
    generic.add(visible).add(hidden);

    po::positional_options_description p;
    p.add("context index filename", 1);
    p.add("output filename", 1);
    p.add("sampling ratio", 1);

    po::variables_map vm;
    try {
        po::command_line_parser parser{argc, argv};
        po::store(parser.options(generic).positional(p).run(), vm);
        po::notify(vm);
    } catch (boost::wrapexcept<boost::program_options::unknown_option> &except) {
        std::cout << except.what() << std::endl;
        print_help(std::cerr, argv[0], visible);
        exit(1);
    }

    if (!vm.count("context index filename")) {
        std::cerr << "The context index filename is mandatory. "
                  << "You can produce it by using `races_build_index`" << std::endl << std::endl;
        print_help(std::cerr, argv[0], visible);
        exit(1);
    }

    if (!std::filesystem::exists(context_index_filename)) {
        std::cerr << "\"" << context_index_filename << "\"  does not exists."
                  << std::endl << std::endl;
        print_help(std::cerr, argv[0], visible);
        exit(1);
    }

    if (!std::filesystem::is_regular_file(context_index_filename)) {
        std::cerr << "\"" << context_index_filename << "\"  is not a regular file."
                  << std::endl << std::endl;
        print_help(std::cerr, argv[0], visible);
        exit(1);
    }

    if (!vm.count("output filename")) {
        std::cerr << "The output context index filename is mandatory." << std::endl << std::endl;
        print_help(std::cerr, argv[0], visible);
        exit(1);
    }

    if (std::filesystem::exists(output_filename)) {
        std::cerr << "\"" << output_filename << "\"  does exists."
                  << std::endl << std::endl;
        print_help(std::cerr, argv[0], visible);
        exit(1);
    }

    if (!vm.count("sampling ratio")) {
        std::cerr << "The sampling ratio is mandatory." << std::endl << std::endl;
        print_help(std::cerr, argv[0], visible);
        exit(1);
    }

    return vm;
}

int main(int argc, char* argv[])
{
    using namespace Races;
    using namespace Races::Passengers;
    namespace po = boost::program_options;

    auto vm = handle_options(argc, argv);

    if (vm.count("quiet")==0) {
        std::cout << "Loading context index..." << std::flush;
    }
    auto context_index = load_context_index(vm["context index filename"].as<std::string>());

    unsigned int sampling_ratio = vm["sampling ratio"].as<unsigned int>();

    size_t total_removal{0}, removed{0}, percentage{0};

    UI::ProgressBar* bar{nullptr};
    if (vm.count("quiet")==0) {
        std::cout << "done" << std::endl << std::flush;

        UI::ProgressBar::hide_console_cursor();

        bar = new UI::ProgressBar();
        bar->set_message("Sampling context index");
    }

    for (auto& [context, positions] : context_index.get_context_positions()) {
        total_removal += positions.size()*(sampling_ratio-1)/sampling_ratio;
    }

    std::mt19937_64 generator(0);
    for (auto& [context, positions] : context_index.get_context_positions()) {
        size_t aimed_size = positions.size()/sampling_ratio;
        while (positions.size() > aimed_size) {
            std::uniform_int_distribution<> u_dist(0,positions.size()-1);

            context_index.extract(context, u_dist(generator));
            ++removed;

            percentage = 100*removed/total_removal;
            if (percentage<99 && percentage > bar->get_progress()) {
                bar->set_progress(percentage);
            }
        }
    }

    if (vm.count("quiet")==0) {
        bar->set_progress(100, "Done");

        UI::ProgressBar::show_console_cursor();

        delete bar;

        std::cout << "Saving context index..." << std::flush;
    }

    {
        Races::Archive::Binary::Out archive(vm["output filename"].as<std::string>());

        archive & context_index;
    }
    
    if (vm.count("quiet")==0) {
        std::cout << "done" << std::endl << std::flush;
    }

    return 0;
}
