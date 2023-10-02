/**
 * @file mutation_engine.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements a class to place passenger mutations on the nodes of a phylogenetic forest
 * @version 0.5
 * @date 2023-10-02
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

#include "mutation_engine.hpp"

namespace Races
{

namespace Passengers
{

MutationStatistics::MutationStatistics():
    num_of_cells(0), num_of_non_filtered_cells(0)
{}

std::ostream& MutationStatistics::write_SNVs_table(std::ostream& os, const char separator)
{
    os << "chr" << separator << "from" << separator << "to" << separator
       << "ref" << separator << "alt" << separator
       << "type" << separator << "context" << separator << "cause"
       << separator << "# of mutated alleles";

    if (num_of_non_filtered_cells != num_of_cells) { 
        os << separator << "# of non-filtered cells" << separator  << "# of cells"
           << separator << "non-filtered rate" << separator << "over-total rate";
    } else {
        os << separator  << "# of cells" << separator << "rate";
    }
    
    os << std::endl;

    for (const auto& [snv, snv_statistics] : SNVs) {

        os << GenomicPosition::chrtos(snv.chr_id) << separator << snv.position 
           << separator << snv.position << separator
           << snv.context.get_central_nucleotide() << separator << snv.mutated_base << separator
           << "SNV" << separator << snv.context << separator
           << snv.cause << separator << snv_statistics.mutated_alleles;

        const auto double_non_filtered = static_cast<double>(snv_statistics.num_of_non_filtered_cells);
        if (num_of_non_filtered_cells != num_of_cells) {
            os << separator << snv_statistics.num_of_non_filtered_cells 
               << separator << snv_statistics.num_of_cells
               << separator << (double_non_filtered/num_of_non_filtered_cells)
               << separator << (double_non_filtered/num_of_cells);
        } else {
            os << separator << snv_statistics.num_of_cells
               << separator << (double_non_filtered/num_of_cells);
        }

        os << std::endl;
    }

    return os;
}

std::ostream& MutationStatistics::write_CNAs_table(std::ostream& os, const char separator)
{
    os << "chr" << separator << "from" << separator << "to" << separator
       << "type" << separator << "major" << separator << "minor" << std::endl;

    for (const auto& cna : CNAs) {

        os << GenomicPosition::chrtos(cna.region.get_chromosome_id()) << separator 
           << cna.region.get_initial_position() << separator 
           << cna.region.get_final_position()  << separator
           << (cna.type == CopyNumberAlteration::Type::AMPLIFICATION?"amplification":"deletion")
           << separator << "???" << separator << "???" << std::endl;
    }

    return os;
}

}   // Passengers

}   // Races
