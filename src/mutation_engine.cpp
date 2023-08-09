/**
 * @file mutation_engine.cpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Implements a class to place passenger mutations on the nodes of a phylogenetic forest
 * @version 0.1
 * @date 2023-08-09
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

MutationStatistics::MutationStatistics()
{}

MutationStatistics& MutationStatistics::record(const GenomeMutations& cell_mutations)
{
    // for all chromosomes
    for (const auto& [chr_id, chromosome]: cell_mutations.get_chromosomes()) {

        for (const auto& cna : chromosome.get_CNAs()) {
            CNAs.push_back(cna);
        }

        // for all chromosome alleles
        for (const auto& [allele_id, allele]: chromosome.get_alleles()) {

            // for all fragments in the allele
            for (const auto& [fragment_pos, fragment]: allele.get_fragments()) {

                // for all SNVs in the fragment
                for (const auto& [snv_pos, snv]: fragment.get_SNVs()) {

                    auto it = SNVs.find(snv);

                    if (it != SNVs.end()) {
                        ++(it->second);
                    } else {
                        SNVs[snv] = 1;   
                    }
                }
            }
        }
    }
 
    return *this;
}

std::ostream& MutationStatistics::write_SNVs_table(std::ostream& os, const char separator)
{
    os << "chr" << separator << "from" << separator << "to" << separator
       << "ref" << separator << "alt" << separator
       << "type" << separator << "context" << separator << "cause" << separator 
       << "multiplicity"<< std::endl;


    for (const auto& [snv, counter] : SNVs) {

        os << GenomicPosition::chrtos(snv.chr_id) << separator << snv.position 
           << separator << snv.position << separator
           << snv.context.get_central_nucleotide() << separator << snv.mutated_base << separator
           << "SNV" << separator << snv.context << separator
           << snv.cause << separator << counter << std::endl;
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
