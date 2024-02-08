/**
 * @file mutation_engine.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements a class to place mutations on a descendants forest
 * @version 0.12
 * @date 2024-02-08
 *
 * @copyright Copyright (c) 2023-2024
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

namespace Mutations
{

MutationStatistics::SNVStatistics::SNVStatistics():
    mutated_alleles(0), num_of_cells(0)
{}

MutationStatistics::SampleMutationStatistics::SampleMutationStatistics():
    num_of_cells(0)
{}

MutationStatistics::MutationStatistics()
{}

MutationStatistics& MutationStatistics::record(const std::string& sample_name,
                                               const CellGenomeMutations& cell_mutations)
{
    std::set<SNV> in_cell;

    auto& sample_stat = sample_statistics[sample_name];

    // for all chromosomes
    for (const auto& [chr_id, chromosome]: cell_mutations.get_chromosomes()) {

        for (const auto& cna : chromosome.get_CNAs()) {
            sample_stat.CNAs.push_back(cna);
            overall_statistics.CNAs.push_back(cna);
        }

        // for all chromosome alleles
        for (const auto& [allele_id, allele]: chromosome.get_alleles()) {

            // for all fragments in the allele
            for (const auto& [fragment_pos, fragment]: allele.get_fragments()) {

                // for all SNVs in the fragment
                for (const auto& [snv_pos, snv]: fragment.get_SNVs()) {

                    ++(sample_stat.SNVs[snv].mutated_alleles);
                    ++(overall_statistics.SNVs[snv].mutated_alleles);

                    in_cell.insert(snv);
                }
            }
        }
    }

    ++(sample_stat.num_of_cells);
    ++(overall_statistics.num_of_cells);
    for (const auto& snv : in_cell) {
        ++(sample_stat.SNVs[snv].num_of_cells);
        ++(overall_statistics.SNVs[snv].num_of_cells);
    }

    return *this;
}

MutationStatistics& MutationStatistics::record(const std::list<Races::Mutations::SampleGenomeMutations>& mutations_list,
                                               UI::ProgressBar* progress_bar)
{
    using namespace Races;
    using namespace Races::Mutations;

    size_t total_steps{0};
    size_t recorded{0};

    if (progress_bar!=nullptr) {
        progress_bar->set_message("Collecting mutation data");

        for (const auto& sample_mutations : mutations_list) {
            total_steps += sample_mutations.mutations.size();
        }
    }

    for (const auto& sample_mutations : mutations_list) {
        const auto& sample_name = sample_mutations.name;
        for (const auto& cell_mutations: sample_mutations.mutations) {
            record(sample_name, cell_mutations);

            if (progress_bar!=nullptr) {
                size_t percentage = 100*(++recorded)/total_steps;
                if (percentage>(progress_bar->get_progress())) {
                    progress_bar->set_progress(percentage);
                }
            }
        }
    }

    if (progress_bar!=nullptr) {
        progress_bar->set_message("Mutation data collected");
    }

    return *this;
}

std::ostream& MutationStatistics::write_SNVs_table(std::ostream& os, const char separator)
{
    os << "chr" << separator << "from" << separator << "to" << separator
       << "ref" << separator << "alt" << separator
       << "type" << separator << "cause" << separator << "mutation type";

    if (sample_statistics.size()>1) {
        for (const auto& [sample_name, sample_stat]: sample_statistics) {
            os << separator << "# of mutated alleles in " << sample_name
               << separator << "# cells in " << sample_name;
        }
        os << separator << "# of mutated alleles in total"
           << separator << "# cells in total";
    } else {
        os << separator << "# of mutated alleles"
           << separator << "# of cells";
    }

    os << std::endl;

    for (const auto& [snv, snv_statistics] : overall_statistics.SNVs) {

        os << GenomicPosition::chrtos(snv.chr_id) << separator << snv.position
           << separator << snv.position << separator
           << snv.ref_base << separator << snv.alt_base << separator
           << "SNV" << separator << snv.cause << separator;

        switch(snv.type) {
            case SNV::Type::DRIVER:
                os << "D";
                break;
            case SNV::Type::PASSENGER:
                os << "P";
                break;
            case SNV::Type::GERMLINE:
                os << "G";
                break;
            default:
                os << "NA";
        };

        if (sample_statistics.size()>1) {
            for (const auto& [sample_name, sample_stat]: sample_statistics) {
                auto found = sample_stat.SNVs.find(snv);
                if (found != sample_stat.SNVs.end()) {
                    const auto &snv_stat = found->second;
                    os << separator << snv_stat.mutated_alleles
                       << separator << snv_stat.num_of_cells;
                } else {
                    os << separator << 0 << separator << 0;
                }
            }
        }

        os << separator << snv_statistics.mutated_alleles
           << separator << snv_statistics.num_of_cells
           << std::endl;
    }

    return os;
}

std::ostream& MutationStatistics::write_CNAs_table(std::ostream& os, const char separator)
{
    os << "chr" << separator << "from" << separator << "to" << separator
       << "type" << separator << "major" << separator << "minor" << std::endl;


    if (sample_statistics.size()>1) {
        for (const auto& cna: overall_statistics.CNAs) {
            os << GenomicPosition::chrtos(cna.region.get_chromosome_id()) << separator
               << cna.region.get_initial_position() << separator
               << cna.region.get_final_position()  << separator
               << (cna.type == CopyNumberAlteration::Type::AMPLIFICATION?"amplification":"deletion")
               << separator << "???" << separator << "???" << std::endl;
        }
    }

    return os;
}

}   // Mutations

}   // Races
