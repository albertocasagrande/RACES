/**
 * @file read_simulator.cpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Implements classes to simulate sequencing
 * @version 0.1
 * @date 2023-09-16
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

#include "read_simulator.hpp"

namespace Races
{

namespace Passengers
{

namespace SequencingSimulations
{

void Statistics::GroupStatistics::increase_coverage(const ChromosomeId& chromosome_id, const ChrPosition& begin_pos, const size_t& read_size)
{
    auto& chr_coverage = coverage[chromosome_id];

    for (auto base=chr_coverage.begin()+begin_pos; base != chr_coverage.begin()+begin_pos+read_size; ++base) {
        ++(*base);
    }
}

void Statistics::GroupStatistics::increase_SNV_occurrences(const SNV& snv)
{
    auto SNV_occ_it = SNV_occurrences.find(snv);

    if (SNV_occ_it == SNV_occurrences.end()) {
        SNV_occurrences.insert({snv, 1});

        return;
    }

    ++(SNV_occ_it->second);
}

Statistics::GroupStatistics Statistics::get_total() const
{
    Statistics::GroupStatistics total;

    for (const auto& [group_name, statistics]: group_statistics) {
        for (const auto& [chr_id, coverage]: statistics.coverage) {
            total.add_chromosome(chr_id, coverage.size());

            auto& total_coverage = total.coverage[chr_id];

            auto total_it = total_coverage.begin();
            auto coverage_it = coverage.begin();

            for (; total_it != total_coverage.end(); ++total_it,++coverage_it) {
                *total_it += *coverage_it;
            }
        }

        for (const auto& [snv, occurrences]: statistics.SNV_occurrences) {
            auto SNV_occ_it = total.SNV_occurrences.find(snv);
            if (SNV_occ_it == total.SNV_occurrences.end()) {
                total.SNV_occurrences.insert({snv, occurrences});
            } else {
                SNV_occ_it->second += occurrences;
            }
        }
    }

    return total;
}

void Statistics::save_VAF_csv(const std::filesystem::path& filename) const
{
    const char separator='\t';

    const auto total = get_total();

    std::ofstream os(filename);

    os << "chr" << separator << "from" << separator << "to" << separator
        << "ref" << separator << "alt" << separator
        << "type" << separator << "context" << separator << "cause" << separator
        << "group coverage" << separator  << "group VAF" << separator
        << "total coverage" << separator << "total VAF" << separator 
        << "total occurrences";
    
    if (group_statistics.size()>1) {
        os << separator << "group" << std::endl;
    } else {
        os << std::endl;
    }

    for (const auto& [group_name, statistics]: group_statistics) {
        for (const auto& [snv, occurrences]: statistics.SNV_occurrences) {

            os << GenomicPosition::chrtos(snv.chr_id) << separator << snv.position 
                << separator << snv.position 
                << separator << snv.context.get_central_nucleotide()
                << separator << snv.mutated_base << separator << "SNV" 
                << separator << snv.context << separator << snv.cause;

            for (const auto coverage_ptr: {&(statistics.coverage), &(total.coverage)}) {
                auto it = coverage_ptr->find(snv.chr_id);
                if (it == coverage_ptr->end()) {
                    throw std::runtime_error("Unknown chromosome "+std::to_string(snv.chr_id));
                }
                const auto& coverage_value = it->second[snv.position];
                os << separator << static_cast<size_t>(coverage_value) << separator 
                    << static_cast<double>(occurrences)/coverage_value;
            }

            auto total_it = total.SNV_occurrences.find(snv);

            if (total_it == total.SNV_occurrences.end()) {
                throw std::runtime_error("Unknown SNV");
            }
            
            os << separator << total_it->second;

            if (group_statistics.size()>1) {
                os << separator << group_name << std::endl;
            } else {
                os << std::endl;
            }
        }
    }
}

} // SequencingSimulator

} // Passengers

} // Races
