/**
 * @file read_simulator.cpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Implements classes to simulate sequencing
 * @version 0.2
 * @date 2023-09-18
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

#include <string>

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

    std::ofstream os(filename);

    os << "chr" << separator << "from" << separator << "to" << separator
        << "ref" << separator << "alt" << separator
        << "type" << separator << "context" << separator << "cause";
    
    Statistics::GroupStatistics total;
    if (group_statistics.size()>1) {
        total = get_total();

        os << separator << "group coverage" << separator  << "group VAF"
           << separator << "total coverage" << separator << "total VAF"
           << separator << "total occurrences" << separator << "group" << std::endl;
    } else {
        os << separator << "coverage" << separator  << "VAF" 
           << separator << "occurrences" << std::endl;
    }

    for (const auto& [group_name, statistics]: group_statistics) {
        for (const auto& [snv, occurrences]: statistics.SNV_occurrences) {

            os << GenomicPosition::chrtos(snv.chr_id) << separator << snv.position 
                << separator << snv.position 
                << separator << snv.context.get_central_nucleotide()
                << separator << snv.mutated_base << separator << "SNV" 
                << separator << snv.context << separator << snv.cause;

            std::vector<const decltype(statistics.coverage)*> coverage_ptr_vector = {&(statistics.coverage)};

            if (group_statistics.size()>1) {
                coverage_ptr_vector.push_back(&(total.coverage));
            }

            for (const auto coverage_ptr: coverage_ptr_vector) {
                auto it = coverage_ptr->find(snv.chr_id);
                if (it == coverage_ptr->end()) {
                    throw std::runtime_error("Unknown chromosome "+std::to_string(snv.chr_id));
                }
                const auto& coverage_value = it->second[snv.position];
                os << separator << static_cast<size_t>(coverage_value) << separator 
                    << static_cast<double>(occurrences)/coverage_value;
            }

            decltype(statistics.SNV_occurrences)::const_iterator total_it;

            if (group_statistics.size()>1) {
                total_it = total.SNV_occurrences.find(snv);
            } else {
                total_it = statistics.SNV_occurrences.find(snv);
            }

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


#ifdef WITH_MATPLOT

std::vector<double> down_sample_coverage(const std::vector<uint8_t>& coverage, const size_t& new_size=300)
{
    std::vector<double> down_sampled;
    
    down_sampled.reserve(new_size);

    auto bin_size = coverage.size()/new_size;
    size_t base_in_bins{0};
    double sum{0};
    for (const auto& base_coverage: coverage) {
        sum += base_coverage;
        ++base_in_bins;

        if (base_in_bins == bin_size) {
            down_sampled.push_back(sum/bin_size);

            base_in_bins = 0;
            sum = 0;
        }
    }

    return down_sampled;
}

void Statistics::save_coverage_images(const std::filesystem::path& filename, const ChromosomeId& chromosome_id) const
{
    using namespace matplot;

    size_t num_of_plots = group_statistics.size();
    const size_t new_size = 300;

    auto f = figure(true);
    f->size(1500, 1500);

    if (num_of_plots==1) {
        const auto& coverage = (group_statistics.begin()->second).coverage.find(chromosome_id)->second;
        
        std::vector<std::vector<double>> d_coverage = {down_sample_coverage(coverage, new_size)};

        area(d_coverage);

        long chr_size = coverage.size();
        xticks({1, new_size/2, new_size});
        xticklabels({"0",std::to_string(chr_size/2),std::to_string(chr_size)});

        xlabel("Pos. in sequence");
        ylabel("Coverage");
    } else {
        std::vector<std::vector<double>> total(1);

        tiledlayout(++num_of_plots, 1);

        long chr_size;
        for (const auto& [group_name, statistics]: group_statistics) {
            nexttile();

            const auto& coverage = statistics.coverage.find(chromosome_id)->second;

            chr_size = coverage.size();
            
            std::vector<std::vector<double>> d_coverage = {down_sample_coverage(coverage, new_size)};

            area(d_coverage);

            title("Read group \""+group_name+"\"");

            xticks({1, new_size/2, new_size});
            xticklabels({"0",std::to_string(chr_size/2),std::to_string(chr_size)});

            ylabel("Coverage");
            if (total[0].size() != d_coverage[0].size()) {
                total[0] = d_coverage[0];
            } else {
                auto total_it = total[0].begin();
                auto d_coverage_it = d_coverage[0].begin();

                for (; total_it != total[0].end(); ++total_it,++d_coverage_it) {
                    *total_it += *d_coverage_it;
                }
            }
        }

        nexttile();

        area(total);
        title("Overall");

        xticks({1, new_size/2, new_size});
        xticklabels({"0",std::to_string(chr_size/2),std::to_string(chr_size)});

        xlabel("Pos. in sequence");
        ylabel("Coverage");
    }

    sgtitle("Chromosome "+GenomicPosition::chrtos(chromosome_id)+"'s coverage");

    f->save(filename);
}

std::map<SNV, std::pair<double, size_t>> get_total_occurrences(const Statistics& statistics, const ChromosomeId& chromosome_id)
{
    std::map<SNV, std::pair<double, size_t>> total_occurrences;

    for (const auto& [group_name, g_statistics]: statistics.group_statistics) {
        const auto& coverage = g_statistics.coverage.find(chromosome_id)->second;
        
        for (const auto& [snv, occurrences]: g_statistics.SNV_occurrences) {
            if (snv.chr_id == chromosome_id) {
                auto total_it = total_occurrences.find(snv);
                if (total_it == total_occurrences.end()) {
                    total_occurrences[snv] = {occurrences, coverage[snv.position]};
                } else {
                    total_it->second.first += occurrences;
                    total_it->second.second += coverage[snv.position];
                };
            }
        }
    }

    return total_occurrences;
}

void Statistics::save_SNV_histogram(const std::filesystem::path& filename, const ChromosomeId& chromosome_id) const
{
    using namespace matplot;

    size_t num_of_bins = 50;
    size_t num_of_plots = group_statistics.size();

    auto f = figure(true);
    f->size(1500, 1500);
    
    if (num_of_plots==1) {
        const auto& statistics = (group_statistics.begin()->second);
        const auto& coverage = statistics.coverage.find(chromosome_id)->second;

        std::vector<double> VAF;

        for (const auto& [snv, occurrences]: statistics.SNV_occurrences) {
            if (snv.chr_id == chromosome_id) {
                VAF.push_back(static_cast<double>(occurrences)/coverage[snv.position]);
            }
        }

        hist(VAF,num_of_bins);

        xlabel("VAF");
        ylabel("Num. of SNVs");
    } else {
        auto total_occurrences = get_total_occurrences(*this, chromosome_id);

        tiledlayout(++num_of_plots, 1);
        for (const auto& [group_name, statistics]: group_statistics) {
            nexttile();

            const auto& coverage = statistics.coverage.find(chromosome_id)->second;

            std::vector<double> VAF;

            for (const auto& [snv, occurrences]: statistics.SNV_occurrences) {
                if (snv.chr_id == chromosome_id) {
                    VAF.push_back(static_cast<double>(occurrences)/coverage[snv.position]);
                }
            }

            hist(VAF,num_of_bins);

            title("Read group \""+group_name+"\"");

            ylabel("Num. of SNVs");
        }

        nexttile();

        std::vector<double> VAF;
        for (const auto& [snv, data]: total_occurrences) {
            VAF.push_back(data.first/data.second);
        }

        hist(VAF,num_of_bins);
        title("Overall");

        xlabel("VAF");
        ylabel("Num. of SNVs");

    }

    sgtitle("Chromosome "+GenomicPosition::chrtos(chromosome_id)+"'s SNV VAF");

    f->save(filename);
}

#endif // WITH_MATPLOT

} // SequencingSimulator

} // Passengers

} // Races
