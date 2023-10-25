/**
 * @file read_simulator.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements classes to simulate sequencing
 * @version 0.6
 * @date 2023-10-25
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
#include <regex>

#include "read_simulator.hpp"

namespace Races
{

namespace Passengers
{

namespace SequencingSimulations
{

void Statistics::SampleStatistics::increase_coverage(const ChromosomeId& chromosome_id, const ChrPosition& begin_pos, const size_t& read_size)
{
    auto& chr_coverage = coverage[chromosome_id];

    for (auto base=chr_coverage.begin()+begin_pos; base != chr_coverage.begin()+begin_pos+read_size; ++base) {
        ++(*base);
    }
}

void Statistics::SampleStatistics::increase_SNV_occurrences(const SNV& snv)
{
    auto SNV_occ_it = SNV_occurrences.find(snv);

    if (SNV_occ_it == SNV_occurrences.end()) {
        SNV_occurrences.insert({snv, 1});

        return;
    }

    ++(SNV_occ_it->second);
}

Statistics::SampleStatistics Statistics::get_total() const
{
    Statistics::SampleStatistics total;

    for (const auto& [sample_name, statistics]: sample_statistics) {
        for (const auto& [chr_id, coverage]: statistics.coverage) {
            total.add_chromosome(chr_id, coverage.size());
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

    for (const auto& [snv, occurrences]: total.SNV_occurrences) {
        for (const auto& [sample_name, statistics]: sample_statistics) {
            auto& total_chr_coverage = total.coverage.at(snv.chr_id);
            const auto& sample_chr_coverage = statistics.coverage.at(snv.chr_id);

            total_chr_coverage[snv.position] += sample_chr_coverage[snv.position];
        }
    }

    return total;
}

const uint16_t& Statistics::SampleStatistics::find_coverage(const GenomicPosition& position) const
{
    auto it = coverage.find(position.chr_id);
    if (it == coverage.end()) {
        throw std::runtime_error("Unknown chromosome "+std::to_string(position.chr_id));
    }
    return it->second[position.position];
}

void print_sample_data_about(std::ofstream& os, const Statistics::SampleStatistics& statistics, 
                             const SNV& snv, const char separator='\t')
{
    const auto& coverage = statistics.find_coverage(snv);

    auto found = statistics.SNV_occurrences.find(snv);
    if (found != statistics.SNV_occurrences.end()) {
        const auto& occurrences = found->second;

        os << separator << occurrences
           << separator << coverage 
           << separator << (static_cast<double>(occurrences)/coverage);
    } else {
        os << separator << "0" << separator << coverage << separator << "0";
    }
}

void Statistics::save_VAF_csv(const std::filesystem::path& filename) const
{
    const char separator='\t';

    std::ofstream os(filename);

    os << "chr" << separator << "from" << separator << "to" << separator
        << "ref" << separator << "alt" << separator
        << "type" << separator << "context" << separator << "cause";
    
    Statistics::SampleStatistics total = get_total();
    if (sample_statistics.size()>1) {
        for (const auto& [sample_name, statistics]: sample_statistics) {
            os << separator << sample_name << " occurrences"
               << separator << sample_name << " coverage" 
               << separator << sample_name << " VAF";
        }
        os << separator << "total occurrences"
           << separator << "total coverage" 
           << separator << "total VAF" << std::endl;
    } else {
        os << separator << "occurrences"
           << separator << "coverage" 
           << separator << "VAF" << std::endl;
    }

    for (const auto& [snv, total_occurrences]: total.SNV_occurrences) {
        os << GenomicPosition::chrtos(snv.chr_id) << separator << snv.position 
           << separator << snv.position 
           << separator << snv.context.get_central_nucleotide()
           << separator << snv.mutated_base << separator << "SNV" 
           << separator << snv.context << separator << snv.cause;

        for (const auto& [sample_name, statistics]: sample_statistics) {
            print_sample_data_about(os, statistics, snv, separator);
        }

        if (sample_statistics.size()>1) {
            print_sample_data_about(os, total, snv, separator);
        }

        os << std::endl;
    }
}


#ifdef WITH_MATPLOT

std::vector<double> down_sample_coverage(const std::vector<uint16_t>& coverage, const size_t& new_size=300)
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

    size_t num_of_plots = sample_statistics.size();
    const size_t new_size = 300;

    auto f = figure(true);
    f->size(1500, 1500);

    if (num_of_plots==1) {
        const auto& coverage = (sample_statistics.begin()->second).coverage.find(chromosome_id)->second;
        
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
        for (const auto& [sample_name, statistics]: sample_statistics) {
            nexttile();

            const auto& coverage = statistics.coverage.find(chromosome_id)->second;

            chr_size = coverage.size();
            
            std::vector<std::vector<double>> d_coverage = {down_sample_coverage(coverage, new_size)};

            area(d_coverage);

            const auto latex_name = std::regex_replace(sample_name, std::regex("_"), "\\\\_");
            title("Sample \""+latex_name+"\"");

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

    for (const auto& [sample_name, g_statistics]: statistics.sample_statistics) {
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
    size_t num_of_plots = sample_statistics.size();

    auto f = figure(true);
    f->size(1500, 1500);
    
    if (num_of_plots==1) {
        const auto& statistics = (sample_statistics.begin()->second);
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
        for (const auto& [sample_name, statistics]: sample_statistics) {
            nexttile();

            const auto& coverage = statistics.coverage.find(chromosome_id)->second;

            std::vector<double> VAF;

            for (const auto& [snv, occurrences]: statistics.SNV_occurrences) {
                if (snv.chr_id == chromosome_id) {
                    VAF.push_back(static_cast<double>(occurrences)/coverage[snv.position]);
                }
            }

            hist(VAF,num_of_bins);

            const auto latex_name = std::regex_replace(sample_name, std::regex("_"), "\\\\_");
            title("Sample \""+latex_name+"\"");

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
