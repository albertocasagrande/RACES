/**
 * @file read_simulator.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements classes to simulate sequencing
 * @version 0.13
 * @date 2024-01-07
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

#include <string>
#include <regex>

#include "read_simulator.hpp"

namespace Races
{

namespace Mutations
{

namespace SequencingSimulations
{

ChrSampleStatistics::ChrSampleStatistics()
{}

ChrSampleStatistics::ChrSampleStatistics(const ChromosomeId& chromosome_id,
                                         const GenomicRegion::Length& size):
    chr_id(chromosome_id), coverage(size)
{}

void ChrSampleStatistics::increase_coverage(const ChrPosition& begin_pos, const size_t& read_size)
{
    if (begin_pos+read_size > coverage.size()) {
        using namespace Races::Mutants;

        throw std::runtime_error("The chromosome " + GenomicPosition::chrtos(chr_id) + " has length "
                                 + std::to_string(coverage.size())+". Coverage cannot be increased "
                                 + "up to position " + std::to_string(begin_pos) + "+"
                                 + std::to_string(read_size) + " as requested.");
    }

    for (auto base=coverage.begin()+begin_pos; base != coverage.begin()+begin_pos+read_size; ++base) {
        ++(*base);
    }
}

void check_in(const SNV& snv, const ChromosomeId& chr_id)
{
    if (snv.chr_id != chr_id) {
        using namespace Races::Mutants;

        std::ostringstream oss;

        oss << snv << " does not lay in the chromosome " 
            << GenomicPosition::chrtos(chr_id) << ".";

        throw std::domain_error(oss.str());
    }
}

void check_in(const GenomicPosition& pos, const ChromosomeId& chr_id)
{
    if (pos.chr_id != chr_id) {
        using namespace Races::Mutants;

        std::ostringstream oss;

        oss << pos << " does not belong to the chromosome " 
            << GenomicPosition::chrtos(chr_id) << ".";

        throw std::domain_error(oss.str());
    }
}

void ChrSampleStatistics::add_occurrence(const SNV& snv)
{
    check_in(snv, chr_id);

    auto SNV_occ_it = SNV_occurrences.find(snv);

    if (SNV_occ_it == SNV_occurrences.end()) {
        SNV_occurrences.insert({snv, 1});
    } else {
        ++(SNV_occ_it->second);
    }
}

const uint16_t& ChrSampleStatistics::get_coverage(const GenomicPosition& position) const
{
    check_in(position, chr_id);
        
    return coverage[position.position];
}

size_t ChrSampleStatistics::number_of_occurrences(const SNV& snv) const
{
    check_in(snv, chr_id);

    auto SNV_occ_it = SNV_occurrences.find(snv);

    if (SNV_occ_it != SNV_occurrences.end()) {
        return SNV_occ_it->second;
    }
    
    return 0;
}

void check_same_chr(const ChrSampleStatistics& a, const ChrSampleStatistics& b)
{
    if (a.get_chr_id() != b.get_chr_id()) {
        using namespace Races::Mutants;

        std::ostringstream oss;

        oss << "The two sequencing chromosome statistics refer "
            << "to different chromosomes: " << GenomicPosition::chrtos(a.get_chr_id())
            << " and " << GenomicPosition::chrtos(b.get_chr_id()) << ".";

        throw std::domain_error(oss.str());
    }

    if (a.get_chr_length() != b.get_chr_length()) {
        using namespace Races::Mutants;

        std::ostringstream oss;

        oss << "The two coverages have different lengths: " << a.get_chr_length()
            << " and " << b.get_chr_length() << ".";

        throw std::domain_error(oss.str());
    }
}

void update_coverage(std::vector<uint16_t>& a, const std::vector<uint16_t>& b)
{
    auto b_it=b.begin();
    for (auto a_it=a.begin(); a_it != a.end(); ++a_it, ++b_it) {
        *a_it += *b_it;
    }
}

std::map<SNV, size_t>::iterator
find_not_before(const SNV& snv, const std::map<SNV, size_t>& occurrences,
                std::map<SNV, size_t>::iterator it)
{
    std::less<SNV> before;
    while (it != occurrences.end() && before(it->first, snv)) {
        ++it;
    }

    return it;
}

ChrSampleStatistics& ChrSampleStatistics::operator+=(const ChrSampleStatistics& chr_stats)
{
    check_same_chr(*this, chr_stats);

    update_coverage(coverage, chr_stats.coverage);

    auto snv_it = SNV_occurrences.begin();
    auto chr_stats_snv_it = chr_stats.SNV_occurrences.begin();

    std::less<SNV> before;
    while (chr_stats_snv_it != chr_stats.SNV_occurrences.end()) {

        // find occurrence of chr_stats_snv_it->first in SNV_occurrences
        snv_it = find_not_before(chr_stats_snv_it->first, SNV_occurrences, snv_it);

        if (snv_it != SNV_occurrences.end()) {
            // found an SNV not laying before chr_stats_snv_it->first

            if (!before(chr_stats_snv_it->first, snv_it->first)) {
                // if the found SNV does not lay after chr_stats_snv_it->first,
                // then it is in the same SNV

                snv_it->second += chr_stats_snv_it->second;
            }
            ++chr_stats_snv_it;
        } else {
            for (; chr_stats_snv_it != chr_stats.SNV_occurrences.end(); ++chr_stats_snv_it) {
                SNV_occurrences.insert(*chr_stats_snv_it);
            }
        }
    }

    return *this;
}

ChrSampleStatistics& ChrSampleStatistics::operator+=(ChrSampleStatistics&& chr_stats)
{
    return this->operator+=(chr_stats);
}

ChrSampleStatistics operator+(const ChrSampleStatistics& a, const ChrSampleStatistics& b)
{
    ChrSampleStatistics result(a);

    result += b;

    return result;
}

ChrSampleStatistics& SampleStatistics::operator[](const ChromosomeId& chr_id)
{
    auto found = chr_stats.find(chr_id);

    if (found == chr_stats.end()) {
        using namespace Races::Mutants;

        throw std::runtime_error("Statistics does not include chromosome " 
                                 + GenomicPosition::chrtos(chr_id) + ".");
    }
    return found->second;
}

const ChrSampleStatistics& SampleStatistics::operator[](const ChromosomeId& chr_id) const
{
    auto found = chr_stats.find(chr_id);

    if (found == chr_stats.end()) {
        using namespace Races::Mutants;

        throw std::runtime_error("Statistics does not include chromosome " 
                                 + GenomicPosition::chrtos(chr_id) + ".");
    }
    return found->second;
}

void SampleStatistics::increase_coverage(const ChromosomeId& chromosome_id, const ChrPosition& begin_pos, 
                                         const size_t& read_size)
{
    this->operator[](chromosome_id).increase_coverage(begin_pos, read_size);
}

void SampleStatistics::add_occurrence(const SNV& snv)
{
    this->operator[](snv.chr_id).add_occurrence(snv);
}

SampleStatistics Statistics::get_total() const
{
    SampleStatistics total;

    auto sample_stats_it = sample_statistics.begin();
    for (const auto& [chr_id, chr_stats]: sample_stats_it->second.get_chr_statistics()) {
        total.add_chromosome(chr_id, chr_stats.get_chr_length());
    }    

    for (const auto& [sample_name, statistics]: sample_statistics) {
        for (const auto& [chr_id, chr_stats]: statistics.get_chr_statistics()) {
            total[chr_id] += chr_stats;
        }
    }

    return total;
}

size_t SampleStatistics::number_of_occurrences(const SNV& snv) const
{
    auto it = chr_stats.find(snv.chr_id);
    if (it == chr_stats.end()) {
        throw std::runtime_error("Unknown chromosome "+std::to_string(snv.chr_id));
    }
    return it->second.number_of_occurrences(snv);
}

const uint16_t& SampleStatistics::get_coverage(const GenomicPosition& position) const
{
    auto it = chr_stats.find(position.chr_id);
    if (it == chr_stats.end()) {
        throw std::runtime_error("Unknown chromosome "+std::to_string(position.chr_id));
    }
    return it->second.get_coverage(position);
}

void print_sample_data_about(std::ofstream& os, const SampleStatistics& statistics, 
                             const SNV& snv, const char separator='\t')
{
    const auto& coverage = statistics.get_coverage(snv);

    auto occurrences = statistics.number_of_occurrences(snv);

    if (occurrences!=0) {
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
    
    SampleStatistics total = get_total();
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

    for (const auto& [chr_id, chr_stats]: total.get_chr_statistics()) {
        for (const auto& [snv, snv_occurrences]: chr_stats.get_SNV_occurrences()) {
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
}


#if WITH_MATPLOT

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
        const auto& coverage = (sample_statistics.begin()->second)[chromosome_id].get_coverage();
        
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

            const auto& coverage = statistics[chromosome_id].get_coverage();

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
        const auto& chr_stats = g_statistics[chromosome_id];
        const auto& coverage = chr_stats.get_coverage();
        
        for (const auto& [snv, occurrences]: chr_stats.get_SNV_occurrences()) {
            auto total_it = total_occurrences.find(snv);
            if (total_it == total_occurrences.end()) {
                total_occurrences[snv] = {occurrences, coverage[snv.position]};
            } else {
                total_it->second.first += occurrences;
                total_it->second.second += coverage[snv.position];
            };
        }
    }

    return total_occurrences;
}

std::vector<double> get_VAF_data(const SampleStatistics& sample_statistics, 
                                 const ChromosomeId& chromosome_id, const double& threshold=0.0)
{
    std::vector<double> VAF;

    const auto found = sample_statistics.get_chr_statistics().find(chromosome_id);
    if (found == sample_statistics.get_chr_statistics().end()) {
        return VAF;
    }

    const auto& coverage = found->second.get_coverage();
    for (const auto& [snv, occurrences]: found->second.get_SNV_occurrences()) {
        auto value = static_cast<double>(occurrences)/coverage[snv.position];

        if (value > threshold) {
            VAF.push_back(value);
        }
    }

    return  VAF;
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
        std::vector<double> VAF = get_VAF_data(statistics, chromosome_id, 0.15);

        if (VAF.size()>0) {
            hist(VAF,num_of_bins);
        }

        xlabel("VAF");
        ylabel("Num. of SNVs");
    } else {
        auto total_occurrences = get_total_occurrences(*this, chromosome_id);

        tiledlayout(++num_of_plots, 1);
        for (const auto& [sample_name, statistics]: sample_statistics) {
            nexttile();

            std::vector<double> VAF = get_VAF_data(statistics, chromosome_id, 0.15);

            if (VAF.size()>0) {
                hist(VAF,num_of_bins);
            }
            
            const auto latex_name = std::regex_replace(sample_name, std::regex("_"), "\\\\_");
            title("Sample \""+latex_name+"\"");

            ylabel("Num. of SNVs");
        }

        nexttile();

        SampleStatistics total = get_total();

        std::vector<double> VAF = get_VAF_data(total, chromosome_id, 0.15);

        if (VAF.size()>0) {
            hist(VAF,num_of_bins);
        }
        title("Overall");

        xlabel("VAF");
        ylabel("Num. of SNVs");
    }

    sgtitle("Chromosome "+GenomicPosition::chrtos(chromosome_id)+"'s SNV VAF");

    f->save(filename);
}

#endif // WITH_MATPLOT

} // SequencingSimulator

} // Mutations

} // Races
