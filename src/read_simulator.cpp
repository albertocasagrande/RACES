/**
 * @file read_simulator.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements classes to simulate sequencing
 * @version 0.14
 * @date 2024-01-09
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
    if (begin_pos+read_size-1 > coverage.size()) {
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

void update_occurrences(std::map<SNV, size_t>& a, const std::map<SNV, size_t>& b)
{
    auto a_it = a.begin();
    auto b_it = b.begin();

    std::less<SNV> before;
    while (b_it != b.end()) {

        // find occurrence of b_it->first in a
        a_it = find_not_before(b_it->first, a, a_it);

        if (a_it != a.end()) {
            // found an SNV not laying before b_it->first

            if (!before(b_it->first, a_it->first)) {
                // if the found SNV does not lay after b_it->first,
                // then it is in the same SNV

                a_it->second += b_it->second;
            }
            ++b_it;
        } else {
            for (; b_it != b.end(); ++b_it) {
                a.insert(*b_it);
            }
        }
    }
}

ChrSampleStatistics& ChrSampleStatistics::operator+=(const ChrSampleStatistics& chr_stats)
{
    if (coverage.size()==0) {
        if (chr_stats.coverage.size()!=0) {
            *this = chr_stats;
        }

        return *this;
    }
    
    if (chr_stats.coverage.size()==0) {
        return *this;
    }

    check_same_chr(*this, chr_stats);

    update_coverage(coverage, chr_stats.coverage);

    update_occurrences(SNV_occurrences, chr_stats.SNV_occurrences);

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

std::filesystem::path SampleStatistics::get_chr_filename(const ChromosomeId& chr_id) const
{
    return get_data_directory()/("chr_"+GenomicPosition::chrtos(chr_id)+".chr_stats");
}

SampleStatistics::SampleStatistics(const std::filesystem::path& data_directory,
                                   const std::string& sample_name):
    data_dir(data_directory/sample_name), sample_name(sample_name)
{}

void create_dir(const std::filesystem::path& directory)
{
    using namespace std::filesystem;

    if (!exists(directory)) {
        create_directory(directory);

        return;
    }

    if (!is_directory(directory)) {
        throw std::runtime_error(std::string(directory)
                                 + " exists, but it is not a directory.");
    }
}

void SampleStatistics::add_chr_statistics(const ChrSampleStatistics& chr_stats)
{
    const auto& chr_id = chr_stats.get_chr_id();

    if (chr_ids.count(chr_id)>0) {
        using namespace Races::Mutants;

        throw std::runtime_error("SampleStatistics " + sample_name + " already "
                                 + "includes statistics about chromosome " 
                                 + GenomicPosition::chrtos(chr_id) + ".");
    }

    chr_ids.insert(chr_id);
    create_dir(get_data_directory());

    const auto chr_filename = get_chr_filename(chr_id);
    Archive::Binary::Out out_archive(chr_filename);
    chr_stats.save(out_archive);
}

ChrSampleStatistics SampleStatistics::get_chr_statistics(const ChromosomeId& chr_id) const
{
    if (chr_ids.count(chr_id)==0) {
        using namespace Races::Mutants;

        throw std::runtime_error("SampleStatistics " + sample_name + " does not "
                                 + "include statistics about chromosome " 
                                 + GenomicPosition::chrtos(chr_id) + ".");
    }

    const auto chr_filename = get_chr_filename(chr_id);
    Archive::Binary::In in_archive(chr_filename);
    return ChrSampleStatistics::load(in_archive);
}

void print_SNV_data(std::ofstream& os, const size_t coverage, 
                    const size_t occurrences, const char separator='\t')
{
    if (occurrences!=0) {
        os << separator << occurrences
           << separator << coverage 
           << separator << (static_cast<double>(occurrences)/coverage);
    } else {
        os << separator << "0" << separator << coverage << separator << "0";
    }
}

std::ofstream& Statistics::stream_VAF_csv_header(std::ofstream& os, const char& separator) const
{
    os << "chr" << separator << "from" << separator << "to" << separator
        << "ref" << separator << "alt" << separator
        << "type" << separator << "context" << separator << "cause";

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

    return os;
}

void check_in(const std::map<std::string, SampleStatistics>& sample_statistics,
              const std::string& sample_name)
{
    if (sample_statistics.count(sample_name)>0) {
        throw std::runtime_error("The statistics of the sample \"" + sample_name 
                                 + "\" are already contained by the object.");
    }   
}

Statistics::Statistics(const std::filesystem::path& data_directory):
    data_dir(data_directory)
{
    create_dir(data_dir);
}

SampleStatistics& Statistics::add_chr_statistics(const std::string& sample_name,
                                                 const ChrSampleStatistics& chr_statistics)
{
    auto found = sample_statistics.find(sample_name);

    if (found == sample_statistics.end()) {
        found = sample_statistics.insert({sample_name, SampleStatistics(data_dir, sample_name)}).first;
    }

    found->second.add_chr_statistics(chr_statistics);

    return found->second;
}

const SampleStatistics& Statistics::operator[](const std::string& sample_name) const
{
    auto found = find(sample_name);
    if (found == end()) {
        throw std::runtime_error("This object does not contain data about of "
                                 "the sample \"" + sample_name + "\".");
    }

    return found->second;
}


void Statistics::save_VAF_CSVs(const std::string& base_name) const
{
    for (const auto& chr_id : repr_chr_ids) {
        save_VAF_CSV(base_name+GenomicPosition::chrtos(chr_id)+".csv", chr_id);
    }
}

void Statistics::save_VAF_CSV(const std::filesystem::path& filename, 
                              const ChromosomeId& chr_id) const
{
    const char separator='\t';

    std::ofstream os(data_dir/filename);

    stream_VAF_csv_header(os, separator);

    if (sample_statistics.size()==0) {
        return;
    }

    std::map<std::string, ChrSampleStatistics> chr_statistics;
    std::map<SNV, size_t> total_SNV_occurrences;
    for (const auto& [sample_name, sample_stats]: sample_statistics) {
        auto chr_stats = sample_stats.get_chr_statistics(chr_id);
        update_occurrences(total_SNV_occurrences, chr_stats.get_SNV_occurrences());

        chr_statistics.insert({sample_name, std::move(chr_stats)});
    }

    for (const auto& [snv, total_snv_occurrences]: total_SNV_occurrences) {
        os << GenomicPosition::chrtos(chr_id) << separator << snv.position 
        << separator << snv.position 
        << separator << snv.context.get_central_nucleotide()
        << separator << snv.mutated_base << separator << "SNV" 
        << separator << snv.context << separator << snv.cause;

        size_t total_snv_coverage{0};
        for (const auto& [sample_name, chr_stats]: chr_statistics) {
            const auto& snv_coverage = chr_stats.get_coverage(snv);
            const auto& snv_occurrences = chr_stats.number_of_occurrences(snv);

            print_SNV_data(os, snv_coverage, snv_occurrences, separator);

            total_snv_coverage += snv_coverage;
        }

        if (sample_statistics.size()>1) {
            print_SNV_data(os, total_snv_coverage, total_snv_occurrences, separator);
        }

        os << std::endl;
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

void Statistics::save_coverage_images(const std::string& base_name) const
{
    for (const auto& chr_id : repr_chr_ids) {
        save_coverage_image(base_name + GenomicPosition::chrtos(chr_id) + "_coverage.jpg", 
                            chr_id);
    }
}

void Statistics::save_coverage_image(const std::filesystem::path& filename, 
                                     const ChromosomeId& chromosome_id) const
{
    using namespace matplot;

    size_t num_of_plots = sample_statistics.size();
    const size_t new_size = 300;

    auto f = figure(true);
    f->size(1500, 1500);

    if (num_of_plots==1) {
        const auto& first_sample_stats = sample_statistics.begin()->second;
        const auto chr_stats = first_sample_stats.get_chr_statistics(chromosome_id);

        const auto& coverage = chr_stats.get_coverage();
        
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
        for (const auto& [sample_name, sample_stats]: sample_statistics) {
            nexttile();

            const auto chr_stats = sample_stats.get_chr_statistics(chromosome_id);

            const auto& coverage = chr_stats.get_coverage();

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

    sgtitle("Chromosome " + GenomicPosition::chrtos(chromosome_id) + "'s coverage");

    f->save(data_dir/filename);
}

std::vector<double> get_VAF_data(const ChrSampleStatistics& chr_statistics, 
                                 const double& threshold=0.0)
{
    std::vector<double> VAF;

    const auto& coverage = chr_statistics.get_coverage();
    for (const auto& [snv, occurrences]: chr_statistics.get_SNV_occurrences()) {
        auto value = static_cast<double>(occurrences)/coverage[snv.position];

        if (value > threshold) {
            VAF.push_back(value);
        }
    }

    return  VAF;
}

void Statistics::save_SNV_histograms(const std::string& base_name) const
{
    for (const auto& chr_id : repr_chr_ids) {
        save_SNV_histogram(base_name + GenomicPosition::chrtos(chr_id) + "_hist.jpg", chr_id);
    }
}

void Statistics::save_SNV_histogram(const std::filesystem::path& filename,
                                    const ChromosomeId& chromosome_id) const
{
    using namespace matplot;

    size_t num_of_bins = 50;
    size_t num_of_plots = sample_statistics.size();

    auto f = figure(true);
    f->size(1500, 1500);
    
    if (num_of_plots==1) {
        const auto& sample_stats = (sample_statistics.begin()->second);

        const auto chr_stats = sample_stats.get_chr_statistics(chromosome_id);
        std::vector<double> VAF = get_VAF_data(chr_stats, 0.15);

        if (VAF.size()>0) {
            hist(VAF,num_of_bins);
        }

        xlabel("VAF");
        ylabel("Num. of SNVs");
    } else {
        ChrSampleStatistics total_stats;

        tiledlayout(++num_of_plots, 1);
        for (const auto& [sample_name, sample_stats]: sample_statistics) {
            nexttile();
            const auto chr_stats = sample_stats.get_chr_statistics(chromosome_id);

            total_stats += chr_stats;

            std::vector<double> VAF = get_VAF_data(chr_stats, 0.15);

            if (VAF.size()>0) {
                hist(VAF,num_of_bins);
            }
            
            const auto latex_name = std::regex_replace(sample_name, std::regex("_"), "\\\\_");
            title("Sample \""+latex_name+"\"");

            ylabel("Num. of SNVs");
        }

        nexttile();

        std::vector<double> VAF = get_VAF_data(total_stats, 0.15);

        if (VAF.size()>0) {
            hist(VAF,num_of_bins);
        }
        title("Overall");

        xlabel("VAF");
        ylabel("Num. of SNVs");
    }

    sgtitle("Chromosome " + GenomicPosition::chrtos(chromosome_id)+"'s SNV VAF");

    f->save(data_dir/filename);
}

#endif // WITH_MATPLOT

} // SequencingSimulator

} // Mutations

} // Races
