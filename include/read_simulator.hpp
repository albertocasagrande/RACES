/**
 * @file read_simulator.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines classes to simulate sequencing
 * @version 0.37
 * @date 2024-04-01
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

#ifndef __RACES_READ_SIMULATOR__
#define __RACES_READ_SIMULATOR__

#include <list>
#include <set>
#include <map>
#include <string>
#include <sstream>
#include <filesystem>
#include <random>

#include "snv.hpp"
#include "genome_mutations.hpp"

#include "fasta_chr_reader.hpp"

#include "progress_bar.hpp"
#include "variables.hpp"

#include "utils.hpp"

#include "sequencer.hpp"

#if WITH_MATPLOT
#include <matplot/matplot.h>
#endif // WITH_MATPLOT

namespace Races
{

namespace Mutations
{

/**
 * @brief The sequencing simulation namespace
 */

namespace SequencingSimulations
{

using BaseCoverage = uint16_t;

/**
 * @brief Simulated sequencing chromosome statistics
 */
class ChrCoverage
{
    ChromosomeId chr_id;                //!< The chromosome id
    std::vector<BaseCoverage> coverage; //!< The chromosome base coverage
public:

    /**
     * @brief The empty constructor
     */
    ChrCoverage();

    /**
     * @brief A constructor
     *
     * @param chromosome_id is the identifier of the chromosome whose coverage refer to
     * @param size is the size of the chromosome whose coverage refer to
     */
    ChrCoverage(const ChromosomeId& chromosome_id, const GenomicRegion::Length& size);

    /**
     * @brief Increase the coverage of a region
     *
     * @param begin_pos is the initial position of the region whose coverage must be increase
     * @param read_size is the size of the region whose coverage must be increase
     */
    void increase_coverage(const ChrPosition& begin_pos, const size_t& read_size);

    /**
     * @brief Get the identifier of the chromosome whose coverage refer to
     *
     * @return a constant reference to the identifier of the chromosome whose
     *          coverage refer to
     */
    inline const ChromosomeId& get_chr_id() const
    {
        return chr_id;
    }

    /**
     * @brief Get the length of the chromosome whose coverage refer to
     *
     * @return the size of the chromosome whose coverage refer to
     */
    inline size_t size() const
    {
        return coverage.size();
    }

    /**
     * @brief Check whether the chromosome coverage is initialized
     *
     * @return `true` if and only if this object was created by using
     *      a non-empty constructor or was copied by an initialized
     *      object
     */
    inline bool is_initialized() const
    {
        return coverage.size()!=0;
    }

    /**
     * @brief Get the coverage of a genomic position
     *
     * @param position is the genomic position whose coverage is aimed
     * @return a constant reference to the coverage of the aimed position
     */
    const BaseCoverage& get_coverage(const GenomicPosition& position) const;

    /**
     * @brief Get the chromosome coverage
     *
     * @return a constant reference to the chromosome coverage vector
     */
    inline const std::vector<BaseCoverage>& get_coverage_vector() const
    {
        return coverage;
    }

    /**
     * @brief Join a chromosome coverage to the current one
     *
     * This method updates the coverage of the current object by
     * adding the coverage of the parameter. Moreover, it adds to
     * the current object SNV occurrences the occurrences of the
     * parameter.
     *
     * @param chr_coverage is the sequencing chromosome coverage to join
     * @return a reference to the updated coverage
     */
    ChrCoverage& operator+=(ChrCoverage&& chr_coverage);

    /**
     * @brief Join a chromosome coverage to the current one
     *
     * This method updates the coverage of the current object by
     * adding the coverage of the parameter. Moreover, it adds to
     * the current object SNV occurrences the occurrences of the
     * parameter.
     *
     * @param chr_coverage is the sequencing chromosome coverage to join
     * @return a reference to the updated coverage
     */
    ChrCoverage& operator+=(const ChrCoverage& chr_coverage);

    /**
     * @brief Save sequencing chromosome coverage in an archive
     *
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        archive & chr_id
                & coverage;
    }

    /**
     * @brief Load sequencing chromosome coverage from an archive
     *
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the sequencing chromosome coverage
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    inline static ChrCoverage load(ARCHIVE& archive)
    {
        ChrCoverage chr_coverage;

        archive & chr_coverage.chr_id
                & chr_coverage.coverage;

        return chr_coverage;
    }
};

/**
 * @brief Sequencing data of any SNV
 */
struct SNVData
{
    BaseCoverage num_of_occurrences;        //!< The number of occurrences
    std::set<std::string> causes;           //!< The causes of the SNV
    std::set<Mutation::Nature> nature_set;  //!< The nature of the SNV

    /**
     * @brief The empty contructor
     */
    SNVData();

    /**
     * @brief Construct a new SNVData object
     *
     * @param snv is the snv whose data refer to
     * @param num_of_occurrences is the number of occurrences of `snv`
     */
    SNVData(const SNV& snv, const BaseCoverage num_of_occurrences);

    /**
     * @brief Update the SNV data
     *
     * This method increases the number of occurrences of the
     * SNV and adds the cause and type of a SNV to the
     * causes and types stored in the SNV data.
     *
     * @param snv is the SNV whose cause and type should be
     *          added to the causes and types stored in the
     *          data
     */
    void account_for(const SNV& snv);

    /**
     * @brief Update the current object
     *
     * This method adds all the data contained in an SNVData
     * object to the current object. It sums the number of
     * occurrences of the parameter to those of the current
     * object. Moreover, it also adds the causes and the
     * types of the parameter to the current object causes
     * and types.
     *
     * @param data is the SNVData object whose data is added
     *          to the current object
     */
    void update(const SNVData& data);
};

/**
 * @brief Simulated sequencing chromosome statistics
 */
class ChrSampleStatistics : public ChrCoverage
{
    std::map<SNV, SNVData> SNV_data;    //!< The SNV data about the simulated reads
public:

    /**
     * @brief The empty constructor
     */
    ChrSampleStatistics();

    /**
     * @brief A constructor
     *
     * @param chromosome_id is the identifier of the chromosome whose sequencing
     *          statistics are collected
     * @param size is the size of the chromosome whose sequencing
     *          statistics are collected
     */
    ChrSampleStatistics(const ChromosomeId& chromosome_id, const GenomicRegion::Length& size);

    /**
     * @brief Update the data associated to an SNV
     *
     * @param snv is the SNV whose data must be updated
     */
    void account_for(const SNV& snv);

    /**
     * @brief Get the collected SNV data
     *
     * @return a constant reference to the collected SNV data
     */
    inline const std::map<SNV, SNVData>& get_SNV_data() const
    {
        return SNV_data;
    }

    /**
     * @brief Get the collected data of an SNV
     *
     * @param snv is an SNV
     * @return the data of `snv`
     */
    SNVData get_SNV_data(const SNV& snv) const;

    /**
     * @brief Join other chromosome statistics to the current one
     *
     * This method updates the coverage of the current object by
     * adding the coverage of the parameter. Moreover, it adds to
     * the current object SNV occurrences the occurrences of the
     * parameter.
     *
     * @param chr_stats is the sequencing chromosome statistics to join
     * @return a reference to the updated objects
     */
    ChrSampleStatistics& operator+=(ChrSampleStatistics&& chr_stats);

    /**
     * @brief Join other chromosome statistics to the current one
     *
     * This method updates the coverage of the current object by
     * adding the coverage of the parameter. Moreover, it adds to
     * the current object SNV occurrences the occurrences of the
     * parameter.
     *
     * @param chr_stats is the sequencing chromosome statistics to join
     * @return a reference to the updated objects
     */
    ChrSampleStatistics& operator+=(const ChrSampleStatistics& chr_stats);

    /**
     * @brief Canonize the SNV data domains of a chromosome sample statistics list
     *
     * This method takes a list of chromosome sample statistics and  alter them
     * so that the domain of map returned by the method `get_SNV_data()` is the
     * union of the domains of the maps returned by the same method on each of
     * them.
     *
     * @param[in, out] chr_stats_list is a list of chromosome sample statistics
     */
    static void canonize_SNV_data(std::list<ChrSampleStatistics>& chr_stats_list);
};

/**
 * @brief Statistics about simulated sequencing on a set of tissue samples
 */
class SampleSetStatistics;

/**
 * @brief Statistics about simulated sequencing of a tissue sample
 */
class SampleStatistics
{
    std::filesystem::path data_dir;     //!< The directory in which the data is saved
    std::string sample_name;            //!< The name of the sample

    std::set<ChromosomeId> chr_ids;     //!< The identifiers of the chromosomes included in the statistics

    std::map<SNV, SNVData> SNV_data;    //!< The SNV data in the sample
    std::map<GenomicPosition, BaseCoverage> locus_coverage;   //!< The coverage of any position hosting an SNV

    bool save_coverage; //!< A flag to enable/disable storage of coverage data

    /**
     * @brief Get the filename of the chromosome coverage file
     *
     * @param chr_id is the identifier of the chromosome whose coverage is
     *          going to be saved/load
     * @return the filename of the chromosome coverage file
     */
    std::filesystem::path get_coverage_filename(const ChromosomeId& chr_id) const;
public:

    /**
     * @brief A constructor
     *
     * @param data_directory is the path of the directory containing the data files
     * @param sample_name is the name of the sample whose statistics refer to
     * @param save_coverage is a flag to enable/disable storage of coverage data
     */
    SampleStatistics(const std::filesystem::path& data_directory, const std::string& sample_name,
                     const bool save_coverage=false);

    /**
     * @brief Add a chromosome data to the statistics
     *
     * @param chromosome_id is the id of the chromosome to be added
     * @param size is the size of the chromosome
     */
    inline void add_chr_statistics(ChrSampleStatistics&& chr_stats)
    {
        add_chr_statistics(chr_stats);
    }

    /**
     * @brief Add a chromosome data to the statistics
     *
     * @param chromosome_id is the id of the chromosome to be added
     * @param size is the size of the chromosome
     */
    void add_chr_statistics(const ChrSampleStatistics& chr_stats);

    /**
     * @brief Get the chromosome coverage by its identifier
     *
     * @param chr_id is the id of the chromosome whose sequencing coverage is aimed
     * @return the sequencing coverage of the chromosome whose identifier is `chr_id`
     */
    ChrCoverage get_chr_coverage(const ChromosomeId& chr_id) const;

    /**
     * @brief Get the chromosome statistics
     *
     * @return a constant reference to the chromosome statistics
     */
    inline const std::set<ChromosomeId>& get_chromosome_ids() const
    {
        return chr_ids;
    }

    /**
     * @brief Get the sample name
     *
     * @return the name of the sample whose data refer to
     */
    inline const std::string& get_sample_name() const
    {
        return sample_name;
    }

    /**
     * @brief Get the SNV data
     *
     * @return a constant reference to the SNV data
     */
    inline const std::map<SNV, SNVData>& get_SNV_data() const
    {
        return SNV_data;
    }

    /**
     * @brief Get the collected data of an SNV
     *
     * @param snv is an SNV
     * @return the data of `snv`
     */
    SNVData get_SNV_data(const SNV& snv) const;

    /**
     * @brief Get the coverage of the positions in which SNVs occurr
     *
     * @return a constant reference to the coverage of the positions in which SNVs
     *      occurr
     */
    inline const std::map<GenomicPosition, BaseCoverage>& get_SNV_coverage() const
    {
        return locus_coverage;
    }

    /**
     * @brief Get the coverage of a position in which an SNV occurs
     *
     * @return a constant reference to the coverage of `pos` assuming that
     *      `pos` is a possition in which an SNV occurs
     */
    const BaseCoverage& get_SNV_coverage(const GenomicPosition& pos) const;

    /**
     * @brief Get the data directory path
     *
     * @return as constant reference to the path of the directory
     *          which contained the object saved data
     */
    inline const std::filesystem::path& get_data_directory() const
    {
        return data_dir;
    }

    /**
     * @brief Save sequencing sample statistics in an archive
     *
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        archive & data_dir
                & sample_name
                & save_coverage
                & chr_ids
                & SNV_data
                & locus_coverage;
    }

    /**
     * @brief Load sequencing sample statistics from an archive
     *
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the sequencing sample statistics
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    inline static SampleStatistics load(ARCHIVE& archive)
    {
        std::string data_dir, sample_name;
        bool save_coverage;

        archive & data_dir
                & sample_name
                & save_coverage;

        SampleStatistics sample_stats(data_dir, sample_name, save_coverage);

        archive & sample_stats.chr_ids
                & sample_stats.SNV_data
                & sample_stats.locus_coverage;

        return sample_stats;
    }

    friend class SampleSetStatistics;
};

/**
 * @brief A read simulator
 *
 * @tparam RANDOM_GENERATOR is the random generator type used to sample genome
 */
template<typename RANDOM_GENERATOR=std::mt19937_64>
class ReadSimulator;

class SampleSetStatistics
{
    std::map<std::string, SampleStatistics> stats_map;  //!< The map of the sample statistics

    std::set<ChromosomeId> repr_chr_ids; //!< The identifiers of the chromosomes whose statistics have been collected

    std::filesystem::path data_dir; //!< The directory in which the data is saved

    /**
     * @brief Stream the header of a VAF csv
     *
     * @param os is the output stream
     * @param separator is the column separator symbol
     * @return the updated stream
     */
    std::ofstream& stream_VAF_csv_header(std::ofstream& os, const char& separator) const;
public:
    /**
     * @brief A constructor
     *
     * @param data_directory is the path of the directory containing the data files
     */
    SampleSetStatistics(const std::filesystem::path& data_directory);

    /**
     * @brief Check whether this object contains statistics about a sample
     *
     * @param sample_name is the sample name
     * @return `true` if and only if this object contains statistics about
     *          a sample named `sample_name`
     */
    inline std::map<std::string, SampleStatistics>::const_iterator find(const std::string& sample_name) const
    {
        return stats_map.find(sample_name);
    }

    /**
     * @brief Get the first constant iterator of a map from the sample names to the corresponding statistics
     *
     * @return the first constant iterator of a map from the sample names to the corresponding statistics
     */
    inline std::map<std::string, SampleStatistics>::const_iterator begin() const
    {
        return stats_map.begin();
    }

    /**
     * @brief Get the last constant iterator of a map from the sample names to the corresponding statistics
     *
     * @return the last constant iterator of a map from the sample names to the corresponding statistics
     */
    inline std::map<std::string, SampleStatistics>::const_iterator end() const
    {
        return stats_map.end();
    }

    /**
     * @brief Get the identifiers of the represented chromosome
     *
     * @return a constant reference to the identifiers of the chromosomes whose statistics have been
     *          collected
     */
    std::set<ChromosomeId> get_represented_chromosomes() const;

    /**
     * @brief Add a new chromosome to the statistics of a sample
     *
     * @param sample_name is the sample name
     * @param chr_statistics are the statistics about a chromosome
     * @param save_coverage is a flag to enable/disable storage of coverage data
     * @return a constant reference to the current statistics of the sample whose name is
     *          `sample_name`
     */
    const SampleStatistics& add_chr_statistics(const std::string& sample_name,
                                               const ChrSampleStatistics& chr_statistics,
                                               const bool save_coverage=false);

    /**
     * @brief Add a new sample statistics
     *
     * @param sample_statistics are the sample statistics to add
     * @return a constant reference to the added sample statistics
     */
    const SampleStatistics& add_statistics(const SampleStatistics& sample_statistics);

    /**
     * @brief Get a sample statistics
     *
     * @param sample_name is the sample name
     * @return a constant reference to the sample statistics
     */
    const SampleStatistics& operator[](const std::string& sample_name) const;

    /**
     * @brief Get the data directory path
     *
     * @return as constant reference to the path of the directory
     *          which contained the object saved data
     */
    inline const std::filesystem::path& get_data_directory() const
    {
        return data_dir;
    }

    /**
     * @brief Get the list of the sample names
     *
     * @return the list of the sample names
     */
    std::list<std::string> get_sample_names() const;

    /**
     * @brief Get the number of samples in the set
     *
     * @return the number of samples in the set
     */
    inline size_t num_of_samples() const
    {
        return stats_map.size();
    }

    /**
     * @brief Check whether the statistics is canonical
     *
     * A statistics is canonical when all its samples stores the occurrences of
     * the very same SNVs.
     *
     * @return `true` if and only if current statistics is canonical
     */
    bool is_canonical() const;

    /**
     * @brief Canonize the statistics
     *
     * A statistics is canonical when all its samples stores the occurrences of
     * the very same SNVs.
     * This method alters the sample statistics of this object and canonizes it.
     */
    void canonize();

    /**
     * @brief Save CSV files reporting the SNV VAFs
     *
     * This method saves one CSV for each chromosome in the reference
     * sequence. Each CSV reports the SNV VAFs in the corresponding
     * chromosome.
     *
     * @param base_name is the prefix of the name of the CSV files.
     *          The filename format is "<base_name><chromosome_name>.csv"
     * @param progress_bar_stream is the output stream for the progress bar
     * @param quiet is a Boolean flag to enable/disable a progress bar
     */
    void save_VAF_CSVs(const std::string& base_name,
                       std::ostream& progress_bar_stream=std::cout,
                       const bool& quiet = false) const;


    /**
     * @brief Save CSV files reporting the SNV VAFs
     *
     * This method saves one CSV for each chromosome in the reference
     * sequence. Each CSV reports the SNV VAFs in the corresponding
     * chromosome. The file names will have the format
     * "chr_<chromosome_name>.csv".
     *
     * @param base_name is the prefix of the name of the CSV files.
     *          The file name format is "<base_name><chromosome_name>.csv"
     * @param progress_bar_stream is the output stream for the progress bar
     * @param quiet is a Boolean flag to enable/disable a progress bar
     */
    inline void save_VAF_CSVs(std::ostream& progress_bar_stream=std::cout,
                              const bool& quiet = false) const
    {
        save_VAF_CSVs("chr_", progress_bar_stream, quiet);
    }

    /**
     * @brief Save a CSV file reporting the SNV VAFs in a chromosome
     *
     * @param filename is the image filename
     * @param chromosome_id is the identifier of the chromosome whose
     *          coverage must be represented
     */
    void save_VAF_CSV(const std::filesystem::path& filename,
                      const ChromosomeId& chromosome_id) const;

#if WITH_MATPLOT
    /**
     * @brief Save images representing the chromosome coverages
     *
     * This method saves one image for each chromosome in the reference
     * sequence. Each image depicts the coverage of the corresponding
     * chromosome.
     *
     * @param base_name is the prefix of the name of the CSV files.
     *          The filename format is "<base_name><chromosome_name>_coverage.jpg"
     * @param progress_bar_stream is the output stream for the progress bar
     * @param quiet is a Boolean flag to enable/disable a progress bar
     */
    void save_coverage_images(const std::string& base_name,
                              std::ostream& progress_bar_stream=std::cout,
                              const bool& quiet = false) const;

    /**
     * @brief Save images representing the chromosome coverages
     *
     * This method saves one image for each chromosome in the reference
     * sequence. Each image depicts the coverage of the corresponding
     * chromosome. The file names will have the format
     * "chr_<chromosome_name>_coverage.jpg".
     *
     * @param progress_bar_stream is the output stream for the progress bar
     * @param quiet is a Boolean flag to enable/disable a progress bar
     */
    inline void save_coverage_images(std::ostream& progress_bar_stream=std::cout,
                                     const bool& quiet = false) const
    {
        save_coverage_images("chr_", progress_bar_stream, quiet);
    }

    /**
     * @brief Save a image representing the coverage of a chromosome
     *
     * @param filename is the image filename
     * @param chromosome_id is the identifier of the chromosome whose
     *          coverage must be represented
     */
    void save_coverage_image(const std::filesystem::path& filename,
                             const ChromosomeId& chromosome_id) const;

    /**
     * @brief Save images representing the histogram of SNVs
     *
     * This method saves one image for each chromosome in the reference
     * sequence. Each image depicts the SNV histogram of the corresponding
     * chromosome.
     *
     * @param base_name is the prefix of the name of the CSV files.
     *          The filename syntax is "<base_name><chromosome_name>_hist.jpg"
     * @param progress_bar_stream is the output stream for the progress bar
     * @param quiet is a Boolean flag to enable/disable a progress bar
     */
    void save_SNV_histograms(const std::string& base_name,
                             std::ostream& progress_bar_stream=std::cout,
                             const bool& quiet = false) const;


    /**
     * @brief Save images representing the histogram of SNVs
     *
     * This method saves one image for each chromosome in the reference
     * sequence. Each image depicts the SNV histogram of the corresponding
     * chromosome. The file names will have the format
     * "chr_<chromosome_name>_hist.jpg".
     *
     * @param progress_bar_stream is the output stream for the progress bar
     * @param quiet is a Boolean flag to enable/disable a progress bar
     */
    inline void save_SNV_histograms(std::ostream& progress_bar_stream=std::cout,
                                    const bool& quiet = false) const
    {
        save_SNV_histograms("chr_", progress_bar_stream, quiet);
    }

    /**
     * @brief Save a image representing the histogram of a chromosome SNVs
     *
     * @param filename is the image filename
     * @param chromosome_id is the identifier of the chromosome whose
     *          SNVs must be represented
     */
    void save_SNV_histogram(const std::filesystem::path& filename,
                            const ChromosomeId& chromosome_id) const;
#endif // WITH_MATPLOT

    /**
     * @brief Save sequencing sample set statistics in an archive
     *
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        archive & data_dir
                & stats_map;
    }

    /**
     * @brief Load sequencing sample set statistics from an archive
     *
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the sequencing sample set statistics
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    inline static SampleSetStatistics load(ARCHIVE& archive)
    {
        std::string data_dir;

        archive & data_dir;

        SampleSetStatistics statistics(data_dir);

        archive & statistics.stats_map;

        for (const auto& [sample_name, sample_stats] : statistics.stats_map ) {
            for (const auto& chr_id : sample_stats.get_chromosome_ids()) {
                statistics.repr_chr_ids.insert(chr_id);
            }
        }

        return statistics;
    }

    template<typename RANDOM_GENERATOR>
    friend class ReadSimulator;
};

template<typename RANDOM_GENERATOR>
class ReadSimulator
{
public:
    /**
     * @brief The type of simulated reads
     */
    enum class ReadType {
        SINGLE_READ,
        PAIRED_READ
    };

    /**
     * @brief SAM generator mode
     */
    enum class Mode {
        CREATE,     // Create
        OVERWRITE,  // Overwrite files
        APPEND      // Append files
    };

    ReadType read_type;     //!< The type of the produced reads

    size_t read_size;       //!< The produced-read size
    size_t insert_size;     //!< The paired-read insert size

    bool write_SAM;         //!< A Boolean flag to write SAM files

private:
    bool save_coverage;     //!< A flag to enable/disable storage of coverage data

    RANDOM_GENERATOR random_generator;  //!< The random generator

    std::filesystem::path output_directory;         //!< The output directory
    std::filesystem::path ref_genome_filename;      //!< The reference genome FASTA filename

    size_t num_of_reads;    //!< Number of already placed reads

    size_t hamming_distance_threshold;    //!< Hamming distance threshold for reads

    /**
     * @brief A structure to maintain read simulation data
     */
    struct ReadSimulationData
    {
        size_t not_covered_allelic_size;   //!< Allelic size to cover yet
        size_t missing_templates;          //!< Number of templates to be placed

        /**
         * @brief The empty constructor
         */
        ReadSimulationData():
            ReadSimulationData(0,0)
        {}

        /**
         * @brief A constructor
         *
         * @param not_covered_allelic_size is the allelic size to cover yet
         * @param missing_templates is the number of templates to be placed
         */
        ReadSimulationData(const size_t& not_covered_allelic_size,
                           const size_t& missing_templates):
            not_covered_allelic_size(not_covered_allelic_size),
            missing_templates(missing_templates)
        {}
    };

    /**
     * @brief Get the overall size of a list of sample mutations
     *
     * This method computes the sum of the sizes of the genomic
     * mutations in a list of sample mutations.
     *
     * @param mutations_list is a list of sample mutations
     * @param progress_bar is a progress bar
     * @return the sum of the sizes of the genomic mutations in `mutations_list`
     */
    static size_t get_overall_size(const std::list<SampleGenomeMutations>& mutations_list,
                                   UI::ProgressBar& progress_bar)
    {
        progress_bar.set_message("Evaluating sample overall genomic sizes");
        size_t overall_size{0};
        for (const auto& sample_mutations: mutations_list) {
            for (const auto& mutations: sample_mutations.mutations) {
                overall_size += mutations->size();
                progress_bar.update_elapsed_time();
            }
        }

        return overall_size;
    }

    /**
     * @brief Check whether all the genomic mutations in a list represent the same genome
     *
     * This method checks whether all the genomic mutations in a list represent genomes having
     * the same number of chromosomes and the chromosome sizes are consistent.
     *
     * @tparam GENOME_MUTATION is the type of genome mutations
     * @param mutations is a list of genome mutations
     * @param progress_bar is a progress bar
     * @return `true` if and only if all the genomic mutations in `mutations` have the same
     *          number of chromosomes and, for all pairs of elements in `mutations`, \f$m_1\f$
     *          and \f$m_2\f$, the $i$-th chromosomes in \f$m_1\f$ and \f$m_2\f$ are equal in
     *          size
     * @todo move this method in SampleGenomicMutations
     */
    template<typename GENOME_MUTATION, std::enable_if_t<std::is_base_of_v<Races::Mutations::GenomeMutations, GENOME_MUTATION>, bool> = true>
    static bool represents_same_genome(const std::list<GENOME_MUTATION>& mutations, UI::ProgressBar& progress_bar)
    {
        if (mutations.size()==0) {
            return true;
        }

        // get the first genome mutations in the list
        const auto& first_mutations = mutations.front();

        progress_bar.set_message("Checking input");
        // for all the genome mutations in the list
        for (const auto& cell_mutations : mutations) {

            // if the first element and the current one differ in the number of chromosome
            if (cell_mutations.get_chromosomes().size()!=first_mutations.get_chromosomes().size()) {
                return false;
            }

            // for all the chromosomes in the two genome mutations
            auto cell_it = cell_mutations.get_chromosomes().begin();
            auto first_it = first_mutations.get_chromosomes().begin();
            for (;cell_it != cell_mutations.get_chromosomes().end(); ++cell_it, ++first_it) {
                // if they differ because of their chromosome identifier or because of their
                // sizes
                if ((cell_it->first != first_it->first)
                        ||(cell_it->second.size()!=first_it->second.size())) {
                    return false;
                }
            }
            progress_bar.update_elapsed_time();
        }

        // otherwise, return true
        return true;
    }

    /**
     * @brief Extract the SNVs of two SNV maps laying on a simulated read
     *
     * @param A the first position-SNV map
     * @param B the second position-SNV map
     * @param genomic_position the first genomic position of the read
     * @param read_size the size of the read
     * @return a position-SNV map containing all the SNVs in `A` and `B`
     *      that lay on the simulated read
     */
    static std::map<GenomicPosition, SNV>
    merge_SNVs_in(const std::map<GenomicPosition, SNV>& A, const std::map<GenomicPosition, SNV>& B,
                  const GenomicPosition& genomic_position, const size_t& read_size)
    {
        using namespace Races::Mutations;

        const ChrPosition last_position = genomic_position.position+read_size-1;

        std::map<GenomicPosition, SNV> merged;

        // for all the SNVs in A on the simulated read
        for (auto snv_it = A.lower_bound(genomic_position);
                snv_it != A.end() && snv_it->second.position<=last_position; ++snv_it) {
            merged.insert(*snv_it);
        }

        // for all the SNVs in B on the simulated read
        for (auto snv_it = B.lower_bound(genomic_position);
                snv_it != B.end() && snv_it->second.position<=last_position; ++snv_it) {
            merged.insert(*snv_it);
        }

        return merged;
    }

    /**
     * @brief Count the number of mismatched in a mismatch vector
     *
     * @param mismatch_vector a mismatch vector
     * @return the number of mismatches represented in the mismatch vector
     */
    static size_t get_hamming_distance(const std::vector<bool>& mismatch_vector)
    {
        size_t total_mismatches{0};
        for (const auto& mismatch: mismatch_vector) {
            if (mismatch) {
                ++total_mismatches;
            }
        }

        return total_mismatches;
    }

    /**
     * @brief Apply mutation to sequencing error vector
     *
     * @param[in,out] mismatch_vector is a Boolean vector representing mismatch with respect
     *      to the reference
     * @param[in] read_SNVs is a map from genomic positions to the corresponding SNVs of the
     *      SNV laying on the read
     * @param[in] genomic_position is the initial position of the simulated read
     * @return a Boolean vector representing the mismatch with respect to the
     *      reference genome of a read. Every `true` in the mismatch vector corresponds to a
     *      mismatch of the read with respect to the reference genome
     */
    static void apply_mutations_to(std::vector<bool>& mismatch_vector,
                                   const std::map<GenomicPosition, SNV>& read_SNVs,
                                   const GenomicPosition& genomic_position)
    {
        // for all the SNVs placed on the simulated read
        for (auto snv_it = read_SNVs.begin(); snv_it != read_SNVs.end(); ++snv_it) {

            const size_t read_pos = snv_it->first.position-genomic_position.position;
            mismatch_vector[read_pos] = true;
        }
    }

    /**
     * @brief Compute the CIGAR string of a genomic region potentially affected by SNVs
     *
     * The CIGAR string a code that represents the alignment of a sequence over a reference
     * in SAM files (see [1]). This method considers a set of SNVs, a genomic position,
     * and a size, and it produces the CIGAR code corresponding to a simulated read whose
     * sequence correspond to that of the reference genome from the specified position with
     * the exception of positions in which the given SNVs were applied.
     *
     * [1] "Sequence Alignment/Map Format Specification", The SAM/BAM Format Specification
     *     Working Group, 22 August 2022, http://samtools.github.io/hts-specs/SAMv1.pdf
     *
     * @param[in] mismatch_vector is a mismatch vector
     * @return the CIGAR code corresponding to the mismatch vector
     */
    static std::string get_CIGAR(const std::vector<bool>& mismatch_vector)
    {
        std::ostringstream oss;

        size_t mismatches=0;
        size_t matches=0;
        for (const auto& mismatch : mismatch_vector) {
            if (mismatch) {
                if (matches>0) {
                    oss << matches << "M";
                    matches = 0;
                }
                ++mismatches;
            } else {
                if (mismatches>0) {
                    oss << mismatches << "X";
                    mismatches = 0;
                }
                ++matches;
            }
        }

        if (mismatches>0) {
            oss << mismatches << "X";
        }
        if (matches>0) {
            oss << matches << "M";
        }

        return oss.str();
    }

    /**
     * @brief Apply a set of SNVs to a fragment of the reference genome
     *
     * @param nucleic_sequence is a fragment of the reference genome
     * @param SNVs is a map from genomic positions to the corresponding SNVs
     * @param genomic_position is the initial position of `nucleic_sequence` in
     *          the reference genome
     * @return a reference to the modified version of `nucleic_sequence` where
     *          the SNVs that were placed in the reference genome region
     *          corresponding to `nucleic_sequence` are applied
     */
    static std::string& apply_SNVs(std::string& nucleic_sequence,
                                   const std::map<GenomicPosition, SNV>& SNVs,
                                   const GenomicPosition& genomic_position)
    {
        using namespace Races::Mutations;

        const ChrPosition last_position = genomic_position.position+nucleic_sequence.size()-1;

        for (auto snv_it = SNVs.lower_bound(genomic_position);
                snv_it != SNVs.end() && snv_it->second.position<=last_position; ++snv_it) {
            const size_t SNV_pos = snv_it->second.position-genomic_position.position;

            {   // check the match between the SNV reference base and the corresponding base
                // in `nucleic_sequence`

                const char snv_ref = snv_it->second.ref_base;

                if (snv_ref != '?') {

                    // if the SNV reference base differs from the corresponding
                    // `nucleic_sequence` base
                    const char ref_base = toupper(nucleic_sequence[SNV_pos]);

                    if (snv_ref!=ref_base) {
                        std::ostringstream oss;

                        oss << "Reference base mismatch in chr. "
                            << GenomicPosition::chrtos(genomic_position.chr_id)
                            << " in position " << snv_it->second.position << ": expected '"
                            << snv_ref << "' and got '"
                            << ref_base << "'."<< std::endl
                            << "Are you using a FASTA file different from that used for "
                            << "building the context index?"
                            << std::endl;
                        throw std::domain_error(oss.str());
                    }
                }
            }

            // change the `nucleic_sequence` nucleotide according to the considered SNV
            nucleic_sequence[SNV_pos] = snv_it->second.alt_base;
        }

        return nucleic_sequence;
    }

    /**
     * @brief Get the template read data
     *
     * @param template_first_position is the template first position
     * @param template_size is the template size
     * @return a vector containing a pair flag-first-base position for
     *       each read in the template
     */
    std::vector<std::pair<int, ChrPosition>>
    get_template_read_data(const ChrPosition& template_first_position,
                           const size_t& template_size)
    {
        std::uniform_int_distribution<int> c_dist(0,1);

        bool rev_comp = c_dist(random_generator)==0;

        std::vector<std::pair<int, ChrPosition>> template_read_data;

        if (read_type == ReadType::PAIRED_READ) {
            int flag = 0x1|0x2;
            int first_flag = flag | 0x20 | 0x40;
            int second_flag = flag | 0x10 | 0x80;

            if (rev_comp) {
                std::swap(first_flag, second_flag);
            }
            template_read_data.push_back({first_flag,
                                          template_first_position});

            template_read_data.push_back({second_flag,
                                          template_first_position+template_size-read_size});

            return template_read_data;
        }

        int flag = 0;
        if (rev_comp) {
            flag = 0x10;
        }

        return {{flag, template_first_position}};
    }

    /**
     * @brief Collect template statistics
     *
     * @param[in] chr_data is the data about the chromosome from which the simulated read come from
     * @param[in] germline_SNVs is a map from genomic positions to the corresponding germline SNVs
     * @param[in] SNVs is a map from genomic positions to the corresponding SNVs
     * @param[in] read_first_position is the first position of the simulated read
     * @param[in] template_size is the size of the template
     * @param[in,out] chr_statistics are the chromosome statistics about a tissue sample
     */
    void collect_read_statistics(const Races::IO::FASTA::ChromosomeData<Races::IO::FASTA::Sequence>& chr_data,
                                 const std::map<GenomicPosition, SNV>& germline_SNVs,
                                 const std::map<GenomicPosition, SNV>& SNVs,
                                 const std::string& read, const ChrPosition& read_position,
                                 ChrSampleStatistics& chr_statistics)
    {
        const auto read_end = read_position+read.size();
        chr_statistics.increase_coverage(read_position, read.size());

        const GenomicPosition genomic_position{chr_data.chr_id, read_position};
        for (const auto snvs : {&SNVs, &germline_SNVs}) {
            for (auto snv_it = snvs->lower_bound(genomic_position);
                    snv_it != snvs->end() && snv_it->second.position<read_end; ++snv_it) {
                const size_t base_pos = snv_it->first.position - read_position;
                if (snv_it->second.ref_base != read[base_pos]) {
                    if (snv_it->second.alt_base != read[base_pos]) {
                        SNV snv = snv_it->second;

                        snv.alt_base = read[base_pos];
                        chr_statistics.account_for(snv);
                    } else {
                        chr_statistics.account_for(snv_it->second);
                    }
                }
            }
        }
    }

    /**
     * @brief Get the name of a template
     *
     *
     * @param template_id is the identifier of the template
     * @return the name of the template whose identifier is
     *      `template_id`.
     */
    static std::string get_template_name(const size_t& template_id)
    {
        std::ostringstream oss_name;

        oss_name << "r" <<  std::setfill('0') << std::setw(11)
                 << template_id;

        return oss_name.str();
    }

    /**
     * @brief Write the SAM alignment line corresponding to a simulated read
     *
     * This method extracts the template from the chromosome sequence, obtains the 
     * template reads, applies germinal and somatic mutations to it, simulates
     * sequencing errors, collects sequencing statistics, and, whenever it is 
     * required, writes the read alignments in the SAM stream.
     *
     * @tparam SEQUENCER is the sequencer model type
     * @param[in,out] sequencer is the sequencer
     * @param[in] chr_data is the data about the chromosome from which the simulated
     *   template come from
     * @param[in] germline_SNVs is a map from genomic positions to the corresponding
     *   germline SNVs
     * @param[in] SNVs is a map from genomic positions to the corresponding SNVs
     * @param[in] template_begin_pos is the first position of the simulated template
     * @param[in] template_size is the size of the template
     * @param[in,out] chr_statistics are the chromosome statistics about a tissue sample
     * @param[out] SAM_stream is a pointer to a SAM stream. If `SAM_stream` is
     *   different from `nullptr` and the hamming distance between the produced read 
     *   and the reference fragment is below the threshold, the read alignments is
     *   written in the SAM stream
     * @param[in] sample_name is the tissue sample name
     */
    template<typename SEQUENCER,
             std::enable_if_t<std::is_base_of_v<Races::Sequencers::BasicSequencer, SEQUENCER>, bool> = true>
    void process_template(SEQUENCER& sequencer,
                          const Races::IO::FASTA::ChromosomeData<Races::IO::FASTA::Sequence>& chr_data,
                          const std::map<GenomicPosition, SNV>& germline_SNVs,
                          const std::map<GenomicPosition, SNV>& SNVs,
                          const ChrPosition& template_begin_pos,
                          const size_t& template_size, ChrSampleStatistics& chr_statistics,
                          std::ostream* SAM_stream,
                          const std::string& sample_name="")
    {
        const auto template_read_data = get_template_read_data(template_begin_pos, template_size);

        const size_t template_id = num_of_reads/template_read_data.size();
        const std::string template_name = get_template_name(template_id);
        for (size_t i=0; i<template_read_data.size(); ++i) {
            const auto& read_first_position = template_read_data[i].second;

            const GenomicPosition genomic_position{chr_data.chr_id, read_first_position};

            std::vector<bool> mismatch_vector(read_size, false);
            std::string read = chr_data.nucleotides.substr(read_first_position-1, read_size);
            apply_SNVs(read, SNVs, genomic_position);
            apply_SNVs(read, germline_SNVs, genomic_position);

            std::string qual = sequencer.simulate_seq(read, genomic_position, i);

            for (size_t i=0; i<read_size; ++i) {
                mismatch_vector[i] = (read[i] != 'N' 
                                        && read[i] != chr_data.nucleotides[read_first_position-1+i]);
            }

            const auto mapq = 33;

            const auto hamming_distance = get_hamming_distance(mismatch_vector);

            if (hamming_distance < hamming_distance_threshold) {
                collect_read_statistics(chr_data, germline_SNVs, SNVs,
                                        read, read_first_position, chr_statistics);

                if (SAM_stream != nullptr) {
                    *SAM_stream << template_name                         // QNAME
                                << '\t' << template_read_data[i].first   // FLAG
                                << '\t' << chr_data.name                 // RNAME
                                << '\t' << read_first_position           // POS
                                << '\t' << mapq                          // MAPQ
                                << '\t' << get_CIGAR(mismatch_vector);   // CIGAR
                    if (read_type == ReadType::PAIRED_READ) {
                        const size_t paired_idx = 1-i;
                        const auto paired_pos = template_read_data[paired_idx].second;

                        *SAM_stream << '\t' << '='               // RNEXT
                                    << '\t' << paired_pos        // PNEXT
                                    << '\t' << (i==0?"":"-")
                                    << template_size;            // TLEN
                    } else {
                        *SAM_stream << '\t' << '*'   // RNEXT
                                    << '\t' << '0'   // PNEXT
                                    << '\t' << '0';  // TLEN
                    }
                    *SAM_stream << '\t' << read  // SEQ
                                << '\t' << qual;  // QUAL

                    // TAGS
                    if (sample_name.size()!=0) {
                        *SAM_stream << "\tRG:Z:" << sample_name;  // The read group
                    }
                    *SAM_stream << "\tNM:i:" << hamming_distance  // The Hamming distance
                                << std::endl;
                }
            }
        }
    }

    /**
     * @brief Generate the reads in an allelic fragment
     *
     * @tparam SEQUENCER is the sequencer model type
     * @param[in,out] sequencer is the sequencer
     * @param[in] chr_data is the data about the chromosome from which the simulated read come from
     * @param[in] fragment is the allele fragment for which reads must be generated
     * @param[in,out] sample_simulation_data are the read simulation data relative to the considered `sample`
     * @param[in] sample_name is the name of the considered tissue sample
     * @param[in,out] chr_statistics are the chromosome statistics of the tissue sample
     * @param[in] total_steps is the total number of steps required to complete the overall procedure
     * @param[in,out] steps is the number of performed steps
     * @param[in,out] progress_bar is the progress bar
     * @param[in,out] SAM_stream is the SAM file output stream
     */
    template<typename SEQUENCER,
             std::enable_if_t<std::is_base_of_v<Races::Sequencers::BasicSequencer, SEQUENCER>, bool> = true>
    void generate_fragment_reads(SEQUENCER& sequencer,
                                 const Races::IO::FASTA::ChromosomeData<Races::IO::FASTA::Sequence>& chr_data,
                                 const AlleleFragment& germline_fragment,
                                 const AlleleFragment& fragment,
                                 ReadSimulationData& sample_simulation_data, const std::string& sample_name,
                                 ChrSampleStatistics& chr_statistics,
                                 const size_t& total_steps, size_t& steps,
                                 Races::UI::ProgressBar& progress_bar, std::ostream* SAM_stream)
    {
        auto template_size = ((read_type==ReadType::PAIRED_READ)
                              ?2*read_size+insert_size:read_size);

        if (fragment.size()>=template_size) {
            const double hit_probability = static_cast<double>(fragment.size())/
                                        sample_simulation_data.not_covered_allelic_size;
            std::binomial_distribution<size_t> b_dist(sample_simulation_data.missing_templates,
                                                        hit_probability);

            size_t num_of_frag_templates = b_dist(random_generator);

            const double fragment_progress_ratio = (100*static_cast<double>(fragment.size())/num_of_frag_templates)/total_steps;

            const double current_progress = 100*static_cast<double>(steps)/total_steps;

            const auto& germline_SNVs = germline_fragment.get_SNVs();
            const auto& SNVs = fragment.get_SNVs();

            auto first_possible_begin = fragment.get_initial_position();
            auto last_possible_begin = static_cast<ChrPosition>(fragment.get_final_position())-template_size+1;

            std::uniform_int_distribution<ChrPosition> dist(first_possible_begin, last_possible_begin);
            size_t simulated_templates = 0;
            while (simulated_templates < num_of_frag_templates) {
                auto begin_pos = dist(random_generator);

                const auto template_read_data = get_template_read_data(begin_pos, template_size);

                process_template(sequencer, chr_data, germline_SNVs, SNVs, begin_pos,
                                 template_size, chr_statistics, SAM_stream, sample_name);

                num_of_reads += ((read_type == ReadType::PAIRED_READ)?2:1);
                ++simulated_templates;
                --sample_simulation_data.missing_templates;

                progress_bar.set_progress(current_progress+simulated_templates*fragment_progress_ratio);
            }
        }

        sample_simulation_data.not_covered_allelic_size -= fragment.size();
        steps += fragment.size();
    }


    /**
     * @brief Generate simulated reads on a chromosome and write their SAM alignments
     *
     * This method takes a list of genome mutations and a reference genome chromosome, it
     * generates simulated reads over it up to a specified coverage of the mutated genomes,
     * and write the corresponding SAM alignments in a stream.
     *
     * @tparam SEQUENCER is the sequencer model type
     * @param[in,out] sequencer is the sequencer
     * @param[in] chr_data is the data about the chromosome from which the simulated read come from
     * @param[in] sample_genome_mutations is a sample genome mutations
     * @param[in,out] sample_simulation_data is the read simulation data of the sample
     * @param[in] total_steps is the total number of steps required to complete the overall procedure
     * @param[in,out] steps is the number of performed steps
     * @param[in,out] progress_bar is the progress bar
     * @param[in,out] SAM_stream is the SAM file output stream
     */
    template<typename SEQUENCER,
             std::enable_if_t<std::is_base_of_v<Races::Sequencers::BasicSequencer, SEQUENCER>, bool> = true>
    ChrSampleStatistics
    generate_chromosome_reads(SEQUENCER& sequencer,
                              const Races::IO::FASTA::ChromosomeData<Races::IO::FASTA::Sequence>& chr_data,
                              const SampleGenomeMutations& sample_genome_mutations,
                              ReadSimulationData& sample_simulation_data,
                              const size_t& total_steps, size_t& steps,
                              Races::UI::ProgressBar& progress_bar,
                              std::ostream* SAM_stream=nullptr)
    {
        ChrSampleStatistics chr_stats{chr_data.chr_id,
				                      static_cast<GenomicRegion::Length>(chr_data.length)};

        for (const auto& cell_mutations: sample_genome_mutations.mutations) {
            const auto& chr_mutations = cell_mutations->get_chromosome(chr_data.chr_id);
            const auto& germline_chr_mut = sample_genome_mutations.germline_mutations.get_chromosome(chr_data.chr_id);

            for (const auto& [allele_id, allele] : chr_mutations.get_alleles()) {
                const auto& germline_allele_id = allele.get_history().front();
                const auto& germline_allele = germline_chr_mut.get_allele(germline_allele_id);

                for (const auto& [position, fragment] : allele.get_fragments()) {
                    // searching for the germline fragment AFTER `fragment`
                    auto germline_it = germline_allele.get_fragments().upper_bound(position);

                    // update `germline_mutations` so to the last germline fragment starting
                    // BEFORE or in the SAME POSITION of `fragment`
                    --germline_it;

                    // since `cell_mutations` derives from `germline_mutations`, `germline_allele`
                    // fully contains `allele`
                    generate_fragment_reads(sequencer, chr_data, germline_it->second, fragment,
                                            sample_simulation_data, sample_genome_mutations.name,
                                            chr_stats, total_steps, steps, progress_bar, SAM_stream);
                }

                progress_bar.set_progress(100*steps/total_steps);
            }
        }

        return chr_stats;
    }

    /**
     * @brief Generate simulated reads on a chromosome and write their SAM alignments
     *
     * This method takes a list of genome mutations and a reference genome chromosome, it
     * generates simulated reads over it up to a specified coverage of the mutated genomes,
     * and write the corresponding SAM alignments in a stream.
     *
     * @tparam SEQUENCER is the sequencer model type
     * @param[in,out] sequencer is the sequencer
     * @param[out] statistics are the statistics about the generated reads
     * @param[in] chr_data is the data about the chromosome from which the simulated read come from
     * @param[in] mutations_list is a list of sample mutations
     * @param[in,out] simulation_data is a list of simulated data corresponding to `mutations_list` elements
     * @param[in] total_steps is the total number of steps required to complete the overall procedure
     * @param[in,out] steps is the number of performed steps
     * @param[in,out] progress_bar is the progress bar
     * @param[in,out] SAM_stream is the SAM file output stream
     */
    template<typename SEQUENCER,
             std::enable_if_t<std::is_base_of_v<Races::Sequencers::BasicSequencer, SEQUENCER>, bool> = true>
    void generate_chromosome_reads(SEQUENCER& sequencer,
                                   SampleSetStatistics& statistics,
                                   const Races::IO::FASTA::ChromosomeData<Races::IO::FASTA::Sequence>& chr_data,
                                   const std::list<SampleGenomeMutations>& mutations_list,
                                   std::list<ReadSimulationData>& simulation_data,
                                   const size_t& total_steps, size_t& steps,
                                   Races::UI::ProgressBar& progress_bar,
                                   std::ostream* SAM_stream=nullptr)
    {
        auto chr_name = GenomicPosition::chrtos(chr_data.chr_id);
        progress_bar.set_message("Processing chr. " + chr_name);

        auto simulation_data_it = simulation_data.begin();

        std::list<ChrSampleStatistics> chr_stats_list;
        for (const auto& mutations : mutations_list) {
            chr_stats_list.push_back(generate_chromosome_reads(sequencer, chr_data, mutations,
                                                               *simulation_data_it, total_steps,
                                                               steps, progress_bar, SAM_stream));
            ++simulation_data_it;
        }

        ChrSampleStatistics::canonize_SNV_data(chr_stats_list);

        auto chr_stats_it = chr_stats_list.cbegin();
        for (const auto& mutations : mutations_list) {
            statistics.add_chr_statistics(mutations.name, *chr_stats_it, save_coverage);

            ++chr_stats_it;
        }
    }

    /**
     * @brief Get the SAM stream
     *
     * @param chr_data in the information available about the chromosome to be considered
     * @param mutations_list is a list of sample mutations
     * @param base_name is the prefix of the filename
     * @return the SAM stream
     */
    std::ofstream
    get_SAM_stream(const Races::IO::FASTA::ChromosomeData<Races::IO::FASTA::Sequence>& chr_data,
                   const std::list<SampleGenomeMutations>& mutations_list,
                   const std::string& base_name="chr_",
                   const std::string& platform="ILLUMINA") const
    {
        auto chr_str = GenomicPosition::chrtos(chr_data.chr_id);

        std::string SAM_filename = base_name + chr_str + ".sam";
        auto SAM_file_path = output_directory/SAM_filename;

        if (std::filesystem::exists(SAM_file_path)) {
            if (!std::filesystem::is_regular_file(SAM_file_path)) {
                throw std::runtime_error("\"" + to_string(SAM_file_path)
                                         + "\" is not a regular file");
            }

            return std::ofstream(SAM_file_path, std::ofstream::app);
        }

        auto SAM_stream = std::ofstream(SAM_file_path);

        // write the SAM header
        SAM_stream << "@HD\tVN:1.6\tSO:unknown" << std::endl;
        SAM_stream << "@SQ\tSN:" << chr_data.name << "\tLN:" << chr_data.length
                << "\tAN:chromosome" << chr_str << ",chr" << chr_str
                << ",chromosome_" << chr_str << ",chr_" << chr_str << std::endl;
        for (const auto& sample_mutations : mutations_list) {
            SAM_stream << "@RG\tID:" << sample_mutations.name
                       << "\tSM:" << sample_mutations.name
                       << "\tPL:" << platform << std::endl;
        }

        return SAM_stream;
    }

    /**
     * @brief Compute the initial read simulation data
     *
     * @param mutations_list is a list of sample mutations
     * @param coverage is the aimed coverage
     * @return the list of the initial sample read simulation data. The order
     *       of the returned list matches that of `mutations_list`
     */
    std::list<ReadSimulationData>
    get_initial_read_simulation_data(const std::list<SampleGenomeMutations>& mutations_list,
                                     const double& coverage) const
    {
        std::list<ReadSimulationData> read_simulation_data;

        size_t total_read_size = (read_type==ReadType::PAIRED_READ?2:1)*read_size;
        for (const auto& sample_mutations : mutations_list) {
            ReadSimulationData sample_data;
            if (sample_mutations.mutations.size()>0) {
                sample_data.missing_templates =
                                 static_cast<size_t>((sample_mutations.mutations.front()->size()*
                                                     coverage)/total_read_size);
            }
            for (const auto& cell_mutations: sample_mutations.mutations) {
                sample_data.not_covered_allelic_size += cell_mutations->allelic_size();
            }

            read_simulation_data.push_back(std::move(sample_data));
        }

        return read_simulation_data;
    }

    /**
     * @brief Get the identifiers of the chromosomes in the genome
     *
     * This method deduces the identifiers of the chromosomes in the
     * considered genome. The method assumes all the samples in the
     * mutations list to have the same chromosomes.
     *
     * @param mutations_list is a list of sample mutations
     * @return the set of identifiers of the genome chromosomes
     */
    static std::set<ChromosomeId>
    get_genome_chromosome_ids(const std::list<SampleGenomeMutations>& mutations_list)
    {
        for (const auto& sample_mutations : mutations_list) {
            for (const auto& mutations : sample_mutations.mutations) {
                std::set<ChromosomeId> ids;
                for (const auto& [chr_id, chr_mutations] : mutations->get_chromosomes()) {
                    ids.insert(chr_id);
                }

                return ids;
            }
        }

        return {};
    }

    /**
     * @brief Generate simulated reads
     *
     * This method takes a list of sample mutations and it generates simulated reads up to a
     * specified coverage of the mutated genomes. The base-coverage and SNV occurrences are
     * collected and, upon request, the SAM alignments corresponding to the produced reads
     * are saved.
     *
     * @tparam SEQUENCER is the sequencer model type
     * @param[in,out] sequencer is the sequencer
     * @param mutations_list is a list of sample mutations
     * @param coverage is the aimed coverage
     * @param base_name is the prefix of the filename
     * @param progress_bar is the progress bar
     * @return the statistics about the generated reads
     */
    template<typename SEQUENCER,
             std::enable_if_t<std::is_base_of_v<Races::Sequencers::BasicSequencer, SEQUENCER>, bool> = true>
    SampleSetStatistics generate_reads(SEQUENCER& sequencer,
                                       const std::list<SampleGenomeMutations>& mutations_list,
                                       const double& coverage, const std::string& base_name,
                                       UI::ProgressBar& progress_bar)
    {
        std::ifstream ref_stream(ref_genome_filename);

        auto read_simulation_data = get_initial_read_simulation_data(mutations_list, coverage);

        size_t steps{0};
        const auto chromosome_ids = get_genome_chromosome_ids(mutations_list);
        size_t total_steps = 2*chromosome_ids.size()*mutations_list.size();
        for (const auto& sample_data: read_simulation_data) {
            total_steps += sample_data.not_covered_allelic_size;
        }

        SampleSetStatistics statistics(output_directory);

        using namespace Races::IO::FASTA;

        ChromosomeData<Sequence> chr_data;
        progress_bar.set_progress((100*steps)/total_steps, "Reading next chromosome");
        while (ChromosomeData<Sequence>::read(ref_stream, chr_data, progress_bar)) {
            if (chromosome_ids.count(chr_data.chr_id)>0) {
                progress_bar.set_progress((100*(++steps))/total_steps);

                if (write_SAM) {
                    std::ofstream SAM_stream = get_SAM_stream(chr_data, mutations_list,
                                                              base_name,
                                                              sequencer.get_platform_name());

                    generate_chromosome_reads(sequencer, statistics, chr_data, mutations_list,
                                              read_simulation_data, total_steps, steps,
                                              progress_bar, &SAM_stream);
                } else {
                    generate_chromosome_reads(sequencer, statistics, chr_data, mutations_list,
                                              read_simulation_data, total_steps,
                                              steps, progress_bar);
                }

                statistics.repr_chr_ids.insert(chr_data.chr_id);
            }
        }

        progress_bar.set_progress(100, "Read simulated");

        return statistics;
    }

    /**
     * @brief A constructor
     *
     * @param output_directory is the output directory
     * @param ref_genome_filename is reference genome filename
     * @param read_type is the type of the produced-read, i.e., single or paired-read
     * @param read_size is the size of the output reads
     * @param insert_size is the size of the insert
     * @param mode is the SAM generator output mode
     * @param save_coverage is a flag to enable/disable storage of coverage data
     * @param seed is the random generator seed
     */
    ReadSimulator(const std::filesystem::path& output_directory,
                  const std::filesystem::path& ref_genome_filename,
                  const ReadType read_type, const size_t& read_size, const size_t& insert_size,
                  const Mode mode, const bool save_coverage, const int& seed):
        read_type(read_type), read_size(read_size), insert_size(insert_size), write_SAM(false),
        save_coverage(save_coverage), random_generator(seed), output_directory(output_directory),
        ref_genome_filename(ref_genome_filename), num_of_reads(0),
        hamming_distance_threshold(30)
    {
        namespace fs = std::filesystem;

        if (fs::exists(this->output_directory)) {
            if (mode==Mode::CREATE) {
                throw std::domain_error("\"" + to_string(output_directory)
                                        + "\" already exists");
            }

            if (!fs::is_directory(this->output_directory)) {
                throw std::domain_error("\"" + to_string(output_directory)
                                        + "\" is not a directory");
            }

            if (mode==Mode::OVERWRITE) {
                fs::remove_all(this->output_directory);
            }
        }

        fs::create_directory(this->output_directory);

        if (!fs::exists(this->ref_genome_filename)) {
            throw std::domain_error("\"" + to_string(ref_genome_filename)
                                    + "\" does not exist");
        }

        if (!fs::is_regular_file(this->ref_genome_filename)) {
            throw std::domain_error("\"" + to_string(ref_genome_filename)
                                    + "\" is not a regular file");
        }

        std::ifstream ref_stream(ref_genome_filename);
        int c = ref_stream.get();

        if (c != EOF && c != '>') {
            throw std::domain_error("\"" + to_string(ref_genome_filename)
                                    + "\" is not a FASTA file");
        }
    }
public:

    /**
     * @brief The empty constructor
     */
    ReadSimulator():
        read_type(ReadType::SINGLE_READ), read_size(0), insert_size(0), write_SAM(false),
        save_coverage(false), random_generator(), output_directory(""), ref_genome_filename(""),
        num_of_reads(0), hamming_distance_threshold(30)
    {}

    /**
     * @brief Create a simulator that produces single reads
     *
     * @param output_directory is the output directory
     * @param ref_genome_filename is reference genome filename
     * @param read_size is the size of the output reads
     * @param mode is the SAM generator output mode
     * @param save_coverage is a flag to enable/disable storage of coverage data
     * @param seed is the random generator seed
     */
    ReadSimulator(const std::filesystem::path& output_directory,
                  const std::filesystem::path& ref_genome_filename,
                  const size_t& read_size, const Mode mode=Mode::CREATE,
                  const bool& save_coverage=false, const int& seed=0):
        ReadSimulator(output_directory, ref_genome_filename, ReadType::SINGLE_READ, read_size,
                      0, mode, save_coverage, seed)
    {}

    /**
     * @brief Create a simulator that produces paired reads
     *
     * @param output_directory is the SAM output directory
     * @param ref_genome_filename is reference genome filename
     * @param read_size is the size of the output reads
     * @param insert_size is the size of the insert
     * @param mode is the SAM generator output mode
     * @param save_coverage is a flag to enable/disable storage of coverage data
     * @param seed is the random generator seed
     */
    ReadSimulator(const std::filesystem::path& output_directory,
                  const std::filesystem::path& ref_genome_filename,
                  const size_t& read_size, const size_t& insert_size, const Mode mode=Mode::CREATE,
                  const bool& save_coverage=false, const int& seed=0):
        ReadSimulator(output_directory, ref_genome_filename, ReadType::PAIRED_READ, read_size,
                      insert_size, mode, save_coverage, seed)
    {}

    /**
     * @brief Generate simulated reads for a list of sample genome mutations
     *
     * This method generates simulated reads for a list of sample genome mutations and,
     * if requested, writes the corresponding SAM alignments. The number of simulated reads
     * depends on the specified coverage.
     *
     * @tparam SEQUENCER is the sequencer model type
     * @param sequencer is the sequencer
     * @param mutations_list is a list of sample mutations
     * @param coverage is the aimed coverage
     * @param purity is ratio between the number of sampled tumoral cells, which are
     *              represented in the mutation list, and that of the overall sampled
     *              cells which contains normal cells too (default: 1.0)
     * @param base_name is the prefix of the filename (default: "chr_")
     * @param progress_bar_stream is the output stream for the progress bar
     *              (default: std::cout)
     * @param quiet is a Boolean flag to avoid progress bar and user messages
     *              (default: false)
     * @return the sample set statistics about the generated reads
     */
    template<typename SEQUENCER,
             std::enable_if_t<std::is_base_of_v<Races::Sequencers::BasicSequencer, SEQUENCER>, bool> = true>
    SampleSetStatistics operator()(SEQUENCER& sequencer,
                                   std::list<Races::Mutations::SampleGenomeMutations> mutations_list,
                                   const double& coverage, const double purity=1.0,
                                   const std::string& base_name="chr_",
                                   std::ostream& progress_bar_stream=std::cout,
                                   const bool quiet=false)
    {
        using namespace Races;

        if (coverage < 0) {
            throw std::domain_error("The coverage value must be non-negative: got "
                                    + std::to_string(coverage) + ".");
        }

        if (purity < 0 || purity > 1) {
            throw std::domain_error("The purity value must be in [0,1]: got "
                                    + std::to_string(purity) + ".");

        }

        if (read_type == ReadType::PAIRED_READ && !sequencer.supports_paired_reads()) {
            throw std::domain_error(sequencer.get_model_name()
                                    + " does not support paired reads.");
        }

        if (mutations_list.size()==0) {
            return SampleSetStatistics(output_directory);
        }

        // if purity<1, then duplicate the structure of the germline genome
        // (no mutations as they will be automatically added by `generate_reads`)
        std::shared_ptr<CellGenomeMutations> normal_structure;

        if (purity<1) {
            const auto& germline_mutations = mutations_list.front().germline_mutations;
            normal_structure = std::make_shared<CellGenomeMutations>(germline_mutations.duplicate_structure());
        }

        for (auto& mutation_list : mutations_list) {
            if (purity>0) {
                size_t tumor_number = mutation_list.mutations.size();
                size_t normal_number = static_cast<size_t>(tumor_number*(1-purity)/purity);

                for (size_t i=0; i<normal_number; ++i) {
                    mutation_list.mutations.push_back(normal_structure);
                }
            } else {
                mutation_list.mutations = {normal_structure};
            }
        }

        UI::ProgressBar progress_bar(progress_bar_stream, quiet);

        return generate_reads<SEQUENCER>(sequencer, mutations_list, coverage,
                                         base_name, progress_bar);
    }

    /**
     * @brief Generate simulated reads for a sample genome mutations
     *
     * This method generates simulated reads for a sample genome mutations and,
     * if requested, writes the corresponding SAM alignments. The number of simulated reads
     * depends on the specified coverage.
     *
     * @tparam SEQUENCER is the sequencer model type
     * @param sequencer is the sequencer
     * @param mutations is a sample genome mutations
     * @param coverage is the aimed coverage
     * @param purity is ratio between the number of sampled tumoral cells, which are
     *              represented in the mutation list, and that of the overall sampled
     *              cells which contains normal cells too (default: 1.0)
     * @param base_name is the prefix of the filename (default: "chr_")
     * @param progress_bar_stream is the output stream for the progress bar
     *              (default: std::cout)
     * @param quiet is a Boolean flag to avoid progress bar and user messages
     *              (default: false)
     * @return the statistics about the generated reads
     */
    template<typename SEQUENCER,
             std::enable_if_t<std::is_base_of_v<Races::Sequencers::BasicSequencer, SEQUENCER>, bool> = true>
    SampleStatistics operator()(SEQUENCER& sequencer,
                                const Races::Mutations::SampleGenomeMutations& mutations,
                                const double& coverage, const double purity=1.0,
                                const std::string& base_name="chr_",
                                std::ostream& progress_bar_stream=std::cout,
                                const bool quiet=false)
    {
        using namespace Races::Mutations;
        std::list<SampleGenomeMutations> mutations_list{mutations};

        auto statistics = operator()(sequencer, mutations_list, coverage, purity,
                                     base_name, progress_bar_stream, quiet);

        return statistics[mutations.name];
    }

    /**
     * @brief Enable SAM file writing
     *
     * @param enable is a flag to enable(true)/disable(false) SAM file writing
     */
    inline void enable_SAM_writing(const bool& enable)
    {
        write_SAM = enable;
    }
};

}   // SequencingSimulations

}   // Mutations

}   // Races

#endif // __RACES_READ_SIMULATOR__
