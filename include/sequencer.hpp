/**
 * @file sequencer.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines sequencer models
 * @version 0.2
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

#ifndef __RACES_SEQUENCER__
#define __RACES_SEQUENCER__

#include <string>
#include <random>
#include <vector>
#include <cmath>

#include "genomic_position.hpp"
#include "genomic_sequence.hpp"

namespace Races
{

/**
 * @brief The namespace of sequencer models
 */
namespace Sequencers
{

/**
 * @brief A basic sequencer model
 */
class BasicSequencer
{
public:
    /**
     * @brief Get the model name
     *
     * @return the model name
     */
    virtual std::string get_model_name() const = 0;

    /**
     * @brief Get the platform name
     *
     * This method returns the sequencer platform as specified by [1].
     *
     * [1] "Sequence Alignment/Map Format Specification", The SAM/BAM Format Specification
     *     Working Group, 22 August 2022, http://samtools.github.io/hts-specs/SAMv1.pdf
     *
     * @return the platform name
     */
    virtual std::string get_platform_name() const = 0;

    /**
     * @brief Simulate the sequencing of a base
     *
     *
     * This method simulates the sequencing of a base. It applies an error to the
     * base according to the sequencer model and it also produces the corresponding
     * quality scores.
     *
     * @param[in,out] base is the base whose sequencing must be simulated
     * @param[in] base_position is the position of the base in the read
     * @param[in] read_position is the position of the first read base in the reference
     *      genome
     * @param[in] read_index is the read index, i.e., 0 for single read and the first
     *      read in paired reads, and 1 for the second read in paired reads
     * @return the quality score of the base in Sanger FASTQ format
     */
    virtual char simulate_seq(char& base, const size_t base_position,
                              const Mutations::GenomicPosition& read_position,
                              const unsigned int read_index) = 0;

    /**
     * @brief Simulate the sequencing of a read
     *
     * This method simulates the sequencing of a read. It applies some errors to the
     * read according to the sequencer model and it also produces the corresponding
     * quality scores.
     *
     * @param[in,out] read is the read whose sequencing must be simulated
     * @param[in] position is the position of the first read base in the reference
     *      genome
     * @param[in] read_index is the read index, i.e., 0 for single read and the first
     *      read in paired reads, and 1 for the second read in paired reads
     * @return the quality scores of the read in Sanger FASTQ format
     */
    virtual std::string simulate_seq(std::string& read,
                                     const Mutations::GenomicPosition& position,
                                     const unsigned int read_index) = 0;

    /**
     * @brief Test whether paired reads are supported by the sequencer
     *
     * @return `true` if and only if the sequencer supports paired reads
     */
    virtual inline bool supports_paired_reads() const
    {
        return false;
    }

    /**
     * @brief The destroyer
     */
    virtual inline ~BasicSequencer()
    {}
};

/**
 * @brief A namespace for Illumina sequencer
 */
namespace Illumina
{

/**
 * @brief A basic error-less Illumina sequencer
 */
class ErrorLessSequencer : public BasicSequencer
{
public:
    /**
     * @brief Construct a new error-less Illumina sequencer object
     */
    inline ErrorLessSequencer()
    {}

    /**
     * @brief Get the model name
     *
     * @return the model name
     */
    inline std::string get_model_name() const override
    {
        return "Error-less Illumina Sequencer";
    }

    /**
     * @brief Get the platform name
     *
     * This method returns the sequencer platform as specified by [1].
     *
     * [1] "Sequence Alignment/Map Format Specification", The SAM/BAM Format Specification
     *     Working Group, 22 August 2022, http://samtools.github.io/hts-specs/SAMv1.pdf
     *
     * @return the platform name
     */
    inline std::string get_platform_name() const override
    {
        return "ILLUMINA";
    }

    /**
     * @brief Simulate the sequencing of a base
     *
     *
     * This method simulates the sequencing of a base. It applies an error to the
     * base according to the sequencer model and it also produces the corresponding
     * quality scores.
     *
     * @param[in,out] base is the base whose sequencing must be simulated
     * @param[in] base_position is the position of the base in the read
     * @param[in] read_position is the position of the first read base in the reference
     *      genome
     * @param[in] read_index is the read index, i.e., 0 for single read and the first
     *      read in paired reads, and 1 for the second read in paired reads
     * @return the quality score of the base in Sanger FASTQ format
     */
    inline char simulate_seq(char& base, const size_t base_position,
                             const Mutations::GenomicPosition& read_position,
                             const unsigned int read_index) override
    {
        (void)base;
        (void)base_position;
        (void)read_position;
        (void)read_index;

        return 33;
    }

    /**
     * @brief Simulate the sequencing of a read
     *
     * This method simulates the sequencing of a read. It applies some errors to the
     * read according to the sequencer model and it also produces the corresponding
     * quality scores.
     *
     * @param[in,out] read is the read whose sequencing must be simulated
     * @param[in] position is the position of the first read base in the reference
     *      genome
     * @param[in] read_index is the read index, i.e., 0 for single read and the first
     *      read in paired reads, and 1 for the second read in paired reads
     * @return the quality scores of the read in Sanger FASTQ format
     */
    inline std::string simulate_seq(std::string& read,
                                    const Mutations::GenomicPosition& position,
                                    const unsigned int read_index) override
    {
        (void)position;
        (void)read_index;

        return std::string(read.size(), 33);
    }

    /**
     * @brief Test whether paired reads are supported by the sequencer
     *
     * @return `true` if and only if the sequencer supports paired reads
     */
    inline bool supports_paired_reads() const override
    {
        return true;
    }

    /**
     * @brief The destroyer
     */
    inline ~ErrorLessSequencer()
    {}
};

/**
 * @brief A basic implementation of Illumina sequencer model
 *
 * This class simulates the sequencing errors by using an error probability
 * that is independent from both the base position in the read and the
 * position in the reference genome. The method uniformly samples a value
 * in the interval [0,1]; whenever the sampled value is smaller than
 * the specified probability threshold a sequencing error is produced.
 * The quality score of a base is determined by using the sample value.
 */
template<typename RANDOM_GENERATOR=std::mt19937_64>
class BasicSequencer : public ErrorLessSequencer
{
    double error_rate;

    RANDOM_GENERATOR number_generator;

    std::uniform_real_distribution<double> error_dist;

    /**
     * @brief Get a wrong base
     *
     * @param orig is the original base
     * @return an uniformly selected base different
     *      from `orig`
     */
    char get_wrong_base(const char& orig)
    {
        std::uniform_int_distribution<size_t> offset_dist(1,3);

        auto pos = (GenomicSequence::get_base_index(orig)
                    + offset_dist(number_generator))%4;

        return GenomicSequence::DNA_bases[pos];
    }
public:
    /**
     * @brief Construct a new basic Illumina sequencer object
     *
     * @param error_rate is the error probability in each base
     * @param seed is the random number generator seed
     */
    explicit BasicSequencer(const double& error_rate, const int& seed=0):
        error_rate(error_rate), number_generator(seed),
        error_dist(0,1)
    {}

    /**
     * @brief Get the model name
     *
     * @return the model name
     */
    inline std::string get_model_name() const override
    {
        return "Basic Illumina Sequencer";
    }

    /**
     * @brief Simulate the sequencing of a base
     *
     *
     * This method simulates the sequencing of a base. It applies an error to the
     * base according to the sequencer model and it also produces the corresponding
     * quality scores.
     *
     * @param[in,out] base is the base whose sequencing must be simulated
     * @param[in] base_position is the position of the base in the read
     * @param[in] read_position is the position of the first read base in the reference
     *      genome
     * @param[in] read_index is the read index, i.e., 0 for single read and the first
     *      read in paired reads, and 1 for the second read in paired reads
     * @return the quality score of the base in Sanger FASTQ format
     */
    char simulate_seq(char& base, const size_t base_position,
                      const Mutations::GenomicPosition& read_position,
                      const unsigned int read_index) override
    {
        (void)base_position;
        (void)read_position;
        (void)read_index;

        if (base == 'N') {
            return 33 + 93;
        }

        auto error_sample = error_dist(number_generator);

        if (error_sample >= error_rate) {
            return 33;
        }

        base = get_wrong_base(base);
        if (error_sample == error_rate) {
            return 33 + 93;
        }

        const char qual = 34 - 10 * std::log10(1-error_sample/error_rate);

        if (qual > 33 + 93) {
            return 33 + 93;
        }

        return qual;
    }

    /**
     * @brief Simulate the sequencing of a read
     *
     * This method simulates the sequencing of a read. It applies some errors to the
     * read according to the sequencer model and it also produces the corresponding
     * quality scores.
     *
     * @param[in,out] read is the read whose sequencing must be simulated
     * @param[in] position is the position of the first read base in the reference
     *      genome
     * @param[in] read_index is the read index, i.e., 0 for single read and the first
     *      read in paired reads, and 1 for the second read in paired reads
     * @return the quality scores of the read in Sanger FASTQ format
     */
    std::string simulate_seq(std::string& read,
                             const Mutations::GenomicPosition& position,
                             const unsigned int read_index) override
    {
        std::string qual(read.size(), 33);

        auto read_it = read.begin();
        auto qual_it = qual.begin();
        size_t base_pos{0};
        while (read_it != read.end()) {
            *qual_it = simulate_seq(*read_it, base_pos, position, read_index);

            ++read_it;
            ++qual_it;
            ++base_pos;
        }

        return qual;
    }

    /**
     * @brief Get the error probability per base
     *
     * @return a constant reference to the error probability per base
     */
    inline const double& get_error_rate() const
    {
        return error_rate;
    }

    /**
     * @brief The destroyer
     */
    ~BasicSequencer()
    {}
};

}   // Illumina

}   // Sequencers

}   // Races

#endif // __RACES_SEQUENCER__