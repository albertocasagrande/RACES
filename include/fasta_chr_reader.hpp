/**
 * @file fasta_chr_reader.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines classes to read chromosomes from FASTA streams
 * @version 0.3
 * @date 2024-03-11
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

#ifndef __RACES_FASTA_CHR_READER__
#define __RACES_FASTA_CHR_READER__

#include <string>
#include <istream>
#include <algorithm>

#include "genomic_position.hpp"

#include "fasta_reader.hpp"
#include "fasta_utils.hpp"

namespace Races
{

namespace IO
{

namespace FASTA
{

/**
 * @brief Filter non-chromosome sequence from a FASTA stream
 */
struct FilterNonChromosomeSequence : public SequenceFilter
{
    Mutations::ChromosomeId last_chr_id;   //!< the identifier of the last processed chromosome header

    /**
     * @brief Check whether an header corresponds to a chromosome
     * 
     * @param header is the header of a FASTA sequence
     * @return `true` if and only if the header does not correspond 
     *          to a chromosome
     */
    inline bool operator()(const std::string& header)
    {
        return !is_chromosome_header(header, last_chr_id);
    } 
};

/**
 * @brief Chromosome data
 * 
 * This template represents chromosome data. The objects of this class store 
 * chromosome information, i.e., name, header, size, and chromosome identifier,
 * and, depending on the parameter, which must be in the hierarchy of the class
 * `Races::IO::FASTA::SequenceInfo`, may also maintain the chromosome nucleic
 * sequence.
 * 
 * @tparam DATA_TYPE is the base type of the template. It must be a inherited
 *      from `Races::IO::FASTA::SequenceInfo`
 */
template<typename DATA_TYPE, std::enable_if_t<std::is_base_of_v<Races::IO::FASTA::SequenceInfo, DATA_TYPE>, bool> = true>
struct ChromosomeData : public DATA_TYPE
{
    Mutations::ChromosomeId chr_id;    //!< the chromosome id

    /**
     * @brief Read the next chromosome data from the FASTA stream
     * 
     * @tparam DATA_TYPE if the type of the data to be read from the stream
     * @param[in,out] FASTA_stream is the FASTA stream
     * @param[out] chr_data is the object that will be filled with the chromosome data
     *          if some chromosome is read from `FASTA_stream`
     * @param[in,out] progress_bar is a progress bar
     * @return `true` if and only if a chromosome sequence is read from `FASTA_stream`
     */
    static bool read(std::istream& FASTA_stream, ChromosomeData<DATA_TYPE>& chr_data, 
                     UI::ProgressBar& progress_bar)
    {
        using namespace Races::IO::FASTA;

        FilterNonChromosomeSequence filter;

        if (DATA_TYPE::read(FASTA_stream, chr_data, filter, progress_bar)) {
            chr_data.chr_id = filter.last_chr_id;

            return true;
        }

        return false;
    }

    /**
     * @brief Read the next chromosome data from the FASTA stream
     * 
     * @tparam DATA_TYPE if the type of the data to be read from the stream
     * @param FASTA_stream[in,out] is the FASTA stream
     * @param chr_data[out] is the object that will be filled with the chromosome data if
     *          some chromosome is read from `FASTA_stream`
     * @param progress_bar_stream is the output stream for the progress bar
     * @return `true` if and only if a chromosome sequence is read from `FASTA_stream`
     */
    inline static bool read(std::istream& FASTA_stream, ChromosomeData<DATA_TYPE>& chr_data,
                            std::ostream& progress_bar_stream=std::cout)
    {
        UI::ProgressBar progress_bar(progress_bar_stream, true);

        return read(FASTA_stream, chr_data, progress_bar);
    }
};

}   // FASTA

}   // IO

}   // Races

#endif // __RACES_FASTA_CHR_READER__