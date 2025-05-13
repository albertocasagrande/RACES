/**
 * @file fasta_chr_reader.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines classes to read chromosomes from FASTA streams
 * @version 1.1
 * @date 2025-05-13
 *
 * @copyright Copyright (c) 2023-2025
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

namespace RACES
{

namespace IO
{

namespace FASTA
{

/**
 * @brief Chromosome data
 *
 * This template represents chromosome data. The objects of this class store
 * chromosome information, i.e., name, header, size, and chromosome identifier,
 * and, depending on the parameter, which must be in the hierarchy of the class
 * `RACES::IO::FASTA::SequenceInfo`, may also maintain the chromosome nucleic
 * sequence.
 *
 * @tparam DATA_TYPE is the base type of the template. It must be a inherited
 *      from `RACES::IO::FASTA::SequenceInfo`
 */
template<typename DATA_TYPE, std::enable_if_t<std::is_base_of_v<RACES::IO::FASTA::SequenceInfo, DATA_TYPE>, bool> = true>
struct ChromosomeData : public DATA_TYPE
{
    Mutations::ChromosomeId chr_id;    //!< the chromosome id

    /**
     * @brief Get the identifier from the sequence header
     *
     * @param header is the header of the sequence
     * @return the identifier associated to the header name
     */
    static std::string get_id(const std::string& header)
    {
        Mutations::ChromosomeId chr_id;

        if (is_chromosome_header(header, chr_id)) {
            return Mutations::GenomicPosition::chrtos(chr_id);
        }

        return "";
    }

    /**
     * @brief Check whether the header is valid
     *
     * @param header is the header of the sequence
     * @return `true` if and only if `header` is a valid
     *      sequence header
     */
    static inline bool is_valid(const std::string& header)
    {
        Mutations::ChromosomeId chr_id;

        return is_chromosome_header(header, chr_id);
    }
};

template<>
bool Reader<ChromosomeData<SequenceInfo>>::read(ChromosomeData<SequenceInfo>& chr_info,
                                                RACES::UI::ProgressBar& progress_bar);

template<>
bool Reader<ChromosomeData<Sequence>>::read(ChromosomeData<Sequence>& chr,
                                            RACES::UI::ProgressBar& progress_bar);

template<>
inline const char* Index<ChromosomeData<Sequence>>::index_extension()
{
    return ".chi";
}

template<>
void Index<ChromosomeData<Sequence>>::save(std::ostream& output_stream) const;

template<>
Index<ChromosomeData<Sequence>>
Index<ChromosomeData<Sequence>>::load(std::istream& input_stream);

template<>
bool IndexedReader<ChromosomeData<Sequence>>::read(ChromosomeData<Sequence>& chr,
                                                   const std::string& chr_name,
                                                   RACES::UI::ProgressBar& progress_bar);

}   // FASTA

}   // IO

}   // RACES

#endif // __RACES_FASTA_CHR_READER__