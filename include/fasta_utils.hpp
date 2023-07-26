/**
 * @file fasta_utils.hpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Defines support utilities for FASTA files
 * @version 0.1
 * @date 2023-07-26
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

#ifndef __RACES_FASTA_UTILS__
#define __RACES_FASTA_UTILS__

#include <string>

#include "fragment.hpp"

namespace Races
{

namespace IO
{

namespace FASTA
{

/**
 * @brief A structure to decode chromosome sequences name
 */
struct SeqNameDecoder
{
    /**
     * @brief Extract the genomic region from a sequence name
     * 
     * @param seq_name is a FASTA sequence name
     * @param chr_region is the variable where the chromosome region will be placed
     * @return `true` if and only if `seq_name` correspond to a DNA chromosome sequence name
     */
    virtual bool get_chromosome_region(const std::string& seq_name, Passengers::GenomicRegion& chr_region) const = 0;
};

/**
 * @brief A structure to decode Ensembl chromosome sequences name 
 */
struct EnsemblSeqNameDecoder : public SeqNameDecoder
{
    /**
     * @brief Extract the genomic region from a sequence name
     * 
     * @param seq_name is a FASTA sequence name
     * @param chr_region is the variable where the chromosome region will be placed
     * @return `true` if and only if `seq_name` correspond to a DNA chromosome sequence 
     *      name in Ensembl format
     */
    bool get_chromosome_region(const std::string& seq_name, Passengers::GenomicRegion& chr_region) const override;
};

}   // FASTA

}   // IO

}   // Races

#endif // __RACES_FASTA_UTILS__