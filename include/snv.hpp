/**
 * @file snv.hpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Defines Single Nucleotide Variation
 * @version 0.1
 * @date 2023-07-22
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

#ifndef __RACES_SNV__
#define __RACES_SNV__

#include <ostream>

#include "genomic_position.hpp"

namespace Races 
{

namespace Passengers
{

/**
 * @brief A class to represent SNV
 */
struct SNV : public GenomicPosition
{
    char orig_base;         //!< the original base
    char mutated_base;      //!< the mutated base

    /**
     * @brief The empty constructor
     */
    SNV();

    /**
     * @brief A constructor
     * 
     * @param chromosome_id is the id of the SNV chromosome
     * @param chromosomic_position is the SNV position in the chromosome
     * @param position is the SNV position
     * @param original_base is the original base
     * @param mutated_base is the mutated base
     * @throw std::domain_error `original_base` and `mutated_base` are the same
     */
    SNV(const ChromosomeId& chromosome_id, const ChrPosition& chromosomic_position, 
        const char original_base, const char mutated_base);

    /**
     * @brief A constructor
     * 
     * @param genomic_position is the SNV genomic_position
     * @param original_base is the original base
     * @param mutated_base is the mutated base
     * @throw std::domain_error `original_base` and `mutated_base` are the same
     */
    SNV(const GenomicPosition& genomic_position, const char original_base, const char mutated_base);
};

}   // Passengers

}   // Races

namespace std
{

/**
 * @brief Write a SNV in a output stream
 * 
 * @param out is the output stream
 * @param snv is the SNV to stream
 * @return a reference of the updated stream
 */
std::ostream& operator<<(std::ostream& out, const Races::Passengers::SNV& snv);

}   // std


#endif // __RACES_SNV__