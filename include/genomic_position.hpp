/**
 * @file genomic_position.hpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Defines genomic position and related functions
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

#ifndef __RACES_GENOMIC_POSITION__
#define __RACES_GENOMIC_POSITION__

#include <cstdint>
#include <functional>   // std::less
#include <ostream>

namespace Races 
{

namespace Passengers
{

/**
 * @brief 
 */
typedef uint8_t ChromosomeId;


typedef uint32_t ChrPosition;

/**
 * @brief A class to represent genomic position
 */
struct GenomicPosition
{
    ChromosomeId chr_id;        //!< the chromosome identifier
    ChrPosition position;   //!< the position in the chromosome

    /**
     * @brief The empty constructor
     * 
     * Build a fake genomic position having 0 as chromosome id
     */
    GenomicPosition();

    /**
     * @brief A constructor
     * 
     * @param chromosome_id is the chromosome identifier
     */
    GenomicPosition(const ChromosomeId& chromosome_id);

        /**
     * @brief A constructor
     * 
     * @param chromosome_id is the chromosome identifier
     * @param position is the position in the chromosome
     */
    GenomicPosition(const ChromosomeId& chromosome_id, const ChrPosition& position);

    /**
     * @brief Test whether two genomic position are the same 
     * 
     * @param genomic_position is the genomic position to test
     * @return `true` if and only if `*this` and `position` have 
     *          the same chromosome id and position 
     */
    inline bool operator==(const GenomicPosition& genomic_position) const
    {
        return ((chr_id == genomic_position.chr_id)
                && (genomic_position.position == position));
    }

    /**
     * @brief Test whether two genomic position differ
     * 
     * @param genomic_position is the genomic position to test
     * @return `true` if and only if `*this` and `position` differ 
     *          in chromosome id or position 
     */
    inline bool operator!=(const GenomicPosition& genomic_position) const
    {
        return !(*this == genomic_position);
    }
};


}   // Passengers

}   // Races

namespace std
{


template<>
struct std::less<Races::Passengers::GenomicPosition>
{
    inline bool operator()(const Races::Passengers::GenomicPosition &lhs,
                           const Races::Passengers::GenomicPosition &rhs) const
    {
        return ((lhs.chr_id<rhs.chr_id) || 
                ((lhs.chr_id==rhs.chr_id) && (lhs.position<rhs.position)));
    }
};

/**
 * @brief Write a genomic position in a output stream
 * 
 * @param out is the output stream
 * @param genomic_position is the position to stream
 * @return a reference of the updated stream
 */
std::ostream& operator<<(std::ostream& out, const Races::Passengers::GenomicPosition& genomic_position);

}   // std


#endif // __RACES_GENOMIC_POSITION__