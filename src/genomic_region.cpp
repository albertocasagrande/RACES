/**
 * @file genomic_region.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements genomic region
 * @version 0.2
 * @date 2023-10-02
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

#include <iostream>

#include "genomic_region.hpp"

namespace Races 
{

namespace Passengers
{

using Length = GenomicRegion::Length;

GenomicRegion::GenomicRegion():
    initial_pos(0), length(0)
{}

GenomicRegion::GenomicRegion(const ChromosomeId chromosome_id, const Length length):
    initial_pos(chromosome_id, 1), length(length)
{
    if (length == 0) {
        throw std::domain_error("the genomic region length must greater than 0");
    }
}

GenomicRegion::GenomicRegion(const GenomicPosition initial_pos, const Length length):
    initial_pos(initial_pos), length(length)
{
    if (length == 0) {
        throw std::domain_error("the genomic region length must greater than 0");
    }
}

bool GenomicRegion::overlaps(const GenomicRegion& genomic_region) const
{
    return (this->contains(genomic_region.get_begin())
            || genomic_region.contains(get_begin())
            || this->contains(genomic_region.get_end())
            || genomic_region.contains(get_end()));
}

GenomicRegion GenomicRegion::split(const GenomicPosition& split_point)
{
    if (!contains(split_point)) {
        throw std::domain_error("the fragment does not contain the split point");
    }

    if (split_point==initial_pos) {
        throw std::domain_error("the fragment initial position and the split point are the same");
    } 

    GenomicRegion new_region(split_point, length-(split_point.position-initial_pos.position));

    length = split_point.position - initial_pos.position;

    return new_region;
}

GenomicRegion& GenomicRegion::join(GenomicRegion& contiguous_region)
{
    if (follows(contiguous_region)) {
        std::swap(contiguous_region, *this);
    }

    if (!precedes(contiguous_region)) {
        throw std::domain_error("the two genomic regions are not contiguous");
    }

    length += contiguous_region.length;

    // reset contiguous_region length
    contiguous_region.length = 0;

    return *this;
}


}   // Passengers

}   // Races

namespace std
{

std::ostream& operator<<(std::ostream& out, const Races::Passengers::GenomicRegion& genomic_region)
{
    auto begin_pos = genomic_region.get_begin();

    out << "GenomicRegion(chr: " << static_cast<int>(begin_pos.chr_id)
        << ", begin: " << begin_pos.position
        << ", length: " << genomic_region.size() << ")";

    return out;
}

}
