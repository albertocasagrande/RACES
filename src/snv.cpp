/**
 * @file snv.cpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Implements Single Nucleotide Variation and related functions
 * @version 0.3
 * @date 2023-07-28
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

#include <ostream>
#include <sstream>

#include "snv.hpp"

#include "context.hpp"

namespace Races 
{

namespace Passengers
{

SNV::SNV():
    GenomicPosition(), orig_base('X'), mutated_base('X')
{}

SNV::SNV(const ChromosomeId& chromosome_id, const ChrPosition& chromosomic_position, 
         const char original_base, const char mutated_base):
    SNV(GenomicPosition(chromosome_id, chromosomic_position), original_base, mutated_base)
{}

SNV::SNV(const GenomicPosition& position, const char original_base, const char mutated_base):
    GenomicPosition(position), orig_base(original_base), mutated_base(mutated_base)
{
    if (!MutationalContext::is_a_base(original_base)) {
        std::ostringstream oss;

        oss << "expected a base; got " << original_base << std::endl;

        throw std::domain_error(oss.str());
    }

    if (!MutationalContext::is_a_base(mutated_base)) {
        std::ostringstream oss;

        oss << "expected a base; got " << mutated_base << std::endl;

        throw std::domain_error(oss.str());
    }

    if (original_base == mutated_base) {
        std::domain_error("the original and the mutated bases are the same");
    }
}

SNV::SNV(GenomicPosition&& position, const char original_base, const char mutated_base):
    SNV(position, original_base, mutated_base)
{}

}   // Passengers

}   // Races

namespace std
{

bool less<Races::Passengers::SNV>::operator()(const Races::Passengers::SNV &lhs,
                                              const Races::Passengers::SNV &rhs) const
{
    less<Races::Passengers::GenomicPosition> op;

    return op(lhs, rhs);
}

std::ostream& operator<<(std::ostream& out, const Races::Passengers::SNV& snv)
{
    out << static_cast<Races::Passengers::GenomicPosition>(snv) << "["
        << snv.orig_base << ">" <<  snv.mutated_base << "]";

    return out;
}

}   // std
