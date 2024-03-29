/**
 * @file snv.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements Single Nucleotide Variation and related functions
 * @version 0.17
 * @date 2024-03-29
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

#include <ostream>
#include <sstream>

#include "snv.hpp"

#include "context.hpp"
#include "genomic_sequence.hpp"

namespace Races
{

namespace Mutations
{

SNV::SNV():
    Mutation(), ref_base('?'), alt_base('?')
{}

SNV::SNV(const ChromosomeId& chr_id, const ChrPosition& chr_position,
         const char alt_base, const Mutation::Nature& nature):
    SNV(chr_id, chr_position, '?', alt_base, "", nature)
{}

SNV::SNV(const ChromosomeId& chr_id, const ChrPosition& chr_position,
         const char ref_base, const char alt_base, const Mutation::Nature& nature):
    SNV(chr_id, chr_position, ref_base, alt_base, "", nature)
{}

SNV::SNV(const GenomicPosition& position, const char alt_base,
         const Mutation::Nature& nature):
    SNV(position.chr_id, position.position, '?', alt_base, "", nature)
{}

SNV::SNV(const GenomicPosition& position, const char ref_base,
         const char alt_base, const Mutation::Nature& nature):
    SNV(position.chr_id, position.position, ref_base, alt_base, "", nature)
{}

SNV::SNV(GenomicPosition&& position, const char alt_base, const Mutation::Nature& nature):
    SNV(position.chr_id, position.position, '?', alt_base, "", nature)
{}

SNV::SNV(GenomicPosition&& position, const char ref_base,
         const char alt_base, const Mutation::Nature& nature):
    SNV(position.chr_id, position.position, ref_base, alt_base, "", nature)
{}

SNV::SNV(const ChromosomeId& chr_id, const ChrPosition& chr_position,
         const char alt_base, const std::string& cause, const Mutation::Nature& nature):
    SNV(chr_id, chr_position, '?', alt_base, cause, nature)
{}

SNV::SNV(const ChromosomeId& chr_id, const ChrPosition& chr_position,
         const char ref_base, const char alt_base,
         const std::string& cause, const Mutation::Nature& nature):
    Mutation(chr_id, chr_position, nature, cause),
    ref_base(ref_base), alt_base(alt_base)
{}

SNV::SNV(const GenomicPosition& position, const char alt_base,
         const std::string& cause, const Mutation::Nature& nature):
    SNV(position.chr_id, position.position, '?',
        alt_base, cause, nature)
{}

SNV::SNV(const GenomicPosition& position, const char ref_base,
         const char alt_base, const std::string& cause,
         const Mutation::Nature& nature):
    Mutation(position, nature, cause), ref_base(ref_base), alt_base(alt_base)
{
    if (!GenomicSequence::is_a_DNA_base(alt_base)) {
        std::ostringstream oss;

        oss << "expected a base; got " << alt_base << std::endl;

        throw std::domain_error(oss.str());
    }

    if (ref_base == alt_base) {
        throw std::domain_error("SNV: the reference and altered bases are the same");
    }
}

SNV::SNV(GenomicPosition&& position, const char alt_base,
         const std::string& cause, const Mutation::Nature& nature):
    SNV(position, alt_base, cause, nature)
{}

SNV::SNV(GenomicPosition&& position, const char ref_base,
         const char alt_base, const std::string& cause,
         const Mutation::Nature& nature):
    SNV(position, ref_base, alt_base, cause, nature)
{}


}   // Mutations

}   // Races

namespace std
{

bool less<Races::Mutations::SNV>::operator()(const Races::Mutations::SNV &lhs,
                                             const Races::Mutations::SNV &rhs) const
{
    {
        less<Races::Mutations::GenomicPosition> gp_op;

        if (gp_op(lhs, rhs)) {
            return true;
        }

        if (gp_op(rhs, lhs)) {
            return false;
        }
    }

    return (lhs.alt_base < rhs.alt_base);
}

std::ostream& operator<<(std::ostream& out, const Races::Mutations::SNV& snv)
{
    out << static_cast<Races::Mutations::GenomicPosition>(snv)
        << "[" << snv.ref_base << ">" <<  snv.alt_base << "]";

    return out;
}

}   // std
