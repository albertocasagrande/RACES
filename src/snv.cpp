/**
 * @file snv.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements Single Nucleotide Variation and related functions
 * @version 0.13
 * @date 2024-01-18
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

namespace Races
{

namespace Mutations
{

SNV::SNV():
    GenomicPosition(), context(), mutated_base('X'), type(UNDEFINED)
{}

SNV::SNV(const ChromosomeId& chromosome_id, const ChrPosition& chromosomic_position,
         const char mutated_base, const Type& type):
    SNV(chromosome_id, chromosomic_position, MutationalContext(), mutated_base, "", type)
{}

SNV::SNV(const ChromosomeId& chromosome_id, const ChrPosition& chromosomic_position,
         const MutationalContext& context, const char mutated_base, const Type& type):
    SNV(chromosome_id, chromosomic_position, context, mutated_base, "", type)
{}

SNV::SNV(const GenomicPosition& position, const char mutated_base, const Type& type):
    SNV(position, MutationalContext(), mutated_base, "", type)
{}

SNV::SNV(const GenomicPosition& position, const MutationalContext& context,
         const char mutated_base, const Type& type):
    SNV(position, context, mutated_base, "", type)
{}

SNV::SNV(GenomicPosition&& position, const char mutated_base, const Type& type):
    SNV(position, MutationalContext(), mutated_base, "", type)
{}

SNV::SNV(GenomicPosition&& position, const MutationalContext& context,
         const char mutated_base, const Type& type):
    SNV(position, context, mutated_base, "", type)
{}

SNV::SNV(const ChromosomeId& chromosome_id, const ChrPosition& chromosomic_position,
         const char mutated_base, const std::string& cause, const Type& type):
    SNV(GenomicPosition(chromosome_id, chromosomic_position), mutated_base,
        cause, type)
{}

SNV::SNV(const ChromosomeId& chromosome_id, const ChrPosition& chromosomic_position,
         const MutationalContext& context, const char mutated_base,
         const std::string& cause, const Type& type):
    SNV(GenomicPosition(chromosome_id, chromosomic_position), context, mutated_base,
        cause, type)
{}

SNV::SNV(const GenomicPosition& position, const char mutated_base,
         const std::string& cause, const Type& type):
    SNV(position, MutationalContext(), mutated_base, cause, type)
{}

SNV::SNV(const GenomicPosition& position, const MutationalContext& context,
         const char mutated_base, const std::string& cause, const Type& type):
    GenomicPosition(position), context(context), mutated_base(mutated_base), cause(cause),
    type(type)
{
    if (!MutationalContext::is_a_base(mutated_base)) {
        std::ostringstream oss;

        oss << "expected a base; got " << mutated_base << std::endl;

        throw std::domain_error(oss.str());
    }

    if (context.is_defined() && context.get_central_nucleotide() == mutated_base) {
        throw std::domain_error("SNV: the original and the mutated bases are the same");
    }
}

SNV::SNV(GenomicPosition&& position, const char mutated_base,
         const std::string& cause, const Type& type):
    SNV(position, mutated_base, cause, type)
{}

SNV::SNV(GenomicPosition&& position, const MutationalContext& context,
         const char mutated_base, const std::string& cause, const Type& type):
    SNV(position, context, mutated_base, cause, type)
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

    return (lhs.mutated_base < rhs.mutated_base);
}

std::ostream& operator<<(std::ostream& out, const Races::Mutations::SNV& snv)
{
    std::string snv_sequence = snv.context.get_sequence();

    out << static_cast<Races::Mutations::GenomicPosition>(snv) << "("
        << snv_sequence[0]
        << "[" << snv_sequence[1] << ">" <<  snv.mutated_base << "]"
        << snv_sequence[2] << ")";

    return out;
}

}   // std
