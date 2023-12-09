/**
 * @file snv.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements Single Nucleotide Variation and related functions
 * @version 0.8
 * @date 2023-12-09
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

namespace Mutations
{

SNV::SNV():
    GenomicPosition(), context(), mutated_base('X')
{}

SNV::SNV(const ChromosomeId& chromosome_id, const ChrPosition& chromosomic_position, 
         const MutationalContext& context, const char mutated_base):
    SNV(chromosome_id, chromosomic_position, context, mutated_base, "")
{}

SNV::SNV(const GenomicPosition& position, const MutationalContext& context,
         const char mutated_base):
    SNV(position, context, mutated_base, "")
{}

SNV::SNV(GenomicPosition&& position, const MutationalContext& context, const char mutated_base):
    SNV(position, context, mutated_base, "")
{}

SNV::SNV(const ChromosomeId& chromosome_id, const ChrPosition& chromosomic_position, 
         const MutationalContext& context, const char mutated_base, const std::string& cause):
    SNV(GenomicPosition(chromosome_id, chromosomic_position), context, mutated_base, cause)
{}

SNV::SNV(const GenomicPosition& position, const MutationalContext& context,
         const char mutated_base, const std::string& cause):
    GenomicPosition(position), context(context), mutated_base(mutated_base), cause(cause)
{
    if (!MutationalContext::is_a_base(mutated_base)) {
        std::ostringstream oss;

        oss << "expected a base; got " << mutated_base << std::endl;

        throw std::domain_error(oss.str());
    }

    if (context.get_central_nucleotide() == mutated_base) {
        throw std::domain_error("the original and the mutated bases are the same");
    }
}

SNV::SNV(GenomicPosition&& position, const MutationalContext& context, const char mutated_base,
         const std::string& cause):
    SNV(position, context, mutated_base, cause)
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
    
    {
        less<Races::Mutations::MutationalContext> mc_op;
        
        if (mc_op(lhs.context, rhs.context)) {
            return true; 
        }

        if (mc_op(rhs.context, lhs.context)) {
            return false; 
        }
    }

    if (lhs.mutated_base < rhs.mutated_base) {
        return true;
    }

    if (lhs.mutated_base > rhs.mutated_base) {
        return false;
    }

    {
        less<std::string> str_op;
        
        return str_op(lhs.cause, rhs.cause);
    }
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
