/**
 * @file snv.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements Single Nucleotide Variation and related functions
 * @version 0.15
 * @date 2024-03-01
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
    GenomicPosition(), ref_base('?'), alt_base('?'), type(UNDEFINED)
{}

SNV::SNV(const ChromosomeId& chromosome_id, const ChrPosition& chromosomic_position,
         const char alt_base, const Type& type):
    SNV(chromosome_id, chromosomic_position, '?', alt_base, "", type)
{}

SNV::SNV(const ChromosomeId& chromosome_id, const ChrPosition& chromosomic_position,
         const char ref_base, const char alt_base, const Type& type):
    SNV(chromosome_id, chromosomic_position, ref_base, alt_base, "", type)
{}

SNV::SNV(const GenomicPosition& position, const char alt_base, const Type& type):
    SNV(position, '?', alt_base, "", type)
{}

SNV::SNV(const GenomicPosition& position, const char ref_base,
         const char alt_base, const Type& type):
    SNV(position, ref_base, alt_base, "", type)
{}

SNV::SNV(GenomicPosition&& position, const char alt_base, const Type& type):
    SNV(position, '?', alt_base, "", type)
{}

SNV::SNV(GenomicPosition&& position, const char ref_base,
         const char alt_base, const Type& type):
    SNV(position, ref_base, alt_base, "", type)
{}

SNV::SNV(const ChromosomeId& chromosome_id, const ChrPosition& chromosomic_position,
         const char alt_base, const std::string& cause, const Type& type):
    SNV(GenomicPosition(chromosome_id, chromosomic_position), alt_base,
        cause, type)
{}

SNV::SNV(const ChromosomeId& chromosome_id, const ChrPosition& chromosomic_position,
         const char ref_base, const char alt_base,
         const std::string& cause, const Type& type):
    SNV(GenomicPosition(chromosome_id, chromosomic_position), ref_base, alt_base,
        cause, type)
{}

SNV::SNV(const GenomicPosition& position, const char alt_base,
         const std::string& cause, const Type& type):
    SNV(position, '?', alt_base, cause, type)
{}

SNV::SNV(const GenomicPosition& position, const char ref_base,
         const char alt_base, const std::string& cause, const Type& type):
    GenomicPosition(position), ref_base(ref_base), alt_base(alt_base), cause(cause),
    type(type)
{
    if (!MutationalContext::is_a_base(alt_base)) {
        std::ostringstream oss;

        oss << "expected a base; got " << alt_base << std::endl;

        throw std::domain_error(oss.str());
    }

    if (ref_base == alt_base) {
        throw std::domain_error("SNV: the reference and altered bases are the same");
    }
}

SNV::SNV(GenomicPosition&& position, const char alt_base,
         const std::string& cause, const Type& type):
    SNV(position, alt_base, cause, type)
{}

SNV::SNV(GenomicPosition&& position, const char ref_base,
         const char alt_base, const std::string& cause, const Type& type):
    SNV(position, ref_base, alt_base, cause, type)
{}

std::string SNV::get_type_description(const SNV::Type& type)
{
    switch(type) {
      case SNV::Type::DRIVER:
        return "driver";
      case SNV::Type::PASSENGER:
        return "passenger";
      case SNV::Type::PRENEOPLASTIC:
        return "pre-neoplastic";
      case SNV::Type::GERMINAL:
        return "germinal";
      default:
        return "unknown";
    }
}


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
