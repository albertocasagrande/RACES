/**
 * @file cna.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements a class for copy number alterations
 * @version 0.9
 * @date 2024-03-12
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

#include "cna.hpp"

namespace Races
{

namespace Mutations
{

CNA::CNA()
{}

CNA::CNA(const GenomicRegion& region, const CNA::Type& type, 
         const AlleleId& source, const AlleleId& destination):
    region(region), source(source), dest(destination), type(type)
{}

}   // Mutations

}   // Races


namespace std
{

bool less<Races::Mutations::CNA>::operator()(const Races::Mutations::CNA &lhs,
                                             const Races::Mutations::CNA &rhs) const
{
    using namespace Races::Mutations;

    // differences in region
    {
        less<GenomicRegion> gr_op;

        if (gr_op(lhs.region, rhs.region)) {
            return true;
        }

        if (gr_op(rhs.region, lhs.region)) {
            return false;
        }
    }

    // differences in source
    {
        if (lhs.source<rhs.source) {
            return true;
        }

        if (rhs.source<lhs.source) {
            return false;
        }
    }

    // differences in destination
    {
        if (lhs.dest<rhs.dest) {
            return true;
        }

        if (rhs.dest<lhs.dest) {
            return false;
        }
    }

    // differences in type
    if ((lhs.type == CNA::Type::AMPLIFICATION)
            && (lhs.type != rhs.type)) {
        return true;
    }

    return false;
}

std::ostream& operator<<(std::ostream& out, const Races::Mutations::CNA& cna)
{
    out << "CNA(";
    using namespace Races::Mutations;
    switch(cna.type) {
        case CNA::Type::AMPLIFICATION:
            out << "\"A\"," << cna.region; 
            if (cna.source != RANDOM_ALLELE) {
                out << ", src: "<< Allele::format_id(cna.source);
            }
            if (cna.dest != RANDOM_ALLELE) {
                out << ", dst: "<< Allele::format_id(cna.dest);
            }
            break;
        case CNA::Type::DELETION:
            out << "\"D\"," << cna.region;
            if (cna.dest != RANDOM_ALLELE) {
                out << ", dst: "<< Allele::format_id(cna.dest);
            }
            break;
        default:
            throw std::runtime_error("Unsupported CNA type "
                                     + std::to_string(static_cast<int>(cna.type)));
    }

    out << ")";

    return out;
}

}   // std
