/**
 * @file genomic_position.cpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Implements genomic position and related functions
 * @version 0.3
 * @date 2023-07-31
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

#include <limits>

#include "genomic_position.hpp"

namespace Races 
{

namespace Passengers
{

GenomicPosition::GenomicPosition():
    GenomicPosition(0)
{}

GenomicPosition::GenomicPosition(const ChromosomeId& chromosome_id):
    GenomicPosition(chromosome_id, 0)
{}

GenomicPosition::GenomicPosition(const ChromosomeId& chromosome_id, const ChrPosition& position):
    chr_id(chromosome_id), position(position)
{}


ChromosomeId GenomicPosition::stochr(const std::string& chr_name)
{
    if (chr_name=="X") {
        return std::numeric_limits<ChromosomeId>::max()-1;
    }

    if (chr_name=="Y") {
        return std::numeric_limits<ChromosomeId>::max();
    }

    return static_cast<Passengers::ChromosomeId>(stoi(chr_name));
}

std::string GenomicPosition::chrtos(const ChromosomeId& chr_id)
{
    if (chr_id == std::numeric_limits<ChromosomeId>::max()-1) {
        return "X";
    }

    if (chr_id == std::numeric_limits<ChromosomeId>::max()) {
        return "Y";
    }

    return std::to_string(chr_id);
}

}   // Passengers

}   // Races

namespace std
{

std::ostream& operator<<(std::ostream& out, const Races::Passengers::GenomicPosition& genomic_position)
{
    out << "chr" <<  Races::Passengers::GenomicPosition::chrtos(genomic_position.chr_id) << "(" 
        << genomic_position.position << ")";

    return out;
}

}   // std
