/**
 * @file fasta_utils.cpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Implements support utilities for FASTA files
 * @version 0.3
 * @date 2023-08-03
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

#include <regex>

#include "fasta_utils.hpp"

namespace Races
{

namespace IO
{

namespace FASTA
{

bool EnsemblSeqNameDecoder::get_chromosome_region(const std::string& seq_name, Passengers::GenomicRegion& chr_region) const
{
    using namespace Passengers;

    const std::regex chr_regex("([0-9]+|X|Y) dna:chromosome chromosome:[a-zA-Z0-9]+:([0-9]+|X|Y):1:([0-9]+):1 .*");

    std::smatch m;

    if (!regex_match(seq_name, m, chr_regex)) {
        return false;
    }

    ChromosomeId chr_id = GenomicPosition::stochr(m[1].str());

    GenomicRegion::Length length = static_cast<GenomicRegion::Length>(std::stoi(m[3].str()));

    chr_region = GenomicRegion(chr_id, length);

    return true;
}

bool NCBISeqNameDecoder::get_chromosome_region(const std::string& seq_name, Passengers::GenomicRegion& chr_region) const
{
    using namespace Passengers;

    const std::regex chr_regex("NC_[0-9]+.[0-9]+ [A-Za-z0-9 ]+ chromosome ([0-9]+|X|Y), [.A-Za-z0-9]+ Primary Assembly");

    std::smatch m;

    if (!regex_match(seq_name, m, chr_regex)) {
        return false;
    }

    ChromosomeId chr_id = GenomicPosition::stochr(m[1].str());

    chr_region = GenomicRegion(chr_id, 1);

    return true;
}

}   // FASTA

}   // IO

}   // Races
