/**
 * @file genomic_sequence.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements a structure to handle genomic sequence
 * @version 0.2
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

#include <sstream>

#include "genomic_sequence.hpp"

namespace Races
{

std::vector<char> GenomicSequence::DNA_bases = {'A', 'C', 'G', 'T'};

bool GenomicSequence::is_a_DNA_base(const char& base)
{
    switch(base) {
        case 'A':
        case 'a':
        case 'C':
        case 'c':
        case 'G':
        case 'g':
        case 'T':
        case 't':
            return true;
        default:
            return false;
    }
}

char GenomicSequence::get_complemented(const char& base)
{
    switch(base) {
        case 'A':
        case 'a':
            return 'T';
        case 'C':
        case 'c':
            return 'G';
        case 'G':
        case 'g':
            return 'C';
        case 'T':
        case 't':
            return 'A';
        default:
        {
            std::ostringstream oss;

            oss << "Unsupported base '" << base << "'";
            throw std::domain_error(oss.str());
        }
    }
}

void GenomicSequence::complement(std::string& sequence)
{
    for (char& nucleotide: sequence) {
        nucleotide = get_complemented(nucleotide);
    }
}

std::string GenomicSequence::get_complemented(const std::string& sequence)
{
    std::string complemented(sequence);

    GenomicSequence::complement(complemented);

    return complemented;
}

std::string GenomicSequence::get_reversed(const std::string& sequence)
{
    std::string reversed(sequence);

    GenomicSequence::reverse(reversed);

    return reversed;
}

std::string GenomicSequence::get_reverse_complemented(const std::string& sequence)
{
    std::string new_seq = GenomicSequence::get_reversed(sequence);

    GenomicSequence::complement(new_seq);

    return new_seq;
}

}   // Races
