/**
 * @file fasta_reader.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements a FASTA file reader and support structures
 * @version 0.7
 * @date 2023-11-01
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

#include <sstream>

#include "fasta_reader.hpp"

namespace Races
{

namespace IO
{

namespace FASTA
{

void SequenceInfo::filter_remaining_sequence(std::istream& fasta_stream, UI::ProgressBar& progress_bar)
{
    size_t counter{0};
    int c = fasta_stream.get();
    while (c != EOF && c != '>') {
        std::string line;
        getline(fasta_stream, line);

        if ((counter = (counter+1)%1000) == 0) {
            // let time elapse in progress bar
            progress_bar.set_progress(progress_bar.get_progress());
        }

        c = fasta_stream.get();
    }

    if (c == '>') {
        fasta_stream.unget();
    }
}

bool SequenceInfo::read(std::istream& fasta_stream, SequenceInfo& seq_info, UI::ProgressBar& progress_bar)
{
    SequenceFilter filter;

    return read(fasta_stream, seq_info, nullptr, filter, progress_bar);
}

bool SequenceInfo::read(std::istream& fasta_stream, SequenceInfo& seq_info)
{
    UI::ProgressBar progress_bar(true);

    return SequenceInfo::read(fasta_stream, seq_info, progress_bar);
}

bool Sequence::read(std::istream& fasta_stream, Sequence& sequence, UI::ProgressBar& progress_bar)
{
    SequenceFilter filter;

    return SequenceInfo::read(fasta_stream, sequence, &(sequence.nucleotides), filter, progress_bar);
}

bool Sequence::read(std::istream& fasta_stream, Sequence& sequence)
{
    UI::ProgressBar progress_bar(true);
    SequenceFilter filter;

    return Sequence::read(fasta_stream, sequence, filter, progress_bar);
}

}   // FASTA

}   // IO

}   // Races
