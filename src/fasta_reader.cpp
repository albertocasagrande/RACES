/**
 * @file fasta_reader.cpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Implementes a FASTA file reader and support structures
 * @version 0.2
 * @date 2023-07-25
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
#include <regex>

#include "fasta_reader.hpp"


namespace Races
{

namespace IO
{

namespace FASTA
{

Passengers::ChromosomeId stochr(const std::string& chr_name)
{
    if (chr_name=="X") {
        return 23;
    }

    if (chr_name=="Y") {
        return 24;
    }

    return static_cast<Passengers::ChromosomeId>(stoi(chr_name));
}

std::string chrtos(const Passengers::ChromosomeId& chr_id)
{
    switch(chr_id) {
        case 23:
            return "X";
        case 24:
            return "Y";
        default:
            return std::to_string(chr_id);
    }
}

bool get_chromosome_region(const std::string& seq_name, Passengers::GenomicRegion& chr_region)
{
    using namespace Passengers;

    const std::regex chr_regex("([0-9]+|X|Y) dna:chromosome chromosome:[a-zA-Z0-9]+:([0-9]+|X|Y):1:([0-9]+):1 .*");

    std::smatch m;

    if (!regex_match(seq_name, m, chr_regex)) {
        return false;
    }

    ChromosomeId chr_id = stochr(m[0].str());

    GenomicRegion::Length length = static_cast<GenomicRegion::Length>(std::stoi(m[2].str()));

    chr_region = GenomicRegion(chr_id, length);

    return true;
}

SequenceReader::nucleotide_iterator::nucleotide_iterator(SequenceReader& reader):
    reader(&reader), valid_info(true), seq_name(reader.get_sequence_name()),
    seq_pos(reader.get_position())
{}

SequenceReader::nucleotide_iterator::nucleotide_iterator():
    reader(nullptr), valid_info(false)
{}

SequenceReader::nucleotide_iterator::reference SequenceReader::nucleotide_iterator::operator*() const
{
    if (reader==nullptr) {
        throw std::domain_error("no more bases in the sequence");
    }

    if (reader->seq_pos != seq_pos) {
        throw std::domain_error("invalid iterator");
    }

    reader->refill_buffer();

    return reader->buffer[reader->buffer_pos];
}

SequenceReader::nucleotide_iterator::pointer SequenceReader::nucleotide_iterator::operator->() const
{
    if (reader==nullptr) {
        throw std::domain_error("no more bases in the sequence");
    }

    if (reader->seq_pos != seq_pos) {
        throw std::domain_error("invalid iterator");
    }

    return &(reader->buffer[reader->buffer_pos]);
}

SequenceReader::nucleotide_iterator& SequenceReader::nucleotide_iterator::operator++()
{
    if (reader==nullptr) {
        throw std::domain_error("no more bases in the sequence");
    }

    if (reader->seq_pos != seq_pos) {
        throw std::domain_error("invalid iterator");
    }

    ++reader->buffer_pos;
    
    if (!reader->refill_buffer()) {
        ++reader->seq_pos;
        ++seq_pos;
    } else {
        reader = nullptr;
    }

    return *this;
}

const std::string& SequenceReader::nucleotide_iterator::get_sequence_name() const
{
    if (!valid_info) {
        throw std::domain_error("no sequence read");
    }

    return seq_name;
}

const uint32_t& SequenceReader::nucleotide_iterator::get_position() const
{
    if (!valid_info) {
        throw std::domain_error("no sequence read");
    }

    return seq_pos;
}

SequenceReader::SequenceReader(std::istream& in):
    in(&in), seq_pos(1), buffer_pos(0), concluded(false)
{
    char tag_char;

    // discharge new lines
    do {
        if (in.eof()) {
            concluded = true;

            return;
        }
        in >> tag_char;
    } while (tag_char == '\n');

    // is a new sequence has NOT been found
    if (tag_char != '>') {

        // putback the last read character 
        in.putback(tag_char);

        // throw an exception
        std::ostringstream oss;
        oss << "SequenceReader: expected '>', got '" << tag_char << "'";
        throw std::domain_error(oss.str());
    }

    // read the sequence name
    getline(in, seq_name);
}

bool SequenceReader::refill_buffer()
{
    // if the sequence has been fully read
    if (concluded) {
        return true;
    }

    // if the buffer is not empty
    if (buffer_pos<buffer.size()) {
        return false;
    }

    // the buffer is empty; reset the buffer position
    buffer_pos = 0;

    char read_char;
    do {
        // if all the character in the file 
        // have been read
        if (in->eof()) {
            concluded = true;

            return true;
        }
        *in >> read_char;
    
    // repeat until the read character is a new line
    // or a space
    } while (read_char == '\n' || read_char == ' ' || read_char == 0);

    // the read character is not '\n'; put it back
    // in the stream
    in->putback(read_char);

    // if the read character is '>'
    if (read_char == '>') {
        // the considered sequence has been fully
        // read
        concluded = true;
    } else {  // otherwise

        // fill the buffer
        *in >> buffer;
    }

    // return `concluded`
    return concluded;
}

void SequenceReader::skip()
{
    do {
        buffer_pos = buffer.size();
    } while (!refill_buffer());
}

Reader::Reader(std::istream& in):
    in(in)
{}

std::istream& skip_to_next_sequence(std::istream& in)
{
    std::string buffer;

    do {
        char read_char;
        do {
            // if all the character in the file 
            // have been read
            if (in.eof()) {
                return in;
            }
            in >> read_char;
        
        // repeat until the read character is a new line
        // or a space
        } while (read_char == '\n' || read_char == ' ' || read_char == 0);

        // if the read character is '>'
        if (read_char == '>') {
            // the considered sequence has been fully
            // read
            in.putback(read_char);
            return in;
        }

        if (in.eof()) {
            return in;
        }
        getline(in, buffer);
    } while (true);
}

SequenceReader Reader::next_sequence_reader()
{
    skip_to_next_sequence(in);

    return SequenceReader(in);
}

}   // FASTA

}   // IO

}   // Races
