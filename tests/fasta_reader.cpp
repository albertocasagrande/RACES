/**
 * @file fasta_reader.cpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Testing Races::IO::FASTA::Reader class
 * @version 0.3
 * @date 2023-08-12
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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE fasta_reader

#include <boost/test/unit_test.hpp>

#include <fstream>

#include "fasta_reader.hpp"


BOOST_AUTO_TEST_CASE(reader_creation)
{
    using namespace Races::IO::FASTA;

    BOOST_CHECK_NO_THROW(Sequence sequence);
    BOOST_CHECK_NO_THROW(SequenceInfo seq_info);
}

struct FASTAFixture
{
    std::map<std::string, std::string> sequences;

    FASTAFixture():
        sequences({
        {"1 test", "AACCCTAACCCTAACCCTA"},
        {"2 another test", "ACCCTAACCCTAAGCCCTACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA"},
        {"Test3", "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA"}
    })
    {}
};

BOOST_FIXTURE_TEST_SUITE( FASTA_test, FASTAFixture )

BOOST_AUTO_TEST_CASE(read_sequence_infos)
{
    using namespace Races::IO::FASTA;

    std::ifstream fasta_stream(FASTA_FILE, std::ios_base::in);

    SequenceInfo seq_info;
    for (const auto& [header, nucleotides]: sequences) {
        BOOST_CHECK(SequenceInfo::read(fasta_stream, seq_info));

        BOOST_CHECK_EQUAL(seq_info.header, header);
        BOOST_CHECK_EQUAL(seq_info.length, nucleotides.size());
    }

    BOOST_CHECK(!SequenceInfo::read(fasta_stream, seq_info));
}

BOOST_AUTO_TEST_CASE(read_sequences)
{
    using namespace Races::IO::FASTA;

    std::ifstream fasta_stream(FASTA_FILE, std::ios_base::in);

    Sequence sequence;
    for (const auto& [header, nucleotides]: sequences) {
        BOOST_CHECK(Sequence::read(fasta_stream, sequence));

        BOOST_CHECK_EQUAL(sequence.header, header);
        BOOST_CHECK_EQUAL(sequence.length, nucleotides.size());
        BOOST_CHECK_EQUAL(sequence.nucleotides, nucleotides);
    }

    BOOST_CHECK(!Sequence::read(fasta_stream, sequence));
}

struct SequenceFilterBySize: public Races::IO::FASTA::SequenceFilter
{
    size_t size;

    SequenceFilterBySize(const size_t size):
        SequenceFilter(), size(size)
    {}

    bool operator()(const std::string& header) const
    {
        return header.size() >= size;
    }
};

BOOST_AUTO_TEST_CASE(read_filtered_sequences)
{
    using namespace Races::IO::FASTA;

    std::ifstream fasta_stream(FASTA_FILE, std::ios_base::in);

    SequenceFilterBySize filter(strlen("2 another test"));
    Sequence sequence;
    
    for (const auto& [header, nucleotides]: sequences) {
        if (!filter(header)) {
            BOOST_CHECK(Sequence::read(fasta_stream, sequence, filter));

            BOOST_CHECK_EQUAL(sequence.header, header);
            BOOST_CHECK_EQUAL(sequence.length, nucleotides.size());
            BOOST_CHECK_EQUAL(sequence.nucleotides, nucleotides);
        } else {
            BOOST_CHECK_EQUAL(header, "2 another test");
        }
    }

    BOOST_CHECK(!Sequence::read(fasta_stream, sequence));
}

BOOST_AUTO_TEST_SUITE_END()
