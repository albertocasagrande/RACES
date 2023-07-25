/**
 * @file fasta_reader.cpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Testing Races::IO::FASTA::Reader class
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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE fasta_reader

#include <boost/test/unit_test.hpp>

#include <fstream>

#include "fasta_reader.hpp"


BOOST_AUTO_TEST_CASE(reader_creation)
{
    using namespace Races::IO::FASTA;

    std::ifstream in(FASTA_FILE, std::ios_base::in);

    BOOST_CHECK_NO_THROW(Reader reader(in));
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

BOOST_AUTO_TEST_CASE(read_sequence_name)
{
    using namespace Races::IO::FASTA;

    std::ifstream in(FASTA_FILE, std::ios_base::in);

    Reader fasta_reader(in);

    auto it = sequences.begin();
    auto seq_reader = fasta_reader.next_sequence_reader();
    while (!fasta_reader.finished()) {

        if (it == sequences.end()) {
            throw std::runtime_error("unexpected sequence");
        }

        BOOST_CHECK_EQUAL(it->first, seq_reader.get_sequence_name());

        ++it;
        seq_reader = fasta_reader.next_sequence_reader();
    }

    BOOST_CHECK(it == sequences.end());
}

BOOST_AUTO_TEST_CASE(read_sequence_bases)
{
    using namespace Races::IO::FASTA;

    std::ifstream in(FASTA_FILE, std::ios_base::in);

    Reader fasta_reader(in);

    auto it = sequences.begin();
    auto seq_reader = fasta_reader.next_sequence_reader();
    while (!fasta_reader.finished()) {
        uint32_t pos=0; 
        const auto& seq = it->second;
        auto base_it = seq_reader.begin();

        while (base_it != seq_reader.end()) {
            BOOST_CHECK_EQUAL(seq[pos], *base_it);
            BOOST_CHECK_EQUAL(pos+1, base_it.get_position());
            BOOST_CHECK_EQUAL(it->first, base_it.get_sequence_name());

            ++base_it;
            ++pos;
        }
        BOOST_CHECK_EQUAL(seq.size(), pos);

        ++it;
        seq_reader = fasta_reader.next_sequence_reader();
    }

    BOOST_CHECK(it == sequences.end());
}


BOOST_AUTO_TEST_SUITE_END()
