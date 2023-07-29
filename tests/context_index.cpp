/**
 * @file context_index.cpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Testing Races::Passengers::ContextIndex class
 * @version 0.1
 * @date 2023-07-29
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
#define BOOST_TEST_MODULE context_index

#include <boost/test/unit_test.hpp>

#include <fstream>

#include "context_index.hpp"


BOOST_AUTO_TEST_CASE(context_index_creation)
{
    using namespace Races::Passengers;

    BOOST_CHECK_NO_THROW(ContextIndex());

    BOOST_CHECK_NO_THROW(ContextIndex<>::build_index(FASTA_FILE));

    std::set<GenomicRegion> regions{{{2,115}, 20},
                                    {{1,5}, 73},
                                    {{2,247}, 11}};

    BOOST_CHECK_NO_THROW(ContextIndex<>::build_index(FASTA_FILE, regions));

    {
        std::ifstream is(FASTA_FILE, std::ios::in);

        BOOST_CHECK_NO_THROW(ContextIndex<>::build_index(is));
    }

    {
        std::ifstream is(FASTA_FILE, std::ios::in);

        BOOST_CHECK_NO_THROW(ContextIndex<>::build_index(is, regions));
    }

    {
        std::ifstream is(FASTA_FILE, std::ios::in);

        BOOST_CHECK_NO_THROW(ContextIndex<>::build_index(is));

        BOOST_CHECK_THROW(ContextIndex<>::build_index(is), std::runtime_error);
    }
}

template<typename ABSOLUTE_GENOMIC_POSITION>
std::set<Races::Passengers::GenomicPosition> get_genomic_positions(const Races::Passengers::ContextIndex<ABSOLUTE_GENOMIC_POSITION>& context_index,
                                                                   const Races::Passengers::MutationalContext& mutational_context)
{
    std::set<Races::Passengers::GenomicPosition> positions;

    for (const auto& abs_pos: context_index[mutational_context]) {
        positions.insert(context_index.get_genomic_position(abs_pos));
    }

    return positions;
}

struct ContextFixture
{
    using MutationalContext = Races::Passengers::MutationalContext;
    using GenomicPosition = Races::Passengers::GenomicPosition;

    std::map<MutationalContext, std::set<GenomicPosition> > test_positions;

    ContextFixture():
        test_positions{
            {"ACT",{{1,76},{2,263},{3,5}}},
            {"GCG",{{1,30},{3,8}}},
            {"TCC",{{1,83},{2,295}}},
            {"TCT",{{1,61},{1,107},{2,163},{2,165}}},
            {"GCT",{{1,81},{2,127},{2,170},{2,293}}},
            {"TCG",{{2,125}}}
        }
    {}
};


BOOST_FIXTURE_TEST_SUITE( context_index_test, ContextFixture )

BOOST_AUTO_TEST_CASE(context_index_whole_genome)
{
    using namespace Races::Passengers;

    auto context_index = ContextIndex<>::build_index(FASTA_FILE);

    for (const auto& [context_test, positions_test]: test_positions) {
        std::set<Races::Passengers::GenomicPosition> positions;

        if (positions_test.size() != 0) {
            BOOST_CHECK_NO_THROW(positions = get_genomic_positions(context_index, context_test));
        } else {
            BOOST_CHECK_THROW(positions = get_genomic_positions(context_index, context_test), std::domain_error);
        }

        BOOST_CHECK_EQUAL(positions.size(), positions_test.size());

        auto it = positions.begin();
        auto it_test = positions_test.begin();
        for (; it != positions.end(); ++it, ++it_test) {
            BOOST_CHECK_EQUAL(*it, *it_test);
        }

    }
}

bool in_regions(const std::set<Races::Passengers::GenomicRegion>& genomic_regions,
               const Races::Passengers::GenomicPosition& genomic_position)
{
    for (const auto& genomic_region: genomic_regions) {
        if (genomic_region.contains(genomic_position)) {
            return true;
        }
    }

    return false;
}

BOOST_AUTO_TEST_CASE(context_index_regions)
{
    using namespace Races::Passengers;

    const std::set<GenomicRegion> regions{{{2,115}, 20}, {{1,5}, 73},
                                          {{2,247}, 11}};

    decltype(test_positions) filtered_test_positions;

    for (const auto& [context_test, positions_test]: test_positions) {
        for (const auto& position_test: positions_test) {
            if (in_regions(regions, position_test)) {
                filtered_test_positions[context_test].insert(position_test);
            }
        }
    }

    auto context_index = ContextIndex<>::build_index(FASTA_FILE, regions);

    for (const auto& [context_test, positions_test]: filtered_test_positions) {
        std::set<Races::Passengers::GenomicPosition> positions;

        if (positions_test.size() != 0) {
            BOOST_CHECK_NO_THROW(positions = get_genomic_positions(context_index, context_test));
        } else {
            BOOST_CHECK_THROW(positions = get_genomic_positions(context_index, context_test), std::domain_error);
        }
        BOOST_CHECK_EQUAL(positions.size(), positions_test.size());

        auto it = positions.begin();
        auto it_test = positions_test.begin();
        for (; it != positions.end(); ++it, ++it_test) {
            BOOST_CHECK_EQUAL(*it, *it_test);
        }

    }
}

BOOST_AUTO_TEST_SUITE_END()