/**
 * @file genome_mutations.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Testing Races::Mutations::GenomeMutations class
 * @version 0.9
 * @date 2024-02-28
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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE genome_mutations

#include <boost/test/unit_test.hpp>

#include <fstream>

#include "genome_mutations.hpp"


BOOST_AUTO_TEST_CASE(chromosome_mutations_creation)
{
    using namespace Races::Mutations;

    BOOST_CHECK_NO_THROW(ChromosomeMutations());
    BOOST_CHECK_NO_THROW(ChromosomeMutations(1, 100, 3));
    BOOST_CHECK_NO_THROW(ChromosomeMutations({{1, 1}, 100}, 3));
}

BOOST_AUTO_TEST_CASE(genome_mutations_creation)
{
    using namespace Races::Mutations;

    BOOST_CHECK_NO_THROW(GenomeMutations());
    BOOST_CHECK_NO_THROW(GenomeMutations({{1, 100, 3}, {2, 100, 2}}));

    auto index = ContextIndex<>::build_index(FASTA_FILE);

    auto chr_regions = index.get_chromosome_regions();

    BOOST_CHECK_NO_THROW(GenomeMutations(chr_regions, 2));

    BOOST_CHECK_THROW(GenomeMutations({{1, 100, 3}, {1, 100, 2}}), std::domain_error);
}

struct GenomeMutationsFixture
{
    Races::Mutations::GenomeMutations genome_mutations;

    GenomeMutationsFixture()
    {
        using namespace Races::Mutations;

        auto index = ContextIndex<>::build_index(FASTA_FILE);

        auto chr_regions = index.get_chromosome_regions();

        genome_mutations = GenomeMutations(chr_regions, 2);
    }
};

BOOST_FIXTURE_TEST_SUITE( genome_mutation_tests, GenomeMutationsFixture )

BOOST_AUTO_TEST_CASE(genome_sizes)
{
    BOOST_CHECK_EQUAL(2*genome_mutations.size(), genome_mutations.allelic_size());
}

BOOST_AUTO_TEST_CASE(genome_insert_SNVs)
{
    using namespace Races::Mutations;

    SNV snv0(1, 32, 'A', 'G');  // chromosome 1, position 32
    
    SNV snv1(snv0.chr_id, snv0.position-1, 'T', 'G'), // at snv0 5'
        snv2(snv0.chr_id, snv0.position+1, 'G', 'C'), // at snv0 3'
        snv3(snv0.chr_id+1, snv0.position, 'T', 'G'); // snv0 and snv3 differ in chromosome

    auto test_genome_mutations = genome_mutations;

    // insert snv0 in allele 0: ok
    BOOST_CHECK(test_genome_mutations.insert(snv0, 0));

    {
        const auto chr = test_genome_mutations.get_chromosomes().at(snv0.chr_id);

        // snv in a wrong position
        SNV snv_err(snv0.chr_id, chr.size()+1, 'A', 'G');

        BOOST_CHECK_THROW(test_genome_mutations.insert(snv_err, 0), std::domain_error);
    }

    // insert snv0 in allele 1: true (free context, non-free context in a different allele)
    BOOST_CHECK(test_genome_mutations.insert(snv0, 1));

    // insert snv0 in allele 1: false (non-free context)
    BOOST_CHECK(!test_genome_mutations.insert(snv0, 0));

    // insert snv1 in allele 1: false (non-free context)
    BOOST_CHECK(!test_genome_mutations.insert(snv1, 1));

    // insert snv2 in allele 1: false (non-free context)
    BOOST_CHECK(!test_genome_mutations.insert(snv2, 1));

    // insert snv3 in allele 7: ko (unknown allele)
    BOOST_CHECK_THROW(test_genome_mutations.insert(snv3, 7), std::out_of_range);

    {
        const auto chr = test_genome_mutations.get_chromosomes().at(snv3.chr_id);

        auto chr_allele = chr.get_alleles();

        BOOST_CHECK_EQUAL(chr_allele.size(), 2);
        BOOST_CHECK_EQUAL(chr_allele.at(0).get_SNVs().size(), 0);
        BOOST_CHECK_EQUAL(chr_allele.at(1).get_SNVs().size(), 0);
    }

    // insert snv3 in allele 1: ok (same position, different chromosome)
    BOOST_CHECK(test_genome_mutations.insert(snv3, 1));

    {
        const auto chr = test_genome_mutations.get_chromosomes().at(snv3.chr_id);

        auto chr_allele = chr.get_alleles();

        BOOST_CHECK_EQUAL(chr_allele.size(), 2);
        BOOST_CHECK_EQUAL(chr_allele.at(0).get_SNVs().size(), 0);
        BOOST_CHECK_EQUAL(chr_allele.at(1).get_SNVs().size(), 1);
    }
}

BOOST_AUTO_TEST_CASE(genome_delete_SNVs)
{
    using namespace Races::Mutations;

    SNV snv0(1, 32, 'A', 'G');
    
    SNV snv1(snv0.chr_id, snv0.position-1, 'T', 'G');

    auto test_genome_mutations = genome_mutations;

    // remove snv0: false (no SNV in the position)
    BOOST_CHECK(!test_genome_mutations.remove_SNV(snv0));

    {
        const auto chr = test_genome_mutations.get_chromosomes().at(snv0.chr_id);

        // snv in a wrong position
        SNV snv_err(snv0.chr_id, chr.size()+1, 'A', 'G');

        BOOST_CHECK_THROW(test_genome_mutations.remove_SNV(snv_err), std::domain_error);
    }

    test_genome_mutations.insert(snv0, 0);

    {
        const auto chr = test_genome_mutations.get_chromosomes().at(snv0.chr_id);

        auto chr_allele = chr.get_alleles();

        BOOST_CHECK_EQUAL(chr_allele.size(), 2);
        BOOST_CHECK_EQUAL(chr_allele.at(0).get_SNVs().size(), 1);
        BOOST_CHECK_EQUAL(chr_allele.at(1).get_SNVs().size(), 0);
    }

    // remove snv1: false (no SNV in the position)
    BOOST_CHECK(!test_genome_mutations.remove_SNV(snv1));

    // remove snv0: true (SNV removed)
    BOOST_CHECK(test_genome_mutations.remove_SNV(snv0));

    {
        const auto chr = test_genome_mutations.get_chromosomes().at(snv0.chr_id);

        auto chr_allele = chr.get_alleles();

        BOOST_CHECK_EQUAL(chr_allele.size(), 2);
        BOOST_CHECK_EQUAL(chr_allele.at(0).get_SNVs().size(), 0);
        BOOST_CHECK_EQUAL(chr_allele.at(1).get_SNVs().size(), 0);
    }
}

BOOST_AUTO_TEST_CASE(genome_amplify_region)
{
    using namespace Races::Mutations;

    SNV snv0(1, 32, 'A', 'G');
    
    SNV snv1(snv0.chr_id, 64, 'A', 'T'),
        snv2(snv0.chr_id, 110, 'T', 'G');

    GenomicRegion gr_32_64(snv0, 33), // from 32 to 64
                  gr_10_110({snv0.chr_id, 10}, 101),  // from 10 to 110
                  gr_60_80({snv0.chr_id, 60}, 21);    // from 60 to 80

    auto test_genome_mutations = genome_mutations;

    // insert snv0 (chromosome1, position 32)
    test_genome_mutations.insert(snv0, 0);

    AlleleId new_allele_id;

    // amplify allele 0 from 32 to 64: ok (allele 0 available)
    BOOST_CHECK(test_genome_mutations.amplify_region(gr_32_64, 0, new_allele_id));

    {
        const auto chr = test_genome_mutations.get_chromosomes().at(gr_32_64.get_chromosome_id());

        auto chr_allele = chr.get_alleles();

        BOOST_CHECK_EQUAL(chr_allele.size(), 3);
        BOOST_CHECK_EQUAL(chr_allele.at(0).get_SNVs().size(), 1);
        BOOST_CHECK_EQUAL(chr_allele.at(1).get_SNVs().size(), 0);
        BOOST_CHECK_EQUAL(chr_allele.at(2).get_SNVs().size(), 1);
        BOOST_CHECK(chr_allele.at(0).get_SNVs()==chr_allele.at(2).get_SNVs());
    }

    {
        const auto chr = test_genome_mutations.get_chromosomes().at(snv0.chr_id);

        // snv in a wrong position
        GenomicRegion gr_err(snv0, chr.size()+1);

        BOOST_CHECK_THROW(test_genome_mutations.amplify_region(gr_err, 0, new_allele_id), 
                          std::domain_error);
    }

    // insert snv1 (chromosome1, position 64)
    test_genome_mutations.insert(snv1, 0);

    // amplify allele 0 from 10 to 110: ok (allele 0 available)
    BOOST_CHECK(test_genome_mutations.amplify_region(gr_10_110, 0, new_allele_id));

    // insert snv2 (chromosome1, position 110)
    test_genome_mutations.insert(snv2, 0);

    // amplify allele 0 from 60 to 80: ok (allele 0 available)
    BOOST_CHECK(test_genome_mutations.amplify_region(gr_60_80, 0, new_allele_id));

    {
        // now the chromosome contains 5 alleles: 
        // - snv0 lays in alleles 0, 2, and 3
        // - snv1 lays in alleles 0, 3
        // - snv2 lays in alleles 0, 4

        const auto chr = test_genome_mutations.get_chromosomes().at(gr_32_64.get_chromosome_id());

        auto chr_allele = chr.get_alleles();

        BOOST_CHECK_EQUAL(chr_allele.size(), 5);
        BOOST_CHECK_EQUAL(chr_allele.at(0).get_SNVs().size(), 3);
        BOOST_CHECK_EQUAL(chr_allele.at(1).get_SNVs().size(), 0);
        BOOST_CHECK_EQUAL(chr_allele.at(2).get_SNVs().size(), 1);
        BOOST_CHECK_EQUAL(chr_allele.at(3).get_SNVs().size(), 2);
        BOOST_CHECK_EQUAL(chr_allele.at(4).get_SNVs().size(), 1);
        BOOST_CHECK(chr_allele.at(2).get_SNVs()!=chr_allele.at(4).get_SNVs());
    }
}

BOOST_AUTO_TEST_CASE(genome_remove_region)
{
    using namespace Races::Mutations;

    SNV snv0(1, 32, 'A', 'G');
    
    SNV snv1(snv0.chr_id, 64, 'A', 'T'),
        snv2(snv0.chr_id, 110, 'T', 'G');

    GenomicRegion gr_32_64(snv0, 33), // from 32 to 64
                  gr_10_110({snv0.chr_id, 10}, 101),  // from 10 to 110
                  gr_60_80({snv0.chr_id, 60}, 21);    // from 60 to 80

    auto test_genome_mutations = genome_mutations;

    // insert snv0 (chromosome1, position 32)
    test_genome_mutations.insert(snv0, 0);

    AlleleId new_allele_id;

    // amplify from 32 to 64: ok (allele 0 available)
    test_genome_mutations.amplify_region(gr_32_64, 0, new_allele_id);

    // insert snv1 (chromosome1, position 64)
    test_genome_mutations.insert(snv1, 0);

    // amplify from 10 to 110: ok (allele 0 available)
    test_genome_mutations.amplify_region(gr_10_110, 0, new_allele_id);
    
    // insert snv2 (chromosome1, position 110)
    test_genome_mutations.insert(snv2, 0);

    // amplify from 60 to 80: ok (allele 0 available)
    test_genome_mutations.amplify_region(gr_60_80, 0, new_allele_id);

    // now the chromosome contains 5 alleles: 
    // - snv0 lays in alleles 0, 2, and 3
    // - snv1 lays in alleles 0, 3
    // - snv2 lays in alleles 0, 4

    // remove allele 0 from 60 to 80: ok
    BOOST_CHECK(test_genome_mutations.remove_region(gr_60_80, 0));

    // now the chromosome contains 5 alleles: 
    // - snv0 lays in alleles 0, 2, and 3
    // - snv1 lays in alleles 3
    // - snv2 lays in alleles 0, 4

    {
        const auto chr = test_genome_mutations.get_chromosomes().at(gr_32_64.get_chromosome_id());

        auto chr_allele = chr.get_alleles();

        BOOST_CHECK_EQUAL(chr_allele.size(), 5);
        BOOST_CHECK_EQUAL(chr_allele.at(0).get_SNVs().size(), 2);
        BOOST_CHECK_EQUAL(chr_allele.at(1).get_SNVs().size(), 0);
        BOOST_CHECK_EQUAL(chr_allele.at(2).get_SNVs().size(), 1);
        BOOST_CHECK_EQUAL(chr_allele.at(3).get_SNVs().size(), 2);
        BOOST_CHECK_EQUAL(chr_allele.at(4).get_SNVs().size(), 1);
        BOOST_CHECK(chr_allele.at(2).get_SNVs()!=chr_allele.at(4).get_SNVs());
    }

    {
        const auto chr = test_genome_mutations.get_chromosomes().at(snv0.chr_id);

        // snv in a wrong position
        GenomicRegion gr_err(snv0, chr.size()+1);

        BOOST_CHECK_THROW(test_genome_mutations.remove_region(gr_err, 0), std::domain_error);
    }

    // amplify allele 0 from 10 to 110: fail (allele 0 does not contain fragment 60-80 anymore)
    BOOST_CHECK(!test_genome_mutations.amplify_region(gr_10_110, 0, new_allele_id));

    // remove allele 0 from 60 to 80: fail (allele 0 does not contain fragment 60-80 anymore)
    BOOST_CHECK(!test_genome_mutations.remove_region(gr_60_80, 0));
}

BOOST_AUTO_TEST_SUITE_END()
