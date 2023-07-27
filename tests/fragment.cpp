/**
 * @file fragment.cpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Some tests for Races::Passengers::Fragment
 * @version 0.3
 * @date 2023-07-27
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
#define BOOST_TEST_MODULE fragment

#include <boost/test/unit_test.hpp>

#include "fragment.hpp"


BOOST_AUTO_TEST_CASE(fragment_creation)
{
    using namespace Races::Passengers;

    GenomicPosition g_pos(12, 235);
    GenomicPosition snv_pos(g_pos.chr_id, g_pos.position+10);

    std::vector<Fragment::Allele> alleles(1);

    alleles[0][snv_pos] = SNV(snv_pos, 'A', 'C');

    Fragment::Length length = snv_pos.position+1000;

    BOOST_CHECK_NO_THROW(Fragment(1, length, 2));
    BOOST_CHECK_NO_THROW(Fragment(snv_pos.chr_id, length, alleles));
    BOOST_CHECK_NO_THROW(Fragment(g_pos, length, 2));
    BOOST_CHECK_NO_THROW(Fragment(g_pos, length, alleles));

    BOOST_CHECK_THROW(Fragment(1, 0, 2), std::domain_error);
    BOOST_CHECK_THROW(Fragment(1, length, alleles), std::domain_error);
    BOOST_CHECK_THROW(Fragment(snv_pos.chr_id, 0, alleles), std::domain_error);
    BOOST_CHECK_THROW(Fragment(snv_pos.chr_id, snv_pos.position-1, alleles), std::domain_error);
}

BOOST_AUTO_TEST_CASE(fragment_contains)
{
    using namespace Races::Passengers;

    GenomicPosition f_pos(12, 235);
    GenomicPosition gen_pos(f_pos.chr_id, f_pos.position+10);

    Fragment::Length length = gen_pos.position+1000;

    Fragment fragment(f_pos, length, 2);

    // ok
    BOOST_CHECK(fragment.contains(gen_pos));

    // wrong chromosome
    BOOST_CHECK(!fragment.contains({static_cast<uint8_t>(gen_pos.chr_id+1), gen_pos.position}));

    // after fragment
    BOOST_CHECK(!fragment.contains({gen_pos.chr_id, f_pos.position+length+10}));

    // before fragment
    BOOST_CHECK(!fragment.contains({gen_pos.chr_id, f_pos.position-10}));
}

BOOST_AUTO_TEST_CASE(fragment_follows)
{
    using namespace Races::Passengers;

    GenomicPosition f_pos_a(12, 235);
    Fragment fragment_a(f_pos_a, 1000, 2);

    auto f_begin_a = fragment_a.get_begin().position;
    auto f_end_a = f_begin_a + fragment_a.size()-1;

    Fragment fragment_b({f_pos_a.chr_id, f_end_a+1}, 1000, 2);

    BOOST_CHECK(fragment_b.follows(fragment_a));
    BOOST_CHECK(!fragment_a.follows(fragment_b));
    
    // changing chromosome
    fragment_b = Fragment({static_cast<uint8_t>(f_pos_a.chr_id+1), f_end_a+1}, 1000, 2);

    BOOST_CHECK(!fragment_b.follows(fragment_a));
    BOOST_CHECK(!fragment_a.follows(fragment_b));

    // non-overlapping / before
    fragment_b = Fragment({f_pos_a.chr_id, f_begin_a-100}, 10);

    BOOST_CHECK(!fragment_b.follows(fragment_a));
    BOOST_CHECK(!fragment_a.follows(fragment_b));

    // overlapping / begins before - ends before
    fragment_b = Fragment(fragment_b.get_begin(), fragment_a.size(), 2);

    BOOST_CHECK(!fragment_b.follows(fragment_a));
    BOOST_CHECK(!fragment_a.follows(fragment_b));

    // overlapping / begins before - ends after
    fragment_b = Fragment({f_pos_a.chr_id, f_begin_a-100}, fragment_a.size()+200, 2);

    BOOST_CHECK(!fragment_b.follows(fragment_a));
    BOOST_CHECK(!fragment_a.follows(fragment_b));

    // overlapping - matching begin / ends before
    fragment_b = Fragment(fragment_a.get_begin(), fragment_a.size()-10, 2);

    BOOST_CHECK(!fragment_b.follows(fragment_a));
    BOOST_CHECK(!fragment_a.follows(fragment_b));

    // overlapping - matching begin / ends after
    fragment_b = Fragment(fragment_a.get_begin(), fragment_a.size()+10, 2);

    BOOST_CHECK(!fragment_b.follows(fragment_b));
    BOOST_CHECK(!fragment_a.follows(fragment_b));

    // overlapping - matching begin / matching end

    BOOST_CHECK(!fragment_a.follows(fragment_a));
    BOOST_CHECK(!fragment_a.follows(fragment_a));
}

BOOST_AUTO_TEST_CASE(fragment_insert)
{
    using namespace Races::Passengers;

    std::vector<Fragment::Allele> alleles(4);

    GenomicPosition f_pos(12, 1000);
    Fragment fragment({12,100}, 1000, alleles);

    SNV snv0({f_pos.chr_id, f_pos.position+10}, 'A', 'C');

    // unknown allele
    BOOST_CHECK_THROW(fragment.insert(snv0, 7), std::out_of_range);

    // not in fragment - wrong chromosome
    snv0.chr_id = f_pos.chr_id + 1;
    BOOST_CHECK_THROW(fragment.insert(snv0, 3), std::domain_error);

    // not in fragment - before fragment
    snv0.chr_id = f_pos.chr_id;
    snv0.position = fragment.get_begin().position-1;
    BOOST_CHECK_THROW(fragment.insert(snv0, 3), std::domain_error);

    // not in fragment - after fragment
    snv0.position = fragment.get_end().position+2;
    BOOST_CHECK_THROW(fragment.insert(snv0, 3), std::domain_error);

    // ok
    snv0.position = f_pos.position+10;
    BOOST_CHECK(fragment.insert(snv0, 3));

    // mutational context not free
    SNV snv1(snv0);
    BOOST_CHECK(!fragment.insert(snv1, 3));
    BOOST_CHECK(!fragment.insert(snv1, 1));

    snv1.position = snv0.position + 1;
    BOOST_CHECK(!fragment.insert(snv1, 3));
    BOOST_CHECK(!fragment.insert(snv1, 1));

    snv1.position = snv0.position - 1;
    BOOST_CHECK(!fragment.insert(snv1, 3));
    BOOST_CHECK(!fragment.insert(snv1, 1));

    // ok
    snv1.position = snv0.position - 2;
    BOOST_CHECK(fragment.insert(snv1, 1));

    // ok
    SNV snv2(snv0);
    snv2.position = snv0.position + 2;
    BOOST_CHECK(fragment.insert(snv2, 1));
}


struct FragmentFixture {

    Races::Passengers::Fragment fragment;
    std::vector<Races::Passengers::SNV> snvs;

    FragmentFixture()
    {
        using namespace Races::Passengers;

        std::vector<Fragment::Allele> alleles(4);

        GenomicPosition f_pos(12, 1000);
        fragment = Fragment({12,100}, 1000, alleles);

        snvs.push_back(SNV({f_pos.chr_id, f_pos.position+10}, 'A', 'C'));
        fragment.insert(snvs.back(), 3);

        // mutational context not free
        SNV snv_tmp(snvs[0]);
        snv_tmp.position = snvs[0].position - 5;
        snvs.push_back(snv_tmp);
        fragment.insert(snvs.back(), 1);

        snv_tmp = snvs[0];
        snv_tmp.position = snvs[0].position + 5;
        snvs.push_back(snv_tmp);
        fragment.insert(snvs.back(), 1);
    }
};


BOOST_FIXTURE_TEST_SUITE( Fragment3SNPs, FragmentFixture )

BOOST_AUTO_TEST_CASE(fragment_remove_SNV)
{
    using namespace Races::Passengers;

    auto fragment_a = fragment;

    GenomicPosition g_pos(snvs[0]);
    // not in fragment - wrong chromosome
    g_pos.chr_id = snvs[0].chr_id + 1;
    BOOST_CHECK_THROW(fragment_a.remove_SNV(g_pos), std::domain_error);

    // not in fragment - before fragment
    g_pos.chr_id = snvs[0].chr_id;
    g_pos.position = fragment_a.get_begin().position-1;
    BOOST_CHECK_THROW(fragment_a.remove_SNV(g_pos), std::domain_error);

    // not in fragment - after fragment
    g_pos.position = fragment_a.get_end().position+2;
    BOOST_CHECK_THROW(fragment_a.remove_SNV(g_pos), std::domain_error);

    BOOST_CHECK(fragment_a.remove_SNV(snvs[0]));
    BOOST_CHECK(!fragment_a.remove_SNV(snvs[0]));
}

BOOST_AUTO_TEST_CASE(fragment_split)
{
    using namespace Races::Passengers;

    std::vector<Fragment::Allele> alleles(4);

    Fragment fragment_a(fragment), fragment_b;

    // wrong chromosome
    BOOST_CHECK_THROW(fragment_a.split({static_cast<uint8_t>(fragment_a.get_chromosome_id()+1), 
                                        fragment_a.get_begin().position}), std::domain_error);

    // before fragment
    BOOST_CHECK_THROW(fragment_a.split({fragment_a.get_chromosome_id(), 
                                        fragment_a.get_begin().position-1}), std::domain_error);

    // at the fragment begin
    BOOST_CHECK_THROW(fragment_a.split(fragment_a.get_begin()), std::domain_error);

    // after the fragment
    BOOST_CHECK_THROW(fragment_a.split({fragment_a.get_chromosome_id(), 
                                        fragment_a.get_end().position+1}), std::domain_error);

    // ok
    BOOST_CHECK_NO_THROW(fragment_b = fragment_a.split(fragment_a.get_end()));

    // check sizes
    BOOST_CHECK_EQUAL(fragment_a.size(),fragment.size()-1);
    BOOST_CHECK_EQUAL(fragment_b.size(),1);

    // fragment_c must follow fragment_b
    BOOST_CHECK(fragment_b.follows(fragment_a));

    // fragment_a, fragment_b, and fragment_c must have the
    // same number of alleles
    BOOST_CHECK_EQUAL(fragment.num_of_alleles(),fragment_a.num_of_alleles());
    BOOST_CHECK_EQUAL(fragment.num_of_alleles(),fragment_b.num_of_alleles());

    // fragment_b should not contain all the
    // fragment_a's SNVs
    for (size_t i=0; i< fragment_a.num_of_alleles(); ++i) {
        BOOST_CHECK_EQUAL(fragment_a[i].size(),fragment[i].size());
    }

    // fragment_c should not contain SVNs
    for (size_t i=0; i< fragment_b.num_of_alleles(); ++i) {
        BOOST_CHECK_EQUAL(fragment_b[i].size(),0);
    }

    fragment_a = fragment;
    BOOST_CHECK_NO_THROW(fragment_b = fragment_a.split({snvs[0].chr_id, snvs[0].position}));

    // check sizes
    BOOST_CHECK_EQUAL(fragment_a.size()+fragment_b.size(),fragment.size());

    // check initial positions
    BOOST_CHECK_EQUAL(fragment_a.get_begin(), fragment.get_begin());
    BOOST_CHECK_EQUAL(fragment_b.get_begin(), snvs[0]);

    // fragment_c must follow fragment_b
    BOOST_CHECK(fragment_b.follows(fragment_a));

    // fragment_b should contain snv1 SVNs
    BOOST_CHECK_EQUAL(fragment_a[0].size(),0);
    BOOST_CHECK_EQUAL(fragment_a[1].size(),1);
    BOOST_CHECK_EQUAL(fragment_a[2].size(),0);
    BOOST_CHECK_EQUAL(fragment_a[3].size(),0);
    BOOST_CHECK_EQUAL(fragment_a[1].at(snvs[1]), snvs[1]);

    // fragment_c should contain snv0 and svn2 SVNs
    BOOST_CHECK_EQUAL(fragment_b[0].size(),0);
    BOOST_CHECK_EQUAL(fragment_b[1].size(),1);
    BOOST_CHECK_EQUAL(fragment_b[2].size(),0);
    BOOST_CHECK_EQUAL(fragment_b[3].size(),1);
    BOOST_CHECK_EQUAL(fragment_b[3].at(snvs[0]), snvs[0]);
    BOOST_CHECK_EQUAL(fragment_b[1].at(snvs[2]), snvs[2]);
}

BOOST_AUTO_TEST_CASE(fragment_join)
{
    using namespace Races::Passengers;

    auto fragment_a = fragment;
    auto fragment_b = fragment_a.split(snvs[1]);
    
    auto f_a_begin = fragment_a.get_begin();

    auto fragment_err = Fragment({static_cast<uint8_t>(f_a_begin.chr_id+1), f_a_begin.position},
                                 fragment_a.size(), fragment_a.get_alleles());

    // not contiguous
    BOOST_CHECK_THROW(fragment_a.join(fragment_a), std::domain_error);
    BOOST_CHECK_THROW(fragment_err.join(fragment_a), std::domain_error);
    BOOST_CHECK_THROW(fragment_a.join(fragment_err), std::domain_error);

    // different number of alleles
    std::vector<Fragment::Allele> err_alleles(5);
    fragment_err = Fragment(f_a_begin, fragment_a.size(), err_alleles);
    BOOST_CHECK_THROW(fragment_b.join(fragment_err), std::domain_error);
    BOOST_CHECK_THROW(fragment_err.join(fragment_b), std::domain_error);

    err_alleles = std::vector<Fragment::Allele>(3);
    fragment_err = Fragment(f_a_begin, fragment_a.size(), err_alleles);
    BOOST_CHECK_THROW(fragment_b.join(fragment_err), std::domain_error);
    BOOST_CHECK_THROW(fragment_err.join(fragment_b), std::domain_error);

    // valid join
    Fragment joint1(fragment_a), joint2(fragment_b);
    BOOST_CHECK_NO_THROW(joint1.join(fragment_b));
    BOOST_CHECK_NO_THROW(joint2.join(fragment_a));

    BOOST_CHECK_EQUAL(joint1.get_begin(), fragment.get_begin());
    BOOST_CHECK_EQUAL(joint2.get_begin(), fragment.get_begin());
    BOOST_CHECK_EQUAL(joint1.size(), fragment.size());
    BOOST_CHECK_EQUAL(joint2.size(), fragment.size());
    BOOST_CHECK_EQUAL(joint1.num_of_alleles(), fragment.num_of_alleles());
    BOOST_CHECK_EQUAL(joint2.num_of_alleles(), fragment.num_of_alleles());

    for (size_t i=0; i<fragment.num_of_alleles(); ++i) {
        auto a_it = fragment[i].begin();
        auto j1_it = joint1[i].begin();
        auto j2_it = joint2[i].begin();
        while (a_it != fragment[i].end() 
                && j1_it != joint1[i].end() 
                && j2_it != joint2[i].end()) {
            BOOST_CHECK_EQUAL(j1_it->first, a_it->first);
            BOOST_CHECK_EQUAL(j1_it->second, a_it->second);
            BOOST_CHECK_EQUAL(j2_it->first, a_it->first);
            BOOST_CHECK_EQUAL(j2_it->second, a_it->second);
            ++a_it;
            ++j1_it;
            ++j2_it;
        }
        BOOST_CHECK(a_it==fragment[i].end());
        BOOST_CHECK(j1_it==joint1[i].end());
        BOOST_CHECK(j2_it==joint2[i].end());
    }
}

BOOST_AUTO_TEST_CASE(fragment_duplicate_allele)
{
    auto fragment_a = fragment;
    
    // wrong allele number
    BOOST_CHECK_THROW(fragment_a.duplicate_allele(fragment_a.num_of_alleles()), std::out_of_range);

    // duplicate
    BOOST_CHECK_NO_THROW(fragment_a.duplicate_allele(1));
    
    // check the existence of a new allele
    BOOST_CHECK_EQUAL(fragment.num_of_alleles()+1, fragment_a.num_of_alleles());

    // check that the old alleles did not change
    for (size_t i=0; i<fragment.num_of_alleles(); ++i) {
        auto it = fragment[i].begin();
        auto a_it = fragment_a[i].begin();
        while (it != fragment[i].end() && a_it != fragment_a[i].end()) {
            BOOST_CHECK_EQUAL(it->first, a_it->first);
            BOOST_CHECK_EQUAL(it->second, a_it->second);
            ++a_it;
            ++it;
        }
        BOOST_CHECK(it==fragment[i].end());
    }

    // check the new allele
    auto it_1 = fragment[1].begin();
    auto it_a_1 = fragment_a[1].begin();
    auto it_a_new = fragment_a[fragment.num_of_alleles()].begin();
    while (it_1 != fragment[1].end() 
            && it_a_1 != fragment_a[1].end()
            && it_a_new != fragment_a[fragment.num_of_alleles()].end()) {
        BOOST_CHECK_EQUAL(it_1->first, it_a_1->first);
        BOOST_CHECK_EQUAL(it_1->second, it_a_1->second);
        BOOST_CHECK_EQUAL(it_1->first, it_a_new->first);
        BOOST_CHECK_EQUAL(it_1->second, it_a_new->second);
        ++it_1;
        ++it_a_1;
        ++it_a_new;
    }
    BOOST_CHECK(it_1==fragment[1].end());
    BOOST_CHECK(it_a_1==fragment_a[1].end());
    BOOST_CHECK(it_a_new==fragment_a[fragment.num_of_alleles()].end());
}

BOOST_AUTO_TEST_CASE(fragment_remove_allele)
{
    auto fragment_a = fragment;

    fragment_a.duplicate_allele(1);

    auto fragment_b = fragment_a;

    // wrong allele number
    BOOST_CHECK_THROW(fragment_b.remove_allele(fragment_b.num_of_alleles()), std::out_of_range);

    // duplicate
    BOOST_CHECK_NO_THROW(fragment_b.remove_allele(1));
    
    // check the existence of a new allele
    BOOST_CHECK_EQUAL(fragment.num_of_alleles(), fragment_b.num_of_alleles());

    // check that the old alleles did not change
    for (size_t i=0; i<fragment.num_of_alleles(); ++i) {
        auto it = fragment[i].begin();
        auto b_it = fragment_b[i].begin();
        while (it != fragment[i].end() && b_it != fragment_a[i].end()) {
            BOOST_CHECK_EQUAL(it->first, b_it->first);
            BOOST_CHECK_EQUAL(it->second, b_it->second);
            ++b_it;
            ++it;
        }
        BOOST_CHECK(it==fragment[i].end());
        BOOST_CHECK(b_it==fragment_b[i].end());
    }
}

BOOST_AUTO_TEST_SUITE_END()