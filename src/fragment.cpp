/**
 * @file fragment.cpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Implements genomic fragment
 * @version 0.1
 * @date 2023-07-22
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

#include <iostream>

#include "fragment.hpp"


namespace Races 
{

namespace Passengers
{

using Allele = Fragment::Allele;
using Length = Fragment::Length;

void validate_fragment(const Fragment& fragment)
{
    if (fragment.size() == 0) {
        throw std::domain_error("the fragment length must greater than 0");
    }

    for (const auto& allele: fragment.get_alleles()) {
        for (const auto& [pos, mutation] : allele) {
            if (!fragment.contains(pos)) {
                throw std::domain_error("the fragment does not contains one of the allele");
            }
        }
    }   
}

Fragment::Fragment():
    initial_pos(0), length(0)
{}

Fragment::Fragment(const ChromosomeId chromosome_id, const Length length):
    Fragment(chromosome_id, length, {})
{}

Fragment::Fragment(const ChromosomeId chromosome_id, const Length length, const std::vector<Allele>& alleles):
    initial_pos(chromosome_id), length(length), alleles(alleles)
{
    validate_fragment(*this);
}

Fragment::Fragment(const GenomicPosition initial_pos, const Length length):
    Fragment(initial_pos, length, {})
{}

Fragment::Fragment(const GenomicPosition initial_pos, const Length length, const std::vector<Allele>& alleles):
    initial_pos(initial_pos), length(length), alleles(alleles)
{
    validate_fragment(*this);
}

Fragment Fragment::split(const GenomicPosition& split_point)
{
    if (!contains(split_point)) {
        throw std::domain_error("the fragment does not contain the split point");
    }

    if (split_point==initial_pos) {
        throw std::domain_error("the fragment initial position and the split point are the same");
    } 

    Fragment new_fragment(split_point, length-(split_point.position-initial_pos.position));
    new_fragment.alleles = std::vector<Allele>(num_of_alleles());

    length = split_point.position - initial_pos.position;

    /*** split the alleles ***/
    auto new_it = new_fragment.alleles.begin();

    // for all the alleles in *this
    for (auto& allele: alleles) {

        // search the iterator to the first mutation in the allele 
        // that lays in the new fragment
        auto in_new_it = allele.lower_bound(split_point);

        // while there are some mutations in allele that should 
        // lays in new_fragment
        while (in_new_it != allele.end()) {
            // extract the mutation from the original allele and
            // insert it in the new allele
            new_it->insert(allele.extract(in_new_it++));
        }

        // move to the next allele in new fragment
        ++new_it;
    }

    return new_fragment;
}

Fragment& Fragment::join(Fragment& contiguous_fragment)
{
    if (num_of_alleles() != contiguous_fragment.num_of_alleles()) {
        throw std::domain_error("the two fragments differ in number of alleles");
    }

    if (follows(contiguous_fragment)) {
        std::swap(contiguous_fragment, *this);
    }

    if (!precedes(contiguous_fragment)) {
        throw std::domain_error("the two fragments are not contiguous");
    }

    auto allele_it = alleles.begin();

    // for all the alleles in contiguous_fragment
    for (auto& old_allele: contiguous_fragment.alleles) {

        auto old_it = old_allele.begin();

        // while there are some mutations in old_allele
        while (old_it != old_allele.end()) {
            // extract the mutation from the old_allele and
            // insert it in *allele_it
            allele_it->insert(old_allele.extract(old_it++));
        }

        // move to the next allele in new fragment
        ++allele_it;
    }

    length += contiguous_fragment.length;

    // reset contiguous_fragment length and alleles
    contiguous_fragment.length = 0;
    contiguous_fragment.alleles.resize(0);

    return *this;
}

Fragment& Fragment::duplicate_allele(const size_t index)
{
    if (num_of_alleles()<=index) {
        throw std::out_of_range("unknown allele");
    }

    alleles.push_back(alleles[index]);

    return *this;
}

Fragment& Fragment::remove_allele(const size_t index)
{
    if (num_of_alleles()<=index) {
        throw std::out_of_range("unknown allele");
    }

    std::swap(alleles[index], alleles[num_of_alleles()-1]);

    alleles.pop_back();

    return *this;
}

bool Fragment::is_mutational_context_free(const GenomicPosition& genomic_position) const
{
    if (!contains(genomic_position)) {
        throw std::domain_error("the fragment does not contains the specified position");
    }
    
    // for all alleles
    for (const auto& allele: alleles) {

        // find the first mutation not occurring before `genomic_position`
        auto it = allele.lower_bound(genomic_position);

        if (it != allele.end()) {
            if (it->first.position - genomic_position.position < 2) {
                return false;
            }
        }

        if (it != allele.begin()) {
            --it;
            if (genomic_position.position - it->first.position < 2) {
                return false;
            }
        }
    }

    return true;
}

bool Fragment::insert(SNV&& snv, const size_t allele_index)
{
    if (!contains(snv)) {
        throw std::domain_error("the fragment does not contains the specified SNV");
    }

    if (num_of_alleles()<=allele_index) {
        throw std::out_of_range("unknown allele");
    }

    if (is_mutational_context_free(snv)) {
        alleles[allele_index][snv] = std::move(snv);

        return true;
    }

    return false;
}

bool Fragment::remove_SNV(const GenomicPosition& genomic_position)
{
    if (!contains(genomic_position)) {
        throw std::domain_error("the fragment does not contains the specified genomic position");
    }

    for (auto& allele: alleles) {
        if (allele.erase(genomic_position)>0) {
            return true;
        }
    }

    return false;
}


}   // Passengers

}   // Races

namespace std
{

std::ostream& operator<<(std::ostream& out, const Races::Passengers::Fragment& fragment)
{
    auto begin_pos = fragment.get_begin();

    out << "Fragment(chr: " << static_cast<int>(begin_pos.chr_id)
        << ", begin: " << begin_pos.position
        << ", length: " << fragment.size()
        << ", # of alleles: " << fragment.num_of_alleles() << ")";

    return out;
}

}
