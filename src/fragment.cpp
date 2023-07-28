/**
 * @file fragment.cpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Implements genomic fragment
 * @version 0.5
 * @date 2023-07-28
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
#include <algorithm>  // std::transform

#include "fragment.hpp"


namespace Races 
{

namespace Passengers
{

using Length = Fragment::Length;

void validate_fragment(const Fragment& fragment)
{
    for (const auto& [allele_id, allele]: fragment.get_alleles()) {
        for (const auto& [pos, mutation] : allele) {
            if (!fragment.contains(pos)) {
                throw std::domain_error("the fragment does not contains one of the allele");
            }
        }
    }   
}

GenomicRegion::GenomicRegion():
    initial_pos(0), length(0)
{}

GenomicRegion::GenomicRegion(const ChromosomeId chromosome_id, const Length length):
    initial_pos(chromosome_id, 1), length(length)
{
    if (length == 0) {
        throw std::domain_error("the genomic region length must greater than 0");
    }
}

GenomicRegion::GenomicRegion(const GenomicPosition initial_pos, const Length length):
    initial_pos(initial_pos), length(length)
{
    if (length == 0) {
        throw std::domain_error("the genomic region length must greater than 0");
    }
}

bool GenomicRegion::overlaps(const GenomicRegion& genomic_region) const
{
    return (this->contains(genomic_region.get_begin())
            || genomic_region.contains(get_begin())
            || this->contains(genomic_region.get_end())
            || genomic_region.contains(get_end()));
}

GenomicRegion GenomicRegion::split(const GenomicPosition& split_point)
{
    if (!contains(split_point)) {
        throw std::domain_error("the fragment does not contain the split point");
    }

    if (split_point==initial_pos) {
        throw std::domain_error("the fragment initial position and the split point are the same");
    } 

    GenomicRegion new_region(split_point, length-(split_point.position-initial_pos.position));

    length = split_point.position - initial_pos.position;

    return new_region;
}

GenomicRegion& GenomicRegion::join(GenomicRegion& contiguous_region)
{
    if (follows(contiguous_region)) {
        std::swap(contiguous_region, *this);
    }

    if (!precedes(contiguous_region)) {
        throw std::domain_error("the two genomic regions are not contiguous");
    }

    length += contiguous_region.length;

    // reset contiguous_region length
    contiguous_region.length = 0;

    return *this;
}

Fragment::Fragment():
    GenomicRegion()
{}

Fragment::Fragment(const GenomicRegion& genomic_region, const size_t& num_of_alleles):
    GenomicRegion(genomic_region)
{
    for (AlleleId i = 0; i < num_of_alleles; ++i) {
        alleles[i] = Allele();
    }
}

Fragment::Fragment(const ChromosomeId chromosome_id, const Length length, const size_t& num_of_alleles):
    GenomicRegion(chromosome_id, length)
{
    for (AlleleId i = 0; i < num_of_alleles; ++i) {
        alleles[i] = Allele();
    }
}

Fragment::Fragment(const ChromosomeId chromosome_id, const Length length, const std::vector<Allele>& alleles):
    GenomicRegion(chromosome_id, length)
{
    AlleleId i = 0;
    for (const auto& allele : alleles) {
        (this->alleles)[i] = allele;
        ++i;
    }

    validate_fragment(*this);
}

Fragment::Fragment(const GenomicPosition initial_pos, const Length length, const size_t& num_of_alleles):
    GenomicRegion(initial_pos, length)
{
    for (size_t i = 0; i < num_of_alleles; ++i) {
        alleles[i] = Allele();
    }
}

Fragment::Fragment(const GenomicPosition initial_pos, const Length length, const std::vector<Allele>& alleles):
    GenomicRegion(initial_pos, length)
{
    AlleleId i = 0;
    for (const auto& allele : alleles) {
        (this->alleles)[i] = allele;
        ++i;
    }

    validate_fragment(*this);
}

bool Fragment::has_the_same_allele_ids(const Fragment& fragment) const
{
    if (num_of_alleles() != fragment.num_of_alleles()) {
        return false;
    }

    auto t_it = alleles.begin();
    auto f_it = fragment.alleles.begin();

    for (; t_it != alleles.end(); ++t_it, ++f_it) {
        if (t_it->first != f_it->first) {
            return false;
        }
    }

    return true;
}

std::set<AlleleId> Fragment::get_allele_ids() const
{
    std::set<AlleleId> allele_ids;

    std::transform(alleles.begin(), alleles.end(), std::inserter(allele_ids, allele_ids.end()),
        [](decltype(alleles)::value_type const &pair) {
            return pair.first;
        });
    
    return allele_ids;
}


Fragment Fragment::split(const GenomicPosition& split_point)
{
    if (!contains(split_point)) {
        throw std::domain_error("the fragment does not contain the split point");
    }

    if (split_point==initial_pos) {
        throw std::domain_error("the fragment initial position and the split point are the same");
    } 

    Fragment new_fragment(split_point, length-(split_point.position-initial_pos.position), 0);

    length = split_point.position - initial_pos.position;

    // for all the alleles in *this
    for (auto& [allele_id, allele]: alleles) {
        // search the iterator to the first mutation in the allele 
        // that lays in the new fragment
        auto in_new_it = allele.lower_bound(split_point);

        // create the new allele and get a reference to it
        auto& new_allele = new_fragment.alleles[allele_id];

        // while there are some mutations in allele that should 
        // lays in new_fragment
        while (in_new_it != allele.end()) {
            // extract the mutation from the original allele and
            // insert it in the new allele
            new_allele.insert(allele.extract(in_new_it++));
        }
    }

    return new_fragment;
}

Fragment& Fragment::join(Fragment& contiguous_fragment)
{
    if (!has_the_same_allele_ids(contiguous_fragment)) {
        throw std::domain_error("the two fragments differ in alleles");
    }

    if (follows(contiguous_fragment)) {
        std::swap(contiguous_fragment, *this);
    }

    if (!precedes(contiguous_fragment)) {
        throw std::domain_error("the two fragments are not contiguous");
    }

    auto allele_it = alleles.begin();

    // for all the alleles in contiguous_fragment
    for (auto& [allele_id, old_allele]: contiguous_fragment.alleles) {

        auto old_it = old_allele.begin();

        // while there are some mutations in old_allele
        while (old_it != old_allele.end()) {
            // extract the mutation from the old_allele and
            // insert it in *allele_it
            allele_it->second.insert(old_allele.extract(old_it++));
        }

        // move to the next allele in new fragment
        ++allele_it;
    }

    length += contiguous_fragment.length;

    // reset contiguous_fragment length and alleles
    contiguous_fragment.length = 0;
    contiguous_fragment.alleles.clear();

    return *this;
}

Fragment& Fragment::duplicate_allele(const AlleleId& allele_id, const AlleleId& new_allele_id)
{
    auto it = alleles.find(allele_id);
    if (it == alleles.end()) {
        throw std::out_of_range("unknown allele");
    }

    if (alleles.count(new_allele_id)>0) {
        throw std::runtime_error("the new allele is already present");
    }

    alleles[new_allele_id] = it->second;

    return *this;
}

Fragment& Fragment::remove_allele(const AlleleId& allele_id)
{
    auto it = alleles.find(allele_id);
    if (it == alleles.end()) {
        throw std::out_of_range("unknown allele");
    }

    alleles.erase(it);

    return *this;
}

bool Fragment::is_mutational_context_free(const GenomicPosition& genomic_position) const
{
    if (!contains(genomic_position)) {
        throw std::domain_error("the fragment does not contains the specified position");
    }
    
    // for all alleles
    for (const auto& [allele_id, allele]: alleles) {

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

bool Fragment::insert(SNV&& snv, const AlleleId allele_id)
{
    if (!contains(snv)) {
        throw std::domain_error("the fragment does not contains the specified SNV");
    }

    auto it = alleles.find(allele_id);
    if (it == alleles.end()) {
        throw std::out_of_range("unknown allele");
    }

    if (is_mutational_context_free(snv)) {
        (it->second)[snv] = std::move(snv);

        return true;
    }

    return false;
}

bool Fragment::remove_SNV(const GenomicPosition& genomic_position)
{
    if (!contains(genomic_position)) {
        throw std::domain_error("the fragment does not contains the specified genomic position");
    }

    for (auto& [allele_id, allele]: alleles) {
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

std::ostream& operator<<(std::ostream& out, const Races::Passengers::GenomicRegion& genomic_region)
{
    auto begin_pos = genomic_region.get_begin();

    out << "GenomicRegion(chr: " << static_cast<int>(begin_pos.chr_id)
        << ", begin: " << begin_pos.position
        << ", length: " << genomic_region.size() << ")";

    return out;
}

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
