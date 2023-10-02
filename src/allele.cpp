/**
 * @file allele.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements allele representation
 * @version 0.2
 * @date 2023-10-02
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

#include "allele.hpp"

namespace Races
{

namespace Passengers
{

AlleleFragment::AlleleFragment():
    GenomicRegion()
{}

AlleleFragment::AlleleFragment(const ChromosomeId& chromosome_id,
                               const ChrPosition& begin, const ChrPosition& end):
    GenomicRegion({chromosome_id, begin}, end-begin+1)
{}

AlleleFragment::AlleleFragment(const GenomicRegion& genomic_region):
    GenomicRegion(genomic_region)
{}

bool AlleleFragment::has_context_free(const GenomicPosition& genomic_position) const
{
    GenomicPosition g_pos(genomic_position);
    if (g_pos.position>0) {
        g_pos.position -= 1;
    }

    auto it = snvs.lower_bound(g_pos);

    if (it != snvs.end() && it->first.position < g_pos.position+3) {
        return false;
    }

    return true;
}

bool AlleleFragment::insert(const SNV& snv)
{
    if (!contains(snv)) {
        throw std::out_of_range("The allele fragment does not contain the SNV");
    }

    if (has_context_free(snv)) {
        snvs[snv] = snv;

        return true;
    }

    return false;
}

bool AlleleFragment::remove_SNV(const GenomicPosition& genomic_position)
{
    if (!contains(genomic_position)) {
        throw std::out_of_range("The allele fragment does not contain the SNV");
    }

    auto it = snvs.find(genomic_position);

    if (it == snvs.end()) {
        return false;
    }

    snvs.extract(it);

    return true;
}

AlleleFragment AlleleFragment::split(const GenomicPosition& split_point)
{
    if (!contains(split_point)) {
        throw std::out_of_range("The allele fragment does not contain the SNV");
    }

    if (get_initial_position() == split_point.position) {
        throw std::domain_error("The split point cannot be the initial position of the fragment");
    }

    AlleleFragment new_fragment(GenomicRegion::split(split_point));
    
    auto it = snvs.lower_bound(split_point);
    while (it != snvs.end()) {
        new_fragment.snvs.insert(snvs.extract(it++));
    }

    return new_fragment;
}

AlleleFragment AlleleFragment::copy(const GenomicRegion& genomic_region) const
{
    AlleleFragment new_fragment(genomic_region);

    auto final_position = genomic_region.get_final_position();
    for (auto it = snvs.lower_bound(genomic_region.get_begin());
            it != snvs.end() && it->first.position <= final_position; ++it) {
        new_fragment.snvs.insert(*it);
    }

    return new_fragment; 
}

Allele::Allele()
{}

Allele::Allele(const ChromosomeId& chromosome_id, const ChrPosition& begin, const ChrPosition& end):
    fragments{{{chromosome_id, begin}, AlleleFragment(chromosome_id, begin, end)}}
{}

Allele::Allele(const GenomicRegion& genomic_region):
    fragments{{genomic_region.get_begin(), AlleleFragment(genomic_region)}}
{}

bool Allele::contains(const GenomicPosition& genomic_position) const
{
    auto it = fragments.upper_bound(genomic_position);

    if (it != fragments.begin()) {
        --it;
    }

    if (it != fragments.end() && it->second.contains(genomic_position)) {
        return true;
    }

    return false;
}


bool Allele::contains(const GenomicRegion& genomic_region) const
{
    auto it = fragments.upper_bound(genomic_region.get_begin());

    if (it != fragments.begin()) {
        --it;
    }

    if (it != fragments.end() && it->second.contains(genomic_region)) {
        return true;
    }

    return false;
}

bool Allele::has_context_free(const GenomicPosition& genomic_position) const
{
    GenomicPosition g_pos(genomic_position);
    g_pos.position -=1;

    auto it = fragments.upper_bound(g_pos);

    if (it!=fragments.begin()) {
        --it;
    }

    return it->second.has_context_free(genomic_position);
}

bool Allele::insert(const SNV& snv)
{
    if (!has_context_free(snv)) {
        return false;
    }

    auto it = fragments.upper_bound(snv);
    if (it != fragments.begin()) {
        --it;
    }
    
    if (it != fragments.end() && it->second.contains(snv)) {
        return it->second.insert(snv);
    }

    return false;
}

bool Allele::remove_SNV(const GenomicPosition& genomic_position)
{
    auto it = fragments.upper_bound(genomic_position);
    if (it != fragments.begin()) {
        --it;
    }

    if (it != fragments.end() && it->second.contains(genomic_position)) {
        return it->second.remove_SNV(genomic_position);
    }

    return false;
}

Allele Allele::copy(const GenomicRegion& genomic_region) const
{
    Allele new_sequence;

    auto it = fragments.upper_bound(genomic_region.get_begin());
    if (it != fragments.begin()) {
        --it;
    }

    if (it == fragments.end() || !it->second.contains(genomic_region)) {
        throw std::domain_error("The allele does not fully contain the genomic region.");
    }

    auto new_fragment = it->second.copy(genomic_region);

    new_sequence.fragments[new_fragment.get_begin()] = std::move(new_fragment);

    return new_sequence;
}

bool Allele::remove(const GenomicRegion& genomic_region)
{
    auto it = fragments.upper_bound(genomic_region.get_begin());
    if (it != fragments.begin()) {
        --it;
    }

    if (it == fragments.end() || !it->second.contains(genomic_region)) {
        return false;
    }

    if (genomic_region.get_initial_position() > it->first.position) {
        auto new_fragment = it->second.split(genomic_region.get_begin());

        if (new_fragment.get_final_position()>genomic_region.get_final_position()) {
            GenomicPosition g_pos(genomic_region.get_end());
            g_pos.position += 1;

            new_fragment = new_fragment.split(g_pos);
        }

        fragments[new_fragment.get_begin()] = std::move(new_fragment);
    } else {
        if (it->second.get_final_position()>genomic_region.get_final_position()) {
            GenomicPosition g_pos(genomic_region.get_end());
            g_pos.position += 1;

            auto new_fragment = it->second.split(g_pos);

            fragments.extract(it);

            fragments[new_fragment.get_begin()] = std::move(new_fragment);
        }
    }

    return true;
}

std::map<GenomicPosition, SNV> Allele::get_SNVs() const
{
    std::map<GenomicPosition, SNV> snvs;

    for (const auto& [pos, fragment]: fragments) {
        for (const auto& [snv_pos, snv]: fragment.get_SNVs()) {
            snvs.insert({snv_pos, snv});
        }
    }

    return snvs;
}


Allele::Length Allele::size() const
{
    Length total_size{0};

    for (const auto& [pos, fragment]: fragments) {
        total_size += fragment.size();
    }

    return total_size;
}

}   // Passengers

}   // Races

namespace std
{

std::ostream& operator<<(std::ostream& os, const Races::Passengers::AlleleFragment& allele_fragment)
{
    os << "["<< allele_fragment.get_initial_position() << "-" 
       << allele_fragment.get_final_position() << "]{";

    std::string sep;
    for (const auto& [pos, snv]: allele_fragment.get_SNVs()) {
        os << sep << snv;
        sep = ",";
    }

    os << "}";

    return os; 
}

std::ostream& operator<<(std::ostream& os, const Races::Passengers::Allele& allele)
{
    std::string sep;
    for (const auto& [pos, fragment]: allele.get_fragments()) {
        os << sep << fragment;

        sep = ",";
    }

    return os;
}

} // std
