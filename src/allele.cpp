/**
 * @file allele.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements allele representation
 * @version 1.0
 * @date 2024-06-10
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

#include "allele.hpp"

namespace RACES
{

namespace Mutations
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

bool AlleleFragment::has_context_free(const SID& mutation) const
{
    if (!strictly_contains(mutation)) {
        return false;
    }

    GenomicPosition g_pos{mutation};
    if (g_pos.position>0) {
        g_pos.position -= 1;
    }

    auto end = g_pos.position + mutation.ref.size() + 3;
    auto it = mutations.lower_bound(g_pos);
    if (it != mutations.end() && it->first.position < end) {
        return false;
    }

    if (it != mutations.begin()) {
        --it;

        return (it->first.position+it->second.ref.size() < g_pos.position+1);
    }

    return true;
}

bool AlleleFragment::insert(const SID& mutation)
{
    if (!contains(mutation)) {
        throw std::out_of_range("The allele fragment does not contain the SNV");
    }

    if (has_context_free(mutation)) {
        mutations[mutation] = mutation;

        return true;
    }

    return false;
}

bool AlleleFragment::remove_mutation(const GenomicPosition& genomic_position)
{
    if (!contains(genomic_position)) {
        throw std::out_of_range("The allele fragment does not contain the SNV");
    }

    auto it = mutations.find(genomic_position);

    if (it == mutations.end()) {
        return false;
    }

    mutations.extract(it);

    return true;
}

bool AlleleFragment::includes(const SID& mutation) const
{
    if (!contains(mutation)) {
        return false;
    }

    auto it = mutations.find(mutation);

    if (it == mutations.end()) {
        return false;
    }

    return (it->second.alt == mutation.alt)
            && (it->second.ref == mutation.ref);
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

    auto it = mutations.lower_bound(split_point);
    while (it != mutations.end()) {
        new_fragment.mutations.insert(mutations.extract(it++));
    }

    return new_fragment;
}

AlleleFragment AlleleFragment::copy(const GenomicRegion& genomic_region) const
{
    AlleleFragment new_fragment(genomic_region);

    auto final_position = genomic_region.get_final_position();
    for (auto it = mutations.lower_bound(genomic_region.get_begin());
            it != mutations.end() && it->first.position <= final_position; ++it) {
        new_fragment.mutations.insert(*it);
    }

    return new_fragment;
}

bool AlleleFragment::has_driver_mutations_in(const GenomicRegion& genomic_region) const
{
    auto final_position = genomic_region.get_final_position();
    for (auto it = mutations.lower_bound(genomic_region.get_begin());
            it != mutations.end() && it->first.position <= final_position; ++it) {
        if (it->second.nature == Mutation::DRIVER) {
            return true;
        }
    }

    return false;
}

Allele::Allele()
{}

Allele::Allele(const AlleleId& identifier, const std::list<AlleleId>& history):
    history(history)
{
    this->history.push_back(identifier);
}

Allele::Allele(const AlleleId& identifier, const ChromosomeId& chromosome_id,
               const ChrPosition& begin, const ChrPosition& end,
               const std::list<AlleleId>& history):
    fragments{{{chromosome_id, begin}, AlleleFragment(chromosome_id, begin, end)}}, history(history)
{
    this->history.push_back(identifier);
}

Allele::Allele(const AlleleId& identifier, const GenomicRegion& genomic_region,
               const std::list<AlleleId>& history):
    fragments{{genomic_region.get_begin(), AlleleFragment(genomic_region)}}, history(history)
{
    this->history.push_back(identifier);
}

template<typename KEY, typename VALUE>
typename std::map<KEY, VALUE>::const_iterator
find_not_after(const std::map<KEY, VALUE>& value_map, const KEY& key)
{
    auto it = value_map.upper_bound(key);

    if (it != value_map.begin()) {
        --it;
    }

    return it;
}

template<typename KEY, typename VALUE>
typename std::map<KEY, VALUE>::iterator
find_not_after(std::map<KEY, VALUE>& value_map, const KEY& key)
{
    auto it = value_map.upper_bound(key);

    if (it != value_map.begin()) {
        --it;
    }

    return it;
}

bool Allele::strictly_contains(const GenomicPosition& genomic_position) const
{
    auto it = find_not_after(fragments, genomic_position);

    return (it != fragments.end() && it->second.strictly_contains(genomic_position));
}

bool Allele::contains(const GenomicPosition& genomic_position) const
{
    auto it = find_not_after(fragments, genomic_position);

    return (it != fragments.end() && it->second.contains(genomic_position));
}

bool Allele::contains(const GenomicRegion& genomic_region) const
{
    auto it = find_not_after(fragments, genomic_region.get_begin());

    return (it != fragments.end() && it->second.contains(genomic_region));
}

bool Allele::includes(const SID& mutation) const
{
    auto it = find_not_after(fragments, static_cast<const GenomicPosition&>(mutation));

    return (it != fragments.end() && it->second.includes(mutation));
}

bool Allele::has_context_free(const SID& mutation) const
{
    GenomicPosition g_pos{mutation};
    g_pos.position -=1;

    auto it = find_not_after(fragments, g_pos);

    return it->second.has_context_free(mutation);
}

bool Allele::insert(const SID& mutation)
{
    if (!has_context_free(mutation)) {
        return false;
    }

    auto it = find_not_after(fragments, static_cast<const GenomicPosition&>(mutation));

    if (it != fragments.end() && it->second.contains(mutation)) {
        return it->second.insert(mutation);
    }

    return false;
}

bool Allele::remove_mutation(const GenomicPosition& genomic_position)
{
    auto it = find_not_after(fragments, genomic_position);

    if (it != fragments.end() && it->second.contains(genomic_position)) {
        return it->second.remove_mutation(genomic_position);
    }

    return false;
}

Allele Allele::copy(const AlleleId& new_allele_id, const GenomicRegion& genomic_region) const
{
    Allele new_sequence(new_allele_id, history);

    auto it = find_not_after(fragments, genomic_region.get_begin());

    if (it == fragments.end() || !it->second.contains(genomic_region)) {
        throw std::domain_error("The allele does not fully contain the genomic region.");
    }

    auto new_fragment = it->second.copy(genomic_region);

    new_sequence.fragments[new_fragment.get_begin()] = std::move(new_fragment);

    return new_sequence;
}

bool Allele::has_driver_mutations_in(const GenomicRegion& genomic_region) const
{
    auto it = find_not_after(fragments, genomic_region.get_begin());

    while (it != fragments.end() && !it->second.begins_after(genomic_region.get_end())) {
        if (it->second.has_driver_mutations_in(genomic_region)) {
            return true;
        }
        ++it;
    }

    return false;
}


bool Allele::remove(const GenomicRegion& genomic_region)
{
    auto it = find_not_after(fragments, genomic_region.get_begin());

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

std::map<GenomicPosition, SID> Allele::get_mutations() const
{
    std::map<GenomicPosition, SID> mutations;

    for (const auto& [pos, fragment]: fragments) {
        for (const auto& [mutation_pos, mutation]: fragment.get_mutations()) {
            mutations.insert({mutation_pos, mutation});
        }
    }

    return mutations;
}


Allele::Length Allele::size() const
{
    Length total_size{0};

    for (const auto& [pos, fragment]: fragments) {
        total_size += fragment.size();
    }

    return total_size;
}

std::string Allele::format_id(const RACES::Mutations::AlleleId& allele_id)
{
    if (allele_id == RANDOM_ALLELE) {
        return "NA";
    }

    return std::to_string(allele_id);
}

Allele Allele::duplicate_structure() const
{
    Allele duplicate;

    duplicate.history = history;

    for (const auto& [pos, fragment] : fragments) {
        AlleleFragment d_frag(static_cast<const GenomicRegion&>(fragment));
        duplicate.fragments.emplace(pos, std::move(d_frag));
    }

    return duplicate;
}

}   // Mutations

}   // RACES

namespace std
{

std::ostream& operator<<(std::ostream& os, const RACES::Mutations::AlleleFragment& allele_fragment)
{
    os << "["<< allele_fragment.get_initial_position() << "-"
       << allele_fragment.get_final_position() << "]{";

    std::string sep;
    for (const auto& [pos, mutation]: allele_fragment.get_mutations()) {
        os << sep << mutation;
        sep = ",";
    }

    os << "}";

    return os;
}

std::ostream& operator<<(std::ostream& os, const RACES::Mutations::Allele& allele)
{
    std::string sep;
    for (const auto& [pos, fragment]: allele.get_fragments()) {
        os << sep << fragment;

        sep = ",";
    }

    return os;
}

} // std
