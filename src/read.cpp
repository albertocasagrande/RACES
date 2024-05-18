/**
 * @file read.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements sequencing reads
 * @version 0.4
 * @date 2024-05-18
 *
 * @copyright Copyright (c) 2023-2024å
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

#include <map>

#include "read.hpp"

namespace Races
{

namespace Mutations
{

/**
 * @brief The sequencing simulation namespace
 */
namespace SequencingSimulations
{

void Read::MutationIterator::set_current_mutation()
{
    if (direction == Direction::FORWARD) {
        if (p_end) {
            p_it_curr = false;
        } else {
            p_it_curr = (g_end ? true : p_it->first.position<g_it->first.position);
        }
    } else {
        if (p_begin) {
            p_it_curr = false;
        } else {
            p_it_curr = (g_begin ? true : p_it->first.position > g_it->first.position);
        }
    }
}

Read::MutationIterator::MutationIterator(const std::map<GenomicPosition, SID>& germlines,
                            const std::map<GenomicPosition, SID>& passengers,
                            const std::map<GenomicPosition, SID>::const_iterator& germline_it,
                            const std::map<GenomicPosition, SID>::const_iterator& passenger_it):
    passengers{&passengers}, germlines{&germlines},
    p_it{passenger_it}, g_it{germline_it},
    direction{Direction::FORWARD},
    p_begin{passengers.begin() == passenger_it},
    p_end{passengers.end() == passenger_it},
    g_begin{germlines.begin() == germline_it},
    g_end{germlines.end() == germline_it}
{
    set_current_mutation();
}

Read::MutationIterator::MutationIterator():
    passengers{nullptr}, germlines{nullptr}, p_it{}, g_it{},
    direction{Read::Direction::FORWARD},
    p_begin{true}, p_end{true}, g_begin{true}, g_end{true}
{}

Read::MutationIterator
Read::MutationIterator::lower_bound(const std::map<GenomicPosition, SID>& germlines,
                                    const std::map<GenomicPosition, SID>& passengers,
                                    const GenomicPosition& genomic_position)
{
    auto g_it = germlines.lower_bound(genomic_position);
    auto p_it = passengers.lower_bound(genomic_position);

    return {germlines, passengers, g_it, p_it};
}

Read::MutationIterator&
Read::MutationIterator::operator++()
{
    if (direction == Direction::BACKWARD) {
        if (!p_end) {
            ++p_it;

            if (p_it == passengers->end()) {
                p_end = true;
            }
        }
        if (!g_end) {
            ++g_it;

            if (g_it == germlines->end()) {
                g_end = true;
            }
        }
        direction = Direction::FORWARD;

        set_current_mutation();

        return *this;
    }

    if (is_end()) {
        return *this;
    }

    if (p_it_curr) {
        ++p_it;
        p_begin = false;

        if (p_it == passengers->end()) {
            p_it_curr = false;
            p_end = true;

            return *this;
        }

        if (g_end) {
            return *this;
        }
    } else {
        ++g_it;
        g_begin = false;

        if (g_it == germlines->end()) {
            p_it_curr = true;
            g_end = true;

            return *this;
        }

        if (p_end) {
            return *this;
        }
    }

    p_it_curr = (p_it->first.position < g_it->first.position);

    return *this;
}

Read::MutationIterator Read::MutationIterator::operator++(int)
{
    auto curr = *this;

    this->operator++();

    return curr;
}

Read::MutationIterator&
Read::MutationIterator::operator--()
{
    if (direction == Direction::FORWARD) {
        if (!p_begin) {
            --p_it;

            if (p_it == passengers->begin()) {
                p_begin = true;
            }
        }
        if (!g_begin) {
            --g_it;

            if (g_it == germlines->begin()) {
                g_begin = true;
            }
        }
        direction = Direction::BACKWARD;

        set_current_mutation();

        return *this;
    }

    if (is_begin()) {
        return *this;
    }

    if (p_it_curr) {
        if (p_it == passengers->begin()) {
            p_it_curr = false;
            p_begin = true;

            return *this;
        }

        --p_it;
        p_end = false;

        if (g_begin) {
            return *this;
        }
    } else {
        if (g_it == germlines->begin()) {
            p_it_curr = true;
            g_begin = true;

            return *this;
        }

        --g_it;
        g_end = false;

        if (p_end) {
            return *this;
        }
    }

    p_it_curr = (p_it->first.position >= g_it->first.position);

    return *this;
}

Read::MutationIterator Read::MutationIterator::operator--(int)
{
    auto curr = *this;

    this->operator--();

    return curr;
}

Read::Read()
{}

size_t Read::Hamming_distance() const
{
    size_t total_mismatches{0};
    for (const auto& mismatch: alignment) {
        if (mismatch != MatchingType::MATCH) {
            ++total_mismatches;
        }
    }

    return total_mismatches;
}

std::list<GenomicRegion> Read::get_covered_reference_regions() const
{
    std::list<GenomicRegion> covered_list;

    if (alignment.size()==0) {
        return covered_list;
    }

    MatchingType last_seq_type = MatchingType::DELETION;
    GenomicPosition seq_begin(this->genomic_position);
    size_t seq_size{0};
    for (const auto& matching : alignment) {
        // if part of the sequence have been deleted
        if (matching == MatchingType::DELETION) {
            if (last_seq_type == MatchingType::MATCH) {
                covered_list.emplace_back(seq_begin, seq_size);
                seq_begin.position += seq_size;
            }
            ++seq_begin.position;
            seq_size = 0;
            last_seq_type = MatchingType::DELETION;
        } else if (matching != MatchingType::INSERTION) { // i.e., MATCH || MISMATCH
            ++seq_size;
            last_seq_type = MatchingType::MATCH;
        }
    }

    if (last_seq_type == MatchingType::MATCH) {
        covered_list.emplace_back(seq_begin, seq_size);
    }

    return covered_list;
}

std::string Read::get_CIGAR() const
{
    std::ostringstream oss;

    std::map<MatchingType, std::string> matching_str{
        {MatchingType::MATCH, "M"},
        {MatchingType::MISMATCH, "X"},
        {MatchingType::INSERTION, "I"},
        {MatchingType::DELETION, "D"},
    };

    MatchingType last_seq_type = MatchingType::MATCH;
    size_t last_seq=0;
    for (const auto& matching : alignment) {
        if (last_seq_type != matching) {
            if (last_seq>0) {
                oss << last_seq
                    << matching_str.at(last_seq_type);
                last_seq = 0;
            }
            last_seq_type = matching;
        }
        ++last_seq;
    }

    if (last_seq>0) {
        oss << last_seq << matching_str.at(last_seq_type);
    }

    return oss.str();
}

inline
void validate_mutation(const std::string& reference, const SID& mutation)
{
    const auto length = std::min(mutation.ref.size(),
                                 reference.size()-mutation.position-1);

    std::string_view ref_substr(reference.c_str()+mutation.position-1, length);

    if (ref_substr != mutation.ref) {
        std::ostringstream oss;

        oss << "Wrong mutation ref sequence: expecting \""
            << mutation.ref << "\" got \"" << ref_substr << "\".";
        throw std::runtime_error(oss.str());
    }
}

void Read::copy_reference(const std::string& reference, const size_t up_to_index,
                          size_t& read_end, size_t& ref_end)
{
    for (; ref_end <= up_to_index && read_end < nucleotides.size();
            ++ref_end,++read_end) {
        nucleotides[read_end] = reference[ref_end-1];
        alignment_index.push_back(alignment.size());
        if (nucleotides[read_end]=='N') {
            alignment.push_back(MatchingType::MISMATCH);
        } else {
            alignment.push_back(MatchingType::MATCH);
        }
    }
}

size_t get_common_prefix_size(const std::string& a, const std::string& b)
{
    const size_t max_size = std::min(a.size(), b.size());

    auto a_it = a.begin();
    auto b_it = b.begin();

    for (size_t i=0; i < max_size; ++i) {
        if (*a_it != *b_it) {
            return i;
        }
    }

    return max_size;
}

inline
void update_alignment(std::vector<MatchingType>& alignment,
                      std::vector<size_t>& alignment_index,
                      const SID& mutation, const size_t& to_be_copied)
{
    const size_t common_size = get_common_prefix_size(mutation.ref,
                                                      mutation.alt);

    auto num_of_matches = std::min(mutation.ref.size(),
                                   to_be_copied);
    size_t i=0;
    for (; i<common_size; ++i) {
        alignment_index.push_back(alignment.size());
        alignment.push_back(MatchingType::MATCH);
    }
    for (; i< num_of_matches; ++i) {
        alignment_index.push_back(alignment.size());
        alignment.push_back(MatchingType::MISMATCH);
    }
    for (size_t i=num_of_matches; i< mutation.ref.size(); ++i) {
        alignment.push_back(MatchingType::DELETION);
    }
    for (size_t i=num_of_matches; i< to_be_copied; ++i) {
        alignment_index.push_back(alignment.size());
        alignment.push_back(MatchingType::INSERTION);
    }
}

Read::Read(const std::string& reference,
           const std::map<GenomicPosition, SID>& germlines,
           const std::map<GenomicPosition, SID>& passengers,
           const GenomicPosition& genomic_position,
           const size_t& read_size):
    genomic_position{genomic_position}
{
    auto it = MutationIterator::lower_bound(germlines, passengers,
                                            genomic_position);

    // remove initial deletions from the read
    if (!it.is_begin()) {
        auto prev = it;

        --prev;

        if (!prev.is_begin()) {
            auto begin_ref = prev->second.get_region().get_final_position()+1;
            if (begin_ref>this->genomic_position.position) {
                this->genomic_position.position = begin_ref;
            }
        }
    }

    mutations = std::vector<SID>();
    nucleotides = std::string(read_size,'N');
    alignment.reserve(2*read_size);
    alignment_index.reserve(read_size);
    mutation_index = std::vector<MutationBaseIndex>(read_size);

    size_t read_end=0, ref_end=this->genomic_position.position;

    while (!it.is_end() && read_end < read_size
            && ref_end <= reference.size()) {

        // copy from the reference up to the mutation begin
        const auto ref_up_to = std::min(static_cast<size_t>(it->first.position-1),
                                        reference.size());
        copy_reference(reference, ref_up_to, read_end, ref_end);

        if (read_end < read_size) {
            const SID& mutation = it->second;

#ifdef __DEBUG__
            validate_mutation(reference, mutation);
#endif // __DEBUG__

            // remove the reference from the read
            ref_end += mutation.ref.size();

            // add the altered sequence to the read
            auto to_be_copied = std::min(mutation.alt.size(),
                                         read_size-read_end);
            nucleotides.replace(read_end, to_be_copied, mutation.alt);
            for (size_t i=0; i<to_be_copied; ++i) {
                mutation_index[read_end] = {mutations.size(), i};
                ++read_end;
            }

            update_alignment(alignment, alignment_index,
                             mutation, to_be_copied);

            mutations.push_back(mutation);
        }

        ++it;
    }

    // copy from the reference until the mutation
    copy_reference(reference, reference.size(), read_end, ref_end);

    nucleotides.resize(read_end);
}


void Read::alter_base(const size_t pos, const char base)
{
    if (nucleotides[pos] != base) {
        MutationBaseIndex mb_index = mutation_index[pos];

        if (mb_index.is_valid()) {
            auto& mutation = mutations[mb_index.index];

            mutation.alt[mb_index.base_pos] = base;
        }
        auto& match_type = alignment[alignment_index[pos]];
        if (match_type == MatchingType::MATCH) {
            match_type = MatchingType::MISMATCH;
        }
    }

    nucleotides[pos] = base;
}

}   // SequencingSimulations

}   // Mutations

}   // Races
