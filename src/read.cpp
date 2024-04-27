/**
 * @file read.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements sequencing reads
 * @version 0.1
 * @date 2024-04-27
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

Read::MutationIterator::MutationIterator(const std::map<GenomicPosition, SID>& germlines,
                            const std::map<GenomicPosition, SID>& passengers,
                            const std::map<GenomicPosition, SID>::const_iterator& germline_it,
                            const std::map<GenomicPosition, SID>::const_iterator& passenger_it):
    passengers{&passengers}, germlines{&germlines},
    p_it{passenger_it}, g_it{germline_it},
    p_end{passengers.end() != passenger_it},
    g_end{germlines.end() != germline_it}
{
    if (p_end) {
        p_first = (g_end ? p_it->first.position<g_it->first.position : true);
    } else {
        p_first = false;
    }
}

Read::MutationIterator::MutationIterator():
    passengers{nullptr}, germlines{nullptr}, p_it{}, g_it{},
    p_end{false}, g_end{false}
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
    if (is_end()) {
        return *this;
    }

    if (p_first) {
        ++p_it;

        if (p_it == passengers->end()) {
            p_first = false;
            p_end = false;

            return *this;
        }

        if (!g_end) {
            return *this;
        }
    } else {
        ++g_it;

        if (g_it == germlines->end()) {
            p_first = true;
            g_end = false;

            return *this;
        }

        if (!p_end) {
            return *this;
        }
    }

    p_first = (p_it->first.position < g_it->first.position);

    return *this;
}

Read::MutationIterator Read::MutationIterator::operator++(int)
{
    auto curr = *this;

    this->operator++();

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

std::string Read::get_CIGAR() const
{
    std::ostringstream oss;

    std::map<MatchingType, std::string> matching_str{
        {MatchingType::MATCH, "M"},
        {MatchingType::MISMATCH, "X"},
        {MatchingType::INSERTION, "I"},
        {MatchingType::DELETION, "D"},
    };

    MatchingType last_seq_type;
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

void Read::copy_reference(const std::string& reference, const size_t length,
                          size_t& read_idx, size_t& reference_idx)
{
    for (; reference_idx <= length && read_idx < nucleotides.size();
            ++reference_idx,++read_idx) {
        nucleotides[read_idx] = reference[reference_idx-1];
        alignment_index.push_back(alignment.size());
        if (nucleotides[read_idx]=='N') {
            alignment.push_back(MatchingType::MISMATCH);
        } else {
            alignment.push_back(MatchingType::MATCH);
        }
    }
}

inline
void update_alignment(std::vector<MatchingType>& alignment,
                      std::vector<size_t>& alignment_index,
                      const SID& mutation, const size_t& to_be_copied)
{
    auto num_of_matches = std::min(mutation.ref.size(),
                                   to_be_copied);
    size_t i=0;
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

#define DEBUG true

Read::Read(const std::string& reference,
           const std::map<GenomicPosition, SID>& germlines,
           const std::map<GenomicPosition, SID>& passengers,
           const GenomicPosition& genomic_position,
           const size_t& read_size):
    genomic_position{genomic_position}
{
    mutations = std::vector<SID>();
    nucleotides = std::string(read_size,'N');
    alignment.reserve(2*read_size);
    alignment_index.reserve(read_size);
    mutation_index = std::vector<MutationBaseIndex>(read_size);

    size_t read_idx=0, reference_idx=genomic_position.position;

    auto it = MutationIterator::lower_bound(germlines, passengers,
                                            genomic_position);

    while (!it.is_end() && read_idx < read_size
            && reference_idx <= reference.size()) {

        // copy from the reference up to the mutation begin
        const auto ref_up_to = std::min(static_cast<size_t>(it->first.position-1),
                                        reference.size());
        copy_reference(reference, ref_up_to, read_idx, reference_idx);

        if (read_idx < read_size) {
            const SID& mutation = it->second;

            if (DEBUG) {
                validate_mutation(reference, mutation);
            }

            // remove the reference from the read
            reference_idx += mutation.ref.size();

            // add the altered sequence to the read
            auto to_be_copied = std::min(mutation.alt.size(),
                                         read_size-read_idx);
            nucleotides.replace(read_idx, to_be_copied, mutation.alt);
            for (size_t i=0; i<to_be_copied; ++i) {
                mutation_index[read_idx] = {mutations.size(), i};
                ++read_idx;
            }
            mutations.push_back(mutation);

            update_alignment(alignment, alignment_index,
                             mutation, to_be_copied);
        }

        ++it;
    }

    // copy from the reference until the mutation
    copy_reference(reference, reference.size(), read_idx, reference_idx);

    nucleotides.resize(read_idx);
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
