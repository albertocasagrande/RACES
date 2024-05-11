/**
 * @file rs_index.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines a class to compute the repeated substring index
 * @version 0.1
 * @date 2024-05-11
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

#ifndef __RACES_RS_INDEX__
#define __RACES_RS_INDEX__

#include <string>
#include <vector>
#include <tuple>
#include <ostream>

#include "utils.hpp"
#include "archive.hpp"

namespace Races
{

namespace Mutations
{

/**
 * @brief Repeated sequences
 *
 * This class represents repeated sequences.
 *
 * @tparam REPETITION_TYPE is the number of repetitions type
 */
template<typename REPETITION_TYPE=uint8_t>
struct Repetition
{
    GenomicPosition begin;  //!< The genomic position of the repeated sequence
    REPETITION_TYPE num_of_repetitions; //!< The number of repetitions of the unit
    std::string unit;   //!< The repeated subsequence
    char prev_base;     //!< The base preceding the first repetition of the unit

    /**
     * @brief The empty constructor
     */
    Repetition():
        begin(), num_of_repetitions(0), unit(""), prev_base('N')
    {}

    /**
     * @brief A constructor
     *
     * @param chr_id is the identifier of the chromosome containing
     *      the repeated sequence
     * @param begin is the position of the repeated sequence first base
     *      in the chromosome
     * @param num_of_repetitions is the number of repetition of the unit
     * @param unit is the pointer to the unit first base in the sequence
     * @param unit_size is the unit size
     * @param previous_base is base preceding the first repetition of
     *      the unit
     */
    Repetition(const ChromosomeId& chr_id, const ChrPosition begin,
               const REPETITION_TYPE num_of_repetitions,
               const char* unit, const REPETITION_TYPE unit_size,
               const char previous_character):
        begin(chr_id, begin), num_of_repetitions(num_of_repetitions),
        unit(unit, unit+unit_size), prev_base(previous_character)
    {
        if (unit_size==0) {
            throw std::domain_error("Unit size must be greater than 0.");
        }
    }

    /**
     * @brief Save a repetition in an archive
     *
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        archive & begin
                & num_of_repetitions
                & unit
                & prev_base;
    }

    /**
     * @brief Load a repetition from an archive
     *
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the loaded repetition
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    inline static Repetition load(ARCHIVE& archive)
    {
        Repetition repetition;

        archive & repetition.begin
                & repetition.num_of_repetitions
                & repetition.unit
                & repetition.prev_base;

        return repetition;
    }
};

}   // Mutations

}   //  Races


namespace std
{
    template<typename REPETITION_TYPE>
    std::ostream& operator<<(std::ostream& os, const Races::Mutations::Repetition<REPETITION_TYPE>& repetition)
    {
        os << "[" << std::string(1, repetition.prev_base) << "]" << repetition.unit << " ("
            << repetition.begin << ":" << static_cast<size_t>(repetition.num_of_repetitions) << ")";

        return os;
    }
}

namespace Races
{

namespace Mutations
{


/**
 * @brief An index of the genome repeated sequences
 *
 * This index stores the some of the repeated sequences of genomic sequences.
 *
 * @tparam REPETITION_TYPE is the number of repetitions type
 */
template<typename REPETITION_TYPE=uint8_t>
struct RSIndex
{
    /**
     * @brief A repeated sequence storage
     *
     * The objects of this class store a sample of the identified
     * repeated sequences.
     */
    struct RepetitionStorage
    {
        size_t total_number;        //!< The total number of the identified repeated sequences
        size_t stored_repetitions;  //!< The number of repetitions stored in storage
        std::vector<Repetition<REPETITION_TYPE>> repetitions;  //!< The samples repeated sequences

        /**
         * @brief The empty constructor
         */
        RepetitionStorage():
            total_number(0), stored_repetitions(0)
        {}

        /**
         * @brief Extract a repetition
         *
         * @param pos is the position of the repetition to be removed
         * @return a constant reference to the extracted repetition
         */
        const Repetition<REPETITION_TYPE>& extract(const size_t& pos)
        {
            if (pos >= stored_repetitions) {
                throw std::out_of_range("Position not available in the storage");
            }

            if (pos+1 != stored_repetitions) {
                std::swap(repetitions[pos], repetitions[stored_repetitions-1]);
            }

            --stored_repetitions;

            return repetitions[pos];
        }

        /**
         * @brief Insert a repetition
         *
         * @param repetition is the repetition to insert
         */
        void push_back(Repetition<REPETITION_TYPE>&& repetition)
        {
            repetitions.push_back(repetition);

            ++stored_repetitions;
            ++total_number;
        }

        /**
         * @brief Get the `n`-th repetition in the storage
         *
         * @param pos is the position of the aimed repetition in the storage
         * @return a constant reference to the `pos`-th repetition in the
         *      storage
         */
        inline const Repetition<REPETITION_TYPE>& operator[](const size_t& pos) const
        {
            return repetitions[pos];
        }

        /**
         * @brief Get the `n`-th repetition in the storage
         *
         * @param pos is the position of the aimed repetition in the storage
         * @return a reference to the `pos`-th repetition in the
         *      storage
         */
        inline Repetition<REPETITION_TYPE>& operator[](const size_t& pos)
        {
            return repetitions[pos];
        }

        /**
         * @brief Get the number of repetitions currently in the storage
         *
         * @return a constant reference to the number of repetitions currently
         *      in the storage
         */
        inline const size_t& size() const
        {
            return stored_repetitions;
        }

        /**
         * @brief Get the initial number of repetitions in the storage
         *
         * @return the initial number of repetitions in the storage
         */
        inline const size_t& full_size() const
        {
            return repetitions.size();
        }

        /**
         * @brief Restore the last extracted repetitions
         *
         * @param num_of_repetitions is the number of the last repetitions to
         *      be restored
         */
        inline void restore_extracted(const size_t num_of_repetitions=1)
        {
            stored_repetitions += num_of_repetitions;
        }

        /**
         * @brief Restore the extracted repetitions
         */
        inline void restore()
        {
            stored_repetitions = full_size();
        }

        /**
         * @brief Save a repetition storage in an archive
         *
         * @tparam ARCHIVE is the output archive type
         * @param archive is the output archive
         */
        template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
        inline void save(ARCHIVE& archive) const
        {
            archive & total_number
                    & stored_repetitions
                    & repetitions;
        }

        /**
         * @brief Load a repetition storage from an archive
         *
         * @tparam ARCHIVE is the input archive type
         * @param archive is the input archive
         * @return the loaded repetition storage
         */
        template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
        inline static RepetitionStorage load(ARCHIVE& archive)
        {
            RepetitionStorage storage;

            archive & storage.total_number
                    & storage.stored_repetitions
                    & storage.repetitions;

            return storage;
        }
    };

    /**
     * @brief The number-of-repetition to repetition-storage map
     */
    using RepetitionMap = std::map<uint8_t, RepetitionStorage>;

    std::map<uint8_t, RepetitionMap> hetero_map;  //!< The heteropolymer map unit-size to repetition-map
    std::map<char, RepetitionMap> homo_map;       //!< Tha homopolymer map nucleotide to repetition-map

    std::mt19937_64 random_gen;     //!< The index random generator

    size_t max_unit_size;           //!< The maximum considered size of the repetition unit
    size_t max_stored_repetitions;  //!< The maximum number of stored repetitions

    /**
     * @brief Initialize the suffix array vector
     *
     * @param s is the sequence whose suffix array is aimed
     * @param suffix_array is the vector that will contain the suffix array
     * @param classes is the vector that labels each sequence position with the
     *      class of its first nucleotide
     * @return the number of the different nucleotides in the sequence
     */
    static size_t init_suffix_array(const char* s, std::vector<ChrPosition>& suffix_array,
                                    std::vector<ChrPosition>& classes)
    {
        const size_t alphabet_size = 1<<8;

        std::vector<ChrPosition> counter(alphabet_size, 0);

        for (size_t i = 0; i < suffix_array.size(); ++i) {
            ++counter[s[i]];
        }
        for (size_t i = 1; i < alphabet_size; ++i) {
            counter[i] += counter[i-1];
        }
        for (size_t i = 1; i <= suffix_array.size(); ++i) {
            auto rev = suffix_array.size()-i;
            suffix_array[--counter[s[rev]]] = rev;
        }
        classes[suffix_array[0]] = 0;
        size_t num_of_classes = 1;
        for (size_t i = 1; i < suffix_array.size(); ++i) {
            if (s[suffix_array[i]] != s[suffix_array[i-1]]) {
                ++num_of_classes;
            }
            classes[suffix_array[i]] = num_of_classes - 1;
        }

        return num_of_classes;
    }

    /**
     * @brief Upgrade the (h-1)-suffix array to a (h)-suffix array
     *
     * @param h is the length of the prefices sorted in the (h)-suffix array
     * @param h_suffix_array is a (h-1)-suffix array
     * @param h_classes is the vector of the classes of the (h-1)-suffix array
     * @param num_of_classes is the number of classes of the (h-1)-suffix array
     * @param tmp_a is a temporary vector having the same size of `h_suffix_array`
     * @param tmp_b is a temporary vector having the same size of `h_suffix_array`
     */
    static void update_suffix_array(const size_t h,
                                    std::vector<ChrPosition>& h_suffix_array,
                                    std::vector<ChrPosition>& h_classes,
                                    size_t& num_of_classes,
                                    std::vector<ChrPosition>& tmp_a,
                                    std::vector<ChrPosition>& tmp_b)
    {
        for (size_t i = 0; i < h_suffix_array.size(); ++i) {
            if (h_suffix_array[i] >= h) {
                tmp_a[i] = h_suffix_array[i] - h;
            } else {
                tmp_a[i] = h_suffix_array[i] + h_suffix_array.size() - h;
            }
        }
        auto& counter = tmp_b;
        fill(counter.begin(), counter.begin() + num_of_classes, 0);
        for (size_t i = 0; i < h_suffix_array.size(); ++i) {
            ++counter[h_classes[tmp_a[i]]];
        }
        for (size_t i = 1; i < num_of_classes; ++i) {
            counter[i] += counter[i-1];
        }
        for (size_t i = 1; i <= h_suffix_array.size(); ++i) {
            const auto& curr = tmp_a[h_suffix_array.size()-i];
            h_suffix_array[--counter[h_classes[curr]]] = curr;
        }

        auto& new_classes = tmp_b;
        new_classes[h_suffix_array[0]] = 0;
        num_of_classes = 1;
        for (size_t i = 1; i < h_suffix_array.size(); i++) {
            const auto& curr = h_suffix_array[i];
            const auto& prev = h_suffix_array[i-1];
            if ((h_classes[curr] != h_classes[prev])
                    || (h_classes[(curr + h) % h_suffix_array.size()] !=
                        h_classes[(prev + h) % h_suffix_array.size()])) {
                ++num_of_classes;
            }
            new_classes[curr] = num_of_classes - 1;
        }
        std::swap(h_classes, new_classes);
    }

    /**
     * @brief Check whether a nucleotide is "A"
     *
     * @param nucleotide is a nucleotide
     * @return `true` is and only if the nucleotide is "A"
     */
    inline static bool is_A(const char& nucleotide)
    {
        return nucleotide=='A' || nucleotide=='a';
    }

    /**
     * @brief Check whether a nucleotide is "T"
     *
     * @param nucleotide is a nucleotide
     * @return `true` is and only if the nucleotide is "T"
     */
    inline static bool is_T(const char& nucleotide)
    {
        return nucleotide=='T' || nucleotide=='t';
    }

    /**
     * @brief Check whether a nucleotide is "C"
     *
     * @param nucleotide is a nucleotide
     * @return `true` is and only if the nucleotide is "C"
     */
    inline static bool is_C(const char& nucleotide)
    {
        return nucleotide=='C' || nucleotide=='c';
    }

    /**
     * @brief Check whether a nucleotide is "G"
     *
     * @param nucleotide is a nucleotide
     * @return `true` is and only if the nucleotide is "G"
     */
    inline static bool is_G(const char& nucleotide)
    {
        return nucleotide=='G' || nucleotide=='g';
    }

    /**
     * @brief Check whether a nucleotide is either "A" or "T"
     *
     * @param nucleotide is a nucleotide
     * @return `true` is and only if the nucleotide is either
     *      "A" or "T"
     */
    inline static bool is_AT(const char& nucleotide)
    {
        return is_A(nucleotide) || is_T(nucleotide);
    }

    /**
     * @brief Check whether a nucleotide is either "C" or "G"
     *
     * @param nucleotide is a nucleotide
     * @return `true` is and only if the nucleotide is either
     *      "C" or "G"
     */
    inline static bool is_CG(const char& nucleotide)
    {
        return is_C(nucleotide) || is_G(nucleotide);
    }

    /**
     * @brief Add a repeated sequence to a repetition map
     *
     * @param r_map is the repetition map in which the repetition should
     *      be added
     * @param chr_id is the identifier of the chromosome containing
     *      the repeated sequence
     * @param begin is the position of the repeated sequence first base
     *      in the chromosome
     * @param num_of_repetitions is the number of repetition of the unit
     * @param unit is the pointer to the unit first base in the sequence
     * @param unit_size is the unit size
     * @param previous_base is base preceding the first repetition of
     *      the unit
     */
    void add(RepetitionMap& r_map, const ChromosomeId& chr_id,
             const ChrPosition& begin,
             const REPETITION_TYPE& num_of_repetitions,
             const char*& unit, const REPETITION_TYPE& unit_size,
             const char& previous_character)
    {
        const uint8_t unit_index = (num_of_repetitions>6?6:num_of_repetitions);

        auto r_it = r_map.find(unit_index);

        if (r_it == r_map.end()) {
            r_it = r_map.emplace(unit_index, RepetitionStorage()).first;
        }

        auto& r_storage = r_it->second;

        if (r_storage.total_number < max_stored_repetitions) {
            r_storage.push_back({chr_id, begin, num_of_repetitions, unit, unit_size,
                                 previous_character});
        } else {
            std::uniform_int_distribution<size_t> dist(0, ++r_storage.total_number);

            auto pos = dist(random_gen);

            if (pos<max_stored_repetitions) {
                r_storage[pos] = Repetition<REPETITION_TYPE>(chr_id, begin,
                                                             num_of_repetitions,
                                                             unit, unit_size,
                                                             previous_character);
            }
        }
    }

    /**
     * @brief Add a repeated sequence to the index
     *
     * @param s is the considered genomic sequence
     * @param chr_id is the identifier of the chromosome containing the
     *      repeated sequence
     * @param begin is the position of the genomic sequence first base
     *      on the chromosome
     * @param unit_size is the size of the repetition unit
     * @param r_begin is the position of the repeated sequence first base
     *      on the considered genomic sequence
     * @param r_end is the position of the repeated sequence last base
     *      on the considered genomic sequence
     * @param covered is the vector of bases in the considered genomic
     *      sequence that belong to a repeated sequence
     */
    void add_repetition(const char* s, const ChromosomeId& chr_id,
                        const ChrPosition& begin, const size_t unit_size,
                        const ChrPosition& r_begin, const ChrPosition& r_end,
                        std::vector<bool>& covered)
    {
        const auto rep_begin = r_begin+begin;
        if (rep_begin>1) {
            auto num_of_repetitions = 1+(r_end-r_begin)/unit_size;
            add(chr_id, rep_begin, static_cast<REPETITION_TYPE>(num_of_repetitions),
                s+r_begin, unit_size, *(s+r_begin-1));

            fill(covered.begin()+r_begin, covered.begin()+r_end+unit_size, true);
        }
    }

    /**
     * @brief Build a random k-mer
     *
     * @param k is the size of the k-mer to build
     * @return a random k-mer
     */
    std::string build_k_mer(const uint8_t k)
    {
        std::uniform_int_distribution<uint8_t> dist(0, 3);

        std::string k_mer(k, 'T');
        for (auto& nucleotide : k_mer) {
            switch(dist(random_gen)) {
                case 0:
                    nucleotide = 'A';
                    break;
                case 1:
                    nucleotide = 'C';
                    break;
                case 2:
                    nucleotide = 'G';
                    break;
                case 3:
                default:
                    //nucleotide = 'T';
                    break;
            }
        }

        return k_mer;
    }

    //! @private
    inline void add_null_heteropolymer(const char* s, const ChromosomeId& chr_id,
                                       const uint8_t unit_size,
                                       const ChrPosition& begin, const ChrPosition& r_begin)
    {
        const auto rep_begin = r_begin+begin;
        if (rep_begin>1) {
            auto& r_storage = hetero_map[unit_size][0];

            if (r_storage.total_number < max_stored_repetitions) {
                std::string k_mer = build_k_mer(unit_size);
                r_storage.push_back(Repetition<REPETITION_TYPE>(chr_id, rep_begin,
                                                                0, k_mer.c_str(),
                                                                unit_size, *(s+r_begin)));
            } else {
                std::uniform_int_distribution<size_t> dist(0, ++r_storage.total_number);

                auto pos = dist(random_gen);

                if (pos<max_stored_repetitions) {
                    std::string k_mer = build_k_mer(unit_size);
                    r_storage[pos] = Repetition<REPETITION_TYPE>(chr_id, rep_begin,
                                                                 0, k_mer.c_str(),
                                                                 unit_size, *(s+r_begin));
                }
            }
        }
    }


    //! @private
    inline void add_null_homopolymer(const size_t nucleotide_index, const char* s,
                                     const ChromosomeId& chr_id, const ChrPosition& begin,
                                     const ChrPosition& r_begin)
    {
        const auto rep_begin = r_begin+begin;
        if (rep_begin>1) {
            add(chr_id, rep_begin, static_cast<uint8_t>(0), s+nucleotide_index, 1,
                *(s+r_begin));
        }
    }

    /**
     * @brief Collect the candidate repeated sequences whose unit size is in [h, 2*h-1]
     *
     * @param begin is the position of the genomic sequence first base
     *      on the chromosome
     * @param h is the order of the suffix array
     * @param h_suffix_array is the (h)-suffix array
     * @param classes is the class vector of the (h)-suffix array
     * @return the candidate repeated sequences. The candidates are encoded in a map
     *      from the position of the candidate first base in the considered sequence 
     *      to a map from the candidate unit size to the position of the final base.
     *      The returned repeated sequence are not added to the index yet because
     *      some of them may be fully included in other candidates
     */
    static std::map<ChrPosition, std::map<size_t, ChrPosition>>
    collect_candidates(const ChrPosition& begin, const size_t& h,
                       std::vector<ChrPosition>& h_suffix_array,
                       std::vector<ChrPosition>& classes)
    {
        ChrPosition next_h = (h>std::numeric_limits<ChrPosition>::max()/2?
                              std::numeric_limits<ChrPosition>::max():2*h);

        std::map<ChrPosition, std::map<size_t, ChrPosition>> candidates;
        ChrPosition r_begin=0, r_end=0, curr_delta = next_h;
        for (size_t i = 1; i < h_suffix_array.size(); ++i) {
            const auto& curr = h_suffix_array[i];
            const auto& prev = h_suffix_array[i-1];
            const auto delta = curr-prev-h;

            if (classes[curr] == classes[prev] && curr >= h + prev
                    && curr < next_h + prev && curr+delta < h_suffix_array.size()
                    && classes[curr+delta] == classes[prev+delta]) {
                if (delta != curr_delta && curr_delta != next_h) {
                    if (r_begin<r_end) {
                        if (begin+r_begin>1) {
                            candidates[r_begin][h+curr_delta] = r_end;
                        }

                        r_begin = curr;
                    }
                }

                curr_delta = delta;
                r_end = curr;
            } else {
                if (r_begin<r_end && begin+r_begin>1) {
                    candidates[r_begin][h+curr_delta] = r_end;
                }

                r_begin = curr;
                r_end = curr;
                curr_delta = next_h;
            }
        }
        if (r_begin<r_end && begin+r_begin>1) {
            candidates[r_begin][h+curr_delta] = r_end;
        }

        return candidates;
    }

    /**
     * @brief Collect the repeated sequences whose unit size is in [h, 2*h-1]
     *
     * @param s is the considered genomic sequence
     * @param chr_id is the identifier of the chromosome containing the
     *      repeated sequence
     * @param begin is the position of the genomic sequence first base
     *      on the chromosome
     * @param h is the order of the suffix array
     * @param h_suffix_array is the (h)-suffix array
     * @param classes is the class vector of the (h)-suffix array
     * @param covered is the vector of bases in the considered genomic
     *      sequence that belong to a repeated sequence
     */
    void collect_repetitions(const char* s, const ChromosomeId& chr_id,
                             const ChrPosition& begin, const size_t& h,
                             std::vector<ChrPosition>& h_suffix_array,
                             std::vector<ChrPosition>& classes,
                             std::vector<bool>& covered)
    {
        auto candidates = collect_candidates(begin, h, h_suffix_array, classes);

        ChrPosition r_begin=0;
        std::map<size_t, ChrPosition> r_endings;
        for (auto c_it=candidates.begin(); c_it != candidates.end(); ++c_it, ++r_begin) {
            for (const auto& [unit_size, r_end]: c_it->second) {
                const auto& r_begin = c_it->first;
                auto e_it = r_endings.find(unit_size);
                if (e_it != r_endings.end()) {
                    if (e_it->second < r_end) {
                        e_it->second = r_end;

                        add_repetition(s, chr_id, begin, unit_size, r_begin, r_end, covered);
                    }
                } else {
                    r_endings.insert({unit_size, r_end});

                    add_repetition(s, chr_id, begin, unit_size, r_begin, r_end, covered);
                }
            }
        }
    }

    /**
     * @brief Analyze non-repeated sequence
     *
     * @param s is the considered genomic sequence
     * @param chr_id is the identifier of the chromosome containing the
     *      repeated sequence
     * @param begin is the position of the genomic sequence first base
     *      on the chromosome
     * @param covered is the vector of bases in the considered genomic
     *      sequence that belong to a repeated sequence
     */
    void analyze_non_repeated_seq(const char* s, const ChromosomeId& chr_id,
                                  const ChrPosition& begin, std::vector<bool>& covered)
    {
        ChrPosition begin_uncovered{0};
        std::vector<ChrPosition> last_char(1<<8, 0);
        for (ChrPosition i=0; i<covered.size(); ++i) {
            if (covered[i]) {
                if (begin_uncovered != i) {
                    for (ChrPosition unit_size=2; unit_size<6; ++unit_size) {
                        for (ChrPosition j=begin_uncovered; j+unit_size<i; ++j) {
                            add_repetition(s, chr_id, begin, unit_size, j, j, covered);
                            add_null_heteropolymer(s, chr_id, unit_size, begin, j);
                        }
                    }
                }
                begin_uncovered = i+1;
            } else {
                if (begin_uncovered == i) {
                    last_char['A'] = last_char['C'] = last_char['G']
                                   = last_char['T'] = i;
                }

                const char& curr_char = *(s+i);
                if (last_char[curr_char]+4<i) {
                    for (ChrPosition j=last_char[curr_char]+2; j<i-2; ++j) {
                        add_null_homopolymer(i, s, chr_id, begin, j);
                    }
                }
                last_char[curr_char] = i;

                add_repetition(s, chr_id, begin, 1, i, i, covered);
            }
        }
    }

    /**
     * @brief Add the repeated sequences of a genomic sequence
     *
     * @param sequence is a genomic sequence
     * @param chr_id is the identifier of the chromosome containing the
     *      genomic sequence
     * @param begin is the position of the genomic sequence first base
     *      on the chromosome
     * @param length is the length of the considered sequence
     */
    void add_repetitions_in(const std::string& sequence, const ChromosomeId& chr_id,
                            const ChrPosition begin, size_t length)
    {
        if (length == 0) {
            return;
        }

        length = std::min(sequence.size()-begin+1, length);
        const char *subseq = sequence.c_str()+begin-1;

        std::vector<ChrPosition> suffix_array(length), classes(length),
                                 tmp_a(length), tmp_b(length);

        size_t num_of_classes = init_suffix_array(subseq, suffix_array, classes);

        std::vector<bool> covered(length, false);

        size_t h;
        size_t next_h;
        const auto h_max = std::min(ceil_div(max_unit_size, static_cast<size_t>(2)), length);
        for (h = 1; h < h_max; h=next_h) {
            next_h = (h>std::numeric_limits<size_t>::max()/2?
                      std::numeric_limits<size_t>::max():2*h);

            collect_repetitions(subseq, chr_id, begin, h, suffix_array, classes, covered);
            update_suffix_array(h, suffix_array, classes, num_of_classes, tmp_a, tmp_b);
        }
        collect_repetitions(subseq, chr_id, begin, h, suffix_array, classes, covered);

        analyze_non_repeated_seq(subseq, chr_id, begin, covered);
    }

    /**
     * @brief Get a polymer reference
     *
     * @tparam INDEX_TYPE is the polymer map index type
     * @param polymer_map is the polymer map (i.e., `hetero_map` or `homo_map`)
     * @param polymer_type is a description of the polymer map
     * @param index is the polymer map index (i.e., unit size or homopolymer nucleotide)
     * @param num_of_repetitions is the number of repetitions of the searched polymer
     * @param for_insertion is a Boolean flag declaring whether the polymer will be
     *      targeted for insertion
     * @return a pair repetition storage-position in the storage refering to a polymer
     *      satisfing the argument specification
     */
    template<typename INDEX_TYPE>
    std::pair<RepetitionStorage&, size_t>
    find_a_polymer(std::map<INDEX_TYPE, RepetitionMap>& polymer_map, const std::string& polymer_type,
                   const INDEX_TYPE index, const uint8_t& num_of_repetitions, const bool& for_insertion)
    {
        const uint8_t repeated_unit = std::min(static_cast<uint8_t>(for_insertion?5:6),
                                               num_of_repetitions);

        auto p_it = polymer_map.find(index);

        std::vector<RepetitionStorage*> r_vectors;

        if (p_it != polymer_map.end()) {
            auto &p_map = p_it->second;
            auto r_it = p_map.find(repeated_unit);
            if (r_it != p_map.end()) {
                r_vectors.push_back(&(r_it->second));
                auto total_size = r_vectors.back()->size();

                if (for_insertion && repeated_unit==5) {
                    r_it = p_map.find(6);

                    if (r_it != p_map.end()) {
                        r_vectors.push_back(&(r_it->second));
                        total_size += r_vectors.back()->size();
                    }
                }

                if (total_size>0) {

                    std::uniform_int_distribution<size_t> dist(0, total_size-1);

                    auto pos = dist(random_gen);

                    for (const auto& r_pointer : r_vectors) {
                        if (pos < r_pointer->size()) {

                            return {*r_pointer, pos};
                        }

                        pos -= r_pointer->size();
                    }
                }
            }
        }
        throw std::runtime_error("No such a " + polymer_type + " available.");
    }

    /**
     * @brief Find a heteropolymer
     *
     * @param unit_size is the unit size of the aimed heteropolymer
     * @param num_of_repetitions is the number of repetitions of the aimed polymer
     * @param for_insertion is a Boolean flag declaring whether the polymer will be
     *      targeted for insertion
     * @return a pair repetition storage-position refering to a heteropolymer having
     *      unit size `unit_size` and whose unit is repeated `num_of_repetitions` times,
     *      when `for_insertion` is `false` or `num_of_repetitions` is smaller than 6,
     *      or either 5 or 6 times, when `for_insertion` is `true` and
     *      `num_of_repetitions` is 5
     */
    inline std::pair<RepetitionStorage&, size_t>
    find_a_heteropolymer(const uint8_t& unit_size, const uint8_t& num_of_repetitions,
                         const bool& for_insertion)
    {
        const uint8_t unit_index = (unit_size>5?5:unit_size);

        return find_a_polymer(hetero_map, "heteropolymer", unit_index, num_of_repetitions,
                              for_insertion);
    }

    /**
     * @brief Find a homopolymer
     *
     * @param unit_nucleotide is the unit nucleotide of the aimed homopolymer
     * @param num_of_repetitions is the number of repetitions of the aimed polymer
     * @param for_insertion is a Boolean flag declaring whether the polymer will be
     *      targeted for insertion or deletion
     * @return a pair repetition storage-position refering to a homopolymer whose
     *      unit is `unit_nucleotide` and it is repeated `num_of_repetitions`
     *      times, when `for_insertion` is `false` or `num_of_repetitions` is
     *      smaller than 6, or either 5 or 6 times, when `for_insertion` is `true`
     *      and `num_of_repetitions is 5
     */
    inline std::pair<RepetitionStorage&, size_t>
    find_a_homopolymer(const char& unit_nucleotide, const uint8_t& num_of_repetitions,
                       const bool& for_insertion)
    {
        return find_a_polymer(homo_map, "homopolymer", canonize(unit_nucleotide),
                              num_of_repetitions, for_insertion);
    }

    /**
     * @brief Canonize a nucleotide
     *
     * @param nucleotide is a nucleotide among 'A', 'C', 'G', and 'T'
     * @return if `nucleotide` is either 'A' or 'T', then 'T'; if
     *      `nucleotide` is either 'A' or 'T', then 'C'
     */
    static char canonize(const char& nucleotide)
    {
        if (is_AT(nucleotide)) {
            return 'T';
        }

        return 'C';
    }

    /**
     * @brief Add a repetition to the index
     *
     * @param chr_id is the identifier of the chromosome containing
     *      the repeated sequence
     * @param begin is the position of the repeated sequence first base
     *      in the chromosome
     * @param num_of_repetitions is the number of repetition of the unit
     * @param unit is the pointer to the unit first base in the sequence
     * @param unit_size is the unit size
     * @param previous_base is base preceding the first repetition of
     *      the unit
     */
    void add(const ChromosomeId& chr_id, const ChrPosition& begin,
             const REPETITION_TYPE& num_of_repetitions,
             const char* unit, const size_t& unit_size,
             const char& previous_base)
    {
        if (unit_size==0) {
            throw std::domain_error("Only initialized repetitions can be added.");
        }

        RepetitionMap* r_map;
        if (unit_size == 1) {
            r_map = &(homo_map[canonize(unit[0])]);
        } else {
            const uint8_t unit_index = (unit_size>5?5:unit_size);

            r_map = &(hetero_map[unit_index]);
        }

        add(*r_map, chr_id, begin, num_of_repetitions, unit, unit_size,
            previous_base);
    }

public:
    /**
     * @brief A constructor
     *
     * @param max_unit_size is the maximum considered size of the repetition unit
     * @param max_stored_repetitions is the maximum number of stored repetitions
     * @param seed is the random generator seed
     */
    RSIndex(const size_t max_unit_size = std::numeric_limits<size_t>::max(),
            const size_t max_stored_repetitions = std::numeric_limits<size_t>::max(),
            const int seed = 0):
        random_gen(seed), max_unit_size(max_unit_size),
        max_stored_repetitions(max_stored_repetitions)
    {}

    /**
     * @brief Add to the index all the repetitions in a chromosome
     *
     * @param chr_id is the identifier of a chrosomome
     * @param chr_sequence is the nucleotide sequence of a chrosomome
     */
    void add_repetitions_in(const ChromosomeId& chr_id, const std::string& chr_sequence)
    {
        size_t begin=1, length=0;
        for (size_t i=0; i<chr_sequence.size(); ++i) {
            if (chr_sequence[i] != 'N' && chr_sequence[i] != 'n') {
                if (length == 0) {
                    begin = i+1;
                }
                ++length;
            } else {
                if (length > 0) {
                    add_repetitions_in(chr_sequence, chr_id,
                                       begin, length);
                    length = 0;
                }
            }
        }
        add_repetitions_in(chr_sequence, chr_id, begin, length);

        std::cout << "Homopolymer" << std::endl;
        for (const auto& [type, repeat_map]: homo_map) {
            std::cout << "  Type: " << std::string(1, type) << std::endl;
            for (const auto& [num_of_repeats, r_storage]: repeat_map) {
                std::cout << "    # " << static_cast<size_t>(num_of_repeats) << ": "
                          << r_storage.total_number << " stored: " 
                          << r_storage.size() << std::endl;
            }
        }

        std::cout << std::endl << "Heteropolymer" << std::endl;
        for (const auto& [unit_size, r_map]: hetero_map) {
            std::cout << "  Unit size: " << static_cast<size_t>(unit_size) << std::endl;

            for (const auto& [num_of_repeats, r_storage]: r_map) {
                std::cout << "    # " << static_cast<size_t>(num_of_repeats) << ": "
                          << r_storage.total_number << " stored: " 
                          << r_storage.size()  << std::endl;
            }
        }
    }

    /**
     * @brief Add a repetition to the index
     *
     * @param repetition is the repetition to be added to the index
     */
    inline void add(const Repetition<REPETITION_TYPE>& repetition)
    {
        add(repetition.begin.chr_id, repetition.begin.position,
            repetition.num_of_repetitions, repetition.unit.c_str(), repetition.unit.size(),
            repetition.prev_base);
    }

    /**
     * @brief Randomly select a heteropolymer
     *
     * @param unit_size is the unit size of the aimed heteropolymer
     * @param num_of_repetitions is the number of repetitions of the aimed polymer
     * @param for_insertion is a Boolean flag declaring whether the polymer will be
     *      targeted for insertion
     * @return a constant reference to a heteropolymer having unit size `unit_size` and
     *      whose unit is repeated `num_of_repetitions` times, when `for_insertion` is
     *      `false` or `num_of_repetitions` is smaller than 6, or either 5 or 6 times,
     *      when `for_insertion` is `true` and `num_of_repetitions` is 6
     */
    inline const Repetition<REPETITION_TYPE>&
    select_a_heteropolymer(const uint8_t& unit_size, const uint8_t& num_of_repetitions,
                           const bool& for_insertion)
    {
        const auto rep_reference = find_a_heteropolymer(unit_size, num_of_repetitions,
                                                        for_insertion);

        return rep_reference.first[rep_reference.second];
    }

    /**
     * @brief Randomly select a heteropolymer and extract it from the index
     *
     * @param unit_size is the unit size of the aimed heteropolymer
     * @param num_of_repetitions is the number of repetitions of the aimed polymer
     * @param for_insertion is a Boolean flag declaring whether the polymer will be
     *      targeted for insertion
     * @return a constant reference to a heteropolymer having unit size `unit_size` and
     *      whose unit is repeated `num_of_repetitions` times, when `for_insertion` is
     *      `false` or `num_of_repetitions` is smaller than 6, or either 5 or 6 times,
     *      when `for_insertion` is `true` and `num_of_repetitions` is 6
     */
    inline const Repetition<REPETITION_TYPE>&
    extract_a_heteropolymer(const uint8_t& unit_size, const uint8_t& num_of_repetitions,
                            const bool& for_insertion)
    {
        auto rep_reference = find_a_heteropolymer(unit_size, num_of_repetitions,
                                                  for_insertion);

        return rep_reference.first.extract(rep_reference.second);
    }

    /**
     * @brief Randomly select a homopolymer
     *
     * @param unit_nucleotide is the unit nucleotide of the aimed homopolymer
     * @param num_of_repetitions is the number of repetitions of the aimed polymer
     * @param for_insertion is a Boolean flag declaring whether the polymer will be
     *      targeted for insertion or deletion
     * @return a constant reference to a homopolymer whose unit is `unit_nucleotide`
     *      and it is repeated `num_of_repetitions` times, when `for_insertion` is
     *      `false` or `num_of_repetitions` is smaller than 6, or either 5 or 6
     *      times, when `for_insertion` is `true` and `num_of_repetitions is 6
     */
    inline const Repetition<REPETITION_TYPE>&
    select_a_homopolymer(const char& unit_nucleotide, const uint8_t& num_of_repetitions,
                         const bool& for_insertion)
    {
        const auto rep_reference = find_a_homopolymer(unit_nucleotide, num_of_repetitions,
                                                      for_insertion);

        return rep_reference.first[rep_reference.second];
    }

    /**
     * @brief Randomly select a homopolymer and extract it from the index
     *
     * @param unit_nucleotide is the unit nucleotide of the aimed homopolymer
     * @param num_of_repetitions is the number of repetitions of the aimed polymer
     * @param for_insertion is a Boolean flag declaring whether the polymer will be
     *      targeted for insertion or deletion
     * @return a constant reference to a homopolymer whose unit is `unit_nucleotide`
     *      and it is repeated `num_of_repetitions` times, when `for_insertion` is
     *      `false` or `num_of_repetitions` is smaller than 6, or either 5 or 6
     *      times, when `for_insertion` is `true` and `num_of_repetitions is 6
     */
    inline const Repetition<REPETITION_TYPE>&
    extract_a_homopolymer(const char& unit_nucleotide, const uint8_t& num_of_repetitions,
                          const bool& for_insertion)
    {
        auto rep_reference = find_a_homopolymer(unit_nucleotide, num_of_repetitions,
                                                for_insertion);

        return rep_reference.first.extract(rep_reference.second);
    }

    /**
     * @brief Restore the last extracted heteropolymer
     *
     * @param unit_size is the unit size of the aimed heteropolymer
     * @param num_of_repetitions is the number of repetitions of the aimed heteropolymer
     */
    inline void restore_heteropolymer(const uint8_t& unit_size, const uint8_t& num_of_repetitions)
    {
        hetero_map.at((unit_size>5?5:unit_size)).at(num_of_repetitions).restore_extracted();
    }

    /**
     * @brief Restore the last extracted homopolymer
     *
     * @param unit_nucleotide is the unit nucleotide of the aimed homopolymer
     * @param num_of_repetitions is the number of repetitions of the aimed polymer
     */
    inline void restore_homopolymer(const char& unit_nucleotide, const uint8_t& num_of_repetitions)
    {
        homo_map.at(canonize(unit_nucleotide)).at(num_of_repetitions).restore_extracted();
    }

    /**
     * @brief Restore the extracted polymers
     */
    void restore()
    {
        for (const auto& [unit_nucleotipe, r_map]: homo_map) {
            for (const auto& [num_of_repeats, r_storage]: r_map) {
                r_map.restore();
            }
        }

        std::cout << std::endl << "Heteropolymer" << std::endl;
        for (const auto& [unit_size, r_map]: hetero_map) {
            for (const auto& [num_of_repeats, r_storage]: r_map) {
                r_map.restore();
            }
        }
    }

    /**
     * @brief Save a repeated sequence index in an archive
     *
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        archive & hetero_map
                & homo_map
                & max_unit_size
                & max_stored_repetitions;
    }

    /**
     * @brief Load a repeated sequence index from an archive
     *
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the loaded repeated sequence index
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    inline static RSIndex<REPETITION_TYPE> load(ARCHIVE& archive)
    {
        RSIndex<REPETITION_TYPE> rs_index;

        archive & rs_index.hetero_map
                & rs_index.homo_map
                & rs_index.max_unit_size
                & rs_index.max_stored_repetitions;

        return rs_index;
    }
};

}   // Mutations

}   // Races

#endif // __RACES_RS_INDEX__