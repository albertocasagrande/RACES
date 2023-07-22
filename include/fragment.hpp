/**
 * @file fragment.hpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Defines genomic fragment
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

#ifndef __RACES_FRAGMENT__
#define __RACES_FRAGMENT__

#include <map>
#include <vector>
#include <ostream>

#include "snv.hpp"

namespace Races 
{

namespace Passengers
{

/**
 * @brief A genomic fragment
 * 
 * This class represents a genomic fragment, its 
 * alleles, and all the SNV that occur on it.
 */
struct Fragment
{
    using Length = ChrPosition;

    /**
     * @brief A type for alleles
     * 
     * An allele is the set of SNVs occurring 
     * on the allele itself.
     */
    typedef std::map<GenomicPosition, SNV> Allele;

protected:
    GenomicPosition initial_pos;    //!< the fragment initial position
    Length length;                  //!< the fragment length

    std::vector<Allele> alleles;    //!< the fragment alleles
public:
    /**
     * @brief The empty constructor
     * 
     */
    Fragment();

    /**
     * @brief A constructor
     * 
     * @param chromosome_id is the identifier of the chromosome on which the fragment lays
     * @param length is the length of the fragment
     * @throw std::domain_error `length` is 0
     */
    Fragment(const ChromosomeId chromosome_id, const Length length);

    /**
     * @brief A constructor
     * 
     * @param chromosome_id is the identifier of the chromosome on which the fragment lays
     * @param length is the length of the fragment
     * @param alleles is the allele vector
     * @throw std::domain_error `length` is 0
     */
    Fragment(const ChromosomeId chromosome_id, const Length length, const std::vector<Allele>& alleles);

    /**
     * @brief A constructor
     * 
     * @param initial_pos is the fragment initial position
     * @param length is the fragment length
     * @throw std::domain_error `length` is 0
     */
    Fragment(const GenomicPosition initial_pos, const Length length);

    /**
     * @brief A constructor
     * 
     * @param initial_pos is the fragment initial position
     * @param length is the fragment length
     * @param alleles is the allele vector 
     * @throw std::domain_error `length` is 0
     */
    Fragment(const GenomicPosition initial_pos, const Length length, const std::vector<Allele>& alleles);

    /**
     * @brief Get the identifier of the chromosome containing the fragment
     * 
     * @return a constant reference to the identifier of the chromosome 
     *          containing the fragment
     */
    inline const ChromosomeId& get_chromosome_id() const
    {
        return initial_pos.chr_id;
    }

    /**
     * @brief Get the fragment initial position
     * 
     * @return a constant reference to the fragment initial position
     */
    inline const GenomicPosition& get_begin() const
    {
        return initial_pos;
    }

    /**
     * @brief Get the fragment final position
     * 
     * @return the fragment final position
     */
    inline GenomicPosition get_end() const
    {
        return GenomicPosition(initial_pos.chr_id, initial_pos.position+length-1);
    }

    /**
     * @brief Get the fragment length
     * 
     * @return a constant reference to the fragment length
     */
    inline const Length& size() const
    {
        return length;
    }

    /**
     * @brief Get the fragment allele vector
     * 
     * @return a constant reference to the allele vector
     */
    inline const std::vector<Allele>& get_alleles() const
    {
        return alleles;
    }

    /**
     * @brief Get the number of fragment alleles
     * 
     * @return the number of fragment alleles
     */
    inline size_t num_of_alleles() const
    {
        return alleles.size();
    }

    /**
     * @brief Get the i-th allele
     * 
     * @param index is the aimed allele index
     * @return a constant reference to the `index`-th allele
     *      in the fragment
     */
    inline const Allele& operator[](const size_t index) const
    {
        return alleles.at(index);
    }

    /**
     * @brief Test whether a fragment follows another one
     * 
     * @param fragment is a fragment
     * @return `true` if and only if the current fragment begins just 
     *          after the end of `fragment`
     */
    inline bool follows(const Fragment& fragment) const
    {
        return ((initial_pos.chr_id == fragment.initial_pos.chr_id)
                && (fragment.initial_pos.position + fragment.length == initial_pos.position));
    }

    /**
     * @brief Test whether a fragment follows another one
     * 
     * @param fragment is a fragment
     * @return `true` if and only if the current fragment begins just 
     *          after the end of `fragment`
     */
    inline bool follows(Fragment&& fragment) const
    {
        return follows(fragment);
    }

    /**
     * @brief Test whether a fragment precedes another one
     * 
     * @param fragment is a fragment
     * @return `true` if and only if the current `fragment` begins just 
     *          after the end of the current fragment
     */
    inline bool precedes(const Fragment& fragment) const
    {
        return fragment.follows(*this);
    }

    /**
     * @brief Test whether a fragment precedes another one
     * 
     * @param fragment is a fragment
     * @return `true` if and only if the current `fragment` begins just 
     *          after the end of the current fragment
     */
    inline bool precedes(Fragment&& fragment) const
    {
        return fragment.follows(*this);
    }

    /**
     * @brief Check whether a position is contained by the fragment
     * 
     * @param genomic_position is the genomic position to check
     * @return `true` if and only if `position` is contained by the
     *      fragment
     */
    inline bool contains(const GenomicPosition& genomic_position) const
    {
        return ((genomic_position.chr_id==initial_pos.chr_id)
                && (initial_pos.position <= genomic_position.position)
                && (genomic_position.position+1 <= initial_pos.position+length));
    }

    /**
     * @brief Check whether a position is contained by the fragment
     * 
     * @param genomic_position is the genomic position to check
     * @return `true` if and only if `position` is contained by the
     *      fragment
     */
    inline bool contains(GenomicPosition&& genomic_position) const
    {
        return contains(genomic_position);
    }

    /**
     * @brief Check whether a position is strictly contained by the fragment
     * 
     * @param genomic_position is the genomic position to check
     * @return `true` if and only if `position` is strictly contained by the
     *      fragment
     */
    inline bool strictly_contains(const GenomicPosition& genomic_position) const
    {
        return ((genomic_position.chr_id==initial_pos.chr_id)
                && (initial_pos.position < genomic_position.position)
                && (genomic_position.position+1 < initial_pos.position+length));
    }

    /**
     * @brief Check whether a position is strictly contained by the fragment
     * 
     * @param genomic_position is the genomic position to check
     * @return `true` if and only if `position` is strictly contained by the
     *      fragment
     */
    inline bool strictly_contains(GenomicPosition&& genomic_position) const
    {
        return strictly_contains(genomic_position);
    }

    /**
     * @brief Split a fragment
     * 
     * This method cuts a fragment in a position, generates a 
     * new fragment beginning at the specified position, and 
     * updates the length of the considered fragment so that 
     * the two fragments are contiguous, i.e., one of the two
     * follows the other one. 
     * The number of alleles the new fragment is the same of 
     * the old one. The SNVs of the original fragment are 
     * distributed on the two fragments according to their 
     * genomic position; their allele indices do not change.
     * 
     * @param split_point is the position of the new fragment
     * @return the fragment originated by the split
     * @throw std::domain_error the fragment does not contain 
     *          `split_point` or `split_point` and the fragment
     *          initial point are the same
     */
    Fragment split(const GenomicPosition& split_point);

    /**
     * @brief Split a fragment
     * 
     * This method cuts a fragment in a position, generates a 
     * new fragment beginning at the specified position, and 
     * updates the length of the considered fragment so that 
     * the two fragments are contiguous, i.e., one of the two
     * follows the other one. 
     * The number of alleles the new fragment is the same of 
     * the old one. The SNVs of the original fragment are 
     * distributed on the two fragments according to their 
     * genomic position; their allele indices do not change.
     * 
     * @param split_point is the position of the new fragment
     * @return the fragment originated by the split
     * @throw std::domain_error the fragment does not contain 
     *          `split_point` or `split_point` and the fragment
     *          initial point are the same
     */
    inline Fragment split(GenomicPosition&& split_point)
    {
        return split(split_point);
    }

    /**
     * @brief Join two contiguous `n`-ploid fragments
     * 
     * Two fragments are contiguous when the inital position 
     * of one of the two is the initial position of the other 
     * one plus the latter's length plus one. 
     * This method joins the current fragment and a contiguous 
     * one. It sets the initial position of the current object 
     * to the first of the two initial positions, and the its 
     * length to the sum of the two lengths. The alleles of 
     * the resulting fragment are obtained by considering the 
     * SNVs in the corresponding alleles of the original 
     * fragments.
     * If the two fragment differ in the number of alleles or 
     * they are not contiguous, a `std::domain_error` is thrown.
     * 
     * @param contiguous_fragment is a fragment contiguous to 
     *          the current one
     * @return a reference to the updated fragment
     * @throw std::domain_error the current fragment and 
     *      `contiguous_fragment` differ in the number of 
     *      alleles or they are not contiguous.
     */
    Fragment& join(Fragment& contiguous_fragment);

    /**
     * @brief Duplicate an allele
     * 
     * @param index is the index of the allele to be duplicated
     * @return a reference to the updated fragment
     * @throw std::out_of_range `allele_index` is not a valid allele index
     */
    Fragment& duplicate_allele(const size_t index);

    /**
     * @brief Remove an allele
     * 
     * @param index is the index of the allele to be removed
     * @return a reference to the updated fragment
     * @throw std::out_of_range `allele_index` is not a valid allele index
     */
    Fragment& remove_allele(const size_t index);

    /**
     * @brief Check whether any SNV occurs a possible mutational context
     * 
     * @param genomic_position is the central position of the mutational context
     * @return `true` if and only if no SNVs 
     * @throw std::domain_error `position` does not lays in the fragment
     */
    bool is_mutational_context_free(const GenomicPosition& genomic_position) const;

    /**
     * @brief Check whether any SNV occurs a possible mutational context
     * 
     * @param genomic_position is the central position of the mutational context
     * @return `true` if and only if no SNVs 
     * @throw std::domain_error `position` does not lays in the fragment
     */
    inline bool is_mutational_context_free(GenomicPosition&& genomic_position) const
    {
        return is_mutational_context_free(genomic_position);
    }

    /**
     * @brief Insert a SNV in one of the fragment alleles
     * 
     * This method tries to insert a SNV in one of the fragment alleles. 
     * If any of fragment alleles contains another SNV occurring in 
     * the mutational context of the specified SNV, then the insertion 
     * fails.
     * 
     * @param snv is the SNV to be inserted
     * @param allele_index is the index of the allele in which inserting 
     *      the SNV
     * @return `true` if and only if the fragment alleles do not 
     *      contain any other SNVs in `snv`'s mutational context
     * @throw std::domain_error `SNV` does not lays in the fragment
     * @throw std::out_of_range `allele_index` is not a valid allele index
     */
    inline bool insert(const SNV& snv, const size_t allele_index)
    {
        return insert(SNV(snv), allele_index);
    }

    /**
     * @brief Insert a SNV in one of the fragment alleles
     * 
     * This method tries to insert a SNV in one of the fragment alleles. 
     * If any of fragment alleles contains another SNV occurring in 
     * the mutational context of the specified SNV, then the insertion 
     * fails.
     * 
     * @param snv is the SNV to be inserted
     * @param allele_index is the index of the allele in which inserting 
     *      the SNV
     * @return `true` if and only if the fragment alleles do not 
     *      contain any other SNVs in `snv`'s mutational context
     * @throw std::domain_error `SNV` does not lays in the fragment
     * @throw std::out_of_range `allele_index` is not a valid allele index
     */
    bool insert(SNV&& snv, const size_t allele_index);

    /**
     * @brief Remove a SNV from the fragment
     * 
     * This method tries to remove a SNV from the fragment.
     * 
     * @param genomic_position is the position of the SNV to be removed
     * @return `true` if and only if a SNV occurred at `position`
     * @throw std::domain_error `position` does not lays in the fragment
     */
    bool remove_SNV(const GenomicPosition& genomic_position);

    /**
     * @brief Remove a SNV from the fragment
     * 
     * This method tries to remove a SNV from the fragment.
     * 
     * @param genomic_position is the position of the SNV to be removed
     * @return `true` if and only if a SNV occurred at `position`
     * @throw std::domain_error `position` does not lays in the fragment
     */
    inline bool remove_SNV(GenomicPosition&& genomic_position)
    {
        return remove_SNV(genomic_position);
    }
};


}   // Passengers

}   // Races

namespace std
{

/**
 * @brief Write a fragment in a output stream
 * 
 * @param out is the output stream
 * @param fragment is the fragment to stream
 * @return a reference of the updated stream
 */
std::ostream& operator<<(std::ostream& out, const Races::Passengers::Fragment& fragment);

}

#endif // __RACES_FRAGMENT__
