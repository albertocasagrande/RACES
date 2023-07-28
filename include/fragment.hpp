/**
 * @file fragment.hpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Defines genomic fragment
 * @version 0.6
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

#ifndef __RACES_FRAGMENT__
#define __RACES_FRAGMENT__

#include <map>
#include <set>
#include <vector>
#include <ostream>

#include "snv.hpp"

namespace Races 
{

namespace Passengers
{

/**
 * @brief A genomic region is a interval of genomic positions
 * 
 */
struct GenomicRegion
{
    using Length = ChrPosition;
protected:
    GenomicPosition initial_pos;    //!< the genomic region initial position
    Length length;                  //!< the genomic region length
public:
    /**
     * @brief The empty constructor
     */
    GenomicRegion();

    /**
     * @brief A constructor
     * 
     * @param chromosome_id is the identifier of the chromosome on which the region lays
     * @param length is the length of the genomic region
     * @throw std::domain_error `length` is 0
     */
    GenomicRegion(const ChromosomeId chromosome_id, const Length length);

    /**
     * @brief A constructor
     * 
     * @param initial_pos is the genomic region initial position
     * @param length is the genomic region length
     * @throw std::domain_error `length` is 0
     */
    GenomicRegion(const GenomicPosition initial_pos, const Length length);

    /**
     * @brief Get the identifier of the chromosome containing the genomic region
     * 
     * @return a constant reference to the identifier of the chromosome 
     *          containing the genomic region
     */
    inline const ChromosomeId& get_chromosome_id() const
    {
        return initial_pos.chr_id;
    }

    /**
     * @brief Get the genomic region initial position
     * 
     * @return a constant reference to the genomic region initial position
     */
    inline const GenomicPosition& get_begin() const
    {
        return initial_pos;
    }

    /**
     * @brief Get the genomic region final position
     * 
     * @return the genomic region final position
     */
    inline GenomicPosition get_end() const
    {
        return GenomicPosition(initial_pos.chr_id, initial_pos.position+length-1);
    }

    /**
     * @brief Get the genomic region initial position in the chromosome
     * 
     * @return the genomic region initial position in the chromosome
     */
    inline const ChrPosition& get_initial_position() const
    {
        return initial_pos.position;
    }

    /**
     * @brief Get the genomic region final position in the chromosome
     * 
     * @return the genomic region final position in the chromosome
     */
    inline ChrPosition get_final_position() const
    {
        return initial_pos.position+length-1;
    }

    /**
     * @brief Get the genomic region length
     * 
     * @return a constant reference to the genomic region length
     */
    inline const Length& size() const
    {
        return length;
    }

    /**
     * @brief Test whether a region follows another region
     * 
     * Let `A` and `B` two regions . "`A` follows `B`" when the initial position 
     * of `A` immediately follows the final position of `B`.  
     * 
     * @param genomic_region is a genomic region
     * @return `true` if and only if the current genomic region begins just 
     *          after the end of `genomic_region`
     */
    inline bool follows(const GenomicRegion& genomic_region) const
    {
        return ((initial_pos.chr_id == genomic_region.initial_pos.chr_id)
                && (genomic_region.initial_pos.position + genomic_region.length == initial_pos.position));
    }

    /**
     * @brief Test whether a region follows another region
     * 
     * Let `A` and `B` two regions . "`A` follows `B`" when the initial position 
     * of `A` immediately follows the final position of `B`.  
     * 
     * @param genomic_region is a genomic region
     * @return `true` if and only if the current genomic region begins just 
     *          after the end of `genomic_region`
     */
    inline bool follows(GenomicRegion&& genomic_region) const
    {
        return follows(genomic_region);
    }

    /**
     * @brief Test whether a region precedes another region
     * 
     * Let `A` and `B` two regions . "`A` precedes `B`" when the final position 
     * of `A` immediately follows the initial position of `B`.  
     * 
     * @param genomic_region is a genomic region
     * @return `true` if and only if the current `genomic_region` ends just 
     *          before the beginning of the current genomic region
     */
    inline bool precedes(const GenomicRegion& genomic_region) const
    {
        return genomic_region.follows(*this);
    }

    /**
     * @brief Test whether a region precedes another region
     * 
     * Let `A` and `B` two regions . "`A` precedes `B`" when the final position 
     * of `A` immediately follows the initial position of `B`.  
     * 
     * @param genomic_region is a genomic region
     * @return `true` if and only if the current `genomic_region` ends just 
     *          before the beginning of the current genomic region
     */
    inline bool precedes(GenomicRegion&& genomic_region) const
    {
        return genomic_region.follows(*this);
    }

    /**
     * @brief Check whether two genomic region overlap
     * 
     * @param genomic_region is a genomic region
     * @return `true` if and only if the current region and 
     *       `genomic_region` 
     */
    bool overlaps(const GenomicRegion& genomic_region) const;

    /**
     * @brief Check whether two genomic region overlap
     * 
     * @param genomic_region is a genomic region
     * @return `true` if and only if the current region and 
     *       `genomic_region` 
     */
    inline bool overlaps(GenomicRegion&& genomic_region) const
    {
        return this->overlaps(genomic_region);
    }

    /**
     * @brief Check whether a region ends before another one begin
     * 
     * @param genomic_region is a genomic region
     * @return `true` if and only if the current region ends before 
     *          `genomic_region`'s initial position
     */
    inline bool ends_before(const GenomicRegion& genomic_position) const
    {
        return ends_before(genomic_position.get_begin());
    }

    /**
     * @brief Check whether a region ends before another one begin
     * 
     * @param genomic_region is a genomic region
     * @return `true` if and only if the current region ends before 
     *          `genomic_region`'s initial position
     */
    inline bool ends_before(GenomicRegion&& genomic_position) const
    {
        return ends_before(genomic_position);
    }

    /**
     * @brief Check whether a region begin after another one end
     * 
     * @param genomic_region is a genomic region
     * @return `true` if and only if the current region ends before 
     *          `genomic_position`
     */
    inline bool begins_after(const GenomicRegion& genomic_position) const
    {
        return genomic_position.ends_before(*this);
    }

    /**
     * @brief Check whether a position is contained in the genomic region
     * 
     * @param genomic_position is the genomic position to check
     * @return `true` if and only if `position` is contained in the
     *      genomic region
     */
    inline bool contains(const GenomicPosition& genomic_position) const
    {
        return ((genomic_position.chr_id==initial_pos.chr_id)
                && (initial_pos.position <= genomic_position.position)
                && (genomic_position.position+1 <= initial_pos.position+length));
    }

    /**
     * @brief Check whether a position is contained in the genomic region
     * 
     * @param genomic_position is the genomic position to check
     * @return `true` if and only if `position` is contained in the
     *      genomic region
     */
    inline bool contains(GenomicPosition&& genomic_position) const
    {
        return contains(genomic_position);
    }

    /**
     * @brief Check whether a position is strictly contained in the genomic region
     * 
     * @param genomic_position is the genomic position to check
     * @return `true` if and only if `position` is strictly contained in the
     *      genomic region
     */
    inline bool strictly_contains(const GenomicPosition& genomic_position) const
    {
        return ((genomic_position.chr_id==initial_pos.chr_id)
                && (initial_pos.position < genomic_position.position)
                && (genomic_position.position+1 < initial_pos.position+length));
    }

    /**
     * @brief Check whether a position is strictly contained in the genomic region
     * 
     * @param genomic_position is the genomic position to check
     * @return `true` if and only if `position` is strictly contained in the
     *      genomic region
     */
    inline bool strictly_contains(GenomicPosition&& genomic_position) const
    {
        return strictly_contains(genomic_position);
    }

    /**
     * @brief Check whether a region ends before a position
     * 
     * @param genomic_position is a genomic position
     * @return `true` if and only if the current region ends before 
     *          `genomic_position`
     */
    inline bool ends_before(const GenomicPosition& genomic_position) const
    {
        return ((genomic_position.chr_id==initial_pos.chr_id)
                && (initial_pos.position+length < genomic_position.position+1));
    }

    /**
     * @brief Check whether a region ends before a position
     * 
     * @param genomic_position is a genomic position
     * @return `true` if and only if the current region ends before 
     *          `genomic_position`
     */
    inline bool ends_before(GenomicPosition&& genomic_position) const
    {
        return ends_before(genomic_position);
    }

    /**
     * @brief Check whether a region begins after a position
     * 
     * @param genomic_position is a genomic position
     * @return `true` if and only if the current region begins after 
     *          `genomic_position`
     */
    inline bool begins_after(const GenomicPosition& genomic_position) const
    {
        return ((genomic_position.chr_id==initial_pos.chr_id)
                && (initial_pos.position > genomic_position.position));
    }

    /**
     * @brief Check whether a region begins after a position
     * 
     * @param genomic_position is a genomic position
     * @return `true` if and only if the current region begins after 
     *          `genomic_position`
     */
    inline bool begins_after(GenomicPosition&& genomic_position) const
    {
        return begins_after(genomic_position);
    }

    /**
     * @brief Split a genomic region
     * 
     * This method cuts a genomic region in a position, generates a 
     * new genomic region beginning at the specified position, and 
     * updates the length of the considered genomic region so that 
     * the two fragments are contiguous, i.e., one of the two
     * follows the other one.
     * 
     * @param split_point is the position of the new genomic region
     * @return the fragment originated by the split
     * @throw std::domain_error the genomic region does not contain 
     *          `split_point` or `split_point` and the genomic region
     *          initial point are the same
     */
    GenomicRegion split(const GenomicPosition& split_point);

    /**
     * @brief Split a genomic region
     * 
     * This method cuts a genomic region in a position, generates a 
     * new genomic region beginning at the specified position, and 
     * updates the length of the considered genomic region so that 
     * the two fragments are contiguous, i.e., one of the two
     * follows the other one.
     * 
     * @param split_point is the position of the new genomic region
     * @return the fragment originated by the split
     * @throw std::domain_error the genomic region does not contain 
     *          `split_point` or `split_point` and the genomic region
     *          initial point are the same
     */
    inline GenomicRegion split(GenomicPosition&& split_point)
    {
        return split(split_point);
    }

    /**
     * @brief Join two contiguous genomic regions
     * 
     * Two genomic regions are contiguous when the inital position 
     * of one of the two is the initial position of the other 
     * one plus the latter's length plus one. 
     * This method joins the current genomic region and a contiguous 
     * one. It sets the initial position of the current object 
     * to the first of the two initial positions, and the its 
     * length to the sum of the two lengths.
     * 
     * @param contiguous_region is a genomic region contiguous to 
     *          the current one
     * @return a reference to the updated genomic region
     */
    GenomicRegion& join(GenomicRegion& contiguous_region);
};

/**
 * @brief A type for alleles
 * 
 * An allele is the set of SNVs occurring 
 * on the allele itself.
 */
typedef std::map<GenomicPosition, SNV> Allele;


/**
 * @brief A identifier type for alleles
 */
typedef size_t AlleleId;

/**
 * @brief A genomic fragment
 * 
 * This class represents a genomic fragment, its 
 * alleles, and all the SNV that occur on it.
 */
struct Fragment : public GenomicRegion
{
    using Length = ChrPosition;

protected:
    std::map<AlleleId, Allele> alleles;    //!< the fragment alleles
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
     * @param num_of_alleles is the initial number of alleles
     * @throw std::domain_error `length` is 0
     */
    Fragment(const GenomicRegion& genomic_region, const size_t& num_of_alleles);    

    /**
     * @brief A constructor
     * 
     * @param chromosome_id is the identifier of the chromosome on which the fragment lays
     * @param length is the length of the fragment
     * @param num_of_alleles is the initial number of alleles
     * @throw std::domain_error `length` is 0
     */
    Fragment(const ChromosomeId chromosome_id, const Length length, const size_t& num_of_alleles);

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
     * @param num_of_alleles is the initial number of alleles
     * @throw std::domain_error `length` is 0
     */
    Fragment(const GenomicPosition initial_pos, const Length length, const size_t& num_of_alleles);

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
     * @brief Get the fragment allele map
     * 
     * @return a constant reference to the allele map
     */
    inline const std::map<AlleleId, Allele>& get_alleles() const
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
     * @brief Check whether two fragments have the same allele ids
     * 
     * @param fragment is a fragment
     * @return `true` if and only if `fragment` and the current 
     *       object have the same allele ids.
     */
    bool has_the_same_allele_ids(const Fragment& fragment) const;

    /**
     * @brief Get the identificators of the alleles
     * 
     * @return the set of the allele identificators available 
     *      in this fragment
     */
    std::set<AlleleId> get_allele_ids() const;

    /**
     * @brief Check whether the fragment contains an allele
     * 
     * @param allele_id is an allele identifier
     * @return `true` if and only if the fragment contains 
     *      the the allele having `allele_id` as id.
     */
    inline bool has_allele(const AlleleId& allele_id) const
    {
        return alleles.count(allele_id)>0;
    }

    /**
     * @brief Get the i-th allele
     * 
     * @param allele_id is the aimed allele identifier
     * @return a constant reference to the aimed allele 
     *      in the fragment
     */
    inline const Allele& operator[](const AlleleId& allele_id) const
    {
        return alleles.at(allele_id);
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
     * @param allele_id is the identifier of the allele to be duplicated
     * @param new_allele_id is the identifier of the new allele
     * @return a reference to the updated fragment
     * @throw std::out_of_range `allele_id` is not a valid allele index
     * @throw std::domain_error `new_allele_id` is already present in the
     *          fragment
     */
    Fragment& duplicate_allele(const AlleleId& allele_id, const AlleleId& new_allele_id);

    /**
     * @brief Remove an allele
     * 
     * @param allele_id is the identifier of the allele to be removed
     * @return a reference to the updated fragment
     * @throw std::out_of_range `allele_id` is not a valid allele index
     */
    Fragment& remove_allele(const AlleleId& allele_id);

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
     * @param allele_id is the identifier of the allele in which 
     *      inserting the SNV
     * @return `true` if and only if the fragment alleles do not 
     *      contain any other SNVs in `snv`'s mutational context
     * @throw std::domain_error `SNV` does not lays in the fragment
     * @throw std::out_of_range `allele_id` is not a valid allele
     *          identifier
     */
    inline bool insert(const SNV& snv, const AlleleId allele_id)
    {
        return insert(SNV(snv), allele_id);
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
     * @param allele_id is the identifier of the allele in which 
     *      inserting the SNV
     * @return `true` if and only if the fragment alleles do not 
     *      contain any other SNVs in `snv`'s mutational context
     * @throw std::domain_error `SNV` does not lays in the fragment
     * @throw std::out_of_range `allele_id` is not a valid allele
     *          identifier
     */
    bool insert(SNV&& snv, const AlleleId allele_index);

    /**
     * @brief Remove a SNV from the fragment
     * 
     * This method tries to remove a SNV from the fragment.
     * 
     * @param genomic_position is the position of the SNV to be removed
     * @return `true` if and only if a SNV occurred at `genomic_position`
     * @throw std::domain_error `genomic_position` does not lays in the
     *      fragment
     */
    bool remove_SNV(const GenomicPosition& genomic_position);

    /**
     * @brief Remove a SNV from the fragment
     * 
     * This method tries to remove a SNV from the fragment.
     * 
     * @param genomic_position is the position of the SNV to be removed
     * @return `true` if and only if a SNV occurred at `genomic_position`
     * @throw std::domain_error `genomic_position` does not lays in the
     *      fragment
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

template<>
struct less<Races::Passengers::GenomicRegion>
{
    inline bool operator()(const Races::Passengers::GenomicRegion &lhs,
                           const Races::Passengers::GenomicRegion &rhs) const
    {
        return ((lhs.get_chromosome_id()<rhs.get_chromosome_id())
                || ((lhs.get_chromosome_id()==rhs.get_chromosome_id()) 
                    && (lhs.get_initial_position()<rhs.get_initial_position()))
                || ((lhs.get_chromosome_id()==rhs.get_chromosome_id()) 
                    && (lhs.get_initial_position()==rhs.get_initial_position()) 
                    && (lhs.size()<rhs.size())));
    }
};

/**
 * @brief Write a genomic region in a output stream
 * 
 * @param out is the output stream
 * @param genomic_region is the genomic region to stream
 * @return a reference of the updated stream
 */
std::ostream& operator<<(std::ostream& out, const Races::Passengers::GenomicRegion& genomic_region);

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
