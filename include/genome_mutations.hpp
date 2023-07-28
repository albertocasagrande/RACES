/**
 * @file genome_mutations.hpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Defines genome and chromosome data structures
 * @version 0.1
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

#ifndef __RACES_GENOME_MUTATIONS__
#define __RACES_GENOME_MUTATIONS__

#include <map>
#include <fstream>

#include "snv_signature.hpp"
#include "context_positions.hpp"

#include "fragment.hpp"
#include "progress_bar.hpp"

namespace Races 
{

namespace Passengers
{

/**
 * @brief A class to represent the passenger mutations in a chromosome
 */
class ChromosomeMutations
{
    using Length = GenomicRegion::Length;
private:
    ChromosomeId identifier;    //!< The chromosome identifier
    Length length;              //!< The chromosome length
    Length allelic_length;      //!< The sum of the lengths of all the alleles

    std::map<GenomicPosition, Fragment> fragments;  //!< The chromosome fragments

    AlleleId next_allele_id;   //!< the identifier of the next allele

public:
    /**
     * @brief The empty constructor
     */
    ChromosomeMutations();

    /**
     * @brief A constructor
     * 
     * @param identifier is the chromosome identifier
     * @param size is the chromosome size
     * @param num_of_alleles is the initial number of alleles
     */
    ChromosomeMutations(const ChromosomeId& identifier, const ChromosomeMutations::Length& size, const size_t& num_of_alleles);

    /**
     * @brief A constructor
     * 
     * @param chromosome_region is the chromosome region
     * @param num_of_alleles is the initial number of alleles
     */
    ChromosomeMutations(const GenomicRegion& chromosome_region, const size_t& num_of_alleles);

    /**
     * @brief Get chromosome identifier
     * 
     * @return the chromosome identifier
     */
    inline ChromosomeId id() const
    {
        return identifier;
    }

    /**
     * @brief Get chromosome size
     * 
     * @return the chromosome size
     */
    inline ChromosomeMutations::Length size() const
    {
        return length;
    }

    /**
     * @brief Get chromosome allelic size
     * 
     * @return the chromosome allelic size
     */
    inline ChromosomeMutations::Length allelic_size() const
    {
        return allelic_length;
    }

    /**
     * @brief Get the chromosome fragments
     * 
     * @return the chromosome fragments
     */
    inline const std::map<GenomicPosition, Fragment>& get_fragments() const
    {
        return fragments;
    } 

    /**
     * @brief Get the identifiers of the alleles that completely include a region
     * 
     * @param genomic_region is a region of the chromosome
     * @return the identifiers of the allele that completely include `genomic_region`
     */
    std::set<AlleleId> get_allele_ids_in(const GenomicRegion& genomic_region) const;

    /**
     * @brief Get the sets of SNVs for all the alleles
     * 
     * This methods returns the SNVs occurring in each allele on the chromosome even 
     * if the allele does not cover the full chromosome.
     * 
     * @return a map that associates the allele id to the SNVs present on that allele
     */
    std::map<AlleleId, std::set<SNV>> get_SNVs_per_allele() const;

    /**
     * @brief Check whether a chromosomic region has an allele 
     * 
     * @param allele_id is the identifier of the allele whose presence must be checked
     * @param genomic_region is the region on which the check is performed 
     * @return `true` if and only if the whole `genomic_region` contains the allele
     *          having `allele_id` as allele identifier
     * @throw std::domain_error `genomic_region` does not lays in this chromosome
     */
    bool has_allele_on(const AlleleId& allele_id, const GenomicRegion& genomic_region) const;

    /**
     * @brief Check whether a genomic position lays in a chromosome
     * 
     * @param genomic_position is a genomic position
     * @return `true` if and only if `genomic_position` lays in the current chromosome
     */
    bool contains(const GenomicPosition& genomic_position) const;

    /**
     * @brief Check whether a whole genomic region is mapped in a chromosome
     * 
     * @param genomic_region is a genomic region
     * @return `true` if and only if `genomic_region` lays in the current chromosome
     */
    bool contains(const GenomicRegion& genomic_region) const;

    /**
     * @brief Check whether any SNV occurs in a possible mutational context
     * 
     * @param genomic_position is the central position of the mutational context
     * @return `true` if and only if no SNVs occurred in the context centered 
     *          in `genomic_position`
     * @throw std::domain_error `genomic_position` does not lays in the fragment
     */
    bool is_mutational_context_free(GenomicPosition genomic_position) const;

    /**
     * @brief Increase the copy number of a genomic region
     * 
     * This method increases the number of copies in a genomic region by selecting 
     * the same allele in each fragment of the region. If not all the fragments touching 
     * the genomic region have the specified allele, the computation ends and the 
     * method returns `false`.
     * 
     * @param genomic_region is the genomic region whose copy number must be increased
     * @param allele_id is the identifier of the allele that must be amplified
     * @return `true` if and only if the all the fragments touching `genomic_region` 
     *      have `allele_id` among their allele ids and the region has been amplified
     * @throw std::domain_error `genomic_region` does not lays in this chromosome
     */
    bool increase_copy_number(const GenomicRegion& genomic_region, const AlleleId& allele_id);

    /**
     * @brief Decrease the copy number of a genomic region removing the last allele
     * 
     * This method decreases the number of copies in a genomic region by selecting 
     * the same allele in each fragment of the region. If not all the fragments touching 
     * the genomic region have the specified allele, the computation ends and the 
     * method returns `false`.
     * 
     * @param genomic_region is the genomic region whose copy number must be decrease
     * @param allele_id is the identifier of the allele that must be deleted
     * @return `true` if and only if the all the fragments touching `genomic_region` 
     *      have `allele_id` among their allele ids and the allele has been removed
     * @throw std::domain_error `genomic_region` does not lays in this chromosome
     */
    bool decrease_copy_number(const GenomicRegion& genomic_region, const AlleleId& allele_id);

    /**
     * @brief Insert a SNV in a context free position
     * 
     * This method tries to insert a SNV. If the chromosome contains another SNV 
     * occurring in the mutational context of the specified SNV, then the insertion 
     * fails.
     * 
     * @param snv is the SNV to be inserted in the chromosome
     * @param allele_id is the identifier of the allele in which the SNV must be placed
     * @return `true` if and only if the chromosome do not 
     *      contain any other SNVs in `snv`'s mutational context
     * @throw std::domain_error `SNV` does not lays in the chromosome
     */
    bool insert(const SNV& snv, const AlleleId allele_id);

    /**
     * @brief Remove a SNV
     * 
     * This method tries to remove a SNV from the chromosome.
     * 
     * @param genomic_position is the position of the SNV to be removed
     * @return `true` if and only if a SNV occurred at `genomic_position`
     * @throw std::domain_error `genomic_position` does not lays in the 
     *      chromosome
     */
    bool remove_SNV(const GenomicPosition& genomic_position);
};

/**
 * @brief A class to represent the passenger mutations of a genome
 */
struct GenomeMutations
{
    using Length = size_t;

    std::map<ChromosomeId, ChromosomeMutations> chromosomes;     //!< the chromosomes

    /**
     * @brief The empty constructor
     */
    GenomeMutations();

    /**
     * @brief A constructor
     * 
     * @param chromosomes is the vector of genome chromosomes
     */
    GenomeMutations(const std::vector<ChromosomeMutations>& chromosomes);    

    /**
     * @brief A constructor
     * 
     * This constructor builds the genome mutations model corresponding 
     * to a set of context positions.
     * 
     * @param context_positions is a context position
     * @param num_of_alleles is the initial number of alleles
     */
    template<typename ABSOLUTE_GENOMIC_POSITIONS>
    GenomeMutations(const ContextPositions<ABSOLUTE_GENOMIC_POSITIONS>& context_positions, const size_t& num_of_alleles)
    {
        std::vector<GenomicRegion> chr_regions = context_positions.get_chromosome_regions();

        for (const auto& chr_region : chr_regions) {
            chromosomes[chr_region.get_chromosome_id()] = ChromosomeMutations(chr_region, num_of_alleles);
        }
    }

    /**
     * @brief Get genome size
     * 
     * @return the genome size
     */
    GenomeMutations::Length size() const;

    /**
     * @brief Get genome allelic size
     * 
     * @return the genome allelic size
     */
    GenomeMutations::Length allelic_size() const;

    /**
     * @brief Get the genome chromosomes
     * 
     * @return the genome chromosomes 
     */
    inline std::map<ChromosomeId, ChromosomeMutations> get_chromosomes() const
    {
        return chromosomes;
    }

    /**
     * @brief Get the identifiers of the alleles that completely include a region
     * 
     * @param genomic_region is a region of the chromosome
     * @return the identifiers of the allele that completely include `genomic_region`
     */
    std::set<AlleleId> get_allele_ids_in(const GenomicRegion& genomic_region) const;

    /**
     * @brief Check whether a chromosomic region has an allele 
     * 
     * @param allele_id is the identifier of the allele whose presence must be checked
     * @param genomic_region is the region on which the check is performed 
     * @return `true` if and only if the whole `genomic_region` contains the allele
     *          having `allele_id` as allele identifier
     * @throw std::domain_error `genomic_region` does not lays in this chromosome
     */
    bool has_allele_on(const AlleleId& allele_id, const GenomicRegion& genomic_region) const;

    /**
     * @brief Increase the copy number of a genomic region
     * 
     * This method increases the number of copies in a genomic region by selecting 
     * the same allele in each fragment of the region. If not all the fragments touching 
     * the genomic region have the specified allele, the computation ends and the 
     * method returns `false`
     * 
     * @param genomic_region is the genomic region whose copy number must be increased
     * @return `true` if and only if the all the fragments touching `genomic_region` 
     *      have `allele_id` among their allele ids and the region has been amplified
     */
    bool increase_copy_number(const GenomicRegion& genomic_region, const AlleleId& allele_id);

    /**
     * @brief Decrease the copy number of a genomic region removing the last allele
     * 
     * This method decreases the number of copies in a genomic region by selecting 
     * the same allele in each fragment of the region. If not all the fragments touching 
     * the genomic region have the specified allele, the computation ends and the 
     * method returns `false`
     * 
     * @param genomic_region is the genomic region whose copy number must be decrease
     * @return the id of the allele to be removed
     * @throw std::domain_error no allele can be removed
     */
    bool decrease_copy_number(const GenomicRegion& genomic_region, const AlleleId& allele_id);

    /**
     * @brief Insert a SNV in a context free position
     * 
     * This method tries to insert a SNV. If the chromosome contains another SNV 
     * occurring in the mutational context of the specified SNV, then the insertion 
     * fails.
     * 
     * @param snv is the SNV to be inserted in the chromosome
     * @return `true` if and only if the chromosome do not 
     *      contain any other SNVs in `snv`'s mutational context
     */
    bool insert(const SNV& snv, const AlleleId& allele_id);

    /**
     * @brief Remove a SNV
     * 
     * This method tries to remove a SNV from the chromosome.
     * 
     * @param genomic_position is the position of the SNV to be removed
     * @return `true` if and only if a SNV occurred at `genomic_position`
     */
    bool remove_SNV(const GenomicPosition& genomic_position);
};

}   // Passengers

}   // Races


#endif // __RACES_GENOME_MUTATIONS__
