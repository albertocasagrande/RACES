/**
 * @file genome_mutations.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines genome and chromosome data structures
 * @version 0.14
 * @date 2023-12-09
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
#include "context_index.hpp"

#include "allele.hpp"
#include "cna.hpp"

#include "cell.hpp"
#include "tissue_sample.hpp"

#include "progress_bar.hpp"

namespace Races 
{

/**
 * @brief A namespace for classes related to mutations
 */
namespace Mutations
{

/**
 * @brief A class to represent the mutations in a chromosome
 */
class ChromosomeMutations
{
    /**
     * @brief The chromosome length type
     */
    using Length = GenomicRegion::Length;
private:
    ChromosomeId identifier;    //!< the chromosome identifier
    Length length;              //!< the chromosome length
    Length allelic_length;      //!< the sum of the lengths of all the alleles

    std::map<AlleleId, Allele> alleles;  //!< the chromosome alleles

    std::list<CopyNumberAlteration> CNAs;   //!< the occurred CNAs

    AlleleId next_allele_id;   //!< the identifier of the next allele

    /**
     * @brief Find an allele in the chromosome
     * 
     * @param allele_id is the identifier of the allele to find
     * @return a constant reference to the allele
     * @throw std::out_of_range `allele_id` is not a valid allele identifier for the 
     *          chromosome
     */
    const Allele& find_allele(const AlleleId& allele_id) const;

    /**
     * @brief Find an allele in the chromosome
     * 
     * @param allele_id is the identifier of the allele to find
     * @return a non-constant reference to the allele
     * @throw std::out_of_range `allele_id` is not a valid allele identifier for the 
     *          chromosome
     */
    Allele& find_allele(const AlleleId& allele_id);

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
     * @brief Get the chromosome alleles
     * 
     * @return a constant reference to chromosome alleles
     */
    inline const std::map<AlleleId, Allele>& get_alleles() const
    {
        return alleles;
    }

    /**
     * @brief Get the CNAs occurred in the chromosome
     * 
     * @return the list of the CNAs occurred in the chromosome
     */
    inline const std::list<CopyNumberAlteration>& get_CNAs() const
    {
        return CNAs;
    }

    /**
     * @brief Get the identifiers of the alleles containing a genomic position 
     * 
     * @param genomic_position is a position of the chromosome
     * @return the identifiers of the alleles containing `genomic_position`
     */
    std::set<AlleleId> get_alleles_containing(const GenomicPosition& genomic_position) const;

    /**
     * @brief Get the identifiers of the alleles that completely include a region
     * 
     * @param genomic_region is a region of the chromosome
     * @return the identifiers of the alleles that completely include `genomic_region`
     */
    std::set<AlleleId> get_alleles_containing(const GenomicRegion& genomic_region) const;

    /**
     * @brief Check whether an allele contains a chromosomic region
     * 
     * @param allele_id is the identifier of an allele
     * @param genomic_region is the region on which the check is performed 
     * @return `true` if and only if the chromosome have an allele having `allele_id` 
     *          as allele identifier and the whole `genomic_region` is contained
     *          in this allele
     * @throw std::domain_error `genomic_region` does not lays in this chromosome
     * @throw std::out_of_range the chromosome has not the allele `allele_id`
     */
    bool allele_contains(const AlleleId& allele_id, const GenomicRegion& genomic_region) const;

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
    bool has_context_free(const GenomicPosition& genomic_position) const;

    /**
     * @brief Amplify a genomic region
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
    bool amplify_region(const GenomicRegion& genomic_region, const AlleleId& allele_id);

    /**
     * @brief Remove a genomic region
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
    bool remove_region(const GenomicRegion& genomic_region, const AlleleId& allele_id);

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
     * @throw std::out_of_range the chromosome has not the allele 
     *      `allele_id` or `SNV` does not lay in the allele
     */
    bool insert(const SNV& snv, const AlleleId& allele_id);

    /**
     * @brief Remove a SNV
     * 
     * This method tries to remove a SNV from the chromosome.
     * 
     * @param genomic_position is the position of the SNV to be removed
     * @return `true` if and only if a SNV occurred at `genomic_position`
     * @throw std::domain_error `genomic_position` does not lays in the 
     *      chromosome
     * @throw std::out_of_range the chromosome has not the allele 
     *      `allele_id` or `genomic_position` does not lay in the allele
     */
    bool remove_SNV(const GenomicPosition& genomic_position);
};

/**
 * @brief A class to represent the mutations of a genome
 */
class GenomeMutations
{
    std::map<ChromosomeId, ChromosomeMutations> chromosomes;     //!< the chromosomes

public:
    /**
     * @brief The genome length type
     */
    using Length = size_t;

    /**
     * @brief The empty constructor
     */
    GenomeMutations();

    /**
     * @brief A constructor
     * 
     * @param chromosomes is the vector of genome chromosomes
     */
    explicit GenomeMutations(const std::vector<ChromosomeMutations>& chromosomes);    

    /**
     * @brief A constructor
     * 
     * This constructor builds the genome mutations model corresponding 
     * to a set of context positions.
     * 
     * @param context_index is a context index
     * @param num_of_alleles is the initial number of alleles
     */
    template<typename GENOME_WIDE_POSITIONS>
    GenomeMutations(const ContextIndex<GENOME_WIDE_POSITIONS>& context_index, const size_t& num_of_alleles)
    {
        std::vector<GenomicRegion> chr_regions = context_index.get_chromosome_regions();

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
     * @brief Get the mutations of genome chromosomes
     * 
     * @return the mutations of genome chromosomes
     */
    inline const std::map<ChromosomeId, ChromosomeMutations>& get_chromosomes() const
    {
        return chromosomes;
    }

    /**
     * @brief Get the mutations of a genome chromosome
     * 
     * @param chromosome_id is the identifier of the aimed chromosome mutations
     * @return the mutations of a genome chromosome
     */
    inline const ChromosomeMutations& get_chromosome(const ChromosomeId& chromosome_id) const
    {
        return chromosomes.at(chromosome_id);
    }

    /**
     * @brief Get the identifiers of the alleles containing a genomic position 
     * 
     * @param genomic_position is a position of the genome
     * @return the identifiers of the alleles containing `genomic_position`
     */
    std::set<AlleleId> get_alleles_containing(const GenomicPosition& genomic_position) const;

    /**
     * @brief Get the identifiers of the alleles that completely include a region
     * 
     * @param genomic_region is a region of the chromosome
     * @return the identifiers of the allele that completely include `genomic_region`
     */
    std::set<AlleleId> get_alleles_containing(const GenomicRegion& genomic_region) const;

    /**
     * @brief Check whether an allele contains a chromosomic region
     * 
     * @param allele_id is the identifier of an allele
     * @param genomic_region is the region on which the check is performed 
     * @return `true` if and only if the whole `genomic_region` is contained
     *          in the allele `allele_id` of the corresponding chromosome
     * @throw std::domain_error `genomic_region` does not lays in this chromosome
     * @throw std::out_of_range the chromosome has not the allele `allele_id`
     */
    bool allele_contains(const AlleleId& allele_id, const GenomicRegion& genomic_region) const;

    /**
     * @brief Amplify a genomic region
     * 
     * This method increases the number of copies in a genomic region by selecting 
     * the same allele in each fragment of the region. If not all the fragments touching 
     * the genomic region have the specified allele, the computation ends and the 
     * method returns `false`
     * 
     * @param genomic_region is the genomic region whose copy number must be increased
     * @param allele_id is the identifier of the allele to be amplified
     * @return `true` if and only if the all the fragments touching `genomic_region` 
     *      have `allele_id` among their allele ids and the region has been amplified
     */
    bool amplify_region(const GenomicRegion& genomic_region, const AlleleId& allele_id);

    /**
     * @brief Remove a genomic region
     * 
     * This method decreases the number of copies in a genomic region by selecting 
     * the same allele in each fragment of the region. If not all the fragments touching 
     * the genomic region have the specified allele, the computation ends and the 
     * method returns `false`
     * 
     * @param genomic_region is the genomic region whose copy number must be decrease
     * @param allele_id is the identifier of the allele from which the region must be removed
     * @return `true` if and only if deletion has been successful 
     * @throw std::domain_error no allele can be removed
     */
    bool remove_region(const GenomicRegion& genomic_region, const AlleleId& allele_id);

    /**
     * @brief Insert a SNV in a context free position
     * 
     * This method tries to insert a SNV. If the chromosome contains another SNV 
     * occurring in the mutational context of the specified SNV, then the insertion 
     * fails.
     * 
     * @param snv is the SNV to be inserted in the chromosome
     * @param allele_id is the identifier of the allele in which the SNV must be applied
     * @return `true` if and only if the chromosome do not contain any other SNVs in 
     *      `snv`'s mutational context
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


/**
 * @brief A class to represent the mutations of a specific cell
 */
struct CellGenomeMutations : public Clones::Cell, public GenomeMutations
{
    /**
     * @brief The genome length type
     */
    using Length = size_t;

    /**
     * @brief The empty constructor
     */
    CellGenomeMutations();

    /**
     * @brief A constructor
     * 
     * @param chromosomes is the vector of genome chromosomes
     */
    explicit CellGenomeMutations(const Clones::Cell& cell, const GenomeMutations& genome_mutations);    
};

/**
 * @brief A structure representing the genome mutations of a tissue sample
 */
struct SampleGenomeMutations : public Clones::Evolutions::TissueSample
{
    std::list<CellGenomeMutations> mutations;   //!< The list of cell genome mutations

    /**
     * @brief A constructor
     * 
     * @param sample is the tissue sample whose cell mutations is represented
     */
    SampleGenomeMutations(const Clones::Evolutions::TissueSample& sample);
};

}   // Mutations

}   // Races


namespace std
{

/**
 * @brief Write chromosome mutations data in a stream
 * 
 * @param os is the output stream
 * @param chromosome_mutations is a chromosome mutations
 * @return a reference to output stream
 */
std::ostream& operator<<(std::ostream& os, const Races::Mutations::ChromosomeMutations& chromosome_mutations);

/**
 * @brief Write genome mutations data in a stream
 * 
 * @param os is the output stream
 * @param genome_mutations is a genome mutations
 * @return a reference to output stream
 */
std::ostream& operator<<(std::ostream& os, const Races::Mutations::GenomeMutations& genome_mutations);

}   // std

#endif // __RACES_GENOME_MUTATIONS__
