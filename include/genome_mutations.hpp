/**
 * @file genome_mutations.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines genome and chromosome data structures
 * @version 0.35
 * @date 2024-05-31
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

#ifndef __RACES_GENOME_MUTATIONS__
#define __RACES_GENOME_MUTATIONS__

#include <map>
#include <list>
#include <fstream>

#include "mutation.hpp"
#include "context_index.hpp"

#include "allele.hpp"
#include "cna.hpp"

#include "cell.hpp"

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

    std::list<CNA> CNAs;        //!< the occurred CNAs

    AlleleId next_allele_id;    //!< the identifier of the next allele

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
    inline const ChromosomeMutations::Length& size() const
    {
        return length;
    }

    /**
     * @brief Get chromosome allelic size
     *
     * @return the chromosome allelic size
     */
    inline const ChromosomeMutations::Length& allelic_size() const
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
     * @brief Get an allele in the chromosome
     *
     * @param allele_id is the identifier of the allele to find
     * @return a constant reference to the allele
     * @throw std::out_of_range `allele_id` is not a valid allele identifier for the
     *          chromosome
     */
    const Allele& get_allele(const AlleleId& allele_id) const;

    /**
     * @brief Get the CNAs occurred in the chromosome
     *
     * @return the list of the CNAs occurred in the chromosome
     */
    inline const std::list<CNA>& get_CNAs() const
    {
        return CNAs;
    }

    /**
     * @brief Get the identifiers of the alleles containing a genomic position
     *
     * @param genomic_position is a position of the chromosome
     * @return a list of the identifiers of the alleles containing `genomic_position`
     */
    std::list<AlleleId> get_alleles_containing(const GenomicPosition& genomic_position) const;

    /**
     * @brief Get the identifiers of the alleles that completely include a region
     *
     * @param genomic_region is a region of the chromosome
     * @return a list of the identifiers of the alleles that completely include `genomic_region`
     */
    std::list<AlleleId> get_alleles_containing(const GenomicRegion& genomic_region) const;

    /**
     * @brief Get the alleles with context free for a mutation
     *
     * @param mutation is a mutation
     * @return a list of the identifiers of the allele in which the context of `mutation`
     *      is free
     */
    std::list<AlleleId> get_alleles_with_context_free_for(const SID& mutation) const;

    /**
     * @brief Get an allele in the chromosome
     *
     * @param allele_id is the identifier of the allele to find
     * @return a non-constant reference to the allele
     * @throw std::out_of_range `allele_id` is not a valid allele identifier for the
     *          chromosome
     */
    Allele& get_allele(const AlleleId& allele_id);

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
    inline bool contains(const GenomicPosition& genomic_position) const
    {
        return (genomic_position.chr_id==id()       // same chromosome
                && genomic_position.position>0      // lays in [1, chromosome size]
                && genomic_position.position<=size());
    }

    /**
     * @brief Check whether a whole mutation lays in a chromosome
     *
     * @param mutation is a SID mutation
     * @return `true` if and only if `mutation` completely lays in the current
     *      chromosome
     */
    bool contains(const SID& mutation) const;

    /**
     * @brief Check whether a whole genomic region is mapped in a chromosome
     *
     * @param genomic_region is a genomic region
     * @return `true` if and only if `genomic_region` lays in the current chromosome
     */
    inline bool contains(const GenomicRegion& genomic_region) const
    {
        return (genomic_region.get_chromosome_id()==id()     // same chromosome
                && genomic_region.get_initial_position()>0   // the lays in [1, chromosome size]
                && genomic_region.get_final_position()<=size());
    }

    /**
     * @brief Check whether a mutation context is free
     *
     * @param mutation is a mutation
     * @return `true` if and only if the allele does not contains
     *      SID mutations in the one-base neighborhood of `mutation`
     * @throw std::domain_error `mutation` does not lays in the chromosome
     */
    bool has_context_free(const SID& mutation) const;

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
     * @param new_allele_id is the identifier of the new allele if amplification
     *      succeeded
     * @param nature is the nature of the amplification
     * @return `true` if and only if the all the fragments touching `genomic_region`
     *      have `allele_id` among their allele ids and the region has been amplified
     * @throw std::domain_error `genomic_region` does not lays in this chromosome
     */
    bool amplify_region(const GenomicRegion& genomic_region, const AlleleId& allele_id,
                        AlleleId& new_allele_id,
                        const Mutation::Nature& nature=Mutation::UNDEFINED);

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
     * @param nature is the nature of the deletion
     * @return `true` if and only if the all the fragments touching `genomic_region`
     *      have `allele_id` among their allele ids and the allele has been removed
     * @throw std::domain_error `genomic_region` does not lays in this chromosome
     */
    bool remove_region(const GenomicRegion& genomic_region, const AlleleId& allele_id,
                       const Mutation::Nature& nature=Mutation::UNDEFINED);

    /**
     * @brief Insert a SID mutation in a context free position
     *
     * This method tries to insert a SID mutation. If the chromosome contains another
     * SID mutation occurring in the mutational context of the specified SID, then
     * the insertion fails.
     *
     * @param mutation is the SID mutation to be inserted in the chromosome
     * @param allele_id is the identifier of the allele in which the SID mutation
     *      must be placed
     * @return `true` if and only if the chromosome do not contain any other
     *      SID mutations in `mutation`'s mutational context
     * @throw std::domain_error `mutation` does not lays in the chromosome
     * @throw std::out_of_range the chromosome has not the allele
     *      `allele_id` or `mutation` does not lay in the allele
     */
    bool insert(const SID& mutation, const AlleleId& allele_id);

    /**
     * @brief Duplicate genomic structure
     *
     * This method duplicates the genomic structure of the current objects.
     * It returns a `ChromosomeMutations` object that has the same alleles of
     * the current objects, but misses the original SID mutations.
     *
     * @return a `ChromosomeMutations` object that has the same alleles of
     *      the current objects, but misses the original SID mutations
     */
    ChromosomeMutations duplicate_structure() const;

    /**
     * @brief Remove a SID mutation
     *
     * This method tries to remove a SID mutation from the chromosome.
     *
     * @param genomic_position is the position of the SID mutation to be
     *      removed
     * @return `true` if and only if a SID occurred at `genomic_position`
     * @throw std::domain_error `genomic_position` does not lays in the
     *      chromosome
     * @throw std::out_of_range the chromosome has not the allele
     *      `allele_id` or `genomic_position` does not lay in the allele
     */
    bool remove_mutation(const GenomicPosition& genomic_position);

    /**
     * @brief Check whether a SID mutation is included among the chromosome mutations
     *
     * @param mutation is the SID mutation whose inclusion among the chromosome
     *      mutations is tested
     * @return `true` if and only if `mutation` is contained among the chromosome
     *      mutations
     */
    bool includes(const SID& mutation) const;

    /**
     * @brief Save chromosome mutations in an archive
     *
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        archive & identifier
                & length
                & allelic_length
                & alleles
                & CNAs
                & next_allele_id;
    }

    /**
     * @brief Load chromosome mutations from an archive
     *
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the load chromosome mutations
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    inline static ChromosomeMutations load(ARCHIVE& archive)
    {
        ChromosomeMutations chr_mutations;

        archive & chr_mutations.identifier
                & chr_mutations.length
                & chr_mutations.allelic_length
                & chr_mutations.alleles
                & chr_mutations.CNAs
                & chr_mutations.next_allele_id;

        return chr_mutations;
    }
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
     * @param chromosome_regions is a list of the chromosome regions
     * @param num_of_alleles is the initial number of alleles
     */
    GenomeMutations(const std::list<GenomicRegion>& chromosome_regions,
                    const size_t& num_of_alleles);

    /**
     * @brief A constructor
     *
     * This constructor builds the genome mutations model corresponding
     * to a set of context positions.
     *
     * @param chromosome_regions is a list of the chromosome regions
     * @param alleles_per_chromosome is the initial number of alleles per chromosome
     */
    GenomeMutations(const std::list<GenomicRegion>& chromosome_regions,
                    const std::map<ChromosomeId, size_t>& alleles_per_chromosome);

    /**
     * @brief A constructor
     *
     * This constructor builds the genome mutations model corresponding
     * to a set of context positions.
     *
     * @param chromosome_regions is a list of the chromosome regions
     * @param num_of_alleles is the initial number of alleles
     */
    inline GenomeMutations(std::list<GenomicRegion>&& chromosome_regions,
                           const size_t& num_of_alleles):
        GenomeMutations(chromosome_regions, num_of_alleles)
    {}

    /**
     * @brief A constructor
     *
     * This constructor builds the genome mutations model corresponding
     * to a set of context positions.
     *
     * @param chromosome_regions is a list of the chromosome regions
     * @param alleles_per_chromosome is the initial number of alleles per chromosome
     */
    inline GenomeMutations(std::list<GenomicRegion>&& chromosome_regions,
                           const std::map<ChromosomeId, size_t>& alleles_per_chromosome):
        GenomeMutations(chromosome_regions, alleles_per_chromosome)
    {}

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
     * @brief Get the absolute chromosome positions
     *
     * @return a map from chromosome identifier to the position of the first base
     *      of the chromosome
     */
    std::map<ChromosomeId, size_t> get_absolute_chromosome_positions() const;

    /**
     * @brief Get the mutations of a genome chromosome
     *
     * @param chromosome_id is the identifier of the aimed chromosome mutations
     * @return the mutations of a genome chromosome
     */
    const ChromosomeMutations& get_chromosome(const ChromosomeId& chromosome_id) const;

    /**
     * @brief Get the mutations of a genome chromosome
     *
     * @param chromosome_id is the identifier of the aimed chromosome mutations
     * @return the mutations of a genome chromosome
     */
    ChromosomeMutations& get_chromosome(const ChromosomeId& chromosome_id);

    /**
     * @brief Get the identifiers of the alleles containing a genomic position
     *
     * @param genomic_position is a position of the genome
     * @return a list of the identifiers of the alleles containing `genomic_position`
     */
    std::list<AlleleId> get_alleles_containing(const GenomicPosition& genomic_position) const;

    /**
     * @brief Get the identifiers of the alleles that completely include a region
     *
     * @param genomic_region is a region of the chromosome
     * @return a list of the identifiers of the allele that completely include `genomic_region`
     */
    std::list<AlleleId> get_alleles_containing(const GenomicRegion& genomic_region) const;

    /**
     * @brief Get the alleles with context free for a mutation
     *
     * @param mutation is a mutation
     * @return a list of the identifiers of the allele in which the context of `mutation`
     *      is free
     */
    std::list<AlleleId> get_alleles_with_context_free_for(const SID& mutation) const;

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
     * method returns `false`.
     *
     * @param genomic_region is the genomic region whose copy number must be increased
     * @param allele_id is the identifier of the allele that must be amplified
     * @param new_allele_id is the identifier of the new allele if amplification
     *      succeeded
     * @param nature is the nature of the amplification
     * @return `true` if and only if the all the fragments touching `genomic_region`
     *      have `allele_id` among their allele ids and the region has been amplified
     * @throw std::domain_error `genomic_region` does not lays in this chromosome
     */
    bool amplify_region(const GenomicRegion& genomic_region, const AlleleId& allele_id,
                        AlleleId& new_allele_id,
                        const Mutation::Nature& nature=Mutation::UNDEFINED);

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
     * @param nature is the nature of the deletion
     * @return `true` if and only if deletion has been successful
     * @throw std::domain_error no allele can be removed
     */
    bool remove_region(const GenomicRegion& genomic_region, const AlleleId& allele_id,
                       const Mutation::Nature& nature=Mutation::UNDEFINED);

    /**
     * @brief Insert a SID mutation in a context free position
     *
     * This method tries to insert a SID mutation. If the chromosome contains another
     * SID mutation occurring in the mutational context of the specified SID, then
     * the insertion fails.
     *
     * @param mutation is the SID mutation to be inserted in the chromosome
     * @param allele_id is the identifier of the allele in which the SID mutation
     *      must be placed
     * @return `true` if and only if the chromosome do not contain any other
     *      SID mutations in `mutation`'s mutational context
     * @throw std::domain_error `mutation` does not lays in the chromosome
     * @throw std::out_of_range the chromosome has not the allele
     *      `allele_id` or `mutation` does not lay in the allele
     */
    bool insert(const SID& mutation, const AlleleId& allele_id);

    /**
     * @brief Remove a SID mutation
     *
     * This method tries to remove a SID mutation from the chromosome.
     *
     * @param genomic_position is the position of the SID to be removed
     * @return `true` if and only if a SID occurred at `genomic_position`
     */
    bool remove_mutation(const GenomicPosition& genomic_position);

    /**
     * @brief Check whether a mutation context is free
     *
     * @param mutation is a mutation
     * @return `true` if and only if the allele does not contains
     *      SID mutations in the one-base neighborhood of `mutation`
     * @throw std::domain_error `mutation` does not lays in the genome
     */
    bool has_context_free(const SID& mutation) const;

    /**
     * @brief Check whether a SID mutation is included among the genome mutations
     *
     * @param mutation is the SID mutation whose inclusion among the genome
     *      mutations is tested
     * @return `true` if and only if `mutation` is contained among the genome
     *      mutations
     */
    bool includes(const SID& mutation) const;

    /**
     * @brief Duplicate genomic structure
     *
     * This method duplicates the genomic structure of the current objects.
     * It returns a `GenomeMutations` object that has the same chromsomes and
     * alleles of the current objects, but misses the original SNVs and indels.
     *
     * @return a `GenomeMutations` object that has the same chromsomes and
     *      alleles of the current objects, but misses the original SNVs and
     *      indels
     */
    GenomeMutations duplicate_structure() const;

    /**
     * @brief Save genome mutations in an archive
     *
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        archive & chromosomes;
    }

    /**
     * @brief Load genome mutations from an archive
     *
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the load genome mutations
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    inline static GenomeMutations load(ARCHIVE& archive)
    {
        GenomeMutations g_mutations;

        archive & g_mutations.chromosomes;

        return g_mutations;
    }
};

/**
 * @brief A class to represent the mutations of a specific cell
 */
struct CellGenomeMutations : public Mutants::Cell, public GenomeMutations
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
     * @brief A constructor for wild-type cells
     *
     * @param germline_mutations are the germline mutations
     */
    explicit CellGenomeMutations(const GenomeMutations& germline_mutations);

    /**
     * @brief A constructor
     *
     * @param cell is the cell
     * @param genome_mutations are the genome mutations
     */
    explicit CellGenomeMutations(const Mutants::Cell& cell, const GenomeMutations& genome_mutations);

    /**
     * @brief Save cell genome mutations in an archive
     *
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        ARCHIVE::write_header(archive, "RACES Genome Mutations", 0);

        archive & static_cast<const Mutants::Cell&>(*this)
                & static_cast<const GenomeMutations&>(*this);
    }

    /**
     * @brief Load cell genome mutations from an archive
     *
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the load cell genome mutations
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    inline static CellGenomeMutations load(ARCHIVE& archive)
    {
        ARCHIVE::read_header(archive, "RACES Genome Mutations", 0);

        CellGenomeMutations cg_mutations;

        archive & static_cast<Mutants::Cell&>(cg_mutations)
                & static_cast<GenomeMutations&>(cg_mutations);

        return cg_mutations;
    }
};

/**
 * @brief A structure representing the genome mutations of a tissue sample
 */
struct SampleGenomeMutations
{
    using CellMutationsPtr = std::shared_ptr<CellGenomeMutations>;

    GenomeMutations germline_mutations;      //!< The germline mutations
    std::list<CellMutationsPtr> mutations;   //!< The list of cell genome mutations

    const std::string name;                  //!< The sample name

    /**
     * @brief The constructor
     *
     * @param name is the sample name
     * @param germline_mutations are the germline mutations
     */
    SampleGenomeMutations(const std::string& name,
                          const GenomeMutations& germline_mutations);
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
