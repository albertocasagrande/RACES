/**
 * @file mutation_engine.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines a class to place mutations on a descendants forest
 * @version 0.58
 * @date 2024-04-23
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

#ifndef __RACES_MUTATION_ENGINE__
#define __RACES_MUTATION_ENGINE__

#include <map>
#include <stack>
#include <random>
#include <ostream>

#include "phylogenetic_forest.hpp"
#include "mutant_id.hpp"

#include "context_index.hpp"
#include "genome_mutations.hpp"
#include "snv_signature.hpp"
#include "mutational_properties.hpp"
#include "driver_storage.hpp"

#include "filter.hpp"


namespace Races
{

namespace Mutations
{

/**
 * @brief Mutation Statistics
 */
class MutationStatistics
{
    /**
     * @brief SID statistics
     */
    struct SIDStatistics
    {
        size_t mutated_alleles;            //!< Number of mutated alleles
        size_t num_of_cells;               //!< Total number of cells containing the SID

        /**
         * @brief The empty constructor
         */
        SIDStatistics();
    };

    struct SampleMutationStatistics
    {
        size_t num_of_cells;          //!< Number of recorded cells

        std::map<SID, SIDStatistics> SIDs;      //!< SIDs
        std::list<CNA> CNAs;   //!< CNAs

        /**
         * @brief The empty constructor
         */
        SampleMutationStatistics();
    };

    SampleMutationStatistics overall_statistics;    //!< Overall statistics
    std::map<std::string, SampleMutationStatistics> sample_statistics;  //!< Statistics per sample

public:
    /**
     * @brief The empty constructor
     */
    MutationStatistics();

    /**
     * @brief Record the SIDs of a cell genome
     *
     * @param sample_name is the name of the sample from which the cell has
     *          been extracted
     * @param cell_mutations are the mutations of the cell whose statistics
     *          is about to be recorded
     * @return a reference to the updated object
     */
    MutationStatistics& record(const std::string& sample_name,
                               const CellGenomeMutations& cell_mutations);

    /**
     * @brief Record the SIDs in a list of sample genomic mutations
     *
     * @param[in] mutations_list is a list of sample mutations
     * @param[in,out] progress_bar is a progress bar pointer
     * @return a reference to the updated object
     */
    MutationStatistics&
    record(const std::list<Races::Mutations::SampleGenomeMutations>& mutations_list,
           UI::ProgressBar* progress_bar=nullptr);

    /**
     * @brief Record the SIDs in a list of sample genomic mutations
     *
     * @param[in] mutations_list is a list of sample mutations
     * @param[in,out] progress_bar is a progress bar pointer
     * @return a reference to the updated object
     */
    inline MutationStatistics&
    record(const std::list<Races::Mutations::SampleGenomeMutations>& mutations_list,
           UI::ProgressBar& progress_bar)
    {
        return record(mutations_list, &progress_bar);
    }

    /**
     * @brief Write in a stream a table summarizing SID statistics
     *
     * @param os is the output stream
     * @param separator is the column separator
     * @return a reference to the updated stream
     */
    std::ostream& write_SIDs_table(std::ostream& os, const char separator='\t');

    /**
     * @brief Write in a stream a table summarizing CNA statistics
     *
     * @param os is the output stream
     * @param separator is the column separator
     * @return a reference to the updated stream
     */
    std::ostream& write_CNAs_table(std::ostream& os, const char separator='\t');
};

/**
 * @brief The SBS exposure type
 *
 * A SBS is a single base mutation substituion signature that provides for any context
 * (i.e., a triplet of bases) the probability for a SID to occur on that context.
 * The SBSs depends on the environmental context and, because of that, more
 * than one SBSs may be active at the same time with different probabilities.
 * An exposure is a discrete probability distribution among SBS signatures.
 */
using Exposure = std::map<std::string, double>;

/**
 * @brief The mutation engine
 *
 * The objects of this class place mutations on the cell genome according
 * to a descendants forest
 *
 * @tparam GENOME_WIDE_POSITION is the type used to represent genome-wise position
 * @tparam RANDOM_GENERATOR is the type of random generator
 */
template<typename GENOME_WIDE_POSITION = uint32_t, typename RANDOM_GENERATOR = std::mt19937_64>
class MutationEngine
{
    /**
     * @brief The SBS inverse cumulative distributions
     *
     * SBS's are probability distributions for mutational types. If we arbitrary assign
     * an order to mutational type, given a SBS and a mutational type T, we can compute
     * the probability of randomly choosing a mutational type smaller or equal to T.
     * The inverse cumulative SBS type associated such a probability to T itself for
     * all the mutational type T.
     */
    using InverseCumulativeSBS = std::map<double, Races::Mutations::MutationalType>;

    /**
     * @brief The stack of contexts removed from the context index
     */
    using ContextStack = std::stack<std::pair<MutationalContext, GENOME_WIDE_POSITION>>;

    RANDOM_GENERATOR generator; //!< the random generator

    ContextIndex<GENOME_WIDE_POSITION> context_index;  //!< the genome context index
    ContextStack context_stack;                        //!< the stack of the contexts removed from `context_index`

    std::map<Time, Exposure> timed_exposures;          //!< the timed exposures
    std::map<std::string, InverseCumulativeSBS> inv_cumulative_SBSs;        //!< the inverse cumulative SBS

    MutationalProperties mutational_properties;  //!< the species mutational properties

    GenomeMutations germline_mutations; //!< the germline mutations

    DriverStorage driver_storage;  //!< the driver storage

    std::vector<CNA> passenger_CNAs;   //!< the admissible passenger CNAs

    /**
     * @brief Select a random value in a set
     *
     * @tparam T is the type of objects in the set
     * @param values is a list of values from which a random value must be selected
     * @return a random value among those in `values`
     */
    template<typename T>
    T select_random_value(const std::list<T>& values)
    {
        std::uniform_int_distribution<> u_dist(0,values.size()-1);
        size_t pos = u_dist(generator);

        size_t index=0;
        auto it = values.begin();
        while (index!=pos) {
            ++it;
            ++index;
        }

        return *it;
    }

    /**
     * @brief Compute the inverse cumulative SBS

     * SBS's are probability distributions for mutational types. If we arbitrary assign
     * an order to mutational type, given a SBS and a mutational type T, we can compute
     * the probability of randomly choosing a mutational type smaller or equal to T.
     * The inverse cumulative SBS type associated such a probability to T itself for
     * all the mutational type T.
     *
     * @param m_signature is a mutational signature (SBS)
     * @return the corresponding mutational signature
     */
    static InverseCumulativeSBS get_inv_cumulative_SBS(const MutationalSignature& m_signature)
    {
        InverseCumulativeSBS inv_cumulative_SBS;
        double cumulative_prob = 0;
        for (const auto& [m_type, prob]: m_signature) {
            if (prob != 0) {
                cumulative_prob += prob;

                inv_cumulative_SBS.insert({cumulative_prob, m_type});
            }
        }

        return inv_cumulative_SBS;
    }

    /**
     * @brief Get the exposure according to a cell birth time
     *
     * The active exposure depend on the simulation time. It must be selected in
     * agreement with cell birth times.
     *
     * @param cell is the cell whose birth time select the exposure
     * @return a constant reference to exposure that are associated
     *         to `cell` birth time
     */
    const Exposure&
    get_active_exposure(const Mutants::Cell& cell) const
    {
        auto exposure_it = timed_exposures.upper_bound(cell.get_birth_time());

        return (--exposure_it)->second;
    }

    /**
     * @brief Select a random passenger SID
     *
     * This method selects a random passenger SID among whose having a
     * specified mutational type and available in a context index.
     * The selected SID is extracted from the context index and inserted
     * into a stack to revert the selection.
     *
     * @param[in] m_type is the mutational type of the SID to be selected
     * @return a passenger SID whose type is `m_type` and which was
     *          available in `context_index`
     */
    SID select_SID(const MutationalType& m_type)
    {
        using namespace Races::Mutations;

        MutationalContext context = m_type.get_context();
        MutationalContext compl_context = context.get_complemented();

        size_t total_pos = context_index[context].size();
        total_pos += context_index[compl_context].size();

        size_t index;
        {
            std::uniform_int_distribution<> u_dist(0,total_pos-1);
            index = u_dist(generator);
        }

        char alt_base = m_type.get_replace_base();
        if (index >= context_index[context].size()) {
            index -= context_index[context].size();
            context = compl_context;

            alt_base = GenomicSequence::get_complemented(alt_base);
        }

        GENOME_WIDE_POSITION pos;

        if (infinite_sites_model) {
            pos = context_index.extract(context, index);

            context_stack.push({context, pos});
        } else {
            pos = context_index[context][index];
        }
        auto genomic_pos = context_index.get_genomic_position(pos);

        return {genomic_pos.chr_id, genomic_pos.position,
                context.get_central_nucleotide(), alt_base, Mutation::UNDEFINED};
    }

    /**
     * @brief Randomly select the number of mutations according to a Poisson distribution
     *
     * @param genome_size is the genome size
     * @param passenger_mutation_rate is the rate of the passanger mutations
     * @return is the randomly selected number of mutations
     */
    size_t number_of_mutations(GenomeMutations::Length genome_size,
                               const double& passenger_mutation_rate)
    {
        auto mean = static_cast<int>(genome_size*passenger_mutation_rate);

        std::poisson_distribution<> p_dist(mean);

        return p_dist(generator);
    }

    /**
     * @brief Place the CNAs in a cell genome in the phylogenetic forest
     *
     * @param node is a phylogenetic forest node representing a cell
     * @param cell_mutations are the cell mutations
     * @param num_of_mutations is the number of CNA to apply
     * @param nature is the nature of the mutations
     */
    void place_CNAs(PhylogeneticForest::node* node, GenomeMutations& cell_mutations,
                    const size_t& num_of_mutations, const Mutation::Nature& nature)
    {
        // while some of the requested CNAs have not been set
        size_t new_CNAs{0};
        size_t available{passenger_CNAs.size()};
        while (new_CNAs < num_of_mutations) {
            const size_t last_pos = available-1;
            std::uniform_int_distribution<size_t> u_dist(0, last_pos);

            // randomly select a CNA
            size_t index = u_dist(generator);

            bool to_be_placed = true;
            while (to_be_placed) {
                CNA cna = passenger_CNAs[index];

                cna.nature = nature;
                if (place_CNA(cna, cell_mutations)) {
                    // if the CNA has been successfully applied,
                    // then increase the number of new CNAs
                    ++new_CNAs;

                    if (node != nullptr) {
                        node->add_new_mutation(cna);
                    }

                    to_be_placed = false;
                } else {
                    if (--available == 0) {
                        throw std::runtime_error("None of the "
                                                + std::to_string(passenger_CNAs.size())
                                                + " admitted CNAs can be applied.");
                    }

                    if (index<available) {
                        // here passenger_CNAs is shuffled and it will not be
                        // resorted as at the beginning of the computation, but
                        // this is still acceptable because the CNAs are
                        // randomly selected
                        std::swap(passenger_CNAs[index], passenger_CNAs[available]);
                    } else {
                        index = available-1;
                    }
                }
            }
        }
    }

    /**
     * @brief Place the CNAs in a cell genome in the phylogenetic forest
     *
     * @param node is a phylogenetic forest node representing a cell
     * @param cell_mutations are the cell mutations
     * @param CNA_rate is the rate of passegers CNAs
     */
    void place_CNAs(PhylogeneticForest::node& node, GenomeMutations& cell_mutations,
                    const double& CNA_rate)
    {
        const auto num_of_CNAs = number_of_mutations(cell_mutations.allelic_size(), CNA_rate);

        place_CNAs(&node, cell_mutations, num_of_CNAs, Mutation::PASSENGER);
    }

    /**
     * @brief Get the available alleles for placing a SID
     *
     * @param cell_mutations are the non-germinal mutations of a cell
     * @param genomic_position is the position in which we want to place the SID
     * @param infinite_sites_model is a Boolean flag to enable/disable the infinite
     *      sites model
     * @return a list of the identifiers of the available alleles
     */
    std::list<AlleleId> get_available_alleles_for_SID(const GenomeMutations& cell_mutations,
                                                      const GenomicPosition& genomic_position,
                                                      const bool& infinite_sites_model) const
    {
        std::list<AlleleId> available_alleles;

        if (infinite_sites_model) {
            const auto& chr_mutations = cell_mutations.get_chromosome(genomic_position.chr_id);
            if (germline_mutations.has_context_free(genomic_position) &&
                    chr_mutations.has_context_free(genomic_position)) {

                for (const auto& [allele_id, allele] : chr_mutations.get_alleles()) {
                    available_alleles.push_back(allele_id);
                }
            }

            return available_alleles;
        }

        auto germline_alleles = germline_mutations.get_alleles_with_context_free_for(genomic_position);

        if (germline_alleles.size()==0) {
            return available_alleles;
        }

        std::set<AlleleId> germline_set(germline_alleles.begin(), germline_alleles.end());

        const auto& chr_mutations = cell_mutations.get_chromosome(genomic_position.chr_id);

        for (const auto& [allele_id, allele] : chr_mutations.get_alleles()) {
            if (allele.has_context_free(genomic_position)
                    && germline_set.count(allele.get_history().front())>0) {
                available_alleles.push_back(allele_id);
            }
        }

        return available_alleles;
    }

    /**
     * @brief Place a passenger SID
     *
     * @param snv is the SID to place
     * @param cell_mutations are the cell mutations
     * @param infinite_sites_model is a Boolean flag to enable/disable the infinite
     *      sites model
     * @return `true` if and only if the SID placement has succeed. This
     *          method returns `false` when the SID cannot be placed
     *          because its context is not free or no allele are available
     *          for it.
     */
    bool place_SID(const SID& snv, GenomeMutations& cell_mutations,
                   const bool& infinite_sites_model)
    {
        auto candidate_alleles = get_available_alleles_for_SID(cell_mutations, snv,
                                                               infinite_sites_model);

        if (candidate_alleles.size()==0) {
            return false;
        }

        AlleleId allele_id = select_random_value(candidate_alleles);

        // try to apply the selected SID
        return cell_mutations.insert(snv, allele_id);
    }

    /**
     * @brief Place a passenger SID
     *
     * @param snv is the SID to place
     * @param cell_mutations are the cell mutations
     * @return `true` if and only if the SID placement has succeed. This
     *          method returns `false` when the SID cannot be placed
     *          because its context is not free or no allele are available
     *          for it.
     */
    inline bool place_SID(const SID& snv, GenomeMutations& cell_mutations)
    {
        return place_SID(snv, cell_mutations, infinite_sites_model);
    }

    /**
     * @brief Place a driver SID with an allele specification
     *
     * @param snv is the SID to place
     * @param cell_mutations are the cell mutations
     * @return `true` if and only if the SID placement has succeed. This
     *          method returns `false` when the SID cannot be placed
     *          because its context is not free or no allele are available
     *          for it.
     */
    bool place_SID(const MutationSpec<SID>& snv,
                   GenomeMutations& cell_mutations)
    {
        if (snv.allele_id == RANDOM_ALLELE) {
            return place_SID(static_cast<const SID&>(snv), cell_mutations,
                             false);
        }

        return cell_mutations.insert(snv, snv.allele_id);
    }

    /**
     * @brief Place the SIDs in a cell genome in the phylogenetic forest
     *
     * @param node is a phylogenetic forest node representing a cell
     * @param cell_mutations are the cell mutations
     * @param SBS_name is the name of the SBS used to find the mutations
     * @param num_of_mutations is the number of SID to apply
     * @param nature is the nature of the mutations
     */
    void place_SIDs(PhylogeneticForest::node* node, GenomeMutations& cell_mutations,
                    const std::string& SBS_name, const size_t& num_of_mutations,
                    const Mutation::Nature& nature)
    {
        // get the inverse cumulative SBS
        const auto& inv_cumulative_SBS = inv_cumulative_SBSs.at(SBS_name);

        std::uniform_real_distribution<> u_dist(0,1);

        // while some of the requested SIDs have not been set
        size_t new_SIDs{0}, context_misses{0};
        while (new_SIDs < num_of_mutations) {
            // randomly pick a mutational type according to the current SBS
            auto it = inv_cumulative_SBS.lower_bound(u_dist(generator));

            if (it != inv_cumulative_SBS.end() &&
                    context_index[it->second.get_context()].size()>0) {

                context_misses = 0;

                // select a SID among those having the picked mutational type
                // and that available in the context index
                SID snv = select_SID(it->second);

                snv.cause = SBS_name;
                snv.nature = nature;

                if (place_SID(snv, cell_mutations)) {
                    // if the SID has been successfully applied,
                    // then increase the number of new SIDs
                    ++new_SIDs;

                    if (node != nullptr) {
                        node->add_new_mutation(snv);
                    }
                }
            } else {
                const size_t max_context_misses{1000};

                if (++context_misses >= max_context_misses) {
                    throw std::runtime_error("Missed available context for \""
                                             + SBS_name + "\" for "
                                             + std::to_string(max_context_misses)
                                             + " in row. Try to reduce the "
                                             + "context sampling.");
                }
            }
        }
    }

    /**
     * @brief Place the SIDs in a cell genome in the phylogenetic forest
     *
     * @tparam GENOME_WIDE_POSITION is the type used to represent genome-wise position
     * @param node is a phylogenetic forest node representing a cell
     * @param cell_mutations are the cell mutations
     * @param SID_rate is the rate of passegers SIDs
     */
    void place_SIDs(PhylogeneticForest::node& node, GenomeMutations& cell_mutations,
                    const double& SID_rate)
    {
        const auto num_of_SIDs = number_of_mutations(cell_mutations.allelic_size(), SID_rate);

        // for each active SBS probability in the active exposure
        for (const auto& [SBS_name, probability]: get_active_exposure(node)) {

            // evaluate how many of the SIDs due to the current SBS
            const size_t SBS_num_of_SIDs = probability*num_of_SIDs;

            place_SIDs(&node, cell_mutations, SBS_name, SBS_num_of_SIDs,
                       Mutation::PASSENGER);
        }
    }

    /**
     * @brief Seach for a random allele in a chromosome containing a region
     *
     * @param[out] allele_id is the identifier of the found allele
     * @param[in] chr_mutations is the chromosome in which the allele is searched
     * @param[in] region is the region that must be contained by the found allele
     * @param[in] no_driver_mutations is a Boolean flag to establish whether
     *          alleles containing driver mutations in the region are discharged
     * @return `true` if and only if one of the alleles of `chr_mutations`
     *          contains the region `region`. In this case, the parameter
     *          `allele_id` is set to the identifier of such an allele.
     *          When `no_driver_mutations` is set to `true`, the search
     *          must avoid alleles containing driver mutations in `region`.
     */
    bool select_allele_containing(AlleleId& allele_id,
                                  const ChromosomeMutations& chr_mutations,
                                  const GenomicRegion& region,
                                  const bool& no_driver_mutations)
    {
        const auto& chr_alleles = chr_mutations.get_alleles();

        if (chr_alleles.size()==0) {
            return false;
        }

        size_t first_idx_to_test = (generator() % chr_alleles.size());
        size_t idx{0};
        for (const auto& [a_id, allele]: chr_alleles) {
            if (idx >= first_idx_to_test) {
                if (allele.contains(region)
                    && !(no_driver_mutations
                         && allele.has_driver_mutations_in(region))) {
                    allele_id = a_id;
                    return true;
                }
            }
            ++idx;
        }
        idx = 0;
        for (const auto& [a_id, allele]: chr_alleles) {
            if (idx >= first_idx_to_test) {
                return false;
            }
            if (allele.contains(region)
                    && !(no_driver_mutations
                         && allele.has_driver_mutations_in(region))) {
                allele_id = a_id;
                return true;
            }
            ++idx;
        }

        return false;
    }

    /**
     * @brief Seach for an allele in which the CNA can be applied
     *
     * @param[in,out] cna is the CNA which should be applied
     * @param[in] chr_mutations is the chromosome in which the allele is searched
     * @param[in] no_driver_mutations is a Boolean flag to establish whether
     *          alleles containing driver mutations in the region are discharged
     * @return `true` if and only if one of the alleles of `chr_mutations`
     *          admits the application of `cna`.
     *          When `no_driver_mutations` is set to `true`, the search
     *          must avoid alleles containing driver mutations in `region`.
     *          If `cna` is an amplification and the source is `RANDOM_ALLELE` or
     *          the `cna` is a deletion and the destination is `RANDOM_ALLELE`,
     *          then, when `true` is returned the source or the destination of
     *          `cna`, respectively, are updated to the identifier of the found
     *          allele.
     */
    bool find_allele_for(CNA& cna, ChromosomeMutations& chr_mutations,
                         const bool& no_driver_mutations)
    {
        const auto cna_region = cna.get_region();
        switch(cna.type) {
            case CNA::Type::AMPLIFICATION:
                if (cna.source == RANDOM_ALLELE) {
                    return select_allele_containing(cna.source, chr_mutations,
                                                    cna_region, no_driver_mutations);
                }
                return true;
            case CNA::Type::DELETION:
                if (cna.dest == RANDOM_ALLELE) {
                    return select_allele_containing(cna.dest, chr_mutations,
                                                    cna_region, no_driver_mutations);
                }
                return true;
            default:
                throw std::runtime_error("Unsupported CNA type");
        }
    }

    /**
     * @brief Try to place a CNA
     *
     * @param[in,out] cna is the CNA to place
     * @param[in] cell_mutations are the cell mutations
     * @return `true` if and only if the CNA placement has succeed. This method
     *          returns `false` when the CNA cannot be placed because no allele
     *          are available for it.
     *          When `no_driver_mutations` is set to `true`, the search
     *          must avoid alleles containing driver mutations in `region`.
     *          If `cna` is an amplification and the source is `RANDOM_ALLELE` or
     *          the `cna` is a deletion and the destination is `RANDOM_ALLELE`,
     *          then, when `true` is returned the source or the destination of
     *          `cna`, respectively, are updated to the identifier of the found
     *          allele.
     */
    bool place_CNA(CNA& cna, GenomeMutations& cell_mutations)
    {
        auto& chr_mutations = cell_mutations.get_chromosome(cna.chr_id);
        if (!find_allele_for(cna, chr_mutations, (cna.nature != Mutation::DRIVER))) {
            return false;
        }

        const auto cna_region = cna.get_region();
        switch(cna.type) {
            case CNA::Type::AMPLIFICATION:
                return chr_mutations.amplify_region(cna_region, cna.source,
                                                    cna.dest, cna.nature);
            case CNA::Type::DELETION:
                return chr_mutations.remove_region(cna_region, cna.dest,
                                                   cna.nature);
            default:
                throw std::runtime_error("Unsupported CNA type");
        }
    }

    /**
     * @brief Place the driver mutations
     *
     * @param node is a phylogenetic forest node representing a cell
     * @param cell_mutations are the cell mutations
     * @param cell_statistics are the statistics of the cell mutations
     * @param driver_mutations is the map associating mutant ids and mutations
     */
    void place_driver_mutations(PhylogeneticForest::node& node,
                                GenomeMutations& cell_mutations,
                                const std::map<Mutants::MutantId, DriverMutations>& driver_mutations)
    {
        if (node.is_root() || node.get_mutant_id()!=node.parent().get_mutant_id()) {
            const auto& mutant_mp = driver_mutations.at(node.get_mutant_id());

            for (auto snv : mutant_mp.SIDs) {
                snv.cause = mutant_mp.name;
                snv.nature = Mutation::DRIVER;

                if (place_SID(snv, cell_mutations)) {
                    node.add_new_mutation(snv);
                } else {
                    std::ostringstream oss;

                    oss << mutant_mp.name << "'s driver mutation "
                        << snv << " cannot be placed.";

                    throw std::runtime_error(oss.str());
                }
            }

            for (auto cna : mutant_mp.CNAs) {
                cna.nature = Mutation::DRIVER;
                if (place_CNA(cna, cell_mutations)) {
                    node.add_new_mutation(cna);
                } else {
                    std::ostringstream oss;

                    oss << mutant_mp.name << "'s driver mutation "
                        << cna << " cannot be placed.";

                    throw std::runtime_error(oss.str());
                }
            }
        }
    }

    /**
     * @brief Place the mutations on the genomes of a descendants forest node
     *
     * This method recursively places the mutations on the nodes of a descendants
     * forest. All the genome mutations associated to the forest leaves are
     * collected and saved in a container that partition them according
     * to the cell sample.
     *
     * @param[in] node is a phylogenetic forest node representing a cell
     * @param[in] ancestor_mutations is the genomic mutation of the ancestor
     * @param[in] passenger_rates is the map associating species ids to their passenger rates
     * @param[in] driver_mutations is the map associating mutant ids and mutations
     * @param[in] ancestor_mutations is the genomic mutation of the ancestor
     * @param[in,out] visited_nodes is the number of visited nodes
     * @param[in,out] progress_bar is a progress bar pointer
     */
    void place_mutations(PhylogeneticForest::node& node, const GenomeMutations& ancestor_mutations,
                         const std::map<Mutants::SpeciesId, PassengerRates>& passenger_rates,
                         const std::map<Mutants::MutantId, DriverMutations>& driver_mutations,
                         size_t& visited_nodes, UI::ProgressBar *progress_bar)
    {
        using namespace Races::Mutations;

        const size_t context_stack_size = context_stack.size();

        CellGenomeMutations cell_mutations(static_cast<const Mutants::Cell&>(node),
                                           ancestor_mutations);

        place_driver_mutations(node, cell_mutations, driver_mutations);

        auto rates_it = passenger_rates.find(node.get_species_id());
        if (rates_it == passenger_rates.end()) {
            const auto& species_name = node.get_forest().get_species_name(node.get_species_id());

            throw std::runtime_error("Unknown species \""+ species_name +"\"");
        }
        place_SIDs(node, cell_mutations, rates_it->second.snv);
        place_CNAs(node, cell_mutations, rates_it->second.cna);

        ++visited_nodes;
        if (progress_bar != nullptr) {

            size_t percentage = (100*visited_nodes)/(node.get_forest()).num_of_nodes();
            if (percentage>progress_bar->get_progress()) {
                progress_bar->set_progress(percentage);
            }
        }

        if (node.is_leaf()) {
            auto mut_ptr = std::make_shared<CellGenomeMutations>();

            std::swap(cell_mutations, *mut_ptr);

            auto& phylo_forest = node.get_forest();
            auto cell_id = static_cast<const Mutants::Cell&>(node).get_id();

            phylo_forest.leaves_mutations[cell_id] = mut_ptr;
        } else {
            for (auto child: node.children()) {
                place_mutations(child, cell_mutations, passenger_rates, driver_mutations,
                                visited_nodes, progress_bar);
            }
        }

        if (!infinite_sites_model) {
            // reverse context index extractions
            while (context_stack.size() > context_stack_size) {
                reinsert_last_context();
            }
        }
    }

    /**
     * @brief Remove the context stack top and insert it in the context index
     */
    inline void reinsert_last_context()
    {
        const auto& top_stack = context_stack.top();
        context_index.insert(top_stack.first, top_stack.second);
        context_stack.pop();
    }

    /**
     * @brief Get the species rate map
     *
     * This method returns a map from species identifier to the corresponding
     * passenger rates.
     *
     * @param descendants_forest is a descendent forest
     * @return a map associating a species in `descendants_forest` to its
     *          passenger rates
     */
    std::map<Mutants::SpeciesId, PassengerRates>
    get_species_rate_map(const Mutants::DescendantsForest& descendants_forest) const
    {
        using namespace Races::Mutants;

        std::map<SpeciesId, PassengerRates> species_rates;

        for (const auto& [species_id, species_data] : descendants_forest.get_species_data()) {
            const auto mutant_name = descendants_forest.get_mutant_name(species_data.mutant_id);
            const auto species_name = mutant_name
                                      + MutantProperties::signature_to_string(species_data.signature);

            auto passenger_it = mutational_properties.get_passenger_rates().find(species_name);

            if (passenger_it == mutational_properties.get_passenger_rates().end()) {
                throw std::runtime_error("\""+species_name+"\" has no mutational rate");
            }

            species_rates[species_id] = passenger_it->second;
        }

        return species_rates;
    }

    /**
     * @brief Get the driver mutation map
     *
     * This method returns a map from the mutant identifier to the correspoding
     * driver mutations.
     *
     * @param descendants_forest is a descendent forest
     * @return a map associating the mutants in `descendants_forest` to the
     *              correspoding driver mutations.
     */
    std::map<Mutants::MutantId, DriverMutations>
    get_driver_mutation_map(const Mutants::DescendantsForest& descendants_forest) const
    {
        using namespace Races::Mutants;

        std::map<MutantId, DriverMutations> driver_mutations;

        for (const auto& [species_id, species_data] : descendants_forest.get_species_data()) {
            const auto mutant_name = descendants_forest.get_mutant_name(species_data.mutant_id);

            auto driver_it = mutational_properties.get_driver_mutations().find(mutant_name);

            if (driver_it == mutational_properties.get_driver_mutations().end()) {
                throw std::runtime_error("\"" + mutant_name + "\" is unknown");
            }
            driver_mutations[species_data.mutant_id] = driver_it->second;
        }

        return driver_mutations;
    }

    /**
     * @brief Filter a list of CNAs according the chromosomes in the context index
     *
     * This method returns the vector of CNAs that are contained in the parameter
     * vector and lay in one of the chromosomes mentioned in the context index.
     *
     * @param CNAs is a list of CNAs
     * @return the vector of CNAs that are contained in `CNAs` and lay in one of
     *          the chromosomes mentioned in the context index.
     */
    std::vector<CNA>
    filter_CNA_by_chromosome_ids(const std::vector<CNA>& CNAs) const
    {
        auto chr_ids = context_index.get_chromosome_ids();

        std::set<ChromosomeId> chr_id_set(chr_ids.begin(),chr_ids.end());

        std::vector<CNA> filtered;
        for (const auto& cna : CNAs) {
            if (chr_id_set.count(cna.chr_id)>0) {
                filtered.push_back(cna);
            }
        }

        return filtered;
    }

    /**
     * @brief Check germline-reference sequence consistency
     *
     * This method checks whether germline contains the reference sequence chromosomes
     * and throws a `std::out_of_range` exception this is not the case.
     *
     * @param context_index is a context index
     * @param germline_mutations are the germline mutations
     */
    static void check_genomes_consistency(const ContextIndex<GENOME_WIDE_POSITION>& context_index,
                                          const GenomeMutations& germline_mutations)
    {
        const auto context_chr_ids = context_index.get_chromosome_ids();
        const auto& germline_chromosomes = germline_mutations.get_chromosomes();

        for (const auto chr_id : context_chr_ids) {
            if (germline_chromosomes.find(chr_id) == germline_chromosomes.end()) {
                throw std::out_of_range("The germline does not contain Chr. "
                                        + GenomicPosition::chrtos(chr_id)
                                        + " which is present in the reference.");
            }
        }
    }
public:
    bool infinite_sites_model;   //!< a flag to enable/disable infinite sites model

    /**
     * @brief The empty constructor
     */
    MutationEngine():
        infinite_sites_model(true)
    {}

    /**
     * @brief A constructor
     *
     * @param context_index is the genome context index
     * @param mutational_signatures is the map of the mutational signatures
     * @param germline_mutations are the germline mutations
     * @param driver_storage is the storage of driver mutations
     * @param passenger_CNAs is the vector of the admissible passenger CNAs
     */
    MutationEngine(ContextIndex<GENOME_WIDE_POSITION>& context_index,
                   const std::map<std::string, MutationalSignature>& mutational_signatures,
                   const GenomeMutations& germline_mutations,
                   const DriverStorage& driver_storage,
                   const std::vector<CNA>& passenger_CNAs={}):
        MutationEngine(context_index, mutational_signatures, MutationalProperties(),
                       germline_mutations, driver_storage, passenger_CNAs)
    {}

    /**
     * @brief A constructor
     *
     * @param context_index is the genome context index
     * @param mutational_signatures is the map of the mutational signatures
     * @param mutational_properties are the mutational properties of all the species
     * @param germline_mutations are the germline mutations
     * @param driver_storage is the storage of driver mutations
     * @param passenger_CNAs is the vector of the admissible passenger CNAs
     */
    MutationEngine(ContextIndex<GENOME_WIDE_POSITION>& context_index,
                   const std::map<std::string, MutationalSignature>& mutational_signatures,
                   const MutationalProperties& mutational_properties,
                   const GenomeMutations& germline_mutations,
                   const DriverStorage& driver_storage,
                   const std::vector<CNA>& passenger_CNAs={}):
        generator(), context_index(context_index), mutational_properties(mutational_properties),
        germline_mutations(germline_mutations), driver_storage(driver_storage),
        infinite_sites_model(true)
    {
        MutationEngine::check_genomes_consistency(context_index, germline_mutations);

        for (const auto& [name, mutational_signature]: mutational_signatures) {
            inv_cumulative_SBSs[name] = get_inv_cumulative_SBS(mutational_signature);
        }

        this->passenger_CNAs = filter_CNA_by_chromosome_ids(passenger_CNAs);
    }

    /**
     * @brief Add the properties of a mutant
     *
     * This method add the properties of a mutant (i.e., its passenger rates,
     * its driver mutations) and all its species to the mutations engine.
     *
     * @param name is the name of the mutant
     * @param epistate_passenger_rates is a map from epigenomic status to
     *          passenger rate
     * @param driver_SIDs is a list of driver SIDs
     * @param driver_CNAs is a list of driver CNAs
     * @return a reference to the updated object
     */
    inline
    MutationEngine& add_mutant(const std::string& name,
                               const std::map<std::string, PassengerRates>& epistate_passenger_rates,
                               const std::list<MutationSpec<SID>>& driver_SIDs={},
                               const std::list<CNA>& driver_CNAs={})
    {
        mutational_properties.add_mutant(name, epistate_passenger_rates, driver_SIDs,
                                         driver_CNAs);

        return *this;
    }

    /**
     * @brief Add a timed exposure
     *
     * This method add a exposure that will be applied from the specified simulation
     * time on.
     *
     * @param time is the simulation time from which the exposure will be applied
     * @param exposure is a SBS exposure
     * @return a reference to the updated object
     */
    MutationEngine& add(const Time& time, Exposure&& exposure)
    {
        if (time<0) {
            throw std::domain_error("Simulation time is a non-negative value");
        }

        if (timed_exposures.count(time)>0) {
            throw std::runtime_error("Another exposure has "
                                     "been set for the specified time");
        }

        for (const auto& [name, coeff]: exposure) {
            if (inv_cumulative_SBSs.count(name)==0) {
                throw std::runtime_error("Unknown signature "+name);
            }
        }

        timed_exposures.emplace(Time(time), std::move(exposure));

        return *this;
    }

    /**
     * @brief Add a timed exposure
     *
     * This method add a exposure that will be applied from the specified simulation
     * time on.
     *
     * @param time is the simulation time from which the exposure will be applied
     * @param exposure is an SBS exposure
     * @return a reference to the updated object
     */
    MutationEngine& add(const Time& time, const Exposure& exposure)
    {
        return add(time, Exposure(exposure));
    }

    /**
     * @brief Add the default exposure
     *
     * This method add a default exposure.
     *
     * @param default_exposure is the exposure at time 0
     * @return a reference to the updated object
     */
    MutationEngine& add(const Exposure& default_exposure)
    {
        return add(0, default_exposure);
    }

    /**
     * @brief Place genomic mutations on a descendants forest
     *
     * @param descendants_forest is a descendants forest
     * @param num_of_preneoplatic_mutations is the number of preneoplastic mutations
     * @param seed is the random generator seed
     * @return a phylogenetic forest having the structure of `descendants_forest`
     */
    inline
    PhylogeneticForest place_mutations(const Mutants::DescendantsForest& descendants_forest,
                                       const size_t& num_of_preneoplatic_mutations,
                                       const int& seed=0)
    {
        return place_mutations(descendants_forest, num_of_preneoplatic_mutations, nullptr, seed);
    }

    /**
     * @brief Place genomic mutations on a descendants forest
     *
     * @param descendants_forest is a descendants forest
     * @param num_of_preneoplatic_mutations is the number of preneoplastic mutations
     * @param progress_bar is a progress bar pointer
     * @param seed is the random generator seed
     * @return a phylogenetic forest having the structure of `descendants_forest`
     */
    PhylogeneticForest place_mutations(const Mutants::DescendantsForest& descendants_forest,
                                       const size_t& num_of_preneoplatic_mutations,
                                       UI::ProgressBar *progress_bar, const int& seed=0)
    {
        using namespace Races::Mutants;
        using namespace Races::Mutants::Evolutions;
        using namespace Races::Mutations;

        if (!has_default_exposure()) {
            throw std::runtime_error("The mutation engine default exposure "
                                     "has not been set yet.");
        }

        generator.seed(seed);

        PhylogeneticForest forest;

        forest.germline_mutations = germline_mutations;

        static_cast<DescendantsForest&>(forest) = descendants_forest;

        auto species_rates = get_species_rate_map(descendants_forest);
        auto driver_mutations = get_driver_mutation_map(descendants_forest);

        auto chr_regions = context_index.get_chromosome_regions();

        auto wild_type_structure = germline_mutations.duplicate_structure();

        size_t visited_node = 0;
        for (auto& root: forest.get_roots()) {
            GenomeMutations mutations = wild_type_structure;

            // place preneoplastic mutations
            place_SIDs(&root, mutations, "SBS1", num_of_preneoplatic_mutations,
                       Mutation::PRENEOPLASTIC);

            place_mutations(root, mutations, species_rates, driver_mutations,
                            visited_node, progress_bar);
        }

        // reverse context index extractions
        while (context_stack.size() > 0) {
            reinsert_last_context();
        }

        return forest;
    }

    /**
     * @brief Place genomic mutations on a descendants forest
     *
     * @param descendants_forest is a descendants forest
     * @param num_of_preneoplatic_mutations is the number of preneoplastic mutations
     * @param progress_bar is a progress bar pointer
     * @param seed is the random generator seed
     * @return a phylogenetic forest having the structure of `descendants_forest`
     */
    inline PhylogeneticForest place_mutations(const Mutants::DescendantsForest& descendants_forest,
                                              const size_t& num_of_preneoplatic_mutations,
                                              UI::ProgressBar &progress_bar, const int& seed=0)
    {
        return place_mutations(descendants_forest, num_of_preneoplatic_mutations, &progress_bar, seed);
    }

    /**
     * @brief Get the a sample of wild-type cells
     *
     * @param num_of_cells is the number of cells in the wild-type sample
     * @param num_of_somatic_mutations is the number of somatic mutations
     * @return a sample of `num_of_cells` wild-type cells
     */
    SampleGenomeMutations get_wild_type_sample(const size_t& num_of_cells,
                                               const size_t& num_of_somatic_mutations)
    {
        using namespace Mutants::Evolutions;

        SampleGenomeMutations sample_mutations("wild-type", germline_mutations);

        auto germline_structure = germline_mutations.duplicate_structure();
        for (size_t i=0; i<num_of_cells; ++i) {
            auto mutations = std::make_shared<CellGenomeMutations>(germline_structure);

            // place preneoplastic mutations
            place_SIDs(nullptr, *mutations, "SBS1", num_of_somatic_mutations,
                       Mutation::PRENEOPLASTIC);

            sample_mutations.mutations.push_back(mutations);
        }

        return sample_mutations;
    }

    /**
     * @brief Check whether a default exposure has been set
     *
     * @return `true` if and only if a default exposure has been set
     */
    inline bool has_default_exposure() const
    {
        return timed_exposures.find(0.0) != timed_exposures.end();
    }

    /**
     * @brief Get the mutational properties
     *
     * @return a constant reference to the mutational properties
     */
    inline const MutationalProperties& get_mutational_properties() const
    {
        return mutational_properties;
    }

    /**
     * @brief Get the exposures
     *
     * @return a constant reference to the exposures
     */
    inline const std::map<Time, Exposure>&
    get_timed_exposures() const
    {
        return timed_exposures;
    }

    /**
     * @brief Get the driver storage
     *
     * @return a constant reference to the driver storage
     */
    inline const DriverStorage& get_driver_storage() const
    {
        return driver_storage;
    }

    /**
     * @brief Bring back the context index to its original state
     *
     * This method may be used after an exception to bring back the
     * context index to its original state.
     */
    void reset_context_index()
    {
        while (!context_stack.empty()) {
            reinsert_last_context();
        }
    }
};

}   // Mutations

}   // Races

#endif // __RACES_MUTATION_ENGINE__
