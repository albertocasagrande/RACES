/**
 * @file mutation_engine.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines a class to place mutations on a descendants forest
 * @version 0.32
 * @date 2023-12-21
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
     * @brief SNV statistics
     */
    struct SNVStatistics
    {
        size_t mutated_alleles;            //!< Number of mutated alleles
        size_t num_of_cells;               //!< Total number of cells containing the SNV

        /**
         * @brief The empty constructor
         */
        SNVStatistics();
    };

    struct SampleMutationStatistics
    {
        size_t num_of_cells;          //!< Number of recorded cells

        std::map<SNV, SNVStatistics> SNVs;      //!< SNVs
        std::list<CopyNumberAlteration> CNAs;   //!< CNAs

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
     * @brief Record the SNVs of a cell genome
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
     * @brief Record the SNVs in a list of sample genomic mutations
     *
     * @param[in] mutations_list is a list of sample mutations
     * @param[in,out] progress_bar is a progress bar pointer
     * @return a reference to the updated object
     */
    MutationStatistics&
    record(const std::list<Races::Mutations::SampleGenomeMutations>& mutations_list,
           UI::ProgressBar* progress_bar=nullptr);

    /**
     * @brief Record the SNVs in a list of sample genomic mutations
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
     * @brief Write in a stream a table summarizing SNV statistics
     *
     * @param os is the output stream
     * @param separator is the column separator
     * @return a reference to the updated stream
     */
    std::ostream& write_SNVs_table(std::ostream& os, const char separator='\t');

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
 * (i.e., a triplet of bases) the probability for a SNV to occur on that context. 
 * The SBSs depends on the environmental context and, because of that, more 
 * than one SBSs may be active at the same time with different probabilities.
 * An exposure is a discrite probability distribution among SBS signatures.
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
     * @brief The stack of contexts removed from a context index
     */
    using ContextStack = std::stack<std::pair<MutationalContext, GENOME_WIDE_POSITION>>;

    RANDOM_GENERATOR generator; //!< the random generator

    ContextIndex<GENOME_WIDE_POSITION> context_index;  //!< the genome context index
    ContextStack context_stack;                        //!< the stack of the contexts removed from `context_index`

    std::map<ChromosomeId, size_t> alleles_per_chromosome;   //!< number of initial alleles per chromosome

    std::map<Time, Exposure> timed_exposures;          //!< the timed exposures
    std::map<std::string, InverseCumulativeSBS> inv_cumulative_SBSs;        //!< the inverse cumulative SBS

    MutationalProperties mutational_properties;  //!< the species mutational properties

    /**
     * @brief Select a random value in a set
     *
     * @tparam T is the type of objects in the set
     * @param values is the set of values from which a random value must be selected
     * @return a random value among those in `values`
     */
    template<typename T>
    T select_random_value(const std::set<T>& values)
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
     * @brief Select a random SNV
     *
     * This method selects a random SNV among whose having a specified
     * mutational type and available in a context index. The selected
     * SNV is extracted from the context index and inserted into a stack
     * to revert the selection.
     *
     * @param[in] m_type is the mutational type of the SNV to be selected
     * @param[in] cause is the SNV cause
     * @return a SNV whose type is `m_type` and which was available in
     *          `context_index`
     */
    SNV select_SNV(const MutationalType& m_type, const std::string& cause)
    {
        using namespace Races::Mutations;

        MutationalContext context = m_type.get_context();
        MutationalContext complement_context = context.get_complement();

        size_t total_pos = context_index[context].size();
            total_pos += context_index[complement_context].size();

        size_t index;
        {
            std::uniform_int_distribution<> u_dist(0,total_pos-1);
            index = u_dist(generator);
        }

        char mutated_base = m_type.get_replace_base();
        if (index >= context_index[context].size()) {
            index -= context_index[context].size();
            context = complement_context;

            mutated_base = MutationalContext::get_complement(mutated_base);
        }

        GENOME_WIDE_POSITION pos = context_index.extract(context, index);

        context_stack.push({context, pos});

        auto genomic_pos = context_index.get_genomic_position(pos);

        return {genomic_pos, context, mutated_base, cause};
    }

    /**
     * @brief Randomly select the number of SNVs according to a Poisson distribution
     *
     * @param genome_size is the genome size
     * @param mean_passenger_mutations is the mean mutations per duplication
     * @return is the randomly selected number of SNVs
     */
    size_t number_of_SNVs(GenomeMutations::Length genome_size,
                          const double& mean_passenger_mutations)
    {
        auto SNVs_mean = static_cast<int>(genome_size*mean_passenger_mutations);

        std::poisson_distribution<> p_dist(SNVs_mean);

        return p_dist(generator);
    }

    /**
     * @brief Try to place a SNV
     *
     * @param snv is the SNV to place
     * @param cell_mutations are the cell mutations
     * @return `true` if and only if the SNV placement has succeed. This
     *          method returns `false` when the SNV cannot be placed
     *          because its context is not free or no allele are available
     *          for it.
     */
    bool place_SNV(const SNV& snv, GenomeMutations& cell_mutations)
    {
        auto allele_having_pos = cell_mutations.get_alleles_containing(snv);
        if (allele_having_pos.size()>0) {
            AlleleId allele_id = select_random_value(allele_having_pos);

            // try to apply the selected SNV
            if (cell_mutations.insert(snv, allele_id)) {
                return true;
            }
        }

        return false;
    }

    /**
     * @brief Place the mutations associated to a cell in the phylogenetic forest
     *
     * @tparam GENOME_WIDE_POSITION is the type used to represent genome-wise position
     * @param node is a phylogenetic forest node representing a cell
     * @param cell_mutations are the cell mutations
     * @param species_rates is the map associating species ids and mutational rates
     */
    void place_SNVs(PhylogeneticForest::node& node, GenomeMutations& cell_mutations,
                    const std::map<Mutants::SpeciesId, double>& species_rates)
    {
        auto num_of_SNVs = number_of_SNVs(cell_mutations.allelic_size(),
                                          species_rates.at(node.get_species_id()));

        // get the active SBS exposure
        const auto& exposure = get_active_exposure(node);

        // for each active SBS probability
        for (const auto& [SBS_name, probability]: exposure) {

            // get the inverse cumulative SBS
            const auto& inv_cumulative_SBS = inv_cumulative_SBSs.at(SBS_name);

            // evaluate how many of the SNVs occurred to the considered cells
            // are due to the current SBS
            size_t SBS_num_of_SNVs = probability*num_of_SNVs;

            std::uniform_real_distribution<> u_dist(0,1);

            // while some of the requested SNVs have not been set
            unsigned int new_SNVs = 0;
            while (new_SNVs < SBS_num_of_SNVs) {
                // randomly pick a mutational type according to the current SBS
                auto it = inv_cumulative_SBS.lower_bound(u_dist(generator));

                if (it != inv_cumulative_SBS.end() &&
                        context_index[it->second.get_context()].size()>0) {

                    // select a SNV among those having the picked mutational type
                    // and that available in the context index
                    SNV snv = select_SNV(it->second, SBS_name);

                    if (place_SNV(snv, cell_mutations)) {
                        // if the SNV has been successfully applied,
                        // then increase the number of new SNVs
                        ++new_SNVs;

                        node.add_new_mutation(snv);
                    }
                }
            }
        }
    }

    /**
     * @brief Try to place a CNA
     *
     * @param CNA is the CNA to place
     * @param cell_mutations are the cell mutations
     * @return `true` if and only if the CNA placement has succeed. This
     *          method returns `false` when the CNA cannot be placed
     *          because no allele are available for it.
     */
    bool place_CNA(const CopyNumberAlteration& CNA, GenomeMutations& cell_mutations)
    {
        switch(CNA.type) {
            case CopyNumberAlteration::Type::AMPLIFICATION:
                return cell_mutations.amplify_region(CNA.region, CNA.source);
            case CopyNumberAlteration::Type::DELETION:
                return cell_mutations.remove_region(CNA.region, CNA.source);
            default:
                throw std::runtime_error("Unsupported CNA type");
        }
    }

    /**
     * @brief Place the driver mutations
     *
     * @param node is a phylogenetic forest node representing a cell
     * @param mutations are the cell mutations
     * @param cell_statistics are the statistics of the cell mutations
     * @param driver_mutations is the map associating mutant ids and mutations
     */
    void place_driver_mutations(PhylogeneticForest::node& node,
                                GenomeMutations& cell_mutations,
                                const std::map<Mutants::MutantId, DriverMutations>& driver_mutations)
    {
        if (node.is_root() || node.get_mutant_id()!=node.parent().get_mutant_id()) {
            const auto& mutant_mp = driver_mutations.at(node.get_mutant_id());

            for (auto snv : mutant_mp.SNVs) {
                snv.cause = mutant_mp.name;

                if (place_SNV(snv, cell_mutations)) {
                    node.add_new_mutation(snv);
                }
            }

            for (auto CNA : mutant_mp.CNAs) {
                if (place_CNA(CNA, cell_mutations)) {
                    node.add_new_mutation(CNA);
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
     * @param[in] species_rates is the map associating species ids and mutational rates
     * @param[in] driver_mutations is the map associating mutant ids and mutations
     * @param[in] ancestor_mutations is the genomic mutation of the ancestor
     * @param[in,out] visited_nodes is the number of visited nodes
     * @param[in,out] progress_bar is a progress bar pointer
     */
    void place_mutations(PhylogeneticForest::node& node, const GenomeMutations& ancestor_mutations,
                         const std::map<Mutants::SpeciesId, double>& species_rates,
                         const std::map<Mutants::MutantId, DriverMutations>& driver_mutations,
                         size_t& visited_nodes, UI::ProgressBar *progress_bar)
    {
        using namespace Races::Mutations;

        size_t context_stack_size = context_stack.size();

        CellGenomeMutations cell_mutations(static_cast<const Mutants::Cell&>(node),
                                           ancestor_mutations);

        place_driver_mutations(node, cell_mutations, driver_mutations);
        place_SNVs(node, cell_mutations, species_rates);
        //place_CNAs(node, mutations);

        ++visited_nodes;
        if (progress_bar != nullptr) {

            size_t percentage = (100*visited_nodes)/(node.get_forest()).num_of_nodes();
            if (percentage>progress_bar->get_progress()) {
                progress_bar->set_progress(percentage);
            }
        }

        if (node.is_leaf()) {
            auto& phylo_forest = node.get_forest();
            auto cell_id = static_cast<const Mutants::Cell&>(node).get_id();
            phylo_forest.leaves_mutations[cell_id] = cell_mutations;
        } else {
            for (auto child: node.children()) {
                place_mutations(child, cell_mutations, species_rates, driver_mutations,
                                visited_nodes, progress_bar);
            }
        }

        // reverse context index extractions
        while (context_stack.size() > context_stack_size) {
            reinsert_last_context();
        }
    }

    /**
     * @brief Remove the context stack top and insert it in the context index
     */
    void reinsert_last_context()
    {
        const auto& top_stack = context_stack.top();
        context_index.insert(top_stack.first, top_stack.second);
        context_stack.pop();
    }

public:
    /**
     * @brief The empty constructor
     */
    MutationEngine()
    {}

    /**
     * @brief A constructor
     *
     * @param context_index is the genome context index
     * @param alleles_per_chromosome is the number of alleles in wild-type cells
     * @param mutational_signatures is the map of the mutational signatures
     */
    MutationEngine(ContextIndex<GENOME_WIDE_POSITION>& context_index,
                   const std::map<ChromosomeId, size_t>& alleles_per_chromosome,
                   const std::map<std::string, MutationalSignature>& mutational_signatures):
        MutationEngine(context_index, alleles_per_chromosome, mutational_signatures,
                       MutationalProperties())
    {}

    /**
     * @brief A constructor
     *
     * @param context_index is the genome context index
     * @param alleles_per_chromosome is the number of alleles in wild-type cells
     * @param mutational_signatures is the map of the mutational signatures
     * @param mutational_properties are the mutational properties of all the species
     */
    MutationEngine(ContextIndex<GENOME_WIDE_POSITION>& context_index,
                   const std::map<ChromosomeId, size_t>& alleles_per_chromosome,
                   const std::map<std::string, MutationalSignature>& mutational_signatures,
                   const MutationalProperties& mutational_properties):
        generator(), context_index(context_index), alleles_per_chromosome(alleles_per_chromosome),
        mutational_properties(mutational_properties)
    {
        for (const auto& [name, mutational_signature]: mutational_signatures) {
            inv_cumulative_SBSs[name] = get_inv_cumulative_SBS(mutational_signature);
        }
    }

    /**
     * @brief Add the properties of a mutant
     *
     * This method add the properties of a mutant (i.e., its mutation rates, its SNVs and
     * CNAs) and all its species to the mutations engine.
     *
     * @param name is the name of the mutant
     * @param epistate_mutation_rates is a map from epigenomic status to
     *          mutational rate
     * @param mutant_SNVs is a list of SNVs characterizing the mutant
     * @param mutant_CNAs is a list of CNAs characterizing the mutant
     * @return a reference to the updated object
     */
    inline
    MutationEngine& add_mutant(const std::string& name,
                               const std::map<std::string, double>& epistate_mutation_rates,
                               const std::list<SNV>& mutant_SNVs={},
                               const std::list<CopyNumberAlteration>& mutant_CNAs={})
    {
        mutational_properties.add_mutant(name, epistate_mutation_rates, mutant_SNVs,
                                         mutant_CNAs);

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
     * @param seed is the random generator seed
     * @return a phylogenetic forest having the structure of `descendants_forest`
     */
    inline
    PhylogeneticForest place_mutations(const Mutants::DescendantsForest& descendants_forest,
                                       const int& seed=0)
    {
        return place_mutations(descendants_forest, nullptr, seed);
    }

    /**
     * @brief Place genomic mutations on a descendants forest
     *
     * @param descendants_forest is a descendants forest
     * @param progress_bar is a progress bar pointer
     * @param seed is the random generator seed
     * @return a phylogenetic forest having the structure of `descendants_forest`
     */
    PhylogeneticForest place_mutations(const Mutants::DescendantsForest& descendants_forest,
                                       UI::ProgressBar *progress_bar, const int& seed=0)
    {
        using namespace Races::Mutants;
        using namespace Races::Mutants::Evolutions;
        using namespace Races::Mutations;

        generator.seed(seed);

        PhylogeneticForest forest;

        static_cast<DescendantsForest&>(forest) = descendants_forest;

        GenomeMutations mutations(context_index, alleles_per_chromosome);

        std::map<SpeciesId, double> species_rates;
        std::map<MutantId, DriverMutations> driver_mutations;

        for (const auto& [species_id, species_data] : descendants_forest.get_species_data()) {
            const auto mutant_name = descendants_forest.get_mutant_name(species_data.mutant_id);
            const auto species_name = mutant_name
                                      + MutantProperties::signature_to_string(species_data.signature);

            auto species_it = mutational_properties.get_species_rates().find(species_name);

            if (species_it == mutational_properties.get_species_rates().end()) {
                throw std::runtime_error("\""+species_name+"\" has no mutational rate");
            }

            species_rates[species_id] = species_it->second;

            driver_mutations[species_data.mutant_id] = mutational_properties.get_driver_mutations().at(mutant_name);
        }

        size_t visited_node = 0;
        for (auto& root: forest.get_roots()) {
            place_mutations(root, mutations, species_rates, driver_mutations,
                            visited_node, progress_bar);
        }

        return forest;
    }

    /**
     * @brief Place genomic mutations on a descendants forest
     *
     * @param descendants_forest is a descendants forest
     * @param progress_bar is a progress bar pointer
     * @param seed is the random generator seed
     * @return a phylogenetic forest having the structure of `descendants_forest`
     */
    inline PhylogeneticForest place_mutations(const Mutants::DescendantsForest& descendants_forest,
                                              UI::ProgressBar &progress_bar, const int& seed=0)
    {
        return place_mutations(descendants_forest, &progress_bar, seed);
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
