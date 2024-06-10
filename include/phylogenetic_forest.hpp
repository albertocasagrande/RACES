/**
 * @file phylogenetic_forest.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines classes and function for phylogenetic forests
 * @version 1.0
 * @date 2024-06-10
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

#ifndef __RACES_PHYLOGENETIC_FOREST__
#define __RACES_PHYLOGENETIC_FOREST__

#include <memory>

#include "descendant_forest.hpp"

#include "genome_mutations.hpp"

namespace RACES
{

namespace Mutations
{

/**
 * @brief A class representing descendants forests
 */
class PhylogeneticForest : public Mutants::DescendantsForest
{
public:
    /**
     * @brief Mutations introduced in a cell
     *
     * This structure contains the mutations that are present
     * in one of the cells represented in the phylogenetic forest,
     * but not in its parent.
     */
    struct NovelMutations
    {
        std::set<SID> SIDs;    //!< The newly introduced SID mutations
        std::set<CNA> CNAs;    //!< The newly introduced CNAs

        /**
         * @brief Save novel mutations in an archive
         *
         * @tparam ARCHIVE is the output archive type
         * @param archive is the output archive
         */
        template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
        inline void save(ARCHIVE& archive) const
        {
            archive & SIDs
                    & CNAs;
        }

        /**
         * @brief Load novel mutations from an archive
         *
         * @tparam ARCHIVE is the input archive type
         * @param archive is the input archive
         * @return the load novel mutations
         */
        template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
        inline static NovelMutations load(ARCHIVE& archive)
        {
            NovelMutations mutations;

            archive & mutations.SIDs
                    & mutations.CNAs;

            return mutations;
        }
    };
private:
    using CellIdSet = std::set<Mutants::CellId>;

    using CellMutationsPtr = std::shared_ptr<CellGenomeMutations>;

    std::map<Mutants::CellId, CellMutationsPtr> leaves_mutations;   //!< The mutations of each cells represented as leaves in the forest
    std::map<Mutants::CellId, NovelMutations> novel_mutations;      //!< The mutations introduces by each cell in the forest
    std::map<SID, CellIdSet> SID_first_cells;                       //!< A map associating each SID to the first cells in which it occured
    std::map<CNA, CellIdSet> CNA_first_cells;      //!< A map associating each CNA to the first cells in which it occured

    GenomeMutations germline_mutations; //!< The germline mutations

public:
    /**
     * @brief A constant node of the forest
     */
    class const_node : public Mutants::DescendantsForest::_const_node<PhylogeneticForest>
    {
    public:
        /**
         * @brief A constructor for a constant node
         *
         * @param forest is the forest of the node
         * @param cell_id is the cell id of a cell in the forest
         */
        const_node(const PhylogeneticForest* forest, const Mutants::CellId cell_id);

        /**
         * @brief Get the genome mutations of the cell represented by the node
         *
         * @return a constant reference to the genome mutations of the cell
         *      represented by the node
         */
        inline const CellGenomeMutations& cell_mutations() const
        {
            if (is_leaf()) {
                return *(forest->leaves_mutations.at(cell_id));
            }

            throw std::domain_error("The node is not a leaf");
        }

        /**
         * @brief Get the parent node
         *
         * @return the parent node of this node
         */
        inline const_node parent() const
        {
            return Mutants::DescendantsForest::_const_node<PhylogeneticForest>::parent<const_node>();
        }

        /**
         * @brief Get the node direct descendants
         *
         * @return a vector containing the direct descendants of
         *         this node
         */
        inline std::vector<const_node> children() const
        {
            return Mutants::DescendantsForest::_const_node<PhylogeneticForest>::children<const_node>();
        }

        friend class PhylogeneticForest;
    };

    class node : public Mutants::DescendantsForest::_node<PhylogeneticForest>
    {
    public:
        /**
         * @brief A constructor for a constant node
         *
         * @param forest is the forest of the node
         * @param cell_id is the cell id of a cell in the forest
         */
        node(PhylogeneticForest* forest, const Mutants::CellId cell_id);

        /**
         * @brief Get the parent node
         *
         * @return the parent node of this node
         */
        inline node parent()
        {
            return Mutants::DescendantsForest::_node<PhylogeneticForest>::parent<node>();
        }

        /**
         * @brief Get the node direct descendants
         *
         * @return a vector containing the direct descendants of
         *         this node
         */
        inline std::vector<node> children()
        {
            return Mutants::DescendantsForest::_node<PhylogeneticForest>::children<node>();
        }

        /**
         * @brief Add a newly introduced mutation
         *
         * @param mutation is a SID mutation that was introduced in the corresponding
         *      cell and was not present in the cell parent
         */
        void add_new_mutation(const SID& mutation);

        /**
         * @brief Add a newly introduced mutation
         *
         * @param cna is a CNA that was introduced in the corresponding
         *      cell and was not present in the cell parent
         */
        void add_new_mutation(const CNA& cna);

        /**
         * @brief Get the genome mutations of the cell represented by the node
         *
         * @return a reference to the genome mutations of the cell
         *      represented by the node
         */
        inline CellGenomeMutations& cell_mutations()
        {
            if (is_leaf()) {
                return *(forest->leaves_mutations[cell_id]);
            }

            throw std::domain_error("The node is not a leaf");
        }

        friend class PhylogeneticForest;
    };

    /**
     * @brief The empty constructor
     */
    PhylogeneticForest();

    /**
     * @brief Get a constant node with the given id
     *
     * @param cell_id is the id of the aimed cell node
     * @return the corresponding constant cell node
     */
    const_node get_node(const Mutants::CellId& cell_id) const;

    /**
     * @brief Get a node with the given id
     *
     * @param cell_id is the id of the aimed cell node
     * @return the corresponding cell node
     */
    node get_node(const Mutants::CellId& cell_id);

    /**
     * @brief Get the forest roots
     *
     * @return std::vector<const_node>
     */
    std::vector<const_node> get_roots() const;

    /**
     * @brief Get the forest roots
     *
     * @return std::vector<const_node>
     */
    std::vector<node> get_roots();

    /**
     * @brief Get the forest for a subset of the tissue samples
     *
     * @param sample_names are the names of the samples to be considered
     * @return the descendants forest for the tissue samples whose name is
     *          in `sample_names`
     */
    PhylogeneticForest get_subforest_for(const std::vector<std::string>& sample_names) const;

    /**
     * @brief Get the genome mutations of a cell represented as a leaf
     *
     * @param cell_id is a cell identifier of a cell represented as a leaf
     * @return a constant reference to the genome mutations of the cell
     *          with identifier `cell_id` provided that this cell is represented
     *          by a forest leaf
     */
    const CellGenomeMutations& get_leaf_mutations(const Mutants::CellId& cell_id) const;

    /**
     * @brief Get the genome mutations of all the forest leaves
     *
     * @return a constant reference to a map associating all the leaves cell
     *      identifiers to the corresponding genome mutations
     */
    inline const std::map<Mutants::CellId, std::shared_ptr<CellGenomeMutations>>&
    get_leaves_mutations() const
    {
        return leaves_mutations;
    }

    /**
     * @brief Get map associating each SID to the first cell in which it occurs
     *
     * @return a constant reference to a map associating each SID in the phylogenetic
     *         forest to the identifier of the first cell in which the SID occured
     */
    inline const std::map<SID, std::set<Mutants::CellId>>& get_mutation_first_cells() const
    {
        return SID_first_cells;
    }

    /**
     * @brief Get map associating each CNA to the first cell in which it occurs
     *
     * @return a constant reference to a map associating each CNA in the phylogenetic
     *         forest to the identifier of the first cell in which the CNA occured
     */
    inline const std::map<CNA, std::set<Mutants::CellId>>& get_CNA_first_cells() const
    {
        return CNA_first_cells;
    }

    /**
     * @brief Get the list of sample mutations
     *
     * @return the list of sample mutations
     */
    std::list<SampleGenomeMutations> get_sample_mutations_list() const;

    /**
     * @brief Get the list of sample mutations
     *
     * @param sample_name is the name of the sample whose mutations are aimed
     * @return the mutations of the sample whose name is `sample_name`
     */
    SampleGenomeMutations get_sample_mutations(const std::string& sample_name) const;

    /**
     * @brief Get the germline mutations
     *
     * @return A constant reference to germline mutations
     */
    inline const GenomeMutations& get_germline_mutations() const
    {
        return germline_mutations;
    }

    /**
     * @brief Clear the forest
     */
    void clear();

    /**
     * @brief Save a phylogenetic forest in an archive
     *
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        ARCHIVE::write_header(archive, "RACES Phylogenetic Forest", 0);

        archive & static_cast<const Mutants::DescendantsForest&>(*this)
                & novel_mutations
                & SID_first_cells
                & CNA_first_cells
                & germline_mutations;

        archive & leaves_mutations.size();
        for (const auto& [cell_id, ptr]: leaves_mutations) {
            archive & cell_id & *ptr;
        }
    }

    /**
     * @brief Load a phylogenetic forest from an archive
     *
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the load phylogenetic forest
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    inline static PhylogeneticForest load(ARCHIVE& archive)
    {
        ARCHIVE::read_header(archive, "RACES Phylogenetic Forest", 0);

        PhylogeneticForest forest;

        archive & static_cast<Mutants::DescendantsForest&>(forest)
                & forest.novel_mutations
                & forest.SID_first_cells
                & forest.CNA_first_cells
                & forest.germline_mutations;

        size_t num_of_leaves;
        archive & num_of_leaves;
        for (size_t i{0}; i<num_of_leaves; ++i) {
            Mutants::CellId cell_id;
            archive & cell_id;

            auto mut_ptr = std::make_shared<CellGenomeMutations>();
            archive & *mut_ptr;
            forest.leaves_mutations[cell_id] = mut_ptr;
        }

        return forest;
    }

    template<typename GENOME_WIDE_POSITION, typename RANDOM_GENERATOR>
    friend class MutationEngine;
};

}   // Mutations

}   // RACES

#endif // __RACES_PHYLOGENETIC_TREE__
