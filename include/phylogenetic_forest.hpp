/**
 * @file phylogenetic_forest.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines classes and function for phylogenetic forests
 * @version 0.3
 * @date 2023-12-17
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

#ifndef __RACES_PHYLOGENETIC_FOREST__
#define __RACES_PHYLOGENETIC_FOREST__

#include "descendant_forest.hpp"

#include "genome_mutations.hpp"

namespace Races
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
        std::set<SNV> SNVs;                     //!< The newly introduced SNVs
        std::set<CopyNumberAlteration> CNAs;    //!< The newly introduced CNAs

        /**
         * @brief Save novel mutations in an archive
         *
         * @tparam ARCHIVE is the output archive type
         * @param archive is the output archive
         */
        template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
        inline void save(ARCHIVE& archive) const
        {
            archive & SNVs
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

            archive & mutations.SNVs
                    & mutations.CNAs;

            return mutations;
        }
    };
private:
    std::map<Mutants::CellId, CellGenomeMutations> leaves_mutations;  //!< The mutations of each cells represented as leaves in the forest

    std::map<Mutants::CellId, NovelMutations> novel_mutations;  //!< The mutations introduces by each cell in the forest

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
                return forest->leaves_mutations.at(cell_id);
            }

            throw std::domain_error("The node is not a leaf");
        }

        /**
         * @brief Get the parent node
         *
         * @return the parent node of this node
         */
        const_node parent() const
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
        node parent()
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
         * @param snv is a SVN that was introduced in the corresponding
         *      cell and was not present in the cell parent
         */
        inline void add_new_mutation(const SNV& snv)
        {
            get_forest().novel_mutations[cell_id].SNVs.insert(snv);
        }

        /**
         * @brief Add a newly introduced mutation
         *
         * @param cna is a CNA that was introduced in the corresponding
         *      cell and was not present in the cell parent
         */
        inline void add_new_mutation(const CopyNumberAlteration& cna)
        {
            get_forest().novel_mutations[cell_id].CNAs.insert(cna);
        }

        /**
         * @brief Get the genome mutations of the cell represented by the node
         *
         * @return a reference to the genome mutations of the cell
         *      represented by the node
         */
        inline CellGenomeMutations& cell_mutations()
        {
            if (is_leaf()) {
                return forest->leaves_mutations[cell_id];
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
    inline const std::map<Mutants::CellId, CellGenomeMutations>&
    get_leaves_mutations() const
    {
        return leaves_mutations;
    }

    std::list<SampleGenomeMutations> get_samples_mutations() const;

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
        archive & static_cast<const Mutants::DescendantsForest&>(*this)
                & leaves_mutations
                & novel_mutations;
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
        PhylogeneticForest forest;

        archive & static_cast<Mutants::DescendantsForest&>(forest)
                & forest.leaves_mutations
                & forest.novel_mutations;

        return forest;
    }

    template<typename GENOME_WIDE_POSITION, typename RANDOM_GENERATOR>
    friend class MutationEngine;
};

}   // Mutations

}   // Races

#endif // __RACES_PHYLOGENETIC_TREE__
