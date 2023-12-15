/**
 * @file phylogenetic_forest.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines classes and function for phylogenetic forests
 * @version 0.1
 * @date 2023-12-15
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
    std::map<Mutants::CellId, CellGenomeMutations> cell_mutations;  //!< The mutations of each cell represented in the forest

public:
    /**
     * @brief A constant node of the forest
     */
    class const_node : public Mutants::DescendantsForest::_const_node<PhylogeneticForest>
    {
    protected:
        /**
         * @brief A constructor for a constant node
         *
         * @param forest is the forest of the node
         * @param cell_id is the cell id of a cell in the forest
         */
        const_node(const PhylogeneticForest* forest, const Mutants::CellId cell_id);

    public:
        /**
         * @brief Get the genome mutations of the cell represented by the node
         *
         * @return a constant reference to the genome mutations of the cell
         *      represented by the node
         */
        inline const CellGenomeMutations& get_mutations() const
        {
            return forest->cell_mutations.at(cell_id);
        }

        friend class PhylogeneticForest;
    };

    class node : public Mutants::DescendantsForest::_node<PhylogeneticForest>, public const_node
    {
    protected:
        /**
         * @brief A constructor for a constant node
         *
         * @param forest is the forest of the node
         * @param cell_id is the cell id of a cell in the forest
         */
        node(PhylogeneticForest* forest, const Mutants::CellId cell_id);

        friend class PhylogeneticForest;
    };

    /**
     * @brief The empty constructor
     */
    PhylogeneticForest();

    /**
     * @brief Construct a phylogentic forest
     *
     * This method builds a phylogentic forest by labelling each cell represented
     * in a descendant forest by using the set of its genome mutations.
     *
     * @param descendant_forest is a descendant forest
     * @param cell_mutations is a map labelling each cell in the forest by its genome mutations
     */
    PhylogeneticForest(const Mutants::DescendantsForest& descendant_forest,
                       std::map<Mutants::CellId, CellGenomeMutations>&& cell_mutations);

    /**
     * @brief Construct a phylogentic forest
     *
     * This method builds a phylogentic forest by labelling each cell represented
     * in a descendant forest by using the set of its genome mutations.
     *
     * @param descendant_forest is a descendant forest
     * @param cell_mutations is a map labelling each cell in the forest by its genome mutations
     */
    PhylogeneticForest(const Mutants::DescendantsForest& descendant_forest,
                       const std::map<Mutants::CellId, CellGenomeMutations>& cell_mutations);

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
     * @brief Clear the forest
     */
    void clear();
};

}   // Mutations

}   // Races

#endif // __RACES_PHYLOGENETIC_TREE__
