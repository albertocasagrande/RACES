/**
 * @file phylogenetic_forest.hpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Define classes and function for phylogenetic trees
 * @version 0.1
 * @date 2023-07-13
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

#include <map>
#include <set>
#include <vector>
#include <queue>

#include "cell.hpp"
#include "binary_logger.hpp"

namespace Races
{

class PhylogeneticForest
{
public:
    using CellDataType = LabelledCell<Time>;
    
private:
    std::map<CellId, CellDataType> cells;          //!< The forest cell id-cell map
 
    std::set<CellId> roots;                        //!< The cell ids of the forest roots

    std::map<CellId, std::set<CellId>> branches;   //!< The descendant branches
public:

    /**
     * @brief A constant node of the forest
     */
    class const_node {
        PhylogeneticForest const* forest;   //!< A pointer to the forest
        CellId cell_id;                     //!< The cell id of a cell in the forest

        /**
         * @brief A constructor for a constant node
         * 
         * @param forest is the forest of the node
         * @param cell_id is the cell id of a cell in the forest
         */
        const_node(const PhylogeneticForest* forest, const CellId cell_id);
    public:
        /**
         * @brief Cast to Cell
         * 
         * This method returns a constant reference to a cell. 
         * Notice that the cells in the forest should be modified 
         * exclusively by using `PhylogeneticForest::node` methods
         * 
         * @return a constant reference to a cell
         */
        operator const PhylogeneticForest::CellDataType&() const;

        /**
         * @brief Get the parent node
         * 
         * @return the parent node of this node
         */
        const_node parent() const;

        /**
         * @brief Get the node direct descendants
         * 
         * @return a vector containing the direct descendants of
         *         this node
         */
        std::vector<const_node> children() const;

        friend class PhylogeneticForest;
    };

    class node {
        PhylogeneticForest* forest;   //!< A pointer to the forest
        CellId cell_id;               //!< The cell id of a cell in the forest

        node(PhylogeneticForest* forest, const CellId cell_id);
    public:
        /**
         * @brief Cast to Cell
         * 
         * This method returns a constant reference to a cell. 
         * Notice that the cells in the forest should be modified 
         * exclusively by using `PhylogeneticForest::node` methods
         * 
         * @return a constant reference to a cell
         */
        operator const PhylogeneticForest::CellDataType&() const;

        /**
         * @brief Get the parent node
         * 
         * @return the parent node of this node
         */
        node parent();

        /**
         * @brief Get the node direct descendants
         * 
         * @return a vector containing the direct descendants of
         *         this node
         */
        std::vector<node> children();

        /**
         * @brief Get the parent node
         * 
         * @return the parent node of this node
         */
        const_node parent() const;

        /**
         * @brief Get the node direct descendants
         * 
         * @return a vector containing the direct descendants of
         *         this node
         */
        std::vector<const_node> children() const;

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
    const_node get_node(const CellId& cell_id) const;

    /**
     * @brief Get a node with the given id
     * 
     * @param cell_id is the id of the aimed cell node
     * @return the corresponding cell node
     */
    node get_node(const CellId& cell_id);

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

    template<typename SAMPLER, typename CELL_STORAGE>
    friend PhylogeneticForest grow_forest_from(SAMPLER& sampler, CELL_STORAGE& cell_storage);
};

/**
 * @brief Grow a forest from a sample of cells
 * 
 * This function grows a forest from a sample of cells.
 * The sample cells are the leaves of the forest and 
 * their ancestors are loaded from a cell storage.
 * 
 * @tparam CELL_SAMPLE is the cell sample type 
 * @tparam CELL_STORAGE is the type of the cell storage
 * @param sample is the cell sample
 * @param cell_storage is the cell storage
 * @return the phylogenetic forest obtained by using the 
 *       cells in `sample` as leafs and recovering 
 */
template<typename CELL_SAMPLE, typename CELL_STORAGE>
PhylogeneticForest grow_forest_from(CELL_SAMPLE& sample, CELL_STORAGE& cell_storage)
{
    PhylogeneticForest forest;

    std::set<CellId> parent_ids;
    for (const auto& cell: sample) {
        parent_ids.insert(cell.get_id());
    }

    std::priority_queue<CellId> queue(parent_ids.begin(), parent_ids.end());

    while (!queue.empty()) {
        auto cell = cell_storage[queue.top()];

        queue.pop();
        parent_ids.erase(cell.get_id());

        // if the cell is not an initial cell
        if (cell.get_id()!=cell.get_parent_id()) {

            // the cell id is not in the queue
            if (parent_ids.count(cell.get_parent_id())==0) {

                // add its id to the queue
                queue.push(cell.get_parent_id());
                parent_ids.insert(cell.get_parent_id());
            }
        } else {
            // it is a root
            forest.roots.insert(cell.get_id());
        }

        forest.cells.insert(std::make_pair(cell.get_id(),cell));
        forest.branches[cell.get_parent_id()].insert(cell.get_id());
    }

    return forest;
}

/* Inline implementations */
inline PhylogeneticForest::const_node::operator const PhylogeneticForest::CellDataType&() const
{
    return forest->cells.at(cell_id);
}

inline PhylogeneticForest::node::operator const PhylogeneticForest::CellDataType&() const
{
    return forest->cells.at(cell_id);
}

inline PhylogeneticForest::const_node PhylogeneticForest::get_node(const CellId& cell_id) const
{
    if (cells.count(cell_id)==0) {
        throw std::runtime_error("The forest does not contain the cell "
                                 "having the specified identifier");
    }
        
    return PhylogeneticForest::const_node(this, cell_id);
}

inline PhylogeneticForest::node PhylogeneticForest::get_node(const CellId& cell_id)
{
    if (cells.count(cell_id)==0) {
        throw std::runtime_error("The forest does not contain the cell "
                                 "having the specified identifier");
    }

    return PhylogeneticForest::node(this, cell_id);
}

}

#endif // __RACES_PHYLOGENETIC_TREE__