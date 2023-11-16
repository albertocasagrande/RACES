/**
 * @file phylogenetic_forest.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines classes and function for phylogenetic forests
 * @version 0.20
 * @date 2023-11-16
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
#include <string>

#include "cell.hpp"
#include "binary_logger.hpp"
#include "tissue_sample.hpp"

namespace Races
{

namespace Drivers 
{

/**
 * @brief A class representing descendants forests
 */
class DescendantsForest
{
    /**
     * @brief Species data
     * 
     * This structure privatly stores species data.
     */
    struct SpeciesData
    {
        GenotypeId genotype_id;             //!< The species genotype id
        MethylationSignature signature;     //!< The species signature

        /**
         * @brief The constructor
         */
        SpeciesData(const GenotypeId& genotype_id, const MethylationSignature& signature);
    }; 

    using EdgeTail = std::set<CellId>;
 
    std::set<CellId> roots;                 //!< The cell ids of the forest roots
    std::map<CellId, Cell> cells;          //!< The forest cell id-cell map
    std::map<CellId, EdgeTail> branches;    //!< The descendant branches

    std::map<SpeciesId, SpeciesData> species_data;          //!< The species id to data map
    std::map<GenotypeId, std::string> genotype_names;       //!< The genotype id to genotype name map

    std::vector<Simulation::TissueSample> samples;  //!< The vector of the samples that produced the forest
    std::map<CellId, uint16_t> coming_from;  //!< The map associating each leaf to the sample which it comes from

    /**
     * @brief Grow a forest from a sample of cells
     * 
     * This function grows a forest from a sample of cells.
     * The sample cells are the leaves of the forest and 
     * their ancestors are loaded from a cell storage.
     * 
     * @tparam CELL_STORAGE is the type of the cell storage
     * @param sample_ids is the cell id sample
     * @param cell_storage is the cell storage
     */
    template<typename CELL_STORAGE>
    void grow_from(const std::list<CellId>& sample, CELL_STORAGE& cell_storage)
    {
        std::set<CellId> parent_ids;
        for (const auto& cell_id: sample) {
            parent_ids.insert(cell_id);

            // record leaves children, i.e., none
            branches[cell_id] = std::set<CellId>();
        }

        std::priority_queue<CellId> queue(parent_ids.begin(), parent_ids.end());

        while (!queue.empty()) {
            auto cell = cell_storage.at(queue.top());

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
                branches[cell.get_parent_id()].insert(cell.get_id());
            } else {
                // it is a root
                roots.insert(cell.get_id());
            }

            cells.insert(std::make_pair(cell.get_id(),cell));
        }
    }

    /**
     * @brief Grow a forest from a list of tissue samples
     * 
     * This method resets a forest and grows it from a list of tissue samples and 
     * a cell storage. The cells in the samples will be the leaves of the forest 
     * and their ancestors are loaded from a cell storage.
     * 
     * @tparam CELL_STORAGE is the type of the cell storage
     * @param tissue_samples is a list of tissue samples coming from the simulation
     * @param cell_storage is the cell storage
     */
    template<typename CELL_STORAGE>
    void grow_from(const std::list<Simulation::TissueSample>& tissue_samples, 
                   CELL_STORAGE& cell_storage)
    {
        clear();
 
        samples = std::vector<Simulation::TissueSample>(tissue_samples.begin(), 
                                                        tissue_samples.end());
 
        uint16_t i{0};
        for (const auto& sample: samples) {
            for (const auto& cell_id: sample.get_cell_ids()) {
                coming_from.insert(std::make_pair(cell_id, i));
            }
            ++i;
        }
 
        std::list<CellId> cell_ids;
        for (auto sample_it = tissue_samples.begin(); sample_it != tissue_samples.end(); ++sample_it) {
            cell_ids.insert(cell_ids.end(), sample_it->get_cell_ids().begin(), 
                            sample_it->get_cell_ids().end());
        }
         
        grow_from(cell_ids, cell_storage);
    }

protected:
    /**
     * @brief Get the cell ids to cells maps
     * 
     * @return a constant reference to the cell ids to cells maps
     */
    inline const std::map<CellId, Cell>& get_cells() const
    {
        return cells;
    }

    /**
     * @brief Get the name of a genotype
     * 
     * @param genotype_id is the identifier of the genotype whose name is aimed
     * @return a constant reference to the name of the genotype having 
     *          `genotype_id` as identifier
     */
    inline const std::string& get_genotype_name(const GenotypeId& genotype_id) const
    {
        return genotype_names.at(genotype_id);
    }

    /**
     * @brief Get the map associating each leaf to the sample which it comes from
     * 
     * @return a constant reference to the map associating each leaf to the 
     *      sample which it comes from
     */
    inline const std::map<CellId, uint16_t>& get_coming_from() const
    {
        return coming_from;
    }

    /**
     * @brief Get the species data
     * 
     * @return a constant reference to the species data.
     */
    inline const std::map<SpeciesId, SpeciesData>& get_species_data() const
    {
        return species_data;
    }
public:

    /**
     * @brief A constant node of the forest
     */
    class const_node 
    {
        DescendantsForest* forest;          //!< A pointer to the forest
        CellId cell_id;                     //!< The cell id of a cell in the forest

        /**
         * @brief A constructor for a constant node
         * 
         * @param forest is the forest of the node
         * @param cell_id is the cell id of a cell in the forest
         */
        const_node(const DescendantsForest* forest, const CellId cell_id);
    public:
        /**
         * @brief Cast to Cell
         * 
         * This method returns a constant reference to a cell. 
         * Notice that the cells in the forest should be modified 
         * exclusively by using `DescendantsForest::node` methods
         * 
         * @return a constant reference to a cell
         */
        inline operator const Cell&() const
        {
            return forest->cells.at(cell_id);
        }

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

        /**
         * @brief Test whether this node is a leaf
         * 
         * @return `true` if and only if this node is a leaf
         */
        inline bool is_leaf() const
        {
            return forest->branches.at(cell_id).size()==0;
        }

        /**
         * @brief Test whether this node is a root
         * 
         * @return `true` if and only if this node is a root
         */
        inline bool is_root() const
        {
            const auto& cell = forest->cells.at(cell_id);

            return cell.get_id() == cell.get_parent_id();
        }

        /**
         * @brief Get the cell id of the node
         * 
         * @return a constant reference to the cell id of the node 
         */
        inline const CellId& get_id() const
        {
            return cell_id;
        }

        /**
         * @brief Get the sample that collected the cell 
         * 
         * @return the sample that collected the identifier 
         *      of this node in the tissue
         * @throw `std::domain_error` when the node is not 
         *      a leaf
         */
        const Simulation::TissueSample& get_sample() const;

        /**
         * @brief Get the node species id
         * 
         * @return a constant reference to the node species id
         */
        inline const SpeciesId& get_species_id() const
        {
            return forest->cells.at(cell_id).get_species_id();
        }

        /**
         * @brief Get the node genotype id
         * 
         * @return a constant reference to the node genotype id
         */
        inline const GenotypeId& get_genotype_id() const
        {
            return forest->species_data.at(get_species_id()).genotype_id;
        }

        /**
         * @brief Get the node genotype name
         * 
         * @return a constant reference to the node genotype name
         */
        inline const std::string& get_genotype_name() const
        {
            return forest->genotype_names.at(get_genotype_id());
        }

        /**
         * @brief Get the node methylation signature
         * 
         * @return a constant reference to the node methylation signature
         */
        inline const MethylationSignature& get_methylation_signature() const
        {
            return forest->species_data.at(get_species_id()).signature;
        }

        /**
         * @brief Get the node forest
         * 
         * @return a constant reference to the node forest
         */
        inline const DescendantsForest& get_forest() const
        {
            return *forest;
        }

        friend class DescendantsForest;
    };

    /**
     * @brief A non-constant node of the forest
     */
    class node : public const_node
    {

        node(DescendantsForest* forest, const CellId cell_id);
    public:
        /**
         * @brief Cast to Cell
         * 
         * @return a non-constant reference to a cell
         */
        inline operator Cell&()
        {
            return forest->cells.at(cell_id);
        }

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

        friend class DescendantsForest;
    };

    /**
     * @brief The empty constructor
     */
    DescendantsForest();

    /**
     * @brief Construct a descendants forest by using a tissue sample of a simulation
     * 
     * This method builds a descendants forest by using clone simulation 
     * pre-sampled cells as leaves.
     * 
     * @param simulation is a simulation
     */
    DescendantsForest(const Simulation::Simulation& simulation);

    /**
     * @brief Construct a descendants forest by using a tissue sample of a simulation
     * 
     * @param simulation is a simulation
     * @param tissue_samples is a list of tissue samples coming from the simulation
     */
    DescendantsForest(const Simulation::Simulation& simulation,
                      const std::list<Simulation::TissueSample>& tissue_samples);

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

    /**
     * @brief Get the number of forest nodes
     * 
     * @return the number of forest nodes 
     */
    inline size_t num_of_nodes() const
    {
        return cells.size();
    }

    /**
     * @brief Get the tissue samples that produced the forest
     * 
     * @return a constant reference to a vector of tissue samples 
     *      that produced the descendants forest
     */
    inline const std::vector<Simulation::TissueSample>& get_samples() const
    {
        return samples;
    }

    /**
     * @brief Get the most recent common ancestors of the forest leaves
     * 
     * @return the most recent common ancestors of the forest leaves
     */
    std::vector<CellId> get_coalescent_cells() const;

    /**
     * @brief Get the most recent common ancestors
     * 
     * @param cell_ids is the list of the cells whose most recent common
     *          ancestors are searched
     * @return the most recent common ancestors of the cells having 
     *          the identifier among those in `cell_ids`
     */
    std::vector<CellId> get_coalescent_cells(const std::list<CellId>& cell_ids) const;

    /**
     * @brief Get the forest for a subset of the tissue samples
     * 
     * @param sample_names are the names of the samples to be considered
     * @return the descendants forest for the tissue samples whose name is 
     *          in `sample_names`
     */
    DescendantsForest get_subforest_for(const std::vector<std::string>& sample_names) const;

    /**
     * @brief Clear the forest
     */
    void clear();
};

}   // Drivers

}   // Races

#endif // __RACES_PHYLOGENETIC_TREE__