/**
 * @file phylogenetic_forest.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements classes and function for phylogenetic forests
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

#include <set>
#include <algorithm>

#include "phylogenetic_forest.hpp"

#include "simulation.hpp"

namespace Races
{

namespace Mutations
{

PhylogeneticForest::PhylogeneticForest():
    DescendantsForest()
{}

PhylogeneticForest::PhylogeneticForest(const Mutants::DescendantsForest& descendant_forest,
                                       std::map<Mutants::CellId, CellGenomeMutations>&& cell_mutations):
    Mutants::DescendantsForest(descendant_forest)
{
    std::swap(this->cell_mutations, cell_mutations);
}

PhylogeneticForest::PhylogeneticForest(const Mutants::DescendantsForest& descendant_forest,
                                       const std::map<Mutants::CellId, CellGenomeMutations>& cell_mutations):
    Mutants::DescendantsForest(descendant_forest), cell_mutations(cell_mutations)
{}

PhylogeneticForest::const_node::const_node(const PhylogeneticForest* forest, const Mutants::CellId cell_id):
    Mutants::DescendantsForest::_const_node<PhylogeneticForest>(forest, cell_id)
{}

PhylogeneticForest::node::node(PhylogeneticForest* forest, const Mutants::CellId cell_id):
    Mutants::DescendantsForest::_node<PhylogeneticForest>(forest, cell_id),
    const_node(forest, cell_id)
{}


PhylogeneticForest::const_node PhylogeneticForest::get_node(const Mutants::CellId& cell_id) const
{
    if (get_cells().count(cell_id)==0) {
        throw std::runtime_error("The forest does not contain the cell "
                                 "having the specified identifier");
    }

    return PhylogeneticForest::const_node(this, cell_id);
}

PhylogeneticForest::node PhylogeneticForest::get_node(const Mutants::CellId& cell_id)
{
    if (get_cells().count(cell_id)==0) {
        throw std::runtime_error("The forest does not contain the cell "
                                 "having the specified identifier");
    }

    return PhylogeneticForest::node(this, cell_id);
}

std::vector<PhylogeneticForest::const_node> PhylogeneticForest::get_roots() const
{
    std::vector<PhylogeneticForest::const_node> nodes;

    for (const auto& root_id : get_root_cells()) {
        nodes.push_back(PhylogeneticForest::const_node(this, root_id));
    }

    return nodes;
}

std::vector<PhylogeneticForest::node> PhylogeneticForest::get_roots()
{
    std::vector<PhylogeneticForest::node> nodes;

    for (const auto& root_id : get_root_cells()) {
        nodes.push_back(PhylogeneticForest::node(this, root_id));
    }

    return nodes;
}

PhylogeneticForest
PhylogeneticForest::get_subforest_for(const std::vector<std::string>& sample_names) const
{
    PhylogeneticForest forest;

    static_cast<Mutants::DescendantsForest&>(forest) = Mutants::DescendantsForest::get_subforest_for(sample_names);

    for (const auto& [cell_id, cell]: forest.get_cells()) {
        forest.cell_mutations[cell_id] = cell_mutations.at(cell_id);
    }

    return forest;
}

void PhylogeneticForest::clear()
{
    DescendantsForest::clear();
    cell_mutations.clear();
}

}   // Mutants

}   // Races
