/**
 * @file phylogenetic_forest.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements classes and function for phylogenetic forests
 * @version 0.10
 * @date 2023-11-11
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

#include "phylogenetic_forest.hpp"

#include "simulation.hpp"

namespace Races
{

namespace Drivers 
{

DescendantsForest::const_node::const_node(const DescendantsForest* forest, const CellId cell_id):
    forest(const_cast<DescendantsForest*>(forest)), cell_id(cell_id)
{}

DescendantsForest::const_node DescendantsForest::get_node(const CellId& cell_id) const
{
    if (cells.count(cell_id)==0) {
        throw std::runtime_error("The forest does not contain the cell "
                                 "having the specified identifier");
    }
        
    return DescendantsForest::const_node(this, cell_id);
}

DescendantsForest::node::node(DescendantsForest* forest, const CellId cell_id):
    const_node(forest, cell_id)
{}

DescendantsForest::node DescendantsForest::get_node(const CellId& cell_id)
{
    if (cells.count(cell_id)==0) {
        throw std::runtime_error("The forest does not contain the cell "
                                 "having the specified identifier");
    }

    return DescendantsForest::node(this, cell_id);
}

DescendantsForest::DescendantsForest()
{}

DescendantsForest::DescendantsForest(const Simulation::Simulation& simulation):
    DescendantsForest(simulation, simulation.get_logger().load_samples())
{}

DescendantsForest::DescendantsForest(const Simulation::Simulation& simulation,
                                       const std::list<Simulation::TissueSample>& tissue_samples):
    samples(tissue_samples.begin(),tissue_samples.end())
{
    for (const auto& sample: samples) {
        for (const auto& cell_id: sample.get_cell_ids()) {
            coming_from.insert(std::make_pair(cell_id, &sample));
        }
    }

    Simulation::BinaryLogger::CellReader reader(simulation.get_logger().get_directory());

    auto genotypes = simulation.tissue().get_species_properties();

    std::list<CellId> cell_ids;
    for (auto sample_it = tissue_samples.begin(); sample_it != tissue_samples.end(); ++sample_it) {
        cell_ids.insert(cell_ids.end(), sample_it->get_cell_ids().begin(), 
                        sample_it->get_cell_ids().end());
    }

    grow_forest_from(cell_ids, reader, genotypes);
}

DescendantsForest::const_node DescendantsForest::const_node::parent() const
{
    if (forest==nullptr) {
        throw std::runtime_error("The forest node has not been initialized");
    }

    return DescendantsForest::const_node(forest, forest->cells.at(cell_id).get_parent_id());
}

const Simulation::TissueSample& DescendantsForest::const_node::get_sample() const
{
    if (forest==nullptr) {
        throw std::runtime_error("The forest node has not been initialized");
    }

    auto found = forest->coming_from.find(cell_id);

    if (found == forest->coming_from.end()) {
        throw std::domain_error("The node does not correspond to a sampled cell");
    }

    return *(found->second);
}

std::vector<DescendantsForest::const_node> DescendantsForest::const_node::children() const
{
    if (forest==nullptr) {
        throw std::runtime_error("The forest node has not been initialized");
    }

    std::vector<DescendantsForest::const_node> nodes;

    for (const auto& child_id: forest->branches.at(cell_id)) {
        nodes.push_back(DescendantsForest::const_node(forest, child_id));
    }

    return nodes;
}

GenotypeId DescendantsForest::const_node::get_genotype_id() const
{
    auto genotype_it = forest->genotype_map.find(get_species_id());
    
    if (genotype_it == forest->genotype_map.end()) {
        throw std::runtime_error("unknown species");
    }

    return genotype_it->second;
}

DescendantsForest::node DescendantsForest::node::parent()
{
    if (forest==nullptr) {
       throw std::runtime_error("The forest node has not been initialized");
    }

    return DescendantsForest::node(forest, forest->cells.at(cell_id).get_parent_id());
}

std::vector<DescendantsForest::node> DescendantsForest::node::children()
{
    if (forest==nullptr) {
        throw std::runtime_error("The forest node has not been initialized");
    }

    std::vector<DescendantsForest::node> nodes;

    for (const auto& child_id: forest->branches.at(cell_id)) {
        nodes.push_back(DescendantsForest::node(forest, child_id));
    }

    return nodes;
}

std::vector<DescendantsForest::const_node> DescendantsForest::get_roots() const
{
    std::vector<DescendantsForest::const_node> nodes;

    for (const auto& root_id : roots) {
        nodes.push_back(DescendantsForest::const_node(this, root_id));
    }

    return nodes;
}

std::vector<DescendantsForest::node> DescendantsForest::get_roots()
{
    std::vector<DescendantsForest::node> nodes;

    for (const auto& root_id : roots) {
        nodes.push_back(DescendantsForest::node(this, root_id));
    }

    return nodes;
}

void DescendantsForest::clear()
{
    cells.clear();
    roots.clear();
    branches.clear();
    genotype_map.clear();
}

}   // Drivers

}   // Races
