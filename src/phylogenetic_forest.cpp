/**
 * @file phylogenetic_forest.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements classes and function for phylogenetic forests
 * @version 0.14
 * @date 2023-11-15
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

namespace Drivers 
{

DescendantsForest::SpeciesData::SpeciesData(const GenotypeId& genotype_id, const MethylationSignature& signature):
    genotype_id(genotype_id), signature(signature)
{}

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
    DescendantsForest(simulation, simulation.get_tissue_samples())
{}

DescendantsForest::DescendantsForest(const Simulation::Simulation& simulation,
                                     const std::list<Simulation::TissueSample>& tissue_samples)
{
    std::set<GenotypeId> genotype_ids;

    auto species_properties = simulation.tissue().get_species_properties();
    for (const auto& s_propeties: species_properties) {
        const auto& genotype_id = s_propeties.get_genotype_id();
        genotype_ids.insert(genotype_id);

        const auto& signature = s_propeties.get_methylation_signature();
        species_data.insert({s_propeties.get_id(), {genotype_id, signature}});
    }

    for (const auto& genotype_id : genotype_ids) {
        genotype_names[genotype_id] = simulation.find_genotype_name(genotype_id);
    }

    Simulation::BinaryLogger::CellReader reader(simulation.get_logger().get_directory());

    grow_from(tissue_samples, reader);
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

    return forest->samples[found->second];
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

const DescendantsForest::SpeciesData& 
DescendantsForest::get_species_data(const SpeciesId& species_id) const
{
    auto data_it = species_data.find(species_id);
    
    if (data_it == species_data.end()) {
        throw std::runtime_error("Unknown species");
    }

    return data_it->second;
}

std::map<CellId, size_t> 
count_descendants_in_subtree(const DescendantsForest& forest,
                             const std::list<CellId>& leaf_ids)
{
    std::map<CellId, size_t> counter;

    for (const auto& cell_id: leaf_ids) {
        auto found = counter.find(cell_id);
        if (found == counter.end()) {
            counter.insert({cell_id, 1});
        } else {
            ++found->second;
        }

        auto cell_node = forest.get_node(cell_id);
        while (!cell_node.is_root()) {
            cell_node = cell_node.parent();

            ++counter[cell_node.get_id()];
        }
    }
    
    return counter;
}

CellId find_tree_coalescent(const DescendantsForest& forest, const CellId& root,
                            const std::map<CellId, size_t>& counter)
{
    const auto& num_of_leaves = counter.at(root);
    auto cell_node = forest.get_node(root);

    while (true) {
        auto children = cell_node.children();

        if (children.size()>0) {
            size_t i{0};
            auto counter_it = counter.find(children[i].get_id());
            while (counter_it==counter.end()) {
                counter_it = counter.find(children[++i].get_id());
            }

            if (counter_it->second != num_of_leaves) {
                return cell_node.get_id(); 
            }

            cell_node = children[i];
        } else {
            return cell_node.get_id();
        }
    }
}

std::vector<CellId> 
DescendantsForest::get_coalescent_cells(const std::list<CellId>& cell_ids) const
{
    const auto counter = count_descendants_in_subtree(*this, cell_ids);

    std::vector<CellId> coalescent_nodes;
    coalescent_nodes.reserve(roots.size());

    for (const auto& root: roots) {
        if (counter.find(root) != counter.end()) {
            coalescent_nodes.push_back(find_tree_coalescent(*this, root, counter));
        }
    }

    return coalescent_nodes;
}

std::vector<CellId> 
DescendantsForest::get_coalescent_cells() const
{
    std::list<CellId> leaf_ids;

    for (const auto& sample: get_samples()) {
        for (const auto& cell_id: sample.get_cell_ids()) {
            leaf_ids.push_back(cell_id);
        }
    }

    return get_coalescent_cells(leaf_ids);
}

DescendantsForest 
DescendantsForest::get_subforest_for(const std::vector<std::string>& sample_names) const
{
    std::set<std::string> names(sample_names.begin(), sample_names.end());
    std::set<std::string> found_names, missing_names;

    std::list<Simulation::TissueSample> tissue_samples;
    for (const auto& sample : samples) {
        if (names.count(sample.get_name())>0) {
            tissue_samples.push_back(sample);
            found_names.insert(sample.get_name());
        }
    }

    if (names.size() != found_names.size()) {
        std::set_difference(names.begin(), names.end(),
                            found_names.begin(), found_names.end(),
                            std::inserter(missing_names, missing_names.end()));

        std::ostringstream oss;
        oss << "Unknown sample names:";
        std::string sep = " ";
        for (const auto& missing_name : missing_names) {
            oss << sep << missing_name;
            sep = ", ";
        }

        throw std::domain_error(oss.str());
    }

    DescendantsForest subforest = *this;
    subforest.grow_from(tissue_samples, cells);

    return subforest;
}

void DescendantsForest::clear()
{
    cells.clear();
    roots.clear();
    branches.clear();
    samples.clear();
    coming_from.clear();
}

}   // Drivers

}   // Races
