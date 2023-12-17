/**
 * @file descendant_forest.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements classes and function for descendant forests
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

#include <set>
#include <algorithm>

#include "descendant_forest.hpp"

#include "simulation.hpp"

namespace Races
{

namespace Mutants
{

DescendantsForest::SpeciesData::SpeciesData()
{}

DescendantsForest::SpeciesData::SpeciesData(const MutantId& mutant_id, const MethylationSignature& signature):
    mutant_id(mutant_id), signature(signature)
{}

DescendantsForest::const_node DescendantsForest::get_node(const CellId& cell_id) const
{
    if (cells.count(cell_id)==0) {
        throw std::runtime_error("The forest does not contain the cell "
                                 "having the specified identifier");
    }

    return DescendantsForest::const_node(this, cell_id);
}

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

DescendantsForest::DescendantsForest(const Evolutions::Simulation& simulation):
    DescendantsForest(simulation, simulation.get_tissue_samples())
{}

DescendantsForest::DescendantsForest(const Evolutions::Simulation& simulation,
                                     const std::list<Evolutions::TissueSample>& tissue_samples)
{
    std::set<MutantId> mutant_ids;

    auto species_properties = simulation.tissue().get_species_properties();
    for (const auto& s_propeties: species_properties) {
        const auto& mutant_id = s_propeties.get_mutant_id();
        mutant_ids.insert(mutant_id);

        const auto& signature = s_propeties.get_methylation_signature();
        species_data.insert({s_propeties.get_id(), {mutant_id, signature}});
    }

    for (const auto& mutant_id : mutant_ids) {
        mutant_names[mutant_id] = simulation.find_mutant_name(mutant_id);
    }

    Evolutions::BinaryLogger::CellReader reader(simulation.get_logger().get_directory());

    grow_from(tissue_samples, reader);
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
    std::set<std::string> found_names;

    std::list<Evolutions::TissueSample> tissue_samples;
    for (const auto& sample : samples) {
        if (names.count(sample.get_name())>0) {
            tissue_samples.push_back(sample);
            found_names.insert(sample.get_name());
        }
    }

    if (names.size() != found_names.size()) {
        std::set<std::string> missing_names;
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

bool DescendantsForest::is_leaf(const CellId& cell_id) const
{
    auto branches_it = branches.find(cell_id);

    if (branches_it == branches.end()) {
        throw std::domain_error(std::to_string(cell_id)+" is not a forest cell");
    }
    return branches_it->second.size()==0;
}

void DescendantsForest::clear()
{
    cells.clear();
    roots.clear();
    branches.clear();
    samples.clear();
    coming_from.clear();
}

}   // Mutants

}   // Races
