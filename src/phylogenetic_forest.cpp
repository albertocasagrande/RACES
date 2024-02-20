/**
 * @file phylogenetic_forest.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements classes and function for phylogenetic forests
 * @version 0.8
 * @date 2024-02-20
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

PhylogeneticForest::const_node::const_node(const PhylogeneticForest* forest, const Mutants::CellId cell_id):
    Mutants::DescendantsForest::_const_node<PhylogeneticForest>(forest, cell_id)
{}

PhylogeneticForest::node::node(PhylogeneticForest* forest, const Mutants::CellId cell_id):
    Mutants::DescendantsForest::_node<PhylogeneticForest>(forest, cell_id)
{}

void PhylogeneticForest::node::add_new_mutation(const SNV& snv)
{
    get_forest().novel_mutations[cell_id].SNVs.insert(snv);
    get_forest().SNV_first_cells[snv].insert(cell_id);
}

void PhylogeneticForest::node::add_new_mutation(const CopyNumberAlteration& cna)
{
    get_forest().novel_mutations[cell_id].CNAs.insert(cna);
    get_forest().CNA_first_cells[cna].insert(cell_id);
}

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

template<typename MUTATION>
std::map<MUTATION, std::set<Mutants::CellId>>
filter_by_cells_in(const std::map<MUTATION, std::set<Mutants::CellId>>& mutation_map,
                   const Mutants::DescendantsForest& forest)
{
    std::map<MUTATION, std::set<Mutants::CellId>> filtered_map;

    for (const auto& [mutation, cell_ids]: mutation_map) {
        auto& in_forest = filtered_map[mutation];
        for (const auto& cell_id : cell_ids) {
            if (forest.get_cells().count(cell_id)>0) {
                in_forest.insert(cell_id);
            }
        }
    }

    return filtered_map;
}

PhylogeneticForest
PhylogeneticForest::get_subforest_for(const std::vector<std::string>& sample_names) const
{
    PhylogeneticForest forest;

    static_cast<Mutants::DescendantsForest&>(forest) = Mutants::DescendantsForest::get_subforest_for(sample_names);

    for (const auto& [cell_id, cell]: forest.get_cells()) {
        forest.novel_mutations[cell_id] = novel_mutations.at(cell_id);
        if (forest.is_leaf(cell_id)) {
            forest.leaves_mutations[cell_id] = leaves_mutations.at(cell_id);
        }
    }

    forest.SNV_first_cells = filter_by_cells_in(SNV_first_cells, forest);
    forest.CNA_first_cells = filter_by_cells_in(CNA_first_cells, forest);

    return forest;
}

const CellGenomeMutations&
PhylogeneticForest::get_leaf_mutations(const Mutants::CellId& cell_id) const
{
    if (is_leaf(cell_id)) {
        return *(leaves_mutations.at(cell_id));
    }

    throw std::domain_error(std::to_string(cell_id)+" is not a leaf of the forest");
}

std::list<SampleGenomeMutations> PhylogeneticForest::get_sample_mutations_list() const
{
    using namespace Mutants::Evolutions;

    std::list<SampleGenomeMutations> sample_mutations_list;
    std::map<TissueSampleId, SampleGenomeMutations*> sample_mutation_map;
    for (const auto& sample: get_samples()) {
        sample_mutations_list.emplace_back(sample.get_name(), germline_mutations);
        sample_mutation_map.insert({sample.get_id(), &(sample_mutations_list.back())});
    }

    for (const auto& [cell_id, mutations_ptr] : leaves_mutations) {
        const auto& sample = get_samples()[get_coming_from().at(cell_id)];
        auto& sample_genomic_mutations = *(sample_mutation_map.at(sample.get_id()));
        sample_genomic_mutations.mutations.push_back(mutations_ptr);
    }

    return sample_mutations_list;
}

size_t
find_sample_index(const std::vector<Races::Mutants::Evolutions::TissueSample>& sample_list,
                 const std::string& sample_name)
{
    size_t idx=0;
    for (const auto& sample: sample_list) {
        if (sample.get_name() == sample_name) {
            return idx;
        }
        ++idx;
    }

    throw std::runtime_error("Unknown sample \""+sample_name+"\".");
}

SampleGenomeMutations PhylogeneticForest::get_sample_mutations(const std::string& sample_name) const
{
    using namespace Mutants::Evolutions;

    SampleGenomeMutations sample_mutations(sample_name, germline_mutations);

    const size_t sample_idx = find_sample_index(get_samples(), sample_name);

    for (const auto& cell_id : get_samples()[sample_idx].get_cell_ids()) {
        sample_mutations.mutations.push_back(leaves_mutations.at(cell_id));
    }

    return sample_mutations;
}

void PhylogeneticForest::clear()
{
    DescendantsForest::clear();
    leaves_mutations.clear();
    novel_mutations.clear();
    SNV_first_cells.clear();
    CNA_first_cells.clear();
}

}   // Mutants

}   // Races
