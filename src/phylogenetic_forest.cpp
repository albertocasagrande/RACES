/**
 * @file phylogenetic_forest.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements classes and function for phylogenetic forests
 * @version 1.10
 * @date 2025-09-24
 *
 * @copyright Copyright (c) 2023-2025
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

namespace RACES
{

namespace Mutations
{

PhylogeneticForest::SampleStatistics::SampleStatistics():
    total_allelic_size{0}, number_of_cells{0}
{}

PhylogeneticForest::PhylogeneticForest():
    DescendantsForest()
{}

PhylogeneticForest::const_node::const_node(const PhylogeneticForest* forest, const Mutants::CellId cell_id):
    Mutants::DescendantsForest::_const_node<PhylogeneticForest>(forest, cell_id)
{}

PhylogeneticForest::node::node(PhylogeneticForest* forest, const Mutants::CellId cell_id):
    Mutants::DescendantsForest::_node<PhylogeneticForest>(forest, cell_id)
{}

const MutationList& PhylogeneticForest::const_node::pre_neoplastic_mutations() const
{
    auto found = forest->get_pre_neoplastic_mutations().find(get_id());

    if (found != forest->get_pre_neoplastic_mutations().end()) {
        return found->second;
    }

    throw std::runtime_error("const_node::pre_neoplastic_mutations(): The "
                             "node is not a forest root.");
}

PhylogeneticForest::GenomeMutationTour::GenomeMutationTour(const PhylogeneticForest *forest,
                                                           const bool only_leaves,
                                                           const bool with_pre_neoplastic,
                                                           const bool with_germinal):
    forest{forest}, only_leaves{only_leaves}, with_pre_neoplastic{with_pre_neoplastic},
    with_germinal{with_germinal}
{}

PhylogeneticForest::GenomeMutationTour::GenomeMutationTour():
    forest{nullptr}, with_pre_neoplastic{true}, with_germinal{false}
{}

PhylogeneticForest::GenomeMutationTour::GenomeMutationTour(const PhylogeneticForest& forest,
                                                           const bool only_leaves,
                                                           const bool with_pre_neoplastic,
                                                           const bool with_germinal):
    GenomeMutationTour(&forest, only_leaves, with_pre_neoplastic, with_germinal)
{}

PhylogeneticForest::GenomeMutationTour::const_iterator::const_iterator(const PhylogeneticForest* forest,
                                                                       const bool& only_leaves,
                                                                       const bool& with_pre_neoplastic,
                                                                       const bool& with_germinal,
                                                                       const bool& begin):
    forest{forest}, only_leaves{only_leaves}, tour_end{false}
{
    if (forest != nullptr && begin) {
        auto forest_roots = forest->get_roots();
        for (auto root_it = forest_roots.rbegin();
                root_it != forest_roots.rend(); ++root_it) {
            const auto& germline_mutations = forest->get_germline_mutations();
            GenomeMutations mutations;

            if (with_germinal) {
                mutations = germline_mutations;
            } else {
                mutations = germline_mutations.copy_structure();
            }
            mutations.apply(root_it->arising_mutations());

            if (with_pre_neoplastic) {
                mutations.apply(root_it->pre_neoplastic_mutations());
            }

            iterator_stack.emplace(static_cast<const Mutants::Cell&>(*root_it),
                                   std::move(mutations));
        }

        std::swap(node_mutations, iterator_stack.top());

        iterator_stack.pop();

        if (only_leaves) {
            this->operator++();
        }
    } else {
        tour_end = true;
    }
}

PhylogeneticForest::GenomeMutationTour::const_iterator::const_iterator():
    forest{nullptr}, only_leaves{false}, tour_end{true}
{}

PhylogeneticForest::GenomeMutationTour::const_iterator&
PhylogeneticForest::GenomeMutationTour::const_iterator::operator++()
{
    const_node node{forest, node_mutations.get_id()};

    if (node.is_leaf()) {
        if (iterator_stack.empty()) {
            tour_end = true;

            return *this;
        }

        // take a new cell genome mutations object from the stack
        std::swap(node_mutations, iterator_stack.top());
        iterator_stack.pop();

        if (!only_leaves) {
            return *this;
        }

        // update the node
        node = const_node{forest, node_mutations.get_id()};
    }

    bool next_node_found = node.is_leaf();
    while (!next_node_found) {
        auto children = node.children();

        // place all children's cell genome mutations, but the first one, into the stack
        for (auto child_it = children.rbegin();
                child_it != children.rend()-1; ++child_it) {

            CellGenomeMutations child_mutations{static_cast<const Mutants::Cell&>(*child_it),
                                                node_mutations};
            child_mutations.apply(child_it->arising_mutations());

            iterator_stack.emplace(std::move(child_mutations));
        }

        // apply the first children mutations to the current cell genome mutations
        std::swap(node, children.front());
        node_mutations.apply(node.arising_mutations());

        next_node_found = node.is_leaf() || !only_leaves;
    }

    // update the cell
    static_cast<Mutants::Cell&>(node_mutations) = static_cast<const Mutants::Cell&>(node);

    return *this;
}

bool PhylogeneticForest::GenomeMutationTour::const_iterator::operator==(const const_iterator& rhs) const
{
    if (this->forest != rhs.forest) {
        return false;
    }

    if (this->iterator_stack.size() != rhs.iterator_stack.size()) {
        return false;
    }

    if (this->iterator_stack.size()==0) {
        if (this->tour_end != rhs.tour_end) {
            return false;
        }

        if (this->tour_end) {
            return true;
        }
    }

    return this->iterator_stack.top() == rhs.iterator_stack.top();
}

PhylogeneticForest::const_node PhylogeneticForest::get_node(const Mutants::CellId& cell_id) const
{
    if (get_cells().count(cell_id)==0) {
        throw std::domain_error("The cell " + std::to_string(cell_id)
                                + " is not in the forest.");
    }

    return PhylogeneticForest::const_node(this, cell_id);
}

PhylogeneticForest::node PhylogeneticForest::get_node(const Mutants::CellId& cell_id)
{
    if (get_cells().count(cell_id)==0) {
        throw std::domain_error("The cell " + std::to_string(cell_id)
                                + " is not in the forest.");
    }

    return PhylogeneticForest::node(this, cell_id);
}

std::vector<PhylogeneticForest::const_node> PhylogeneticForest::get_roots() const
{
    std::vector<PhylogeneticForest::const_node> nodes;

    for (const auto& root_id : get_root_cell_ids()) {
        nodes.push_back(PhylogeneticForest::const_node(this, root_id));
    }

    return nodes;
}

std::vector<PhylogeneticForest::node> PhylogeneticForest::get_roots()
{
    std::vector<PhylogeneticForest::node> nodes;

    for (const auto& root_id : get_root_cell_ids()) {
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
        forest.arising_mutations.emplace(cell_id, arising_mutations.at(cell_id));
    }

    for (const auto& root_id : forest.get_root_cell_ids()) {
        forest.pre_neoplastic_mutations.emplace(root_id, pre_neoplastic_mutations.at(root_id));
    }

    const auto samples = forest.get_samples();
    const auto sample_id = samples.front().get_id();

    forest.sample_statistics.emplace(sample_id, sample_statistics.at(sample_id));

    forest.SID_first_cells = filter_by_cells_in(SID_first_cells, forest);
    forest.CNA_first_cells = filter_by_cells_in(CNA_first_cells, forest);

    return forest;
}

CellGenomeMutations PhylogeneticForest::get_cell_mutations(const Mutants::CellId& cell_id,
                                                           const bool& with_pre_neoplastic,
                                                           const bool& with_germinal) const
{
    if (arising_mutations.count(cell_id) == 0) {
        throw std::domain_error("The cell " + std::to_string(cell_id)
                                + " is not in the forest.");
    }

    auto node = get_node(cell_id);

    const Mutants::Cell cell = static_cast<const Mutants::Cell&>(node);

    // find the branch from root to the node
    std::stack<Mutants::CellId> cell_id_stack;
    do {
        cell_id_stack.push(node.get_id());

        node = node.parent();
    } while (!node.is_root());

    // initialize the node genome mutations to the wild-type with or
    // without germline
    CellGenomeMutations node_mutations;
    if (with_germinal) {
        node_mutations = CellGenomeMutations{cell, *germline_mutations};
    } else {
        node_mutations = CellGenomeMutations{cell, germline_mutations->copy_structure()};
    }

    if (with_pre_neoplastic) {
        node_mutations.apply(pre_neoplastic_mutations.at(node.get_id()));
    }
    node_mutations.apply(arising_mutations.at(node.get_id()));

    // browse the branch down to the node and apply the arising mutations
    while (!cell_id_stack.empty()) {
        node_mutations.apply(arising_mutations.at(cell_id_stack.top()));

        cell_id_stack.pop();
    }

    return node_mutations;
}

std::list<SampleGenomeMutations>
PhylogeneticForest::get_sample_mutations_list(const bool& with_germinal) const
{
    // TODO: For memory reasons, we don't want to store a list of the genome
    // mutations of all samples. This method will became deprecated, and it
    // will be replaced by a less demanding alternative.

    using namespace Mutants::Evolutions;

    std::list<SampleGenomeMutations> sample_mutations_list;
    std::map<TissueSampleId, SampleGenomeMutations*> sample_mutation_map;
    for (const auto& sample: get_samples()) {
        sample_mutations_list.emplace_back(sample.get_name(), germline_mutations);
        sample_mutation_map.insert({sample.get_id(), &(sample_mutations_list.back())});
    }

    for (const auto& leaf_mutations : get_leaf_mutation_tour(with_germinal)) {
        const auto& sample = get_samples()[get_coming_from().at(leaf_mutations.get_id())];
        auto& sample_genomic_mutations = *(sample_mutation_map.at(sample.get_id()));
        sample_genomic_mutations.mutations.push_back(std::make_shared<CellGenomeMutations>(leaf_mutations));
    }

    return sample_mutations_list;
}

size_t
find_sample_index(const std::vector<RACES::Mutants::Evolutions::TissueSample>& sample_list,
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

SampleGenomeMutations
PhylogeneticForest::get_sample_mutations(const std::string& sample_name,
                                         const bool& with_germinal) const
{
    // TODO: For memory reasons, we don't want to store a list of the genome
    // mutations of all samples. This method will became deprecated, and it
    // will be replaced by a less demanding alternative.

    using namespace Mutants::Evolutions;

    SampleGenomeMutations sample_mutations(sample_name, germline_mutations);

    const size_t sample_idx = find_sample_index(get_samples(), sample_name);

    for (const auto& cell_id : get_samples()[sample_idx].get_cell_ids()) {
        auto cell_mutations = std::make_shared<CellGenomeMutations>();
        *cell_mutations = get_cell_mutations(cell_id, with_germinal);
        sample_mutations.mutations.push_back(cell_mutations);
    }

    return sample_mutations;
}

SampleGenomeMutations PhylogeneticForest::get_normal_sample(const std::string& name,
                                                            const bool& with_pre_neoplastic) const
{
    SampleGenomeMutations sample(name, germline_mutations);

    if (with_pre_neoplastic) {
        for (auto& [cell_id, cell_genome] : get_normal_genomes(with_pre_neoplastic)) {
            auto genome = std::make_shared<CellGenomeMutations>();

            std::swap(cell_genome, *genome);

            sample.mutations.push_back(genome);
        }
    } else {
        sample.mutations.push_back(std::make_shared<CellGenomeMutations>(germline_mutations->copy_structure()));
    }

    return sample;
}

std::map<Mutants::CellId, CellGenomeMutations>
PhylogeneticForest::get_normal_genomes(const bool& with_pre_neoplastic) const
{
    std::map<Mutants::CellId, CellGenomeMutations> normal_genomes;

    if (with_pre_neoplastic) {
        for (const auto& [cell_id, pnp_mutations] : get_pre_neoplastic_mutations()) {
            auto cell_genome = CellGenomeMutations(germline_mutations->copy_structure());
            cell_genome.apply(pnp_mutations);

            normal_genomes.insert({cell_id, std::move(cell_genome)});
        }
    } else {
        for (const auto& root : get_roots()) {
            normal_genomes.insert({root.get_id(),
                                   CellGenomeMutations(germline_mutations->copy_structure())});
        }
    }

    return normal_genomes;
}

std::map<ChromosomeId, std::set<ChrPosition>>
PhylogeneticForest::get_CNA_break_points() const
{
    std::map<ChromosomeId, std::set<ChrPosition>> b_points;

    for (const auto& leaf_mutations : get_leaf_mutation_tour()) {
        for (const auto& [chr_id, cb_points] : leaf_mutations.get_CNA_break_points()) {
            auto& chr_b_points = b_points[chr_id];
            for (const auto& cb_point : cb_points) {
                chr_b_points.insert(cb_point);
            }
        }
    }

    return b_points;
}

void update_allelic_count_on(PhylogeneticForest::AllelicCount& allelic_count,
                             const GenomeMutations::AllelicMap& allelic_map)
{
    std::greater<uint32_t> cmp;

    for (const auto& [chr_id, chr_allelic_map] : allelic_map) {
        auto& chr_allelic_count = allelic_count[chr_id];
        for (auto [pos, atype] : chr_allelic_map) {
            std::sort(atype.begin(), atype.end(), cmp);
            ++(chr_allelic_count[pos][atype]);
        }
    }
}

PhylogeneticForest::AllelicCount
PhylogeneticForest::get_allelic_count(const size_t& min_allelic_size) const
{
    auto b_points = get_CNA_break_points();

    PhylogeneticForest::AllelicCount allelic_count;

    for (const auto& leaf_mutations : get_leaf_mutation_tour()) {
        auto allelic_map = leaf_mutations.get_allelic_map(b_points, min_allelic_size);

        update_allelic_count_on(allelic_count, allelic_map);
    }

    return allelic_count;
}

PhylogeneticForest::AllelicCount
PhylogeneticForest::get_allelic_count(const std::list<Mutants::CellId>& cell_ids,
                                      const size_t& min_allelic_size) const
{
    auto b_points = get_CNA_break_points();

    PhylogeneticForest::AllelicCount allelic_count;

    for (const auto& cell_id : cell_ids) {
        const auto cell_mutations = get_cell_mutations(cell_id);

        auto allelic_map = cell_mutations.get_allelic_map(b_points, min_allelic_size);

        update_allelic_count_on(allelic_count, allelic_map);
    }

    return allelic_count;
}


PhylogeneticForest::AllelicCount
PhylogeneticForest::get_allelic_count(const std::string& sample_name,
                                      const size_t& min_allelic_size) const
{
    for (const auto& sample : get_samples()) {
        if (sample.get_name() == sample_name) {
            return get_allelic_count(sample.get_cell_ids(), min_allelic_size);
        }
    }

    throw std::domain_error("\"" + sample_name +"\" is not a sample of the "
                            + "phylogenetic forest");
}

void PhylogeneticForest::clear()
{
    arising_mutations.clear();
    SID_first_cells.clear();
    CNA_first_cells.clear();
    sample_statistics.clear();
    DescendantsForest::clear();
}

}   // Mutants

}   // RACES
