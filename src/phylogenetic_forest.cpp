/**
 * @file phylogenetic_forest.cpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Implement classes and function for phylogenetic trees
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

#include "phylogenetic_forest.hpp"

namespace Races
{

PhylogeneticForest::const_node::const_node(const PhylogeneticForest* forest, const CellId cell_id):
    forest(forest), cell_id(cell_id)
{}

PhylogeneticForest::node::node(PhylogeneticForest* forest, const CellId cell_id):
    forest(forest), cell_id(cell_id)
{}

PhylogeneticForest::PhylogeneticForest()
{}

PhylogeneticForest::const_node PhylogeneticForest::const_node::parent() const
{
    if (forest==nullptr) {
        std::runtime_error("The forest node has not been initialized");
    }

    return PhylogeneticForest::const_node(forest, forest->cells.at(cell_id).get_parent_id());
}

std::vector<PhylogeneticForest::const_node> PhylogeneticForest::const_node::children() const
{
    if (forest==nullptr) {
        std::runtime_error("The forest node has not been initialized");
    }

    std::vector<PhylogeneticForest::const_node> nodes;

    for (const auto& child_id: forest->branches.at(cell_id)) {
        nodes.push_back(PhylogeneticForest::const_node(forest, child_id));
    }

    return nodes;
}

PhylogeneticForest::const_node PhylogeneticForest::node::parent() const
{
    if (forest==nullptr) {
        std::runtime_error("The forest node has not been initialized");
    }

    return PhylogeneticForest::const_node(forest, forest->cells.at(cell_id).get_parent_id());
}

std::vector<PhylogeneticForest::const_node> PhylogeneticForest::node::children() const
{
    if (forest==nullptr) {
        std::runtime_error("The forest node has not been initialized");
    }

    std::vector<PhylogeneticForest::const_node> nodes;

    for (const auto& child_id: forest->branches.at(cell_id)) {
        nodes.push_back(PhylogeneticForest::const_node(forest, child_id));
    }

    return nodes;
}

PhylogeneticForest::node PhylogeneticForest::node::parent()
{
    if (forest==nullptr) {
        std::runtime_error("The forest node has not been initialized");
    }

    return PhylogeneticForest::node(forest, forest->cells.at(cell_id).get_parent_id());
}

std::vector<PhylogeneticForest::node> PhylogeneticForest::node::children()
{
    if (forest==nullptr) {
        std::runtime_error("The forest node has not been initialized");
    }

    std::vector<PhylogeneticForest::node> nodes;

    for (const auto& child_id: forest->branches.at(cell_id)) {
        nodes.push_back(PhylogeneticForest::node(forest, child_id));
    }

    return nodes;
}


std::vector<PhylogeneticForest::const_node> PhylogeneticForest::get_roots() const
{
    std::vector<PhylogeneticForest::const_node> nodes;

    for (const auto& root_id : roots) {
        nodes.push_back(PhylogeneticForest::const_node(this, root_id));
    }

    return nodes;
}

std::vector<PhylogeneticForest::node> PhylogeneticForest::get_roots()
{
    std::vector<PhylogeneticForest::node> nodes;

    for (const auto& root_id : roots) {
        nodes.push_back(PhylogeneticForest::node(this, root_id));
    }

    return nodes;
}

}
