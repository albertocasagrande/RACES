/**
 * @file tissue.cpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Define tissue class
 * @version 0.5
 * @date 2023-07-05
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

#include <random>

#include "tissue.hpp"
#include "driver_genotype.hpp"

namespace Races {

Tissue::CellInTissueConstantProxy::CellInTissueConstantProxy(const Tissue &tissue, const PositionInTissue position):
    BaseCellInTissueProxy<const Tissue>(tissue, position)
{
}

Tissue::CellInTissueProxy::CellInTissueProxy(Tissue &tissue, const PositionInTissue position):
    BaseCellInTissueProxy<Tissue>(tissue, position)
{}

Tissue::CellInTissueProxy& Tissue::CellInTissueProxy::operator=(const Cell& cell)
{
    // kill the cell in the position
    kill();

    Species& species = tissue.get_species(cell.get_driver_genotype());

    // if the new cell is already in its species
    if (species.contains(cell.get_id())) {
        auto cell_ptr = &(species(cell.get_id()));

        // reset the corresponding tissue position
        tissue.cell_pointer(*cell_ptr) = nullptr;

        // update the cell position in the species
        *cell_ptr = position;
        tissue.cell_pointer(position) = cell_ptr;
    } else { // if, otherwise, the new cell is not in its species
        // add the cell to the species
        tissue.cell_pointer(position) = species.add(CellInTissue(cell, position));
    }

    return *this;
}

void Tissue::CellInTissueProxy::kill()
{
    // if the position already contains a cell with driver mutation
    if (has_driver_mutations()) {
        CellInTissue*& space_ptr = tissue.cell_pointer(position);

        // remove the cell from its species
        Species& former_species = tissue.get_species(space_ptr->get_driver_genotype());
        former_species.remove(space_ptr->get_id());

        space_ptr = nullptr;
    }
}

CellInTissue Tissue::CellInTissueProxy::copy_and_kill()
{
    CellInTissue copy = *this;

    kill();

    return copy;
}

Tissue::Tissue(const std::string name, const AxisSize x_size, const AxisSize  y_size, const AxisSize  z_size):
    name(name)
{
    std::vector<CellInTissue *> z_vector(z_size, nullptr);
    std::vector<std::vector<CellInTissue *>> y_vector(y_size, z_vector);
    std::vector<std::vector<std::vector<CellInTissue *>>> x_vector(x_size, y_vector);

    std::swap(x_vector, space);
}

Tissue::Tissue(const std::string name, const std::vector<DriverGenotype> genotypes, const AxisSize  x_size, const AxisSize  y_size, const AxisSize  z_size):
    Tissue(name, x_size, y_size, z_size)
{
    for (const auto& genotype: genotypes) {
        add_species(genotype);
    }
}

Tissue::Tissue(const AxisSize  x_size, const AxisSize  y_size, const AxisSize  z_size):
    Tissue("", x_size, y_size, z_size)
{
}

Tissue::Tissue(const std::vector<DriverGenotype> genotypes, const AxisSize  x_size, const AxisSize  y_size, const AxisSize  z_size):
    Tissue("", genotypes, x_size, y_size, z_size)
{
}

size_t Tissue::num_of_mutated_cells() const
{
    size_t mutated{0};

    for (const auto& S: species) {
        mutated += S.num_of_cells();
    }

    return mutated;
}

Tissue& Tissue::add(const DriverGenotypeId genotype, const PositionInTissue position, const unsigned int passenger_mutations)
{
    auto*& cell = space[position.x][position.y][position.z];

    if (cell!=nullptr) {
        throw std::runtime_error("The position has been already taken");
    }

    Species& species = get_species(genotype);

    cell = species.add(CellInTissue(genotype, passenger_mutations, position));

    return *this;
}

Tissue& Tissue::add_species(const DriverGenotype& genotype)
{
    // check whether the genotype is already in the tissue
    if (pos_map.count(genotype.get_id())>0) {
        throw std::runtime_error("Driver mutation already in the tissue");
    }

    // place the new species at the end of the species vector
    pos_map[genotype.get_id()] = species.size();
    species.push_back(genotype);

    return *this;
}

Tissue& Tissue::add_driver_somatic_mutation(const DriverGenotypeId& src, const DriverGenotypeId& dst, const Time delay)
{
    somatic_graph.add_edge(pos_map.at(src), pos_map.at(dst), delay);

    return *this;
}

Tissue& Tissue::add_driver_epigenetic_mutation(const DriverGenotypeId& src, const DriverGenotypeId& dst, const double probability)
{
    epigenetic_graph.add_edge(pos_map.at(src), pos_map.at(dst), probability);

    return *this;
}

Tissue::CellInTissueProxy Tissue::operator()(const PositionInTissue& position)
{
    return CellInTissueProxy(*this, position);
}


const Tissue::CellInTissueConstantProxy Tissue::operator()(const PositionInTissue& position) const
{
    return CellInTissueConstantProxy(*this, position);
}

PositionInTissue Tissue::get_non_driver_in_direction(PositionInTissue position, const Direction& direction) const
{
    if (!is_valid(position)) {
        return position;
    }
 
    PositionDelta delta(direction);
    while ((*this)(position).has_driver_mutations()) {
        position += delta;

        if (!is_valid(position)) {
            return position;
        }
    }

    return position;
}

std::list<Cell> Tissue::push_cells(const PositionInTissue from_position, const Direction& direction)
{
    PositionDelta delta(direction);
    PositionInTissue to_position(from_position+delta);
    CellInTissue *to_be_moved = cell_pointer(from_position);
    while (to_be_moved!=nullptr&&is_valid(to_position)) {
        CellInTissue* &dest_ptr = cell_pointer(to_position);

        std::swap(dest_ptr, to_be_moved);

        *dest_ptr = to_position;

        to_position += delta;
    }

    std::list<Cell> lost_cell;

    if (to_be_moved!=nullptr) {
        lost_cell.push_back(*to_be_moved);

        Species &species = get_species(to_be_moved->get_driver_genotype());

        species.remove(to_be_moved->get_id());
    }

    cell_pointer(from_position) = nullptr;

    return lost_cell;
}

};