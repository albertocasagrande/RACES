/**
 * @file tissue.cpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Define tissue class
 * @version 0.15
 * @date 2023-09-17
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

namespace Races
{

namespace Drivers
{

namespace Simulation
{

Tissue::SpeciesView::SpeciesView(const std::vector<Species>& species, const std::vector<size_t>& species_pos):
    species(species), species_pos(species_pos)
{}

Tissue::SpeciesView::const_iterator::const_iterator(const std::vector<Species>& species, const std::vector<size_t>::const_iterator it):
    species(&species), it(std::make_shared<std::vector<size_t>::const_iterator>(it))
{}

Tissue::SpeciesView::const_iterator Tissue::SpeciesView::const_iterator::operator++(int)
{
    Tissue::SpeciesView::const_iterator copy(*this);

    this->operator++();

    return copy;
}

Tissue::SpeciesView::const_iterator Tissue::SpeciesView::const_iterator::operator--(int)
{
    Tissue::SpeciesView::const_iterator copy(*this);

    this->operator--();

    return copy;
}

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

    Species& species = tissue.get_species(cell.get_epigenetic_id());

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
        Species& former_species = tissue.get_species(space_ptr->get_epigenetic_id());
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

Tissue::CellInTissueProxy::operator CellInTissue&()
{
    const auto ptr = tissue.cell_pointer(position);

    if (ptr!=nullptr) {
        return *ptr;
    }

    throw std::runtime_error("Non-driver mutated cell");
}

Tissue::Tissue(const std::string& name, const std::vector<AxisSize>& sizes):
    name(name), dimensions(sizes.size())
{
    AxisSize z_size{1};
    if (dimensions==3) {
        z_size = sizes[2];
    } else {
        if (dimensions!=2) {
            throw std::domain_error("The tissue must be a either a 3D or 2D shape");
        }
    }

    std::vector<CellInTissue *> z_vector(z_size, nullptr);
    std::vector<std::vector<CellInTissue *>> y_vector(sizes[1], z_vector);
    space = std::vector<std::vector<std::vector<CellInTissue *>>>(sizes[0], y_vector);
}

Tissue::Tissue(const std::string& name, const AxisSize x_size,const AxisSize y_size, const AxisSize z_size):
    Tissue(name, {x_size, y_size, z_size})
{
}

Tissue::Tissue(const std::string& name, const AxisSize x_size, const AxisSize y_size):
    Tissue(name, {x_size, y_size})
{
}

Tissue::Tissue(const std::string& name, const std::vector<Genotype>& genotypes,
               const AxisSize  x_size, const AxisSize y_size, const AxisSize z_size):
    Tissue(name, {x_size, y_size, z_size})
{
    for (const auto& genotype: genotypes) {
        add_species(genotype);
    }
}

Tissue::Tissue(const std::string& name, const std::vector<Genotype>& genotypes, 
               const AxisSize x_size, const AxisSize y_size):
    Tissue(name, {x_size, y_size})
{
    for (const auto& genotype: genotypes) {
        add_species(genotype);
    }
}

Tissue::Tissue(const std::vector<AxisSize>& sizes):
    Tissue("", sizes)
{
}

Tissue::Tissue(const AxisSize x_size, const AxisSize y_size, const AxisSize z_size):
    Tissue("", {x_size, y_size, z_size})
{
}

Tissue::Tissue(const std::vector<Genotype>& genotypes,
               const AxisSize x_size, const AxisSize y_size, const AxisSize z_size):
    Tissue("", genotypes, x_size, y_size, z_size)
{
}

Tissue::Tissue(const std::vector<Genotype>& genotypes,
               const AxisSize  x_size, const AxisSize  y_size):
    Tissue("", genotypes, x_size, y_size)
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

std::vector<EpigeneticGenotype> Tissue::get_genotypes() const
{
    std::vector<EpigeneticGenotype> genotypes;

    for (const auto& species: *this) {
        genotypes.push_back(species);
    }

    return genotypes;
}

const Species& Tissue::get_species(const EpigeneticGenotypeId& genotype_id) const
{
    const auto pos_it = pos_map.find(genotype_id);
    if (pos_it == pos_map.end()) {
        throw std::out_of_range("the species has not been included in the tissue");
    }

    return species[pos_it->second];
}

Species& Tissue::get_species(const EpigeneticGenotypeId& genotype_id)
{
    const auto pos_it = pos_map.find(genotype_id);
    if (pos_it == pos_map.end()) {
        throw std::out_of_range("the species has not been included in the tissue");
    }

    return species[pos_it->second];
}

Tissue& Tissue::add_cell(const EpigeneticGenotypeId& genotype_id, const PositionInTissue position)
{
    if (!is_valid(position)) {
        throw std::runtime_error("The position is not in the tissue");
    }
 
    auto*& cell = space[position.x][position.y][position.z];

    if (cell!=nullptr) {
        throw std::runtime_error("The position has been already taken");
    }

    Species& species = get_species(genotype_id);

    cell = species.add(CellInTissue(genotype_id, position));

    return *this;
}

Tissue& Tissue::add_species(const Genotype& genotype)
{
    // check whether the genotype is already in the tissue
    if (genotope_pos.count(genotype.get_id())>0) {
        throw std::runtime_error("Genomic genotype already in the tissue");
    }

    // check whether any of the epigenetic genotypes is already in the tissue
    for (const auto& e_genotype: genotype.epigenetic_genotypes()) {
        if (pos_map.count(e_genotype.get_id())>0) {
            throw std::runtime_error("Epigenetic genotype already in the tissue");
        }
    }

    // insert the genotypes in the tissue 
    auto& pos = genotope_pos[genotype.get_id()];
    for (const auto& e_genotype: genotype.epigenetic_genotypes()) {
        // place the new species at the end of the species vector
        pos.push_back(species.size());

        pos_map[e_genotype.get_id()] = species.size();
        species.push_back(Species(e_genotype));
    }

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

        Species &species = get_species(to_be_moved->get_epigenetic_id());

        species.remove(to_be_moved->get_id());
    }

    cell_pointer(from_position) = nullptr;

    return lost_cell;
}

std::vector<AxisSize> Tissue::size() const
{
    std::vector<AxisSize> sizes(dimensions);

    sizes[0] = static_cast<AxisSize>(space.size());
    sizes[1] = static_cast<AxisSize>(space[0].size());

    if (dimensions==3) {
        sizes[2] = static_cast<AxisSize>(space[0][0].size());
    }
    return sizes;
}

}   // Simulation

}   // Drivers

}   // Races
