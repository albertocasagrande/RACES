/**
 * @file tissue.cpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Define tissue class
 * @version 0.7
 * @date 2023-07-10
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

    Species& species = tissue.get_species(cell.get_genotype_id());

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
        Species& former_species = tissue.get_species(space_ptr->get_genotype_id());
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

Tissue::Tissue(const std::string name, const std::vector<SomaticGenotype> genotypes, const AxisSize  x_size, const AxisSize  y_size, const AxisSize  z_size):
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

Tissue::Tissue(const std::vector<SomaticGenotype> genotypes, const AxisSize  x_size, const AxisSize  y_size, const AxisSize  z_size):
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

Tissue& Tissue::add(const EpigeneticGenotypeId genotype, const PositionInTissue position)
{
    auto*& cell = space[position.x][position.y][position.z];

    if (cell!=nullptr) {
        throw std::runtime_error("The position has been already taken");
    }

    Species& species = get_species(genotype);

    cell = species.add(CellInTissue(genotype, position));

    return *this;
}

Tissue& Tissue::add_species(const SomaticGenotype& somatic_genotype)
{
    // check whether the somatic genotype is already in the tissue
    if (somatic_genotope_pos.count(somatic_genotype.get_id())>0) {
        throw std::runtime_error("Somatic genotype already in the tissue");
    }

    // check whether any of the epigenetic genotypes is already in the tissue
    for (const auto& genotype: somatic_genotype.epigenetic_genotypes()) {
        if (pos_map.count(genotype.get_id())>0) {
            throw std::runtime_error("Epigenetic genotype already in the tissue");
        }
    }

    // insert the genotypes in the tissue 
    auto& pos = somatic_genotope_pos[somatic_genotype.get_id()];
    for (const auto& genotype: somatic_genotype.epigenetic_genotypes()) {
        // place the new species at the end of the species vector
        pos.push_back(species.size());

        pos_map[genotype.get_id()] = species.size();
        species.push_back(genotype);
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

        Species &species = get_species(to_be_moved->get_genotype_id());

        species.remove(to_be_moved->get_id());
    }

    cell_pointer(from_position) = nullptr;

    return lost_cell;
}

};