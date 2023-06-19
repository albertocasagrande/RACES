/**
 * @file tissue.cpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Define tissue class
 * @version 0.1
 * @date 2023-05-30
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

Tissue::Tissue(const std::string name, const unsigned int x_size, const unsigned int  y_size, const unsigned int  z_size):
    name(name)
{
    CellInTissue template_cell(*this, 0, 0, 0);
    std::vector<CellInTissue> z_vector(z_size, template_cell);
    std::vector<std::vector<CellInTissue>> y_vector(y_size, z_vector);
    std::vector<std::vector<std::vector<CellInTissue>>> x_vector(x_size, y_vector);

    unsigned int x{0};
    for (auto x_it=std::begin(x_vector); x_it!=std::end(x_vector); ++x_it, ++x) {
        unsigned int y{0};
        for (auto y_it=std::begin(*x_it); y_it!=std::end(*x_it); ++y_it, ++y) {
        unsigned int z{0};
            for (auto cell_it=std::begin(*y_it); cell_it!=std::end(*y_it); ++cell_it, ++z) {
                cell_it->x = x;
                cell_it->y = y;
                cell_it->z = z;
            }
        }
    }

    std::swap(x_vector, space);
}

Tissue::Tissue(const std::string name, const std::vector<DriverGenotype> genotypes, const unsigned int  x_size, const unsigned int  y_size, const unsigned int  z_size):
    Tissue(name, x_size, y_size, z_size)
{
    for (const auto& genotype: genotypes) {
        add_species(genotype);
    }
}

Tissue::Tissue(const unsigned int  x_size, const unsigned int  y_size, const unsigned int  z_size):
    Tissue("", x_size, y_size, z_size)
{
}

Tissue::Tissue(const std::vector<DriverGenotype> genotypes, const unsigned int  x_size, const unsigned int  y_size, const unsigned int  z_size):
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

Tissue& Tissue::add(const DriverGenotypeId genotype, const PositionInTissue position, const Time time)
{
    auto& cell = (*this)(position);

    if (cell.has_driver_genotype()) {
        throw std::runtime_error("The position has been already taken");
    }

    Species& local_species = get_species(genotype);

    local_species.add(cell, time);

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

CellInTissue& Tissue::operator()(const PositionInTissue& position)
{
    if (!is_valid(position)) {
        throw std::out_of_range("The position does not belong to the tissue");
    }

    return space[position.x][position.y][position.z];
}

const CellInTissue& Tissue::operator()(const PositionInTissue& position) const
{
    if (!is_valid(position)) {
        throw std::out_of_range("The position does not belong to the tissue");
    }

    return space[position.x][position.y][position.z];
}

PositionInTissue Tissue::get_non_driver_in_direction(PositionInTissue position, const Direction& direction) const
{
    if (!is_valid(position)) {
        return position;
    }
 
    PositionDelta delta(direction);
    while ((*this)(position).has_driver_genotype()) {
        position += delta;

        if (!is_valid(position)) {
            return position;
        }
    }

    return position;
}

Tissue& Tissue::push_cells(PositionInTissue position, const Direction& direction)
{

    // search the first empty space in the tissue 
    auto free_cell_pos = get_non_driver_in_direction(position, direction);
    if (!is_valid(free_cell_pos)) {
        // if there is no available space in that direction
        throw std::domain_error("Trying to push outside borders");
    }

    PositionDelta delta(direction);
    PositionInTissue from_pos(free_cell_pos);
    while (free_cell_pos != position) {
        from_pos -= delta;
        Cell& dest_cell = (*this)(free_cell_pos);
        dest_cell.clone((*this)(from_pos));
        Species &species = this->get_species(dest_cell.get_driver_genotype());

        species.update_cell((*this)(free_cell_pos));

        free_cell_pos = from_pos;
    }

    return *this;
}


};