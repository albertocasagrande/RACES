/**
 * @file species.cpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Implements species representation methods
 * @version 0.9
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

#include "species.hpp"

namespace Races 
{

namespace Drivers
{

namespace Simulation
{

Species::Species():
    EpigeneticGenotype()
{}

void Species::reset()
{
    for (auto& cell: cells) {
        delete cell;
    }

    cells.clear();
    pos_map.clear();
}

Species::Species(const EpigeneticGenotype& genotype):
    EpigeneticGenotype(genotype), last_insertion_time(0)
{}

Species::Species(const Species& orig):
    EpigeneticGenotype(orig), last_insertion_time(0)
{
    for (const auto& cell: orig.cells) {
        this->add(*cell);
    }
}

Species& Species::operator=(const Species& orig)
{
    reset();

    static_cast<EpigeneticGenotype&>(*this) = static_cast<EpigeneticGenotype>(orig);

    for (const auto& cell: orig.cells) {
        this->add(*cell);
    }

    return *this;
}

Species& Species::operator=(Species&& orig)
{
    reset();

    std::swap(*this, orig);

    return *this;
}

void Species::remove(const CellId& cell_id)
{
    // find the cell to be removed
    auto pos = pos_map.at(cell_id);

    // delete it from the position map
    pos_map.erase(cell_id);

    // delete the cell
    delete *pos;

    // remove the cell pointer
    cells.erase(pos);
}

CellInTissue* Species::add(CellInTissue* cell)
{   
    // update the genotype id
    cell->epigenetic_id = get_id();

    // update `last_insertion_time`
    last_insertion_time = cell->get_birth_time();

    cells.push_back(cell);

    // find a position for the new cell
    auto final_pos = cells.end();

    // set it position in the position map
    pos_map[cell->get_id()] = --final_pos;

    return cell;
}

CellInTissue* Species::add(CellInTissue&& cell)
{
    cell.epigenetic_id = get_id();

    return add(new CellInTissue(cell));
}

CellInTissue* Species::add(CellInTissue& cell)
{
    cell.epigenetic_id = get_id();

    return add(new CellInTissue(cell));
}

Species::~Species()
{
    reset();
}

Species::const_iterator Species::begin() const
{
    return const_iterator(std::begin(cells));
}

Species::const_iterator Species::end() const
{
    return const_iterator(std::end(cells));
}

Species::const_iterator::const_iterator(const std::list<CellInTissue*>::const_iterator it): it(it)
{}

Species::const_iterator::const_iterator(): it()
{}

Species::const_iterator Species::const_iterator::operator++(int) 
{
    const_iterator copy = *this;
    ++(*this);
    return copy;
}

Species::const_iterator Species::const_iterator::operator--(int) 
{
    const_iterator copy = *this;
    --(*this);
    return copy;
}


void swap(Species& a, Species& b)
{
    std::swap(static_cast<EpigeneticGenotype&>(a),
              static_cast<EpigeneticGenotype&>(b));
    std::swap(a.cells,b.cells);
    std::swap(a.pos_map,b.pos_map);
    std::swap(a.last_insertion_time,b.last_insertion_time);
}

}   // Simulation

}   // Drivers

}   // Races

namespace std
{

std::ostream& operator<<(std::ostream& out, const Races::Drivers::Simulation::Species& species)
{
    out << "{genotype: " << static_cast<Races::Drivers::EpigeneticGenotype>(species) 
        << ", cells: {";
    std::string sep{""};
    for (const auto& cell: species) {
        out << sep << cell;
        sep = ",";
    }
    out <<"}";

    return out;
}

}   // std
