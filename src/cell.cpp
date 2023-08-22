/**
 * @file cell.cpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Cell representation
 * @version 0.10
 * @date 2023-08-22
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

#include "cell.hpp"
#include "species.hpp"
#include "tissue.hpp"
#include "driver_genotype.hpp"

namespace Races
{

namespace Drivers
{

uint64_t Cell::counter = 0;

Cell::Cell():
    id(0), parent(0), birth_time(0), genotype(NON_DRIVER_GENOTYPE)
{
}

Cell::Cell(const EpigeneticGenotypeId genotype):
    id(Cell::counter), parent(Cell::counter), birth_time(0), genotype(genotype)
{
    ++Cell::counter;
}

Cell::Cell(const EpigeneticGenotypeId genotype, const CellId parent_id):
    Cell(genotype,parent_id,0)
{}

Cell::Cell(const EpigeneticGenotypeId genotype, const CellId parent_id, const Time birth_time):
    id(Cell::counter), parent(parent_id), birth_time(birth_time), genotype(genotype)
{
    ++Cell::counter;
}

void swap(Cell& a, Cell &b)
{
    std::swap(a.id,b.id);
    std::swap(a.parent, b.parent);
    std::swap(a.birth_time, b.birth_time);
    std::swap(a.genotype, b.genotype);
}

namespace Simulation
{


CellInTissue::CellInTissue():
    Cell(), PositionInTissue()
{}

CellInTissue::CellInTissue(const EpigeneticGenotypeId genotype, const PositionInTissue& position):
    Cell(genotype), PositionInTissue(position)
{
}

CellInTissue::CellInTissue(const PositionInTissue& position):
    Cell(), PositionInTissue(position)
{
}

CellInTissue::CellInTissue(const AxisValue& x, const AxisValue& y, const AxisValue& z):
    Cell(), PositionInTissue(x, y, z)
{
}

CellInTissue::CellInTissue(const Cell& cell, const PositionInTissue& position):
    Cell(cell), PositionInTissue(position)
{
}

CellInTissue& CellInTissue::operator=(const Cell& cell)
{
    Cell::operator=(cell);

    return *this;
}


CellInTissue& CellInTissue::operator=(const PositionInTissue& position)
{
    x = position.x;
    y = position.y;
    z = position.z;

    return *this;
}

}   // Simulation

}   // Drivers

}   // Races

namespace std
{

std::ostream& operator<<(std::ostream& os, const Races::Drivers::Cell& cell)
{
    os << "Cell{id: "<< cell.get_id()
       << ", parent_id: "<< cell.get_parent_id()
       << ", birth_time: "<< cell.get_birth_time()
       << ", driver_genotype: " << cell.get_genotype_id()
       << "}";

    return os;
}

std::ostream& operator<<(std::ostream& os, const Races::Drivers::Simulation::CellInTissue& cell)
{
    os << "Cell{id: "<< cell.get_id()
       << ", parent_id: "<< cell.get_parent_id()
       << ", birth_time: "<< cell.get_birth_time()
       << ", driver_genotype: " << cell.get_genotype_id()
       << ", position: " << static_cast<Races::Drivers::Simulation::PositionInTissue>(cell)
       << "}";

    return os;
}

}  // std