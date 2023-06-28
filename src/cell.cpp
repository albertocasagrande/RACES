/**
 * @file cell.cpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Cell representation
 * @version 0.2
 * @date 2023-06-25
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

namespace Races {

uint64_t Cell::counter = 0;

Cell::Cell():
    id(0), parent(0), genotype(NON_DRIVER_GENOTYPE)
{
}

Cell::Cell(const DriverGenotypeId genotype, const CellId parent_id):
    id(++Cell::counter), parent(parent_id), genotype(genotype)
{
}

void swap(Cell& a, Cell &b)
{
    std::swap(a.id,b.id);
    std::swap(a.parent, b.parent);
    std::swap(a.genotype, b.genotype);
}

Cell& Cell::clone(const Cell& cell)
{
    id = cell.id;
    parent = cell.parent;
    genotype = cell.genotype;

    return *this;
}

std::ostream& operator<<(std::ostream& os, const Cell& cell)
{
    os << "Cell{id: "<< cell.get_id() 
       << ", parent_id: "<< cell.get_parent_id()
       << ", driver_genotype: " << cell.get_driver_genotype() 
       << "}";

    return os;
}

CellInTissue::CellInTissue(const AxisValue& x, const AxisValue& y, const AxisValue& z):
    Cell(), PositionInTissue(x, y, z)
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

std::ostream& operator<<(std::ostream& os, const CellInTissue& cell)
{
    os << "Cell{id: "<< cell.get_id() 
       << ", parent_id: "<< cell.get_parent_id()
       << ", driver_genotype: " << cell.get_driver_genotype() 
       << ", position: " << static_cast<PositionInTissue>(cell)
       << "}";

    return os;
}

};
