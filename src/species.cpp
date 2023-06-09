/**
 * @file species.cpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Cell representation
 * @version 0.6
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

#include "species.hpp"

namespace Races {

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
    EpigeneticGenotype(genotype)
{}

Species::Species(const Species& orig):
    EpigeneticGenotype(orig)
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
    const size_t pos = pos_map.at(cell_id);

    // delete it from the position map
    pos_map.erase(cell_id);

    // swap the to-be-removed cell and the last one
    // in the cell vector
    const size_t last_pos = cells.size()-1;
    if (last_pos>0) {
        if (pos != last_pos) {
            std::swap(cells[pos], cells[last_pos]);

            // update the former last cell position
            pos_map[cells[pos]->get_id()] = pos;
        }
    }

    delete cells[last_pos];

    // remove the cell to be removed 
    cells.pop_back();
}

CellInTissue* Species::add(CellInTissue* cell)
{
    // find a position for the new cell
    const size_t new_pos = cells.size();

    // set it position in the position map
    pos_map[cell->get_id()] = new_pos;
    
    cell->genotype = get_id();

    cells.push_back(cell);

    return cell;
}

CellInTissue* Species::add(CellInTissue&& cell)
{
    cell.genotype = get_id();

    return add(new CellInTissue(cell));
}

CellInTissue* Species::add(CellInTissue& cell)
{
    cell.genotype = get_id();

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

Species::const_iterator::const_iterator(const std::vector<CellInTissue*>::const_iterator it): it(it)
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

std::ostream& operator<<(std::ostream& out, const Species& species)
{
    out << "{genotype: " << static_cast<EpigeneticGenotype>(species) 
        << ", cells: {";
    std::string sep{""};
    for (const auto& cell: species) {
        out << sep << cell;
        sep = ",";
    }
    out <<"}";

    return out;
}

}

namespace std
{
void swap(Races::Species& a, Races::Species& b)
{
    std::swap(static_cast<Races::EpigeneticGenotype&>(a),static_cast<Races::EpigeneticGenotype&>(b));
    std::swap(a.cells,b.cells);
    std::swap(a.pos_map,b.pos_map);
}
}