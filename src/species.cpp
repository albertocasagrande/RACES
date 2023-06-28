/**
 * @file species.cpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Cell representation
 * @version 0.2
 * @date 2023-06-28
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

Species::Species(const std::string& name,
                 const std::map<CellEventType, double>& rates,
                 const bool methylated):
    DriverGenotype(name, rates, methylated)
{}

Species::Species(const Species& species):
    DriverGenotype(species), cells(species.cells), pos_map(species.pos_map)
{}

Species::Species(const DriverGenotype& driver):
    DriverGenotype(driver)
{}

Species& Species::remove(const CellId& cell_id)
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

    cells[last_pos]->genotype = NON_DRIVER_GENOTYPE;

    // remove the cell to be removed 
    cells.pop_back();

    return *this;
}

Species& Species::add(CellInTissue& cell)
{
    // find a position for the new cell
    const size_t new_pos = num_of_cells();

    // set it position in the position map
    pos_map[cell.get_id()] = new_pos;

    cells.push_back(&cell);

    cell.genotype = get_id();

    return *this;
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
    out << "Species {genotype: " << static_cast<DriverGenotype>(species) 
        << ", cells: ";
    char sep{'{'};
    for (const auto& cell: species) {
        out << sep << " {id: " << cell.get_id() 
            << ", parent: " << cell.get_parent_id()
            << ", position: " << static_cast<PositionInTissue>(cell) 
            << "}";
        sep = ',';
    }
    out <<"}}";

    return out;
}

};