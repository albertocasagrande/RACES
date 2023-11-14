/**
 * @file tissue_sample.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements tissue samples
 * @version 0.5
 * @date 2023-11-14
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

#include "tissue_sample.hpp"

namespace Races
{

namespace Drivers
{

namespace Simulation
{

TissueSampleId TissueSample::counter = 0;

TissueSample::TissueSample():
    id(0), time(0), region({{0,0},{0,0}}), name("")
{}

TissueSample::TissueSample(const Time& time, const RectangleSet& region):
    TissueSample(time, region, {})
{}

TissueSample::TissueSample(const std::string& name, const Races::Time& time,
                           const Races::Drivers::RectangleSet& region):
    TissueSample(name, time, region, {})
{}

TissueSample::TissueSample(const Time& time, const RectangleSet& region,
                           const std::list<Races::Drivers::CellId>& cell_ids):
    TissueSample("", time, region, cell_ids)
{
    name = "S_"+std::to_string(id);
}

TissueSample::TissueSample(const std::string& name, const Time& time, const RectangleSet& region,
                           const std::list<Races::Drivers::CellId>& cell_ids):
    id(counter++), time(time), region(region), cell_ids(cell_ids), name(name)
{}

void TissueSample::add_cell_id(const Races::Drivers::CellId& cell_id)
{
    if (region.size() == cell_ids.size()) {
        throw std::domain_error("The sample already contains all the cell ids");
    }

    cell_ids.push_back(cell_id);
}

}   // Simulation

}   // Drivers

}   // Races

std::ostream& operator<<(std::ostream& os, const Races::Drivers::Simulation::TissueSample& tissue_sample)
{
    const Races::Drivers::RectangleSet& region =  tissue_sample.get_region();

    os << "# " << tissue_sample.get_time()
       << " " << tissue_sample.get_cell_ids().size()
       << " " << region.size()
       << " "<< region.lower_corner
       << " "<< region.upper_corner << std::endl;

    for (const auto& cell_id: tissue_sample.get_cell_ids()) {
        os << cell_id << std::endl;
    }

    return os;
}
