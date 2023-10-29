/**
 * @file cell_event.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements cell events
 * @version 0.6
 * @date 2023-10-29
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

#include <string>
#include <map>

#include "cell_event.hpp"

namespace Races 
{

namespace Drivers 
{

std::map<CellEventType, std::string> cell_event_names = {
    {CellEventType::DEATH, "death"},
    {CellEventType::DUPLICATION, "duplication"},
    {CellEventType::EPIGENETIC_EVENT, "epigenetic mutation"},
    {CellEventType::DUPLICATION_AND_EPIGENETIC_EVENT, 
        "duplication and epigenetic mutation"},
    {CellEventType::DRIVER_GENETIC_MUTATION, 
        "driver genomic mutation"}
};

namespace Simulation
{

CellEvent::CellEvent():
    type(), position(), initial_genotype(), final_genotype(), delay()
{}

}

}   // Drivers

}   // Races
