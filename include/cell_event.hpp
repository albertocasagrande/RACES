/**
 * @file cell_event.hpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Cell events
 * @version 0.4
 * @date 2023-07-08
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

#ifndef __RACES_CELL_EVENT__
#define __RACES_CELL_EVENT__

#include <string>
#include <map>

#include "driver_genotype_id.hpp"
#include "position.hpp"
#include "time.hpp"

namespace Races {

/**
 * @brief Supported event type
 * 
 */
enum class CellEventType {
    DIE,
    DUPLICATE,
    EPIGENETIC_EVENT,
    DUPLICATION_AND_EPIGENETIC_EVENT,
    PASSENGER_MUTATION,
    DRIVER_SOMATIC_MUTATION
};

/**
 * @brief A map associating a name to each cell event
 */
extern std::map<CellEventType, std::string> cell_event_names;

/**
 * @brief A description of a simulation event
 */
struct CellEvent
{
    CellEventType type;                     //!< event type
    Position position;                      //!< event position
    EpigeneticGenotypeId initial_genotype;  //!< genotype of the cell on which event occurs
    EpigeneticGenotypeId final_genotype;    //!< final genotype for mutational events
    Time delay;                             //!< event delay
};

};

#endif // __RACES_CELL_EVENT__