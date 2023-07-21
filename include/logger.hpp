/**
 * @file logger.hpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Define simulation logger
 * @version 0.6
 * @date 2023-07-21
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

#ifndef __RACES_LOGGER__
#define __RACES_LOGGER__

#include "cell.hpp"
#include "cell_event.hpp"
#include "time.hpp"

namespace Races 
{

namespace Drivers 
{

namespace Simulation 
{

class Simulation;

/**
 * @brief The simulator logger concept
 */
struct BasicLogger
{
    /**
     * @brief The empty constructor
     */
    BasicLogger();

    /**
     * @brief Record an event
     * 
     * @param type is the event type
     * @param cell is the cell on which event has been occurred
     * @param time it the event time
     */
    void record(const CellEventType& type, const CellInTissue& cell, const Time& time);

    /**
     * @brief Record an initial cell
     * 
     * @param cell is the initial cell to record 
     */
    void record_initial_cell(const CellInTissue& cell);

    /**
     * @brief Save a simulation snapshot
     * 
     * @param simulation is the simulation whose snapshot is requested
     */
    void snapshot(const Simulation& simulation);

    /**
     * @brief Flush archive data
     */
    inline void flush_archives()
    {}

    /**
     * @brief Close open archives
     */
    inline void close()
    {}

    /**
     * @brief Reset the logger
     */
    inline void reset()
    {}
};

struct JSONLogger : public BasicLogger
{
    std::ostream& os;  //!< The output stream

    /**
     * @brief The empty constructor
     */
    JSONLogger();

    /**
     * @brief Record an event
     * 
     * @param type is the event type
     * @param cell is the cell on which event has been occurred
     * @param time it the event time
     */
    void record(const CellEventType& type, const CellInTissue& cell, const Time& time);

    /**
     * @brief Record an initial cell
     * 
     * @param cell is the initial cell to record 
     */
    void record_initial_cell(const CellInTissue& cell);

    /**
     * @brief Save a simulation snapshot
     * 
     * @param simulation is the simulation whose snapshot is requested
     */
    void snapshot(const Simulation& simulation);

    /**
     * @brief Flush archive data
     */
    inline void flush_archives()
    {
        os.flush();
    }
};

}   // Simulation

}   // Drivers

}   // Races

#endif // __RACES_LOGGER__