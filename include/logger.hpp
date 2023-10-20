/**
 * @file logger.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines simulation loggers
 * @version 0.9
 * @date 2023-10-20
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

#include <string>
#include <filesystem>

#include "cell.hpp"
#include "cell_event.hpp"
#include "time.hpp"
#include "position_set.hpp"

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
class BasicLogger
{
protected:
    std::filesystem::path directory;   //!< the log directory

    /**
     * @brief Get a string representing the current time
     * 
     * This method produces a string in the format
     *    {year}{month}{day}_{hour}{minutes}{seconds}
     * 
     * @return a string representing the current time 
     */
    static std::string get_time_string();
public:
    /**
     * @brief The empty constructor
     */
    BasicLogger();

    /**
     * @brief A constructor
     * 
     * @param simulation_dir is the new simulation logging directory
     */
    explicit BasicLogger(const std::filesystem::path simulation_dir);

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
     * @brief Save sampled cell ids
     * 
     * @param simulation_dir is the path of the simulation directory
     * @param sampled_cell_ids is the list of sampled cell ids
     * @param time it the sampling time
     * @param sampled_region is the sampled region
     * @return the updated output stream
     */
    static void save_sampled_ids(const std::filesystem::path simulation_dir,
                                 const std::list<Races::Drivers::CellId>& sampled_cell_ids,
                                 const Races::Time& time,
                                 const Races::Drivers::RectangleSet& sampled_region);

    /**
     * @brief Load the sampled cell ids
     * 
     * @param simulation_dir is the path of the simulation directory
     * @return the list of the sampled cell identifiers
     */
    static std::list<Races::Drivers::CellId>
    load_sampled_ids(const std::filesystem::path simulation_dir);

    /**
     * @brief Close open archives
     */
    inline void close()
    {}

    /**
     * @brief Reset the logger
     * 
     * @param simulation_dir is the new simulation logging directory
     */
    inline void reset(const std::filesystem::path simulation_dir)
    {
        close();

        directory = simulation_dir;
    }

    /**
     * @brief Get the log directory
     * 
     * @return the directory containing all the logs
     */
    inline const std::filesystem::path& get_directory() const
    {
        return directory;
    }
    /**
     * @brief Save sampled cell ids
     * 
     * @param sampled_cell_ids is the list of sampled cell ids
     * @param time it the sampling time
     * @param sampled_region is the sampled region
     * @return the updated output stream
     */
    inline void save_sampled_ids(const std::list<Races::Drivers::CellId>& sampled_cell_ids,
                                 const Races::Time& time,
                                 const Races::Drivers::RectangleSet& sampled_region) const
    {
        BasicLogger::save_sampled_ids(directory, sampled_cell_ids, time, sampled_region);
    }

    /**
     * @brief Load the sampled cell ids
     * 
     * @return the list of the sampled cell identifiers
     */
    inline std::list<Races::Drivers::CellId> load_sampled_ids() const
    {
        return BasicLogger::load_sampled_ids(directory);
    }
};

}   // Simulation

}   // Drivers

}   // Races

#endif // __RACES_LOGGER__