/**
 * @file sqlite_logger.hpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Define a binary simulation logger
 * @version 0.1
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

#ifndef __RACES_BINARY_LOGGER__
#define __RACES_BINARY_LOGGER__

#include <fstream>
#include <filesystem>

#include "logger.hpp"

namespace Races {

/**
 * @brief A binary logger
 * 
 * The objects of this class record cell duplications and 
 * tissue status in binary files 
 */
struct BinaryLogger : public BasicLogger
{
    std::filesystem::path directory;   //!< the log directory 
    std::ofstream cell_of;             //!< the current cell output file
    const size_t cells_per_file;       //!< the number of cells to be saved before changing file
    size_t cell_in_current_file;       //!< cells saved in current file
    uint16_t file_number;                //!< the current file number

    std::filesystem::path get_next_cell_path() const;

    std::filesystem::path get_next_snapshot_path(const Time& time) const;
public:
    /**
     * @brief The empty constructor
     */
    BinaryLogger();

    /**
     * @brief A constructor
     * 
     * @param name is the prefix of the filename
     */
    BinaryLogger(const std::string prefix_name);

    /**
     * @brief A constructor
     * 
     * @param prefix_name is the prefix of the filename
     * @param cells_per_file is the number of cells per file
     */
    BinaryLogger(const std::string prefix_name, const size_t cells_per_file);
    
    /**
     * @brief Record an event
     * 
     * @param type is the event type
     * @param cell is the cell on which event has been occurred
     * @param time it the event time
     */
    void record(const CellEventType& type, const CellInTissue& cell, const Time& time);

    /**
     * @brief Save a tissue snapshot
     * 
     * @param tissue is the tissue whose snapshot is requested
     * @param time is the snapshot time
     */
    void snapshot(const Tissue& tissue, const Time& time);

    /**
     * @brief The destructor
     */
    ~BinaryLogger();
};

}

#endif // __RACES_BINARY_LOGGER__