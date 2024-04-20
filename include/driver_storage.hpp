/**
 * @file driver_storage.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines class to load and store driver mutations
 * @version 0.2
 * @date 2024-04-20
 *
 * @copyright Copyright (c) 2023-2024
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

#ifndef __RACES_DRIVER_STORAGE__
#define __RACES_DRIVER_STORAGE__

#include <map>
#include <list>
#include <string>
#include <filesystem>

#include "snv.hpp"

namespace Races
{

namespace Mutations
{

class DriverStorage
{
    std::map<std::string, SNV> SNVs;    //!< The code-SNV map

    std::filesystem::path source_path;  //!< The driver filename
public:
    /**
     * @brief The empty constructor
     */
    DriverStorage();

    /**
     * @brief Get the driver SNVs
     *
     * @return a constant reference to a map associating a driver mutation
     *      code to the correspoding driver SNV
     */
    inline const std::map<std::string, SNV>& get_SNVs() const
    {
        return SNVs;
    }

    /**
     * @brief Get the SNV positions
     *
     * @return A list containing the genomic position of the SNVs
     */
    std::list<GenomicPosition> get_SNV_positions() const;

    /**
     * @brief Get the source path of the drivers
     *
     * @return a constant reference to the driver source path
     */
    const std::filesystem::path& get_source_path() const
    {
        return source_path;
    }

    /**
     * @brief Load the driver SNVs
     *
     * @param filename is the driver mutation filename
     * @return a map associating a driver code to the corresponding driver SNV
     */
    static DriverStorage load(const std::filesystem::path& filename);
};

}   // Mutations

}   // Races

#endif // __RACES_DRIVER_MUTATION_STORAGE__