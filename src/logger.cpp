/**
 * @file logger.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements simulation loggers
 * @version 0.11
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

#include <iostream>
#include <chrono>

#include "logger.hpp"
#include "tissue.hpp"

namespace Races 
{

namespace Drivers 
{

namespace Simulation 
{

std::string BasicLogger::get_time_string()
{
    std::time_t time;
    std::tm* info;
    char buffer[81];

    std::time(&time);
    info = std::localtime(&time);

    std::strftime(buffer,80,"%Y%m%d-%H%M%S",info);

    return buffer;
}

BasicLogger::BasicLogger():
    BasicLogger("races_"+get_time_string())
{}

BasicLogger::BasicLogger(const std::filesystem::path simulation_dir):
    directory(simulation_dir)
{
}

void BasicLogger::record(const CellEventType& type, const CellInTissue& cell, const Time& time)
{
    (void)type;
    (void)cell;
    (void)time;
}

void BasicLogger::record_initial_cell(const CellInTissue& cell)
{
    if (cell.get_id() != cell.get_parent_id()) {
        throw std::domain_error("The provided cell is not an initial cell");
    }

    (void)cell;
}

void BasicLogger::snapshot(const Simulation& simulation)
{
    (void)simulation;
}


void
BasicLogger::save_sampled_ids(const std::filesystem::path simulation_dir,
                              const std::list<Races::Drivers::CellId>& sampled_cell_ids,
                              const Races::Time& time,
                              const Races::Drivers::RectangleSet& sampled_region)
{

    std::ofstream os(simulation_dir/"sampled_ids.list", std::ofstream::app);

    os << "# " << time
       << " " << sampled_cell_ids.size()
       << " " << sampled_region.size()
       << " "<< sampled_region.lower_corner
       << " "<< sampled_region.upper_corner << std::endl;

    for (const auto& cell_id: sampled_cell_ids) {
        os << cell_id << std::endl;
    }
}

std::list<Races::Drivers::CellId>
BasicLogger::load_sampled_ids(const std::filesystem::path simulation_dir)
{
    using namespace Races::Drivers;

    std::list<CellId> sample;

    std::ifstream sample_is(simulation_dir/"sampled_ids.list");
    while (!sample_is.eof()) {
        std::string line;

        std::getline(sample_is, line);
        if (line[0] != '#' && line[0] != '\n') {
            std::istringstream iss(line);

            CellId cell_id;
            iss >> cell_id;

            sample.push_back(cell_id);
        }
    }

    return sample;
}

}   // Simulation

}   // Drivers

}   // Races
