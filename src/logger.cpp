/**
 * @file logger.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements simulation loggers
 * @version 0.12
 * @date 2023-10-23
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
BasicLogger::save_sample(const std::filesystem::path simulation_dir,
                         const Races::Drivers::Simulation::TissueSample& tissue_sample)
{

    std::ofstream os(simulation_dir/"samples.list", std::ofstream::app);

    os << tissue_sample;
}

}   // Simulation

}   // Drivers

}   // Races
