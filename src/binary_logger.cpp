/**
 * @file binary_logger.cpp
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

#include <sstream>
#include <ctime>
#include <limits>
#include <cmath> // log10 and ceil

#include "binary_logger.hpp"
#include "tissue.hpp"

namespace Races {

BinaryLogger::BinaryLogger():
    BinaryLogger("races")
{
}

BinaryLogger::BinaryLogger(const std::string prefix_name):
    BinaryLogger(prefix_name, std::numeric_limits<size_t>::max())
{
}

std::filesystem::path BinaryLogger::get_next_cell_path() const
{
    using namespace std;

    std::ostringstream oss;

    size_t digits = static_cast<size_t>(ceil(log10(numeric_limits<uint16_t>::max())));

    oss << "cells_" << std::setfill('0') << std::setw(digits) << file_number;

    return directory / (oss.str()+".dat") ;
}

std::filesystem::path BinaryLogger::get_next_snapshot_path(const Time& time) const
{
    std::ostringstream oss;

    oss << "snapshot_" << time << ".dat";

    return directory / oss.str();
}

std::string get_directory_name(const std::string& prefix_name)
{
    std::time_t time;
    std::tm* info;
    char buffer [80];

    std::time(&time);
    info = std::localtime(&time);

    std::strftime(buffer,80,"%Y%m%d-%H%M%S",info);

    std::ostringstream oss;

    oss << prefix_name << "_" << buffer;

    return oss.str();
}

BinaryLogger::BinaryLogger(const std::string prefix_name, const size_t cells_per_file):
    BasicLogger(), directory(get_directory_name(prefix_name)), cell_of(), 
    cells_per_file(cells_per_file), cell_in_current_file(0), file_number(0)
{
    std::filesystem::create_directory(directory);

    auto filename = get_next_cell_path();

    cell_of.open(filename, std::ios::out | std::ios::binary);
}

BinaryLogger::~BinaryLogger()
{
    cell_of.close();
}

void BinaryLogger::record(const CellEventType& type, const CellInTissue& cell, const Time& time)
{
    if (type==CellEventType::DUPLICATE || type==CellEventType::DUPLICATION_AND_EPIGENETIC_EVENT) {

        cell_of.write((char*)(&(cell.get_parent_id())), sizeof(CellId));
        cell_of.write((char*)(&(cell.get_driver_genotype())), sizeof(DriverGenotypeId));
        cell_of.write((char*)(&time), sizeof(Time));

        if (++cell_in_current_file>=cells_per_file) {
            cell_of.close();

            cell_in_current_file=0;
            ++file_number;

            auto filename = get_next_cell_path();

            cell_of.open(filename, std::ios::out | std::ios::binary);
        }
    }
}

void BinaryLogger::snapshot(const Tissue& tissue, const Time& time)
{
    auto filename = get_next_snapshot_path(time);

    std::ofstream ofs(filename, std::ios::out | std::ios::binary);

    for (const auto& species:tissue) {
        for (const auto& cell:species) {
            ofs.write((char*)(&(cell.get_id())), sizeof(CellId));
            ofs.write((char*)(&(time)), sizeof(Time));
            ofs.write((char*)(&(cell.x)), sizeof(int16_t));
            ofs.write((char*)(&(cell.y)), sizeof(int16_t));
            ofs.write((char*)(&(cell.z)), sizeof(int16_t));
        }
    }

    ofs.close();
}

}