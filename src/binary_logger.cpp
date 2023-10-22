/**
 * @file binary_logger.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements a binary simulation logger
 * @version 0.19
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

#include <sstream>
#include <ctime>
#include <limits>
#include <cmath> // log10 and ceil

#include "binary_logger.hpp"
#include "tissue.hpp"

#include "simulation.hpp"

namespace Races 
{

namespace Drivers 
{

namespace Simulation 
{

BinaryLogger::BinaryLogger():
    BinaryLogger("races_"+BasicLogger::get_time_string())
{
}

BinaryLogger::BinaryLogger(const std::filesystem::path simulation_dir, const size_t cells_per_file):
    BasicLogger(simulation_dir), cell_archive(), 
    cells_per_file(cells_per_file), cell_in_current_file(0), current_file_number(0)
{
}

BinaryLogger::BinaryLogger(const BinaryLogger& orig):
    BasicLogger(orig),  cell_archive(),
    cells_per_file(orig.cells_per_file), cell_in_current_file(orig.cell_in_current_file), 
    current_file_number(orig.current_file_number)
{
}

std::filesystem::path BinaryLogger::get_cell_archive_path(const std::filesystem::path& directory, const uint16_t& file_number)
{
    using namespace std;

    std::ostringstream oss;

    size_t digits = static_cast<size_t>(ceil(log10(numeric_limits<uint16_t>::max())));

    oss << "cells_" << std::setfill('0') << std::setw(digits) << file_number;

    return directory / (oss.str()+".dat") ;
}

std::filesystem::path BinaryLogger::get_snapshot_path() const
{
    std::ostringstream oss;

    oss << "snapshot_" << BasicLogger::get_time_string() << ".dat";

    return directory / oss.str();
}

BinaryLogger& BinaryLogger::operator=(const BinaryLogger& orig)
{
    close();

    directory = orig.directory;
    cells_per_file = orig.cells_per_file;
    cell_in_current_file = orig.cell_in_current_file;
    current_file_number = orig.current_file_number;

    return *this;
}

std::streampos compute_bytes_per_cell(const std::filesystem::path& directory)
{ 
    const auto first_archive_name = BinaryLogger::get_cell_archive_path(directory, 0);
    Races::Archive::Binary::In cell_archive(first_archive_name);

    const auto cell = Cell::load(cell_archive); 
    return Archive::Binary::ByteCounter::bytes_required_by(cell);
}

BinaryLogger::CellReader::CellReader(const std::filesystem::path& directory):
    directory(directory), bytes_per_cell(1), number_of_cells(0), cell_archive_id(0), 
    cell_archive(get_cell_archive_path(directory, 0))
{
    std::ifstream cell_file(directory / "cell_file_info.txt", std::fstream::in);
    cell_file >> cells_per_archive >> last_file_number;

    // if no cell was saved 
    if (cell_archive.size()==0) {
        return;
    }

    bytes_per_cell = compute_bytes_per_cell(directory);

    {   // compute the total number of cells
        const auto filename = get_cell_archive_path(directory, last_file_number);

        const auto lastfile_bytes = Archive::Binary::In(filename).size();

        number_of_cells = static_cast<uint64_t>(lastfile_bytes/bytes_per_cell)
                            + cells_per_archive * last_file_number;
    }
}

Cell BinaryLogger::CellReader::operator[](const CellId& cell_id) 
{
    if (cell_id >= number_of_cells) {
        throw std::out_of_range("The cell has not been created");
    }

    uint16_t aimed_id = static_cast<uint16_t>(cell_id / cells_per_archive);

    if (!cell_archive.is_open() || cell_archive_id!=aimed_id) {
        if (cell_archive_id!=aimed_id) {
            cell_archive.close();
        }
        
        auto filename = get_cell_archive_path(directory, last_file_number);

        cell_archive.open(filename);
        cell_archive_id = aimed_id;
    }

    cell_archive.seekg(get_cell_pos(cell_id));

    auto cell = Cell::load(cell_archive);

    if (cell.get_id()!=cell_id) {
        std::cerr << "Aiming cell id " << cell_id << " read " 
                  << cell << std::endl;
        throw std::runtime_error("Wrong cell file format.");
    }

    return cell;
}

void BinaryLogger::rotate_cell_file()
{
    if (cell_archive.is_open()) {
        cell_archive.close();
    }

    {
        // save the number of cells per file and the next file number to be used
        std::ofstream cell_file(directory / "cell_file_info.txt", std::fstream::out);
        cell_file << cells_per_file << " " << ++current_file_number << std::endl; 
    }

    auto filename = get_cell_archive_path(directory, current_file_number);

    cell_in_current_file=0;
    cell_archive.open(filename, std::fstream::binary);
}

void BinaryLogger::record_cell(const CellInTissue& cell)
{
    if (!cell_archive.is_open()) {
        if (std::filesystem::exists(directory)) {
            auto filename = get_cell_archive_path(directory, current_file_number);

            cell_archive.open(filename, std::fstream::binary | std::fstream::app);
        } else {
            std::filesystem::create_directory(directory);

            rotate_cell_file();
        }
    }

    cell_archive & static_cast<const Cell&>(cell);

    if (++cell_in_current_file>=cells_per_file) {
        rotate_cell_file();
    } 
}

void BinaryLogger::record(const CellEventType& type, const CellInTissue& cell, const Time& time)
{
    (void)time;

    if (type==CellEventType::DUPLICATE || 
            type==CellEventType::DUPLICATION_AND_EPIGENETIC_EVENT ||
            type==CellEventType::DRIVER_GENETIC_MUTATION) {

        record_cell(cell);
    }
}

void BinaryLogger::record_initial_cell(const CellInTissue& cell) 
{
    if (cell.get_id() != cell.get_parent_id()) {
        throw std::domain_error("The provided cell is not an initial cell");
    }

    record_cell(cell);
}


void BinaryLogger::snapshot(const Simulation& simulation)
{
    Archive::Binary::Out archive(get_snapshot_path());

    archive & simulation;
}

void BinaryLogger::reset(const std::string& output_directory)
{
    close();

    directory = output_directory;
    cell_in_current_file = 0;
    current_file_number = 0;
}

void
BinaryLogger::save_sample(const std::filesystem::path simulation_dir,
                          const Races::Drivers::Simulation::TissueSample& tissue_sample)
{

    Races::Archive::Binary::Out archive(simulation_dir/"samples.dat", std::ios_base::app);

    archive & tissue_sample;
}


std::list<Races::Drivers::Simulation::TissueSample>
BinaryLogger::load_samples(const std::filesystem::path simulation_dir)
{

    Races::Archive::Binary::In archive(simulation_dir/"samples.dat");

    std::list<Races::Drivers::Simulation::TissueSample> samples;

    while (!archive.eof()) {
        samples.push_back(Races::Drivers::Simulation::TissueSample::load(archive));
    }

    return samples;
}

BinaryLogger::~BinaryLogger()
{
    close();
}

}   // Simulation

}   // Drivers

}   // Races
