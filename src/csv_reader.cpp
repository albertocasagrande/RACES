/**
 * @file csv_reader.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements a class to read CSV
 * @version 0.3
 * @date 2024-02-14
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

#include <sstream>

#include "csv_reader.hpp"

namespace Races
{

namespace IO
{

CSVReader::CSVRow::CSVRow():
    fields()
{}

CSVReader::CSVRow::CSVRow(const std::string& row, const char& column_separator)
{
    std::string field;
    std::stringstream oss_line(row);

    while(std::getline(oss_line, field, column_separator)) {
        fields.push_back(field);
    }
}

const std::string& CSVReader::CSVRow::get_field(const size_t& index) const
{
    if (index<fields.size()) {
        return fields[index];
    }

    throw std::out_of_range("The CSV row contains " 
                            + std::to_string(fields.size())
                            + " fields. Requested field in position "
                            + std::to_string(index));
}

bool CSVReader::const_iterator::get_row(CSVReader::CSVRow& row)
{
    std::string line;

    if (!getline(ifs, line)) {
        return false;
    }

    row = CSVRow(line, reader->col_sep);

    return true;
}

CSVReader::const_iterator::const_iterator(const CSVReader* reader, const std::streamoff& off, 
                                          const std::ios_base::seekdir& way):
    reader(reader), ifs(reader->filename)
{
    ifs.seekg(off, way);

    eof = !get_row(curr_row);
}

CSVReader::const_iterator& CSVReader::const_iterator::operator++()
{
    eof = !get_row(curr_row);

    return *this;
}

CSVReader::CSVReader(const std::filesystem::path& filename, const bool has_header,
                     const char column_separator):
    filename(filename), col_sep(column_separator), first_row_pos(0)
{
    if (has_header) {
        std::ifstream ifs(filename);
        std::string line;
        
        if (!getline(ifs, line)) {
            throw std::runtime_error("The file does not have a header row.");
        }

        std::string column_name;
        std::stringstream oss_line(line);

        size_t i{0};
        while(std::getline(oss_line, column_name, col_sep)) {
            header.push_back(column_name);
            columns[column_name] = i;

            ++i;
        }

        first_row_pos = ifs.tellg();
    }
}

}   // IO

}   // Races