/**
 * @file driver_storage.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements class to load and store driver mutations
 * @version 1.1
 * @date 2025-07-10
 *
 * @copyright Copyright (c) 2023-2025
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

#include "driver_storage.hpp"

#include "csv_reader.hpp"

namespace RACES
{

namespace Mutations
{

DriverStorage::DriverStorage()
{}

std::map<SID, std::string> DriverStorage::get_reverse_map() const
{
    std::map<SID, std::string> reverse_map;

    for (const auto& [code, snv] : mutations) {
        reverse_map[snv] = code;
    }

    return reverse_map;
}

std::list<GenomicPosition> DriverStorage::get_mutation_positions() const
{
    std::list<GenomicPosition> genomic_positions;

    for (const auto& [code, mutation] : mutations) {
        genomic_positions.push_back(mutation);
    }

    return genomic_positions;
}

DriverStorage DriverStorage::load(const std::filesystem::path& filename)
{
    DriverStorage mutations;

    RACES::IO::CSVReader csv_reader(filename, true, '\t');

    size_t chr_col = csv_reader.get_column_position("chr");
    size_t pos_col = csv_reader.get_column_position("from");
    size_t alt_col = csv_reader.get_column_position("alt");
    size_t ref_col = csv_reader.get_column_position("ref");
    size_t type_col = csv_reader.get_column_position("mutation_type");
    size_t code_col = csv_reader.get_column_position("driver_code");

    for (const auto& row : csv_reader) {
        const auto type = row.get_field(type_col);
        if (type=="SNV" || type=="indel") {
            auto chr_str = row.get_field(chr_col);

            auto chr_id = GenomicPosition::stochr(chr_str);
            auto pos = static_cast<ChrPosition>(stoul(row.get_field(pos_col)));

            mutations.mutations.insert({row.get_field(code_col),
                                        {chr_id, pos, row.get_field(ref_col),
                                         row.get_field(alt_col),
                                         Mutations::Mutation::DRIVER}});
        }
    }

    mutations.source_path = filename;

    return mutations;
}

}   // Mutations

}   // RACES
