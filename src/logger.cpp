/**
 * @file logger.cpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Define simulation logger
 * @version 0.2
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

#include <iostream>

#include "logger.hpp"
#include "tissue.hpp"

namespace Races {

BasicLogger::BasicLogger()
{}

void BasicLogger::record(const CellEventType& type, const CellInTissue& cell, const Time& time)
{
    (void)type;
    (void)cell;
    (void)time;
}

void BasicLogger::snapshot(const Tissue& tissue, const Time& time)
{
    (void)tissue;
    (void)time;
}

JSONLogger::JSONLogger():
    BasicLogger(), os(std::cout)
{}

void JSONLogger::record(const CellEventType& type, const CellInTissue& cell, const Time& time)
{
    switch(type) {
        case CellEventType::DIE:
            os << "K " << cell.get_id();
            break;
        case CellEventType::DUPLICATE:
            os << "D " << cell.get_id() << " " << cell.get_parent_id() << " " << cell.get_driver_genotype();
            break;
        case CellEventType::DUPLICATION_AND_EPIGENETIC_EVENT:
            os << "F " << cell.get_id() << " " << cell.get_parent_id() << " " << cell.get_driver_genotype();
            break;
        case CellEventType::EPIGENETIC_EVENT:
            os << "E " << cell.get_id() << " " << cell.get_parent_id() << " " << cell.get_driver_genotype();
            break;
        default:
            throw std::domain_error("Unhandled event");
    }

    os << " " << time << std::endl;
}

/*
void JSONLogger::snapshot(const CellInTissue& cell)
{
    out << " {id: " << cell.get_id() 
            << ", parent: " << cell.get_parent_id()
            << ", position: " << static_cast<PositionInTissue>(cell) 
            << "}";
        sep = ',';
}

void JSONLogger::snapshot(const Species& species)
{
    out << "{name: " << species.get_name()
        << ", id: " << species.get_id() 
        << ", methylated: " << (species.is_methylated()?1:0)
        << ", cells: ";
    char sep{'{'};
    for (const auto& cell: species) {
        out << sep << " {id: " << cell.get_id() 
            << ", parent: " << cell.get_parent_id()
            << ", position: " << static_cast<Position>(cell) 
            << "}";
        sep = ',';
    }
    out <<"}}";

    os << "{time: "<< time <<", {name: \"" << tissue.get_name() << "\"";
    for (const auto& species: tissue) {
        os << "," << species;
    }

    os << "}";
}
*/

void JSONLogger::snapshot(const Tissue& tissue, const Time& time)
{
    os << "{time: "<< time <<", {name: \"" << tissue.get_name() << "\"";
    for (const auto& species: tissue) {
        os << "," << species;
    }

    os << "}" << std::endl;
}

}
