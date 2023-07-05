/**
 * @file driver.cpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Driver genotype representation
 * @version 0.2
 * @date 2023-07-05
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
#include <sstream>

#include "driver_genotype.hpp"

namespace Races {

unsigned int Races::DriverGenotype::counter = 0;

DriverGenotype::DriverGenotype(const std::string& name,
                               const std::map<CellEventType, double>& rates,
                               const bool methylated): id(counter++), name(name), 
                               rates(rates), methylated(methylated)
{
    std::map<CellEventType, std::string> event_types{
        {{CellEventType::DIE}, "death rate"},
        {{CellEventType::DUPLICATE}, "duplicate rate"},
        {{CellEventType::PASSENGER_MUTATION}, "passenger mutation rate"},
    };

    for (const auto& [event_type, name]: event_types) {
        if (rates.find(event_type)==rates.end()) {
            std::ostringstream oss;

            oss << "Missing " << name << " for driver \""
                << get_epigenetic_name() << "\".";

            throw std::domain_error(oss.str());
        }
    }
}

DriverGenotype::DriverGenotype(const DriverGenotype& orig): id(orig.id), name(orig.name), 
                               rates(orig.rates), methylated(orig.methylated)
{}

DriverGenotype& DriverGenotype::operator=(const DriverGenotype& orig)
{
    id = orig.id;
    name = orig.name;
    rates = orig.rates;
    methylated = orig.methylated;

    return *this;
}

DriverGenotype& DriverGenotype::operator=(DriverGenotype&& orig)
{
    id = orig.id;
    std::swap(name,orig.name);
    std::swap(rates,orig.rates);
    methylated = orig.methylated;

    return *this;
}

std::ostream& operator<<(std::ostream& out, const DriverGenotype& driver)
{
    out << driver.get_epigenetic_name() 
        << "(die rate: " << driver.get_rate(CellEventType::DIE)<< ", "
        << "duplication rate: " << driver.get_rate(CellEventType::DUPLICATE)
        << ")";
    return out; 
}

};