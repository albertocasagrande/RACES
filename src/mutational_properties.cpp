/**
 * @file mutational_properties.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements a class to represent the mutational properties
 * @version 0.17
 * @date 2024-04-19
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

#include "mutational_properties.hpp"

namespace Races
{

namespace Mutations
{

PassengerRates::PassengerRates():
    snv{0}, cna{0}
{}

PassengerRates::PassengerRates(const double& SNV_rate, const double& CNA_rate):
    snv{SNV_rate}, cna{CNA_rate}
{}

DriverMutations::DriverMutations()
{}

DriverMutations::DriverMutations(const std::string& mutant_name,
                                 const std::list<MutationSpec<SNV>>& SNVs,
                                 const std::list<CNA>& CNAs):
    name(mutant_name), SNVs(SNVs.begin(), SNVs.end()), CNAs(CNAs.begin(), CNAs.end())
{}

MutationalProperties::MutationalProperties()
{}

MutationalProperties&
MutationalProperties::add_mutant(const std::string& mutant_name,
                                 const std::map<std::string, PassengerRates>& epistate_passenger_rates,
                                 const std::list<MutationSpec<SNV>>& driver_SNVs,
                                 const std::list<CNA>& driver_CNAs)
{
    if (driver_mutations.count(mutant_name)>0) {
        throw std::domain_error("The mutational properties of mutant \""+mutant_name+
                                "\" has been already added.");
    }

    driver_mutations.insert({mutant_name, {mutant_name, driver_SNVs, driver_CNAs}});


    for (const auto& [epistate, passenger_rate] : epistate_passenger_rates) {
        auto species_name = mutant_name+epistate;

        passenger_rates[species_name] = passenger_rate;
    }

    return *this;
}

}   // Mutations

}   // Races
