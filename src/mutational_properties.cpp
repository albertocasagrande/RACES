/**
 * @file mutational_properties.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements a class to represent the mutational properties
 * @version 1.1
 * @date 2024-08-04
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

#include "mutational_properties.hpp"

namespace RACES
{

namespace Mutations
{

PassengerRates::PassengerRates():
    indel{0}, snv{0}, cna{0}
{}

PassengerRates::PassengerRates(const double& indel_rate, const double& SNV_rate,
                               const double& CNA_rate):
    indel{indel_rate}, snv{SNV_rate}, cna{CNA_rate}
{}

DriverMutations::DriverMutations()
{}

//! @private
template<typename MUTATION>
std::set<MUTATION> build_mutation_set(const std::string& mutant_name,
                                      const std::list<MUTATION>& mutations)
{
    std::set<MUTATION> mutation_set;

    for (const auto& mutation : mutations) {
        if (mutation_set.count(mutation)>0) {
            std::ostringstream oss;

            oss << mutation << " added twice among "
                << mutant_name << "'s driver mutations.";

            throw std::domain_error(oss.str());
        }

        mutation_set.insert(mutation);
    }

    return mutation_set;
}

std::list<DriverMutations::MutationType>
DriverMutations::get_default_order(const std::list<MutationSpec<SID>>& SIDs,
                                   const std::list<CNA>& CNAs,
                                   const bool& wg_doubling)
{
    std::list<DriverMutations::MutationType> application_order;

    for (size_t i=0; i<SIDs.size(); ++i) {
        application_order.push_back(DriverMutations::MutationType::SID_TURN);
    }

    for (size_t i=0; i<CNAs.size(); ++i) {
        application_order.push_back(DriverMutations::MutationType::CNA_TURN);
    }

    if (wg_doubling) {
        application_order.push_back(DriverMutations::MutationType::WGD_TURN);
    }

    return application_order;
}

DriverMutations::DriverMutations(const std::string& mutant_name,
                                 const std::list<MutationSpec<SID>>& SIDs,
                                 const std::list<CNA>& CNAs,
                                 const bool& wg_doubling):
    name{mutant_name}, SIDs{SIDs}, CNAs{CNAs},
    application_order{DriverMutations::get_default_order(SIDs,CNAs,wg_doubling)}
{}

DriverMutations::DriverMutations(const std::string& mutant_name,
                                 const std::list<MutationSpec<SID>>& SIDs,
                                 const std::list<CNA>& CNAs,
                                 const std::list<MutationType>& application_order):
    name{mutant_name}, SIDs{SIDs}, CNAs{CNAs},
    application_order{application_order}
{
    size_t SID_count{0}, CNA_count{0};
    for (auto order_it = application_order.begin();
         order_it != application_order.end(); ++order_it) {

        switch(*order_it) {
            case SID_TURN:
                ++SID_count;
                break;
            case CNA_TURN:
                ++CNA_count;
                break;
            case WGD_TURN:
                break;
            default:
                throw std::domain_error("Unsupported driver mutation type.");
        }
    }

    if (SID_count != SIDs.size()) {
        throw std::domain_error("The number of SNVs/indels differs from "
                                "that of the same kind of mutations "
                                "in the application order list.");
    }

    if (CNA_count != CNAs.size()) {
        throw std::domain_error("The number of CNAs differs from "
                                "that of the same kind of mutations "
                                "in the application order list.");
    }
}

MutationalProperties::MutationalProperties()
{}

MutationalProperties&
MutationalProperties::add_mutant(const std::string& mutant_name,
                                 const std::map<std::string, PassengerRates>& epistate_passenger_rates,
                                 const std::list<MutationSpec<SID>>& driver_SIDs,
                                 const std::list<CNA>& driver_CNAs,
                                 const bool& wg_doubling)
{
    auto application_order = DriverMutations::get_default_order(driver_SIDs, driver_CNAs, wg_doubling);

    return add_mutant(mutant_name, epistate_passenger_rates, driver_SIDs,
                      driver_CNAs, application_order);
}

MutationalProperties&
MutationalProperties::add_mutant(const std::string& mutant_name,
                                 const std::map<std::string, PassengerRates>& epistate_passenger_rates,
                                 const std::list<MutationSpec<SID>>& driver_SIDs,
                                 const std::list<CNA>& driver_CNAs,
                                 const std::list<DriverMutations::MutationType>& application_order)
{
    if (driver_mutations.count(mutant_name)>0) {
        throw std::domain_error("The mutational properties of mutant \""+mutant_name+
                                "\" has been already added.");
    }

    driver_mutations.insert({mutant_name, {mutant_name, driver_SIDs,
                                           driver_CNAs, application_order}});

    for (const auto& [epistate, passenger_rate] : epistate_passenger_rates) {
        auto species_name = mutant_name+epistate;

        passenger_rates[species_name] = passenger_rate;
    }

    return *this;
}

}   // Mutations

}   // RACES
