/**
 * @file mutational_properties.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines a class to represent the mutational properties
 * @version 1.2
 * @date 2024-08-09
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

#ifndef __RACES_MUTATIONAL_PROPERTIES__
#define __RACES_MUTATIONAL_PROPERTIES__

#include <map>
#include <list>
#include <string>

#include "sid.hpp"
#include "cna.hpp"
#include "mutation_spec.hpp"

#include "mutation_list.hpp"
namespace RACES
{

namespace Mutations
{

/**
 * @brief Rates of the passenger mutations for a species
 */
struct PassengerRates
{
    double indel;   //!< The species indel rate
    double snv;     //!< The species SNV rate
    double cna;     //!< The species CNA rate

    /**
     * @brief The empty constructor
     */
    PassengerRates();

    /**
     * @brief A constructor
     *
     * @param indel_rate is the the SNV rate
     * @param SNV_rate is the the SNV rate
     * @param CNA_rate is the the CNA rate
     */
    PassengerRates(const double& indel_rate, const double& SNV_rate,
                   const double& CNA_rate);
};

/**
 * @brief A class representing the mutational properties of all the species
 *
 */
class MutationalProperties
{
    std::map<std::string, PassengerRates> passenger_rates;      //!< The species passenger rates
    std::map<std::string, DriverMutations> driver_mutations;    //!< The mutation per mutant

public:

    /**
     * @brief The empty constructor
     */
    MutationalProperties();

    /**
     * @brief Add the properties of a mutant
     *
     * @param mutant_name is the name of the mutant
     * @param epistate_passenger_rates is a map from epigenomic state to
     *          passenger rates
     * @param driver_SIDs is a list of driver SIDs
     * @param driver_CNAs is a list of driver CNAs
     * @param wg_doubling is a Boolean flag to enable whole genome doubling
     * @return a reference to the updated object
     */
    MutationalProperties& add_mutant(const std::string& mutant_name,
                                     const std::map<std::string, PassengerRates>& epistate_passenger_rates,
                                     const std::list<MutationSpec<SID>>& driver_SIDs={},
                                     const std::list<CNA>& driver_CNAs={},
                                     const bool& wg_doubling=false);

    /**
     * @brief Add the properties of a mutant
     *
     * @param mutant_name is the name of the mutant
     * @param epistate_passenger_rates is a map from epigenomic state to
     *          passenger rates
     * @param driver_SIDs is a list of driver SIDs
     * @param driver_CNAs is a list of driver CNAs
     * @param application_order is the list of mutation application order
     * @return a reference to the updated object
     */
    MutationalProperties& add_mutant(const std::string& mutant_name,
                                     const std::map<std::string, PassengerRates>& epistate_passenger_rates,
                                     const std::list<MutationSpec<SID>>& driver_SIDs,
                                     const std::list<CNA>& driver_CNAs,
                                     const std::list<DriverMutations::MutationType>& application_order);

    /**
     * @brief Get the species passeger rates
     *
     * @return a constant reference to the species passeger rates
     */
    inline const std::map<std::string, PassengerRates>& get_passenger_rates() const
    {
        return passenger_rates;
    };

    /**
     * @brief Get the driver mutation map
     *
     * @return a constant reference to map associating a mutant name to its
     *      driver mutation
     */
    inline const std::map<std::string, DriverMutations>& get_driver_mutations() const
    {
        return driver_mutations;
    };

};

}   // Mutations

}   // RACES

#endif // __RACES_MUTATIONAL_PROPERTIES__
