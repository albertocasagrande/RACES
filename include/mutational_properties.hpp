/**
 * @file mutational_properties.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines a class to represent the mutational properties
 * @version 0.13
 * @date 2023-12-18
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

#ifndef __RACES_MUTATIONAL_PROPERTIES__
#define __RACES_MUTATIONAL_PROPERTIES__

#include <map>
#include <list>
#include <string>

#include "simulation.hpp"
#include "snv.hpp"
#include "cna.hpp"

#include "descendant_forest.hpp"
namespace Races
{

namespace Mutations
{

/**
 * @brief A class representing the mutational properties of all the species
 *
 */
class MutationalProperties
{
    /**
     * @brief The mutational properties of a species
     */
    struct SpeciesMutationalProperties
    {
        std::string name;                       //!< the species name
        double mu;                              //!< the species mutation rate
        std::list<SNV> SNVs;                    //!< the species SNVs
        std::list<CopyNumberAlteration> CNAs;   //!< the species CNAs

        /**
         * @brief The empty constructor
         */
        SpeciesMutationalProperties();

        /**
         * @brief A constructor
         *
         * @param name is the species name
         * @param SNVs is the vector of species specific SNVs
         * @param CNAs is the vector of species specific CNAs
         */
        SpeciesMutationalProperties(const std::string& name, const double& mu,
                                     const std::list<SNV>& SNVs,
                                     const std::list<CopyNumberAlteration>& CNAs);
    };

    std::map<std::string, Mutants::SpeciesId> species_names2ids;           //!< The map from species name to identifiers
    std::map<Mutants::SpeciesId, SpeciesMutationalProperties> properties;  //!< The species property map
public:

    /**
     * @brief The empty constructor
     */
    MutationalProperties();

    /**
     * @brief A constructor
     *
     * @param species_simulation is a species simulation
     */
    MutationalProperties(const Mutants::Evolutions::Simulation& species_simulation);

    /**
     * @brief A constructor
     *
     * @param descendant_forest is a descendant forest
     */
    MutationalProperties(const Mutants::DescendantsForest& descendant_forest);

    /**
     * @brief Add the properties of a mutant
     *
     * @param name is the name of the mutant
     * @param epistate_mutation_rates is a map from epigenomic status to
     *          mutational rate
     * @param mutant_SNVs is a list of SNVs characterizing the mutant
     * @param mutant_CNAs is a list of CNAs characterizing the mutant
     * @return a reference to the updated object
     */
    MutationalProperties& add_mutant(const std::string& name,
                                     const std::map<std::string, double>& epistate_mutation_rates,
                                     const std::list<SNV>& mutant_SNVs={},
                                     const std::list<CopyNumberAlteration>& mutant_CNAs={});

    /**
     * @brief Get the properties of a species
     *
     * @param species_id is the identifier of the species
     * @return the mutational properties of the species having
     *          `species_id` as identifier.
     */
    inline const MutationalProperties::SpeciesMutationalProperties& at(const Mutants::SpeciesId& species_id) const
    {
        return properties.at(species_id);
    }

    /**
     * @brief Get the species mutational properties
     *
     * @return a constant reference to the species mutational properties
     */
    inline const std::map<Mutants::SpeciesId, SpeciesMutationalProperties>&
    get_species_properties() const
    {
        return properties;
    };
};

}   // Mutations

}   // Races

#endif // __RACES_MUTATIONAL_PROPERTIES__
