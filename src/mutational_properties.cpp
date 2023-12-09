/**
 * @file mutational_properties.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements a class to represent the mutational properties
 * @version 0.7
 * @date 2023-12-09
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

#include "mutational_properties.hpp"

namespace Races
{

namespace Mutations
{

MutationalProperties::SpeciesMutationalProperties::SpeciesMutationalProperties()
{}

MutationalProperties::SpeciesMutationalProperties::SpeciesMutationalProperties(const std::string& name, const double& mu,
                                                                                        const std::list<SNV>& SNVs, 
                                                                                        const std::list<CopyNumberAlteration>& CNAs):
    name(name), mu(mu), SNVs(SNVs), CNAs(CNAs)
{}

MutationalProperties::MutationalProperties()
{}

MutationalProperties& 
MutationalProperties::add_genotype(const Clones::Evolutions::Simulation& clones_simulation,  
                                         const std::string& name,
                                         const std::map<std::string, double>& epigenetic_rates,
                                         const std::list<SNV>& species_SNVs,
                                         const std::list<CopyNumberAlteration>& species_CNAs)
{
    std::map<std::string, const Clones::Evolutions::Species*> species_map;

    for (const auto& species : clones_simulation.tissue()) {
        species_map[species.get_name()] = &species;
    }

    for (const auto& [epigenetic_status, epigenetic_rate] : epigenetic_rates) {
        auto full_name = name+epigenetic_status;
        auto species_it = species_map.find(full_name);

        if (species_it == species_map.end()) {
            throw std::domain_error("Unknown species name: \""+full_name+"\"");
        }

        Clones::SpeciesId id = (species_it->second)->get_id();

        auto prop_it = properties.find(id);

        if (prop_it != properties.end()) {
            throw std::domain_error("The mutational properties of species \""+full_name+
                                    "\" has been already added.");
        }

        properties[id] = {full_name, epigenetic_rate, species_SNVs, species_CNAs};
    }

    return *this;
}

}   // Mutations

}   // Races
