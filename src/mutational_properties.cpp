/**
 * @file mutational_properties.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements a class to represent the mutational properties
 * @version 0.13
 * @date 2023-12-19
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

MutationalProperties::MutantMutations::MutantMutations()
{}

MutationalProperties::MutantMutations::MutantMutations(const std::string& mutant_name,
                                                       const std::list<SNV>& SNVs,
                                                       const std::list<CopyNumberAlteration>& CNAs):
    name(mutant_name), SNVs(SNVs.begin(), SNVs.end()), CNAs(CNAs.begin(), CNAs.end())
{}

MutationalProperties::MutationalProperties()
{}

MutationalProperties&
MutationalProperties::add_mutant(const std::string& mutant_name,
                                 const std::map<std::string, double>& epistate_mutation_rates,
                                 const std::list<SNV>& mutant_SNVs,
                                 const std::list<CopyNumberAlteration>& mutant_CNAs)
{
    if (mutant_mutations.count(mutant_name)>0) {
        throw std::domain_error("The mutational properties of mutant \""+mutant_name+
                                "\" has been already added.");
    }

    mutant_mutations.insert({mutant_name, {mutant_name, mutant_SNVs, mutant_CNAs}});


    for (const auto& [epistate, mutation_rate] : epistate_mutation_rates) {
        auto species_name = mutant_name+epistate;

        rates[species_name] = mutation_rate;
    }

    return *this;
}

}   // Mutations

}   // Races
