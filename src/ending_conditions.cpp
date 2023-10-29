/**
 * @file ending_conditions.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements simulation ending conditions
 * @version 0.1
 * @date 2023-10-29
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

#include "ending_conditions.hpp"

namespace Races
{

namespace Drivers
{

namespace Simulation
{

TimeTest::TimeTest(const Time& threshold):
    threshold(threshold)
{}

SpeciesCountTest::SpeciesCountTest(const EpigeneticGenotypeId& species_id, const size_t& threshold):
    species_id(species_id), threshold(threshold)
{}

bool SpeciesCountTest::operator()(const Simulation& simulation) const
{
    const auto& species = simulation.tissue().get_species(species_id);

    return threshold <= species.num_of_cells();
}

uint8_t SpeciesCountTest::percentage(const Simulation& simulation) const
{
    const auto& species = simulation.tissue().get_species(species_id);

    return static_cast<uint8_t>((100*species.num_of_cells()/threshold));
}

GenotypeCountTest::GenotypeCountTest(const GenotypeId& genotype_id, const size_t& threshold):
    genotype_id(genotype_id), threshold(threshold)
{}

bool GenotypeCountTest::operator()(const Simulation& simulation) const
{
    const auto& species = simulation.tissue().get_genotype_species(genotype_id);

    return threshold <= species.num_of_cells();
}

uint8_t GenotypeCountTest::percentage(const Simulation& simulation) const
{
    const auto genotype_species = simulation.tissue().get_genotype_species(genotype_id);

    const size_t num_of_cells = genotype_species.num_of_cells();

    return static_cast<uint8_t>((100*num_of_cells/threshold));
}

}   // Simulation

}   // Drivers

}   // Races

