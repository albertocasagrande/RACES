/**
 * @file ending_conditions.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines simulation ending conditions
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

#ifndef __RACES_ENDING_CONDITIONS__
#define __RACES_ENDING_CONDITIONS__

#include "simulation.hpp"

namespace Races
{

namespace Drivers
{

namespace Simulation
{

/**
 * @brief Simulation time test
 * 
 * The objects of this class have a method testing the time 
 * simulated by a simulation. If it is below a threshold,
 * then the test returns `false`; otherwise, it returns `true`.
 */
struct TimeTest : public Simulation::BasicTest
{
    Time threshold;

    /**
     * @brief A constructor
     * 
     * @param threshold is the threshold for simulation time test
     */
    TimeTest(const Time& threshold);

    /**
     * @brief Test whether the simulated time is not below the test threshold
     * 
     * @param simulation is the considered simulation
     * @return `false` if and only if the time simulated by `simulation`
     *          is lower than the test threshold
     */
    inline bool operator()(const Simulation& simulation) const override
    {
        return threshold <= simulation.get_time();
    }

    /**
     * @brief Return the percentage of the completed simulation
     * 
     * @param simulation is the considered simulation
     * @return the percentage of the completed simulation
     */
    inline uint8_t percentage(const Simulation& simulation) const override
    {
        return static_cast<uint8_t>((100*simulation.get_time()/threshold));
    }
};

/**
 * @brief Species count test
 * 
 * The objects of this class have a method testing the number of 
 * cells in a species. If this number is below a threshold,
 * then the test returns `false`; otherwise, it returns `true`.
 */
struct SpeciesCountTest : public Simulation::BasicTest
{
    EpigeneticGenotypeId species_id;
    size_t threshold;

    /**
     * @brief A constructor
     * 
     * @param species_id is the identifier of the species whose 
     *          number of cells is counted
     * @param threshold is the threshold for the count test
     */
    SpeciesCountTest(const EpigeneticGenotypeId& species_id, const size_t& threshold);

    /**
     * @brief Test whether the number of cells is below the threshold
     * 
     * @param simulation is the considered simulation
     * @return `false` if and only if the number of cells of the 
     *          considered species is below the test threshold
     */
    bool operator()(const Simulation& simulation) const override;

    /**
     * @brief Return the percentage of the completed simulation
     * 
     * @param simulation is the considered simulation
     * @return the percentage of the completed simulation
     */
    uint8_t percentage(const Simulation& simulation) const override;
};

/**
 * @brief Genotype count test
 * 
 * The objects of this class have a method testing the number of 
 * cells of a genotype. If this number is below a threshold,
 * then the test returns `false`; otherwise, it returns `true`.
 */
struct GenotypeCountTest : public Simulation::BasicTest
{
    GenotypeId genotype_id;
    size_t threshold;

    /**
     * @brief A constructor
     * 
     * @param genotype_id is the identifier of the genotype whose 
     *          number of cells is counted
     * @param threshold is the threshold for the count test
     */
    GenotypeCountTest(const GenotypeId& genotype_id, const size_t& threshold);

    /**
     * @brief Test whether the number of cells is below the threshold
     * 
     * @param simulation is the considered simulation
     * @return `false` if and only if the number of cells of the 
     *          considered genotype is below the test threshold
     */
    bool operator()(const Simulation& simulation) const override;

    /**
     * @brief Return the percentage of the completed simulation
     * 
     * @param simulation is the considered simulation
     * @return the percentage of the completed simulation
     */
    uint8_t percentage(const Simulation& simulation) const override;
};

}   // Simulation

}   // Drivers

}   // Races

#endif // __RACES_ENDING_CONDITIONS__
