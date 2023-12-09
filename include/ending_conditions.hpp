/**
 * @file ending_conditions.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines simulation ending conditions
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

#ifndef __RACES_ENDING_CONDITIONS__
#define __RACES_ENDING_CONDITIONS__

#include "simulation.hpp"

namespace Races
{

namespace Clones
{

namespace Evolutions
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
    const Time threshold; //!< The time threshold

    /**
     * @brief A constructor
     * 
     * @param threshold is the threshold for simulation time test
     */
    explicit TimeTest(const Time& threshold);

    /**
     * @brief Test whether the simulated time is not below the test threshold
     * 
     * @param simulation is the considered simulation
     * @return `false` if and only if the time simulated by `simulation`
     *          is lower than the test threshold
     */
    inline bool operator()(const Simulation& simulation) override
    {
        return threshold < simulation.get_time();
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
    const SpeciesId species_id; //!< The species id
    const size_t threshold;  //!< The size threshold

    /**
     * @brief A constructor
     * 
     * @param species_id is the identifier of the species whose 
     *          number of cells is counted
     * @param threshold is the threshold for the count test
     */
    SpeciesCountTest(const SpeciesId& species_id, const size_t& threshold);

    /**
     * @brief Test whether the number of cells is below the threshold
     * 
     * @param simulation is the considered simulation
     * @return `false` if and only if the number of cells of the 
     *          considered species is below the test threshold
     */
    bool operator()(const Simulation& simulation) override;

    /**
     * @brief Return the percentage of the completed simulation
     * 
     * @param simulation is the considered simulation
     * @return the percentage of the completed simulation
     */
    uint8_t percentage(const Simulation& simulation) const override;
};

/**
 * @brief Clone count test
 * 
 * The objects of this class have a method testing the number of 
 * cells of a clone. If this number is below a threshold,
 * then the test returns `false`; otherwise, it returns `true`.
 */
struct CloneCountTest : public Simulation::BasicTest
{
    const CloneId clone_id;  //!< The clone id
    const size_t threshold;  //!< The size threshold

    /**
     * @brief A constructor
     * 
     * @param clone_id is the identifier of the clone whose 
     *          number of cells is counted
     * @param threshold is the threshold for the count test
     */
    CloneCountTest(const CloneId& clone_id, const size_t& threshold);

    /**
     * @brief Test whether the number of cells is below the threshold
     * 
     * @param simulation is the considered simulation
     * @return `false` if and only if the number of cells of the 
     *          considered clone is below the test threshold
     */
    bool operator()(const Simulation& simulation) override;

    /**
     * @brief Return the percentage of the completed simulation
     * 
     * @param simulation is the considered simulation
     * @return the percentage of the completed simulation
     */
    uint8_t percentage(const Simulation& simulation) const override;
};

/**
 * @brief Event count test
 * 
 * The objects of this class have a method testing the number of 
 * cells stochastic events occurred to a species. If this number 
 * is below a threshold, then the test returns `false`; otherwise, 
 * it returns `true`.
 */
class EventCountTest : public Simulation::BasicTest
{
    /**
     * @brief Get the number of events tested by this object
     * 
     * @param simulation is the considered simulation
     * @return the number of events tested by this object
     */
    size_t get_event_number(const Simulation& simulation) const;

    const CellEventType event_type; //!< The cell event
    const SpeciesId species_id;     //!< The id of the species in which the events occurred
    const SpeciesId dst_id;         //!< The id of the destination species in case of epigenetic events
    const size_t threshold;         //!< The threshold on the number of events

public:

    /**
     * @brief A constructor
     * 
     * Whenever `event_type` is `CellEventType::EPIGENETIC_SWITCH`
     * the built object compare the threshold and the sum of all 
     * the epigenetic events.
     * 
     * @param event_type is the aimed event type
     * @param species_id is the identifier of the species whose 
     *          number of events is counted
     * @param threshold is the threshold for the number of events
     */
    EventCountTest(const CellEventType& event_type, const SpeciesId& species_id, 
                   const size_t& threshold);

    /**
     * @brief A constructor
     * 
     * This constructor builds an object specifically designed to 
     * test whether the number of epigenetic event from `src_id` 
     * to `dst_id` has reached a threshold.
     * 
     * @param src_id is the identifier of the source species of 
     *          counted epigenetic events
     * @param dst_id is the identifier of the destination species of 
     *          counted epigenetic events
     * @param threshold is the threshold for the number of events
     */
    EventCountTest(const SpeciesId& src_id, const SpeciesId& dst_id,
                   const size_t& threshold);

    /**
     * @brief Test whether the number of cells is below the threshold
     * 
     * @param simulation is the considered simulation
     * @return `false` if and only if the number of cells of the 
     *          considered species is below the test threshold
     */
    bool operator()(const Simulation& simulation) override;

    /**
     * @brief Return the percentage of the completed simulation
     * 
     * @param simulation is the considered simulation
     * @return the percentage of the completed simulation
     */
    uint8_t percentage(const Simulation& simulation) const override;
};

}   // Evolutions

}   // Clones

}   // Races

#endif // __RACES_ENDING_CONDITIONS__
