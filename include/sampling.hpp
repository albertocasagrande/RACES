/**
 * @file sampling.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines tissue samplings
 * @version 0.2
 * @date 2023-10-25
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

#ifndef __RACES_SAMPLING__
#define __RACES_SAMPLING__

#include "position_set.hpp"
#include "simulation_event.hpp"

namespace Races 
{

namespace Drivers
{

namespace Simulation
{

/**
 * @brief A structure to represent tissue sampling event
 */
struct Sampling : public SimulationEvent
{
    using Type = SimulationEvent::Type;

    std::list<RectangleSet> sample_set; //!< The set of tissue region to be sampled
    bool preserve_tissue;               //!< The flag establishing whether the tissue must remain unchanged

    /**
     * @brief A constructor
     * 
     * @param sample_set is a list of rectangular region to sample
     * @param preserve_tissue is a Boolean flag to establish whether the sample 
     *          must be remove from the tissue
     */
    Sampling(const std::list<RectangleSet>& sample_set, const bool& preserve_tissue=true);

    /**
     * @brief Save a timed genomic mutation in an archive
     * 
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    void save(ARCHIVE& archive) const
    {
        archive & sample_set & preserve_tissue;
    }

    /**
     * @brief Load a timed genomic mutation from an archive
     * 
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the loaded timed genomic mutation
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    static Sampling load(ARCHIVE& archive)
    {
        std::list<RectangleSet> sample_set;
        bool preserve_tissue;

        archive & sample_set & preserve_tissue;

        return {sample_set, preserve_tissue};
    }

    inline Type type() const {
        return Type::SAMPLING;
    } 
};

}   // Simulation

}   // Drivers

}   // Races


/**
 * @brief Test the equivalence between two sampling
 * 
 * @param lhs is the left-hand side of the equivalence
 * @param rhs is the right-hand side of the equivalence
 * @return `true` if and only if the two liveness rate updates represent
 *      the same event
 */
bool operator==(const Races::Drivers::Simulation::Sampling& lhs, 
                const Races::Drivers::Simulation::Sampling& rhs);

#endif // __RACES_RATE_UPDATE__