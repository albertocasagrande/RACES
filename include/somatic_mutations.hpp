/**
 * @file somatic_mutations.hpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Implementing timed somatic mutations
 * @version 0.2
 * @date 2023-07-13
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

#ifndef __RACES_SOMATIC_MUTATIONS__
#define __RACES_SOMATIC_MUTATIONS__

#include <functional>

#include "driver_genotype.hpp"
#include "time.hpp"

namespace Races {

/**
 * @brief A structure to represent timed somatic mutation
 */
struct TimedSomaticMutation {
    SomaticGenotypeId initial_id;   //!< The initial somatic genotype identifier
    SomaticGenotypeId final_id;     //!< The final somatic genotype identifier

    Time time;                      //!< The mutation time

    /**
     * @brief A constructor
     * 
     * @param initial_id is the identifier of the initial somatic genotype
     * @param final_id is the identifier of the final somatic genotype
     * @param time is the simulation time of the somatic mutation
     */
    TimedSomaticMutation(const SomaticGenotypeId& initial_id, const SomaticGenotypeId& final_id, const Time& time);

    /**
     * @brief A constructor
     * 
     * @param initial_id is the initial somatic genotype
     * @param final_id is the final somatic genotype
     * @param time is the simulation time of the somatic mutation
     */
    TimedSomaticMutation(const SomaticGenotype& initial_genotype, const SomaticGenotype& final_genotype, const Time& time);

    /**
     * @brief Save a timed somatic mutation in an archive
     * 
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        archive & initial_id
                & final_id
                & time;
    }

    /**
     * @brief Load a timed somatic mutation from an archive
     * 
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the loaded timed somatic mutation
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    static TimedSomaticMutation load(ARCHIVE& archive)
    {
        SomaticGenotypeId initial_id;
        SomaticGenotypeId final_id;
        Time time;

        archive & initial_id
                & final_id
                & time;

        return TimedSomaticMutation(initial_id,final_id,time);
    }
};

};

template<>
struct std::less<Races::TimedSomaticMutation> {
    inline constexpr bool operator()(const Races::TimedSomaticMutation &lhs, const Races::TimedSomaticMutation &rhs) const 
    {
        return lhs.time < rhs.time;
    }
};

#endif // __RACES_SOMATIC_MUTATIONS__