/**
 * @file timed_mutations.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines timed genomic mutations
 * @version 0.4
 * @date 2023-10-02
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

#ifndef __RACES_TIMED_MUTATIONS__
#define __RACES_TIMED_MUTATIONS__

#include <functional>

#include "driver_genotype.hpp"
#include "time.hpp"

namespace Races 
{

namespace Drivers
{

namespace Simulation
{

/**
 * @brief A structure to represent timed driver genomic mutation
 */
struct TimedGenomicMutation {
    GenotypeId initial_id;   //!< The initial genomic genotype identifier
    GenotypeId final_id;     //!< The final genomic genotype identifier

    Time time;                      //!< The mutation time

    /**
     * @brief A constructor
     * 
     * @param initial_id is the identifier of the initial genomic genotype
     * @param final_id is the identifier of the final genomic genotype
     * @param time is the simulation time of the genomic mutation
     */
    TimedGenomicMutation(const GenotypeId& initial_id, const GenotypeId& final_id, const Time& time);

    /**
     * @brief A constructor
     * 
     * @param initial_genotype is the initial genomic genotype
     * @param final_genotype is the final genomic genotype
     * @param time is the simulation time of the genomic mutation
     */
    TimedGenomicMutation(const Genotype& initial_genotype, const Genotype& final_genotype, const Time& time);

    /**
     * @brief Save a timed genomic mutation in an archive
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
     * @brief Load a timed genomic mutation from an archive
     * 
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the loaded timed genomic mutation
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    static TimedGenomicMutation load(ARCHIVE& archive)
    {
        GenotypeId initial_id;
        GenotypeId final_id;
        Time time;

        archive & initial_id
                & final_id
                & time;

        return TimedGenomicMutation(initial_id,final_id,time);
    }
};

}   // Simulation

}   // Drivers

}   // Races


template<>
struct std::greater<Races::Drivers::Simulation::TimedGenomicMutation> {
    inline constexpr bool operator()(const Races::Drivers::Simulation::TimedGenomicMutation &lhs, 
                                     const Races::Drivers::Simulation::TimedGenomicMutation &rhs) const 
    {
        return lhs.time > rhs.time;
    }
};

#endif // __RACES_TIMED_MUTATIONS__