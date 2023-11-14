/**
 * @file sample_specification.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines sample specification
 * @version 0.2
 * @date 2023-11-14
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

#ifndef __RACES_SAMPLE_SPECIFICATION__
#define __RACES_SAMPLE_SPECIFICATION__

#include "position_set.hpp"

namespace Races 
{

namespace Drivers
{

namespace Simulation
{

/**
 * @brief Representation for a sample specification
 */
class SampleSpecification
{
protected:
    std::string name;           //!< The name of the sample
    RectangleSet region;        //!< The set of tissue region to be sampled

public:
    /**
     * @brief A constructor
     * 
     * @param region is a rectangular region to sample
     */
    SampleSpecification(const RectangleSet& region);

    /**
     * @brief A constructor
     * 
     * @param name is the specification name
     * @param region is a rectangular region to sample
     */
    SampleSpecification(const std::string& name, const RectangleSet& region);

    /**
     * @brief Set the specification name
     * 
     * @param name is the new specification name
     * @return a constant reference to the updated specification name
     */
    inline const std::string& set_name(const std::string& name)
    {
        return (this->name = name);
    }

    /**
     * @brief Get the sample specification name
     * 
     * @param region is a rectangular region to sample
     * @return a constant reference to the sample specification name
     */
    inline const std::string& get_name() const
    {
        return name;
    }

    /**
     * @brief Get the sample specification region
     * 
     * @return the sample specification region
     */
    inline const RectangleSet& get_region() const
    {
        return region;
    }

    /**
     * @brief Get the default sample specification name for a region
     * 
     * @param region is a rectangular region to sample
     * @return the default name for the sample specification
     */
    static std::string get_default_name(const RectangleSet& region);

    /**
     * @brief Save a sample specification in an archive
     * 
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    void save(ARCHIVE& archive) const
    {
        archive & name & region;
    }

    /**
     * @brief Load a sample specification from an archive
     * 
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the loaded sample specification
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    static SampleSpecification load(ARCHIVE& archive)
    {
        std::string name;
        RectangleSet region;

        archive & name & region;

        return {name, region};
    }
};

}   // Simulation

}   // Drivers

}   // Races


/**
 * @brief Test the equivalence between two sample specifications
 * 
 * @param lhs is the left-hand side of the equivalence
 * @param rhs is the right-hand side of the equivalence
 * @return `true` if and only if the two sample specifications
 *      deal with the same region and they have the same 
 *      name. 
 */
inline
bool operator==(const Races::Drivers::Simulation::SampleSpecification& lhs, 
                const Races::Drivers::Simulation::SampleSpecification& rhs)
{
    return (lhs.get_name() == rhs.get_name() && lhs.get_region() == rhs.get_region());
}

#endif // __RACES_RATE_UPDATE__