/**
 * @file driver_genotype.hpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Driver genotype representation
 * @version 0.2
 * @date 2023-06-28
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

#ifndef __RACES_DRIVER_GENOTYPE__
#define __RACES_DRIVER_GENOTYPE__

#include <map>
#include <string>

#include "cell_event.hpp"

namespace Races {

#define NULL_MUTATION std::numeric_limits<DriverGenotypeId>::max() 

/**
 * @brief Driver genotype class
 */
class DriverGenotype {
    static unsigned int counter;   //!< Total number of driver along the computation
    DriverGenotypeId id;           //!< Driver identifier

    std::string name;              //!< Driver genotype name
    std::map<CellEventType, double> rates;    //!< Event rates

    bool methylated;               //!< Methylation flag
public:

    /**
     * @brief Construct a new DriverGenotype
     * 
     * @param name is the name of the driver genotype
     * @param rates is a map from the events to rates
     * @param methylated is a methylation flag
     */
    DriverGenotype(const std::string& name,
           const std::map<CellEventType, double>& rates,
           const bool methylated=false);

    /**
     * @brief Copy constructor
     * 
     * @param orig is the original driver genotype
     */
    DriverGenotype(const DriverGenotype& orig);

    /**
     * @brief Copy operator
     * 
     * @param orig is the original driver genotype
     * @return a reference to the updated object
     */
    DriverGenotype& operator=(const DriverGenotype& orig);

    /**
     * @brief Move operator
     * 
     * @param orig is the original driver genotype
     * @return a reference to the updated object
     */
    DriverGenotype& operator=(DriverGenotype&& orig);

    /**
     * @brief Get the driver genotype identifier
     * 
     * @return a constant reference to the driver identifier
     */
    const DriverGenotypeId& get_id() const;

    /**
     * @brief Get the driver genotype name
     * 
     * @return a constant reference to the driver identifier
     */
    const std::string& get_name() const;

    /**
     * @brief Get the driver genotype event rates
     * 
     * @return a constant reference to the genotype event rates
     */
    const std::map<CellEventType, double>& get_rates() const;

    /**
     * @brief Get the rate of an event 
     * 
     * @param event is the event whose rate is required
     * @return a constant reference to event rate
     */
    const double& get_rate(const CellEventType& event) const;

    /**
     * @brief Test whether this is the methylated version of the genotype
     * 
     * @return true if and only if this object represent 
     *        the methylated version of the driver
     */
    const bool& is_methylated() const;
};

/**
 * @brief Write the information about a driver genotype in an output stream
 * 
 * @param out is the output stream
 * @param driver is the driver genotype to be streamed
 * @return a reference to the output stream
 */
std::ostream& operator<<(std::ostream& out, const DriverGenotype& driver);

/* Inline methods definitions */

inline const DriverGenotypeId& DriverGenotype::get_id() const
{
    return id;
}

inline const std::string& DriverGenotype::get_name() const
{
    return name;
}

inline const std::map<CellEventType, double>& DriverGenotype::get_rates() const
{
    return rates;
}

inline const double& DriverGenotype::get_rate(const CellEventType& event) const
{
    return rates.at(event);
}

inline const bool& DriverGenotype::is_methylated() const
{
    return methylated;
}

};

#endif // __RACES_DRIVER_GENOTYPE__