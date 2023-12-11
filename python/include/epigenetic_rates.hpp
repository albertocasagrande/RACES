/**
 * @file epigenetic_rates.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Define the Python wrapper class and functions for `Races::EpigeneticRates`
 * @version 0.6
 * @date 2023-12-11
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

#ifndef __RACES_PYTHON_EPIGENETIC_RATES__
#define __RACES_PYTHON_EPIGENETIC_RATES__

#include <memory>

#include <boost/python.hpp>

#include "mutant_properties.hpp"


using namespace boost::python;

/**
 * @brief A wrapper structure for `Races::EpigeneticRates`
 */
struct EpigeneticRatesWrapper 
{
    /**
     * @brief Create a new `Races::EpigeneticRates` object
     * 
     * This method creates a new `Races::EpigeneticRates` object
     * by using as parameters the values stored in Python list.
     * If the list contains 1 or 2 real values, then this 
     * method creates the object. Otherwise, it throws a 
     * domain exception.
     * 
     * @param rates_list is a Python list of object
     * @return a shared pointer to the newly created 
     *        `Races::EpigeneticRates` object
     */
    static std::shared_ptr<Races::Mutants::EpigeneticRates>
    create(boost::python::list const& rates_list);

    /**
     * @brief Set the methylation rate of a `Races::EpigeneticRates` object
     * 
     * @param epigenetic_rates is the `Races::EpigeneticRates` object 
     * @param value is the methylation rate value to be set
     */
    inline static
    void set_methylation_rate(Races::Mutants::EpigeneticRates *epigenetic_rates, const double& value)
    {
        epigenetic_rates->set_methylation_rate(value);
    }

    /**
     * @brief Set the demethylation rate of a `Races::EpigeneticRates` object
     * 
     * @param epigenetic_rates is the `Races::EpigeneticRates` object 
     * @param value is the demethylation rate value to be set
     */
    inline static
    void set_demethylation_rate(Races::Mutants::EpigeneticRates *epigenetic_rates, const double& value)
    {
        epigenetic_rates->set_demethylation_rate(value);
    }
};

#endif // __RACES_PYTHON_EPIGENETIC_RATES__
