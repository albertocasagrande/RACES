/**
 * @file clone.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Define the Python wrapper class and functions for `CloneProperties`
 * @version 0.1
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

#ifndef __RACES_PYTHON_CLONE__
#define __RACES_PYTHON_CLONE__

#include <memory>

#include <boost/python.hpp>

#include "clone_properties.hpp"


using namespace boost::python;

/**
 * @brief A wrapper structure for `CloneProperties`
 * 
 * This structure defines all the static methods required by 
 * the Python wrapping of `CloneProperties`.
 */
struct CloneWrapper 
{

    /**
     * @brief Create a new `Clone` object
     * 
     * @param name is the name of the new object
     * @param epigenetic_rates is the Python list of epigenetic rates
     *         being either pairs of doubles of object of the class
     *         `EpigeneticRates`
     * @return a shared pointer of the newly created object
     */
    static std::shared_ptr<Races::Clones::CloneProperties>
    create(std::string const& name, boost::python::list const& epigenetic_rates);

    /**
     * @brief Set the rates of the species
     * 
     * Every clone is associated to many different species.
     * The species is a clone equipped with a epigenetic status
     * (i.e., methylated/demethylated) that is represented by the methylation 
     * signature (i.e., "+"/"-"). This method set the rates of a species 
     * by specifying its clone, its methylation signature, and the new 
     * rates. 
     * 
     * @param clone is the clone whose rate should be set
     * @param methylation_signature is a string representing the methylation signature
     * @param rates are the new rates 
     */
    static
    void set_rates(Races::Clones::CloneProperties *clone, const std::string& methylation_signature, 
                boost::python::dict const& rates);

    /**
     * @brief Get the rates of the species
     * 
     * Every clone is associated to many different species.
     * The species is a clone equipped with a epigenetic status
     * (i.e., methylated/demethylated) that is represented by the methylation 
     * signature (i.e., "+"/"-"). This method get the rates of a species 
     * by specifying its clone, its methylation signature, and the 
     * event type.
     * 
     * @param clone is the clone whose rate should be set
     * @param methylation_signature is a string representing the methylation signature
     * @param 
     */
    inline static
    double get_rate(const Races::Clones::CloneProperties *clone, 
                    const std::string& methylation_signature,
                    const Races::Clones::CellEventType event_type)
    {
        return (*clone)[methylation_signature].get_rate(event_type);
    }

};

#endif // __RACES_PYTHON_CLONE__