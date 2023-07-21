/**
 * @file genotype.hpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Define the Python wrapper class and functions for `Races::Genotype`
 * @version 0.2
 * @date 2023-07-21
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

#ifndef __RACES_PYTHON_GENOMIC_GENOTYPE__
#define __RACES_PYTHON_GENOMIC_GENOTYPE__

#include <memory>

#include <boost/python.hpp>

#include "driver_genotype.hpp"


using namespace boost::python;

/**
 * @brief A wrapper structure for `Races::Genotype`
 * 
 * This structure defines all the static methods required by 
 * the Python wrapping of `Races::Genotype`.
 */
struct GenotypeWrapper 
{

    /**
     * @brief Create a new `Genotype` object
     * 
     * @param name is the name of the new object
     * @param epigenetic_rates is the Python list of epigenetic rates
     *         being either pairs of doubles of object of the class
     *         `EpigeneticRates`
     * @return a shared pointer of the newly created object
     */
    static std::shared_ptr<Races::Drivers::Genotype>
    create(std::string const& name, boost::python::list const& epigenetic_rates);

    /**
     * @brief Set the rates of the epigenetic genotype
     * 
     * Every genomic genotype is associated to many different epigenetic genotype.
     * The epigenetic genotype is a genomic genotype equipped with a epigenetic status
     * (i.e., methylated/demethylated) that is represented by the methylation 
     * signature (i.e., "+"/"-"). This method set the rates of an epigenetic genotype 
     * by specifying its genomic genotype, its methylation signature, and the new 
     * rates. 
     * 
     * @param genotype is the genomic genotype whose rate should be set
     * @param methylation_signature is a string representing the methylation signature
     * @param rates are the new rates 
     */
    static
    void set_rates(Races::Drivers::Genotype *genotype, const std::string& methylation_signature, 
                boost::python::dict const& rates);

    /**
     * @brief Get the rates of the epigenetic genotype
     * 
     * Every genomic genotype is associated to many different epigenetic genotype.
     * The epigenetic genotype is a genomic genotype equipped with a epigenetic status
     * (i.e., methylated/demethylated) that is represented by the methylation 
     * signature (i.e., "+"/"-"). This method get the rates of an epigenetic genotype 
     * by specifying its genomic genotype, its methylation signature, and the 
     * event type.
     * 
     * @param genotype is the genomic genotype whose rate should be set
     * @param methylation_signature is a string representing the methylation signature
     * @param 
     */
    inline static
    double get_rate(const Races::Drivers::Genotype *genotype, 
                    const std::string& methylation_signature,
                    const Races::Drivers::CellEventType event_type)
    {
        return (*genotype)[methylation_signature].get_rate(event_type);
    }

};

#endif // __RACES_PYTHON_GENOMIC_GENOTYPE__