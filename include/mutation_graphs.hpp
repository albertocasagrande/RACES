/**
 * @file mutation_graphs.hpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Represent admissible driver mutations
 * @version 0.1
 * @date 2023-05-30
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

#ifndef __RACES_MUTATION_GRAPHS__
#define __RACES_MUTATION_GRAPHS__

#include "digraph.hpp"
#include "time.hpp"

namespace Races {

/**
 * @brief Driver somatic mutation graph
 * 
 * This graph represents the somatic mutations admitted for driver.
 */
class DriverSomaticGraph : public DiGraph<Time>
{
public:
    /**
     * @brief The empty driver somatic graph constructor
     */
    DriverSomaticGraph();

    /**
     * @brief The copy constructor
     * 
     * @param orig is the template object
     */
    DriverSomaticGraph(const DriverSomaticGraph& orig);
};

/**
 * @brief Driver epigenetic mutation graph
 * 
 * This graph represents the epigenetic mutations admitted for drivers.
 */
class DriverEpigeneticGraph : public DiGraph<double>
{
public:
    /**
     * @brief The empty driver mutation graph constructor
     */
    DriverEpigeneticGraph();

    /**
     * @brief The copy constructor
     * 
     * @param orig is the template object
     */
    DriverEpigeneticGraph(const DriverEpigeneticGraph& orig);
};

};

#endif // __RACES_MUTATION_GRAPHS__