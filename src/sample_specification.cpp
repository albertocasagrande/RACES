/**
 * @file sample_specification.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements sample specification
 * @version 0.4
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

#include <sstream>

#include "sample_specification.hpp"

namespace Races 
{

namespace Mutants
{

namespace Evolutions
{

std::string SampleSpecification::get_default_name(const RectangleSet& region)
{
    std::ostringstream oss;

    oss << "S_" << region.lower_corner << "-" << region.upper_corner; 

    return oss.str();
}

SampleSpecification::SampleSpecification(const RectangleSet& region):
    name(get_default_name(region)), region(region)
{}

SampleSpecification::SampleSpecification(const std::string& name, const RectangleSet& region):
    name(name), region(region)
{}

}   // Evolutions

}   // Mutants

}   // Races
