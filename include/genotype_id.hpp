/**
 * @file genotype_id.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines genotype ans species identifiers
 * @version 0.3
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

#ifndef __RACES_GENOTYPE_ID__
#define __RACES_GENOTYPE_ID__

#include <limits>
#include <cstdint>

namespace Races 
{

/**
 * @brief The namespace of classes handing clones
 */
namespace Clones
{

/**
 * @brief The genotype identifier type
 */
using GenotypeId = uint16_t;

/**
 * @brief The species identifier type
 */
using SpeciesId = uint16_t;


}   // Clones

}   // Races

/**
 * @brief An identifier for wild-type species
 * 
 * This macro represents the species identifier of the wild-type cells.
 */
#define WILD_TYPE_SPECIES std::numeric_limits<Races::Clones::SpeciesId>::max()

#endif // __RACES_GENOTYPE_ID__