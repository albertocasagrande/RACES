/**
 * @file mutation.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implement a general class for mutations
 * @version 1.0
 * @date 2024-06-10
 *
 * @copyright Copyright (c) 2023-2024
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

#include "mutation.hpp"

namespace RACES
{

namespace Mutations
{

Mutation::Mutation():
    GenomicPosition(), nature(UNDEFINED), cause()
{}

Mutation::Mutation(const ChromosomeId& chr_id, const ChrPosition& chr_position,
                   const Nature& nature, const std::string& cause):
    GenomicPosition(chr_id, chr_position), nature(nature), cause(cause)
{}

Mutation::Mutation(const GenomicPosition& position, const Nature& nature,
                   const std::string& cause):
    Mutation(position.chr_id, position.position, nature, cause)
{}

std::string Mutation::get_nature_description(const Mutation::Nature& nature)
{
    switch(nature) {
      case Mutation::Nature::DRIVER:
        return "driver";
      case Mutation::Nature::PASSENGER:
        return "passenger";
      case Mutation::Nature::PRENEOPLASTIC:
        return "pre-neoplastic";
      case Mutation::Nature::GERMINAL:
        return "germinal";
      default:
        return "unknown";
    }
}


}   // Mutations

}   // RACES
