/**
 * @file sampling.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements samplings
 * @version 0.1
 * @date 2023-10-18
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

#include "sampling.hpp"

namespace Races 
{

namespace Drivers
{

namespace Simulation
{

Sampling::Sampling(const std::list<RectangleSet>& sample_set, const bool& remove_sample):
    sample_set(sample_set), remove_sample(remove_sample)
{}

}   // Simulation

}   // Drivers

}   // Races


bool operator==(const Races::Drivers::Simulation::Sampling& lhs, 
                const Races::Drivers::Simulation::Sampling& rhs)
{
    if (lhs.sample_set.size() != rhs.sample_set.size()
        || lhs.remove_sample != rhs.remove_sample) {
        return false;
    }

    auto lhs_it = lhs.sample_set.begin();
    auto rhs_it = rhs.sample_set.begin();

    for (;lhs_it!=lhs.sample_set.end();++lhs_it,++rhs_it) {
        if (!(*lhs_it == *rhs_it)) {
            return false;
        }
    }

    return true;
}
