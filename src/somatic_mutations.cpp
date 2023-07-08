/**
 * @file somatic_mutations.cpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Implementing timed somatic mutations
 * @version 0.1
 * @date 2023-07-08
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

#include "somatic_mutations.hpp"

namespace Races {

TimedSomaticMutation::TimedSomaticMutation(const SomaticGenotypeId& initial_id, const SomaticGenotypeId& final_id, const Time& time):
    initial_id(initial_id), final_id(final_id), time(time)
{}

TimedSomaticMutation::TimedSomaticMutation(const SomaticGenotype& initial_genotype, const SomaticGenotype& final_genotype, const Time& time):
    TimedSomaticMutation(initial_genotype.get_id(), final_genotype.get_id(), time)
{}

}
