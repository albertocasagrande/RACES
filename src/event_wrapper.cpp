/**
 * @file event_wrapper.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implement a simulation event wrapper
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

#include "event_wrapper.hpp"

namespace Races 
{

namespace Clones
{

namespace Evolutions
{

SimulationEventWrapper::SimulationEventWrapper(const GenotypeMutation& genotype_mutation):
    type(SimulationEvent::Type::GENOTYPE_MUTATION),
    event(std::make_shared<GenotypeMutation>(genotype_mutation))
{}

SimulationEventWrapper::SimulationEventWrapper(const RateUpdate& rate_update):
    type(SimulationEvent::Type::LIVENESS_RATE_UPDATE),
    event(std::make_shared<RateUpdate>(rate_update))
{}

SimulationEventWrapper::SimulationEventWrapper(const Sampling& sampling):
    type(SimulationEvent::Type::SAMPLING),
    event(std::make_shared<Sampling>(sampling))
{}

}   // Evolutions

}   // Clones

}   // Races


bool operator==(const Races::Clones::Evolutions::SimulationEventWrapper& lhs, 
                const Races::Clones::Evolutions::SimulationEventWrapper& rhs)
{
    using namespace Races::Clones::Evolutions;

    if (lhs.type != rhs.type) {
        return false;
    }

    switch(lhs.type) {
        case SimulationEvent::Type::GENOTYPE_MUTATION:
            {
                const auto& m_lhs = lhs.get_event<GenotypeMutation>();
                const auto& m_rhs = rhs.get_event<GenotypeMutation>();

                return (m_lhs == m_rhs);
            }
        case SimulationEvent::Type::LIVENESS_RATE_UPDATE:
            {
                const auto& ru_lhs = lhs.get_event<RateUpdate>();
                const auto& ru_rhs = rhs.get_event<RateUpdate>();

                return (ru_lhs == ru_rhs);
            }
        default:
            throw std::runtime_error("Unsupported event type");
    }
}