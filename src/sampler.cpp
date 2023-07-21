/**
 * @file sampler.cpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Tissue sampler class implementation
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

#include "sampler.hpp"

namespace Races
{

namespace Drivers
{

BasicSampler::const_iterator BasicSampler::begin() const
{
    throw std::runtime_error("The BasicSampler class must be inherited");
}

BasicSampler::const_iterator BasicSampler::end() const
{
    throw std::runtime_error("The BasicSampler class must be inherited");
}

RectangleSampler::const_iterator::const_iterator(const RectangleSampler* sampler,
                                                 const Simulation::PositionInTissue& position):
    sampler(sampler), pos(position)
{
    if (!sampler->tissue(pos).has_driver_mutations()) {
        this->operator++();
    }
}

RectangleSampler::const_iterator::const_iterator():
    sampler(nullptr)
{
}

template<typename T>
bool increase_and_update(T& value, const T& min, const T& max)
{
    if (value >= max) {
        value = min;

        return true;
    } 
    
    ++value;
    return false;
}

void set_to_invalid_position(Simulation::PositionInTissue& pos, const Simulation::Tissue& tissue)
{
    auto sizes = tissue.size();
    pos.x = sizes[0];
    pos.y = sizes[1];
    if (sizes.size()==2) {
        pos.z = 0;
    } else {
        pos.z = sizes[2];
    }
}

RectangleSampler::const_iterator& RectangleSampler::const_iterator::operator++()
{
    if (sampler==nullptr) {
        return *this;
    }

    bool found_a_cell{false};
    while (!found_a_cell) {
        if (increase_and_update(pos.z, sampler->lower_corner.z, sampler->upper_corner.z)) {
            if (increase_and_update(pos.y, sampler->lower_corner.y, sampler->upper_corner.y)) {
                if (increase_and_update(pos.x, sampler->lower_corner.x, sampler->upper_corner.x)) {
                    set_to_invalid_position(pos, sampler->tissue);

                    return *this;
                }
            }
        }

        found_a_cell = sampler->tissue(pos).has_driver_mutations();
    }

    return *this;
}

template<typename T>
bool decrease_and_update(T& value, const T& min, const T& max)
{
    if (value <= min) {
        value = max;

        return true;
    } 
    
    --value;
    return false;
}

RectangleSampler::const_iterator& RectangleSampler::const_iterator::operator--()
{
    if (sampler==nullptr) {
        return *this;
    }

    bool found_a_cell{false};
    while (found_a_cell) {
        if (decrease_and_update(pos.z, sampler->lower_corner.z, sampler->upper_corner.z)) {
            if (decrease_and_update(pos.y, sampler->lower_corner.y, sampler->upper_corner.y)) {
                if (decrease_and_update(pos.x, sampler->lower_corner.x, sampler->upper_corner.x)) {
                    set_to_invalid_position(pos, sampler->tissue);
        
                    return *this;
                }
            }
        }

        found_a_cell = sampler->tissue(pos).has_driver_mutations();
    }

    return *this;
}

RectangleSampler::RectangleSampler(const Simulation::Tissue& tissue, 
                                   const Simulation::PositionInTissue& lower_corner, 
                                   const Simulation::PositionInTissue& upper_corner):
    tissue(tissue), lower_corner(lower_corner), upper_corner(upper_corner)
{
    if (tissue.num_of_dimensions()==2 && lower_corner.z!=0) { 
        throw std::domain_error("The tissue is a 2D space. The lower corner must be a 2D position.");
    }
    if (tissue.num_of_dimensions()==2 &&  upper_corner.z!=0) { 
        throw std::domain_error("The tissue is a 2D space. The upper corner must be a 2D position.");
    }

    if (!tissue.is_valid(lower_corner)) {
        throw std::domain_error("The lower corner is not valid for the tissue");
    }

    if (!tissue.is_valid(upper_corner)) {
        throw std::domain_error("The upper corner is not valid for the tissue");
    }
}

RectangleSampler::RectangleSampler(const Simulation::Tissue& tissue,
                                   const Simulation::PositionInTissue& lower_corner,
                                   const Simulation::AxisSize& x_size,
                                   const Simulation::AxisSize& y_size,
                                   const Simulation::AxisSize& z_size):
    tissue(tissue), lower_corner(lower_corner), 
    upper_corner(lower_corner.x+x_size-1,lower_corner.y+y_size-1,lower_corner.z+z_size-1)
{
    if (tissue.num_of_dimensions()==2) { 
        throw std::domain_error("The tissue is a 2D space. The sampler "
                                "hyper-rectangle must be a 2D rectangle.");
    }

    if (!tissue.is_valid(lower_corner)) {
        throw std::domain_error("The lower corner is not valid for the tissue");
    }

    if (!tissue.is_valid(upper_corner)) {
        throw std::domain_error("The rectangle sizes are not compatible with the tissue sizes");
    }
}

RectangleSampler::RectangleSampler(const Simulation::Tissue& tissue,
                                   const Simulation::PositionInTissue& lower_corner,
                                   const Simulation::AxisSize& x_size,
                                   const Simulation::AxisSize& y_size):
    tissue(tissue), lower_corner(lower_corner), 
    upper_corner(lower_corner.x+x_size-1,lower_corner.y+y_size-1)
{
    if (tissue.num_of_dimensions()==3) { 
        throw std::domain_error("The tissue is a 3D space. The sampler "
                                "hyper-rectangle sizes must be 3 in number.");
    }
    if (lower_corner.z!=0) { 
        throw std::domain_error("The lower corner must be a 2D position.");
    }

    if (!tissue.is_valid(lower_corner)) {
        throw std::domain_error("The lower corner is not valid for the tissue");
    }

    if (!tissue.is_valid(upper_corner)) {
        throw std::domain_error("The rectangle sizes are not compatible with the tissue sizes");
    }
}

RectangleSampler::const_iterator RectangleSampler::begin() const
{
    return const_iterator(this, lower_corner);
}

RectangleSampler::const_iterator RectangleSampler::end() const
{
    const_iterator it(this, upper_corner);

    return ++it;
}

}   // Drivers

}   // Races
