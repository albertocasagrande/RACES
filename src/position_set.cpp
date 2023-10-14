/**
 * @file position_set.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements classes to represent tissue position set
 * @version 0.1
 * @date 2023-10-14
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

#include "position_set.hpp"

namespace Races
{

namespace Drivers
{

BasicPositionSet::const_iterator BasicPositionSet::begin() const
{
    throw std::runtime_error("The BasicPositionSet class must be inherited");
}

BasicPositionSet::const_iterator BasicPositionSet::end() const
{
    throw std::runtime_error("The BasicPositionSet class must be inherited");
}

RectangleSet::const_iterator::const_iterator(const RectangleSet* rectangle,
                                             const Simulation::PositionInTissue& position):
    rectangle(rectangle), pos(position)
{
}

RectangleSet::const_iterator::const_iterator():
    rectangle(nullptr)
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

void set_to_invalid_position(Simulation::PositionInTissue& pos, const RectangleSet& rectangle)
{
    pos = rectangle.upper_corner;
    ++pos.x;
}

inline bool is_invalid_in(Simulation::PositionInTissue& pos, const RectangleSet& rectangle)
{
    return ((pos.x<rectangle.lower_corner.x || pos.x>rectangle.upper_corner.x
            || pos.y<rectangle.lower_corner.y || pos.y>rectangle.upper_corner.y
            || pos.z<rectangle.lower_corner.z || pos.z>rectangle.upper_corner.z));
}

RectangleSet::const_iterator& RectangleSet::const_iterator::operator++()
{
    if (rectangle==nullptr) {
        return *this;
    }

    if (is_invalid_in(pos, *rectangle)) {
        throw std::runtime_error("The iterator has already reached the end of the set");
    }

    if (increase_and_update(pos.z, rectangle->lower_corner.z, rectangle->upper_corner.z)) {
        if (increase_and_update(pos.y, rectangle->lower_corner.y, rectangle->upper_corner.y)) {
            if (increase_and_update(pos.x, rectangle->lower_corner.x, rectangle->upper_corner.x)) {
                set_to_invalid_position(pos, *rectangle);

                return *this;
            }
        }
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

RectangleSet::const_iterator& RectangleSet::const_iterator::operator--()
{
    if (rectangle==nullptr) {
        return *this;
    }

    if (is_invalid_in(pos, *rectangle)) {
        throw std::runtime_error("The iterator has already reached the end of the set");
    }

    if (decrease_and_update(pos.z, rectangle->lower_corner.z, rectangle->upper_corner.z)) {
        if (decrease_and_update(pos.y, rectangle->lower_corner.y, rectangle->upper_corner.y)) {
            if (decrease_and_update(pos.x, rectangle->lower_corner.x, rectangle->upper_corner.x)) {
                set_to_invalid_position(pos, *rectangle);

                return *this;
            }
        }
    }

    return *this;
}

RectangleSet::RectangleSet(const Simulation::PositionInTissue& position):
    RectangleSet(position, position)
{}

RectangleSet::RectangleSet(const Simulation::PositionInTissue& lower_corner, 
                           const Simulation::PositionInTissue& upper_corner):
    lower_corner(lower_corner), upper_corner(upper_corner)
{
}

RectangleSet::RectangleSet(const Simulation::PositionInTissue& lower_corner,
                            const Simulation::AxisSize& x_size,
                            const Simulation::AxisSize& y_size,
                            const Simulation::AxisSize& z_size):
    lower_corner(lower_corner), 
    upper_corner(lower_corner.x+x_size-1,lower_corner.y+y_size-1,lower_corner.z+z_size-1)
{
}

RectangleSet::RectangleSet(const Simulation::PositionInTissue& lower_corner,
                           const Simulation::AxisSize& x_size,
                           const Simulation::AxisSize& y_size):
    lower_corner(lower_corner), 
    upper_corner(lower_corner.x+x_size-1,lower_corner.y+y_size-1)
{
    if (lower_corner.z!=0) { 
        throw std::domain_error("The lower corner must be a 2D position.");
    }
}

RectangleSet::const_iterator RectangleSet::begin() const
{
    return const_iterator(this, lower_corner);
}

RectangleSet::const_iterator RectangleSet::end() const
{
    const_iterator it(this, upper_corner);

    return ++it;
}

}   // Drivers

}   // Races
