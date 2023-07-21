/**
 * @file position.cpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Defines a position class in a tissue
 * @version 0.3
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

#include "position.hpp"
#include "tissue.hpp"

namespace Races 
{

namespace Drivers 
{

namespace Simulation 
{

Direction::Direction(const uint8_t value):
    bit_vector(value)
{
    for (size_t i=0; i<3; ++i) {
        const uint8_t axis = (value >> 2*i)&0x03;

        if (axis == 0x03) {
            throw std::domain_error("Unsupported direction");
        }
    }
}


int Direction::get_delta(const size_t& axis_index) const
{
    const uint8_t axis = (bit_vector>>2*axis_index)&0x03;

    if (axis == 0x00) {
        return 0;
    }

    if (axis == 0x01) {
        return 1;
    }

    return -1;
}


PositionDelta::PositionDelta(const int x, const int y, const int z):
    x(x), y(y), z(z)
{
}

PositionDelta::PositionDelta(const Direction& direction):
    PositionDelta(direction.get_delta_x(),direction.get_delta_y(),direction.get_delta_z())
{
}

PositionDelta PositionDelta::operator-() const
{
    return PositionDelta{-this->x,-this->y,-this->z};
}

PositionInTissue::PositionInTissue():
    x(0), y(0), z(0)
{}

PositionInTissue::PositionInTissue(const AxisValue x, const AxisValue y, const AxisValue z):
    x(x), y(y), z(z)
{}

PositionInTissue& PositionInTissue::operator+=(const PositionDelta delta)
{
    x += delta.x;
    y += delta.y;
    z += delta.z;

    return *this;
}

PositionInTissue& PositionInTissue::operator-=(const PositionDelta delta)
{
    x -= delta.x;
    y -= delta.y;
    z -= delta.z;

    return *this;
}


PositionInTissue operator+(const PositionInTissue& position, const PositionDelta delta)
{
    PositionInTissue res{position};

    res += delta;

    return res;
}


PositionInTissue operator-(const PositionInTissue& position, const PositionDelta delta)
{
    PositionInTissue res{position};

    res -= delta;

    return res;  
}

bool operator==(const PositionInTissue& p1, const PositionInTissue& p2)
{
    return p1.x==p2.x && p1.y==p2.y && p1.z==p2.z;
}


bool operator!=(const PositionInTissue& p1, const PositionInTissue& p2)
{
    return !(p1==p2);
}

size_t Manhattan_distance(const PositionInTissue& p1, const PositionInTissue& p2)
{
    size_t distance{0};

    distance += (p1.x>p2.x ? p1.x-p2.x : p2.x-p1.x);
    distance += (p1.y>p2.y ? p1.y-p2.y : p2.y-p1.y);
    distance += (p1.z>p2.z ? p1.z-p2.z : p2.z-p1.z);

    return distance;
}

Position::Position():
    PositionInTissue(), tissue(nullptr)
{}

Position::Position(Tissue& tissue, const AxisValue& x, const AxisValue& y, const AxisValue& z):
    PositionInTissue(x, y, z), tissue(&tissue)
{}

Position::Position(Tissue& tissue, const PositionInTissue& pos):
    PositionInTissue(pos), tissue(&tissue)
{}

}   // Simulation

}   // Drivers

}   // Races

namespace std 
{

std::ostream& operator<<(std::ostream& os, const Races::Drivers::Simulation::Direction& direction)
{
    for (size_t i=0; i<3; ++i) {
        switch(direction.get_delta(i)) {
            case 1:
                os << "U";
                break;
            case -1:
                os << "D";
                break;
            case 0:
                os << "N";
                break;
            default:
                throw std::runtime_error("Unknown direction");
        }
    }

    return os;
}

std::ostream& operator<<(std::ostream& os, const Races::Drivers::Simulation::PositionDelta& delta)
{
    os << "(" << delta.x << "," << delta.y << "," << delta.z << ")";

    return os;
}

std::ostream& operator<<(std::ostream& os, const Races::Drivers::Simulation::PositionInTissue& position)
{
    os << "(" << position.x <<","<< position.y <<","<< position.z <<")";

    return os;
}

std::ostream& operator<<(std::ostream& os, const Races::Drivers::Simulation::Position& position)
{
    os << "\""<< position.tissue->get_name() << "\""
       << static_cast<Races::Drivers::Simulation::PositionInTissue>(position);

    return os;
}

}  // std
