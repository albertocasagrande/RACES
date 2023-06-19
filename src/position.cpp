/**
 * @file position.cpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Defines a position class in a tissue
 * @version 0.1
 * @date 2023-05-31
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

namespace Races {

PositionDelta::PositionDelta(const int x, const int y, const int z):
    x(x), y(y), z(z)
{
}

PositionDelta::PositionDelta(const Direction& direction):
    PositionDelta(0,0,0)
{
    switch(direction) {
        case Direction::X_LEFT:
            x = -1;

            break;
        case Direction::X_RIGHT:
            x = 1;

            break;
        case Direction::Y_LEFT:
            y = -1;

            break;
        case Direction::Y_RIGHT:
            y = 1;

            break;
        case Direction::Z_LEFT:
            z = -1;

            break;
        case Direction::Z_RIGHT:
            z = 1;

            break;
        default:
            throw std::domain_error("Unsupported direction");
    }
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

std::ostream& operator<<(std::ostream& os, const PositionInTissue& position)
{
    os << "(" << position.x <<","<< position.y <<","<< position.z <<")";

    return os;
}

Position::Position():
    PositionInTissue(), tissue(nullptr)
{}

Position::Position(Tissue& tissue, const AxisValue& x, const AxisValue& y, const AxisValue& z):
    PositionInTissue(x, y, z), tissue(&tissue)
{}

std::ostream& operator<<(std::ostream& os, const Position& position)
{
    os << "\""<< position.tissue->get_name() << "\""
       << static_cast<PositionInTissue>(position);

    return os;
}

};
