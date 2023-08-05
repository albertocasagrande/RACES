/**
 * @file position.hpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Defines a position class in a tissue
 * @version 0.4
 * @date 2023-08-05
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

#ifndef __RACES_POSITION__
#define __RACES_POSITION__

#include <iostream>
#include <type_traits>

namespace Races 
{

namespace Drivers 
{

namespace Simulation 
{

class Tissue;

/**
 * @brief Directions
 */
struct Direction {
    uint8_t bit_vector; 

    int get_delta(const size_t& axis_index) const;

public:
    enum Value {
        X_NULL = 0x00,
        X_DOWN = 0x01,
        X_UP = 0x02,
        Y_NULL = 0x00 << 2,
        Y_DOWN = 0x01 << 2,
        Y_UP = 0x02 << 2,
        Z_NULL = 0x00 << 4,
        Z_DOWN = 0x01 << 4,
        Z_UP = 0x02 << 4
    };

    Direction(const uint8_t value);

    int get_delta_x() const;

    int get_delta_y() const;

    int get_delta_z() const;

    friend std::ostream& operator<<(std::ostream& os, const Direction& direction);
};

inline
int Direction::get_delta_x() const
{
    return get_delta(0);
}

inline
int Direction::get_delta_y() const
{
    return get_delta(1);
}

inline
int Direction::get_delta_z() const
{
    return get_delta(2);
}

struct PositionDelta {
    int x;  //!< x axis
    int y;  //!< y axis
    int z;  //!< z axis

    /**
     * @brief Create a new position delta
     * 
     * @param x is the x-axis delta
     * @param y is the y-axis delta
     * @param z is the z-axis delta
     */
    PositionDelta(const int x, const int y, const int z);

    /**
     * @brief Create a new position delta towards a direction 
     * 
     * @param direction 
     */
    explicit PositionDelta(const Direction& direction);

    /**
     * @brief Negate position delta;
     * 
     */
    PositionDelta operator-() const;
};

using AxisValue = uint16_t;

/**
 * @brief A 3D position in a tissues
 * 
 */
struct PositionInTissue {
    int16_t x;  //!< x axis
    int16_t y;  //!< y axis
    int16_t z;  //!< z axis

    /**
     * @brief A constructor
     */
    PositionInTissue();

    /**
     * @brief Construct a new Position In Tissue object
     * 
     * @param x is the x-axis position
     * @param y is the y-axis position
     * @param z is the z-axis position
     */
    PositionInTissue(const AxisValue x, const AxisValue y, const AxisValue z=0);

    /**
     * @brief Add a delta to the position
     * 
     * @param delta is the delta to be added
     * @return a reference to the updated position
     */
    PositionInTissue& operator+=(const PositionDelta delta);

    /**
     * @brief Subtract a delta from a position 
     * 
     * @param delta is the delta to be subtracted to
     * @return a reference to the updated position
     */
    PositionInTissue& operator-=(const PositionDelta delta);
};

/**
 * @brief Add a delta to a position 
 * 
 * @param position is the position to which a delta should be added to
 * @param delta is the delta to be added to
 * @return the position reachable from `position` adding `delta`
 */
PositionInTissue operator+(const PositionInTissue& position, const PositionDelta delta);

/**
 * @brief Subtract a delta from a position 
 * 
 * @param position is the position to which a delta should be subtracted to
 * @param delta is the delta to be subtracted to
 * @return the position reachable from `position` subtracting `delta`
 */
PositionInTissue operator-(const PositionInTissue& position, const PositionDelta delta);

/**
 * @brief Check whether two positions are the same
 * 
 * @param p1 is a position
 * @param p2 is a position
 * @return true if and only if the two positions have the same coordinates
 */
bool operator==(const PositionInTissue& p1, const PositionInTissue& p2);

/**
 * @brief Check whether two positions differ
 * 
 * @param p1 is a position
 * @param p2 is a position
 * @return true if and only if the two positions have different coordinates
 */
bool operator!=(const PositionInTissue& p1, const PositionInTissue& p2);

/**
 * @brief Compute the Manhattan distance between two positions
 * 
 * @param p1 is a position
 * @param p2 is a position
 * @return the Manhattan distance
 */
size_t Manhattan_distance(const PositionInTissue& p1, const PositionInTissue& p2);

/**
 * @brief The class representing cell position
 */
struct Position : public PositionInTissue
{
    Tissue* tissue;         //!< tissue

    /**
     * @brief An empty constructor
     */
    Position();

    /**
     * @brief Construct a new Position object
     * 
     * @param tissue is the tissue referred by the position
     * @param x is the x axis position in the tissue
     * @param y is the y axis position in the tissue
     * @param z is the z axis position in the tissue
     */
    Position(Tissue& tissue, const AxisValue& x, const AxisValue& y, const AxisValue& z);

    /**
     * @brief Construct a new Position object
     * 
     * @param tissue is the tissue referred by the position
     * @param position is the position in the tissue
     */
    Position(Tissue& tissue, const PositionInTissue& position);
};


}   // Simulation

}   // Drivers

}   // Races


namespace std
{

std::ostream& operator<<(std::ostream& os, const Races::Drivers::Simulation::Direction& direction);

std::ostream& operator<<(std::ostream& os, const Races::Drivers::Simulation::PositionDelta& delta);

/**
 * @brief Stream the position in output
 * 
 * @param os is the output stream
 * @param position is the position to be streamed
 * @return a reference to the output stream
 */
std::ostream& operator<<(std::ostream& os, const Races::Drivers::Simulation::PositionInTissue& position);

/**
 * @brief Stream the position in output
 * 
 * @param os is the output stream
 * @param position is the position to be streamed
 * @return a reference to the output stream
 */
std::ostream& operator<<(std::ostream& os, const Races::Drivers::Simulation::Position& position);

}   // std

#endif // __RACES_POSITION__
