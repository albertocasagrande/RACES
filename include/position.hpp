/**
 * @file position.hpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Defines a position class in a tissue
 * @version 0.1
 * @date 2023-05-30
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

namespace Races {

class Tissue;

/**
 * @brief Directions
 */
enum class Direction {
    X_LEFT,
    X_RIGHT,
    Y_LEFT,
    Y_RIGHT,
    Z_LEFT,
    Z_RIGHT
};

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
    PositionDelta(const Direction& direction);

    /**
     * @brief Negate position delta;
     * 
     */
    PositionDelta operator-() const;
};

/**
 * @brief A 3D position in a tissues
 * 
 */
struct PositionInTissue {
    unsigned int x;  //!< x axis
    unsigned int y;  //!< y axis
    unsigned int z;  //!< z axis

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
    PositionInTissue(const unsigned int x, const unsigned int y, const unsigned int z=0);

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
 * @brief Stream the position in output
 * 
 * @param os is the output stream
 * @param position is the position to be streamed
 * @return a reference to the output stream
 */
std::ostream& operator<<(std::ostream& os, const PositionInTissue& position);

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
    Position(Tissue& tissue, const unsigned int& x, const unsigned int& y, const unsigned int& z);

    /**
     * @brief Construct a new Position object
     * 
     * @param tissue is the tissue referred by the position
     * @param position is the position in the tissue
     */
    Position(Tissue& tissue, const PositionInTissue& position);
};

/**
 * @brief Stream the position in output
 * 
 * @param os is the output stream
 * @param position is the position to be streamed
 * @return a reference to the output stream
 */
std::ostream& operator<<(std::ostream& os, const Position& position);

};

#endif // __RACES_POSITION__
