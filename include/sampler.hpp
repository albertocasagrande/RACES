/**
 * @file sampler.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines classes to sample cells in a tissue
 * @version 0.5
 * @date 2023-10-02
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

#ifndef __RACES_SAMPLER__
#define __RACES_SAMPLER__

#include <cstdint>
#include <vector>

#include "tissue.hpp"

namespace Races
{

namespace Drivers
{

/**
 * @brief The basic class for tissue cell samplers 
 */
struct BasicSampler
{
    struct const_iterator {
    };

    /**
     * @brief Get the initial sampled cell iterator
     * 
     * This method must be over-loaded.
     * 
     * @return the initial sampled cell iterator
     * @throws std::runtime_error
     */
    const_iterator begin() const;

    /**
     * @brief Get the final sampled cell iterator
     * 
     * This method must be over-loaded.
     * 
     * @return the final sampled cell iterator
     * @throws std::runtime_error
     */
    const_iterator end() const;
};

/**
 * @brief A class to extract the tissue cells in a hyper-rectangle
 */
class RectangleSampler : public BasicSampler
{
    Simulation::Tissue const& tissue;           //!< The tissue          

    Simulation::PositionInTissue lower_corner;  //!< The lower corner in the sampler hyper-rectangle
    Simulation::PositionInTissue upper_corner;  //!< The upper corner in the sampler hyper-rectangle

public:
    /**
     * @brief Constant iterators for the sample cells
     */
    class const_iterator
    {
        const RectangleSampler* sampler;

        Simulation::PositionInTissue pos;

        /**
         * @brief A constructor
         * 
         * @param sampler is the sampler whose cells are traversed by the iterator
         * @param position is the position in the tissue of the new iterator
         */
        const_iterator(const RectangleSampler* sampler, const Simulation::PositionInTissue& position);
    public:
        using difference_type   =   std::ptrdiff_t;
        using value_type        =   Cell;
        using pointer           =   const Cell*;
        using reference         =   const Cell&;
        using iterator_category =   std::bidirectional_iterator_tag;

        /**
         * @brief An empty constructor
         */
        const_iterator();

        /**
         * @brief Reference operator
         * 
         * @return a reference to the species pointer by the iterator 
         */
        inline reference operator*() const 
        { 
            return sampler->tissue(pos);
        }

        /**
         * @brief Pointer operator
         * 
         * @return a pointer to the species pointer by the iterator 
         */
        inline pointer operator->() 
        {
            return &(static_cast<const Simulation::CellInTissue&>(sampler->tissue(pos)));
        }

        /**
         * @brief The prefix increment
         * 
         * @return a reference to the updated object
         */
        const_iterator& operator++();

        /**
         * @brief The postfix increment
         * 
         * @return a copy of the original object
         */
        inline const_iterator operator++(int)
        {
            RectangleSampler::const_iterator copy(*this);

            this->operator++();

            return copy;
        }

        /**
         * @brief The prefix decrement
         * 
         * @return a reference to the updated object
         */
        const_iterator& operator--();

        /**
         * @brief The postfix decrement
         * 
         * @return a copy of the original object
         */
        inline const_iterator operator--(int)
        {
            RectangleSampler::const_iterator copy(*this);

            this->operator--();

            return copy;
        }

        /**
         * @brief Test whether two iterators are the same
         * 
         * @param a is the first iterator to compare
         * @param b is the second iterator to compare
         * @return `true` if and only if the two iterators 
         *      refer to the same object
         */
        friend inline bool operator==(const const_iterator& a, const const_iterator& b)
        { 
            return (a.pos == b.pos) && (a.sampler == b.sampler); 
        }

        friend class RectangleSampler;
    };

    /**
     * @brief A hyper-rectangle sampler constructor
     * 
     * @param tissue is the tissue to sample
     * @param lower_corner is the hyper-rectangle lower corner
     * @param upper_corner is the hyper-rectangle upper corner
     */
    RectangleSampler(const Simulation::Tissue& tissue, 
                     const Simulation::PositionInTissue& lower_corner, 
                     const Simulation::PositionInTissue& upper_corner);

    /**
     * @brief A cuboid sampler constructor
     * 
     * @param tissue is the tissue to sample
     * @param lower_corner is the hyper-rectangle lower corner 
     * @param x_size is the hyper-rectangle size along the x-axis
     * @param y_size is the hyper-rectangle size along the y-axis
     * @param z_size is the hyper-rectangle size along the z-axis
     */
    RectangleSampler(const Simulation::Tissue& tissue, const Simulation::PositionInTissue& lower_corner, 
                     const Simulation::AxisSize& x_size, const Simulation::AxisSize& y_size, 
                     const Simulation::AxisSize& z_size);

    /**
     * @brief A rectangle sampler constructor
     * 
     * @param tissue is the tissue to sample
     * @param lower_corner is the hyper-rectangle lower corner 
     * @param x_size is the hyper-rectangle size along the x-axis
     * @param y_size is the hyper-rectangle size along the y-axis
     */
    RectangleSampler(const Simulation::Tissue& tissue, const Simulation::PositionInTissue& lower_corner, 
                     const Simulation::AxisSize& x_size, const Simulation::AxisSize& y_size);

    /**
     * @brief Get the initial sampled cell iterator
     * 
     * @return the initial sampled cell iterator
     */
    const_iterator begin() const;

    /**
     * @brief Get the final sampled cell iterator
     * 
     * @return the final sampled cell iterator
     */
    const_iterator end() const;

    /**
     * @brief Test whether two rectangle samplers are the same
     * 
     * @param a is the first rectangle sampler to compare
     * @param b is the second rectangle sampler to compare
     * @return `true` if and only if the two rectangle samplers 
     *      refer to the same object
     */
    friend inline bool operator==(const RectangleSampler& a, const RectangleSampler& b)
    { 
        return ((&a.tissue==&b.tissue)
                && (a.lower_corner==b.lower_corner)
                && (a.upper_corner==b.upper_corner));
    }

    friend RectangleSampler::const_iterator;
};


/**
 * @brief Test whether two iterators differs
 * 
 * @param a is the first iterator to compare
 * @param b is the second iterator to compare
 * @return `true` if and only if the two iterators 
 *      do not refer to the same object
 */
inline bool operator!=(const RectangleSampler::const_iterator& a, const RectangleSampler::const_iterator& b)
{ 
    return !(a==b); 
}

}   // Drivers

}   // Races

#endif // __RACES_SAMPLER__