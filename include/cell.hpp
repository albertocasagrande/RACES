/**
 * @file cell.hpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Cell representation
 * @version 0.3
 * @date 2023-06-28
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

#ifndef __RACES_CELL__
#define __RACES_CELL__

#include <limits>
#include <cstdint>

#include "time.hpp"
#include "position.hpp"
#include "driver_genotype.hpp"

#define NON_DRIVER_GENOTYPE std::numeric_limits<DriverGenotypeId>::max()

namespace Races {

using CellId = uint64_t;

class Cell {
protected:
    static uint64_t counter;   //!< total number of cell along the computation
    CellId id;                 //!< cell identifier
    CellId parent;             //!< parent cell

    DriverGenotypeId genotype;     //!< cell species reference

    /**
     * @brief The empty constructor
     */
    Cell();

public:
    /**
     * @brief Create a new cell
     * 
     * @param genotype is the driver genotype identifier
     * @param parent_id is the parent cell identifier
     */
    Cell(const DriverGenotypeId genotype, const CellId parent_id=0);

    /**
     * @brief Get the cell identifier
     * 
     * @return a constant reference to the cell identifier
     */
    const CellId& get_id() const;

    /**
     * @brief Get the cell parent identifier
     * 
     * @return a constant reference to the cell parent identifier
     */
    const CellId& get_parent_id() const;

    /**
     * @brief Get the cell driver genotype
     * 
     * @return a constant reference to the cell driver genotype
     */
    const DriverGenotypeId& get_driver_genotype() const;

    /**
     * @brief Clone a cell
     * 
     * @param cell is the template cell
     * @return a reference to the updated object
     */
    Cell& clone(const Cell& cell);

    friend class Tissue;
    friend class Species;

    template<typename LOGGER, typename PLOT_WINDOW>
    friend class BasicSimulator; 

    friend void swap(Cell& a, Cell &b);
};

/**
 * @brief Swap two cells
 * 
 * @param a is a cell
 * @param b is a cell
 */
void swap(Cell& a, Cell &b);

/**
 * @brief Write a cell in an output stream
 * 
 * @param os is the output stream
 * @param cell is the cell to be streamed
 * @return a reference to the updated stream
 */
std::ostream& operator<<(std::ostream& os, const Cell& cell);

class Species;

/**
 * @brief A class to represent a cell in a tissue
 * 
 * A cell in a tissue is a cell provided with a position 
 * in the tissue
 */
class CellInTissue : public Cell, public PositionInTissue {

protected:
    /**
     * @brief A cell in tissue constructor
     * 
     * This constructor is meant to build non-driver-genotype cells
     * 
     * @param tissue is the tissue referred by the position
     * @param x is the x axis position in the tissue
     * @param y is the y axis position in the tissue
     * @param z is the z axis position in the tissue
     */
    CellInTissue(const AxisValue& x, const AxisValue& y, const AxisValue& z);

    /**
     * @brief Copy a cell into a cell in a tissue
     * 
     * @param position is the position in the tissue
     * @return a reference to the updated object
     */
    CellInTissue& operator=(const PositionInTissue& position);

public:

    /**
     * @brief Copy a cell into a cell in a tissue
     * 
     * @param cell is the original free cell
     * @return a reference to the updated object
     */
    CellInTissue& operator=(const Cell& cell);

    /**
     * @brief Test whether the cell has a driver genotype
     * 
     * @return true if and only if the cell has a driver genotype
     */
    bool has_driver_genotype() const;

    friend class Tissue;
};


/**
 * @brief Write a cell in a tissue in an output stream
 * 
 * @param os is the output stream
 * @param cell is the cell in the tissue to be streamed
 * @return a reference to the updated stream
 */
std::ostream& operator<<(std::ostream& os, const CellInTissue& cell);

/* Inline Implementation */

inline const CellId& Cell::get_id() const
{
    return id;
}

inline const CellId& Cell::get_parent_id() const
{
    return parent;
}


inline const DriverGenotypeId& Cell::get_driver_genotype() const
{
    return genotype;
}

inline bool CellInTissue::has_driver_genotype() const {
    return get_driver_genotype() != NON_DRIVER_GENOTYPE;
}

}

#endif // __RACES_CELL__