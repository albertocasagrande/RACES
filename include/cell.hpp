/**
 * @file cell.hpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Cell representation
 * @version 0.8
 * @date 2023-07-09
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

#include "archive.hpp"
#include "time.hpp"
#include "position.hpp"
#include "driver_genotype.hpp"

#define NON_DRIVER_GENOTYPE std::numeric_limits<EpigeneticGenotypeId>::max()

namespace Races {

using CellId = uint64_t;

class Cell {
protected:
    static uint64_t counter;            //!< total number of cell along the computation
    CellId id;                          //!< cell identifier
    CellId parent;                      //!< parent cell

    EpigeneticGenotypeId genotype;          //!< cell species reference

    unsigned int passenger_mutations;   //!< number of passenger mutations

    /**
     * @brief The empty constructor
     */
    Cell();

public:
    /**
     * @brief Create a new cell with no passenger mutations
     * 
     * @param genotype is the driver genotype identifier
     */
    Cell(const EpigeneticGenotypeId genotype);

    /**
     * @brief Create a new cell
     * 
     * @param genotype is the driver genotype identifier
     * @param passenger_mutations is the number of passenger mutations
     */
    Cell(const EpigeneticGenotypeId genotype, unsigned int passenger_mutations);

    /**
     * @brief Create a new cell
     * 
     * @param genotype is the driver genotype identifier
     * @param passenger_mutations is the number of passenger mutations
     * @param parent_id is the parent cell identifier
     */
    Cell(const EpigeneticGenotypeId genotype, unsigned int passenger_mutations, const CellId parent_id);

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
    const EpigeneticGenotypeId& get_genotype_id() const;

    /**
     * @brief Get the number of passenger mutations
     * 
     * @return a constant reference to the number of passenger mutations
     */
    const unsigned int& get_passenger_mutations() const;

    /**
     * @brief Generate a descendent cell
     * 
     * @param new_passenger_mutations is the number of new passenger mutations
     * @return the generated cell
     */
    Cell generate_descendent(const unsigned int new_passenger_mutations=0) const;

    /**
     * @brief Save a cell in an archive
     * 
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Out::Basic, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        archive & id 
                & parent 
                & genotype 
                & passenger_mutations;
    }

    /**
     * @brief Load a cell from an archive
     * 
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the load cell
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::In::Basic, ARCHIVE>, bool> = true>
    inline static Cell load(ARCHIVE& archive)
    {
        Cell cell;

        archive & cell.id 
                & cell.parent 
                & cell.genotype 
                & cell.passenger_mutations;

        return cell;
    }

    friend class Tissue;
    friend class Species;

    template<typename LOGGER>
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
     * @brief The empty constructor
     */
    CellInTissue();

    /**
     * @brief A cell in tissue constructor
     * 
     * @param genotype is the cell driver genotype id
     * @param passenger_mutations is the number of passenger mutations
     * @param position is the cell position
     */
    CellInTissue(const EpigeneticGenotypeId genotype, unsigned int passenger_mutations, const PositionInTissue& position);

    /**
     * @brief A cell in tissue constructor
     * 
     * @param cell is the cell
     * @param position is the cell position
     */
    CellInTissue(const Cell& cell, const PositionInTissue& position);

    /**
     * @brief A cell in tissue constructor
     * 
     * This constructor is meant to build non-driver-genotype cells.
     * 
     * @param position is the cell position
     */
    CellInTissue(const PositionInTissue& position);

    /**
     * @brief A cell in tissue constructor
     * 
     * This constructor is meant to build non-driver-genotype cells.
     * 
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
    bool has_driver_mutations() const;

    /**
     * @brief Save a cell in an archive
     * 
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Out::Basic, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        archive & static_cast<const Cell &>(*this) & x & y & z;
    }

    /**
     * @brief Load a cell from an archive
     * 
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the loaded cell
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::In::Basic, ARCHIVE>, bool> = true>
    static CellInTissue load(ARCHIVE& archive)
    {
        CellInTissue cell;

        archive & static_cast<Cell &>(cell) 
                & cell.x & cell.y & cell.z;

        return cell;
    }

    friend class Tissue;
    friend class Species;
    friend class CellInTissueProxy;
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

inline const EpigeneticGenotypeId& Cell::get_genotype_id() const
{
    return genotype;
}

inline const unsigned int& Cell::get_passenger_mutations() const
{
    return passenger_mutations;
}

inline Cell Cell::generate_descendent(const unsigned int new_passenger_mutations) const
{
    return Cell(get_genotype_id(), passenger_mutations+new_passenger_mutations, get_id());
}

inline bool CellInTissue::has_driver_mutations() const {
    return get_genotype_id() != NON_DRIVER_GENOTYPE;
}

}

#endif // __RACES_CELL__