/**
 * @file species.hpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Cell representation
 * @version 0.5
 * @date 2023-07-08
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

#ifndef __RACES_SPECIES__
#define __RACES_SPECIES__

#include <iostream>

#include <map>
#include <vector>
#include <random>
#include <cstddef> // std::ptrdiff_t

#include "time.hpp"
#include "cell.hpp"
#include "driver_genotype.hpp"

namespace Races {

class Tissue;

/**
 * @brief Cell species
 * 
 * This class represents the set of cells having the same driver genotype.
 */
class Species: public EpigeneticGenotype {
    std::vector<CellInTissue *> cells;   //!< species cell
    std::map<CellId, size_t> pos_map;    //!< map from CellId to `cells` position

    /**
     * @brief A constant iterator for the species cells
     */
    class const_iterator {
        std::vector<CellInTissue*>::const_iterator it; //!< the cell pointer vector iterator

        const_iterator(const std::vector<CellInTissue*>::const_iterator it);
    public:
        using difference_type   =   std::ptrdiff_t;
        using value_type        =   CellInTissue;
        using pointer           =   const CellInTissue*;
        using reference         =   const CellInTissue&;
        using iterator_category =   std::random_access_iterator_tag;

        const_iterator();

        inline reference operator*() const 
        { 
            return **it; 
        }

        inline pointer operator->() 
        {
            return *it;
        }

        // Prefix increment
        inline const_iterator& operator++() 
        {
            ++it;
            return *this;
        }  

        // Postfix increment
        const_iterator operator++(int);

        // Prefix decrement
        inline const_iterator& operator--() 
        {
            --it;
            return *this;
        }  

        // Postfix decrement
        const_iterator operator--(int);

        inline const_iterator operator+(const int delta) 
        {
            return const_iterator(it + delta);
        }

        inline const_iterator operator-(const int delta) 
        {
            return const_iterator(it - delta);
        }

        inline const_iterator& operator+=(const int& delta) {
            it += delta;

            return *this;
        }

        inline const_iterator& operator-=(const int& delta) {
            it -= delta;

            return *this;
        }

        inline reference operator[](const int& delta) const 
        {
            return *(it[delta]);
        }

        friend inline bool operator==(const const_iterator& a, const const_iterator& b)
        { 
            return a.it == b.it; 
        }

        friend inline bool operator!=(const const_iterator& a, const const_iterator& b)
        { 
            return a.it != b.it; 
        }

        friend class Species;
    };

public:

    /**
     * @brief A constructor
     * 
     * @param genotype is the epigenetic genotype of the species
     */
    Species(const EpigeneticGenotype& genotype);

    /**
     * @brief Get the initial iterator for the species cells
     * 
     * @return the initial iterator for the species cells
     */
    const_iterator begin() const;

    /**
     * @brief Get the final iterator for the species cells
     * 
     * @return the final iterator for the species cells
     */
    const_iterator end() const;

    /**
     * @brief Get the number of cells in the species
     * 
     * @return the number of cells in the species
     */
    size_t num_of_cells() const;

    /**
     * @brief Randomly select a cell in the species
     * 
     * @tparam GENERATOR is a random number generator type
     * @param generator is the random number generator used to select 
     *      the cell
     * @return a non-constant reference to a randomly selected
     *       cell in the species 
     */
    template<typename GENERATOR>
    const CellInTissue& choose_a_cell(GENERATOR& generator) const
    {
        if (num_of_cells()==0) {
            throw std::domain_error("No cells in the species");
        }

        std::uniform_int_distribution<size_t> distribution(0,num_of_cells()-1);

        size_t cell_pos = distribution(generator);

        return *(cells[cell_pos]);
    }

    /**
     * @brief Remove a cell from the species
     * 
     * @param cell_id is the id of the cell to be removed
     * @return a reference to the updated species
     */
    void remove(const CellId& cell_id);

    /**
     * @brief Add a cell to the species
     * 
     * @param cell is the cell to be added to the species
     * @return a reference to the updated species
     */
    CellInTissue* add(CellInTissue cell);

    /**
     * @brief Get cell by identifier
     * 
     * @param cell_id is the identifier of the aimed cell
     * @return a non-constant reference to the aimed cell
     */
    CellInTissue& operator()(const CellId& cell_id);

    /**
     * @brief Get cell by identifier
     * 
     * @param cell_id is the identifier of the aimed cell
     * @return a constant reference to the aimed cell
     */
    const CellInTissue& operator()(const CellId& cell_id) const;

    /**
     * @brief Test whether a cell is already contained in a species
     * 
     * @param cell_id is the identifier of the cell to be tested
     * @return `true` is an only if the cell is already in the species
     */
    bool contains(const CellId& cell_id) const;

    /**
     * @brief The destroyer
     */
    ~Species();

    template<typename LOGGER, typename PLOT_WINDOW>
    friend class BasicSimulator;
};

/**
 * @brief Write the information about a species in an output stream
 * 
 * @param out is the output stream
 * @param species is the species to be streamed
 * @return a reference to the output stream
 */
std::ostream& operator<<(std::ostream& out, const Species& species);

/* Inline Implementation */

inline bool Species::contains(const CellId& cell_id) const
{
    return pos_map.find(cell_id) != pos_map.end();
}

inline CellInTissue& Species::operator()(const CellId& cell_id)
{
    return *(cells[pos_map.at(cell_id)]);
}

inline const CellInTissue& Species::operator()(const CellId& cell_id) const
{
    return *(cells[pos_map.at(cell_id)]);
}

inline size_t Species::num_of_cells() const
{
    return cells.size();
}

};

#endif // __RACES_SPECIES__