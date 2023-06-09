/**
 * @file tissue.hpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Define tissue class
 * @version 0.11
 * @date 2023-07-12
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

#ifndef __RACES_TISSUE__
#define __RACES_TISSUE__

#include <vector>
#include <list>
#include <map>
#include <set>
#include <string>
#include <memory> // SpeciesView::const_iterator

#include "archive.hpp"
#include "time.hpp"
#include "species.hpp"
#include "cell.hpp"

namespace Races {

using AxisSize = uint16_t;

class Tissue {
    using SomaticGenotypePosition = std::map<SomaticGenotypeId, std::vector<size_t>>;

    std::string name;                               //!< The tissue name
    std::vector<Species> species;                   //!< The species in the tissue
    std::map<EpigeneticGenotypeId, size_t> pos_map; //!< The identifier to position map
    SomaticGenotypePosition somatic_genotope_pos;   //!< The positions of the species associated to the same somatic genotype 
    
    uint8_t dimensions;                             //!< The number of space dimension for the tissue

    std::vector<std::vector<std::vector<CellInTissue *>>> space;     //!< Space in the tissue

    /**
     * @brief Get the pointer to the cell in a position
     * 
     * This method returns a pointer to a cell in a given position.
     * The returned value is undefined when the position is not 
     * a valid position for the tissue.
     * 
     * @param position is the position of the aimed cell pointer
     * @return a non-constant reference to a pointer to the cell 
     *          in `position`
     */
    CellInTissue*& cell_pointer(const PositionInTissue& position);

    /**
     * @brief Get the pointer to the cell in a position
     * 
     * This method returns a pointer to a cell in a given position.
     * The returned value is undefined when the position is not 
     * a valid position for the tissue.
     * 
     * @param position is the position of the aimed cell pointer
     * @return a constant reference to a pointer to the cell 
     *          in `position`
     */
    const CellInTissue* cell_pointer(const PositionInTissue& position) const;

public:

    /**
     * @brief A view class for species
     * 
     * This class allows to have a partial view of the 
     * tissue's species. For instance, it allows to 
     * list all the species having the same somatic 
     * genotype.
     */
    class SpeciesView {
        const std::vector<Species>& species;       //!< The species vector
        const std::vector<size_t>& species_pos;    //!< The species position vector

        /**
         * @brief A constructor
         * 
         * @param species is a reference to the tissue's species vector
         * @param species_pos is a vector of valid positions for `vector`
         */
        SpeciesView(const std::vector<Species>& species, const std::vector<size_t>& species_pos);
    public:

        /**
         * @brief A constant iterator over species views
         */
        class const_iterator {
            const std::vector<Species>* species;    //!< The tissue's species vector
            std::shared_ptr<std::vector<size_t>::const_iterator> it; //!< A constant iterator over a vector of position in `species`

            /**
             * @brief A private constructor
             * 
             * @param species is the vector of the tissue's species
             * @param it is a constant iterator over a vector of the positions in `species`
             */
            const_iterator(const std::vector<Species>& species, const std::vector<size_t>::const_iterator it);
        public:
            using difference_type   =   std::ptrdiff_t;
            using value_type        =   Species;
            using pointer           =   const Species*;
            using reference         =   const Species&;
            using iterator_category =   std::random_access_iterator_tag;

            /**
             * @brief An empty construtor
             */
            const_iterator();

            /**
             * @brief Reference operator
             * 
             * @return a reference to the species pointer by the iterator 
             */
            inline reference operator*() const 
            { 
                return (*species)[**it]; 
            }

            /**
             * @brief Pointer operator
             * 
             * @return a pointer to the species pointer by the iterator 
             */
            inline pointer operator->() 
            {
                return &((*species)[**it]);
            }

            /**
             * @brief The prefix increment
             * 
             * @return a reference to the updated object
             */
            inline const_iterator& operator++() 
            {
                ++(*it);
                return *this;
            }  

            /**
             * @brief The postfix increment
             * 
             * @return a copy of the original object
             */
            const_iterator operator++(int);

            /**
             * @brief The prefix decrement
             * 
             * @return a reference to the updated object
             */
            inline const_iterator& operator--() 
            {
                --(*it);
                return *this;
            }  

            /**
             * @brief The postfix decrement
             * 
             * @return a copy of the original object
             */
            const_iterator operator--(int);

            /**
             * @brief Add operator
             * 
             * @param delta is the value to add
             * @return a new iterator that points `delta` position ahead 
             *      with respect to the original object 
             */
            inline const_iterator operator+(const int delta) 
            {
                return const_iterator(*species, *it + delta);
            }

            /**
             * @brief Subtract operator
             * 
             * @param delta is the value to subtract
             * @return a new iterator that points `delta` position backwards 
             *      with respect to the original object 
             */
            inline const_iterator operator-(const int delta) 
            {
                return const_iterator(*species, *it - delta);
            }

            /**
             * @brief Inplace add operator
             * 
             * @param delta is the value to add
             * @return a reference to the update object  
             */
            inline const_iterator& operator+=(const int& delta) {
                (*it) += delta;

                return *this;
            }

            /**
             * @brief Inplace subtract operator
             * 
             * @param delta is the value to subtract
             * @return a reference to the update object  
             */
            inline const_iterator& operator-=(const int& delta) {
                (*it) -= delta;

                return *this;
            }

            /**
             * @brief Index operator
             * 
             * @param index is the index to access 
             * @return a reference to the `index`-th objects after 
             *      that pointer by the current iterator 
             */
            inline reference operator[](const int& index) const 
            {
                return (*species)[(*it)[index]];
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
                return (*a.it == *b.it) && (a.species == b.species); 
            }

            friend class Tissue::SpeciesView;
        };

        /**
         * @brief Index operator
         * 
         * @param index is the index to access
         * @return a constant reference to the `index`-th species 
         *      in the view
         */
        inline const Species& operator[](const size_t& index) const
        {
            return species[species_pos[index]];
        }

        /**
         * @brief Get the view size
         * 
         * @return the view size
         */
        inline size_t size() const
        {
            return species_pos.size();
        }

        /**
         * @brief Get the view begin
         * 
         * @return a constant iterator to the first element 
         *      in the view 
         */
        inline const_iterator begin() const
        {
            return const_iterator(species, species_pos.begin());
        }

        /**
         * @brief Get the view end
         * 
         * @return a constant iterator to the view end
         */
        inline const_iterator end() const
        {
            return const_iterator(species, species_pos.end());
        }

        friend class Tissue;
    };

    /**
     * @brief This class wraps pointer to constant cells in tissue space
     */
    template<typename TISSUE_TYPE>
    class BaseCellInTissueProxy
    {
    protected:
        TISSUE_TYPE &tissue;                //!< tissue
        const PositionInTissue position;    //!< position of the cell

        BaseCellInTissueProxy(TISSUE_TYPE &tissue, const PositionInTissue position);
    public:

        bool has_driver_mutations() const;

        operator const CellInTissue&() const;

        friend class Tissue;
    };

    /**
     * @brief This class wraps pointer to constant cells in tissue space
     */    
    class CellInTissueConstantProxy : public BaseCellInTissueProxy<const Tissue>
    {
        CellInTissueConstantProxy(const Tissue &tissue, const PositionInTissue position);

        friend class Tissue;
    };

    /**
     * @brief This class wraps pointer to cells in tissue space
     */
    class CellInTissueProxy : public BaseCellInTissueProxy<Tissue>
    {
        CellInTissueProxy(Tissue &tissue, const PositionInTissue position);
    public:
        CellInTissueProxy& operator=(const Cell& cell);

        operator CellInTissue&();

        void kill();

        CellInTissue copy_and_kill();

        friend class Tissue;
    };

    /**
     * @brief A constructor
     * 
     * @param sizes are the sizes of the tissue
     */
    Tissue(const std::vector<AxisSize> sizes);

    /**
     * @brief A constructor for a 3D tissue
     * 
     * @param x_size is the size of the tissue on the x axis
     * @param y_size is the size of the tissue on the y axis
     * @param z_size is the size of the tissue on the z axis
     */
    Tissue(const AxisSize x_size, const AxisSize  y_size, const AxisSize z_size);

    /**
     * @brief A constructor for a 2D tissue
     * 
     * @param x_size is the size of the tissue on the x axis
     * @param y_size is the size of the tissue on the y axis
     */
    Tissue(const AxisSize x_size, const AxisSize y_size);

    /**
     * @brief A constructor for a 3D tissue
     * 
     * @param genotypes is the vector of driver genotypes
     * @param x_size is the size of the tissue on the x axis
     * @param y_size is the size of the tissue on the y axis
     * @param z_size is the size of the tissue on the z axis
     */
    Tissue(const std::vector<SomaticGenotype> genotypes, const AxisSize  x_size, const AxisSize  y_size, const AxisSize z_size);

    /**
     * @brief A constructor for a 2D tissue
     * 
     * @param genotypes is the vector of driver genotypes
     * @param x_size is the size of the tissue on the x axis
     * @param y_size is the size of the tissue on the y axis
     */
    Tissue(const std::vector<SomaticGenotype> genotypes, const AxisSize  x_size, const AxisSize  y_size);

    /**
     * @brief A constructor
     * 
     * @param name is the tissue name
     * @param sizes are the sizes of the tissue
     */
    Tissue(const std::string name, const std::vector<AxisSize> sizes);

    /**
     * @brief A constructor for a 3D tissue
     * 
     * @param name is the tissue name
     * @param x_size is the size of the tissue on the x axis
     * @param y_size is the size of the tissue on the y axis
     * @param z_size is the size of the tissue on the z axis
     */
    Tissue(const std::string name, const AxisSize x_size, const AxisSize  y_size, const AxisSize  z_size);

    /**
     * @brief A constructor for a 2D tissue
     * 
     * @param name is the tissue name
     * @param x_size is the size of the tissue on the x axis
     * @param y_size is the size of the tissue on the y axis
     */
    Tissue(const std::string name, const AxisSize x_size, const AxisSize  y_size);

    /**
     * @brief A constructor for a 3D tissue
     * 
     * @param name is the tissue name
     * @param genotypes is the vector of driver genotypes
     * @param x_size is the size of the tissue on the x axis
     * @param y_size is the size of the tissue on the y axis
     * @param z_size is the size of the tissue on the z axis
     */
    Tissue(const std::string name, const std::vector<SomaticGenotype> genotypes, const AxisSize  x_size, const AxisSize  y_size, const AxisSize z_size);

    /**
     * @brief A constructor for a 2D tissue
     * 
     * @param name is the tissue name
     * @param genotypes is the vector of driver genotypes
     * @param x_size is the size of the tissue on the x axis
     * @param y_size is the size of the tissue on the y axis
     */
    Tissue(const std::string name, const std::vector<SomaticGenotype> genotypes, const AxisSize  x_size, const AxisSize  y_size);

    /**
     * @brief Get the initial iterator for the tissue species
     * 
     * @return the initial iterator for the tissue species
     */
    std::vector<Species>::const_iterator begin() const;

    /**
     * @brief Get the final iterator for the tissue species
     * 
     * @return the final iterator for the tissue species
     */
    std::vector<Species>::const_iterator end() const;

    /**
     * @brief Get a tissue species by driver identifier
     * 
     * @param genotype_id is the driver genotype identifier
     * @return a constant reference to the tissue species
     */
    const Species& get_species(const EpigeneticGenotypeId& genotype_id) const;

    /**
     * @brief Get a tissue species by driver identifier
     * 
     * @param genotype_id is the driver genotype identifier
     * @return a non-constant reference to the tissue species
     */
    Species& get_species(const EpigeneticGenotypeId& genotype_id);

    /**
     * @brief Add cell driver genotype
     * 
     * @param genotype is the driver genotype of the cell
     * @param position is the initial position in the tissue
     * @return a reference to the updated object
     */
    Tissue& add(const EpigeneticGenotypeId genotype, const PositionInTissue position);

    /**
     * @brief Add a new species to the tissue
     * 
     * @param genotype is the driver genotype of the new species
     * @return a reference to the updated object
     */
    Tissue& add_species(const SomaticGenotype& genotype);

    /**
     * @brief Test whether a position is valid in a tissue
     * 
     * @param position is the position to be tested
     * @return true if and only if `position` is a valid 
     *      position for the tissue
     */
    bool is_valid(const PositionInTissue& position) const;

    /**
     * @brief Get the number of species in the tissue
     * 
     * @return the number of tissue species
     */
    size_t num_of_species() const;

    /**
     * @brief Get the number of driver mutated cells in the tissue
     * 
     * @return the number of driver mutated cells in the tissue
     */
    size_t num_of_mutated_cells() const;

    /**
     * @brief Get an iterator over the species having the same somatic genotype
     * 
     * @param genotype_id is the identifier of the somatic genotype whose species are aimed 
     * @return an interator over the tissue's species having `genotype_id` as somatic
     *       genotype identifier
     */
    SpeciesView get_somatic_genotype_species(const SomaticGenotypeId& genotype_id) const;

    /**
     * @brief Get the cell in a position
     * 
     * @param position is the position of the aimed cell
     * @return a cell in tissue proxy
     */
    CellInTissueProxy operator()(const PositionInTissue& position);

    /**
     * @brief Get the cell in a position
     * 
     * @param position is the position of the aimed cell
     * @return a constant cell in tissue proxy
     */
    const CellInTissueConstantProxy operator()(const PositionInTissue& position) const;

    /**
     * @brief Get tissue name
     * 
     * @return a constant reference to the tissue name
     */
    const std::string& get_name() const;

    /**
     * @brief Get tissue name
     * 
     * @return a constant reference to the tissue name
     */
    const uint8_t& num_of_dimensions() const;

    /**
     * @brief Push contiguous driver mutated cells in a direction
     * 
     * This method pushes the contiguous driver mutated cells from a
     * position towards a direction and returns the list of the cells
     * that are pushed outside the tissue. 
     * 
     * @param from_position is the position from which the cells are pushed
     * @param direction is the push direction
     * @return the list of the cells that have been pushed outside the 
     *      tissue border
     */
    std::list<Cell> push_cells(const PositionInTissue from_position, const Direction& direction);

    /**
     * @brief Get the tissue size
     * 
     * @return the tissue size for the 3 dimensions
     */
    std::vector<AxisSize> size() const;

    /**
     * @brief Save a tissue in an archive
     * 
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        archive & size()
                & name
                & species
                & somatic_genotope_pos;
    }

    /**
     * @brief Load a tissue from an archive
     * 
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the loaded tissue
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    static Tissue load(ARCHIVE& archive)
    {
        std::vector<AxisSize> sizes;

        archive & sizes;

        Tissue tissue(sizes);

        archive & tissue.name
                & tissue.species
                & tissue.somatic_genotope_pos;

        size_t i=0;
        for (auto& species : tissue.species) {
            tissue.pos_map[species.get_id()] = ++i;
            for (auto& cell : species) {
                tissue.cell_pointer(cell) = const_cast<CellInTissue*>(&cell);
            }
        }

        return tissue;
    }

    template<typename TISSUE_TYPE>
    friend class BaseCellInTissueProxy;
};

/**
 * @brief Test whether two iterators differs
 * 
 * @param a is the first iterator to compare
 * @param b is the second iterator to compare
 * @return `true` if and only if the two iterators 
 *      do not refer to the same object
 */
inline bool operator!=(const Tissue::SpeciesView::const_iterator& a, const Tissue::SpeciesView::const_iterator& b)
{ 
    return !(a==b); 
}

/* Inline implementation */

inline std::vector<Species>::const_iterator Tissue::begin() const
{
    return std::begin(species);
}

inline std::vector<Species>::const_iterator Tissue::end() const
{
    return std::end(species);
}

inline Tissue::SpeciesView Tissue::get_somatic_genotype_species(const SomaticGenotypeId& genotype_id) const
{
    return Tissue::SpeciesView(species, somatic_genotope_pos.at(genotype_id));
}

inline const Species& Tissue::get_species(const EpigeneticGenotypeId& genotype) const
{
    return species[pos_map.at(genotype)];
}

inline Species& Tissue::get_species(const EpigeneticGenotypeId& genotype)
{
    return species[pos_map.at(genotype)];
}

inline size_t Tissue::num_of_species() const
{
    return species.size();
}

inline bool Tissue::is_valid(const PositionInTissue& position) const
{
    return (static_cast<size_t>(position.x)<space.size() && 
            static_cast<size_t>(position.y)<space[0].size() && 
            static_cast<size_t>(position.z)<space[0][0].size() &&
            position.x>=0 && position.y>=0 && position.z>=0);
}

inline const std::string& Tissue::get_name() const
{
    return name;
}

inline const uint8_t& Tissue::num_of_dimensions() const
{
    return dimensions;
}

inline CellInTissue*& Tissue::cell_pointer(const PositionInTissue& position)
{
    return this->space[position.x][position.y][position.z];
}


inline const CellInTissue* Tissue::cell_pointer(const PositionInTissue& position) const
{
    return this->space[position.x][position.y][position.z];
}

template<typename TISSUE_TYPE>
Tissue::BaseCellInTissueProxy<TISSUE_TYPE>::BaseCellInTissueProxy(TISSUE_TYPE &tissue, const PositionInTissue position):
    tissue(tissue), position(position)
{
    if (!tissue.is_valid(position)) {
        throw std::out_of_range("The position does not belong to the tissue");
    }
}

template<typename TISSUE_TYPE>
inline bool Tissue::BaseCellInTissueProxy<TISSUE_TYPE>::has_driver_mutations() const
{
    return tissue.cell_pointer(position)!=nullptr;
}

template<typename TISSUE_TYPE>
Tissue::BaseCellInTissueProxy<TISSUE_TYPE>::operator const CellInTissue&() const
{
    const auto ptr = tissue.cell_pointer(position);

    if (ptr!=nullptr) {
        return *ptr;
    }

    throw std::runtime_error("Non-driver mutated cell");
}

};
#endif // __RACES_TISSUE__