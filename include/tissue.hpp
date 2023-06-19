/**
 * @file tissue.hpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Define tissue class
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

#ifndef __RACES_TISSUE__
#define __RACES_TISSUE__

#include <vector>
#include <map>
#include <string>

#include "time.hpp"
#include "species.hpp"
#include "cell.hpp"
#include "mutation_graphs.hpp"

namespace Races {

class Tissue {
    std::string name;                               //!< tissue name
    std::vector<Species> species;                   //!< species in the tissue
    std::map<DriverGenotypeId, size_t> pos_map;     //!< identifier to position map

    DriverSomaticGraph somatic_graph;               //!< somatic mutation graph
    DriverEpigeneticGraph epigenetic_graph;         //!< epigenetic mutation graph
    
    std::vector<std::vector<std::vector<CellInTissue>>> space;     //!< Space in the tissue

public:

    /**
     * @brief A constructor
     * 
     * @param x_size is the size of the tissue on the x axis
     * @param y_size is the size of the tissue on the y axis
     * @param z_size is the size of the tissue on the z axis
     */
    Tissue(const unsigned int x_size, const unsigned int  y_size, const unsigned int  z_size=1);

    /**
     * @brief A constructor
     * 
     * @param genotypes is the vector of driver genotypes
     * @param x_size is the size of the tissue on the x axis
     * @param y_size is the size of the tissue on the y axis
     * @param z_size is the size of the tissue on the z axis
     */
    Tissue(const std::vector<DriverGenotype> genotypes, const unsigned int  x_size, const unsigned int  y_size, const unsigned int z_size=1);

    /**
     * @brief A constructor
     * 
     * @param name is the tissue name
     * @param x_size is the size of the tissue on the x axis
     * @param y_size is the size of the tissue on the y axis
     * @param z_size is the size of the tissue on the z axis
     */
    Tissue(const std::string name, const unsigned int x_size, const unsigned int  y_size, const unsigned int  z_size=1);

    /**
     * @brief A constructor
     * 
     * @param name is the tissue name
     * @param genotypes is the vector of driver genotypes
     * @param x_size is the size of the tissue on the x axis
     * @param y_size is the size of the tissue on the y axis
     * @param z_size is the size of the tissue on the z axis
     */
    Tissue(const std::string name, const std::vector<DriverGenotype> genotypes, const unsigned int  x_size, const unsigned int  y_size, const unsigned int z_size=1);

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
     * @return a non-constant reference to the tissue species
     */
    Species& get_species(const DriverGenotypeId& genotype_id);

    /**
     * @brief Add cell driver genotype
     * 
     * @param genotype is the driver genotype of the cell
     * @param position is the initial position in the tissue
     * @param time is the insertion time
     * @return a reference to the updated object
     */
    Tissue& add(const DriverGenotypeId genotype, const PositionInTissue position, const Time time);

    /**
     * @brief Add a new species to the tissue
     * 
     * @param genotype is the driver genotype of the new species
     * @return a reference to the updated object
     */
    Tissue& add_species(const DriverGenotype& genotype);

    /**
     * @brief Add a delayed somatic mutation between driver genotypes
     * 
     * @param src is the original driver genotype
     * @param dst is the final driver genotype
     * @param delay is the time delay for the occurrence of the mutation
     * @return the updated tissue
     */
    Tissue& add_driver_somatic_mutation(const DriverGenotypeId& src, const DriverGenotypeId& dst, const Time delay);

    /**
     * @brief Add a delayed epigenetic mutation between driver genotypes
     * 
     * @param src is the original driver genotype
     * @param dst is the final driver genotype
     * @param probability is the methylation/demethylation probability 
     * @return the updated tissue
     */
    Tissue& add_driver_epigenetic_mutation(const DriverGenotypeId& src, const DriverGenotypeId& dst, const double probability);

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
     * @brief Get the position of the first non-driver cell in a direction
     * 
     * @param position is a position
     * @param direction is the direction along which non-driver cells are searched
     * @return get the position of the first non-driver cell in a direction when 
     *     available. When such a cell is not available, a non valid position.
     */
    PositionInTissue get_non_driver_in_direction(PositionInTissue position, const Direction& direction) const;

    /**
     * @brief Push a cell in a direction
     * 
     * @param position is the position of the cell to push
     * @param direction is the direction of the push
     * @return a reference to the updated object
     */
    Tissue& push_cells(PositionInTissue position, const Direction& direction);

    /**
     * @brief Get the cell in a position
     * 
     * @param position is the position of the aimed cell
     * @return a non-constant reference to the aimed cell
     */
    CellInTissue& operator()(const PositionInTissue& position);

    /**
     * @brief Get the cell in a position
     * 
     * @param position is the position of the aimed cell
     * @return a constant reference to the aimed cell
     */
    const CellInTissue& operator()(const PositionInTissue& position) const;

    /**
     * @brief Get tissue name
     * 
     * @return a constant reference to the tissue name
     */
    const std::string& get_name() const;

    /**
     * @brief Get the tissue somatic graph
     * 
     * @return a constant reference to the tissue somatic graph
     */
    const DriverSomaticGraph& get_somatic_graph() const;

    /**
     * @brief Get the tissue epigenetic graph
     * 
     * @return a constant reference to the tissue epigenetic graph
     */
    const DriverEpigeneticGraph& get_epigenetic_graph() const;

    /**
     * @brief Get the tissue size
     * 
     * @return the tissue size for the 3 dimensions
     */
    std::vector<size_t> size() const;
};

/* Inline implementation */

inline std::vector<Species>::const_iterator Tissue::begin() const
{
    return std::begin(species);
}

inline std::vector<Species>::const_iterator Tissue::end() const
{
    return std::end(species);
}

inline Species& Tissue::get_species(const DriverGenotypeId& genotype)
{
    return species[pos_map.at(genotype)];
}

inline size_t Tissue::num_of_species() const
{
    return species.size();
}

inline bool Tissue::is_valid(const PositionInTissue& position) const
{
    return position.x<space.size() && position.y<space[0].size() && position.z<space[0][0].size();
}

inline std::vector<size_t> Tissue::size() const
{
    return std::vector<size_t>{space.size(), space[0].size(), space[0][0].size()};
}

inline const std::string& Tissue::get_name() const
{
    return name;
}

inline const DriverSomaticGraph& Tissue::get_somatic_graph() const
{
    return somatic_graph;
}

inline const DriverEpigeneticGraph& Tissue::get_epigenetic_graph() const
{
    return epigenetic_graph;
}

};
#endif // __RACES_TISSUE__