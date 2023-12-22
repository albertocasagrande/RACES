/**
 * @file snv.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines Single Nucleotide Variation
 * @version 0.10
 * @date 2023-12-22
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

#ifndef __RACES_SNV__
#define __RACES_SNV__

#include <ostream>

#include "genomic_position.hpp"
#include "context.hpp"

namespace Races
{

namespace Mutations
{

/**
 * @brief A class to represent SNV
 */
struct SNV : public GenomicPosition
{
    /**
     * @brief The type of the SNV
     *
     * This enumeration defines the types of an SNV.
     * An object of this type can either be DRIVER,
     * PASSENGER, or UNDEFINED.
     */
    enum Type {
        DRIVER,
        PASSENGER,
        UNDEFINED
    };

    MutationalContext context;  //!< The context base
    char mutated_base;          //!< The mutated base
    std::string cause;          //!< The cause of the SNV
    Type type;                  //!< The type of the SNV

    /**
     * @brief The empty constructor
     */
    SNV();

    /**
     * @brief A constructor
     *
     * @param chromosome_id is the id of the SNV chromosome
     * @param chromosomic_position is the SNV position in the chromosome
     * @param context is the SNV context
     * @param mutated_base is the mutated base
     * @param type is the type of the SNV
     * @throw std::domain_error `original_base` and `mutated_base` are the same
     */
    SNV(const ChromosomeId& chromosome_id, const ChrPosition& chromosomic_position,
        const MutationalContext& context, const char mutated_base, const Type& type=UNDEFINED);

    /**
     * @brief A constructor
     *
     * @param genomic_position is the SNV genomic_position
     * @param context is the SNV context
     * @param mutated_base is the mutated base
     * @param type is the type of the SNV
     * @throw std::domain_error `original_base` and `mutated_base` are the same
     */
    SNV(const GenomicPosition& genomic_position, const MutationalContext& context,
        const char mutated_base, const Type& type=UNDEFINED);

    /**
     * @brief A constructor
     *
     * @param genomic_position is the SNV genomic_position
     * @param context is the SNV context
     * @param mutated_base is the mutated base
     * @param type is the type of the SNV
     * @throw std::domain_error `original_base` and `mutated_base` are the same
     */
    SNV(GenomicPosition&& genomic_position, const MutationalContext& context, const char mutated_base,
        const Type& type=UNDEFINED);

    /**
     * @brief A constructor
     *
     * @param chromosome_id is the id of the SNV chromosome
     * @param chromosomic_position is the SNV position in the chromosome
     * @param context is the SNV context
     * @param mutated_base is the mutated base
     * @param cause is the cause of the SNV
     * @param type is the type of the SNV
     * @throw std::domain_error `original_base` and `mutated_base` are the same
     */
    SNV(const ChromosomeId& chromosome_id, const ChrPosition& chromosomic_position,
        const MutationalContext& context, const char mutated_base, const std::string& cause,
        const Type& type=UNDEFINED);

    /**
     * @brief A constructor
     *
     * @param genomic_position is the SNV genomic_position
     * @param context is the SNV context
     * @param mutated_base is the mutated base
     * @param cause is the cause of the SNV
     * @param type is the type of the SNV
     * @throw std::domain_error `original_base` and `mutated_base` are the same
     */
    SNV(const GenomicPosition& genomic_position, const MutationalContext& context,
        const char mutated_base, const std::string& cause, const Type& type=UNDEFINED);

    /**
     * @brief A constructor
     *
     * @param genomic_position is the SNV genomic_position
     * @param context is the SNV context
     * @param mutated_base is the mutated base
     * @param cause is the cause of the SNV
     * @param type is the type of the SNV
     * @throw std::domain_error `original_base` and `mutated_base` are the same
     */
    SNV(GenomicPosition&& genomic_position, const MutationalContext& context,
        const char mutated_base, const std::string& cause, const Type& type=UNDEFINED);

    /**
     * @brief Save a SNV in an archive
     *
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        archive & static_cast<const GenomicPosition&>(*this)
                & context
                & mutated_base
                & cause
                & type;
    }


    /**
     * @brief Load a SNV from an archive
     *
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the load SNV
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    inline static SNV load(ARCHIVE& archive)
    {
        SNV snv;

        archive & static_cast<GenomicPosition&>(snv)
                & snv.context
                & snv.mutated_base
                & snv.cause
                & snv.type;

        return snv;
    }
};

}   // Mutations

}   // Races

namespace std
{

template<>
struct less<Races::Mutations::SNV>
{
    bool operator()(const Races::Mutations::SNV &lhs,
                    const Races::Mutations::SNV &rhs) const;
};

/**
 * @brief Write a SNV in a output stream
 *
 * @param out is the output stream
 * @param snv is the SNV to stream
 * @return a reference of the updated stream
 */
std::ostream& operator<<(std::ostream& out, const Races::Mutations::SNV& snv);

}   // std


#endif // __RACES_SNV__