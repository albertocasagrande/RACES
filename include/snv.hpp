/**
 * @file snv.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines Single Nucleotide Variation
 * @version 0.15
 * @date 2024-03-12
 * 
 * @copyright Copyright (c) 2023-2024
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

#include "mutation.hpp"

namespace Races
{

namespace Mutations
{

/**
 * @brief A class to represent SNV
 */
struct SNV : public Mutation
{
    char ref_base;      //!< The ref_base base
    char alt_base;      //!< The altered base

    /**
     * @brief The empty constructor
     */
    SNV();

    /**
     * @brief A constructor
     *
     * @param chromosome_id is the id of the SNV chromosome
     * @param chromosomic_position is the SNV position in the chromosome
     * @param alt_base is the altered base
     * @param nature is the nature of the SNV
     * @throw std::domain_error reference and altered bases are the same
     */
    SNV(const ChromosomeId& chromosome_id, const ChrPosition& chromosomic_position,
        const char alt_base, const Mutation::Nature& nature=Mutation::UNDEFINED);

    /**
     * @brief A constructor
     *
     * @param chromosome_id is the id of the SNV chromosome
     * @param chromosomic_position is the SNV position in the chromosome
     * @param ref_base is the reference base
     * @param alt_base is the altered base
     * @param nature is the nature of the SNV
     * @throw std::domain_error reference and altered bases are the same
     */
    SNV(const ChromosomeId& chromosome_id, const ChrPosition& chromosomic_position,
        const char ref_base, const char alt_base,
        const Mutation::Nature& nature=Mutation::UNDEFINED);

    /**
     * @brief A constructor
     *
     * @param genomic_position is the SNV genomic_position
     * @param alt_base is the altered base
     * @param nature is the nature of the SNV
     * @throw std::domain_error reference and altered bases are the same
     */
    SNV(const GenomicPosition& genomic_position, const char alt_base,
        const Mutation::Nature& nature=Mutation::UNDEFINED);

    /**
     * @brief A constructor
     *
     * @param genomic_position is the SNV genomic_position
     * @param ref_base is the reference base
     * @param alt_base is the altered base
     * @param nature is the nature of the SNV
     * @throw std::domain_error reference and altered bases are the same
     */
    SNV(const GenomicPosition& genomic_position, const char ref_base,
        const char alt_base, const Mutation::Nature& nature=Mutation::UNDEFINED);

    /**
     * @brief A constructor
     *
     * @param genomic_position is the SNV genomic_position
     * @param alt_base is the altered base
     * @param nature is the nature of the SNV
     * @throw std::domain_error reference and altered bases are the same
     */
    SNV(GenomicPosition&& genomic_position, const char alt_base,
        const Mutation::Nature& nature=Mutation::UNDEFINED);

    /**
     * @brief A constructor
     *
     * @param genomic_position is the SNV genomic_position
     * @param ref_base is the reference base
     * @param alt_base is the altered base
     * @param nature is the nature of the SNV
     * @throw std::domain_error reference and altered bases are the same
     */
    SNV(GenomicPosition&& genomic_position, const char ref_base,
        const char alt_base, const Mutation::Nature& nature=Mutation::UNDEFINED);

    /**
     * @brief A constructor
     *
     * @param chromosome_id is the id of the SNV chromosome
     * @param chromosomic_position is the SNV position in the chromosome
     * @param alt_base is the altered base
     * @param cause is the cause of the SNV
     * @param nature is the nature of the SNV
     * @throw std::domain_error reference and altered bases are the same
     */
    SNV(const ChromosomeId& chromosome_id, const ChrPosition& chromosomic_position,
        const char alt_base, const std::string& cause,
        const Mutation::Nature& nature=Mutation::UNDEFINED);

    /**
     * @brief A constructor
     *
     * @param chromosome_id is the id of the SNV chromosome
     * @param chromosomic_position is the SNV position in the chromosome
     * @param ref_base is the reference base
     * @param alt_base is the altered base
     * @param cause is the cause of the SNV
     * @param nature is the nature of the SNV
     * @throw std::domain_error reference and altered bases are the same
     */
    SNV(const ChromosomeId& chromosome_id, const ChrPosition& chromosomic_position,
        const char ref_base, const char alt_base, const std::string& cause,
        const Mutation::Nature& nature=Mutation::UNDEFINED);

    /**
     * @brief A constructor
     *
     * @param genomic_position is the SNV genomic_position
     * @param alt_base is the altered base
     * @param cause is the cause of the SNV
     * @param nature is the nature of the SNV
     * @throw std::domain_error reference and altered bases are the same
     */
    SNV(const GenomicPosition& genomic_position, const char alt_base,
        const std::string& cause, const Mutation::Nature& nature=Mutation::UNDEFINED);

    /**
     * @brief A constructor
     *
     * @param genomic_position is the SNV genomic_position
     * @param ref_base is the reference base
     * @param alt_base is the altered base
     * @param cause is the cause of the SNV
     * @param nature is the nature of the SNV
     * @throw std::domain_error reference and altered bases are the same
     */
    SNV(const GenomicPosition& genomic_position, const char ref_base,
        const char alt_base, const std::string& cause,
        const Mutation::Nature& nature=Mutation::UNDEFINED);

    /**
     * @brief A constructor
     *
     * @param genomic_position is the SNV genomic_position
     * @param alt_base is the altered base
     * @param cause is the cause of the SNV
     * @param nature is the nature of the SNV
     * @throw std::domain_error reference and altered bases are the same
     */
    SNV(GenomicPosition&& genomic_position, const char alt_base, 
        const std::string& cause,
        const Mutation::Nature& nature=Mutation::UNDEFINED);

    /**
     * @brief A constructor
     *
     * @param genomic_position is the SNV genomic_position
     * @param ref_base is the reference base
     * @param alt_base is the altered base
     * @param cause is the cause of the SNV
     * @param nature is the nature of the SNV
     * @throw std::domain_error reference and altered bases are the same
     */
    SNV(GenomicPosition&& genomic_position, const char ref_base,
        const char alt_base, const std::string& cause,
        const Mutation::Nature& nature=Mutation::UNDEFINED);

    /**
     * @brief Save a SNV in an archive
     *
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        archive & static_cast<const Mutation&>(*this)
                & ref_base
                & alt_base;
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

        archive & static_cast<Mutation&>(snv)
                & snv.ref_base
                & snv.alt_base;

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