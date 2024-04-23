/**
 * @file sid.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines SNV, Insertion, and Deletion mutations
 * @version 0.1
 * @date 2024-04-23
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

#ifndef __RACES_SID__
#define __RACES_SID__

#include <ostream>

#include "mutation.hpp"

namespace Races
{

namespace Mutations
{

/**
 * @brief A class to represent SNVs, Insertions and Deletions
 */
struct SID : public Mutation
{
    std::string ref;      //!< The reference sequence
    std::string alt;      //!< The altered sequence

    /**
     * @brief The empty constructor
     */
    SID();

    /**
     * @brief A constructor
     *
     * @param chromosome_id is the id of the chromosome in which the mutation occurs
     * @param chromosomic_position is the position in the chromosome of the mutation
     * @param ref is the reference sequence
     * @param alt is the altered sequence
     * @param nature is the nature of the mutation
     * @throw std::domain_error reference and altered sequence are the same
     */
    SID(const ChromosomeId& chromosome_id, const ChrPosition& chromosomic_position,
        const char ref_base, const char alt_base,
        const Mutation::Nature& nature=Mutation::UNDEFINED);

    /**
     * @brief A constructor
     *
     * @param chromosome_id is the id of the chromosome in which the mutation occurs
     * @param chromosomic_position is the position in the chromosome of the mutation
     * @param ref is the reference sequence
     * @param alt is the altered sequence
     * @param cause is the cause of the mutation
     * @param nature is the nature of the mutation
     * @throw std::domain_error reference and altered sequence are the same
     */
    SID(const ChromosomeId& chromosome_id, const ChrPosition& chromosomic_position,
        const char ref_base, const char alt_base, const std::string& cause,
        const Mutation::Nature& nature=Mutation::UNDEFINED);

    /**
     * @brief A constructor
     *
     * @param chromosome_id is the id of the chromosome in which the mutation occurs
     * @param chromosomic_position is the position in the chromosome of the mutation
     * @param ref is the reference sequence
     * @param alt is the altered sequence
     * @param nature is the nature of the mutation
     * @throw std::domain_error reference and altered sequence are the same
     */
    SID(const ChromosomeId& chromosome_id, const ChrPosition& chromosomic_position,
        const std::string& ref, const std::string& alt,
        const Mutation::Nature& nature=Mutation::UNDEFINED);

    /**
     * @brief A constructor
     *
     * @param genomic_position is the mutation genomic position
     * @param ref is the reference sequence
     * @param alt is the altered sequence
     * @param nature is the nature of the mutation
     * @throw std::domain_error reference and altered sequence are the same
     */
    SID(const GenomicPosition& genomic_position,
        const std::string& ref, const std::string& alt,
        const Mutation::Nature& nature=Mutation::UNDEFINED);

    /**
     * @brief A constructor
     *
     * @param genomic_position is the mutation genomic position
     * @param ref is the reference sequence
     * @param alt is the altered sequence
     * @param nature is the nature of the mutation
     * @throw std::domain_error reference and altered sequence are the same
     */
    SID(GenomicPosition&& genomic_position,
        const std::string& ref, const std::string& alt,
        const Mutation::Nature& nature=Mutation::UNDEFINED);


    /**
     * @brief A constructor
     *
     * @param chromosome_id is the id of the chromosome in which the mutation occurs
     * @param chromosomic_position is the position in the chromosome of the mutation
     * @param ref is the reference sequence
     * @param alt is the altered sequence
     * @param cause is the cause of the mutation
     * @param nature is the nature of the mutation
     * @throw std::domain_error reference and altered sequence are the same
     */
    SID(const ChromosomeId& chromosome_id, const ChrPosition& chromosomic_position,
        const std::string& ref, const std::string& alt,
        const std::string& cause,
        const Mutation::Nature& nature=Mutation::UNDEFINED);

    /**
     * @brief A constructor
     *
     * @param genomic_position is the mutation genomic position
     * @param ref is the reference sequence
     * @param alt is the altered sequence
     * @param cause is the cause of the mutation
     * @param nature is the nature of the mutation
     * @throw std::domain_error reference and altered sequence are the same
     */
    SID(const GenomicPosition& genomic_position,
        const std::string& ref, const std::string& alt,
        const std::string& cause,
        const Mutation::Nature& nature=Mutation::UNDEFINED);

    /**
     * @brief A constructor
     *
     * @param genomic_position is the mutation genomic position
     * @param ref is the reference sequence
     * @param alt is the altered sequence
     * @param cause is the cause of the mutation
     * @param nature is the nature of the mutation
     * @throw std::domain_error reference and altered sequence are the same
     */
    SID(GenomicPosition&& genomic_position,
        const std::string& ref, const std::string& alt,
        const std::string& cause,
        const Mutation::Nature& nature=Mutation::UNDEFINED);

    /**
     * @brief Test whether the mutation is a SNV
     * 
     * @return `true` if and only if both the reference and the alternative
     *      sequences have size 1
     */
    inline bool is_SNV() const
    {
        return ref.size()==1 && alt.size()==1;
    }

    /**
     * @brief Save an SID in an archive
     *
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        archive & static_cast<const Mutation&>(*this)
                & ref
                & alt;
    }

    /**
     * @brief Load an SID from an archive
     *
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the loaded SID mutation
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    inline static SID load(ARCHIVE& archive)
    {
        SID sid;

        archive & static_cast<Mutation&>(sid)
                & sid.ref
                & sid.alt;

        return sid;
    }
};

}   // Mutations

}   // Races

namespace std
{

template<>
struct less<Races::Mutations::SID>
{
    bool operator()(const Races::Mutations::SID &lhs,
                    const Races::Mutations::SID &rhs) const;
};

/**
 * @brief Write an SID in a output stream
 *
 * @param out is the output stream
 * @param sid is the SID to stream
 * @return a reference of the updated stream
 */
std::ostream& operator<<(std::ostream& out, const Races::Mutations::SID& sid);

}   // std


#endif // __RACES_SID__