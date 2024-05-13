/**
 * @file sbs_signature.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines SBS signature
 * @version 0.16
 * @date 2024-05-13
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

#ifndef __RACES_SBS_SIGNATURE__
#define __RACES_SBS_SIGNATURE__

#include <string>
#include <functional> // std::less
#include <iostream>
#include <sstream>

#include "sbs_context.hpp"
#include "genomic_sequence.hpp"

#include "signature.hpp"

namespace Races
{

namespace Mutations
{

/**
 * @brief A class to represent SBS type
 *
 * A SBS type is a SBS context and the
 * replacing nucleic base for the central nucleotide.
 * For consistency with the literature, the SBS
 * context will always be either a 'C' or a 'T'. We
 * call it normal context. If the provided context
 * is not a normal context, the complement context
 * and replacing nucleic base are used.
 */
class SBSType
{
    SBSContext context;    //!< the normal context
    char replace_base;     //!< the replace base
public:
    /**
     * @brief The empty constructor
     */
    SBSType();

    /**
     * @brief A constructor
     *
     * @param context is the SBS type context
     * @param replace_base is the base replacing the context central nucleotide
     */
    SBSType(const SBSContext& context, const char& replace_base);

    /**
     * @brief A constructor
     *
     * @param context is the SBS type context
     * @param replace_base is the base replacing the context central nucleotide
     */
    SBSType(const std::string& context, const char& replace_base);

    /**
     * @brief A constructor
     *
     * A SBS type is conventionally represented by a string in the
     * form `X[Y>W]K` where `X` and `K` and the bases on 5' and 3', respectively,
     * and `Y` and `W` are the central nucleotide before and after the mutation.
     * This constructor takes as parameter a string in this format and build
     * the corresponding `SBSType` object.
     *
     * @param type is the textual representation of a SBS type
     */
    explicit SBSType(const std::string& type);

    /**
     * @brief Get the mutational context
     *
     * @return the mutational context
     */
    inline const SBSContext& get_context() const
    {
        return context;
    }

    /**
     * @brief Get the replace base
     *
     * @return the replace base
     */
    inline const char& get_replace_base() const
    {
        return replace_base;
    }

    /**
     * @brief Get complement base of the replace base
     *
     * @return the replace base
     */
    inline char get_complement_replace_base() const
    {
        return GenomicSequence::get_complemented(replace_base);
    }

    /**
     * @brief Test whether two SBS types are equivalent
     *
     * @param type is the SBS type to compare
     * @return `true` if and only if the two SBS types refer
     *      to equivalent contexts and the same replace base
     */
    inline bool operator==(const SBSType& type) const
    {
        return (context == type.context) && (replace_base == type.replace_base);
    }

    /**
     * @brief Test whether two SBS types differ
     *
     * @param type is the SBS type to compare
     * @return `true` if and only if the two SBS types refer
     *      to different contexts or different replace bases
     */
    inline bool operator!=(const SBSType& type) const
    {
        return !(*this == type);
    }
};

}   // Mutations

}   // Races


namespace std
{

template<>
struct less<Races::Mutations::SBSType>
{
    bool operator()(const Races::Mutations::SBSType &lhs,
                    const Races::Mutations::SBSType &rhs) const;
};

/**
 * @brief Stream the SBS type in a stream
 *
 * @param out is the output stream
 * @param type is the SBS type to stream
 * @return a reference to the output stream
 */
std::ostream& operator<<(std::ostream& out, const Races::Mutations::SBSType& type);

/**
 * @brief Stream the SBS type from a stream
 *
 * @param in is the input stream
 * @param type is the object where the streamed SBS type will be placed
 * @return a reference to the input stream
 */
std::istream& operator>>(std::istream& in, Races::Mutations::SBSType& type);

}  // std

namespace Races
{

namespace Mutations
{

using SBSSignature = Signature<SBSType>;

}   // Mutations

}   // Races


#endif // __RACES_SBS_SIGNATURE__