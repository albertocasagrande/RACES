/**
 * @file sbs_signature.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines SBS signature
 * @version 0.15
 * @date 2024-05-11
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
#include <map>
#include <set>
#include <functional> // std::less
#include <iostream>
#include <sstream>

#include "sbs_context.hpp"
#include "genomic_sequence.hpp"

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

template<>
struct std::less<Races::Mutations::SBSType>
{
    bool operator()(const Races::Mutations::SBSType &lhs,
                    const Races::Mutations::SBSType &rhs) const;
};

namespace Races
{

namespace Mutations
{


class SBSSignature;

/**
 * @brief A class to represent the result of a SBS signature expression
 *
 * This class is meant to represent temporary object that are evaluated during
 * the computation of expressions of the kind:
 *
 * \f$\alpha_1 * \beta_1 + alpha_2 * beta_2 + \ldots\f$
 *
 * where the \f$\alpha_i\f$'s are real values in the interval \f$[0,1]\f$ and
 * the \f$\beta_i\f$'s are `SBSSignature` objects.
 *
 * Even if the final result of above expression is a SBS signature, the
 * partial results may be different from a probability distribution and that
 * is why this class is needed.
 */
class SBSSignatureExprResult
{
    std::map<SBSType, double> value_map; //!< the SBS type-value map

    /**
     * @brief The constructor
     *
     * This constructor is private and it is meant to be exclusively called by
     * `SBSSignature`'s methods.
     *
     * @param value_map is a SBS type-value map
     */
    SBSSignatureExprResult(const std::map<SBSType, double>& value_map);
public:
    /**
     * @brief The empty constructor
     */
    SBSSignatureExprResult();

    /**
     * @brief Cast to `SBSSignature`
     *
     * This method tries to cast a SBS signature expression to a
     * SBS signature. When the expression does not represent a
     * probability distribution a `std::domain_error` is thrown.
     *
     * @return the corresponding `SBSSignature` object
     */
    operator SBSSignature();

    /**
     * @brief Inplace multiply by an arithmetic value
     *
     * @tparam T is the type of the multiplicand
     * @param value is the multiplicand
     * @return a reference to the updated object
     */
    template<typename T, std::enable_if_t<std::is_arithmetic_v<T>, bool> = true>
    SBSSignatureExprResult& operator*(const T& value)
    {
        if (value>1 || value<0) {
            std::ostringstream oss;

            oss << "the multiplying value must be in the real interval [0,1]. "
                << std::to_string(value) << " passed.";

            throw std::domain_error(oss.str());
        }

        for (auto& [type, type_value]: value_map) {
            type_value *= value;
        }

        return *this;
    }

    /**
     * @brief Inplace add a SBS signature expression value
     *
     * @param expression_value is a SBS signature expression value
     * @return a reference to the updated object
     */
    SBSSignatureExprResult& operator+(SBSSignatureExprResult&& expression_value);

    /**
     * @brief Inplace add a signature
     *
     * @param signature is a SBS signature
     * @return a reference to the updated object
     */
    SBSSignatureExprResult& operator+(const SBSSignature& signature);

    friend class SBSSignature;
};

/**
 * @brief A class to represent a SBS signature
 *
 * A SBS signature is a probability distribution on
 * the set of SBS types.
 */
class SBSSignature
{
    std::map<SBSType, double> dist_map; //!< the signature probability distribution map
public:
    using const_iterator = std::map<SBSType, double>::const_iterator;

    /**
     * @brief The empty constructor
     */
    SBSSignature();

    /**
     * @brief A constructor
     *
     * @param distribution is a mutation type-value map representing a distribution
     */
    explicit SBSSignature(const std::map<SBSType, double>& distribution);

    /**
     * @brief Get the initial constant iterator
     *
     * @return the initial constant iterator
     */
    inline const_iterator begin() const
    {
        return dist_map.begin();
    }

    /**
     * @brief Get the final constant iterator
     *
     * @return the final constant iterator
     */
    inline const_iterator end() const
    {
        return dist_map.end();
    }

    /**
     * @brief Get the probability associated to a SBS type
     *
     * @param type is the SBS type whose probability is aimed
     * @return the probability of `type`
     */
    double operator()(const SBSType& type) const;

    /**
     * @brief Multiply by an arithmetic value
     *
     * @tparam T is the type of the multiplicand
     * @param value is the multiplicand
     * @return the resulting SBS signature expression value
     */
    template<typename T, std::enable_if_t<std::is_arithmetic_v<T>, bool> = true>
    inline SBSSignatureExprResult operator*(const T& value) const
    {
        return SBSSignatureExprResult(dist_map) * value;
    }

    /**
     * @brief Add a signature
     *
     * @param signature is a SBS signature
     * @return the resulting SBS signature expression value
     */
    inline SBSSignatureExprResult operator+(const SBSSignature& signature) const
    {
        return SBSSignatureExprResult(dist_map) + signature;
    }

    /**
     * @brief Read SBS signature from a input stream
     *
     * @param in is the input stream
     * @return a map that associates the name of the signatures in the file
     *         and the corresponding signature.
     */
    static std::map<std::string, SBSSignature> read_from_stream(std::istream& in);

    /**
     * @brief Read SBS signature from a input stream
     *
     * @param in is the input stream
     * @param signature_names is the set of the requested signature
     * @return a map that associates the name of the signatures in the file that
     *         match `signature_names` and the corresponding signature.
     */
    static std::map<std::string, SBSSignature> read_from_stream(std::istream& in, const std::set<std::string>& signature_names);
};

/**
 * @brief Multiply an arithmetic value and a SBS signature
 *
 * @tparam T is the type of the arithmetic value
 * @param value is the arithmetic value
 * @param signature
 * @return a `SBSSignatureExprResult` object representing the
 *         the multiplication result
 */
template<typename T, std::enable_if_t<std::is_arithmetic_v<T>, bool> = true>
inline SBSSignatureExprResult operator*(const T& value, const SBSSignature& signature)
{
    return signature * value;
}

}   // Mutations

}   // Races


namespace std
{
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

#endif // __RACES_SBS_SIGNATURE__