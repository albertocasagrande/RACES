/**
 * @file snv_signature.hpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Define Single Variation Mutation mutational signature
 * @version 0.4
 * @date 2023-07-21
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

#ifndef __RACES_SNV_DISTRIBUTION__
#define __RACES_SNV_DISTRIBUTION__

#include <string>
#include <cstdint>
#include <map>
#include <set>
#include <functional> // std::less
#include <iostream>
#include <sstream>


namespace Races
{

namespace Passengers
{

namespace SNV
{

/**
 * @brief A class to represent mutational context
 * 
 * A context is a nucleic triplet where a mutation may 
 * occur.
 * Every context univocally corresponds to a code that 
 * represents the triplet. This class guarantees the 
 * mutational context code to be in the interval 
 * natural [0,63]: 5', 3', and the central nucleotide 
 * have 4 possible values, i.e., A, C, G, and T.
 * Moreover, all the codes for contexts having either 
 * a 'C' and a 'T' as central nucleotide are guaranteed
 * to be in the interval [0,31].
 */
class MutationalContext
{
    uint8_t code;   //!< the code associated to the context
public:
    /**
     * @brief The empty constructor
     */
    MutationalContext();

    /**
     * @brief A constructor
     * 
     * @param nucleic_triplet is the string of the nucleic triplet
     */
    MutationalContext(const std::string& nucleic_triplet);

    /**
     * @brief A constructor
     * 
     * @param code is the context code
     */
    MutationalContext(const uint8_t code);

    /**
     * @brief Get the mutational context code 
     * 
     * @return a constant reference to the mutational context code 
     */
    inline const uint8_t& get_code() const
    {
        return code;
    }

    /**
     * @brief Get the nucleic sequence of the context
     * 
     * @return the nucleic sequence of the context
     */
    std::string get_sequence() const;

    /**
     * @brief Get the context central nucleotide
     * 
     * @return the context central nucleotide
     */
    char get_central_nucleotide() const;

    /**
     * @brief Get the complement mutational context
     * 
     * @return the complement mutational context 
     */
    inline MutationalContext get_complement() const
    {
        return MutationalContext(get_complement(code));
    }

    /**
     * @brief Test whether two mutational contexts are equivalent
     * 
     * @param context is the mutational context to compare
     * @return `true` if and only if the two mutational contexts represent
     *      the same nucleic triplet
     */
    inline bool operator==(const MutationalContext& context) const
    {
        return get_code() == context.get_code();
    }

    /**
     * @brief Test whether two mutational contexts differ
     * 
     * @param context is the mutational context to compare
     * @return `true` if and only if the two mutational contexts represent
     *      different nucleic triplets
     */
    inline bool operator!=(const MutationalContext& context) const
    {
        return get_code() != context.get_code();
    }

    /**
     * @brief Get the code of the complement mutational context
     * 
     * @param code is the mutational context code whose complement is request 
     * @return the code of the complement context of the 
     *      context whose code is `code` 
     */
    static uint8_t get_complement(const uint8_t& code);

    /**
     * @brief Get the complement base
     * 
     * @param base is the base whose complement is request 
     * @return the complement base of `base` 
     */
    static char get_complement(const char& base);

    /**
     * @brief Get the complement sequence
     * 
     * @param sequence is the sequence whose complement is request 
     * @return the complement sequence of `sequence` 
     */
    static std::string get_complement(const std::string& sequence);

    /**
     * @brief Check whether a symbol represents a nucleic base
     * 
     * @param symbol is the symbol to be tested 
     * @return `true` if and only if symbol represents a base, i.e., 
     *       it is one among 'A', 'C', 'G', 'T', 'a', 'c', 'g', or 't'.
     */
    static bool is_a_base(const char& symbol);
};

}   // SNV

}   // Passengers

}   // Races


template<>
struct std::less<Races::Passengers::SNV::MutationalContext>
{
    inline constexpr bool operator()(const Races::Passengers::SNV::MutationalContext &lhs,
                                     const Races::Passengers::SNV::MutationalContext &rhs) const
    {
        return lhs.get_code() < rhs.get_code();
    }
};

namespace Races 
{

namespace Passengers
{

namespace SNV
{
/**
 * @brief A class to represent mutational type
 * 
 * A mutational type is a mutational context and the 
 * replacing nucleic base for the central nucleotide.
 * For consistency with the literature, the mutational 
 * context will always be either a 'C' or a 'T'. We 
 * call it normal context. If the provided context 
 * is not a normal context, the complement context 
 * and replacing nucleic base are used.
 */
class MutationalType
{
    MutationalContext context;    //!< the normal context
    char replace_base;            //!< the replace base
public:
    /**
     * @brief The empty constructor
     */
    MutationalType();

    /**
     * @brief A constructor
     * 
     * @param context is the mutational type context
     * @param substituting_base is the base replacing the context central nucleotide
     */
    MutationalType(const MutationalContext& context, const char& replace_base);

    /**
     * @brief A constructor
     * 
     * @param context is the mutational type context
     * @param substituting_base is the base replacing the context central nucleotide
     */
    MutationalType(const std::string& context, const char& replace_base);

    /**
     * @brief A constructor
     * 
     * A mutational type is conventionally represented by a string in the 
     * form `X[Y>W]K` where `X` and `K` and the bases on 5' and 3', respectively,  
     * and `Y` and `W` are the central nucleotide before and after the mutation.
     * This constructor takes as parameter a string in this format and build 
     * the corresponding `MutationalType` object.
     * 
     * @param type is the textual representation of a mutational type
     */
    MutationalType(const std::string& type);

    /**
     * @brief Get the mutational context
     * 
     * @return the mutational context
     */
    inline const MutationalContext& get_context() const
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
        return MutationalContext::get_complement(replace_base);
    }

    /**
     * @brief Test whether two mutational types are equivalent
     * 
     * @param type is the mutational type to compare
     * @return `true` if and only if the two mutational types refer
     *      to equivalent contexts and the same replace base
     */
    inline bool operator==(const MutationalType& type) const
    {
        return (context == type.context) && (replace_base == type.replace_base);
    }

    /**
     * @brief Test whether two mutational types differ
     * 
     * @param type is the mutational type to compare
     * @return `true` if and only if the two mutational types refer
     *      to different contexts or different replace bases
     */
    inline bool operator!=(const MutationalType& type) const
    {
        return !(*this == type);
    }
};

}   // SNV

}   // Passengers

}   // Races

template<>
struct std::less<Races::Passengers::SNV::MutationalType>
{
    bool operator()(const Races::Passengers::SNV::MutationalType &lhs,
                    const Races::Passengers::SNV::MutationalType &rhs) const
    {
        const auto& lhs_code = lhs.get_context().get_code();
        const auto& rhs_code = rhs.get_context().get_code();

        return ((lhs_code < rhs_code) || 
                ((lhs_code == rhs_code) && (lhs.get_replace_base()<rhs.get_replace_base())));
    }
};

namespace Races 
{

namespace Passengers
{

namespace SNV
{

class MutationalSignature;

/**
 * @brief A class to represent the result of a mutational signature expression
 * 
 * This class is meant to represent temporary object that are evaluated during 
 * the computation of expressions of the kind:
 *  
 * \f$\alpha_1 * \beta_1 + alpha_2 * beta_2 + \ldots\f$
 * 
 * where the \f$\alpha_i\f$'s are real values in the interval \f$[0,1]\f$ and 
 * the \f$\beta_i\f$'s are `MutationalSignature` objects. 
 * 
 * Even if the final result of above expression is a mutational signature, the
 * partial results may be different from a probability distribution and that 
 * is why this class is needed.
 */
class MutationalSignatureExprResult
{
    std::map<MutationalType, double> value_map; //!< the mutational type-value map

    /**
     * @brief The constructor
     * 
     * This constructor is private and it is meant to be exclusively called by 
     * `MutationalSignature`'s methods.
     * 
     * @param value_map is a mutational type-value map
     */
    MutationalSignatureExprResult(const std::map<MutationalType, double>& value_map);
public:
    /**
     * @brief The empty constructor
     */
    MutationalSignatureExprResult();

    /**
     * @brief Cast to `MutationalSignature`
     * 
     * This method tries to cast a mutational signature expression to a 
     * mutational signature. When the expression does not represent a 
     * probability distribution a `std::domain_error` is thrown.
     * 
     * @return the corresponding `MutationalSignature` object
     */
    operator MutationalSignature();

    /**
     * @brief Inplace multiply by an arithmetic value
     * 
     * @tparam T is the type of the multiplicand
     * @param value is the multiplicand
     * @return a reference to the updated object
     */
    template<typename T, std::enable_if_t<std::is_arithmetic_v<T>, bool> = true>
    MutationalSignatureExprResult& operator*(const T& value)
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
     * @brief Inplace add a mutational signature expression value
     * 
     * @param expression_value is a mutational signature expression value
     * @return a reference to the updated object
     */
    MutationalSignatureExprResult& operator+(MutationalSignatureExprResult&& expression_value);

    /**
     * @brief Inplace add a signature
     * 
     * @param signature is a mutational signature
     * @return a reference to the updated object
     */
    MutationalSignatureExprResult& operator+(const MutationalSignature& signature);

    friend class MutationalSignature;
};

/**
 * @brief A class to represent a mutational signature
 * 
 * A mutational signature is a probability distribution on 
 * the set of mutational types.
 */
class MutationalSignature
{
    std::map<MutationalType, double> dist_map; //!< the signature probability distribution map
public:
    using const_iterator = std::map<MutationalType, double>::const_iterator;

    /**
     * @brief The empty constructor
     */
    MutationalSignature();

    /**
     * @brief A constructor
     * 
     * @param distribution is a mutation type-value map representing a distribution
     */
    MutationalSignature(const std::map<MutationalType, double>& distribution);

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
     * @brief Get the probability associated to a mutational type
     * 
     * @param type is the mutational type whose probability is aimed
     * @return the probability of `type`
     */
    inline double operator()(const MutationalType& type) const
    {
        auto it = dist_map.find(type);
        if (it != dist_map.end()) {
            return it->second;
        } 
    
        return 0;
    }

    /**
     * @brief Multiply by an arithmetic value
     * 
     * @tparam T is the type of the multiplicand
     * @param value is the multiplicand
     * @return the resulting mutational signature expression value
     */
    template<typename T, std::enable_if_t<std::is_arithmetic_v<T>, bool> = true>
    inline MutationalSignatureExprResult operator*(const T& value) const
    {
        return MutationalSignatureExprResult(dist_map) * value;
    }

    /**
     * @brief Add a signature
     * 
     * @param signature is a mutational signature
     * @return the resulting mutational signature expression value
     */
    inline MutationalSignatureExprResult operator+(const MutationalSignature& signature) const
    {
        return MutationalSignatureExprResult(dist_map) + signature;
    }

    /**
     * @brief Read mutational signature from a input stream
     * 
     * @param in is the input stream
     * @return a map that associates the name of the signatures in the file 
     *         and the corresponding signature.
     */
    static std::map<std::string, MutationalSignature> read_from_stream(std::istream& in);

    /**
     * @brief Read mutational signature from a input stream
     * 
     * @param in is the input stream
     * @param signature_names is the set of the requested signature
     * @return a map that associates the name of the signatures in the file that
     *         match `signature_names` and the corresponding signature.
     */
    static std::map<std::string, MutationalSignature> read_from_stream(std::istream& in, const std::set<std::string>& signature_names);
};

/**
 * @brief Multiply an arithmetic value and a mutational signature 
 * 
 * @tparam T is the type of the arithmetic value
 * @param value is the arithmetic value
 * @param signature 
 * @return a `MutationalSignatureExprResult` object representing the
 *         the multiplication result 
 */
template<typename T, std::enable_if_t<std::is_arithmetic_v<T>, bool> = true>
inline MutationalSignatureExprResult operator*(const T& value, const MutationalSignature& signature)
{
    return signature * value;
}

}   // SNV

}   // Passengers

}   // Races


namespace std 
{
/**
 * @brief Stream the mutational context in a stream
 * 
 * @param out is the output stream
 * @param context is the mutational context to stream
 * @return a reference to the output stream
 */
inline std::ostream& operator<<(std::ostream& out, const Races::Passengers::SNV::MutationalContext& context)
{
    return (out << context.get_sequence());
}

/**
 * @brief Stream the mutational type in a stream
 * 
 * @param out is the output stream
 * @param type is the mutational type to stream
 * @return a reference to the output stream
 */
std::ostream& operator<<(std::ostream& out, const Races::Passengers::SNV::MutationalType& type);

/**
 * @brief Stream the mutational type from a stream
 * 
 * @param in is the input stream
 * @param type is the object where the streamed mutational type will be placed
 * @return a reference to the input stream
 */
std::istream& operator>>(std::istream& in, Races::Passengers::SNV::MutationalType& type);

}  // std

#endif // __RACES_SNV_DISTRIBUTION__