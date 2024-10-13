/**
 * @file utils.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines utility functions
 * @version 1.1
 * @date 2024-10-13
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

#include <string>
#include <filesystem>
#include <map>
#include <set>
#include <algorithm>
#include <ostream>

#ifndef __RACES_UTILS__
#define __RACES_UTILS__

/**
 * @brief Convert a path to a string
 * 
 * @param fs_path is the path to be converted
 * @return a string representation of the path
 */
inline std::string to_string(const std::filesystem::path& fs_path)
{
#if defined(__WIN32__) || defined(__WIN64__)
    std::wstring w_path(fs_path);

    return std::string(w_path.begin(), w_path.end());
#else
    return fs_path;
#endif
}

/**
 * @brief Compute the ceiling of a division
 *
 * @tparam TYPE is the type of the operands
 * @param x is the dividend
 * @param y is the divisor
 * @return the ceiling of the division `x/y`
 */
template<typename TYPE>
inline static TYPE ceil_div(const TYPE& x, const TYPE& y)
{
    return 1 + ((x - 1) / y);
}

namespace std
{

/**
 * @brief Test whether a string is a suffix of another string
 * 
 * @param suffix is the suffix to be tested
 * @param text is the test whose suffix is tested
 * @return `true` if and only if `suffix` is a suffix of `text` 
 */
inline bool is_suffix_of(const std::string& suffix, const std::string& text) {
    return ((text.size()>=suffix.size()) 
            && (text.compare(text.size()-suffix.size(), suffix.size(), suffix)));
}

/**
 * @brief Get the union of two sets
 * 
 * @tparam T is the type of the set values
 * @param A is a set of T values
 * @param B is a set of T values
 * @return the set \$A \cup B\$
 */
template<typename T>
std::set<T> get_union(const std::set<T>& A, const std::set<T>& B)
{
    std::set<T> C;

    std::set_union(A.begin(), A.end(), B.begin(), B.end(),
                   std::inserter(C, C.end()));

    return C;
}

/**
 * @brief Stream a map
 * 
 * @tparam KEYS is the type of the map keys
 * @tparam VALUES is the type of the map values
 * @param os is the output stream
 * @param m is the map to be streamed
 * @return a reference to the updated stream
 */
template<typename KEYS, typename VALUES>
std::ostream& operator<<(std::ostream& os, const std::map<KEYS, VALUES>& m)
{
  os << "{";
  std::string sep = "";
  for (const auto& [key, value]: m) {
    os << sep << "\"" << key << "\": " << value;
    sep = ", ";
  }
  os << "}";

  return os;
}

}   // std

#endif // __RACES_UTILS__
