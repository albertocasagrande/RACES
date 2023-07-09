/**
 * @file archive.hpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Define some archive classes and their methods
 * @version 0.1
 * @date 2023-07-09
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

#ifndef __RACES_ARCHIVE__
#define __RACES_ARCHIVE__

#include <map>
#include <vector>
#include <queue>
#include <fstream>
#include <filesystem>
#include <type_traits>

namespace Races {

/**
 * @brief A structure to test the presence of a save method
 * 
 * This structure is meant to detect whether a class is 
 * equipped with a save method. The idea at the basis of 
 * this code was taken from https://stackoverflow.com/a/16824239
 * 
 * @tparam C is the class whose save method is aimed
 * @tparam ARCHIVE is the type of the save parameter
 */
template<typename C, typename ARCHIVE>
struct has_save {
private:
    template<typename T>
    static constexpr auto check(T*)
    -> typename
        std::is_same<
            decltype( std::declval<T>().save(std::declval<ARCHIVE&>()) ),
            void    
        >::type;

    template<typename>
    static constexpr std::false_type check(...);

    typedef decltype(check<C>(0)) type;

public:
    static constexpr bool value = type::value;
};

/**
 * @brief A structure to test the presence of a load method
 * 
 * This structure is meant to detect whether a class is 
 * equipped with a static load method. The idea at the basis of 
 * this code was taken from https://stackoverflow.com/a/16824239
 * 
 * @tparam C is the class whose save method is aimed
 * @tparam ARCHIVE is the type of the load parameter
 */
template<typename C, typename ARCHIVE>
struct has_load {
private:
    template<typename T>
    static constexpr auto check(T*)
    -> typename
        std::is_same<
            decltype( T::load(std::declval<ARCHIVE&>()) ),
            T    
        >::type;

    template<typename>
    static constexpr std::false_type check(...);

    typedef decltype(check<C>(0)) type;

public:
    static constexpr bool value = type::value;
};

namespace Archive {

namespace Out {

class Basic {};

class Binary : public Basic 
{
    std::ofstream ofs;
public:
    Binary();

    Binary(std::filesystem::path position);

    inline void open(std::filesystem::path position)
    {
        ofs.open(position, std::fstream::binary | std::fstream::out);
    }

    inline void close()
    {
        ofs.close();
    }

    inline bool is_open() const
    {
        return ofs.is_open();
    }

    Binary& operator&(const std::string& text);

    template<typename ARITHMETIC_TYPE, std::enable_if_t<std::is_arithmetic_v<ARITHMETIC_TYPE>, bool> = true>
    Binary& operator&(const ARITHMETIC_TYPE& value)
    {
        ofs.write((char const*)(&value), sizeof(ARITHMETIC_TYPE));

        return *this;
    }

    template<class T, class Alloc>
    Binary& operator&(const std::vector<T,Alloc>& v)
    {
        *this & v.size();
        for (const T& value : v) {
            *this & value;
        }

        return *this;
    }

    template<class Key, class T, class Compare, class Allocator>
    Binary& operator&(const std::map<Key,T,Compare,Allocator>& m)
    {
        *this & m.size();
        for (const auto& [key,value] : m) {
            *this & key & value;
        }

        return *this;
    }

    template<class T, class Compare>
    Binary& operator&(const std::priority_queue<T,std::vector<T>,Compare>& queue)
    {
        *this & queue.size();
        for (const T& value : queue) {
            *this & value;
        }

        return *this;
    }

    ~Binary();
};

}

namespace In 
{

class Basic {};

class Binary : public Basic 
{
    std::ifstream ifs;
public:
    Binary(std::filesystem::path position);

    inline void close()
    {
        ifs.close();
    }

    inline bool is_open() const
    {
        return ifs.is_open();
    }

    template<typename ARITHMETIC_TYPE, std::enable_if_t<std::is_arithmetic_v<ARITHMETIC_TYPE>, bool> = true>
    Binary& operator&(ARITHMETIC_TYPE& value)
    {
        ifs.read((char *)(&value), sizeof(ARITHMETIC_TYPE));

        return *this;
    }

    Binary& operator&(std::string& text);

    template<class T, class Alloc, std::enable_if_t<has_load<T, Races::Archive::In::Binary>::value, bool> = true>
    Binary& operator&(std::vector<T,Alloc>& v)
    {
        size_t size;

        *this & size;

        v.clear();

        for (size_t i=0; i<size; ++i) {
            v.push_back(T::load(*this));
        }

        return *this;
    }

/*
    template<class Key, class T, class Compare, class Allocator, 
             std::enable_if_t<has_load<Key, Races::Archive::In::Binary>::value, bool> = true>
    Binary& operator&(std::map<Key,T,Compare,Allocator>& m)
    {
        size_t size;

        *this & size;

        m.clear();

        for (size_t i=0; i<size; ++i) {

            m.emplace(key, value);
        }

        return *this;
    }
*/

    template<class T, class Alloc, std::enable_if_t<!has_load<T, Races::Archive::In::Binary>::value, bool> = true>
    Binary& operator&(std::vector<T,Alloc>& v)
    {
        size_t size;

        *this & size;

        v.clear();

        for (size_t i=0; i<size; ++i) {
            T value;

            *this & value;

            v.push_back(std::move(value));
        }

        return *this;
    }

    template<class Key, class T, class Compare, class Allocator, std::enable_if_t<!has_load<T, Races::Archive::In::Binary>::value, bool> = true>
    Binary& operator&(std::map<Key,T,Compare,Allocator>& m)
    {
        size_t size;

        *this & size;

        m.clear();

        for (size_t i=0; i<size; ++i) {
            Key key;
            T value;

            *this & key & value;

            m.emplace(key, value);
        }

        return *this;
    }

    template<class T, class Compare>
    Binary& operator&(std::priority_queue<T,std::vector<T>,Compare>& queue)
    {
        std::vector<T> queue_v;

        *this & queue_v;

        Compare compare;

        queue = std::priority_queue<T,std::vector<T>,Compare>(compare, std::move(queue_v));

        return *this;
    }

    ~Binary();
};

}

}

template<class ARCHIVE, class VALUE_TYPE, std::enable_if_t<has_save<VALUE_TYPE, ARCHIVE>::value &&
                                                           std::is_base_of_v<Archive::Out::Basic, ARCHIVE>, bool> = true>
ARCHIVE& operator&(ARCHIVE& archive, const VALUE_TYPE& value)
{
    value.save(archive);

    return archive;
}

template<class ARCHIVE, class VALUE_TYPE, std::enable_if_t<has_load<VALUE_TYPE, ARCHIVE>::value &&
                                                           std::is_base_of_v<Archive::In::Basic, ARCHIVE>, bool> = true>
ARCHIVE& operator&(ARCHIVE& archive, VALUE_TYPE& value)
{
    value = VALUE_TYPE::load(archive);

    return archive;
}

template<class ARCHIVE, class VALUE_TYPE, std::enable_if_t<std::is_enum_v<VALUE_TYPE> && 
                                                           std::is_base_of_v<Archive::Out::Basic, ARCHIVE>, bool> = true>
ARCHIVE& operator&(ARCHIVE& archive, const VALUE_TYPE& type)
{
    archive & static_cast<typename std::underlying_type<VALUE_TYPE>::type>(type);

    return archive;
}

template<class ARCHIVE, class VALUE_TYPE, std::enable_if_t<std::is_enum_v<VALUE_TYPE> && 
                                                           std::is_base_of_v<Archive::In::Basic, ARCHIVE>, bool> = true>
ARCHIVE& operator&(ARCHIVE& archive, VALUE_TYPE& type)
{
    typename std::underlying_type<VALUE_TYPE>::type value;
    
    archive & value;

    type = static_cast<VALUE_TYPE>(value);

    return archive;
}

}

#endif // __RACES_ARCHIVE__