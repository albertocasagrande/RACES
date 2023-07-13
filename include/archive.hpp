/**
 * @file archive.hpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Define some archive classes and their methods
 * @version 0.6
 * @date 2023-07-13
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
#include <set>
#include <list>
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

/**
 * @brief The Archive module namespace
 */
namespace Archive 
{

/**
 * @brief The namespace for basic classes
 */
namespace Basic
{

/**
 * @brief A basic archive class
 */
struct Basic
{
    std::fstream fs;  //!< The archive file stream

    /**
     * @brief The empty constructor
     */
    Basic();

    /**
     * @brief A constructor
     * 
     * This method creates a basic archive and opens the corresponding stream.
     * 
     * @param filename is the archive file name
     * @param mode is the archive mode 
     */
    Basic(std::filesystem::path filename, std::ios_base::openmode mode);

    /**
     * @brief Open an archive file
     * 
     * @param filename is the archive file name
     * @param mode is the archive mode 
     */
    inline void open(std::filesystem::path filename, std::ios_base::openmode mode)
    {
        fs.open(filename, mode);
    }

    /**
     * @brief Close the archive
     * 
     */
    inline void close()
    {
        fs.close();
    }

    /**
     * @brief Test whether the archive file is open
     * 
     * @return `true` if and only if the archive is open
     */
    inline bool is_open() const
    {
        return fs.is_open();
    }

    /**
     * @brief The destructor
     */
    ~Basic();
};

/**
 * @brief The basic output archive
 */
struct Out : public Basic
{
    /**
     * @brief The empty constructor
     */
    Out();

    /**
     * @brief A constructor
     * 
     * This method creates an output archive and opens the corresponding stream 
     * in output.
     * 
     * @param filename is the archive file name
     */
    Out(std::filesystem::path filename);

    /**
     * @brief A constructor
     * 
     * This method creates an output archive and opens the corresponding stream 
     * in output.
     * 
     * @param filename is the archive file name
     * @param mode is the archive mode 
     */
    Out(std::filesystem::path filename, std::ios_base::openmode mode);
    
    /**
     * @brief Open an archive file
     * 
     * This method opens the corresponding stream in output.
     * 
     * @param filename is the archive file name
     * @param mode is the archive mode 
     */
    inline void open(std::filesystem::path filename, std::ios_base::openmode mode)
    {
        Basic::open(filename, mode | std::fstream::out);
    }

    /**
     * @brief Flush the archive buffer
     */
    inline Out& flush()
    {
        fs.flush();

        return *this;
    }
};

/**
 * @brief The basic input archive
 */
struct In : public Basic
{
    /**
     * @brief A constructor
     * 
     * This method creates an input archive and opens the corresponding stream 
     * in input.
     * 
     * @param filename is the archive file name
     */
    In(std::filesystem::path filename);

    /**
     * @brief A constructor
     * 
     * This method creates an input archive and opens the corresponding stream 
     * in input.
     * 
     * @param filename is the archive file name
     * @param mode is the archive mode 
     */
    In(std::filesystem::path filename, std::ios_base::openmode mode);
        
    /**
     * @brief Open an archive file
     * 
     * This method opens the corresponding stream in input.
     * 
     * @param filename is the archive file name
     * @param mode is the archive mode 
     */
    inline void open(std::filesystem::path filename, std::ios_base::openmode mode)
    {
        Basic::open(filename, mode | std::fstream::in);
    }

    /**
     * @brief Set the archive position
     * 
     * This method sets the position of the archive stream to the 
     * specified position
     * 
     * @param pos is the aimed position in the archive stream
     * @return a reference to the updated object
     */
    inline In& seekg(std::streampos pos)
    {
        fs.seekg(pos);

        return *this;
    }

    /**
     * @brief Get the archive file size
     * 
     * @return the archive file size
     */
    std::streampos size();
};

}  // Basic

}  // Archive

}  // Races

/**
 * @brief Load an object implementing the static method `load(ARCHIVE&)`
 * 
 * This method overloads `operator&` for the objects implementing 
 * the static method `load(ARCHIVE&)` where `ARCHIVE` is derived 
 * from `Archive::Basic::In`.
 * 
 * @tparam ARCHIVE is the archive type and it is derived from `Archive::Basic::In`
 * @tparam VALUE_TYPE is the type of the value to load
 * @param archive is an output archive
 * @param value is the object in which the value is loaded
 * @return a reference to the updated archive 
 */
template<class ARCHIVE, class VALUE_TYPE, std::enable_if_t<Races::has_load<VALUE_TYPE, ARCHIVE>::value &&
                                                           std::is_base_of_v<Races::Archive::Basic::In, ARCHIVE>, bool> = true>
inline ARCHIVE& operator&(ARCHIVE& archive, VALUE_TYPE& value)
{
    value = VALUE_TYPE::load(archive);

    return archive;
}

/**
 * @brief Save an object implementing the method `save(ARCHIVE&)`
 * 
 * This method overloads `operator&` for the objects implementing 
 * the method `save(ARCHIVE&)` where `ARCHIVE` is derived from 
 * `Archive::Basic::Out`.
 * 
 * @tparam ARCHIVE is the archive type and it is derived from `Archive::Basic::Out`
 * @tparam VALUE_TYPE is the type of the value to save
 * @param archive is an output archive
 * @param value is the value to save
 * @return a reference to the updated archive 
 */
template<class ARCHIVE, class VALUE_TYPE, std::enable_if_t<Races::has_save<VALUE_TYPE, ARCHIVE>::value &&
                                                           std::is_base_of_v<Races::Archive::Basic::Out, ARCHIVE>, bool> = true>
inline ARCHIVE& operator&(ARCHIVE& archive, const VALUE_TYPE& value)
{
    value.save(archive);

    return archive;
}

/**
 * @brief Save an enum class object
 * 
 * @tparam ARCHIVE is the archive type and it is derived from `Archive::Basic::Out`
 * @tparam VALUE_TYPE is the enum class type of the value to save
 * @param archive is an output archive
 * @param value is the value to save
 * @return a reference to the updated archive 
 */
template<class ARCHIVE, class VALUE_TYPE, std::enable_if_t<std::is_enum_v<VALUE_TYPE> && 
                                                           std::is_base_of_v<Races::Archive::Basic::Out, ARCHIVE>, bool> = true>
inline ARCHIVE& operator&(ARCHIVE& archive, const VALUE_TYPE& type)
{
    archive & static_cast<typename std::underlying_type<VALUE_TYPE>::type>(type);

    return archive;
}

/**
 * @brief Load an enum class object
 * 
 * @tparam ARCHIVE is the archive type and it is derived from `Archive::Basic::In`
 * @tparam VALUE_TYPE is the enum class type of the value to load
 * @param archive is an output archive
 * @param value is the object in which the value is loaded
 * @return a reference to the updated archive 
 */
template<class ARCHIVE, class VALUE_TYPE, std::enable_if_t<std::is_enum_v<VALUE_TYPE> && 
                                                           std::is_base_of_v<Races::Archive::Basic::In, ARCHIVE>, bool> = true>
inline ARCHIVE& operator&(ARCHIVE& archive, VALUE_TYPE& type)
{
    typename std::underlying_type<VALUE_TYPE>::type value;
    
    archive & value;

    type = static_cast<VALUE_TYPE>(value);

    return archive;
}

/**
 * @brief Load a vector or a list from the archive 
 * 
 * This method is used only when `T` implements the 
 * static method `load(Archive::Binary::In&)`. 
 * 
 * @tparam ARCHIVE is the input archive type
 * @tparam CONTAINER is the container type, i.e., either `std::list` or `std::vector`
 * @tparam T is the type of the container elements and it implements the static method `load`
 * @tparam Alloc is the type of the container allocator
 * @param archive is the input archive
 * @param container is the object in which the container is loaded
 * @return a reference to the updated archive 
 */
template<class ARCHIVE, template<class,class> class CONTAINER, class T, class Alloc,
         std::enable_if_t<std::is_base_of_v<Races::Archive::Basic::In, ARCHIVE> &&
                            (std::is_same_v<CONTAINER<T,Alloc>, std::list<T,Alloc>> ||
                             std::is_same_v<CONTAINER<T,Alloc>, std::vector<T,Alloc>>) && 
                          Races::has_load<T, ARCHIVE>::value, bool> = true>
ARCHIVE& operator&(ARCHIVE& archive, CONTAINER<T,Alloc>& container)
{
    size_t size;

    archive & size;

    container.clear();

    for (size_t i=0; i<size; ++i) {
        container.push_back(T::load(archive));
    }

    return archive;
}

/**
 * @brief Load a vector or a list from the archive
 * 
 * This method is used only when `T` does NOT implement the 
 * static method `load(Archive::Binary::In&)`.
 * 
 * @tparam ARCHIVE is the input archive type
 * @tparam CONTAINER is the container type, i.e., either `std::list` or `std::vector`
 * @tparam T is the type of the container elements and it does not implement the static method `load`
 * @tparam Alloc is the type of the container allocator
 * @param archive is the input archive
 * @param container is the object in which the container is loaded
 * @return a reference to the updated archive 
 */
template<class ARCHIVE, template<typename,typename> class CONTAINER, class T, class Alloc, 
         std::enable_if_t<std::is_base_of_v<Races::Archive::Basic::In, ARCHIVE> &&
                            (std::is_same_v<CONTAINER<T,Alloc>, std::list<T,Alloc>> ||
                             std::is_same_v<CONTAINER<T,Alloc>, std::vector<T,Alloc>>) && 
                          !Races::has_load<T, ARCHIVE>::value, bool> = true>
ARCHIVE& operator&(ARCHIVE& archive, CONTAINER<T,Alloc>& container)
{
    size_t size;

    archive & size;

    container.clear();

    for (size_t i=0; i<size; ++i) {
        T value;

        archive & value;

        container.push_back(std::move(value));
    }

    return archive;
}

/**
 * @brief Load a set from an input archive
 * 
 * This method is used only when `T` implements the 
 * static method `load(Archive::Binary::In&)`.
 * 
 * @tparam ARCHIVE is the input archive type
 * @tparam T is the type of the set elements and it implements the static method `load`
 * @tparam Compare is the comparator type for the set elements
 * @tparam Alloc is the type of the set allocator
 * @param archive is the input archive
 * @param S is the object in which the set is loaded
 * @return a reference to the updated archive 
 */
template<class ARCHIVE, class T, class Compare, class Alloc, 
         std::enable_if_t<std::is_base_of_v<Races::Archive::Basic::In, ARCHIVE> && 
                          Races::has_load<T, ARCHIVE>::value, bool> = true>
ARCHIVE& operator&(ARCHIVE& archive, std::set<T,Compare,Alloc>& S)
{
    size_t size;

    archive & size;

    S = std::set<T,Compare,Alloc>();

    for (size_t i=0; i<size; ++i) {
        S.insert(T::load(archive));
    }

    return archive;
}

/**
 * @brief Load a set from an input archive
 * 
 * This method is used only when `T` does NOT implement the 
 * static method `load(Archive::Binary::In&)`.
 * 
 * @tparam ARCHIVE is the input archive type
 * @tparam T is the type of the set elements and it does not implement the static method `load`
 * @tparam Compare is the comparator type for the set elements
 * @tparam Alloc is the type of the set allocator
 * @param archive is the input archive
 * @param S is the object in which the set is loaded
 * @return a reference to the updated archive 
 */
template<class ARCHIVE, class T, class Compare, class Alloc, 
         std::enable_if_t<std::is_base_of_v<Races::Archive::Basic::In, ARCHIVE> && 
                          !Races::has_load<T, ARCHIVE>::value, bool> = true>
ARCHIVE& operator&(ARCHIVE& archive, std::set<T,Compare,Alloc>& S)
{
    size_t size;

    archive & size;

    S = std::set<T,Compare,Alloc>();

    for (size_t i=0; i<size; ++i) {
        T value;

        archive & value;

        S.insert(std::move(value));
    }

    return archive;
}

/**
 * @brief Load a map from the archive
 * 
 * @tparam ARCHIVE is the input archive type
 * @tparam Key is the type of the map keys and it does not implement the static method `load`
 * @tparam T is the type of the map values and it does not implement the static method `load`
 * @tparam Compare is the type of the key comparator 
 * @tparam Allocator is the type of the map allocator
 * @param archive is the input archive
 * @param m is the object in which the map is loaded
 * @return a reference to the updated archive 
 */
template<class ARCHIVE, class Key, class T, class Compare, class Allocator, 
            std::enable_if_t<std::is_base_of_v<Races::Archive::Basic::In, ARCHIVE> &&
                             !Races::has_load<Key, ARCHIVE>::value && 
                             !Races::has_load<T, ARCHIVE>::value, bool> = true>
ARCHIVE& operator&(ARCHIVE& archive, std::map<Key,T,Compare,Allocator>& m)
{
    size_t size;

    archive & size;

    m.clear();

    for (size_t i=0; i<size; ++i) {
        Key key;
        T value;

        archive & key & value;

        m.emplace(key, value);
    }

    return archive;
}

/**
 * @brief Load a map from the archive
 * 
 * @tparam ARCHIVE is the input archive type
 * @tparam Key is the type of the map keys and it does not implement the static method `load`
 * @tparam T is the type of the map values and it implements the static method `load`
 * @tparam Compare is the type of the key comparator
 * @tparam Allocator is the type of the map allocator
 * @param archive is the input archive
 * @param m is the object in which the map is loaded
 * @return a reference to the updated archive 
 */
template<class ARCHIVE, class Key, class T, class Compare, class Allocator, 
            std::enable_if_t<std::is_base_of_v<Races::Archive::Basic::In, ARCHIVE> &&
                             !Races::has_load<Key, ARCHIVE>::value && 
                             Races::has_load<T, ARCHIVE>::value, bool> = true>
ARCHIVE& operator&(ARCHIVE& archive, std::map<Key,T,Compare,Allocator>& m)
{
    size_t size;

    archive & size;

    m.clear();

    for (size_t i=0; i<size; ++i) {
        Key key;
        T value;

        archive & key;

        m.emplace(key, T::load(archive));
    }

    return archive;
}

/**
 * @brief Load a priority queue from an input archive
 * 
 * @tparam ARCHIVE is the input archive type
 * @tparam T is the type of values in the priority queue
 * @tparam Compare is the type of the comparator
 * @param archive is the output archive
 * @param queue is the object in which the priority queue is loaded
 * @return a reference to the updated archive 
 */
template<class ARCHIVE, class T, class Compare, std::enable_if_t<std::is_base_of_v<Races::Archive::Basic::In, ARCHIVE>, bool> = true>
ARCHIVE& operator&(ARCHIVE& archive, std::priority_queue<T,std::vector<T>,Compare>& queue)
{
    std::vector<T> queue_v;

    archive & queue_v;

    Compare compare;

    queue = std::priority_queue<T,std::vector<T>,Compare>(compare, std::move(queue_v));

    return archive;
}

/**
 * @brief Save a vector or a list in an output archive
 * 
 * @tparam ARCHIVE is the output archive type
 * @tparam CONTAINER is the container type, i.e., either `std::list` or `std::vector`
 * @tparam T is the type of the container elements
 * @tparam Alloc is the type of the container allocator
 * @param archive is the output archive
 * @param container is the container to save
 * @return a reference to the updated archive 
 */
template<class ARCHIVE, template<class,class> class CONTAINER, class T, class Alloc,
         std::enable_if_t<std::is_base_of_v<Races::Archive::Basic::Out, ARCHIVE> && 
                            (std::is_same_v<CONTAINER<T,Alloc>, std::list<T,Alloc>> ||
                             std::is_same_v<CONTAINER<T,Alloc>, std::vector<T,Alloc>>), bool> = true>
ARCHIVE& operator&(ARCHIVE& archive, const CONTAINER<T,Alloc>& container)
{
    archive & container.size();
    for (const T& value : container) {
        archive & value;
    }

    return archive;
}


/**
 * @brief Save a set in an output archive
 * 
 * @tparam ARCHIVE is the output archive type
 * @tparam T is the type of the set elements
 * @tparam Compare is the comparator type for the set elements
 * @tparam Alloc is the type of the set allocator
 * @param archive is the output archive
 * @param S is the set to save
 * @return a reference to the updated archive 
 */
template<class ARCHIVE, class T, class Compare, class Alloc,
         std::enable_if_t<std::is_base_of_v<Races::Archive::Basic::Out, ARCHIVE>, bool> = true>
ARCHIVE& operator&(ARCHIVE& archive, const std::set<T,Compare,Alloc>& S)
{
    archive & S.size();
    for (const T& value : S) {
        archive & value;
    }

    return archive;
}

/**
 * @brief Save a map in an output archive
 * 
 * @tparam ARCHIVE is the output archive type
 * @tparam Key is the type of the map keys
 * @tparam T is the type of the map values
 * @tparam Compare is the type of the key comparator
 * @tparam Allocator is the type of the map allocator
 * @param archive is the output archive
 * @param m is the map to save
 * @return a reference to the updated archive 
 */
template<class ARCHIVE, class Key, class T, class Compare, class Allocator, 
         std::enable_if_t<std::is_base_of_v<Races::Archive::Basic::Out, ARCHIVE>, bool> = true>
ARCHIVE& operator&(ARCHIVE& archive, const std::map<Key,T,Compare,Allocator>& m)
{
    archive & m.size();
    for (const auto& [key,value] : m) {
        archive & key & value;
    }

    return archive;
}

/**
 * @brief Save a priority queue in an output archive
 * 
 * @tparam ARCHIVE is the output archive type
 * @tparam T is the type of values in the priority queue
 * @tparam Compare is the type of the comparator
 * @param archive is the output archive
 * @param queue is the priority queue to save
 * @return a reference to the updated archive 
 */
template<class ARCHIVE, class T, class Compare, std::enable_if_t<std::is_base_of_v<Races::Archive::Basic::Out, ARCHIVE>, bool> = true>
ARCHIVE& operator&(ARCHIVE& archive, const std::priority_queue<T,std::vector<T>,Compare>& queue)
{
    std::priority_queue<T,std::vector<T>,Compare> queue_copy(queue);
    
    archive & queue_copy.size();
    while (!queue_copy.empty()) {
        archive & queue_copy.top();

        queue_copy.pop();
    }

    return archive;
}

namespace Races
{

/**
 * @brief The Archive module namespace
 */
namespace Archive 
{

/**
 * @brief The binary sub-module namespace
 */
namespace Binary 
{

/**
 * @brief The binary output archive
 */
struct Out : public Archive::Basic::Out
{
    /**
     * @brief The empty constructor
     */
    Out();

    /**
     * @brief A constructor
     * 
     * This method creates an output archive and opens the corresponding 
     * binary stream in output.
     * 
     * @param filename is the archive file name
     */
    Out(std::filesystem::path filename);

    /**
     * @brief A constructor
     * 
     * This method creates an output archive and opens the corresponding 
     * binary stream in output.
     * 
     * @param filename is the archive file name
     * @param mode is the archive mode 
     */
    Out(std::filesystem::path filename, std::ios_base::openmode mode);

    /**
     * @brief Open an archive file
     * 
     * This method opens the corresponding binary stream in output.
     * 
     * @param filename is the archive file name
     */
    inline void open(std::filesystem::path filename)
    {
        open(filename, std::fstream::binary);
    }

    /**
     * @brief Open an archive file
     * 
     * This method opens the corresponding binary stream in output.
     * 
     * @param filename is the archive file name
     * @param mode is the archive mode 
     */
    inline void open(std::filesystem::path filename, std::ios_base::openmode mode)
    {
        Archive::Basic::Out::open(filename, std::fstream::binary | mode);
    }

    /**
     * @brief Save a string in the archive
     * 
     * @param text is the string to save
     * @return a reference to the updated archive 
     */
    Out& operator&(const std::string& text);

    /**
     * @brief Save an arithmetic value in the archive
     * 
     * @tparam ARITHMETIC_TYPE is the type of the value to save
     * @param value is the value to save
     * @return a reference to the updated archive 
     */
    template<typename ARITHMETIC_TYPE, std::enable_if_t<std::is_arithmetic_v<ARITHMETIC_TYPE>, bool> = true>
    inline Out& operator&(const ARITHMETIC_TYPE& value)
    {
        fs.write((char const*)(&value), sizeof(ARITHMETIC_TYPE));

        return *this;
    }
};

/**
 * @brief A class to count the bytes required to store an object
 */
class ByteCounter : public Out
{
    std::streampos byte_counter;     //!< the byte counter 

    /**
     * @brief The constructor
     * 
     * This is the only constructor and it is meant to be 
     * exclusively used by the static method `bytes_in_archive`.
     */
    ByteCounter();

public:
    /**
     * @brief Measure the space required to store a string
     * 
     * @param text is the string whose archive space is required
     * @return a reference to the updated archive 
     */
    ByteCounter& operator&(const std::string& text);

    /**
     * @brief Measure the space required by an arithmetic value
     * 
     * @tparam ARITHMETIC_TYPE is the type of the value
     * @param value is the value whore archive required space is aimed
     * @return a reference to the updated archive
     */
    template<typename ARITHMETIC_TYPE, std::enable_if_t<std::is_arithmetic_v<ARITHMETIC_TYPE>, bool> = true>
    inline ByteCounter& operator&(const ARITHMETIC_TYPE& value)
    {
        (void)value;

        byte_counter += sizeof(ARITHMETIC_TYPE);

        return *this;
    }

    template<typename T>
    static std::streampos bytes_in_archive(const T& object) {
        ByteCounter bc;

        bc & object;

        return bc.byte_counter;
    }
};

struct In : public Archive::Basic::In 
{
    /**
     * @brief A constructor
     * 
     * This method creates a input archive and opens the corresponding 
     * binary stream in input.
     * 
     * @param filename is the archive file name
     */
    In(std::filesystem::path filename);

    /**
     * @brief A constructor
     * 
     * This method creates an input archive and opens the corresponding 
     * binary stream in input.
     * 
     * @param filename is the archive file name
     * @param mode is the archive mode 
     */
    In(std::filesystem::path filename, std::ios_base::openmode mode);
        
    /**
     * @brief Open an input archive file
     * 
     * This method opens the corresponding binary stream in input.
     * 
     * @param filename is the archive file name
     * @param mode is the archive mode 
     */
    inline void open(std::filesystem::path filename)
    {
        open(filename, std::fstream::binary );
    }
        
    /**
     * @brief Open an input archive file
     * 
     * This method opens the corresponding binary stream in input.
     * 
     * @param filename is the archive file name
     */
    inline void open(std::filesystem::path filename, std::ios_base::openmode mode)
    {
        Archive::Basic::In::open(filename, std::fstream::binary | mode);
    }

    /**
     * @brief Load an arithmetic value from the archive
     * 
     * @tparam ARITHMETIC_TYPE is the type of the value to load
     * @param value is the object in which the value is loaded
     * @return a reference to the updated archive 
     */
    template<typename ARITHMETIC_TYPE, std::enable_if_t<std::is_arithmetic_v<ARITHMETIC_TYPE>, bool> = true>
    inline In& operator&(ARITHMETIC_TYPE& value)
    {
        fs.read((char *)(&value), sizeof(ARITHMETIC_TYPE));

        return *this;
    }

    /**
     * @brief Load a string from the archive
     * 
     * @param text is the object in which the string is loaded
     * @return a reference to the updated archive 
     */
    In& operator&(std::string& text);
};

}   // Binary

}   // Archive

}   // Races

#endif // __RACES_ARCHIVE__