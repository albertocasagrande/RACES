/**
 * @file archive.cpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Define some archive classes and their methods
 * @version 0.3
 * @date 2023-07-11
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

#include "archive.hpp"

namespace Races {

namespace Archive {

namespace Basic {

Basic::Basic():
    fs()
{}

Basic::Basic(std::filesystem::path position, std::ios_base::openmode mode):
    fs(position, mode)
{}

Basic::~Basic()
{
    if (is_open()) {
        close();
    }
}

Out::Out():
    Basic()
{}

Out::Out(std::filesystem::path position):
    Basic(position, std::fstream::out)
{}

Out::Out(std::filesystem::path position, std::ios_base::openmode mode):
    Basic(position, mode | std::fstream::out)
{}

In::In(std::filesystem::path position):
    Basic(position, std::fstream::in)
{}

In::In(std::filesystem::path position, std::ios_base::openmode mode):
    Basic(position, mode | std::fstream::in)
{}

std::streampos In::size()
{
    // store current position
    const auto pos = fs.tellg();

    // seek the end
    fs.seekg(0, std::ios::end);

    // read the number of bytes from the beginning
    const auto file_bytes = fs.tellg();

    // move to the original position
    fs.seekg(pos, std::ios::beg);

    return file_bytes;
}

}  // Basic

namespace Binary {

Out::Out():
    Archive::Basic::Out()
{}

Out::Out(std::filesystem::path position):
    Archive::Basic::Out(position, std::fstream::binary)
{}

Out::Out(std::filesystem::path position, std::ios_base::openmode mode):
    Archive::Basic::Out(position, mode | std::fstream::binary)
{}

Out& Out::operator&(const std::string& text)
{
    size_t size = text.size();
    fs.write((char const*)(&size), sizeof(size_t));
    fs.write(text.c_str(), text.size()*sizeof(char));

    return *this;
}

ByteCounter::ByteCounter():
    Out(), byte_counter(0)
{}

ByteCounter& ByteCounter::operator&(const std::string& text)
{
    byte_counter += sizeof(size_t);
    byte_counter += text.size()*sizeof(char);

    return *this;
}

In::In(std::filesystem::path position):
    Archive::Basic::In(position, std::fstream::binary)
{}

In::In(std::filesystem::path position, std::ios_base::openmode mode):
    Archive::Basic::In(position, std::fstream::binary | mode)
{}

In& In::operator&(std::string& text)
{
    size_t size;
    fs.read((char *)(&size), sizeof(size_t));

    char* buffer = new char[size+1];

    fs.read(buffer, size*sizeof(char));

    buffer[size] = '\0';

    text = std::string(buffer);

    delete[] buffer;

    return *this;
}

}

}

}
