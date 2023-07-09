/**
 * @file archive.cpp
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

#include "archive.hpp"

namespace Races {

namespace Archive {

namespace Out {

Binary::Binary()
{}

Binary::Binary(std::filesystem::path position):
    Binary()
{
    open(position);
}

Binary& Binary::operator&(const std::string& text)
{
    size_t size = text.size();
    ofs.write((char const*)(&size), sizeof(size_t));
    ofs.write(text.c_str(), text.size()*sizeof(char));

    return *this;
}

Binary::~Binary()
{
    if (is_open()) {
        close();
    }
}

}

namespace In 
{

Binary::Binary(std::filesystem::path position):
    ifs(position, std::fstream::binary | std::fstream::in)
{}

Binary& Binary::operator&(std::string& text)
{
    size_t size;
    ifs.read((char *)(&size), sizeof(size_t));

    char* buffer = new char[size+1];

    ifs.read(buffer, size*sizeof(char));

    buffer[size] = '\0';

    text = std::string(buffer);

    delete[] buffer;

    return *this;
}

Binary::~Binary()
{
    if (is_open()) {
        close();
    }
}

}

}

}
