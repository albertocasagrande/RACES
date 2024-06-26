/**
 * @file palette.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements the Color class and defines the TissuePlotter species palette
 * @version 1.0
 * @date 2024-06-10
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

#include <vector>

#include "plot_2D.hpp"

namespace RACES
{

namespace UI
{


Color::Color():
	Color(0x00,0x00,0x00,0xFF)
{}

Color::Color(const uint8_t red, const uint8_t green, const uint8_t blue):
	Color(red,green,blue,0xFF)
{}

Color::Color(const uint8_t red, const uint8_t green, const uint8_t blue, const uint8_t alpha):
	red(red), green(green), blue(blue), alpha(alpha)
{}

std::vector<Color> palette{
/*
	{0x00,0x00,0xFF},
	{0xFF,0x00,0x00},
	{0x01,0xFF,0xFE},
	{0xFF,0xA6,0xFE},
*/
	{0x00,0x73,0xC2},
	{0xEF,0xC0,0x00},
	{0x86,0x86,0x86},
	{0xCD,0x53,0x4C},

	{0xFF,0xDB,0x66},
	{0x00,0x64,0x01},
	{0x01,0x00,0x67},
	{0x95,0x00,0x3A},
	{0x00,0x7D,0xB5},
	{0xFF,0x00,0xF6},
	{0xFF,0xEE,0xE8},
	{0x77,0x4D,0x00},
	{0x90,0xFB,0x92},
	{0x00,0x76,0xFF},
	{0xD5,0xFF,0x00},
	{0xFF,0x93,0x7E},
	{0x6A,0x82,0x6C},
	{0xFF,0x02,0x9D},
	{0xFE,0x89,0x00},
	{0x7A,0x47,0x82},
	{0x7E,0x2D,0xD2},
	{0x85,0xA9,0x00},
	{0xFF,0x00,0x56},
	{0xA4,0x24,0x00},
	{0x00,0xAE,0x7E},
	{0x68,0x3D,0x3B},
	{0xBD,0xC6,0xFF},
	{0x26,0x34,0x00},
	{0xBD,0xD3,0x93},
	{0x00,0xB9,0x17},
	{0x9E,0x00,0x8E},
	{0x00,0x15,0x44},
	{0xC2,0x8C,0x9F},
	{0xFF,0x74,0xA3},
	{0x01,0xD0,0xFF},
	{0x00,0x47,0x54},
	{0xE5,0x6F,0xFE},
	{0x78,0x82,0x31},
	{0x0E,0x4C,0xA1},
	{0x91,0xD0,0xCB},
	{0xBE,0x99,0x70},
	{0x96,0x8A,0xE8},
	{0xBB,0x88,0x00},
	{0x43,0x00,0x2C},
	{0xDE,0xFF,0x74},
	{0x00,0xFF,0xC6},
	{0xFF,0xE5,0x02},
	{0x62,0x0E,0x00},
	{0x00,0x8F,0x9C},
	{0x98,0xFF,0x52},
	{0x75,0x44,0xB1},
	{0xB5,0x00,0xFF},
	{0x00,0xFF,0x78},
	{0xFF,0x6E,0x41},
	{0x00,0x5F,0x39},
	{0x6B,0x68,0x82},
	{0x5F,0xAD,0x4E},
	{0xA7,0x57,0x40},
	{0xA5,0xFF,0xD2},
	{0xFF,0xB1,0x67},
	{0x00,0x9B,0xFF},
	{0xE8,0x5E,0xBE}
};

}

}
