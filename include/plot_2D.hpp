/**
 * @file plot_2D.hpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Implement a 2D plot window
 * @version 0.1
 * @date 2023-06-14
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

#ifndef __RACES_PLOT_2D__
#define __RACES_PLOT_2D__

#include <string>
#include <vector>
#include <cstdint>

namespace Races 
{

namespace UI 
{

struct Color {
	uint8_t red;
	uint8_t green;
	uint8_t blue;
	uint8_t alpha;

	Color();

	Color(const uint8_t red, const uint8_t green, const uint8_t blue);

	Color(const uint8_t red, const uint8_t green, const uint8_t blue, const uint8_t alpha);
};

class Plot2DWindow {
protected:
	Color background_color;
public:
	Plot2DWindow(const unsigned int width, const unsigned int height, const std::string name);

	void clear();

	Color get_color() const;

	void set_color(const Color& color);

	Color get_background_color() const;

	void set_background_color(const Color& color);

	void draw_point(const unsigned int x, const unsigned int y);

	void draw_rectangle(const unsigned int upper_left_x, const unsigned int upper_left_y,
						const unsigned int width, const unsigned int height, const unsigned int thickness=1);

	void draw_filled_rectangle(const unsigned int upper_left_x, const unsigned int upper_left_y,
							   const unsigned int width, const unsigned int height);

	void get_text_size(const std::string& text, unsigned int& width, unsigned int& height);

	void draw_text(const std::string& text, const unsigned int upper_left_x, const unsigned int upper_left_y);

	void delete_point(const unsigned int x, const unsigned int y);

	bool update();

	static constexpr uint8_t dimensions();
};

inline constexpr uint8_t Plot2DWindow::dimensions()
{
	return 2;
}

inline void Plot2DWindow::clear()
{
}
	
inline void Plot2DWindow::set_color(const Color& color)
{
	(void)color;
}

inline void Plot2DWindow::set_background_color(const Color& color)
{
	background_color = color;
}

	
inline Color Plot2DWindow::get_color() const
{
	return Color();
}

inline Color Plot2DWindow::get_background_color() const
{
	return background_color;
}

inline void Plot2DWindow::draw_point(const unsigned int x, const unsigned int y)
{
	(void)x;
	(void)y;
}

inline void Plot2DWindow::draw_rectangle(const unsigned int upper_left_x, const unsigned int upper_left_y,
						   				 const unsigned int width, const unsigned int height, const unsigned int thickness)
{
	(void)upper_left_x;
	(void)upper_left_y;
	(void)width;
	(void)height;
	(void)thickness;
}

inline void Plot2DWindow::draw_filled_rectangle(const unsigned int upper_left_x, const unsigned int upper_left_y,
								  				const unsigned int width, const unsigned int height)
{
	(void)upper_left_x;
	(void)upper_left_y;
	(void)width;
	(void)height;
}

inline void Plot2DWindow::get_text_size(const std::string& text, unsigned int& width, unsigned int& height)
{
	(void)text;
	(void)width;
	(void)height;
}

inline void Plot2DWindow::draw_text(const std::string& text, const unsigned int upper_left_x, const unsigned int upper_left_y)
{
	(void)upper_left_x;
	(void)upper_left_y;
	(void)text;
}

inline void Plot2DWindow::delete_point(const unsigned int x, const unsigned int y)
{
	(void)x;
	(void)y;
}

inline bool Plot2DWindow::update()
{
	return true;
}
}

}

#endif // __RACES_PLOT_2D__