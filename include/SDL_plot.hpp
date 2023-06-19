/**
 * @file SDL_plot.hpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Implement a 2D plot window by using SDL2
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

#ifndef __RACES_SDL_PLOT__
#define __RACES_SDL_PLOT__

#include <SDL.h>
#include <SDL_image.h>
#include <SDL_ttf.h>

#include <iostream>

#include "plot_2D.hpp"

namespace Races 
{

namespace UI 
{

class SDLWindow : public Plot2DWindow
{
	enum {
		RED = 0,
		GREEN = 1,
		BLUE = 2,
		ALPHA
	};

	SDL_Window* window;         //!< graphical window
	SDL_Renderer* renderer;     //!< renderer

	TTF_Font* font;				//!< text font
public:
	SDLWindow(const unsigned int width, const unsigned int height, const std::string name);

	void clear();

	void set_color(const Color& color);

	Color get_color() const;

	void draw_point(const unsigned int x, const unsigned int y);

	void draw_rectangle(const unsigned int upper_left_x, const unsigned int upper_left_y,
						const unsigned int width, const unsigned int height, const unsigned int thickness=1);

	void draw_filled_rectangle(const unsigned int upper_left_x, const unsigned int upper_left_y,
							   const unsigned int width, const unsigned int height);

	void get_text_size(const std::string& text, unsigned int& width, unsigned int& height);

	void draw_text(const std::string& text, const unsigned int upper_left_x, const unsigned int upper_left_y);

	void delete_point(const unsigned int x, const unsigned int y);

	void update();

	~SDLWindow();
};


inline Color SDLWindow::get_color() const
{
	Color color;
	SDL_GetRenderDrawColor( renderer, &color.red, &color.green, 
		      						  &color.blue, &color.alpha );

	return color;
}


inline void SDLWindow::set_color(const Color& color)
{
	SDL_SetRenderDrawColor( renderer, color.red, color.green, color.blue, color.alpha );
}

inline void SDLWindow::draw_point(const unsigned int x, const unsigned int y)
{
	SDL_RenderDrawPoint( renderer, x, y );
}

inline void SDLWindow::clear()
{
	set_color(this->background_color);
	SDL_RenderClear( renderer );
}

}

}

#endif // __RACES_SDL_PLOT__