/**
 * @file SDL_plot.hpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Implement a 2D plot window by using SDL2
 * @version 0.1
 * @date 2023-05-30
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

#include <iostream>

#include "SDL_plot.hpp"

namespace Races 
{

namespace UI 
{

SDLWindow::SDLWindow(const unsigned int width, const unsigned int height, const std::string name):
	Plot2DWindow(width,height,name)
{
	if(SDL_Init(SDL_INIT_VIDEO) < 0) {
		throw std::runtime_error( "SDL initialization error" );
	}

	if (TTF_Init() < 0 ) {
		throw std::runtime_error("SDL_ttf initialization error");
	}

	font = TTF_OpenFont("fonts/Roboto-Regular.ttf", 24);
	if ( !font ) {
		throw std::runtime_error("Failed to load font.");
	}

	SDL_SetHint( SDL_HINT_RENDER_SCALE_QUALITY, "1");

	//Create window
	window = SDL_CreateWindow(name.c_str(), SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED,
								width, height, SDL_WINDOW_SHOWN);
	if (window == NULL) {
		throw std::runtime_error( "Window could not be created" );
	}

	//Create renderer for window
	renderer = SDL_CreateRenderer( window, -1, SDL_RENDERER_ACCELERATED );
	if( renderer == NULL ) {
		throw std::runtime_error( "Renderer could not be created" );
	}
	//Initialize renderer color
	set_color(this->background_color);
}


void SDLWindow::draw_rectangle(const unsigned int upper_left_x, const unsigned int upper_left_y,
					           const unsigned int width, const unsigned int height, const unsigned int thickness)
{

	SDL_Rect outer  = { static_cast<int>(upper_left_x), 
						static_cast<int>(upper_left_y), 
					    static_cast<int>(width), 
					    static_cast<int>(height) };

	SDL_RenderFillRect( renderer, &outer );

	if (width>2*thickness &&  height>2*thickness) {
		Color color = get_color();

		set_color(this->background_color);
		SDL_Rect inner = { static_cast<int>(upper_left_x+thickness), 
						   static_cast<int>(upper_left_y+thickness), 
						   static_cast<int>(width-2*thickness), 
						   static_cast<int>(height-2*thickness) };

		SDL_RenderFillRect( renderer, &inner );
		set_color(color);
	}
}

void SDLWindow::draw_filled_rectangle(const unsigned int upper_left_x, const unsigned int upper_left_y,
							          const unsigned int width, const unsigned int height)
{
	SDL_Rect rect  = { static_cast<int>(upper_left_x), 
						static_cast<int>(upper_left_y), 
					    static_cast<int>(width), 
					    static_cast<int>(height) };

	SDL_RenderFillRect( renderer, &rect );
}

void SDLWindow::delete_point(const unsigned int x, const unsigned int y)
{
	Color color = get_color();
	
	// set background color 
	set_color(this->background_color);

	// plot the point
	draw_point(x,y);

	set_color(color);
}

void SDLWindow::get_text_size(const std::string& text, unsigned int& width, unsigned int& height)
{
	int w, h;
	TTF_SizeText(font, text.c_str(), &w, &h);

	width = static_cast<unsigned int>(w);
	height = static_cast<unsigned int>(h);
}

void SDLWindow::draw_text(const std::string& text, 
						  const unsigned int upper_left_x, const unsigned int upper_left_y)
{
	Color color = get_color();

	SDL_Surface* msg_sfc = TTF_RenderText_Blended(font, text.c_str(),
												  {color.red, color.green, 
												   color.blue, color.alpha}); 

	// now you can convert it into a texture
	SDL_Texture* msg = SDL_CreateTextureFromSurface(renderer, msg_sfc);

	unsigned int width, height;

	get_text_size(text, width, height);

	SDL_Rect msg_rect  = { static_cast<int>(upper_left_x), 
						   static_cast<int>(upper_left_y), 
					       static_cast<int>(width), 
					       static_cast<int>(height) };

	SDL_RenderCopy(renderer, msg, NULL, &msg_rect);

	SDL_FreeSurface(msg_sfc);
	SDL_DestroyTexture(msg);
}

void SDLWindow::update()
{
	SDL_RenderPresent( renderer );

	SDL_Event e;
	while(SDL_PollEvent(&e)!=0) {
		//User requests quit
		if( e.type == SDL_QUIT ) {
			this->w_closed = true;

			return;
		}
    }
}

SDLWindow::~SDLWindow()
{
	//Destroy window	
	SDL_DestroyRenderer( renderer );
	SDL_DestroyWindow( window );

	//Quit SDL subsystems
	IMG_Quit();
	TTF_CloseFont(font);
	SDL_Quit();
}

}

}