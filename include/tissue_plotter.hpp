/**
 * @file tissue_plotter.hpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Define a UI window to plot a tissue
 * @version 0.2
 * @date 2023-06-28
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

#ifndef __RACES_TISSUE_PLOTTER__
#define __RACES_TISSUE_PLOTTER__

#include <memory>
#include <sstream>

#include "tissue.hpp"

namespace Races 
{

namespace UI 
{

extern std::vector<Color> driver_palette;

template<typename PLOT_WINDOW>
class TissuePlotter {
	const Tissue& tissue;					//!< the tissue to plot

	std::unique_ptr<PLOT_WINDOW> window;	//!< the plotting window

	const unsigned int frame_border;		//!< frame border width
	const unsigned int frame_thickness;		//!< frame line thickness

	const unsigned int legend_rectangle_width;	//!< legend rectangle width
	const unsigned int legend_rectangle_height;	//!< legend rectangle height
	const unsigned int legend_label_height;		//!< legend label height

	const size_t redraw_interval;	            //!< epochs between two redraws
	size_t epochs_since_last_redraw;            //!< epochs since the last redraw

	/**
	 * @brief Get the total width of the plotting window
	 * 
	 * @return the total width of the plotting window in pixels
	 */
	unsigned int total_width() const
	{
		return static_cast<unsigned int>(tissue.size()[0])
				+3*frame_border+legend_rectangle_width;
	}

	/**
	 * @brief Get the total height of the plotting window
	 * 
	 * @return the total height of the plotting window in pixels
	 */
	unsigned int total_height() const
	{
		unsigned int height =  tissue.num_of_species()*(legend_rectangle_height+frame_border);

		return std::max(height, static_cast<unsigned int>(tissue.size()[1]+2*frame_border+legend_rectangle_height));
	}

	/**
	 * @brief Draw the tissue frame
	 */
	void draw_tissue_frame()
	{
		auto sizes(tissue.size());
		if constexpr(PLOT_WINDOW::dimensions()==2) {
			window->set_color(Color(0xFF,0xFF,0xFF));
			window->draw_rectangle(frame_border-frame_thickness, 
								   legend_rectangle_height+frame_border-frame_thickness,
							       sizes[0]+2*frame_thickness, sizes[1]+2*frame_thickness);
		}
	}

	/**
	 * @brief Draw the elapsed time
	 * 
	 * @param time is the time to be printed
	 */
	void draw_time(const Time& time)
	{
		unsigned int label_width, label_height;

		std::ostringstream oss;
		oss << "Time: " << time << "s";

		window->get_text_size(oss.str(), label_width, label_height);
		const unsigned int label_x_pos = frame_border;
		const unsigned int label_y_pos = (frame_border+legend_rectangle_height-label_height)/2;

		window->set_color(Color(0xFF,0xFF,0xFF));
		window->draw_text(oss.str(), label_x_pos, label_y_pos);
	}

	/**
	 * @brief Draw the tissue
	 */
	void draw_tissue()
	{	
		draw_tissue_frame();
	
		size_t species_idx=0;
		for (const auto& species: tissue) {
			const auto& color = driver_palette[++species_idx];

			window->set_color(color);

			for (const auto& cell: species) {
				if constexpr(PLOT_WINDOW::dimensions()==2) {
					window->draw_point(cell.x+frame_border,
									   cell.y+frame_border+legend_rectangle_height);
				} else {
					window->draw_point(cell.x+frame_border,
									   cell.y+frame_border,
									   cell.z+frame_border);
				}
			}
		}
	}

	/**
	 * @brief Draw the leggend
	 */
	void draw_legend()
	{
		auto sizes(tissue.size());
		const unsigned int x = 2*frame_border + static_cast<unsigned int>(sizes[0]);
		unsigned int y = frame_border;

		size_t species_idx=0;
		for (const auto& species: tissue) {

			(void)species;
	
			window->set_color(driver_palette[++species_idx]);
			window->draw_filled_rectangle(x,y,legend_rectangle_width,legend_rectangle_height);

			window->set_color(window->get_background_color());

			std::string label{species.get_name()+(species.is_methylated()?"+":"-")};

			unsigned int label_width, label_height;

			window->get_text_size(label, label_width, label_height);
			const unsigned int label_x_pos = x+(legend_rectangle_width-label_height)/2;
			const unsigned int label_y_pos = y+(legend_rectangle_height-label_height)/2;
			window->draw_text(label, label_x_pos, label_y_pos);

			y += legend_rectangle_height+frame_border;
		}
	}

public:
	/**
	 * @brief A constructor
	 * 
	 * @param tissue is the tissue to plot
	 */
	TissuePlotter(const Tissue& tissue, const size_t redraw_interval=2000):
		TissuePlotter<PLOT_WINDOW>(tissue, "RACES Simulator"+((tissue.get_name()=="")?"":" - "+tissue.get_name()), redraw_interval)
	{
	}

	/**
	 * @brief A constructor
	 * 
	 * @param tissue is the tissue to plot
	 * @param name is the name of the plotting window
	 */
	TissuePlotter(const Tissue& tissue, const std::string name, const size_t redraw_interval=2000):
		tissue(tissue), frame_border(20), frame_thickness(5), legend_rectangle_width(200),
		legend_rectangle_height(35), legend_label_height(16), redraw_interval(redraw_interval), 
		epochs_since_last_redraw(redraw_interval)
	{
		if (tissue.num_of_species()>driver_palette.size()) {
			throw std::domain_error("The color palette does not support so many species");
		}

		auto sizes(tissue.size());
		if constexpr(PLOT_WINDOW::dimensions()==2) {
			window = std::make_unique<PLOT_WINDOW>(total_width(), total_height(), name);
			plot(0);
			
			return;
		}

		throw std::runtime_error("Unsupported plot window type");
	}

	/**
	 * @brief Plot the tissue
	 * 
	 * @param time is the time to report
	 */
    void plot(const Time& time)
	{
		if (++epochs_since_last_redraw>=redraw_interval) {
			epochs_since_last_redraw = 0;
			window->clear();
			
			draw_time(time);

			draw_legend();
			draw_tissue();
			
			window->update();
		}
	}

	/**
	 * @brief Test whether the plotting window has been closed
	 * 
	 * @return `true` if and only if the plotting window has been closed
	 */
	inline const bool& closed() const
	{
		return window->closed();
	}

	/**
	 * @brief Test whenever the plotting window is waiting for any event
	 * 
	 * @return `true` if and only if the plotting window is waiting 
	 *     for any event
	 */
	inline bool waiting_end() const
	{
		return window->waiting_end();
	}

	/**
	 * @brief The destroyer
	 */
	~TissuePlotter()
	{}
};

}

}

#endif // __RACES_TISSUE_PLOTTER__