/**
 * @file simulator.hpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Define a tumor evolution simulator
 * @version 0.1
 * @date 2023-06-09
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

#ifndef __RACES_SIMULATOR__
#define __RACES_SIMULATOR__

#include <random>

#include "tissue.hpp"
#include "logger.hpp"
#include "plot_2D.hpp"
#include "tissue_plotter.hpp"

namespace Races {

/**
 * @brief A tumor evolution simulator
 * 
 */
template<typename LOGGER=BasicLogger, typename PLOT_WINDOW=UI::Plot2DWindow>
class BasicSimulator 
{
    Time time;  //!< simulation time
    Tissue &tissue; //!< simulated tissue
    std::mt19937_64 random_gen;  //!< pseudo-random generator
    Time last_snapshot_time;    //!< time of the last snapshot
    size_t num_of_mutated_cells; //!< number of mutated cells

    LOGGER logger;                          //!< event logger
    UI::TissuePlotter<PLOT_WINDOW> plotter; //!< tissue plotter UI

    bool quitting;                          //!< quitting flag

	const size_t plot_time_interval;	    //!< epochs between two plots
	size_t epochs_since_last_plot;          //!< epochs since the last plot

    DriverGenotypeId select_epigenetic_clone(const DriverEpigeneticGraph& epigenetic_graph,
                                             const DriverGenotypeId& genotype);

    /**
     * @brief Duplicate the cell in the specified position
     * 
     * This method simulates the duplication of a cell in the 
     * tissue. The duplicating cell is identified by its position. 
     * One of the two sibling cells is placed in the parent cell 
     * position; the other one is randomly placed near to the 
     * former parent. This is done by selecting one direction among 
     * the 6 possible directions (i.e., left or right on x-axis, 
     * y-axis, or z-axis) and by pushing all the cells on that 
     * direction one position await from the parent cell.
     * If the cell in the provided position has non-driver 
     * genotype, then nothing is done.  
     * 
     * @param position is the position of the cell to be duplicate
     * @param time is the cell duplication time
     * @return a pair of pointers to the sibling cells
     */
    BasicSimulator<LOGGER,PLOT_WINDOW>& duplicate_cell(Position& position, const Time& time);

    /**
     * @brief Kill the cell in the specified position
     * 
     * This method simulates the death of a cell in tissue.
     * If the cell in the provided position has non-driver 
     * genotype, then nothing is done.
     * 
     * @param position is the position of the cell to be killed
     * @param time is the cell death time
     * @return a reference to the updated object
     */
    BasicSimulator<LOGGER,PLOT_WINDOW>& kill_cell(Position& position, const Time& time);

public:
    Time snapshot_interval;     //!< time interval between two snapshots
    /**
     * @brief The basic simulator constructor
     * 
     * @param tissue is the tissue on which simulation is performed
     * @param logger is the simulation logger
     * @param random_seed is the simulator random seed
     * @param snapshot_interval is the time interval between two snapshots
     */
    BasicSimulator(Tissue &tissue, LOGGER logger=LOGGER(), Time snapshot_interval=std::numeric_limits<Time>::max(), int random_seed=0);

    /**
     * @brief Select the next event
     * 
     * This method select the next event by using the Gillespie's 
     * first reaction method as detailed at page 42 of: 
     *   Gillespie DT. Stochastic simulation of chemical kinetics. 
     *   Annu Rev Phys Chem. 2007;58:35-55.
     *   doi: 10.1146/annurev.physchem.58.032806.104637.PMID: 17037977.
     * 
     * @return a cell event 
     */
    CellEvent select_next_event();

    /**
     * @brief Simulate up to the next event
     * 
     * @return a reference to the updated simulator
     */
    BasicSimulator<LOGGER,PLOT_WINDOW>& run_up_to_next_event();

    /**
     * @brief Simulate up to 
     * 
     * @return a reference to the updated simulator
     */
    BasicSimulator<LOGGER,PLOT_WINDOW>& run_up_to(const Time& final_time);

    /**
     * @brief Get the current simulator time
     * 
     * @return a constant reference to the simulator time
     */
    const Time& get_time() const;
};

/* Implementation */

template<typename LOGGER, typename PLOT_WINDOW>
inline const Time& BasicSimulator<LOGGER,PLOT_WINDOW>::get_time() const
{
    return time;
}

template<typename LOGGER, typename PLOT_WINDOW>
BasicSimulator<LOGGER,PLOT_WINDOW>::BasicSimulator(Tissue &tissue, LOGGER logger, Time snapshot_interval, int random_seed):
    time(0), tissue(tissue), random_gen(), last_snapshot_time(0), num_of_mutated_cells(0), 
    logger(logger), plotter(tissue), quitting(false), plot_time_interval(1000), 
    epochs_since_last_plot(0), snapshot_interval(snapshot_interval)
{
    random_gen.seed(random_seed);
}

template<typename LOGGER, typename PLOT_WINDOW>
CellEvent BasicSimulator<LOGGER,PLOT_WINDOW>::select_next_event()
{
    std::uniform_real_distribution<double> uni_dist(0.0, 1.0);

    std::vector<CellEventType> event_types{CellEventType::DIE, CellEventType::DUPLICATE};

    CellEvent event;

    event.delay = std::numeric_limits<double>::max();

    for (const auto& species: tissue) {
        const auto num_of_cells{species.num_of_cells()};
        for (const auto& event_type: event_types) {
            const double r_value = uni_dist(random_gen);
            const double event_rate = species.get_rate(event_type);
            const double candidate_delay =  -log(r_value) / (num_of_cells * event_rate);

            if (event.delay>candidate_delay) {
                event.type = event_type;
                event.delay = candidate_delay;
                event.position = species.choose_a_cell(random_gen);
            }
        }
    }

    return event;
}

template<typename LOGGER, typename PLOT_WINDOW>
BasicSimulator<LOGGER,PLOT_WINDOW>& BasicSimulator<LOGGER,PLOT_WINDOW>::kill_cell(Position& position, const Time& time)
{
    CellInTissue& cell = (*(position.tissue))(position);

    if (!cell.has_driver_genotype()) {
        return *this;
    }
    num_of_mutated_cells -= 1;

    logger.record(CellEventType::DIE, cell, time);

    // remove former cell from the species
    Species& species = cell.get_species();
    species.remove(cell.get_id(), time);

    return *this;
}

template<typename LOGGER, typename PLOT_WINDOW>
DriverGenotypeId BasicSimulator<LOGGER,PLOT_WINDOW>::select_epigenetic_clone(const DriverEpigeneticGraph& epigenetic_graph, const DriverGenotypeId& genotype)
{
    const auto& dst_map = epigenetic_graph.get_outgoing_edge_map(genotype);

    std::uniform_real_distribution<double> uni_dist(0,1);

    for (const auto& [new_genotype, probability]: dst_map) {
        if (uni_dist(random_gen)<probability) {
            return new_genotype;
        }
    }

    return genotype;
}

template<typename GENERATOR>
inline Direction select_random_direction(GENERATOR& random_gen, const Position& position)
{
    const Tissue& tissue = *(position.tissue);

    std::vector<Direction> directions{Direction::X_LEFT,
                                      Direction::X_RIGHT,
                                      Direction::Y_LEFT,
                                      Direction::Y_RIGHT,
                                      Direction::Z_LEFT,
                                      Direction::Z_RIGHT};
    
     std::vector<Direction> valid_dirs;
    for (const Direction& dir: directions) {
        auto free_cell_pos = tissue.get_non_driver_in_direction(position, dir);

        if (tissue.is_valid(free_cell_pos)) {
            valid_dirs.push_back(dir);
        }
    };

    if (valid_dirs.size()==0) {
        throw std::runtime_error("No valid direction from position");
    }

    std::uniform_int_distribution<size_t> distribution(0,valid_dirs.size()-1);

    return valid_dirs[distribution(random_gen)];
}  

template<typename LOGGER, typename PLOT_WINDOW>
BasicSimulator<LOGGER,PLOT_WINDOW>& BasicSimulator<LOGGER,PLOT_WINDOW>::duplicate_cell(Position& position, const Time& time)
{
    Tissue& tissue = *(position.tissue);

    CellInTissue& cell1 = tissue(position);

    if (!cell1.has_driver_genotype()) {
        return *this;
    }

    Direction push_dir = select_random_direction(random_gen, position);

    // copy the current cell in push direction
    tissue.push_cells(position, push_dir);

    CellInTissue& cell2 = tissue(position + PositionDelta(push_dir));

    DriverGenotypeId genotype(cell1.get_driver_genotype());

    // remove former cell from the species
    tissue.get_species(genotype).remove(cell1.get_id(), time);

    DriverGenotypeId new_genotype;
    
    // select the epigenetic clone
    new_genotype = select_epigenetic_clone(tissue.get_epigenetic_graph(), genotype);

    // set siblings' parent
    auto parent_id = cell1.get_id();
    cell1 = Cell(new_genotype, parent_id);
    cell2 = Cell(new_genotype, parent_id);

    num_of_mutated_cells +=1;

    // add new cells to the new species
    tissue.get_species(new_genotype).add(cell1, time)
                                    .add(cell2, time);

    logger.record(CellEventType::DUPLICATE, cell1, time);
    logger.record(CellEventType::DUPLICATE, cell2, time);

    return *this;
}


template<typename LOGGER, typename PLOT_WINDOW>
BasicSimulator<LOGGER,PLOT_WINDOW>& BasicSimulator<LOGGER,PLOT_WINDOW>::run_up_to_next_event()
{
    CellEvent event = select_next_event();

    time += event.delay;

    switch(event.type) {
        case CellEventType::DIE:
            kill_cell(event.position, time);
            break;
        case CellEventType::DUPLICATE:
            try {
                duplicate_cell(event.position, time);
            } catch (std::exception&) {
                time -= event.delay;
            }
            break;
        default:
            throw std::runtime_error("Unknown event type");
    }
    
    if (++epochs_since_last_plot>=plot_time_interval) {
        epochs_since_last_plot = 0;
        if (!plotter.closed()) {
            plotter.plot(time);
            quitting = plotter.closed();
        }
    }

    return *this;
}

template<typename LOGGER, typename PLOT_WINDOW>
BasicSimulator<LOGGER,PLOT_WINDOW>& BasicSimulator<LOGGER,PLOT_WINDOW>::run_up_to(const Time& final_time)
{
    size_t total_cells{1};

    for (const auto& axis_size: tissue.size()) {
        total_cells *= axis_size;
    }

    num_of_mutated_cells = tissue.num_of_mutated_cells();
    while (!quitting && num_of_mutated_cells>0 && time < final_time) {
        run_up_to_next_event();

        if (last_snapshot_time+snapshot_interval<time) {
            last_snapshot_time = time;
            logger.snapshot(tissue, last_snapshot_time);
        }
    }

    while (!quitting) {
        plotter.plot(time);
        quitting = plotter.closed();
    }

    return *this;
}

}

#endif // __RACES_SIMULATOR__