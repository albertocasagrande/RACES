/**
 * @file simulator.hpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Define a tumor evolution simulator
 * @version 0.4
 * @date 2023-06-30
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
#include <list>

#include "tissue.hpp"
#include "logger.hpp"
#include "plot_2D.hpp"
#include "tissue_plotter.hpp"
#include "statistics.hpp"

namespace Races {

/**
 * @brief A tumor evolution simulator
 * 
 */
template<typename LOGGER=BasicLogger, typename PLOT_WINDOW=UI::Plot2DWindow>
class BasicSimulator 
{
    struct EventAffectedCells
    {
        std::list<CellInTissue> new_cells;
        std::list<Cell> lost_cells;
    };

    Time time;  //!< simulation time
    Tissue &tissue; //!< simulated tissue
    std::mt19937_64 random_gen;  //!< pseudo-random generator
    Time last_snapshot_time;    //!< time of the last snapshot

    bool logging_enabled;           //!< a flag to pause logging
    LOGGER *logger;                         //!< event logger

    UI::TissuePlotter<PLOT_WINDOW> plotter; //!< tissue plotter UI
    bool quitting;                          //!< quitting flag

    TissueStatistics statistics;    //!< the tissue simulation statistics

    DriverGenotypeId select_epigenetic_clone(const DriverEpigeneticGraph& epigenetic_graph,
                                             const DriverGenotypeId& genotype);

    /**
     * @brief Push cells towards a non-driver mutated cell
     * 
     * @param tissue is the tissue in which cells are pushed
     * @param position is the position from which cells are pushed
     * @param direction is the pushing direction
     * @return the list of the cells that have been pushed outside 
     *      the tissue border. Either this list is empty or 
     *      it contains one single cell.
     */
    std::list<Cell> push_cells(Tissue& tissue, PositionInTissue from_position, const Direction& direction);    

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
     * @return list of the affected cells. Whenever the specified position 
     *      is not in the tissue or does not correspond to a driver 
     *      mutation, the returned list is empty. If one of the two 
     *      siblings is pushed outside the tissue borders, the returned 
     *      list contains only the new cell in position `position`. 
     *      In the remaining case, the returned list contains two cells
     */
    EventAffectedCells duplicate_cell(const Position& position);

    /**
     * @brief Kill the cell in the specified position
     * 
     * This method simulates the death of a cell in tissue.
     * If the cell in the provided position has non-driver 
     * genotype, then nothing is done.
     * 
     * @param position is the position of the cell to be killed
     * @return list of the affected cells, i.e., a list containing 
     *      the former status of the killed cell at most
     */
    EventAffectedCells kill_cell(const Position& position);

    /**
     * @brief Apply an epigenetic event on the cell in the specified position
     * 
     * This method simulates an epigenetic event on a cell in tissue.
     * If the cell in the provided position has non-driver genotype, 
     * then nothing is done.  
     * 
     * @param position is the position of the cell to be methylated/demethylated
     * @param destination_id is the resulting genotype identifier of the cell
     * @return list of the affected cells, i.e., a list containing 
     *      the status of the mutated cell at most
     */
    EventAffectedCells apply_epigenetic_event(const Position& position, const DriverGenotypeId& destination_id);

    /**
     * @brief Duplicate the cell and apply an epigenetic event
     * 
     * This method simulates the duplication of a cell in the 
     * tissue and applies an epigenetic event to one of the siblings.
     * 
     * @param position is the position of the cell to be duplicate
     * @param destination_id is the resulting genotype identifier of the cell
     * @return list of the affected cells. Whenever the specified position 
     *      is not in the tissue or does not correspond to a driver 
     *      mutation, the returned list is empty. If one of the two 
     *      siblings is pushed outside the tissue borders, the returned 
     *      list contains only the new cell in position `position`. 
     *      In the remaining case, the returned list contains two cells
     */
    EventAffectedCells duplicate_cell_and_apply_epigenetic_event(const Position& position, const DriverGenotypeId& destination_id);

    /**
     * @brief Update position of a cell in the species data
     * 
     * @param tissue is the tissue containing the cell
     * @param position is the position of the cell to be updated
     */
    static void update_position_in_species(Tissue& tissue, const PositionInTissue& position);
public:
    Time snapshot_interval;     //!< time interval between two snapshots

    /**
     * @brief The basic simulator constructor
     * 
     * @param tissue is the tissue on which simulation is performed
     * @param logger is a pointer to the simulation logger
     * @param snapshot_interval is the time interval between two snapshots
     * @param random_seed is the simulator random seed
     */
    BasicSimulator(Tissue &tissue, LOGGER *logger, Time snapshot_interval, int random_seed=0);

    /**
     * @brief The basic simulator constructor
     * 
     * @param tissue is the tissue on which simulation is performed
     * @param random_seed is the simulator random seed
     */
    BasicSimulator(Tissue &tissue, int random_seed=0);

    /**
     * @brief The basic simulator constructor
     * 
     * @param tissue is the tissue on which simulation is performed
     * @param logger is a pointer to the simulation logger
     * @param random_seed is the simulator random seed
     */
    BasicSimulator(Tissue &tissue, LOGGER *logger, int random_seed=0);

    /**
     * @brief The basic simulator constructor
     * 
     * @param tissue is the tissue on which simulation is performed
     * @param snapshot_interval is the time interval between two snapshots
     * @param random_seed is the simulator random seed
     */
    BasicSimulator(Tissue &tissue, Time snapshot_interval, int random_seed=0);

    /**
     * @brief Set the simulator logger
     * 
     * @param logger is the simulator logger
     */
    void set_logger(LOGGER *logger);

    /**
     * @brief Enable logging
     */
    void enable_logging();

    /**
     * @brief Disable logging
     */
    void disable_logging();

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
     * @param final_time is the final simulation time
     * @return a reference to the updated simulator
     */
    BasicSimulator<LOGGER,PLOT_WINDOW>& run_up_to(const Time& final_time);

    /**
     * @brief Get the current simulator time
     * 
     * @return a constant reference to the simulator time
     */
    const Time& get_time() const;

    /**
     * @brief Get the simulation plotter 
     * 
     * @return a non-constant reference to the plotter
     */
    UI::TissuePlotter<PLOT_WINDOW>& get_plotter();
};

/* Implementation */

template<typename LOGGER, typename PLOT_WINDOW>
inline const Time& BasicSimulator<LOGGER,PLOT_WINDOW>::get_time() const
{
    return time;
}

template<typename LOGGER, typename PLOT_WINDOW>
BasicSimulator<LOGGER,PLOT_WINDOW>::BasicSimulator(Tissue &tissue, LOGGER *logger, Time snapshot_interval, int random_seed):
    time(0), tissue(tissue), random_gen(), last_snapshot_time(0), 
    logging_enabled(logger!=nullptr), logger(logger), plotter(tissue), 
    quitting(false), statistics(tissue), snapshot_interval(snapshot_interval)
{
    random_gen.seed(random_seed);
}

template<typename LOGGER, typename PLOT_WINDOW>
BasicSimulator<LOGGER,PLOT_WINDOW>::BasicSimulator(Tissue &tissue, int random_seed):
    BasicSimulator(tissue, nullptr, std::numeric_limits<Time>::max(), random_seed)
{
}

template<typename LOGGER, typename PLOT_WINDOW>
BasicSimulator<LOGGER,PLOT_WINDOW>::BasicSimulator(Tissue &tissue, LOGGER *logger, int random_seed):
    BasicSimulator(tissue, logger, std::numeric_limits<Time>::max(), random_seed)
{
}

template<typename LOGGER, typename PLOT_WINDOW>
BasicSimulator<LOGGER,PLOT_WINDOW>::BasicSimulator(Tissue &tissue, Time snapshot_interval, int random_seed):
    BasicSimulator(tissue, nullptr, snapshot_interval, random_seed)
{
}

template<typename LOGGER, typename PLOT_WINDOW>
CellEvent BasicSimulator<LOGGER,PLOT_WINDOW>::select_next_event()
{
    std::uniform_real_distribution<double> uni_dist(0.0, 1.0);

    CellEvent event;
    event.delay = std::numeric_limits<Time>::max();

    for (const auto& species: tissue) {
        const auto num_of_cells{species.num_of_cells()};

        // deal with exclusively somatic events
        for (const auto& event_type: {CellEventType::DIE, CellEventType::DUPLICATE}) {
            const Time r_value = uni_dist(random_gen);
            const Time event_rate = species.get_rate(event_type);
            const Time candidate_delay =  -log(r_value) / (num_of_cells * event_rate);

            if (event.delay>candidate_delay) {
                event.type = event_type;
                event.delay = candidate_delay;
                event.position = Position(tissue, species.choose_a_cell(random_gen));
                event.initial_genotype = species.get_id();
            }
        }

        // Deal with possible epigenetic events whenever admitted
        const auto& epigenetic_graph = tissue.get_epigenetic_graph();
        for (const auto& [dst_id, event_rate]: epigenetic_graph.get_outgoing_edge_map(species.get_id())) {
            const Time r_value = uni_dist(random_gen);
            const Time candidate_delay =  -log(r_value) / (num_of_cells * event_rate);
            
            if (event.delay>candidate_delay) {
                event.type = CellEventType::DUPLICATION_AND_EPIGENETIC_EVENT;
                event.delay = candidate_delay;
                event.position = Position(tissue, species.choose_a_cell(random_gen));
                event.epigenetic_genotype = static_cast<DriverGenotypeId>(dst_id);
                event.initial_genotype = species.get_id();
            }
        }
    }

    return event;
}

template<typename LOGGER, typename PLOT_WINDOW>
typename BasicSimulator<LOGGER,PLOT_WINDOW>::EventAffectedCells 
BasicSimulator<LOGGER,PLOT_WINDOW>::kill_cell(const Position& position)
{
    CellInTissue cell = (*(position.tissue))(position);

    if (!cell.has_driver_genotype()) {
        return {{},{}};
    }

    // remove former cell from the species
    Species& species = position.tissue->get_species(cell.get_driver_genotype());
    species.remove(cell.get_id());

    return {{cell},{}};
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
inline Direction select_2D_random_direction(GENERATOR& random_gen, const Position& position)
{
    std::vector<Direction> directions{Direction::X_UP, Direction::X_DOWN,
                                      Direction::Y_UP, Direction::Y_DOWN,
                                      Direction::X_UP|Direction::Y_UP,
                                      Direction::X_UP|Direction::Y_DOWN,
                                      Direction::X_DOWN|Direction::Y_UP,
                                      Direction::X_DOWN|Direction::Y_DOWN};

    (void)position;

    std::uniform_int_distribution<size_t> distribution(0,directions.size()-1);

    return directions[distribution(random_gen)];
}  

template<typename LOGGER, typename PLOT_WINDOW>
void BasicSimulator<LOGGER,PLOT_WINDOW>::update_position_in_species(Tissue& tissue, const PositionInTissue& position)
{
    CellInTissue& cell = tissue(position);
    Species &species = tissue.get_species(cell.get_driver_genotype());

    species.update_cell(cell);
}

template<typename LOGGER, typename PLOT_WINDOW>
std::list<Cell> BasicSimulator<LOGGER,PLOT_WINDOW>::push_cells(Tissue& tissue, PositionInTissue from_position, const Direction& direction)
{
    PositionDelta delta(direction);
    PositionInTissue to_position(from_position+delta);
    Cell to_be_moved=tissue(from_position);
    while (tissue.is_valid(to_position)&&tissue(to_position).has_driver_genotype()) {   
        swap(tissue(to_position), to_be_moved);

        BasicSimulator<LOGGER,PLOT_WINDOW>::update_position_in_species(tissue, to_position);

        to_position += delta;
    }

    if (!tissue.is_valid(to_position)) {
        if (to_be_moved.get_id() != tissue(from_position).get_id()) {
            Species &species = tissue.get_species(to_be_moved.get_driver_genotype());

            species.remove(to_be_moved.get_id());
        }

        return {to_be_moved};
    }

    swap(tissue(to_position), to_be_moved);

    BasicSimulator<LOGGER,PLOT_WINDOW>::update_position_in_species(tissue, to_position);

    return {};
}

template<typename LOGGER, typename PLOT_WINDOW>
typename BasicSimulator<LOGGER,PLOT_WINDOW>::EventAffectedCells 
BasicSimulator<LOGGER,PLOT_WINDOW>::apply_epigenetic_event(const Position& position, const DriverGenotypeId& destination_id)
{
    const auto& graph = tissue.get_epigenetic_graph();

    CellInTissue& cell = tissue(position);

    DriverGenotypeId genotype(cell.get_driver_genotype());

    const auto& dst_map = graph.get_outgoing_edge_map(genotype);

    Species& species = tissue.get_species(genotype);

    if (dst_map.find(destination_id)==dst_map.end()) {
        std::ostringstream oss;

        oss << "Unsupported epigenetic event from driver genotype " 
            << static_cast<DriverGenotype>(species) << " to " 
            << static_cast<DriverGenotype>(tissue.get_species(destination_id)) << std::endl;
        
        throw std::domain_error(oss.str());
    }

    species.remove(cell.get_id());

    cell.genotype = destination_id;

    tissue.get_species(destination_id).add(cell);

    return {{tissue(position)},{}};
}

template<typename LOGGER, typename PLOT_WINDOW>
typename BasicSimulator<LOGGER,PLOT_WINDOW>::EventAffectedCells 
BasicSimulator<LOGGER,PLOT_WINDOW>::duplicate_cell(const Position& position)
{
    Tissue& tissue = *(position.tissue);
    EventAffectedCells affected;

    CellInTissue& cell1 = tissue(position);

    if (!cell1.has_driver_genotype()) {
        return affected;
    }

    Direction push_dir = select_2D_random_direction(random_gen, position);

    // push the cell in position towards push direction
    affected.lost_cells = push_cells(tissue, position, push_dir);

    DriverGenotypeId genotype(cell1.get_driver_genotype());
    Species& species{tissue.get_species(genotype)};

    auto parent_id = cell1.get_id();

    species.remove(parent_id);

    cell1 = Cell(genotype, parent_id);
    species.add(cell1);

    affected.new_cells.push_back(cell1);

    PositionInTissue new_cell_position{position + PositionDelta(push_dir)};
    if (tissue.is_valid(new_cell_position)) {
        CellInTissue& cell2 = tissue(new_cell_position);
        cell2 = Cell(genotype, parent_id);
        species.add(cell2);

        affected.new_cells.push_back(cell2);
    }

    return affected;
}

template<typename LOGGER, typename PLOT_WINDOW>
typename BasicSimulator<LOGGER,PLOT_WINDOW>::EventAffectedCells 
BasicSimulator<LOGGER,PLOT_WINDOW>::duplicate_cell_and_apply_epigenetic_event(const Position& position, const DriverGenotypeId& destination_id)
{
    auto affected = duplicate_cell(position);

    if (affected.new_cells.size()==0) {
        return affected;
    }

    apply_epigenetic_event(position, destination_id);

    const Tissue& tissue = *position.tissue;

    for (auto& cell : affected.new_cells) {
        cell = tissue(cell);
    }

    return affected;
}

template<typename LOGGER, typename PLOT_WINDOW>
BasicSimulator<LOGGER,PLOT_WINDOW>& BasicSimulator<LOGGER,PLOT_WINDOW>::run_up_to_next_event()
{
    CellEvent event = select_next_event();

    time += event.delay;

    EventAffectedCells affected;

    switch(event.type) {
        case CellEventType::DIE:
            affected = kill_cell(event.position);
            break;
        case CellEventType::DUPLICATE:
            affected = duplicate_cell(event.position);
            break;
        case CellEventType::DUPLICATION_AND_EPIGENETIC_EVENT:
            affected = duplicate_cell_and_apply_epigenetic_event(event.position, event.epigenetic_genotype);
            break;
        default:
            throw std::runtime_error("Unhandled event type");
    }

    if (logging_enabled) {
        for (const auto& cell : affected.new_cells) {
            logger->record(event.type, cell, time);
        }
    }

    statistics.record_event(event, time);

    for (const auto& cell: affected.lost_cells) {
        statistics.record_lost(cell.get_driver_genotype(), time);
    }

    if (!plotter.closed()) {
        plotter.plot(statistics);
    }
    
    quitting = quitting || plotter.closed();

    return *this;
}

template<typename LOGGER, typename PLOT_WINDOW>
BasicSimulator<LOGGER,PLOT_WINDOW>& BasicSimulator<LOGGER,PLOT_WINDOW>::run_up_to(const Time& final_time)
{
    size_t total_cells{1};

    for (const auto& axis_size: tissue.size()) {
        total_cells *= axis_size;
    }

    while (!quitting && tissue.num_of_mutated_cells()>0 && time < final_time) {
        run_up_to_next_event();

        if (last_snapshot_time+snapshot_interval<time) {
            last_snapshot_time = time;
            if (logging_enabled) {
                logger->snapshot(tissue, last_snapshot_time);
            }
        }
    }

    while (!quitting) {
        plotter.plot(statistics);
        quitting = !plotter.waiting_end();
    }

    if (logging_enabled) {
        logger->snapshot(tissue, time);
    }

    return *this;
}

template<typename LOGGER, typename PLOT_WINDOW>
void BasicSimulator<LOGGER,PLOT_WINDOW>::set_logger(LOGGER *logger)
{
    logging_enabled = logging_enabled | (this->logger==nullptr);
    
    this->logger = logger;
}

template<typename LOGGER, typename PLOT_WINDOW>
void BasicSimulator<LOGGER,PLOT_WINDOW>::enable_logging()
{
    logging_enabled = (this->logger!=nullptr);
}

template<typename LOGGER, typename PLOT_WINDOW>
void BasicSimulator<LOGGER,PLOT_WINDOW>::disable_logging()
{
    logging_enabled = false;
}

template<typename LOGGER, typename PLOT_WINDOW>
UI::TissuePlotter<PLOT_WINDOW>& BasicSimulator<LOGGER,PLOT_WINDOW>::get_plotter()
{
    return plotter;
}

}

#endif // __RACES_SIMULATOR__