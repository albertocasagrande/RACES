/**
 * @file simulator.hpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Define a tumor evolution simulator
 * @version 0.10
 * @date 2023-07-08
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
#include <queue>
#include <vector>

#include "tissue.hpp"
#include "logger.hpp"
#include "plot_2D.hpp"
#include "tissue_plotter.hpp"
#include "statistics.hpp"
#include "somatic_mutations.hpp"

namespace Races {

/**
 * @brief A tumor evolution simulator
 */
template<typename LOGGER=BasicLogger>
class BasicSimulator 
{
    struct EventAffectedCells
    {
        std::list<CellInTissue> new_cells;
        std::list<Cell> lost_cells;
    };

    Time time;  //!< Simulation time
    Tissue &tissue; //!< Simulated tissue
    std::map<EpigeneticGenotypeId, bool> death_enabled;  //!< Death enabled flag

    std::mt19937_64 random_gen;     //!< Pseudo-random generator

    std::priority_queue<TimedSomaticMutation> somatic_mutation_queue;   //!< The timed somatic mutation queue 

    std::chrono::system_clock::time_point last_snapshot_time;    //!< Time of the last snapshot
    long secs_between_snapshots;     //!< The number of minutes between two snapshots

    bool logging_enabled;            //!< A flag to pause logging
    LOGGER *logger;                  //!< Event logger

    bool quitting;                          //!< Quitting flag

    TissueStatistics statistics;    //!< The tissue simulation statistics

    /**
     * @brief Simulate a cell duplication
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
    EventAffectedCells simulate_duplication(const Position& position);

    /**
     * @brief Simulate the death of a cell
     * 
     * This method simulates the death of a cell in tissue.
     * If the cell in the provided position has non-driver 
     * genotype, then nothing is done.
     * 
     * @param position is the position of the cell to be killed
     * @return list of the affected cells, i.e., a list containing 
     *      the former status of the killed cell at most
     */
    EventAffectedCells simulate_death(const Position& position);

    /**
     * @brief Simulate a cell mutation
     * 
     * This method simulates a mutation on a cell in tissue.
     * If the cell in the provided position has non-driver genotype, 
     * then nothing is done.  
     * 
     * @param position is the position of the cell that will mutate
     * @param final_id is the resulting genotype identifier of the cell
     * @return list of the affected cells, i.e., a list containing 
     *      the status of the mutated cell at most
     */
    EventAffectedCells simulate_mutation(const Position& position, const EpigeneticGenotypeId& final_id);

    /**
     * @brief Simulate a duplication and an epigenetic event
     * 
     * This method simulates the duplication of a cell in the 
     * tissue and applies an epigenetic event to one of the siblings.
     * 
     * @param position is the position of the cell to be duplicate
     * @param final_id is the resulting genotype identifier of the cell
     * @return list of the affected cells. Whenever the specified position 
     *      is not in the tissue or does not correspond to a driver 
     *      mutation, the returned list is empty. If one of the two 
     *      siblings is pushed outside the tissue borders, the returned 
     *      list contains only the new cell in position `position`. 
     *      In the remaining case, the returned list contains two cells
     */
    EventAffectedCells simulate_duplication_epigenetic_event(const Position& position, const EpigeneticGenotypeId& final_id);

    /**
     * @brief Randomly select a cell among those having a specified somatic genotype
     * 
     * @tparam GENERATOR is the random number generator type
     * @param generator is a random number generator
     * @param tissue is the tissue in which cell must be choosen
     * @param genotype_id is the identifier of the somatic genotype that must have the selected cell
     * @return whenever the targeted somatic species contain at least 
     *        one cell, a pointer to a randomly selected cell whose 
     *        driver somatic genotype identifier is `genotype_id`. 
     *        If otherwise there are no cells whose driver somatic 
     *        genotype identifier is `genotype_id`, this method 
     *        returns `nullptr`
     * 
     */
    template<typename GENERATOR>
    static const CellInTissue* choose_a_cell_in_somatic_species(const Tissue& tissue, const SomaticGenotypeId& genotype_id, GENERATOR& generator);

public:
    size_t death_activation_level;  //!< The minimum number of cells required to activate death

    /**
     * @brief The basic simulator constructor
     * 
     * @param tissue is the tissue on which simulation is performed
     * @param logger is the simulation logger
     * @param random_seed is the simulator random seed
     */
    BasicSimulator(Tissue &tissue, LOGGER &logger, int random_seed=0);

    /**
     * @brief The basic simulator constructor
     * 
     * @param tissue is the tissue on which simulation is performed
     * @param logger is a pointer to the simulation logger
     * @param random_seed is the simulator random seed
     */
    BasicSimulator(Tissue &tissue, LOGGER *logger=nullptr, int random_seed=0);

    /**
     * @brief Set the simulator logger
     * 
     * @param logger is the simulator logger
     */
    void set_logger(LOGGER *logger);

    /**
     * @brief Set the simulator logger
     * 
     * @param logger is the simulator logger
     */
    void set_logger(LOGGER &logger);

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
     * @brief Add a timed driver somatic mutation
     * 
     * @param src is the source driver somatic genotype
     * @param dst is the destination driver somatic genotype
     * @param time is the mutation timing
     * @return a reference to the updated simulator
     */
    BasicSimulator<LOGGER>& add_somatic_mutation(const SomaticGenotype& src, const SomaticGenotype& dst, const Time time);

    /**
     * @brief Simulate up to the next event
     * 
     * This method simulates a tissue up to the next 
     * event. If the user provide a pointer to a 
     * plotter, then the simulation is also plotted
     * in a graphical window.
     * 
     * @tparam PLOT_WINDOW is the plotting window type
     * @param plotter is a tissue plotter pointer
     * @return a reference to the updated simulator
     */
    template<typename PLOT_WINDOW>
    BasicSimulator<LOGGER>& run_up_to_next_event(UI::TissuePlotter<PLOT_WINDOW>* plotter);

    /**
     * @brief Simulate a tissue up to a given time
     * 
     * This method simulates a tissue up to a given 
     * simulated time. If the user provide a pointer 
     * to a plotter, then the simulation is also plotted
     * in a graphical window.
     * 
     * @tparam PLOT_WINDOW is the plotting window type
     * @param final_time is the final simulation time
     * @param plotter is a tissue plotter pointer
     * @return a reference to the updated simulator
     */
    template<typename PLOT_WINDOW>
    BasicSimulator<LOGGER>& run_up_to(const Time& final_time, UI::TissuePlotter<PLOT_WINDOW>* plotter);

    /**
     * @brief Simulate up to the next event
     * 
     * This method simulates a tissue up to the next 
     * event. The simulation is also plotted in a 
     * graphical window.
     * 
     * @tparam PLOT_WINDOW is the plotting window type
     * @param plotter is a tissue plotter
     * @return a reference to the updated simulator
     */
    template<typename PLOT_WINDOW>
    BasicSimulator<LOGGER>& run_up_to_next_event(UI::TissuePlotter<PLOT_WINDOW>& plotter);

    /**
     * @brief Simulate a tissue up to a given time
     * 
     * This method simulates a tissue up to a given 
     * simulated time. The simulation is also plotted
     * in a graphical window.
     * 
     * @tparam PLOT_WINDOW is the plotting window type
     * @param final_time is the final simulation time
     * @param plotter is a tissue plotter
     * @return a reference to the updated simulator
     */
    template<typename PLOT_WINDOW>
    BasicSimulator<LOGGER>& run_up_to(const Time& final_time, UI::TissuePlotter<PLOT_WINDOW>& plotter);

    /**
     * @brief Simulate up to the next event
     * 
     * This method simulates a tissue up to the next 
     * event. If the user provide a pointer to a 
     * plotter, then the simulation is also plotted
     * in a graphical window.
     * 
     * @return a reference to the updated simulator
     */
    BasicSimulator<LOGGER>& run_up_to_next_event();

    /**
     * @brief Simulate a tissue up to a given time
     * 
     * This method simulates a tissue up to a given 
     * simulated time. If the user provide a pointer 
     * to a plotter, then the simulation is also plotted
     * in a graphical window.
     * 
     * @param final_time is the final simulation time
     * @return a reference to the updated simulator
     */
    BasicSimulator<LOGGER>& run_up_to(const Time& final_time);

    /**
     * @brief Get the current simulator time
     * 
     * @return a constant reference to the simulator time
     */
    const Time& get_time() const;

    /**
     * @brief Set the interval between snapshots
     * 
     * @param time_interval is the time interval between two snapshots
     */
    template<class Rep, class Period>
    void set_interval_between_snapshots(const std::chrono::duration<Rep,Period> time_interval);
};

/* Implementation */

template<typename LOGGER>
inline const Time& BasicSimulator<LOGGER>::get_time() const
{
    return time;
}

template<typename LOGGER>
BasicSimulator<LOGGER>::BasicSimulator(Tissue &tissue, LOGGER *logger, int random_seed):
    time(0), tissue(tissue), random_gen(), last_snapshot_time(std::chrono::system_clock::now()), 
    secs_between_snapshots(0), logging_enabled(logger!=nullptr), logger(logger),
    quitting(false), statistics(tissue), death_activation_level(1)
{
    random_gen.seed(random_seed);

    for (const auto& species: tissue) {
        death_enabled[species.get_id()] = false;
    }
}

template<typename LOGGER>
BasicSimulator<LOGGER>::BasicSimulator(Tissue &tissue, LOGGER &logger, int random_seed):
    BasicSimulator(tissue, &logger, random_seed)
{}

template<typename GENERATOR_TYPE>
void select_liveness_event_in_species(CellEvent& event, Tissue& tissue, 
                                      const std::map<EpigeneticGenotypeId, bool>& death_enabled, 
                                      const Species& species,
                                      std::uniform_real_distribution<double>& uni_dist,
                                      GENERATOR_TYPE& random_gen)
{
    const auto& num_of_cells{species.num_of_cells()};

    std::list<CellEventType> event_types{CellEventType::DUPLICATE};

    if (death_enabled.at(species.get_id())) {
        event_types.push_back(CellEventType::DIE);
    }

    // deal with exclusively somatic events
    for (const auto& event_type: event_types) {
        const Time event_rate = species.get_rate(event_type);
        const Time r_value = uni_dist(random_gen);
        const Time candidate_delay =  -log(r_value) / (num_of_cells * event_rate);
        if (event.delay>candidate_delay) {
            event.type = event_type;
            event.delay = candidate_delay;
            event.position = Position(tissue, species.choose_a_cell(random_gen));
            event.initial_genotype = species.get_id();
        }
    }
}

template<typename GENERATOR_TYPE>
void select_epigenetic_event_in_species(CellEvent& event, Tissue& tissue, const Species& species,
                                        std::uniform_real_distribution<double>& uni_dist, GENERATOR_TYPE& random_gen)
{
    const auto& num_of_cells{species.num_of_cells()};

    // Deal with possible epigenetic events whenever admitted
    for (const auto& [dst_id, event_rate]: species.get_epigenentic_mutation_rates()) {
        const Time r_value = uni_dist(random_gen);
        const Time candidate_delay = -log(r_value) / (num_of_cells * event_rate);
        
        if (event.delay>candidate_delay) {
            event.type = CellEventType::DUPLICATION_AND_EPIGENETIC_EVENT;
            event.delay = candidate_delay;
            event.position = Position(tissue, species.choose_a_cell(random_gen));
            event.initial_genotype = species.get_id();
            event.final_genotype = dst_id;
        }
    }
}

template<typename GENERATOR_TYPE>
void select_next_event_in_species(CellEvent& event, Tissue& tissue,
                                  const std::map<EpigeneticGenotypeId, bool>& death_enabled, 
                                  const Species& species, 
                                  std::uniform_real_distribution<double>& uni_dist,
                                  GENERATOR_TYPE& random_gen)
{
    if (species.num_of_cells()==0) {
        return;
    }

    select_liveness_event_in_species(event, tissue, death_enabled, species, uni_dist, random_gen);
    select_epigenetic_event_in_species(event, tissue, species, uni_dist, random_gen);
}

template<typename LOGGER>
CellEvent BasicSimulator<LOGGER>::select_next_event()
{
    std::uniform_real_distribution<double> uni_dist(0.0, 1.0);

    CellEvent event;
    event.delay = std::numeric_limits<Time>::max();

    for (const Species& species: tissue) {
        select_next_event_in_species(event, tissue, death_enabled, species, uni_dist, random_gen);
    }

    // update candidate event with somatic mutation is possible
    while (!somatic_mutation_queue.empty()) {
        const TimedSomaticMutation& somatic_mutation = somatic_mutation_queue.top();

        if (somatic_mutation.time < event.delay + time) {
            const auto* cell = choose_a_cell_in_somatic_species(tissue, somatic_mutation.initial_id, random_gen);
            if (cell != nullptr) {
                event.delay = (somatic_mutation.time >= time? somatic_mutation.time - time:0);
                event.type = CellEventType::DRIVER_SOMATIC_MUTATION;
                event.position = Position(tissue, *cell);
                event.initial_genotype = cell->get_genotype_id();

                const Species& initial_species = tissue.get_species(event.initial_genotype);

                const auto& somatic_species = tissue.get_somatic_genotype_species(somatic_mutation.final_id);
                const size_t index = SomaticGenotype::signature_to_index(initial_species.get_methylation_signature());

                event.final_genotype = somatic_species[index].get_id();

                somatic_mutation_queue.pop();

                return event;
            }

            somatic_mutation_queue.pop();
        } else {
            return event;
        }
    }

    return event;
}

template<typename LOGGER>
typename BasicSimulator<LOGGER>::EventAffectedCells 
BasicSimulator<LOGGER>::simulate_death(const Position& position)
{
    auto cell = (*(position.tissue))(position);

    if (!cell.has_driver_mutations()) {
        return {{},{}};
    }

    return {{cell.copy_and_kill()},{}};
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

template<typename LOGGER>
typename BasicSimulator<LOGGER>::EventAffectedCells 
BasicSimulator<LOGGER>::simulate_mutation(const Position& position, const EpigeneticGenotypeId& final_id)
{
    CellInTissue cell = tissue(position);

    cell.genotype = final_id;

    tissue(position) = cell;

    return {{tissue(position)},{}};
}

template<typename LOGGER>
typename BasicSimulator<LOGGER>::EventAffectedCells 
BasicSimulator<LOGGER>::simulate_duplication(const Position& position)
{
    Tissue& tissue = *(position.tissue);
    EventAffectedCells affected;

    if (!tissue(position).has_driver_mutations()) {
        return affected;
    }

    Cell parent_cell = tissue(position);

    // push the cell in position towards a random direction
    Direction push_dir = select_2D_random_direction(random_gen, position);
    affected.lost_cells = tissue.push_cells(position, push_dir);
    
    tissue(position) = parent_cell.generate_descendent();

    affected.new_cells.push_back(tissue(position));

    PositionInTissue new_cell_position{position + PositionDelta(push_dir)};
    if (tissue.is_valid(new_cell_position)) {
        tissue(new_cell_position) = parent_cell.generate_descendent();

        affected.new_cells.push_back(tissue(new_cell_position));
    }

    return affected;
}

template<typename LOGGER>
typename BasicSimulator<LOGGER>::EventAffectedCells 
BasicSimulator<LOGGER>::simulate_duplication_epigenetic_event(const Position& position, const EpigeneticGenotypeId& final_id)
{
    auto affected = simulate_duplication(position);

    if (affected.new_cells.size()==0) {
        return affected;
    }

    simulate_mutation(position, final_id);

    const Tissue& tissue = *position.tissue;

    for (auto& cell : affected.new_cells) {
        cell = static_cast<Cell>(tissue(cell));
    }

    return affected;
}

template<typename LOGGER>
template<typename GENERATOR>
const CellInTissue*
BasicSimulator<LOGGER>::choose_a_cell_in_somatic_species(const Tissue& tissue, const SomaticGenotypeId& genotype_id, GENERATOR& generator)
{
    std::vector<size_t> num_of_cells;

    size_t total = 0;

    const auto& species = tissue.get_somatic_genotype_species(genotype_id);
    for (const auto& S : species) {
        total += S.num_of_cells();
        num_of_cells.push_back(total);
    }

    if (total == 0) {
        return nullptr;
    }

    std::uniform_int_distribution<size_t> distribution(1,total);
    const size_t value = distribution(generator);

    size_t i=0;
    while (value > num_of_cells[i]) {
        ++i;
    }

    return &(species[i].choose_a_cell(generator));
}

template<typename LOGGER>
BasicSimulator<LOGGER>&
BasicSimulator<LOGGER>::add_somatic_mutation(const SomaticGenotype& src, const SomaticGenotype& dst, const Time time)
{
    if (src.num_of_promoters()>dst.num_of_promoters()) {
        std::ostringstream oss;

        oss << "Incompatible number of methylable promoters: \"" << src.get_name() 
            << "\" must have at least as many methylable promoters as \"" << dst.get_name() 
            << "\".";
        throw std::domain_error(oss.str());
    }

    somatic_mutation_queue.emplace(src, dst, time);

    return *this;
}

template<typename LOGGER>
template<typename PLOT_WINDOW>
BasicSimulator<LOGGER>& BasicSimulator<LOGGER>::run_up_to_next_event(UI::TissuePlotter<PLOT_WINDOW>* plotter)
{
    CellEvent event = select_next_event();

    time += event.delay;

    EventAffectedCells affected;

    switch(event.type) {
        case CellEventType::DIE:
            affected = simulate_death(event.position);
            break;
        case CellEventType::DUPLICATE:
            affected = simulate_duplication(event.position);
            break;
        case CellEventType::DUPLICATION_AND_EPIGENETIC_EVENT:
            affected = simulate_duplication_epigenetic_event(event.position, event.final_genotype);
            break;
        case CellEventType::DRIVER_SOMATIC_MUTATION:
            affected = simulate_mutation(event.position, event.final_genotype);
            break;
        default:
            throw std::runtime_error("Unhandled event type");
    }

    for (const auto& cell : affected.new_cells) {
        if (logging_enabled) {
            logger->record(event.type, cell, time);
        }

        bool& death_enabled_in_species = death_enabled[cell.get_genotype_id()];
        if (!death_enabled_in_species) {
            Species& species = tissue.get_species(cell.get_genotype_id());
            
            death_enabled_in_species = species.num_of_cells()>=death_activation_level;
        }
    }

    statistics.record_event(event, time);
    for (const auto& cell: affected.lost_cells) {
        statistics.record_lost(cell.get_genotype_id(), time);
    }

    if (plotter != nullptr && !plotter->closed()) {
        plotter->plot(statistics);
    }
    
    quitting = quitting || (plotter != nullptr && plotter->closed());

    return *this;
}

template<typename LOGGER>
template<typename PLOT_WINDOW>
BasicSimulator<LOGGER>& BasicSimulator<LOGGER>::run_up_to(const Time& final_time, UI::TissuePlotter<PLOT_WINDOW>* plotter)
{
    size_t total_cells{1};

    for (const auto& axis_size: tissue.size()) {
        total_cells *= axis_size;
    }

    while (!quitting && tissue.num_of_mutated_cells()>0 && time < final_time) {
        using namespace std::chrono;

        run_up_to_next_event(plotter);

        const auto curr_time = system_clock::now();
        const auto from_last_snapshot = duration_cast<seconds>(curr_time-last_snapshot_time);

        if (secs_between_snapshots>0 &&from_last_snapshot.count()>=secs_between_snapshots) {
            last_snapshot_time = curr_time;
            if (logging_enabled) {
                logger->snapshot(tissue);
            }
        }
    }

    if (plotter != nullptr) {
        while (!quitting) {
            plotter->plot(statistics);
            quitting = !plotter->waiting_end();
        }
    }

    if (logging_enabled) {
        logger->snapshot(tissue);
    }

    return *this;
}

template<typename LOGGER>
template<typename PLOT_WINDOW>
inline BasicSimulator<LOGGER>& BasicSimulator<LOGGER>::run_up_to_next_event(UI::TissuePlotter<PLOT_WINDOW>& plotter)
{
    return run_up_to_next_event(&plotter);
}

template<typename LOGGER>
template<typename PLOT_WINDOW>
inline BasicSimulator<LOGGER>& BasicSimulator<LOGGER>::run_up_to(const Time& final_time, UI::TissuePlotter<PLOT_WINDOW>& plotter)
{
    return run_up_to(final_time, &plotter);
}

template<typename LOGGER>
inline BasicSimulator<LOGGER>& BasicSimulator<LOGGER>::run_up_to_next_event()
{
    return run_up_to_next_event<UI::Plot2DWindow>(nullptr);
}

template<typename LOGGER>
inline BasicSimulator<LOGGER>& BasicSimulator<LOGGER>::run_up_to(const Time& final_time)
{
    return run_up_to<UI::Plot2DWindow>(final_time, nullptr);
}

template<typename LOGGER>
inline void BasicSimulator<LOGGER>::set_logger(LOGGER& logger)
{
    this->logger = &logger;

    enable_logging();
}

template<typename LOGGER>
inline void BasicSimulator<LOGGER>::enable_logging()
{
    logging_enabled = (this->logger!=nullptr);
}

template<typename LOGGER>
inline void BasicSimulator<LOGGER>::disable_logging()
{
    logging_enabled = false;
}
 
template<typename LOGGER>
template<class Rep, class Period>
inline void BasicSimulator<LOGGER>::set_interval_between_snapshots(const std::chrono::duration<Rep,Period> time_interval)
{
    using namespace std::chrono;

    secs_between_snapshots = duration_cast<seconds>(time_interval).count();
}

}

#endif // __RACES_SIMULATOR__