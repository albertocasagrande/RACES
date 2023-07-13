/**
 * @file simulation.hpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Define a tumor evolution simulation
 * @version 0.1
 * @date 2023-07-13
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
#include "binary_logger.hpp"
#include "plot_2D.hpp"
#include "tissue_plotter.hpp"
#include "statistics.hpp"
#include "somatic_mutations.hpp"

namespace Races {

/**
 * @brief A tumor evolution simulation
 */
class Simulation 
{
    struct EventAffectedCells
    {
        std::list<CellInTissue> new_cells;
        std::list<Cell> lost_cells;
    };

    using SomaticQueue = std::priority_queue<TimedSomaticMutation>;
    using system_clock = std::chrono::system_clock;

    std::vector<Tissue> tissues;     //!< Simulated tissues
    BinaryLogger *logger;            //!< Event logger

    std::vector<Direction> valid_directions;   //!< valid simulation tissue directions

    system_clock::time_point last_snapshot_time;    //!< Time of the last snapshot
    long secs_between_snapshots;    //!< The number of minutes between two snapshots

    TissueStatistics statistics;     //!< The tissue simulation statistics

    Time time;                       //!< Simulation time
    std::mt19937_64 random_gen;      //!< Pseudo-random generator

    SomaticQueue somatic_mutation_queue;   //!< The timed somatic mutation queue 

    std::set<EpigeneticGenotypeId> death_enabled;  //!< Species having reached the death activation level

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
     *      children is pushed outside the tissue borders, the returned 
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
     * tissue and applies an epigenetic event to one of the children.
     * 
     * @param position is the position of the cell to be duplicate
     * @param final_id is the resulting genotype identifier of the cell
     * @return list of the affected cells. Whenever the specified position 
     *      is not in the tissue or does not correspond to a driver 
     *      mutation, the returned list is empty. If one of the two 
     *      children is pushed outside the tissue borders, the returned 
     *      list contains only the new cell in position `position`. 
     *      In the remaining case, the returned list contains two cells
     */
    EventAffectedCells simulate_duplication_epigenetic_event(const Position& position, const EpigeneticGenotypeId& final_id);

    /**
     * @brief Record tissue initial cells in log
     */
    void log_initial_cells();

    /**
     * @brief Initialize the valid simulation tissue direction vector
     * 
     * This method must be called every time the simulation tissue
     * changes.
     */
    void init_valid_directions();

    /**
     * @brief Randomly select a cell among those having a specified somatic genotype
     * 
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
    static const CellInTissue* choose_a_cell_in_somatic_species(const Tissue& tissue, const SomaticGenotypeId& genotype_id, std::mt19937_64& generator);

public:
    size_t death_activation_level;  //!< The minimum number of cells required to activate death

    /**
     * @brief The basic simulation constructor
     * 
     * @param random_seed is the simulation random seed
     */
    Simulation(int random_seed=0);

    /**
     * @brief A swap constructor
     * 
     * @param orig is the original simulation
     */
    Simulation(Simulation&& orig);

    /**
     * @brief A copy operator
     * 
     * @param orig is the original simulation
     * @return A reference of the updated object
     */
    Simulation& operator=(Simulation&& orig);

    /**
     * @brief Add a timed driver somatic mutation
     * 
     * @param src is the source driver somatic genotype
     * @param dst is the destination driver somatic genotype
     * @param time is the mutation timing
     * @return a reference to the updated simulation
     */
    Simulation& add_somatic_mutation(const SomaticGenotype& src, const SomaticGenotype& dst, const Time time);

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
     * This method simulates a tissue up to the next 
     * event. If the user provide a pointer to a 
     * plotter, then the simulation is also plotted
     * in a graphical window.
     * 
     * @tparam PLOT_WINDOW is the plotting window type
     * @param plotter is a tissue plotter pointer
     * @param logging is a flag to enable/disable logging
     * @return a reference to the updated simulation
     */
    template<typename PLOT_WINDOW>
    Simulation& run_up_to_next_event(UI::TissuePlotter<PLOT_WINDOW>* plotter, const bool logging = true);

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
     * @param logging is a flag to enable/disable logging
     * @return a reference to the updated simulation
     */
    template<typename PLOT_WINDOW>
    Simulation& run_up_to(const Time& final_time, UI::TissuePlotter<PLOT_WINDOW>* plotter, const bool logging = true);

    /**
     * @brief Simulate up to the next event
     * 
     * This method simulates a tissue up to the next 
     * event. The simulation is also plotted in a 
     * graphical window.
     * 
     * @tparam PLOT_WINDOW is the plotting window type
     * @param plotter is a tissue plotter
     * @param logging is a flag to enable/disable logging
     * @return a reference to the updated simulation
     */
    template<typename PLOT_WINDOW>
    Simulation& run_up_to_next_event(UI::TissuePlotter<PLOT_WINDOW>& plotter, const bool logging = true);

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
     * @param logging is a flag to enable/disable logging
     * @return a reference to the updated simulation
     */
    template<typename PLOT_WINDOW>
    Simulation& run_up_to(const Time& final_time, UI::TissuePlotter<PLOT_WINDOW>& plotter, const bool logging = true);

    /**
     * @brief Simulate up to the next event
     * 
     * This method simulates a tissue up to the next 
     * event. If the user provide a pointer to a 
     * plotter, then the simulation is also plotted
     * in a graphical window.
     * 
     * @param logging is a flag to enable/disable logging
     * @return a reference to the updated simulation
     */
    Simulation& run_up_to_next_event(const bool logging = true);

    /**
     * @brief Simulate a tissue up to a given time
     * 
     * This method simulates a tissue up to a given 
     * simulated time. If the user provide a pointer 
     * to a plotter, then the simulation is also plotted
     * in a graphical window.
     * 
     * @param final_time is the final simulation time
     * @param logging is a flag to enable/disable logging
     * @return a reference to the updated simulation
     */
    Simulation& run_up_to(const Time& final_time, const bool logging = true);

    /**
     * @brief Get the current simulation time
     * 
     * @return a constant reference to the simulation time
     */
    const Time& get_time() const;

    /**
     * @brief Get the simulation statistics
     * 
     * @return the simulation statistics
     */
    const TissueStatistics& get_statistics() const;

    /**
     * @brief Set a new simulation tissue
     * 
     * This method resets the simulation and sets a 
     * new simulation tissue.
     * 
     * @param name is the tissue name
     * @param sizes are the sizes of the tissue
     * @return a reference to the new tissue
     */
    Tissue& set_tissue(const std::string name, const std::vector<AxisSize> sizes);

    /**
     * @brief Get the simulation tissue
     * 
     * This method returns a reference to the simulation tissue.
     * A `std::runtime_error` object is throws if no tissue has been 
     * associated to the simulation yet.
     * 
     * @return A constant reference to the associated tissue
     */
    const Tissue& tissue() const;

    /**
     * @brief Get the simulated tissue
     * 
     * This method returns a reference to the simulated tissue.
     * A `std::runtime_error` object is throws if no tissue has been 
     * associated to the simulation yet.
     * 
     * @return A constant reference to the associated tissue
     */
    Tissue& tissue();

    /**
     * @brief Set the interval between snapshots
     * 
     * @param time_interval is the time interval between two snapshots
     */
    template<class Rep, class Period>
    void set_interval_between_snapshots(const std::chrono::duration<Rep,Period> time_interval);

    /**
     * @brief Reset a simulation
     */
    void reset();

    /**
     * @brief Save a simulation in an archive
     * 
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        archive & tissues
                & secs_between_snapshots
                & statistics
                & time
                & somatic_mutation_queue
                & death_enabled
                & death_activation_level;
    }

    /**
     * @brief Load a simulation from an archive
     * 
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the loaded simulation
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    static Simulation load(ARCHIVE& archive)
    {
        Simulation simulation;

        archive & simulation.tissues
                & simulation.secs_between_snapshots
                & simulation.statistics
                & simulation.time
                & simulation.somatic_mutation_queue
                & simulation.death_enabled
                & simulation.death_activation_level;

        simulation.init_valid_directions();

        simulation.last_snapshot_time = system_clock::now();

        return simulation;
    }

    ~Simulation();
};

/* Inline and template implementations */

inline const Time& Simulation::get_time() const
{
    return time;
}

inline const TissueStatistics& Simulation::get_statistics() const
{
    return statistics;
}

template<typename PLOT_WINDOW>
Simulation& Simulation::run_up_to_next_event(UI::TissuePlotter<PLOT_WINDOW>* plotter, const bool logging)
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
        if (logging) {
            logger->record(event.type, cell, time);
        }

        // if death has not been enabled yet
        const auto species_id = cell.get_genotype_id();
        if (death_enabled.count(species_id)==0) {
            Species& species = tissue().get_species(species_id);
            
            // and the death activation level has been reached
            if (species.num_of_cells()>=death_activation_level) {

                // enable death
                death_enabled.insert(species_id);
            }
        }
    }

    statistics.record_event(event, time);
    for (const auto& cell: affected.lost_cells) {
        statistics.record_lost(cell.get_genotype_id(), time);
    }

    if (plotter != nullptr && !plotter->closed()) {
        plotter->plot(statistics);
    }

    return *this;
}

template<typename PLOT_WINDOW>
Simulation& Simulation::run_up_to(const Time& final_time, UI::TissuePlotter<PLOT_WINDOW>* plotter, const bool logging)
{
    // the tissue() call checks whether a tissue has been 
    // associated to the simulation and, if this is not the 
    // case, it throws an std::runtime_error 
    (void)tissue();

    if (logger == nullptr) {
        logger = new BinaryLogger();
    }

    // if we are at the beginning of the computation, log the initial cells
    if (time == 0 && logging) {
        log_initial_cells();
    }

    while ((plotter == nullptr || !plotter->closed()) && 
           tissue().num_of_mutated_cells()>0 && time < final_time) {
        using namespace std::chrono;

        run_up_to_next_event(plotter, logging);

        const auto curr_time = system_clock::now();
        const auto from_last_snapshot = duration_cast<seconds>(curr_time-last_snapshot_time);

        if (secs_between_snapshots>0 &&from_last_snapshot.count()>=secs_between_snapshots) {
            last_snapshot_time = curr_time;
            if (logging) {
                logger->snapshot(*this);
            }
        }
    }

    if (logging) {
        logger->snapshot(*this);
        logger->flush_archives();
    }

    return *this;
}

template<typename PLOT_WINDOW>
inline Simulation& Simulation::run_up_to_next_event(UI::TissuePlotter<PLOT_WINDOW>& plotter, const bool logging)
{
    return run_up_to_next_event(&plotter, logging);
}

template<typename PLOT_WINDOW>
inline Simulation& Simulation::run_up_to(const Time& final_time, UI::TissuePlotter<PLOT_WINDOW>& plotter, const bool logging)
{
    return run_up_to(final_time, &plotter, logging);
}

inline Simulation& Simulation::run_up_to_next_event(const bool logging)
{
    return run_up_to_next_event<UI::Plot2DWindow>(nullptr, logging);
}

inline Simulation& Simulation::run_up_to(const Time& final_time,const bool logging)
{
    return run_up_to<UI::Plot2DWindow>(final_time, nullptr, logging);
}

template<class Rep, class Period>
inline void Simulation::set_interval_between_snapshots(const std::chrono::duration<Rep,Period> time_interval)
{
    using namespace std::chrono;

    secs_between_snapshots = duration_cast<seconds>(time_interval).count();
}

}   // Races

#endif // __RACES_SIMULATOR__