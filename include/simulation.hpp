/**
 * @file simulation.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines a tumor evolution simulation
 * @version 0.30
 * @date 2023-11-03
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
#include <chrono>

#include "tissue.hpp"
#include "binary_logger.hpp"
#include "plot_2D.hpp"
#include "tissue_plotter.hpp"
#include "statistics.hpp"
#include "timed_event.hpp"
#include "tissue_sample.hpp"

#include "phylogenetic_forest.hpp"

#include "progress_bar.hpp"

namespace Races 
{

namespace Drivers
{

/**
 * @brief The namespace of driver simulations
 */
namespace Simulation
{

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

    using TimedEventQueue = std::priority_queue<TimedEvent,
                                                std::vector<TimedEvent>,
                                                std::greater<TimedEvent>>;
    using system_clock = std::chrono::system_clock;

    std::vector<Tissue> tissues;     //!< Simulated tissues
    std::map<std::string, GenotypeId> genotype_name2id; //!< A map associating the genotype name to its id

    BinaryLogger logger;             //!< Event logger

    std::vector<Direction> valid_directions;   //!< valid simulation tissue directions

    system_clock::time_point last_snapshot_time;    //!< Time of the last snapshot
    long secs_between_snapshots;    //!< The number of seconds between two snapshots

    TissueStatistics statistics;     //!< The tissue simulation statistics

    Time time;                       //!< Simulation time
    std::mt19937_64 random_gen;      //!< Pseudo-random generator

    TimedEventQueue timed_event_queue;   //!< The timed event queue 

    std::set<SpeciesId> death_enabled;  //!< Species having reached the death activation level

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
     * If the cell in the provided position has non-genotype, 
     * then nothing is done.  
     * 
     * @param position is the position of the cell that will mutate
     * @param final_id is the resulting genotype identifier of the cell
     * @return list of the affected cells, i.e., a list containing 
     *      the status of the mutated cell at most
     */
    EventAffectedCells simulate_mutation(const Position& position, const SpeciesId& final_id);

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
    EventAffectedCells simulate_duplication_epigenetic_event(const Position& position, const SpeciesId& final_id);

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
     * @brief Handle a driver mutation
     * 
     * This method tries to handle a driver mutation during the cell event selection and 
     * it succeeds if and only if there exists at least one cell in the origin.
     * 
     * @param timed_genotype_mutation is the timed driver mutation to be applied
     * @param candidate_event is the current candidate cell event
     * @return `true` if and only if there exists at least one cell in the origin and 
     *          the candidate cell event has been updated
     */
    bool handle_timed_genotype_mutation(const TimedEvent& timed_genotype_mutation, CellEvent& candidate_event);

    /**
     * @brief Apply a rate update
     * 
     * @param time_rate_update is the timed rate update to be applied
     */
    void handle_timed_rate_update(const TimedEvent& timed_rate_update);

    /**
     * @brief Sample the tissue
     * 
     * This method also update next candidate event.
     * 
     * @param time_sampling is the timed sampling to be applied
     * @param candidate_event is the current candidate cell event
     */
    void handle_timed_sampling(const TimedEvent& time_sampling, CellEvent& candidate_event);

    /**
     * @brief Handle time event queue during cell event selection
     * 
     * This method handles the time event queue during the cell event selection
     * by extracting the timed events that occurs before the candidate cell 
     * event and applying them. If necessary (for instance when the next 
     * time event is a driver mutation), the method also updates the
     * candidate event.
     * 
     * @param candidate_event is a candidate cell event
     */
    void handle_timed_event_queue(CellEvent& candidate_event);

    /**
     * @brief Select a cell event among those due to cell liveness
     * 
     * @return a cell event among those due to cell liveness
     */
    CellEvent select_next_cell_event();

public:
    /**
     * @brief Simulation test
     */
    struct BasicTest
    {
        /**
         * @brief Test a simulation
         * 
         * @param simulation is the considered simulation
         * @return a Boolean value
         */
        virtual bool operator()(const Simulation& simulation) const = 0;

        /**
         * @brief Return the percentage of the completed simulation
         * 
         * @param simulation is the considered simulation
         * @return the percentage of the completed simulation
         */
        virtual uint8_t percentage(const Simulation& simulation) const = 0;

        virtual ~BasicTest()
        {}
    };

    size_t death_activation_level;  //!< The minimum number of cells required to activate death

    bool duplicate_internal_cells; //!< A flag to enable/disable duplication in internal cells
    bool storage_enabled;          //!< A flag to enable/disable storage

    /**
     * @brief The basic simulation constructor
     * 
     * @param random_seed is the simulation random seed
     */
    explicit Simulation(int random_seed=0);

    /**
     * @brief A simulation constructor
     * 
     * @param log_directory is the simulation log directory
     * @param random_seed is the simulation random seed
     */
    explicit Simulation(const std::filesystem::path& log_directory, int random_seed=0);

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
     * @brief Schedule a timed genotype mutation
     * 
     * @param src is the source genotype
     * @param dst is the destination genotype
     * @param time is the mutation timing
     * @return a reference to the updated simulation
     */
    Simulation& schedule_genotype_mutation(const GenotypeProperties& src, 
                                           const GenotypeProperties& dst, const Time time);

    /**
     * @brief Schedule a timed genotype mutation
     * 
     * @param src is the source genotype name
     * @param dst is the destination genotype name
     * @param time is the mutation timing
     * @return a reference to the updated simulation
     */
    Simulation& schedule_genotype_mutation(const std::string& src, const std::string& dst, const Time time);

    /**
     * @brief Schedule a timed event
     * 
     * @param timed_event is the timed event to schedule 
     * @return a reference to the updated simulation
     */
    Simulation& schedule_timed_event(const TimedEvent& timed_event);

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
     * @brief  Simulate a genotype mutation on a position
     * 
     * This method simulates both the duplication of the cell in the 
     * specified position and the birth of a cells of a different 
     * genotype preserving the epigenetic status of the original cell. 
     * 
     * @param position is the position in which the 
     * @param dst_genotype_name is the name of the mutated cell genotype
     * @return a reference to the updated simulation
     */
    Simulation& simulate_genotype_mutation(const PositionInTissue& position,
                                           const std::string& dst_genotype_name);

    /**
     * @brief  Simulate a genotype mutation on a position
     * 
     * This method simulates both the duplication of the cell in the 
     * specified position and the birth of a cells of a different 
     * genotype preserving the epigenetic status of the original cell. 
     * 
     * @param position is the position in which the 
     * @param dst_genotype_id is the identifier of the mutated cell genotype
     * @return a reference to the updated simulation
     */
    Simulation& simulate_genotype_mutation(const PositionInTissue& position,
                                           const GenotypeId& dst_genotype_id);

    /**
     * @brief Simulate an event
     * 
     * This method simulates an event on the tissue. If `plotter` is not
     * `nullptr`, then the simulation is also plotted in a graphical window.
     * 
     * @tparam PLOT_WINDOW is the plotting window type
     * @param event is the event to be simulated
     * @param plotter is a tissue plotter pointer
     * @return a reference to the updated simulation
     */
    template<typename PLOT_WINDOW>
    Simulation& simulate(const CellEvent& event, UI::TissuePlotter<PLOT_WINDOW>* plotter);

    /**
     * @brief Simulate an event
     * 
     * This method simulates an event on the tissue. If `plotter` is not
     * `nullptr`, then the simulation is also plotted in a graphical window.
     * 
     * @tparam PLOT_WINDOW is the plotting window type
     * @param event is the event to be simulated
     * @param plotter is a tissue plotter pointer
     * @return a reference to the updated simulation
     */
    template<typename PLOT_WINDOW>
    inline Simulation& simulate(const CellEvent& event, UI::TissuePlotter<PLOT_WINDOW>& plotter)
    {
        return simulate(event, &plotter);
    }

    /**
     * @brief Simulate an event
     * 
     * @param event is the event to be simulated
     * @return a reference to the updated simulation
     */
    inline Simulation& simulate(const CellEvent& event)
    {
        return simulate<UI::Plot2DWindow>(event, nullptr);
    }

    /**
     * @brief Simulate up to the next event
     * 
     * This method simulates a tissue up to the next event. If `plotter` is not
     * `nullptr`, then the simulation is also plotted in a graphical window.
     * 
     * @tparam PLOT_WINDOW is the plotting window type
     * @param plotter is a tissue plotter pointer
     * @return a reference to the updated simulation
     */
    template<typename PLOT_WINDOW>
    Simulation& run_up_to_next_event(UI::TissuePlotter<PLOT_WINDOW>* plotter);

    /**
     * @brief Simulate up to the next event
     * 
     * This method simulates a tissue up to the next event. The simulation is also 
     * plotted in a graphical window.
     * 
     * @tparam PLOT_WINDOW is the plotting window type
     * @param plotter is a tissue plotter
     * @return a reference to the updated simulation
     */
    template<typename PLOT_WINDOW>
    inline Simulation& run_up_to_next_event(UI::TissuePlotter<PLOT_WINDOW>& plotter)
    {
        return run_up_to_next_event(&plotter);
    }

    /**
     * @brief Simulate up to the next event
     * 
     * This method simulates a tissue up to the next 
     * event.
     * 
     * @return a reference to the updated simulation
     */
    inline Simulation& run_up_to_next_event()
    {
        return run_up_to_next_event<UI::Plot2DWindow>(nullptr);
    }

    /**
     * @brief Simulate a tissue
     * 
     * This method simulates a tissue until a test on the simulation 
     * object holds. If the user provide a pointer to a plotter, then 
     * the simulation is also plotted in a graphical window.
     * 
     * @tparam SIMULATION_TEST is a class whose objects can test whether 
     *              the simulation is concluded
     * @tparam PLOT_WINDOW is the plotting window type
     * @tparam INDICATOR_TYPE is the progress indicator type
     * @param done is a test verifying whether the simulation has come to an end
     * @param plotter is a tissue plotter pointer
     * @param indicator is the progress indicator pointer
     * @return a reference to the updated simulation
     */
    template<class SIMULATION_TEST, typename PLOT_WINDOW, typename INDICATOR_TYPE,
             std::enable_if_t<std::is_base_of_v<BasicTest, SIMULATION_TEST> &&
                              std::is_base_of_v<UI::Plot2DWindow, PLOT_WINDOW>, bool> = true>
    Simulation& run(SIMULATION_TEST& done, UI::TissuePlotter<PLOT_WINDOW>* plotter,
                    INDICATOR_TYPE* indicator)
    {
        // the tissue() call checks whether a tissue has been 
        // associated to the simulation and, if this is not the 
        // case, it throws an std::runtime_error 
        (void)tissue();

        if (storage_enabled) {
            // if we are at the beginning of the computation, 
            // log the initial cells
            if (time == 0) {
                log_initial_cells(); 
            }
        }

        while ((plotter == nullptr || !plotter->closed()) && !done(*this) 
            && tissue().num_of_mutated_cells()>0) {

            run_up_to_next_event(plotter);

            if (storage_enabled) {
                snapshot_on_time(indicator);
            }


            if (indicator != nullptr) {
                indicator->set_progress(done.percentage(*this),
                                        "Cells: " + std::to_string(tissue().num_of_mutated_cells()));
            }
        }

        if (storage_enabled) {
            if (indicator != nullptr) {
                indicator->set_message("Saving snapshot");
            }
            logger.snapshot(*this);
            logger.flush_archives();
        }

        if (indicator != nullptr) {
            indicator->set_message("Cells: " + std::to_string(tissue().num_of_mutated_cells()));
        }

        return *this;
    }

    /**
     * @brief Simulate a tissue
     * 
     * This method simulates a tissue until a test on the simulation 
     * object holds. The simulation is also plotted in a graphical window.
     *
     * @tparam SIMULATION_TEST is a class whose objects can test whether 
     *              the simulation is concluded
     * @tparam PLOT_WINDOW is the plotting window type
     * @tparam INDICATOR_TYPE is the progress indicator type
     * @param done is a test verifying whether the simulation has come to an end
     * @param plotter is a tissue plotter
     * @param indicator is the progress indicator pointer
     * @return a reference to the updated simulation
     */
    template<class SIMULATION_TEST, typename PLOT_WINDOW, typename INDICATOR_TYPE,
             std::enable_if_t<std::is_base_of_v<BasicTest, SIMULATION_TEST> &&
                              std::is_base_of_v<UI::Plot2DWindow, PLOT_WINDOW>, bool> = true>
    inline Simulation& run(SIMULATION_TEST& done, UI::TissuePlotter<PLOT_WINDOW>& plotter,
                           INDICATOR_TYPE& indicator)
    {
        return run(done, &plotter, &indicator);
    }

    /**
     * @brief Simulate a tissue
     * 
     * This method simulates a tissue until a test on the simulation holds.
     * 
     * @tparam SIMULATION_TEST is a class whose objects can test whether 
     *              the simulation is concluded
     * @tparam INDICATOR_TYPE is the progress indicator type
     * @param done is a test verifying whether the simulation has come to an end
     * @param indicator is the progress indicator
     * @return a reference to the updated simulation
     */
    template<class SIMULATION_TEST, typename INDICATOR_TYPE,
             std::enable_if_t<std::is_base_of_v<BasicTest, SIMULATION_TEST>, bool> = true>
    inline Simulation& run(SIMULATION_TEST& done, INDICATOR_TYPE& indicator)
    {
        return run<SIMULATION_TEST, UI::Plot2DWindow, INDICATOR_TYPE>(done, nullptr, &indicator);
    }

    /**
     * @brief Simulate a tissue
     * 
     * This method simulates a tissue until a test on the simulation 
     * object holds. The simulation is also plotted in a graphical window.
     *
     * @tparam SIMULATION_TEST is a class whose objects can test whether 
     *              the simulation is concluded
     * @tparam PLOT_WINDOW is the plotting window type
     * @param done is a test verifying whether the simulation has come to an end
     * @param plotter is a tissue plotter
     * @return a reference to the updated simulation
     */
    template<class SIMULATION_TEST, typename PLOT_WINDOW,
             std::enable_if_t<std::is_base_of_v<BasicTest, SIMULATION_TEST> &&
                              std::is_base_of_v<UI::Plot2DWindow, PLOT_WINDOW>, bool> = true>
    inline Simulation& run(SIMULATION_TEST& done, UI::TissuePlotter<PLOT_WINDOW>& plotter)
    {
        return run<SIMULATION_TEST, PLOT_WINDOW, UI::ProgressBar>(done, &plotter, nullptr);
    }

    /**
     * @brief Simulate a tissue
     * 
     * This method simulates a tissue until a test on the simulation 
     * object holds.
     *
     * @tparam SIMULATION_TEST is a class whose objects can test whether 
     *              the simulation is concluded
     * @param done is a test verifying whether the simulation has come to an end
     * @return a reference to the updated simulation
     */
    template<class SIMULATION_TEST,
             std::enable_if_t<std::is_base_of_v<BasicTest, SIMULATION_TEST>, bool> = true>
    inline Simulation& run(SIMULATION_TEST& done)
    {
        return run<SIMULATION_TEST, UI::Plot2DWindow, UI::ProgressBar>(done, nullptr, nullptr);
    }

    /**
     * @brief Randomly select a cell among those having a specified genotype
     * 
     * @param genotype_id is the identifier of the genotype that must contain 
     *          the selected cell
     * @return whenever the set of cells having `genotype_id` as genotype 
     *         identifier is not empty, a randomly selected cell in it. 
     *         Otherwise, if the set is empty, a domain error is thrown.
     */
    const CellInTissue& choose_cell_in(const GenotypeId& genotype_id);

    /**
     * @brief Randomly select a cell among those having a specified genotype
     * 
     * @param genotype_name is the name of the genotype that must contain 
     *          the selected cell
     * @return whenever the set of cells having `genotype_id` as genotype 
     *         identifier is not empty, a randomly selected cell in it. 
     *         Otherwise, if the set is empty, a domain error is thrown.
     */
    const CellInTissue& choose_cell_in(const std::string& genotype_name);

    /**
     * @brief Randomly select in a tissue rectangle a cell having a specified genotype
     * 
     * @param genotype_id is the identifier of the genotype that must contain 
     *          the selected cell
     * @param rectangle is the tissue rectangle in which the cell must be selected 
     * @return whenever the set of tissue cells that have `genotype_id` as genotype 
     *         and lay in one the positions specified by `rectangle` is not empty, 
     *         a randomly selected cell in it. Otherwise, if the set is empty, a 
     *         domain error is thrown.
     */
    const CellInTissue& choose_cell_in(const GenotypeId& genotype_id,
                                       const RectangleSet& rectangle);

    /**
     * @brief Randomly select in a tissue rectangle a cell having a specified genotype
     * 
     * @param genotype_name is the name of the genotype that must contain 
     *          the selected cell
     * @param rectangle is the tissue rectangle in which the cell must be selected 
     * @return whenever the set of tissue cells that have `genotype_id` as genotype 
     *         and lay in one the positions specified by `rectangle` is not empty, 
     *         a randomly selected cell in it. Otherwise, if the set is empty, a 
     *         domain error is thrown.
     */
    const CellInTissue& choose_cell_in(const std::string& genotype_name,
                                       const RectangleSet& rectangle);

    /**
     * @brief Get the current simulation time
     * 
     * @return a constant reference to the simulation time
     */
    inline const Time& get_time() const
    {
        return time;
    }

    /**
     * @brief Get the simulation statistics
     * 
     * @return the simulation statistics
     */
    inline const TissueStatistics& get_statistics() const
    {
        return statistics;
    }

    /**
     * @brief Add a new species to the tissue
     * 
     * @param genotype_properties is the genotype properties of the new species
     * @return a reference to the updated object
     */
    Simulation& add_species(const GenotypeProperties& genotype_properties);

    /**
     * @brief Add a cell to the simulated tissue
     * 
     * @param species_properties is the species properties of the new cell
     * @param position is the initial position in the tissue
     * @return a reference to the updated object
     */
    inline Simulation& add_cell(const SpeciesProperties& species_properties,
                                const PositionInTissue& position)
    {
        return add_cell(species_properties.get_id(), position);
    }

    /**
     * @brief Add a cell to the simulated tissue
     * 
     * @param species_id is the species identifier of the new cell
     * @param position is the initial position in the tissue
     * @return a reference to the updated object
     */
    inline Simulation& add_cell(const SpeciesId& species_id, const PositionInTissue& position)
    {
        tissue().add_cell(species_id, position);

        return *this;
    }

    /**
     * @brief Set a new simulation tissue
     * 
     * This method resets the simulation and sets a 
     * new simulation tissue.
     * 
     * @param name is the tissue name
     * @param sizes are the sizes of the tissue
     * @return a reference to the updated object
     */
    Simulation& set_tissue(const std::string& name, const std::vector<AxisSize>& sizes);

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
    inline void set_interval_between_snapshots(const std::chrono::duration<Rep,Period> time_interval)
    {
        using namespace std::chrono;

        secs_between_snapshots = duration_cast<seconds>(time_interval).count();
    }

    /**
     * @brief Reset a simulation
     */
    void reset();

    /**
     * @brief Update the log directory
     * 
     * @param log_directory is the new simulation log directory
     */
    inline void rename_log_directory(const std::filesystem::path& log_directory)
    {
        logger.rename_directory(log_directory);
    }

    /**
     * @brief Inject random generator seed
     * 
     * @param random_seed is the simulation random seed
     */
    inline Simulation& random_generator_seed(int random_seed)
    {
        random_gen.seed(random_seed);

        return *this;
    }
    
    /**
     * @brief Performs a simulation snapshot
     * 
     * This method performs a simulation snapshot.
     * 
     * @tparam INDICATOR is the type of progress bar
     * @param indicator is the progress bar
     */
    template<typename INDICATOR>
    void make_snapshot(INDICATOR *indicator);

    /**
     * @brief Performs a simulation snapshot
     * 
     * This method performs a simulation snapshot if the time 
     * elapsed from the last snapshot is greater than that set 
     * in `secs_between_snapshots`.
     * 
     * @tparam INDICATOR is the type of progress bar
     * @param indicator is the progress bar
     */
    template<typename INDICATOR>
    void snapshot_on_time(INDICATOR *indicator);

    /**
     * @brief Sample the simulation tissue
     * 
     * This method samples the simulation tissue in a non-destructive way.
     * 
     * @param rectangle is the rectangle to sample
     * @return the sample of the tissue in `rectangle`
     */
    TissueSample sample_tissue(const RectangleSet& rectangle) const;

    /**
     * @brief Sample the simulation tissue
     * 
     * This method samples the simulation tissue. The sampling 
     * can either preserve the tissue or remove the sampled cells 
     * having driver mutations according to the value of a 
     * Boolean parameter (`preserve_tissue`).
     * 
     * @param rectangle is the rectangle to sample
     * @param preserve_tissue is a Boolean flag to enable/disable 
     *          sample removal from simulation tissue
     * @return the sample of the tissue in `rectangle`
     */
    TissueSample sample_tissue(const RectangleSet& rectangle, const bool& preserve_tissue=false);

    /**
     * @brief Get the simulation logger
     * 
     * @return a constant reference to the simulation logger
     */
    inline const BinaryLogger& get_logger() const
    {
        return logger;
    }

    /**
     * @brief Find a genotype id by name
     * 
     * @param genotype_name is the name of the genotype whose id is aimed
     * @return the genotype id of the genotype called `genotype_name
     */
    GenotypeId find_genotype_id(const std::string& genotype_name) const;

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
                & genotype_name2id
                & logger
                & secs_between_snapshots
                & statistics
                & time
                & timed_event_queue
                & death_enabled
                & death_activation_level
                & duplicate_internal_cells
                & storage_enabled
                & Cell::counter;
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
                & simulation.genotype_name2id
                & simulation.logger
                & simulation.secs_between_snapshots
                & simulation.statistics
                & simulation.time
                & simulation.timed_event_queue
                & simulation.death_enabled
                & simulation.death_activation_level
                & simulation.duplicate_internal_cells
                & simulation.storage_enabled
                & Cell::counter;

        simulation.init_valid_directions();

        simulation.last_snapshot_time = system_clock::now();

        return simulation;
    }

    ~Simulation();
};

/* Template implementations */

template<typename INDICATOR>
void Simulation::snapshot_on_time(INDICATOR *indicator)
{
    using namespace std::chrono;

    const auto from_last_snapshot = duration_cast<seconds>(system_clock::now()-last_snapshot_time);

    if (secs_between_snapshots>0 && from_last_snapshot.count()>=secs_between_snapshots) {
        make_snapshot(indicator);
    }
}

template<typename INDICATOR>
void Simulation::make_snapshot(INDICATOR *indicator)
{
    last_snapshot_time = system_clock::now();
    if (indicator != nullptr) {
        indicator->set_message("Saving snapshot");
    }
    if (storage_enabled) {
        logger.snapshot(*this);
    }
}

template<typename PLOT_WINDOW>
inline
Simulation& Simulation::run_up_to_next_event(UI::TissuePlotter<PLOT_WINDOW>* plotter)
{
    CellEvent event = select_next_event();

    return simulate(event, plotter);
}

template<typename PLOT_WINDOW>
Simulation& Simulation::simulate(const CellEvent& event, UI::TissuePlotter<PLOT_WINDOW>* plotter)
{
    time += event.delay;

    EventAffectedCells affected;

    switch(event.type) {
        case CellEventType::DEATH:
            affected = simulate_death(event.position);
            break;
        case CellEventType::DUPLICATION:
            affected = simulate_duplication(event.position);
            break;
        case CellEventType::DUPLICATION_AND_EPIGENETIC_EVENT:
            affected = simulate_duplication_epigenetic_event(event.position, event.final_species);
            break;
        case CellEventType::DRIVER_GENETIC_MUTATION:
            affected = simulate_mutation(event.position, event.final_species);
            break;
        default:
            throw std::runtime_error("Unhandled event type");
    }

    for (const auto& cell : affected.new_cells) {
        if (storage_enabled) {
            logger.record(event.type, cell, time);
        }

        // if death has not been enabled yet
        const auto species_id = cell.get_species_id();
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
        statistics.record_lost(cell.get_species_id(), time);
    }

    if (plotter != nullptr && !plotter->closed()) {
        plotter->plot(statistics);
    }

    return *this;
}

}   // Simulation

}   // Drivers

}   // Races

#endif // __RACES_SIMULATOR__