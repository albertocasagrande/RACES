/**
 * @file simulation_wrapper.hpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Define the Python wrapper class and functions for `Races::Simulation`
 * @version 0.4
 * @date 2023-08-05
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

#ifndef __RACES_PYTHON_SIMULATION_WRAPPER__
#define __RACES_PYTHON_SIMULATION_WRAPPER__

#include <memory>
#include <shared_mutex>

#include <boost/python.hpp>

#include "driver_genotype.hpp"
#include "simulation.hpp"

using namespace boost::python;

/**
 * @brief A wrapper structure for `Races::Simulation`
 * 
 * The `Races::Simulation` class cannot be used directly 
 * because it cannot be copied by design. This class 
 * wraps `Races::Simulation` by maintaining a shared 
 * pointer to a `Races::Simulation` object and 
 * synchronizing all the accesses to it by using a 
 * mutex.
 */
class SimulationWrapper 
{
    struct _SimulationWrapper
    {
        using Simulation = Races::Drivers::Simulation::Simulation;

        Simulation simulation;         //!< the c++ simulation object

        std::shared_mutex s_mutex;     //!< the simulation mutex

        explicit _SimulationWrapper(int random_seed);
    };

    std::shared_ptr<SimulationWrapper::_SimulationWrapper> obj_ptr; //!< the pointer to the real c++ wrapper

public:
    /**
     * @brief The basic simulation constructor
     * 
     * @param random_seed is the simulation random seed
     */
    explicit SimulationWrapper(int random_seed=0);

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
     * @param quiet is a flag to enable/disable the progress bar
     * @param plot is a flag to enable/disable plotting
     * @return a reference to the updated simulation
     */
    void run_up_to(const Races::Time& final_time, const bool logging = true,
                   const bool quiet = false, const bool plot = false);

    /**
     * @brief Get the current simulation time
     * 
     * @return a constant reference to the simulation time
     */
    const Races::Time& get_time() const;

    /**
     * @brief Add a timed driver genomic mutation
     * 
     * @param src is the source driver genomic genotype
     * @param dst is the destination driver genomic genotype
     * @param time is the mutation timing
     * @return a reference to the updated simulation
     */
    void add_genomic_mutation(const Races::Drivers::Genotype& src, const Races::Drivers::Genotype& dst,
                              const Races::Time time);

    /**
     * @brief Add a new species to the tissue
     * 
     * @param genotype is the driver genotype of the new species
     * @return a reference to the updated object
     */
    void add_species(const Races::Drivers::Genotype& genotype);

    /**
     * @brief Add a cell to the simulated tissue
     * 
     * @param genotype is the driver genotype of the new cell
     * @param methylation_signature is the methylation signature of the new cell
     * @param position is the initial position in the tissue
     * @return a reference to the updated object
     */
    void add_cell(const Races::Drivers::Genotype& genotype, const std::string& methylation_signature, 
                  boost::python::list const& position);

    /**
     * @brief Set a new simulation tissue
     * 
     * This method resets the simulation and sets a 
     * new simulation tissue.
     * 
     * @param name is the tissue name
     * @param sizes_list are the sizes of the tissue
     * @return a reference to the updated object
     */
    void set_tissue(const std::string& name, boost::python::list const& sizes_list);

    /**
     * @brief Get the death activation level
     * 
     * @return the death activation level
     */
    inline size_t get_death_activation_level() const
    {
        return obj_ptr->simulation.death_activation_level;
    }

    /**
     * @brief Set the death activation level
     * 
     * @param death_activation_level is the new death activation level 
     */
    inline void set_death_activation_level(const size_t& death_activation_level)
    {
        obj_ptr->simulation.death_activation_level = death_activation_level;
    }

    /**
     * @brief Create a new `SimulationWrapper` object
     * 
     * @param minutes_between_snapshots is the number of minutes between snapshots
     * @param random_seed is the simulation random seed
     * @return a shared pointer to the newly created 
     *        `SimulationWrapper` object
     */
    inline static std::shared_ptr<SimulationWrapper> create(const unsigned int minutes_between_snapshots=5, int random_seed=0)
    {
        auto wrapper = std::make_shared<SimulationWrapper>(random_seed);

        const std::chrono::minutes snapshot_interval(minutes_between_snapshots);

        wrapper->obj_ptr->simulation.set_interval_between_snapshots(snapshot_interval);

        return wrapper;
    }

    /**
     * @brief Simulate a tissue up to a given time
     * 
     * This method simulates a tissue up to a given 
     * simulated time. If the user provide a pointer 
     * to a plotter, then the simulation is also plotted
     * in a graphical window.
     * 
     * @param wrapper is the simulation wrapper on which the 
     *                  computation is performed
     * @param final_time is the final simulation time
     * @param logging is a flag to enable/disable logging
     * @param quiet is a flag to enable/disable the progress bar
     * @param plot is a flag to enable/disable plotting
     */
    inline static void static_run_up_to(SimulationWrapper *wrapper, const Races::Time& final_time, 
                                 const bool logging = true, const bool quiet = false,
                                 const bool plot = false)
    {
        wrapper->run_up_to(final_time, logging, quiet, plot);
    }
};

#endif // __RACES_PYTHON_SIMULATION_WRAPPER__