/**
 * @file statistics.hpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Define simulation statistics
 * @version 0.1
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

#ifndef __RACES_STATISTICS__
#define __RACES_STATISTICS__

#include <list>

#include "species.hpp"
#include "cell_event.hpp"

namespace Races {

class TissueStatistics;

/**
 * @brief A class for simulation statistics about one species
 */
class SpeciesStatistics 
{
    /**
     * @brief The empty constructor
     * 
     * @param num_of_cells is the number of cells in the species
     */
    SpeciesStatistics();   
public:
    Time rise_time;         //!< the time of the first appearance 
    Time extinction_time;   //!< the time of the last appearance

    size_t total_cells;     //!< the number of cells that have been in the species
    size_t curr_cells;      //!< the number of cells currently in the species
    size_t killed_cells;    //!< the number of species cells that were killed
    size_t lost_cells;      //!< the number of cells that have overcome the tissue border

    /**
     * @brief A constructor
     * 
     * @param num_of_cells is the number of cells in the species
     */
    SpeciesStatistics(const size_t& num_of_cells);

    friend class TissueStatistics;
};

/**
 * @brief A class to collect simulation statistics about a tissue
 */
class TissueStatistics
{
    using time_point = std::chrono::steady_clock::time_point;

    std::map<DriverGenotypeId, SpeciesStatistics> s_statistics;  //!< a map from species id to statistics
    
    std::list<Time> sim_times;             //!< the simulated times of the last recorded events
    std::list<time_point> real_times;      //!< the recording times of the last events 

    const size_t max_stored_times;         //!< the maximum number of times to store

    size_t total_events;                   //!< number of recorded events
    time_point first_event_time;           //!< first recoded event time
public:
    /**
     * @brief A constructor
     * 
     * @param tissue is the tissue whose statistics are collected
     */
    TissueStatistics(const Tissue& tissue);

    /**
     * @brief Get the statistics of a species
     * 
     * @param species is the species whose statistics are aimed
     * @return a non-constant reference to the statistics of `species` 
     */
    SpeciesStatistics& operator[](const Species& species);

    /**
     * @brief Get the statistics of a species
     * 
     * @param species_id is the identifier of the species whose statistics are aimed
     * @return a non-constant reference to the statistics of species having 
     *      `species_id` as identifier 
     */
    SpeciesStatistics& operator[](const DriverGenotypeId& species_id);

    /**
     * @brief Get the statistics of a species
     * 
     * @param species is the species whose statistics are aimed
     * @return a non-constant reference to the statistics of `species` 
     */
    const SpeciesStatistics& at(const Species& species) const;

    /**
     * @brief Get the statistics of a species
     * 
     * @param species_id is the identifier of the species whose statistics are aimed
     * @return a constant reference to the statistics of species having 
     *      `species_id` as identifier 
     */
    const SpeciesStatistics& at(const DriverGenotypeId& species_id) const;

    /**
     * @brief Record a death
     * 
     * @param genotype_id is the genotype id of the dying cell
     * @param time is the death time
     */
    void record_death(const DriverGenotypeId& genotype_id, const Time &time);

    /**
     * @brief Record a lost cell
     * 
     * A cell is lost whenever it is pushed outside the tissue border. 
     * This method record in the statistics that a cell is lost.
     * 
     * @param genotype_id is the genotype id of the lost cell
     * @param time is the time in which the lost cell has been pushed 
     *       outside the tissue border
     */
    void record_lost(const DriverGenotypeId& genotype_id, const Time &time);

    /**
     * @brief Record a cell duplication
     * 
     * @param genotype_id is the driver genotype id of the duplicating cell 
     * @param lost_a_cell is a flag to signal whenever a cell has been lost 
     *     during the event.
     */
    void record_duplication(const DriverGenotypeId& genotype_id, const bool& lost_a_cell);

    /**
     * @brief Record a cell duplication and an epigenetic event
     * 
     * @param genotype_id is the driver genotype id of the duplicating cell 
     * @param epigenetic_genotype is the driver genotype id of the mutated cell
     * @param time is the duplication time
     * @param lost_a_cell is a flag to signal whenever a cell has been lost 
     *     during the event. The lost cell is assumed to always be the 
     *     non-mutated one
     */
    void record_duplication_and_epigenetic_event(const DriverGenotypeId& genotype_id, const DriverGenotypeId& epigenetic_genotype, const Time &time, const bool& lost_a_cell);

    /**
     * @brief Record the last event
     * 
     * This method records a new event and, whenever the 
     * queues of the recorded times overlaps the maximum size, 
     * deletes the oldest recorded times.
     * 
     * @param event is the event to record
     * @param time is the event time
     * @param lost_a_cell is a flag to signal whenever a cell has been lost 
     *     during the event
     */
    void record_event(const CellEvent& event, const Time &time, const bool lost_a_cell);

    /**
     * @brief Get the number of recorded times
     * 
     * @return the number of recorded times
     */
    const size_t& get_recorded_time_number() const;

    /**
     * @brief Get the elapsed time
     * 
     * @return the elapsed time
     */
    std::chrono::steady_clock::duration get_elapsed_time() const;

    /**
     * @brief Get the simulated time
     * 
     * @return the simulated time
     */
    Time get_simulated_time() const;

    /**
     * @brief Get the number of the last recorded events over time
     * 
     * @return the number of the last recorded events over time
     */
    template<class toDuration>
    double get_last_recorded_events_over_time() const;
};

/* Inline implementations */

inline SpeciesStatistics& TissueStatistics::operator[](const Species& species)
{
    return s_statistics.at(species.get_id());
}

inline SpeciesStatistics& TissueStatistics::operator[](const DriverGenotypeId& species_id)
{
    return s_statistics.at(species_id);
}

inline const SpeciesStatistics& TissueStatistics::at(const Species& species) const
{
    return s_statistics.at(species.get_id());
}

inline const SpeciesStatistics& TissueStatistics::at(const DriverGenotypeId& species_id) const
{
    return s_statistics.at(species_id);
}

inline std::chrono::steady_clock::duration TissueStatistics::get_elapsed_time() const
{
    return std::chrono::steady_clock::now()-first_event_time;
}

inline Time TissueStatistics::get_simulated_time() const
{
    return sim_times.back();
}

template <class Rep, std::intmax_t num, std::intmax_t denom>
auto extract_time_units(std::chrono::duration<Rep, std::ratio<num, denom>> duration)
{
    using namespace std::chrono;

    const auto hrs = duration_cast<hours>(duration);
    duration -= hrs;

    const auto mins = duration_cast<minutes>(duration);
    duration -= mins;

    const auto secs = duration_cast<seconds>(duration);

    return std::make_tuple(hrs, mins, secs);
}

template<class ToDuration> 
double TissueStatistics::get_last_recorded_events_over_time() const
{
    using namespace std::chrono;

    auto time = duration_cast<ToDuration>(real_times.back()-real_times.front()).count();

    return static_cast<double>(real_times.size())/time;
}

}

#endif // __RACES_STATISTICS__