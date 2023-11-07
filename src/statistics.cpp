/**
 * @file statistics.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Define simulation statistics
 * @version 0.13
 * @date 2023-11-07
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

#include <cassert>

#include "statistics.hpp"
#include "tissue.hpp"

namespace Races 
{

namespace Drivers
{

namespace Simulation
{

SpeciesStatistics::SpeciesStatistics():
    rise_time(0), extinction_time(0), total_cells(0),
    curr_cells(0), killed_cells(0), lost_cells(0), num_duplications(0)
{}

SpeciesStatistics::SpeciesStatistics(const size_t& num_of_cells):
    rise_time(0), extinction_time(0), total_cells(0),
    curr_cells(num_of_cells), killed_cells(0), lost_cells(0), 
    num_duplications(0)
{}

size_t SpeciesStatistics::num_of_epigenetic_events() const
{
    size_t num_of_events{0};

    for (const auto& [id, num_of_dst_events]: epigenetic_events) {
        num_of_events += num_of_dst_events;
    }

    return num_of_events;
}

size_t SpeciesStatistics::num_of_epigenetic_events(const SpeciesId& dst_id) const
{
    if (epigenetic_events.count(dst_id)==0) {
        return 0;
    }
    return epigenetic_events.at(dst_id);
}

TissueStatistics::TissueStatistics():
    s_statistics(), sim_times(), max_stored_times(100), total_events(0), time_series_sampling_time(1)
{
    assert(max_stored_times>0);
}

TissueStatistics::TissueStatistics(const Time& delta):
    s_statistics(), sim_times(), max_stored_times(100), total_events(0), time_series_sampling_time(delta)
{
    assert(max_stored_times>0);
}


Time TissueStatistics::time_elapsed_from_last_time_series_sample(const Time &time) const
{
    if (time_series.size()==0) {
        return time;
    }

    if (time < time_series.rbegin()->first) {
        throw std::domain_error("The provided time is earlier than the last "
                                "time series sample time");
    }
    return time - time_series.rbegin()->first;
}

void TissueStatistics::sample_time_series_if_needed(const Time &time)
{
    if (time_elapsed_from_last_time_series_sample(time) >= time_series_sampling_time) {
        auto last_time_slot = static_cast<size_t>(trunc(time/time_series_sampling_time))*
                                    time_series_sampling_time;
        if (last_time_slot == time && time > time_series_sampling_time) {
            last_time_slot -= time_series_sampling_time;
        }
        store_current_data_in_time_series(last_time_slot);
    }
}

void TissueStatistics::record_death(const SpeciesId& species_id, const Time &time)
{
    sample_time_series_if_needed(time);

    auto& s_stats = s_statistics[species_id];

    s_stats.killed_cells += 1;
    s_stats.curr_cells -= 1;
    if (s_stats.curr_cells==0) {
        s_stats.extinction_time = time;
    }
}

void TissueStatistics::record_lost(const SpeciesId& species_id, const Time &time)
{
    sample_time_series_if_needed(time);

    auto& s_stats = s_statistics[species_id];

    s_stats.lost_cells += 1;
    s_stats.curr_cells -= 1;
    if (s_stats.curr_cells==0) {
        s_stats.extinction_time = time;
    }
}

void TissueStatistics::record_duplication(const SpeciesId& species_id, const Time &time)
{
    sample_time_series_if_needed(time);

    auto& s_stats = s_statistics[species_id];

    s_stats.total_cells += 2;
    s_stats.curr_cells += 1;
    s_stats.num_duplications += 1;
}

void TissueStatistics::record_species_change(const SpeciesId& src_species, const SpeciesId& dst_species, const Time &time)
{
    sample_time_series_if_needed(time);

    auto& s_stats = s_statistics[src_species];
    s_stats.curr_cells -= 1;
    s_stats.total_cells -= 1;
    if (s_stats.curr_cells==0) {
        s_stats.extinction_time = time;
    }

    auto& s_stats2 = s_statistics[dst_species];
    s_stats2.total_cells += 1;
    s_stats2.curr_cells += 1;
    if (s_stats2.total_cells==1) {
        s_stats2.rise_time = time;
    }
}

void TissueStatistics::record_duplication_and_species_change(const SpeciesId& src_species,
                                                             const SpeciesId& dst_species,
                                                             const Time &time)
{
    sample_time_series_if_needed(time);

    record_duplication(src_species, time);
    record_species_change(src_species, dst_species, time);
}

void TissueStatistics::record_epigenetic_switch(const SpeciesId& src_species, 
                                                const SpeciesId& dst_species,
                                                const Time &time)
{
    record_duplication_and_species_change(src_species, dst_species, time);

    auto& src_stats = s_statistics[src_species];

    ++(src_stats.epigenetic_events[dst_species]);
}

void TissueStatistics::record_event(const CellEvent& event, const Time &time)
{
    ++total_events;

    real_times.push_back(std::chrono::steady_clock::now());
    sim_times.push_back(time);

    if (real_times.size()==1) {
        first_event_time = real_times.back();
    }

    if (sim_times.size()>max_stored_times) {
        sim_times.pop_front();
        real_times.pop_front();
    }

    switch(event.type) {
        case CellEventType::DEATH:
            record_death(event.initial_species, time);
            break;
        case CellEventType::DUPLICATION:
            record_duplication(event.initial_species, time);
            break;
        case CellEventType::EPIGENETIC_SWITCH:
            record_epigenetic_switch(event.initial_species, 
                                     event.final_species, time);
            break;
        case CellEventType::GENOTYPE_MUTATION:
            record_genotype_mutation(event.initial_species, 
                                     event.final_species, time);
            break;
        default:
            throw std::runtime_error("Unsupported event type");
    }
}

}   // Simulation

}   // Drivers

}   // Races
