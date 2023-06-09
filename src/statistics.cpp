/**
 * @file statistics.cpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Define simulation statistics
 * @version 0.4
 * @date 2023-07-12
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

namespace Races {

SpeciesStatistics::SpeciesStatistics():
    rise_time(0), extinction_time(0), total_cells(0),
    curr_cells(0), killed_cells(0), lost_cells(0)
{}

SpeciesStatistics::SpeciesStatistics(const size_t& num_of_cells):
    rise_time(0), extinction_time(0), total_cells(0),
    curr_cells(num_of_cells), killed_cells(0), lost_cells(0)
{}

TissueStatistics::TissueStatistics(const Tissue& tissue):
    s_statistics(), sim_times(), max_stored_times(100), total_events(0)
{
    assert(max_stored_times>0);

    for (const auto& species: tissue) {
        s_statistics.insert({species.get_id(), SpeciesStatistics(species.num_of_cells())});
    }
}

void TissueStatistics::record_death(const EpigeneticGenotypeId& genotype_id, const Time &time)
{
    auto& s_stats = s_statistics.at(genotype_id);

    s_stats.killed_cells += 1;
    s_stats.curr_cells -= 1;
    if (s_stats.curr_cells==0) {
        s_stats.extinction_time = time;
    }
}

void TissueStatistics::record_lost(const EpigeneticGenotypeId& genotype_id, const Time &time)
{
    auto& s_stats = s_statistics.at(genotype_id);

    s_stats.lost_cells += 1;
    s_stats.curr_cells -= 1;
    if (s_stats.curr_cells==0) {
        s_stats.extinction_time = time;
    }
}

void TissueStatistics::record_duplication(const EpigeneticGenotypeId& genotype_id)
{
    auto& s_stats = s_statistics.at(genotype_id);

    s_stats.total_cells += 2;
    s_stats.curr_cells += 1;
}

void TissueStatistics::record_mutation(const EpigeneticGenotypeId& initial_id, const EpigeneticGenotypeId& final_id, const Time &time)
{
    auto& s_stats = s_statistics.at(initial_id);
    s_stats.curr_cells -= 1;
    s_stats.total_cells -= 1;
    if (s_stats.curr_cells==0) {
        s_stats.extinction_time = time;
    }

    auto& s_stats2 = s_statistics.at(final_id);
    s_stats2.total_cells += 1;
    s_stats2.curr_cells += 1;
    if (s_stats2.total_cells==1) {
        s_stats2.rise_time = time;
    }
}

void TissueStatistics::record_duplication_epigenetic_event(const EpigeneticGenotypeId& genotype_id, const EpigeneticGenotypeId& epigenetic_genotype, const Time &time)
{
    record_duplication(genotype_id);
    record_mutation(genotype_id,epigenetic_genotype,time);
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
        case CellEventType::DIE:
            record_death(event.initial_genotype, time);
            break;
        case CellEventType::DUPLICATE:
            record_duplication(event.initial_genotype);
            break;
        case CellEventType::DUPLICATION_AND_EPIGENETIC_EVENT:
            record_duplication(event.initial_genotype);
            record_mutation(event.initial_genotype, 
                            event.final_genotype, time);
            break;
        case CellEventType::DRIVER_SOMATIC_MUTATION:
            record_mutation(event.initial_genotype, 
                            event.final_genotype, time);
            break;
        default:
            throw std::runtime_error("Unsupported event type");
    }
}

}
