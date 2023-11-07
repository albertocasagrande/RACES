/**
 * @file ending_conditions.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements simulation ending conditions
 * @version 0.6
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

#include "ending_conditions.hpp"

namespace Races
{

namespace Drivers
{

namespace Simulation
{

TimeTest::TimeTest(const Time& threshold):
    threshold(threshold)
{}

SpeciesCountTest::SpeciesCountTest(const SpeciesId& species_id, const size_t& threshold):
    species_id(species_id), threshold(threshold)
{}

bool SpeciesCountTest::operator()(const Simulation& simulation) const
{
    const auto& species = simulation.tissue().get_species(species_id);

    return threshold < species.num_of_cells();
}

uint8_t SpeciesCountTest::percentage(const Simulation& simulation) const
{
    const auto& species = simulation.tissue().get_species(species_id);

    return static_cast<uint8_t>((100*species.num_of_cells()/threshold));
}

GenotypeCountTest::GenotypeCountTest(const GenotypeId& genotype_id, const size_t& threshold):
    genotype_id(genotype_id), threshold(threshold)
{}

bool GenotypeCountTest::operator()(const Simulation& simulation) const
{
    const auto& genotype_species = simulation.tissue().get_genotype_species(genotype_id);

    return threshold < genotype_species.num_of_cells();
}

uint8_t GenotypeCountTest::percentage(const Simulation& simulation) const
{
    const auto genotype_species = simulation.tissue().get_genotype_species(genotype_id);

    return static_cast<uint8_t>((100*genotype_species.num_of_cells()/threshold));
}

EventCountTest::EventCountTest(const CellEventType& event_type, 
                               const SpeciesId& species_id, const size_t& threshold):
    event_type(event_type), species_id(species_id), 
    dst_id(std::numeric_limits<SpeciesId>::max()), threshold(threshold)
{
    switch(event_type) {
        case CellEventType::DEATH:
        case CellEventType::DUPLICATION:
        case CellEventType::EPIGENETIC_SWITCH:
            return;
        default:
            throw std::domain_error("EventCountTest does not support event "+
                                    cell_event_names[event_type]);
    }
}

EventCountTest::EventCountTest(const SpeciesId& src_id, const SpeciesId& dst_id,
                               const size_t& threshold):
    event_type(CellEventType::EPIGENETIC_SWITCH), species_id(src_id), dst_id(dst_id), 
    threshold(threshold)
{}

//@private
size_t get_epigenetic_events_to(const SpeciesStatistics& s_stats, const SpeciesId& dst_id)
{
    if (dst_id != std::numeric_limits<SpeciesId>::max()) {
        return s_stats.num_of_epigenetic_events(dst_id);
    }

    return s_stats.num_of_epigenetic_events();
}

size_t EventCountTest::get_event_number(const Simulation& simulation) const
{
    const auto& t_stats = simulation.get_statistics();

    if (!t_stats.contains_data_for(species_id)) {
        return false;
    }

    const auto& s_stats = t_stats.at(species_id);

    switch(event_type) {
        case CellEventType::DEATH:
            return s_stats.killed_cells;
        case CellEventType::DUPLICATION:
            return s_stats.num_duplications;
        case CellEventType::EPIGENETIC_SWITCH:
            return get_epigenetic_events_to(s_stats, dst_id);
        default:
            throw std::domain_error("EventCountTest does not support event "+
                                    cell_event_names[event_type]);
    }
}

bool EventCountTest::operator()(const Simulation& simulation) const
{
    return threshold < get_event_number(simulation);
}

uint8_t EventCountTest::percentage(const Simulation& simulation) const
{
    return static_cast<uint8_t>((100*get_event_number(simulation)/threshold));
}

}   // Simulation

}   // Drivers

}   // Races

