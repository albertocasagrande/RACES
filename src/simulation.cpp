/**
 * @file simulation.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Define a tumor evolution simulation
 * @version 0.23
 * @date 2023-10-25
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

#include "simulation.hpp"

namespace Races 
{

namespace Drivers
{

namespace Simulation
{


Simulation::Simulation(int random_seed):
    logger(), last_snapshot_time(system_clock::now()), secs_between_snapshots(0), 
    time(0), death_activation_level(1), duplicate_internal_cells(true), storage_enabled(true)
{
    random_gen.seed(random_seed);
}

Simulation::Simulation(Simulation&& orig):
    Simulation()
{
    *this = std::move(orig);
}

Simulation& Simulation::operator=(Simulation&& orig)
{
    std::swap(logger, orig.logger);
    std::swap(tissues, orig.tissues);
    std::swap(valid_directions, orig.valid_directions);
    std::swap(last_snapshot_time, orig.last_snapshot_time);
    std::swap(secs_between_snapshots, orig.secs_between_snapshots);
    std::swap(statistics, orig.statistics);
    std::swap(time, orig.time);
    std::swap(random_gen, orig.random_gen);
    std::swap(timed_event_queue, orig.timed_event_queue);
    std::swap(death_enabled, orig.death_enabled);
    std::swap(death_activation_level, orig.death_activation_level);
    std::swap(duplicate_internal_cells, orig.duplicate_internal_cells);

    return *this;
}

const Tissue& Simulation::tissue() const
{
    if (tissues.size()==0) {
        throw std::runtime_error("No tissue has been associated to the simulation yet.");
    }
    return tissues[0];
}

Tissue& Simulation::tissue()
{
    if (tissues.size()==0) {
        throw std::runtime_error("No tissue has been associated to the simulation yet.");
    }
    return tissues[0];
}

template<typename GENERATOR_TYPE>
void select_liveness_event_in_species(CellEvent& event, Tissue& tissue, 
                                      const std::set<EpigeneticGenotypeId>& death_enabled, 
                                      const Species& species,
                                      std::uniform_real_distribution<double>& uni_dist,
                                      GENERATOR_TYPE& random_gen)
{
    std::list<CellEventType> event_types{CellEventType::DUPLICATE};

    if (death_enabled.count(species.get_id())>0) {
        event_types.push_back(CellEventType::DIE);
    }

    bool selected_new_event = false;
    // deal with exclusively genomic events
    for (const auto& event_type: event_types) {
        const Time event_rate = species.get_rate(event_type);
        const Time r_value = uni_dist(random_gen);
        size_t num_of_cells;
        
        num_of_cells = species.num_of_cells_available_for(event_type);
        //num_of_cells = species.num_of_cells();
        const Time candidate_delay =  -log(r_value) / (num_of_cells * event_rate);
        if (event.delay>candidate_delay) {
            event.type = event_type;
            event.delay = candidate_delay;
            selected_new_event = true;
        }
    }

    if (selected_new_event) {
        event.position = Position(tissue, species.choose_a_cell(random_gen, event.type));
        event.initial_genotype = species.get_id();
    }
}

template<typename GENERATOR_TYPE>
void select_epigenetic_event_in_species(CellEvent& event, Tissue& tissue, const Species& species,
                                        std::uniform_real_distribution<double>& uni_dist, GENERATOR_TYPE& random_gen)
{
    const auto& num_of_cells{species.num_of_cells()};

    bool selected_new_event = false;
    // Deal with possible epigenetic events whenever admitted
    for (const auto& [dst_id, event_rate]: species.get_epigenetic_mutation_rates()) {
        const Time r_value = uni_dist(random_gen);
        const Time candidate_delay = -log(r_value) / (num_of_cells * event_rate);
        
        if (event.delay>candidate_delay) {
            event.delay = candidate_delay;
            event.final_genotype = dst_id;
            selected_new_event = true;
        }
    }

    if (selected_new_event) {
        event.type = CellEventType::DUPLICATION_AND_EPIGENETIC_EVENT;
        event.position = Position(tissue, species.choose_a_cell(random_gen, event.type));
        event.initial_genotype = species.get_id();
    }
}

template<typename GENERATOR_TYPE>
void select_next_event_in_species(CellEvent& event, Tissue& tissue,
                                  const std::set<EpigeneticGenotypeId>& death_enabled, 
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

bool Simulation::handle_timed_driver_mutation(const TimedEvent& timed_driver_mutation, CellEvent& candidate_event)
{
    const auto& driver_mutation = timed_driver_mutation.get_event<DriverMutation>();
    const auto* cell = choose_a_cell_in_genomic_species(random_gen, tissue(), driver_mutation.initial_id);
    if (cell != nullptr) {
        candidate_event.delay = (timed_driver_mutation.time >= time? timed_driver_mutation.time - time:0);
        candidate_event.type = CellEventType::DRIVER_GENETIC_MUTATION;
        candidate_event.position = Position(tissue(), *cell);
        candidate_event.initial_genotype = cell->get_epigenetic_id();

        const Species& initial_species = tissue().get_species(candidate_event.initial_genotype);

        const auto& genomic_species = tissue().get_genotype_species(driver_mutation.final_id);
        const size_t index = Genotype::signature_to_index(initial_species.get_methylation_signature());

        candidate_event.final_genotype = genomic_species[index].get_id();

        return true;
    }

    return false;
}

void Simulation::handle_timed_rate_update(const TimedEvent& timed_rate_update)
{
    const auto& rate_update = timed_rate_update.get_event<RateUpdate>();
    Species& species = tissue().get_species(rate_update.species_id);
    
    species.set_rate(rate_update.event_type, rate_update.new_rate);
}

void Simulation::handle_timed_sampling(const TimedEvent& timed_sampling, CellEvent& candidate_event)
{
    const auto& sampling = timed_sampling.get_event<Sampling>();
    
    // sample tissue
    for (const auto& region : sampling.sample_set) {
        TissueSample sample = sample_tissue(region, sampling.preserve_tissue);

        if (storage_enabled) {
            logger.save_sample(sample);
        }
    }

    // update candidate event to avoid removed regions
    candidate_event = select_next_cell_event();

    if (candidate_event.delay+time < timed_sampling.time) {
        candidate_event.delay = 0;
    } else {
        candidate_event.delay -= (time - timed_sampling.time);
    }

    time = timed_sampling.time;
}

void Simulation::handle_timed_event_queue(CellEvent& candidate_event)
{
    // if the timed event queue is not empty and the next time event occurs before the candidate cell event
    while (!timed_event_queue.empty() && (timed_event_queue.top().time <= candidate_event.delay + time)) {
        TimedEvent timed_event = timed_event_queue.top();

        timed_event_queue.pop();

        if (storage_enabled) {
            if (timed_event_queue.top().time != timed_event.time 
                || timed_event_queue.top().type != timed_event.type) {
                logger.snapshot(*this);
            }
        }

        switch(timed_event.type) {
            case TimedEvent::Type::DRIVER_MUTATION:
                {
                    auto candidate_updated = handle_timed_driver_mutation(timed_event, candidate_event);

                    if (candidate_updated) {
                        return;
                    }
                }
                break;
            case TimedEvent::Type::LIVENESS_RATE_UPDATE:
                {
                    handle_timed_rate_update(timed_event);
                }
                break;
            case TimedEvent::Type::SAMPLING:
                {
                    handle_timed_sampling(timed_event, candidate_event);
                }
                break;
            default:
                std::cout << "Unsupported timed event" << std::endl;
                exit(1);
        }
    }
}

CellEvent Simulation::select_next_cell_event()
{
    std::uniform_real_distribution<double> uni_dist(0.0, 1.0);

    CellEvent event;
    event.delay = std::numeric_limits<Time>::max();

    for (const Species& species: tissue()) {
        select_next_event_in_species(event, tissue(), death_enabled, species, uni_dist, random_gen);
    }

    return event;
}

CellEvent Simulation::select_next_event()
{
    // the tissue() call checks whether a tissue has been 
    // associated to the simulation and, if this is not the 
    // case, it throws an std::runtime_error 
    (void)tissue();

    CellEvent event = select_next_cell_event();

    handle_timed_event_queue(event);

    return event;
}

void enable_duplication_on_neighborhood_externals(Tissue& tissue, const PositionInTissue& position)
{
    auto sizes = tissue.size();
    PositionInTissue pos;
    pos.x = (position.x>0?position.x-1:0);
    for (; (pos.x < position.x+2 &&  pos.x < sizes[0]); ++pos.x) {
        pos.y = (position.y>0?position.y-1:0);
        for (; (pos.y < position.y+2 &&  pos.y < sizes[1]); ++pos.y) {
            pos.z = (position.z>0?position.z-1:0);
            for (; ((pos.z < position.z+2 &&  pos.z < sizes[2])
                    || (sizes[2]==0 && pos.z==0)); ++pos.z) {
                Tissue::CellInTissueProxy cell_in_tissue = tissue(pos);
                if (cell_in_tissue.has_driver_mutations() && cell_in_tissue.is_on_border()) {
                    cell_in_tissue.enable_duplication();
                }
            }                     
        }   
    }
}

typename Simulation::EventAffectedCells 
Simulation::simulate_death(const Position& position)
{
    auto cell = (*(position.tissue))(position);

    if (!cell.has_driver_mutations()) {
        return {{},{}};
    }

    Simulation::EventAffectedCells affected = {{cell.copy_and_erase()},{}};

    enable_duplication_on_neighborhood_externals(*(position.tissue), position);

    return affected;
}

template<typename GENERATOR>
inline const Direction& select_push_direction(GENERATOR& random_gen,
                                              const Tissue& tissue,
                                              const PositionInTissue& position,
                                              const std::vector<Direction>& directions)
{
    std::vector<double> cells_to_push;
    double total{0.0};
    for (const auto& direction : directions) {
        size_t num_of_cells = tissue.count_driver_cells_from(position, direction);
        cells_to_push.push_back(static_cast<double>(1.0)/num_of_cells);
        total += cells_to_push.back(); 
    }

    std::uniform_real_distribution<double> distribution(0,total);

    auto selected = distribution(random_gen);

    total = 0.0;
    auto it = cells_to_push.begin();
    for (const auto& direction : directions) {
        total += *it;
        ++it;

        if (total >= selected) {
            return direction;
        }
    }

    return directions.back();
}


template<typename GENERATOR>
inline const Direction& select_min_push_direction(GENERATOR& random_gen,
                                                  const Tissue& tissue,
                                                  const PositionInTissue& position,
                                                  const std::vector<Direction>& directions)
{
    std::list<size_t> cells_to_push;
    size_t min_cells_to_push = std::numeric_limits<size_t>::max();
    uint8_t num_of_mins = 0;
    for (const auto& direction : directions) {
        cells_to_push.push_back(tissue.count_driver_cells_from(position, direction));
        if (min_cells_to_push>cells_to_push.back()) {
            num_of_mins = 1;
            min_cells_to_push = cells_to_push.back();
        } else {
            if (min_cells_to_push==cells_to_push.back()) {
                ++num_of_mins;
            }
        }
    }

    std::uniform_int_distribution<uint8_t> distribution(1,num_of_mins);

    auto selected = distribution(random_gen);

    auto it = cells_to_push.begin();
    uint8_t counter = 0;
    for (const auto& direction : directions) {
        if (*it==min_cells_to_push) {
            if (++counter==selected) {
                return direction;
            }
        }
        ++it;
    }

    throw std::runtime_error("The direction has not been selected");
}

template<typename GENERATOR, typename T>
inline const T& select_random_value(GENERATOR& random_gen, const std::vector<T>& values)
{
    std::uniform_int_distribution<size_t> distribution(0,values.size()-1);

    return values[distribution(random_gen)];
}

typename Simulation::EventAffectedCells 
Simulation::simulate_mutation(const Position& position, const EpigeneticGenotypeId& final_id)
{
    Tissue& tissue = *(position.tissue);

    Cell parent_cell = tissue(position);

    parent_cell.epigenetic_id = final_id;

    tissue(position) = parent_cell.generate_descendent(time);

    return {{tissue(position)},{}};
}


void disable_duplication_on_neighborhood_internals(Tissue& tissue, const PositionInTissue& position)
{
    auto sizes = tissue.size();
    PositionInTissue pos;
    pos.x = (position.x>0?position.x-1:0);
    for (; (pos.x < position.x+2 &&  pos.x < sizes[0]); ++pos.x) {
        pos.y = (position.y>0?position.y-1:0);
        for (; (pos.y < position.y+2 &&  pos.y < sizes[1]); ++pos.y) {
            pos.z = (position.z>0?position.z-1:0);
            for (; ((pos.z < position.z+2 &&  pos.z < sizes[2])
                    || (sizes[2]==0 && pos.z==0)); ++pos.z) {
                Tissue::CellInTissueProxy cell_in_tissue = tissue(pos);
                if (cell_in_tissue.has_driver_mutations() && !cell_in_tissue.is_on_border()) {
                    cell_in_tissue.disable_duplication();
                }
            }                     
        }   
    }
}

typename Simulation::EventAffectedCells 
Simulation::simulate_duplication(const Position& position)
{
    Tissue& tissue = *(position.tissue);
    EventAffectedCells affected;

    if (!tissue(position).has_driver_mutations()) {
        return affected;
    }

    Cell parent_cell = tissue(position);

    // push the cell in position towards a random direction
    const Direction& push_dir = select_random_value(random_gen, valid_directions);
    //const Direction& push_dir = select_min_push_direction(random_gen, tissue, position, 
    //                                                      valid_directions);
    
    affected.lost_cells = tissue.push_cells(position, push_dir);
    
    Tissue::CellInTissueProxy cell_in_tissue = tissue(position);

    cell_in_tissue = parent_cell.generate_descendent(time);

    affected.new_cells.push_back(cell_in_tissue);

    PositionInTissue new_cell_position{position + PositionDelta(push_dir)};
    if (tissue.is_valid(new_cell_position)) {
        Tissue::CellInTissueProxy cell_in_tissue = tissue(new_cell_position);

        cell_in_tissue = parent_cell.generate_descendent(time);

        if (!duplicate_internal_cells) {
            disable_duplication_on_neighborhood_internals(tissue, new_cell_position);
        }

        affected.new_cells.push_back(cell_in_tissue);
    }

    return affected;
}

typename Simulation::EventAffectedCells 
Simulation::simulate_duplication_epigenetic_event(const Position& position, const EpigeneticGenotypeId& final_id)
{
    Tissue& tissue = *(position.tissue);

    auto affected = simulate_duplication(position);

    if (affected.new_cells.size()==0) {
        return affected;
    }

    Cell cell = tissue(position);

    cell.epigenetic_id = final_id;

    tissue(position) = cell;

    for (auto& cell : affected.new_cells) {
        cell = static_cast<Cell>(tissue(cell));
    }

    return affected;
}

const CellInTissue*
Simulation::choose_a_cell_in_genomic_species(std::mt19937_64& generator, const Tissue& tissue, const GenotypeId& genotype_id)
{
    std::vector<size_t> num_of_cells;

    size_t total = 0;

    const auto& species = tissue.get_genotype_species(genotype_id);
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

Simulation&
Simulation::add_driver_mutation(const Genotype& src, const Genotype& dst, const Time time)
{
    // the tissue() call checks whether a tissue has been 
    // associated to the simulation and, if this is not the 
    // case, it throws an std::runtime_error 
    (void)tissue();

    if (src.num_of_promoters()>dst.num_of_promoters()) {
        std::ostringstream oss;

        oss << "Incompatible number of methylable promoters: \"" << src.get_name() 
            << "\" must have at least as many methylable promoters as \"" << dst.get_name() 
            << "\".";
        throw std::domain_error(oss.str());
    }

    DriverMutation mutation(src, dst);

    timed_event_queue.emplace(time, mutation);

    return *this;
}

Simulation&
Simulation::add_timed_event(const TimedEvent& timed_event)
{
    // the tissue() call checks whether a tissue has been 
    // associated to the simulation and, if this is not the 
    // case, it throws an std::runtime_error 
    (void)tissue();

    timed_event_queue.push(timed_event);

    return *this;  
}

void Simulation::init_valid_directions()
{
    if (tissues.size()==0) {
        return;
    }

    valid_directions.clear();

    for (const auto &x_move : {Direction::X_UP, Direction::X_DOWN, Direction::X_NULL}) {
        for (const auto &y_move : {Direction::Y_UP, Direction::Y_DOWN, Direction::Y_NULL}) {
            if (tissue().size().size()==3) {
                for (const auto &z_move : {Direction::Z_UP, Direction::Z_DOWN, Direction::Z_NULL}) {
                    valid_directions.push_back(x_move|y_move|z_move);
                }
            } else {
                valid_directions.push_back(x_move|y_move);
            }
        }
    }

    // remove null move
    valid_directions.pop_back();
}

void Simulation::reset()
{
    tissues.clear();

    statistics = TissueStatistics();
    time = 0; 
    
    death_enabled.clear();

    logger = BinaryLogger();
    
    //death_activation_level = 1;
}

Simulation& Simulation::set_tissue(const std::string& name, const std::vector<AxisSize>& sizes)
{
    if (tissues.size()>0) {
        reset();
    }

    Cell::counter = 0;

    tissues.push_back(Tissue(name, sizes));

    init_valid_directions();

    return *this;
}

void Simulation::log_initial_cells()
{
    if (storage_enabled) {
        for (const auto& species: tissue()) {
            for (const auto& cell: species) {
                logger.record_initial_cell(cell);
            }
        }
    }
}

TissueSample
Simulation::sample_tissue(const RectangleSet& rectangle) const
{
    TissueSample sample(time, rectangle);

    for (const auto& position: rectangle) {
        if (tissue().is_valid(position)) {
            auto cell = tissue()(position);
            if (cell.has_driver_mutations()) {
                sample.add_cell_id(static_cast<const CellInTissue&>(cell).get_id());
            }
        }
    }

    return sample;
}

TissueSample
Simulation::sample_tissue(const RectangleSet& rectangle, const bool& preserve_tissue)
{
    TissueSample sample(time, rectangle);

    for (const auto& position: rectangle) {
        if (tissue().is_valid(position)) {
            auto cell = tissue()(position);
            if (cell.has_driver_mutations()) {
                auto cell_in_tissue = static_cast<const CellInTissue&>(cell);

                sample.add_cell_id(cell_in_tissue.get_id());
                if (!preserve_tissue) {
                    ++(statistics[cell_in_tissue.get_epigenetic_id()].lost_cells);
                    cell.erase();
                }
            }
        }
    }

    return sample;
}

Simulation::~Simulation()
{
}

}   // Simulation

}   // Drivers

}   // Races
