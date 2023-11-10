/**
 * @file simulation.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Define a tumor evolution simulation
 * @version 0.41
 * @date 2023-11-10
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

#include <set>
#include <list>

namespace Races 
{

namespace Drivers
{

namespace Simulation
{

Simulation::AddedCell::AddedCell():
    species_id(WILD_TYPE_SPECIES)
{}

Simulation::AddedCell::AddedCell(const SpeciesId& species, const PositionInTissue& position, 
                                 const Time& time):
    PositionInTissue(position), species_id{species}, time{time}
{}

Simulation::Simulation(int random_seed):
    logger(), last_snapshot_time(system_clock::now()), secs_between_snapshots(0), 
    time(0), death_activation_level(1), duplicate_internal_cells(true), storage_enabled(true)
{
    random_gen.seed(random_seed);

    // Create a default tissue
    set_tissue("A tissue", {1000, 1000});
}

Simulation::Simulation(const std::filesystem::path& log_directory, int random_seed):
    logger(log_directory), last_snapshot_time(system_clock::now()), secs_between_snapshots(0), 
    time(0), death_activation_level(1), duplicate_internal_cells(true), storage_enabled(true)
{
    random_gen.seed(random_seed);

    // Create a default tissue
    set_tissue("A tissue", {1000, 1000});
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
                                      const std::set<SpeciesId>& death_enabled, 
                                      const Species& species,
                                      std::uniform_real_distribution<double>& uni_dist,
                                      GENERATOR_TYPE& random_gen)
{
    std::list<CellEventType> event_types{CellEventType::DUPLICATION};

    if (death_enabled.count(species.get_id())>0) {
        event_types.push_back(CellEventType::DEATH);
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
        event.initial_species = species.get_id();
    }
}

template<typename GENERATOR_TYPE>
void select_epigenetic_event_in_species(CellEvent& event, Tissue& tissue, const Species& species,
                                        std::uniform_real_distribution<double>& uni_dist, GENERATOR_TYPE& random_gen)
{
    const auto& num_of_cells{species.num_of_cells()};

    bool selected_new_event = false;
    // Deal with possible epigenetic events whenever admitted
    for (const auto& [dst_id, event_rate]: species.get_epigenetic_switch_rates()) {
        const Time r_value = uni_dist(random_gen);
        const Time candidate_delay = -log(r_value) / (num_of_cells * event_rate);
        
        if (event.delay>candidate_delay) {
            event.delay = candidate_delay;
            event.final_species = dst_id;
            selected_new_event = true;
        }
    }

    if (selected_new_event) {
        event.type = CellEventType::EPIGENETIC_SWITCH;
        event.position = Position(tissue, species.choose_a_cell(random_gen, event.type));
        event.initial_species = species.get_id();
    }
}

template<typename GENERATOR_TYPE>
void select_next_event_in_species(CellEvent& event, Tissue& tissue,
                                  const std::set<SpeciesId>& death_enabled, 
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

const CellInTissue&
Simulation::choose_cell_in(const GenotypeId& genotype_id)
{
    std::vector<size_t> num_of_cells;

    size_t total = 0;

    const auto& species = tissue().get_genotype_species(genotype_id);
    for (const auto& S : species) {
        total += S.num_of_cells();
        num_of_cells.push_back(total);
    }

    if (total == 0) {
        throw std::runtime_error("No cell available for the genotype "+std::to_string(genotype_id));
    }

    std::uniform_int_distribution<size_t> distribution(1,total);
    const size_t value = distribution(random_gen);

    size_t i=0;
    while (value > num_of_cells[i]) {
        ++i;
    }

    return species[i].choose_a_cell(random_gen);
}

const CellInTissue& Simulation::choose_cell_in(const std::string& genotype_name)
{
    auto genotype_id = find_genotype_id(genotype_name);

    return choose_cell_in(genotype_id);
}

const CellInTissue& Simulation::choose_cell_in(const GenotypeId& genotype_id,
                                               const RectangleSet& rectangle)
{
    std::set<SpeciesId> species_ids;
    auto genotype_species = tissue().get_genotype_species(genotype_id);

    for (const auto& species : genotype_species) {
        species_ids.insert(species.get_id());
    }

    std::list<PositionInTissue> position_vector;
    for (const auto& position : rectangle) {
        if (!tissue()(position).is_wild_type()) {
            const CellInTissue& cell = tissue()(position);

            if (species_ids.count(cell.get_species_id())>0) {
                position_vector.push_back(position);
            }
        }
    }

    std::uniform_int_distribution<size_t> distribution(1,position_vector.size());

    auto selected = distribution(random_gen);

    size_t i{0};
    for (const auto& position: position_vector) {
        if (++i==selected) {
            return tissue()(position);
        }
    }

    std::ostringstream oss;

    oss <<"No cell having the specified genotype is available in " << rectangle;
    throw std::domain_error(oss.str());
}

const CellInTissue& Simulation::choose_cell_in(const std::string& genotype_name,
                                               const RectangleSet& rectangle)
{
    auto genotype_id = find_genotype_id(genotype_name);

    return choose_cell_in(genotype_id, rectangle);
}

/**
 * @brief Create a genotype mutation event
 * 
 * @param tissue is the tissue in which the event will occurs
 * @param position is the position of the parent cell that will give birth to the mutated cell
 * @param final_id is the genotype identifier of the mutated cell
 * @param delay is the delay of the mutation with respect to the current simulation clock
 * @return the created cell event 
 */
CellEvent create_genotype_mutation_event(Tissue& tissue, const PositionInTissue& position,
                                         const GenotypeId& final_id, const Time& delay)
{
    CellEvent event;

    event.delay = delay;
    event.type = CellEventType::GENOTYPE_MUTATION;
    event.position = Position(tissue, position);
    event.initial_species = static_cast<const CellInTissue&>(tissue(position)).get_species_id();

    const Species& initial_species = tissue.get_species(event.initial_species);

    const auto& genotype_species = tissue.get_genotype_species(final_id);
    const size_t index = GenotypeProperties::signature_to_index(initial_species.get_methylation_signature());

    event.final_species = genotype_species[index].get_id();

    return event;
}

bool Simulation::handle_timed_genotype_mutation(const TimedEvent& timed_genotype_mutation, CellEvent& candidate_event)
{
    const auto& genotype_mutation = timed_genotype_mutation.get_event<GenotypeMutation>();

    try {
        const auto& cell = choose_cell_in(genotype_mutation.initial_id);

        auto delay = (timed_genotype_mutation.time >= time ? 
                        timed_genotype_mutation.time-time : 0);
        candidate_event = create_genotype_mutation_event(tissue(), cell, genotype_mutation.final_id, delay);

        return true;
    } catch (std::runtime_error&) {
        return false;
    }
}

Simulation& Simulation::simulate_genotype_mutation(const PositionInTissue& position, 
                                                   const std::string& dst_genotype_name)
{
    auto dst_genotype_id = find_genotype_id(dst_genotype_name);

    return simulate_genotype_mutation(position, dst_genotype_id);
}

Simulation& Simulation::simulate_genotype_mutation(const PositionInTissue& position,
                                                   const GenotypeId& dst_genotype_id)
{
    auto mutation_event = create_genotype_mutation_event(tissue(), position, 
                                                         dst_genotype_id, 0);

    simulate(mutation_event);

    if (storage_enabled) {
        logger.snapshot(*this);
        logger.flush_archives();
    }

    return *this;
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
    
    const Time previous_time = time;

    time = timed_sampling.time;

    // sample tissue
    TissueSample sample = sample_tissue(sampling.get_region(), sampling.is_preserving_tissue());

    if (sampling.get_name()=="") {
        sample.set_name(sample.get_default_name());
    } else {
        sample.set_name(sampling.get_name());
    }

    if (storage_enabled) {
        logger.save_sample(sample);
    }

    // update candidate event to avoid removed regions
    candidate_event = select_next_cell_event();

    if (candidate_event.delay+previous_time < time) {
        candidate_event.delay = 0;
    } else {
        candidate_event.delay -= (time - previous_time);
    }
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
            case TimedEvent::Type::GENOTYPE_MUTATION:
                {
                    auto candidate_updated = handle_timed_genotype_mutation(timed_event, candidate_event);

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

template<typename T>
inline
T get_min(const T v1, const size_t& v2)
{
    return std::min(v1, static_cast<T>(v2));
}

void enable_duplication_on_neighborhood_externals(Tissue& tissue, const PositionInTissue& position)
{
    auto sizes = tissue.size();
    PositionInTissue pos;
    pos.x = (position.x>0?position.x-1:0);
    auto max_x = get_min(position.x+2, sizes[0]);
    for (; pos.x < max_x; ++pos.x) {
        pos.y = (position.y>0?position.y-1:0);
        auto max_y = get_min(position.y+2, sizes[1]);
        for (; pos.y < max_y; ++pos.y) {
            pos.z = (position.z>0?position.z-1:0);
            auto max_z = (sizes.size()==2?1:get_min(position.z+2, sizes[2]));
            for (; pos.z < max_z; ++pos.z) {
                Tissue::CellInTissueProxy cell_in_tissue = tissue(pos);
                if (!cell_in_tissue.is_wild_type() && cell_in_tissue.is_on_border()) {
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

    if (cell.is_wild_type()) {
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
Simulation::simulate_mutation(const Position& position, const SpeciesId& final_id)
{
    Tissue& tissue = *(position.tissue);

    Cell parent_cell = tissue(position);

    parent_cell.species_id = final_id;

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
                if (!cell_in_tissue.is_wild_type() && !cell_in_tissue.is_on_border()) {
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

    if (tissue(position).is_wild_type()) {
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
Simulation::simulate_duplication_and_mutation_event(const Position& position, const SpeciesId& final_id)
{
    Tissue& tissue = *(position.tissue);

    auto affected = simulate_duplication(position);

    if (affected.new_cells.size()==0) {
        return affected;
    }

    Cell cell = tissue(position);

    if (!lineage_graph.has_edge(cell.get_species_id(), final_id)) {
        lineage_graph.add_edge(cell.get_species_id(), final_id, time);
    }

    cell.species_id = final_id;

    tissue(position) = cell;

    for (auto& cell : affected.new_cells) {
        cell = static_cast<Cell>(tissue(cell));
    }

    return affected;
}

Simulation&
Simulation::schedule_genotype_mutation(const GenotypeProperties& src, const GenotypeProperties& dst, const Time time)
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

    GenotypeMutation mutation(src, dst);

    timed_event_queue.emplace(time, mutation);

    return *this;
}

const GenotypeId& Simulation::find_genotype_id(const std::string& genotype_name) const
{
    auto found_genotype = genotype_name2id.find(genotype_name);
    if (found_genotype == genotype_name2id.end()) {
        throw std::out_of_range("Unknown genotype name \""+genotype_name+"\"");
    }

    return found_genotype->second;
}

const std::string& Simulation::find_genotype_name(const GenotypeId& genotype_id) const
{
    for (const auto& [name, id] : genotype_name2id) {
        if (id == genotype_id) {
            return name;
        }
    }

    throw std::out_of_range("Unknown genotype id "+std::to_string(genotype_id));
}

Simulation&
Simulation::schedule_genotype_mutation(const std::string& src, const std::string& dst, const Time time)
{
    // the tissue() call checks whether a tissue has been 
    // associated to the simulation and, if this is not the 
    // case, it throws an std::runtime_error 
    (void)tissue();

    auto src_id = find_genotype_id(src);
    auto dst_id = find_genotype_id(dst);

    GenotypeMutation mutation(src_id, dst_id);

    timed_event_queue.emplace(time, mutation);

    return *this;
}

Simulation&
Simulation::schedule_timed_event(const TimedEvent& timed_event)
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
    added_cells.clear();

    statistics = TissueStatistics();
    time = 0; 
    
    death_enabled.clear();

    logger = BinaryLogger(logger.get_directory());
    
    //death_activation_level = 1;
}

Simulation& Simulation::add_genotype(const GenotypeProperties& genotype)
{
    if (genotype_name2id.count(genotype.get_name())>0) {
        throw std::out_of_range("Genotype \""+genotype.get_name()+"\" already in the simulation");
    }

    tissue().add_genotype(genotype);

    genotype_name2id[genotype.get_name()] = genotype.get_id();

    return *this;
}

Simulation& Simulation::add_cell(const SpeciesId& species_id, const PositionInTissue& position)
{
    tissue().add_cell(species_id, position);

    added_cells.emplace_back(species_id, position, time);

    LineageEdge edge{WILD_TYPE_SPECIES, species_id};
    if (!lineage_graph.has_edge(edge)) {
        lineage_graph.add_edge(edge, time);
    }

    return *this;
}

Simulation& Simulation::set_tissue(const std::string& name, const std::vector<AxisSize>& sizes)
{
    if (time>0) {
        std::runtime_error("The tissue properties can be exclusively changed before"
                           "the simulation beginning");
    }

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
            if (!cell.is_wild_type()) {
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
            if (!cell.is_wild_type()) {
                auto cell_in_tissue = static_cast<const CellInTissue&>(cell);

                sample.add_cell_id(cell_in_tissue.get_id());
                if (!preserve_tissue) {
                    ++(statistics[cell_in_tissue.get_species_id()].lost_cells);
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
