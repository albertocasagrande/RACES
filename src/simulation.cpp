/**
 * @file simulation.cpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Define a tumor evolution simulation
 * @version 0.10
 * @date 2023-09-21
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
    logger(nullptr), last_snapshot_time(system_clock::now()), secs_between_snapshots(0), 
    time(0), death_activation_level(1)
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
    std::swap(genomic_mutation_queue, orig.genomic_mutation_queue);
    std::swap(death_enabled, orig.death_enabled);
    std::swap(death_activation_level, orig.death_activation_level);

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
    const auto& num_of_cells{species.num_of_cells()};

    std::list<CellEventType> event_types{CellEventType::DUPLICATE};

    if (death_enabled.count(species.get_id())>0) {
        event_types.push_back(CellEventType::DIE);
    }

    bool selected_new_event = false;
    // deal with exclusively genomic events
    for (const auto& event_type: event_types) {
        const Time event_rate = species.get_rate(event_type);
        const Time r_value = uni_dist(random_gen);
        const Time candidate_delay =  -log(r_value) / (num_of_cells * event_rate);
        if (event.delay>candidate_delay) {
            event.type = event_type;
            event.delay = candidate_delay;
            selected_new_event = true;
        }
    }

    if (selected_new_event) {
        event.position = Position(tissue, species.choose_a_cell(random_gen));
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
        event.position = Position(tissue, species.choose_a_cell(random_gen));
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

CellEvent Simulation::select_next_event()
{
    // the tissue() call checks whether a tissue has been 
    // associated to the simulation and, if this is not the 
    // case, it throws an std::runtime_error 
    (void)tissue();

    std::uniform_real_distribution<double> uni_dist(0.0, 1.0);

    CellEvent event;
    event.delay = std::numeric_limits<Time>::max();

    for (const Species& species: tissue()) {
        select_next_event_in_species(event, tissue(), death_enabled, species, uni_dist, random_gen);
    }

    // update candidate event when genomic mutation is possible
    while (!genomic_mutation_queue.empty()) {
        const TimedGenomicMutation& genomic_mutation = genomic_mutation_queue.top();

        if (genomic_mutation.time < event.delay + time) {
            const auto* cell = choose_a_cell_in_genomic_species(tissue(), genomic_mutation.initial_id, random_gen);
            if (cell != nullptr) {
                event.delay = (genomic_mutation.time >= time? genomic_mutation.time - time:0);
                event.type = CellEventType::DRIVER_GENETIC_MUTATION;
                event.position = Position(tissue(), *cell);
                event.initial_genotype = cell->get_epigenetic_id();

                const Species& initial_species = tissue().get_species(event.initial_genotype);

                const auto& genomic_species = tissue().get_genotype_species(genomic_mutation.final_id);
                const size_t index = Genotype::signature_to_index(initial_species.get_methylation_signature());

                event.final_genotype = genomic_species[index].get_id();

                genomic_mutation_queue.pop();

                return event;
            }
        } else {
            return event;
        }

        genomic_mutation_queue.pop();
    }

    return event;
}

typename Simulation::EventAffectedCells 
Simulation::simulate_death(const Position& position)
{
    auto cell = (*(position.tissue))(position);

    if (!cell.has_driver_mutations()) {
        return {{},{}};
    }

    return {{cell.copy_and_kill()},{}};
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
    //const Direction& push_dir = select_random_value(random_gen, valid_directions);
    const Direction& push_dir = select_push_direction(random_gen, tissue, position, 
                                                      valid_directions);
    affected.lost_cells = tissue.push_cells(position, push_dir);
    
    tissue(position) = parent_cell.generate_descendent(time);

    affected.new_cells.push_back(tissue(position));

    PositionInTissue new_cell_position{position + PositionDelta(push_dir)};
    if (tissue.is_valid(new_cell_position)) {
        tissue(new_cell_position) = parent_cell.generate_descendent(time);

        affected.new_cells.push_back(tissue(new_cell_position));
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
Simulation::choose_a_cell_in_genomic_species(const Tissue& tissue, const GenotypeId& genotype_id, std::mt19937_64& generator)
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
Simulation::add_genomic_mutation(const Genotype& src, const Genotype& dst, const Time time)
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

    genomic_mutation_queue.emplace(src, dst, time);

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

    if (logger != nullptr) {
        delete logger;
    }

    logger = new BinaryLogger();
    
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
    for (const auto& species: tissue()) {
        for (const auto& cell: species) {
            logger->record_initial_cell(cell);
        }
    }
}

Simulation::~Simulation()
{
    if (logger != nullptr) {
        delete logger;
    }
}

}   // Simulation

}   // Drivers

}   // Races
