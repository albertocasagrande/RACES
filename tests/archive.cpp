/**
 * @file archive.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Some archive tests
 * @version 0.19
 * @date 2023-11-05
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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE archive

#include <sstream>
#include <filesystem>

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include "archive.hpp"
#include "simulation.hpp"
#include "ending_conditions.hpp"
#include "logger.hpp"


struct ArchiveFixture {
    long double time_horizon;

    Races::Drivers::Simulation::Simulation simulation;

    ArchiveFixture():
        time_horizon(70), simulation()
    {
        using namespace Races::Drivers;
        
        GenotypeProperties A("A",{{0.01,0.01}});
        A["-"].set_rates({{CellEventType::DEATH, 0.1},
                          {CellEventType::DUPLICATION, 0.3}});
        A["+"].set_rates({{CellEventType::DEATH, 0.1},
                          {CellEventType::DUPLICATION, 0.45}});
    
        GenotypeProperties B("B",{{0.01,0.01}});
        B["-"].set_rates({{CellEventType::DEATH, 0.1},
                          {CellEventType::DUPLICATION, 0.2}});
        B["+"].set_rates({{CellEventType::DEATH, 0.01},
                          {CellEventType::DUPLICATION, 0.02}});

        simulation.set_tissue("Liver", {1000,1000})
                  .add_genotype(A)
                  .add_genotype(B)
                  .add_cell(B["-"], {250, 500})
                  .schedule_genotype_mutation(B,A,50);

        simulation.death_activation_level = 100;
        simulation.storage_enabled = false;

        Races::Drivers::Simulation::TimeTest done(time_horizon);

        simulation.run(done);
    }

    ~ArchiveFixture()
    {}
};

std::filesystem::path get_a_temporary_path()
{
    namespace fs = std::filesystem;

    fs::path tmp_dir = fs::temp_directory_path();

    std::string basename{"file_"};
    std::string name = basename + "0";
    size_t counter = 0;

    while (fs::exists(tmp_dir / name)) {
        std::ostringstream oss;

        oss << basename << ++counter;

        name = oss.str();
    }

    return tmp_dir / name;
}

template<typename T>
void basic_type_test(const std::vector<T>& to_save)
{
    auto filename = get_a_temporary_path();
    {
        Races::Archive::Binary::Out o_archive(filename);

        for (const auto& value : to_save) {
            o_archive & value;
        }
    }

    {
        Races::Archive::Binary::In i_archive(filename);

        for (const auto& value : to_save) {
            T read_value;

            i_archive & read_value;

            BOOST_CHECK(read_value==value);
        }
    }

    std::filesystem::remove(filename);
}

template<typename T>
bool operator==(const std::vector<T>& a, const std::vector<T>& b)
{
    if (a.size()!=b.size()) {
        return false;
    }

    auto a_it=a.begin(), b_it=b.begin();
    for (; a_it!=a.end(); ++a_it, ++b_it) {
        if (*a_it!=*b_it) {
            return false;
        }
    }

    return true;
}

template<typename KEY, typename T>
bool operator==(const std::map<KEY,T>& a, const std::map<KEY,T>& b)
{
    if (a.size()!=b.size()) {
        return false;
    }

    auto a_it=a.begin(), b_it=b.begin();
    for (; a_it!=a.end(); ++a_it, ++b_it) {
        if (a_it->first!=b_it->first || a_it->second!=b_it->second) {
            return false;
        }
    }

    return true;
}

bool operator==(const Races::Drivers::Cell& a, const Races::Drivers::Cell& b)
{
    return (a.get_id()==b.get_id() && 
            a.get_parent_id()==b.get_parent_id() &&
            a.get_species_id()==b.get_species_id());
}

bool operator==(const Races::Drivers::Simulation::CellInTissue& a, const Races::Drivers::Simulation::CellInTissue& b)
{
    using namespace  Races::Drivers;

    return (static_cast<const Cell&>(a)==static_cast<const Cell&>(b) &&
            static_cast<const Simulation::PositionInTissue&>(a)==static_cast<const Simulation::PositionInTissue&>(b));
}

bool operator==(const Races::Drivers::SpeciesProperties& a, const Races::Drivers::SpeciesProperties& b)
{
    return (a.get_name()==b.get_name() &&
            a.get_id()==b.get_id() &&
            a.get_genotype_id()==b.get_genotype_id() &&
            a.get_rates()==b.get_rates() &&
            a.get_epigenetic_switch_rates()==b.get_epigenetic_switch_rates());
}

inline bool operator!=(const Races::Drivers::SpeciesProperties& a, const Races::Drivers::SpeciesProperties& b)
{
    return !(a==b);
}

bool operator==(const Races::Drivers::Simulation::Species& a, const Races::Drivers::Simulation::Species& b)
{
    using namespace Races::Drivers;

    if (static_cast<const SpeciesProperties&>(a)!=static_cast<const SpeciesProperties&>(b)) {
        return false;
    }

    auto cell_a_it = a.begin(), cell_b_it = b.begin();

    for (; cell_a_it != a.end(); ++cell_a_it, ++cell_b_it) {
        if (*cell_a_it != *cell_b_it) {
            return false;
        }
    }

    return true;
}

bool operator==(const Races::Drivers::Simulation::Tissue& a, const Races::Drivers::Simulation::Tissue& b)
{
    using namespace Races::Drivers;

    if (a.get_name()!=b.get_name()) {
        return false;
    }

    std::set<GenotypeId> genotype_ids;

    auto a_it = a.begin(), b_it = b.begin();

    for (; a_it != a.end(); ++a_it, ++b_it) {

        if (*a_it != *b_it) {
            return false;
        }
        genotype_ids.insert(a_it->get_genotype_id());
    }

    for (const auto& genotype_id : genotype_ids) {
        auto a_view = a.get_genotype_species(genotype_id);
        auto b_view = b.get_genotype_species(genotype_id);

        auto a_it = a_view.begin(), b_it = b_view.begin();

        for (; a_it != a_view.end(); ++a_it, ++b_it) {

            if (a_it->get_id() != b_it->get_id()) {
                return false;
            }
        }
    }

    return true;
}

bool operator==(const Races::Drivers::Simulation::Simulation& a, const Races::Drivers::Simulation::Simulation& b)
{
    return a.get_time()==b.get_time() &&
           a.tissue()==b.tissue() &&
           a.death_activation_level==b.death_activation_level;
}

BOOST_AUTO_TEST_CASE(binary_size_t)
{
    basic_type_test<size_t>({3, 4, 0, 5});
}

BOOST_AUTO_TEST_CASE(binary_double)
{
    basic_type_test<double>({-3.3, 1/4, -0, 2/3});
}

BOOST_AUTO_TEST_CASE(binary_char)
{
    basic_type_test<char>({'a', '.', '\n', '\0'});
}

BOOST_AUTO_TEST_CASE(binary_string)
{
    basic_type_test<std::string>({"string", "", "\n", "\"ci\0ao\""});
}

BOOST_AUTO_TEST_CASE(binary_vector)
{
    basic_type_test<std::vector<size_t>>({{},{3, 4}, {}, {0}, {5}});
    basic_type_test<std::vector<double>>({{-3.3}, {}, {1/4, -0},{}, {2/3}});
    basic_type_test<std::vector<char>>({{},{},{'a', '.'}, {}, {'\n', '\0'},{}});
    basic_type_test<std::vector<std::string>>({{},{"string", ""}, {}, {"\n"}, {"\"ci\0ao\""}});
}

BOOST_AUTO_TEST_CASE(binary_list)
{
    basic_type_test<std::list<size_t>>({{},{3, 4}, {}, {0}, {5}});
    basic_type_test<std::list<double>>({{-3.3}, {}, {1/4, -0},{}, {2/3}});
    basic_type_test<std::list<char>>({{},{},{'a', '.'}, {}, {'\n', '\0'},{}});
    basic_type_test<std::list<std::string>>({{},{"string", ""}, {}, {"\n"}, {"\"ci\0ao\""}});
}

BOOST_AUTO_TEST_CASE(binary_map)
{
    std::map<std::string, std::vector<std::vector<double>>> to_save{{"ciao", {{}, {-0.3}, {5,6}}},
                                                                    {"", {}},
                                                                    {"ult\0imo", {{-3/7},{}}}};
    auto filename = get_a_temporary_path();
    {
        Races::Archive::Binary::Out o_archive(filename);

        o_archive & to_save;
    }

    {
        Races::Archive::Binary::In i_archive(filename);

        std::map<std::string, std::vector<std::vector<double>>> read_value;

        i_archive & read_value;

        BOOST_CHECK(read_value==to_save);
    }

    std::filesystem::remove(filename);
}

BOOST_AUTO_TEST_CASE(binary_cell)
{
    using namespace Races::Drivers;

    std::vector<Cell> to_save{Cell(0), Cell(1,300), Cell(2, 200)};

    auto filename = get_a_temporary_path();
    {
        Races::Archive::Binary::Out o_archive(filename);

        for (const auto& value : to_save) {
            o_archive & value;
        }
    }

    {
        Races::Archive::Binary::In i_archive(filename);

        for (auto& value : to_save) {
            Cell read_value = Cell::load(i_archive);

            BOOST_CHECK(read_value.get_id()==value.get_id());
            BOOST_CHECK(read_value.get_parent_id()==value.get_parent_id());
            BOOST_CHECK(read_value.get_species_id()==value.get_species_id());
        }
    }

    std::filesystem::remove(filename);
}

BOOST_AUTO_TEST_CASE(binary_timed_genotype_mutation)
{
    using namespace Races::Drivers::Simulation;

    std::vector<TimedEvent> to_save{{5,SimulationEventWrapper(GenotypeMutation(0,1))},
                                    {3.5,SimulationEventWrapper(GenotypeMutation(1,7))},
                                    {8.1,SimulationEventWrapper(GenotypeMutation(2,1))}};

    auto filename = get_a_temporary_path();
    {
        Races::Archive::Binary::Out o_archive(filename);

        for (const auto& value : to_save) {
            o_archive & value;
        }
    }

    {
        Races::Archive::Binary::In i_archive(filename);

        for (auto& value : to_save) {
            TimedEvent read_value = TimedEvent::load(i_archive);

            BOOST_CHECK(read_value==value);
        }
    }

    std::filesystem::remove(filename);
}


BOOST_AUTO_TEST_CASE(binary_timed_genotype_mutation_queue)
{
    using namespace Races::Drivers::Simulation;

    std::vector<TimedEvent> to_save{{5,SimulationEventWrapper(GenotypeMutation(0,1))},
                                    {3.5,SimulationEventWrapper(GenotypeMutation(1,7))},
                                    {8.1,SimulationEventWrapper(GenotypeMutation(2,1))}};

    using PriorityQueue = std::priority_queue<TimedEvent,
                                              std::vector<TimedEvent>,
                                              std::greater<TimedEvent>>;

    PriorityQueue queue(to_save.begin(), to_save.end());

    auto filename = get_a_temporary_path();
    {
        Races::Archive::Binary::Out o_archive(filename);

        o_archive & queue;
    }

    {
        Races::Archive::Binary::In i_archive(filename);
    
        PriorityQueue i_queue;

        i_archive & i_queue;
    
        BOOST_CHECK(i_queue.size()==queue.size());

        while (!i_queue.empty() && !queue.empty()) {

            BOOST_CHECK(i_queue.top()==queue.top());
            i_queue.pop();
            queue.pop();
        }

        BOOST_CHECK(i_queue.size()==queue.size());
    }

    std::filesystem::remove(filename);
}

BOOST_AUTO_TEST_CASE(binary_epigenetic_genotype)
{
    using namespace Races::Drivers;
    
    GenotypeProperties to_save("A",{{0.01,0.01},{0.01,0.01}});
    to_save["--"].set_rates({{CellEventType::DEATH, 0.1},
                             {CellEventType::DUPLICATION, 0.3}});
    to_save["+-"].set_rates({{CellEventType::DEATH, 0.1},
                             {CellEventType::DUPLICATION, 0.45}});
    to_save["-+"].set_rates({{CellEventType::DEATH, 0.1},
                             {CellEventType::DUPLICATION, 0.2}});
    to_save["+-"].set_rates({{CellEventType::DEATH, 0.01},
                             {CellEventType::DUPLICATION, 0.02}});

    auto filename = get_a_temporary_path();
    {
        Races::Archive::Binary::Out o_archive(filename);

        for (const auto& species : to_save.get_species()) {
            o_archive & species;
        }
    }

    {
        Races::Archive::Binary::In i_archive(filename);

        for (const auto& species : to_save.get_species()) {
            SpeciesProperties in_species =  SpeciesProperties::load(i_archive);

            BOOST_CHECK(species==in_species);
        }
    }
    std::filesystem::remove(filename);
}

BOOST_FIXTURE_TEST_SUITE( simulatedData, ArchiveFixture )

BOOST_AUTO_TEST_CASE(binary_tissue)
{
    using namespace Races::Drivers::Simulation;

    auto filename = get_a_temporary_path();

    {
        Races::Archive::Binary::Out o_archive(filename);

        o_archive & simulation.tissue();
    }

    {
        Races::Archive::Binary::In i_archive(filename);

        Tissue in_tissue = Tissue::load(i_archive);

        BOOST_CHECK(simulation.tissue()==in_tissue);
    }

    std::filesystem::remove(filename);
}

BOOST_AUTO_TEST_CASE(simulation_statistics)
{
    using namespace Races::Drivers::Simulation;

    auto filename = get_a_temporary_path();

    {
        Races::Archive::Binary::Out o_archive(filename);

        o_archive & simulation.get_statistics();
    }

    {
        Races::Archive::Binary::In i_archive(filename);

        TissueStatistics statistics = TissueStatistics::load(i_archive);
    }

    std::filesystem::remove(filename);
}

BOOST_AUTO_TEST_CASE(simulation_tissue)
{
    using namespace Races::Drivers::Simulation;

    auto filename = get_a_temporary_path();

    {
        Races::Archive::Binary::Out o_archive(filename);

        o_archive & simulation;
    }

    {
        Races::Archive::Binary::In i_archive(filename);

        Simulation in_simulation = Simulation::load(i_archive);

        BOOST_CHECK(simulation==in_simulation);
    }

    std::filesystem::remove(filename);
}

BOOST_AUTO_TEST_SUITE_END()