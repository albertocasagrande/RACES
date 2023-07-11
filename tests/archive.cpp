/**
 * @file archive.cpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Some archive tests
 * @version 0.3
 * @date 2023-07-11
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
#include "simulator.hpp"

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

bool operator==(const Races::Cell& a, const Races::Cell& b)
{
    return (a.get_id()==b.get_id() && 
            a.get_parent_id()==b.get_parent_id() &&
            a.get_genotype_id()==b.get_genotype_id());
}

bool operator==(const Races::CellInTissue& a, const Races::CellInTissue& b)
{
    return (static_cast<const Races::Cell&>(a)==static_cast<const Races::Cell&>(b) &&
            static_cast<const Races::PositionInTissue&>(a)==static_cast<const Races::PositionInTissue&>(b));
}

bool operator==(const Races::EpigeneticGenotype& a, const Races::EpigeneticGenotype& b)
{
    return (a.get_name()==b.get_name() &&
            a.get_id()==b.get_id() &&
            a.get_somatic_id()==b.get_somatic_id() &&
            a.get_rates()==b.get_rates() &&
            a.get_epigenentic_mutation_rates()==b.get_epigenentic_mutation_rates());
}

inline bool operator!=(const Races::EpigeneticGenotype& a, const Races::EpigeneticGenotype& b)
{
    return !(a==b);
}


bool operator==(const Races::Species& a, const Races::Species& b)
{
    if (static_cast<const Races::EpigeneticGenotype&>(a)!=static_cast<const Races::EpigeneticGenotype&>(b)) {
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

bool operator==(const Races::Tissue& a, const Races::Tissue& b)
{
    if (a.get_name()!=b.get_name()) {
        return false;
    }

    std::set<Races::SomaticGenotypeId> somatic_ids;

    auto a_it = a.begin(), b_it = b.begin();

    for (; a_it != a.end(); ++a_it, ++b_it) {

        if (*a_it != *b_it) {
            return false;
        }
        somatic_ids.insert(a_it->get_somatic_id());
    }

    for (const auto& somatic_id : somatic_ids) {
        auto a_view = a.get_somatic_genotype_species(somatic_id);
        auto b_view = b.get_somatic_genotype_species(somatic_id);

        auto a_it = a_view.begin(), b_it = b_view.begin();

        for (; a_it != a_view.end(); ++a_it, ++b_it) {

            if (a_it->get_id() != b_it->get_id()) {
                return false;
            }
        }
    }

    return true;
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
    std::vector<Races::Cell> to_save{{0}, {1,300}, {2, 200}};

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
            Races::Cell read_value = Races::Cell::load(i_archive);

            BOOST_CHECK(read_value.get_id()==value.get_id());
            BOOST_CHECK(read_value.get_parent_id()==value.get_parent_id());
            BOOST_CHECK(read_value.get_genotype_id()==value.get_genotype_id());
        }
    }

    std::filesystem::remove(filename);
}

BOOST_AUTO_TEST_CASE(binary_epigenetic_genotype)
{
    using namespace Races;
    
    SomaticGenotype to_save("A",{{0.01,0.01},{0.01,0.01}});
    to_save["--"].set_rates({{CellEventType::DIE, 0.1},
                             {CellEventType::DUPLICATE, 0.3},
                             {CellEventType::PASSENGER_MUTATION, 20}});
    to_save["+-"].set_rates({{CellEventType::DIE, 0.1},
                             {CellEventType::DUPLICATE, 0.45},
                             {CellEventType::PASSENGER_MUTATION, 2}});
    to_save["-+"].set_rates({{CellEventType::DIE, 0.1},
                             {CellEventType::DUPLICATE, 0.2},
                             {CellEventType::PASSENGER_MUTATION, 1}});
    to_save["+-"].set_rates({{CellEventType::DIE, 0.01},
                             {CellEventType::DUPLICATE, 0.02},
                             {CellEventType::PASSENGER_MUTATION, 20}});

    auto filename = get_a_temporary_path();
    {
        Races::Archive::Binary::Out o_archive(filename);

        for (const auto& e_genotype : to_save.epigenetic_genotypes()) {
            o_archive & e_genotype;
        }
    }

    {
        Races::Archive::Binary::In i_archive(filename);

        for (const auto& e_genotype : to_save.epigenetic_genotypes()) {
            Races::EpigeneticGenotype in_genotype =  Races::EpigeneticGenotype::load(i_archive);

            BOOST_CHECK(e_genotype==in_genotype);
        }
    }
    std::filesystem::remove(filename);
}

BOOST_AUTO_TEST_CASE(binary_tissue)
{
    using namespace Races;
    
    long double time_horizon = 70;

    SomaticGenotype A("A",{{0.01,0.01}});
    A["-"].set_rates({{CellEventType::DIE, 0.1},
                      {CellEventType::DUPLICATE, 0.3},
                      {CellEventType::PASSENGER_MUTATION, 20}});
    A["+"].set_rates({{CellEventType::DIE, 0.1},
                      {CellEventType::DUPLICATE, 0.45},
                      {CellEventType::PASSENGER_MUTATION, 2}});

    SomaticGenotype B("B",{{0.01,0.01}});
    B["-"].set_rates({{CellEventType::DIE, 0.1},
                      {CellEventType::DUPLICATE, 0.2},
                      {CellEventType::PASSENGER_MUTATION, 1}});
    B["+"].set_rates({{CellEventType::DIE, 0.01},
                      {CellEventType::DUPLICATE, 0.02},
                      {CellEventType::PASSENGER_MUTATION, 20}});

    Tissue tissue("Liver", 1000,1000);

    tissue.add_species(A);
    tissue.add_species(B);

    tissue.add(B["-"].get_id(), {250, 500});
    {
        BasicSimulator<> simulator(tissue);

        simulator.death_activation_level = 100;

        simulator.add_somatic_mutation(B,A,50);

        simulator.run_up_to(time_horizon);
    }

    auto filename = get_a_temporary_path();

    {
        Races::Archive::Binary::Out o_archive(filename);

        o_archive & tissue;
    }

    {
        Races::Archive::Binary::In i_archive(filename);

        Tissue in_tissue = Tissue::load(i_archive);

        BOOST_CHECK(tissue==in_tissue);
    }

    std::filesystem::remove(filename);
}