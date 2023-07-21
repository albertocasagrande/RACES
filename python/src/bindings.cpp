/**
 * @file bindings.cpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Implement Python bindings
 * @version 0.2
 * @date 2023-07-21
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

#include <memory>

#include <boost/python.hpp>

#include "position.hpp" 
#include "driver_genotype.hpp"

#include "simulation_wrapper.hpp"
#include "genotype.hpp"
#include "epigenetic_rates.hpp"

using namespace boost::python;

namespace RacesSim = Races::Drivers::Simulation;
namespace RacesDrv = Races::Drivers;


BOOST_PYTHON_FUNCTION_OVERLOADS(run_up_to_overloads, SimulationWrapper::static_run_up_to, 2, 5)

BOOST_PYTHON_MODULE(RACES)
{
    class_<RacesSim::PositionInTissue>("Position", init<RacesSim::AxisValue, RacesSim::AxisValue>())
        .def(init<RacesSim::AxisValue, RacesSim::AxisValue, RacesSim::AxisValue>())
        .def(self_ns::str(self))
        .def_readwrite("x", &RacesSim::PositionInTissue::x)
        .def_readwrite("y", &RacesSim::PositionInTissue::y)
        .def_readwrite("z", &RacesSim::PositionInTissue::z)
    ;

    enum_<RacesDrv::CellEventType>("CellEventType")
        .value("DIE", RacesDrv::CellEventType::DIE)
        .value("DUPLICATE", RacesDrv::CellEventType::DUPLICATE)
        .value("EPIGENETIC_EVENT", RacesDrv::CellEventType::EPIGENETIC_EVENT)
        .value("DUPLICATION_AND_EPIGENETIC_EVENT", RacesDrv::CellEventType::DUPLICATION_AND_EPIGENETIC_EVENT)
        .value("DRIVER_GENETIC_MUTATION", RacesDrv::CellEventType::DRIVER_GENETIC_MUTATION)
    ;

    class_<RacesDrv::EpigeneticRates, std::shared_ptr<RacesDrv::EpigeneticRates>>("EpigeneticRates", init<double, double>())
        .def("__init__", make_constructor(EpigeneticRatesWrapper::create))
        .def(self_ns::str(self))
        .def("get_methylation_rate", make_function(&RacesDrv::EpigeneticRates::get_methylation_rate, return_value_policy<copy_const_reference>()))
        .def("set_methylation_rate", &EpigeneticRatesWrapper::set_methylation_rate)
        .def("get_demethylation_rate", make_function(&RacesDrv::EpigeneticRates::get_demethylation_rate, return_value_policy<copy_const_reference>()))
        .def("set_demethylation_rate", &EpigeneticRatesWrapper::set_demethylation_rate)
    ;

    class_<RacesDrv::Genotype, std::shared_ptr<RacesDrv::Genotype>>("DriverGenotype", no_init)
        .def("__init__", make_constructor(GenotypeWrapper::create))
        .def("set_rates", &GenotypeWrapper::set_rates)
        .def("get_rate", &GenotypeWrapper::get_rate)
        .add_property("num_of_promoters", &RacesDrv::Genotype::num_of_promoters)
        .add_property("name", make_function(&RacesDrv::Genotype::get_name, return_value_policy<copy_const_reference>()))
        .add_property("id", make_function(&RacesDrv::Genotype::get_id, return_value_policy<copy_const_reference>()))
    ;

    class_<SimulationWrapper, std::shared_ptr<SimulationWrapper>>("Simulation", no_init)
        .def("__init__", make_constructor(SimulationWrapper::create, default_call_policies(),
                                          (arg("minutes_between_snapshot")=5, arg("random_seed")=0)))
        .def("run_up_to", &SimulationWrapper::static_run_up_to,
             run_up_to_overloads((arg("wrapper"), arg("time"), arg("logging")=true,
                                  arg("quiet")=false, arg("plot")=false)))
        .def("get_time", make_function(&SimulationWrapper::get_time, return_value_policy<copy_const_reference>()))
        .def("add_species", &SimulationWrapper::add_species)
        .def("add_genomic_mutation", &SimulationWrapper::add_genomic_mutation)
        .def("add_cell", &SimulationWrapper::add_cell)
        .def("set_tissue", &SimulationWrapper::set_tissue)
        .add_property("death_activation_level", &SimulationWrapper::get_death_activation_level, 
                      &SimulationWrapper::set_death_activation_level)
    ;
}