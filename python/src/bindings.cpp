/**
 * @file bindings.cpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Implement Python bindings
 * @version 0.1
 * @date 2023-07-19
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
#include "somatic_genotype.hpp"
#include "epigenetic_rates.hpp"

using namespace boost::python;


BOOST_PYTHON_FUNCTION_OVERLOADS(run_up_to_overloads, SimulationWrapper::static_run_up_to, 2, 5)

BOOST_PYTHON_MODULE(RACES)
{
    class_<Races::PositionInTissue>("Position", init<Races::AxisValue, Races::AxisValue>())
        .def(init<Races::AxisValue, Races::AxisValue, Races::AxisValue>())
        .def(self_ns::str(self))
        .def_readwrite("x", &Races::PositionInTissue::x)
        .def_readwrite("y", &Races::PositionInTissue::y)
        .def_readwrite("z", &Races::PositionInTissue::z)
    ;

    enum_<Races::CellEventType>("CellEventType")
        .value("DIE", Races::CellEventType::DIE)
        .value("DUPLICATE", Races::CellEventType::DUPLICATE)
        .value("EPIGENETIC_EVENT", Races::CellEventType::EPIGENETIC_EVENT)
        .value("DUPLICATION_AND_EPIGENETIC_EVENT", Races::CellEventType::DUPLICATION_AND_EPIGENETIC_EVENT)
        .value("PASSENGER_MUTATION", Races::CellEventType::PASSENGER_MUTATION)
        .value("DRIVER_SOMATIC_MUTATION", Races::CellEventType::DRIVER_SOMATIC_MUTATION)
    ;

    class_<Races::EpigeneticRates, std::shared_ptr<Races::EpigeneticRates>>("EpigeneticRates", init<double, double>())
        .def("__init__", make_constructor(EpigeneticRatesWrapper::create))
        .def(self_ns::str(self))
        .def("get_methylation_rate", make_function(&Races::EpigeneticRates::get_methylation_rate, return_value_policy<copy_const_reference>()))
        .def("set_methylation_rate", &EpigeneticRatesWrapper::set_methylation_rate)
        .def("get_demethylation_rate", make_function(&Races::EpigeneticRates::get_demethylation_rate, return_value_policy<copy_const_reference>()))
        .def("set_demethylation_rate", &EpigeneticRatesWrapper::set_demethylation_rate)
    ;

    class_<Races::SomaticGenotype, std::shared_ptr<Races::SomaticGenotype>>("SomaticGenotype", no_init)
        .def("__init__", make_constructor(SomaticGenotypeWrapper::create))
        .def("set_rates", &SomaticGenotypeWrapper::set_rates)
        .def("get_rate", &SomaticGenotypeWrapper::get_rate)
        .add_property("num_of_promoters", &Races::SomaticGenotype::num_of_promoters)
        .add_property("name", make_function(&Races::SomaticGenotype::get_name, return_value_policy<copy_const_reference>()))
        .add_property("id", make_function(&Races::SomaticGenotype::get_id, return_value_policy<copy_const_reference>()))
    ;

    class_<SimulationWrapper, std::shared_ptr<SimulationWrapper>>("Simulation", no_init)
        .def("__init__", make_constructor(SimulationWrapper::create, default_call_policies(),
                                          (arg("minutes_between_snapshot")=5, arg("random_seed")=0)))
        .def("run_up_to", &SimulationWrapper::static_run_up_to,
             run_up_to_overloads((arg("wrapper"), arg("time"), arg("logging")=true,
                                  arg("quiet")=false, arg("plot")=false)))
        .def("get_time", make_function(&SimulationWrapper::get_time, return_value_policy<copy_const_reference>()))
        .def("add_species", &SimulationWrapper::add_species)
        .def("add_somatic_mutation", &SimulationWrapper::add_somatic_mutation)
        .def("add_cell", &SimulationWrapper::add_cell)
        .def("set_tissue", &SimulationWrapper::set_tissue)
        .add_property("death_activation_level", &SimulationWrapper::get_death_activation_level, 
                      &SimulationWrapper::set_death_activation_level)
    ;
}