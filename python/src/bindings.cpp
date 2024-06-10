/**
 * @file bindings.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements Python bindings
 * @version 0.18
 * @date 2024-06-10
 *
 * @copyright Copyright (c) 2023-2024
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
#include "mutant_properties.hpp"

#include "simulation_wrapper.hpp"
#include "mutant.hpp"
#include "epigenetic_rates.hpp"

using namespace boost::python;

namespace RACESSim = RACES::Mutants::Evolutions;
namespace RACESDrv = RACES::Mutants;


BOOST_PYTHON_FUNCTION_OVERLOADS(run_up_to_overloads, SimulationWrapper::static_run_up_to, 2, 4)

BOOST_PYTHON_MODULE(RACES)
{
    class_<RACESSim::PositionInTissue>("Position", init<RACESSim::AxisPosition, RACESSim::AxisPosition>())
        .def(init<RACESSim::AxisPosition, RACESSim::AxisPosition, RACESSim::AxisPosition>())
        .def(self_ns::str(self))
        .def_readwrite("x", &RACESSim::PositionInTissue::x)
        .def_readwrite("y", &RACESSim::PositionInTissue::y)
        .def_readwrite("z", &RACESSim::PositionInTissue::z)
    ;

    enum_<RACESDrv::CellEventType>("CellEventType")
        .value("DEATH", RACESDrv::CellEventType::DEATH)
        .value("DUPLICATION", RACESDrv::CellEventType::DUPLICATION)
        .value("EPIGENETIC_SWITCH", RACESDrv::CellEventType::EPIGENETIC_SWITCH)
        .value("MUTATION", RACESDrv::CellEventType::MUTATION)
    ;

    class_<RACESDrv::EpigeneticRates, std::shared_ptr<RACESDrv::EpigeneticRates>>("EpigeneticRates", init<double, double>())
        .def("__init__", make_constructor(EpigeneticRatesWrapper::create))
        .def(self_ns::str(self))
        .def("get_methylation_rate", make_function(&RACESDrv::EpigeneticRates::get_methylation_rate, return_value_policy<copy_const_reference>()))
        .def("set_methylation_rate", &EpigeneticRatesWrapper::set_methylation_rate)
        .def("get_demethylation_rate", make_function(&RACESDrv::EpigeneticRates::get_demethylation_rate, return_value_policy<copy_const_reference>()))
        .def("set_demethylation_rate", &EpigeneticRatesWrapper::set_demethylation_rate)
    ;

    class_<RACESDrv::MutantProperties, std::shared_ptr<RACESDrv::MutantProperties>>("Clone", no_init)
        .def("__init__", make_constructor(CloneWrapper::create))
        .def("set_rates", &CloneWrapper::set_rates)
        .def("get_rate", &CloneWrapper::get_rate)
        .add_property("num_of_promoters", &RACESDrv::MutantProperties::num_of_promoters)
        .add_property("name", make_function(&RACESDrv::MutantProperties::get_name, return_value_policy<copy_const_reference>()))
        .add_property("id", make_function(&RACESDrv::MutantProperties::get_id, return_value_policy<copy_const_reference>()))
    ;

    class_<SimulationWrapper, std::shared_ptr<SimulationWrapper>>("Simulation", no_init)
        .def("__init__", make_constructor(SimulationWrapper::create, default_call_policies(),
                                          (arg("minutes_between_snapshot")=5, arg("random_seed")=0)))
        .def("run_up_to", &SimulationWrapper::static_run_up_to,
             run_up_to_overloads((arg("wrapper"), arg("time"), arg("quiet")=false, arg("plot")=false)))
        .def("get_time", make_function(&SimulationWrapper::get_time, return_value_policy<copy_const_reference>()))
        .def("add_mutant", &SimulationWrapper::add_mutant)
        .def("schedule_mutation", &SimulationWrapper::schedule_mutation)
        .def("place_cell", &SimulationWrapper::place_cell)
        .def("set_tissue", &SimulationWrapper::set_tissue)
        .def("rename_log_directory", &SimulationWrapper::rename_log_directory)
        .add_property("death_activation_level", &SimulationWrapper::get_death_activation_level,
                      &SimulationWrapper::set_death_activation_level)
        .add_property("storage_enabled", &SimulationWrapper::get_storage_enabled,
                      &SimulationWrapper::set_storage_enabled)
    ;
}
