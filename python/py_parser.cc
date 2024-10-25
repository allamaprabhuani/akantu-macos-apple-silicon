/**
 * Copyright (©) 2020-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * This file is part of Akantu
 *
 * Akantu is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 */

/* -------------------------------------------------------------------------- */
#include "py_aka_array.hh"
/* -------------------------------------------------------------------------- */
#include <aka_common.hh>
#include <parameter_registry.hh>
#include <parsable.hh>
#include <parser.hh>
/* -------------------------------------------------------------------------- */
#include <map>
#include <pybind11/pybind11.h>
/* -------------------------------------------------------------------------- */
namespace py = pybind11;
/* -------------------------------------------------------------------------- */

namespace akantu {
std::map<void *, std::map<std::string, void *>> map_params;

void register_parser(py::module & mod) {
  py::enum_<ParameterAccessType>(mod, "ParameterAccessType", py::arithmetic())
      .value("_pat_internal", _pat_internal)
      .value("_pat_writable", _pat_writable)
      .value("_pat_readable", _pat_readable)
      .value("_pat_modifiable", _pat_modifiable)
      .value("_pat_parsable", _pat_parsable)
      .value("_pat_parsmod", _pat_parsmod)
      .export_values();

  py::class_<ParameterRegistry>(mod, "ParameterRegistry",
                                py::multiple_inheritance())
      .def("registerParamReal",
           [](ParameterRegistry & self, const std::string & name, UInt type,
              const std::string & description) {
             Real * p = new Real;
             map_params[&self][name] = p;
             self.registerParam<Real>(name, *p, ParameterAccessType(type),
                                      description);
           })
      .def("registerParamReal",
           [](ParameterRegistry & self, const Real & _default,
              const std::string & name, UInt type,
              const std::string & description) {
             Real * p = new Real;
             map_params[&self][name] = p;
             self.registerParam<Real>(name, *p, _default,
                                      ParameterAccessType(type), description);
           })

      .def("setReal", [](ParameterRegistry & self, const std::string & name,
                         const Real value) { self.set(name, value); })
      .def("getReal",
           [](ParameterRegistry & self, const std::string & name) -> Real {
             return self.get(name);
           })
      .def("setRandom",
           [](ParameterRegistry & self, const std::string & name, Real value,
              const RandomDistributionType & distribution, Real a, Real b) {
             auto param = make_random_parameter(value, distribution, a, b);
             self.set(name, value);
           })

      .def("setBool", [](ParameterRegistry & self, const std::string & name,
                         const bool value) { self.set(name, value); })
      .def("getBool",
           [](ParameterRegistry & self, const std::string & name) -> bool {
             return self.get(name);
           })

      .def("setString",
           [](ParameterRegistry & self, const std::string & name,
              const std::string & value) { self.set(name, value); })
      .def("getString",
           [](ParameterRegistry & self,
              const std::string & name) -> std::string {
             std::string tmp = self.get(name);
             return tmp;
           })

      .def("setInt", [](ParameterRegistry & self, const std::string & name,
                        const Int value) { self.set(name, value); })
      .def("getInt",
           [](ParameterRegistry & self, const std::string & name) -> Int {
             return self.get(name);
           })

      .def(
          "getMatrix",
          [](ParameterRegistry & self, const std::string & name) {
            const Matrix<Real> & res =
                static_cast<const Matrix<Real> &>(self.get(name));
            return res;
          },
          py::return_value_policy::copy);

  py::class_<Parsable, ParameterRegistry>(mod, "Parsable",
                                          py::multiple_inheritance())
      .def(py::init<const ParserType &, const ID &>())
      .def("set", [](Parsable & self, std::string param, std::string value) {
        std::any val = value.c_str();
        self.set(param, val);
      });

  mod.def(
      "parseInput",
      [](const std::string & input_file) {
        getStaticParser().parse(input_file);
      },
      "Parse an Akantu input file");
}
} // namespace akantu
