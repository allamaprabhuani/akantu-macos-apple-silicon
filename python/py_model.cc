/**
 * Copyright (©) 2019-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include <constitutive_laws_handler.hh>
#include <model.hh>
#include <non_linear_solver.hh>
#include <solver_callback.hh>
#include <sparse_matrix_aij.hh>
#include <time_step_solver.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
namespace py = pybind11;
/* -------------------------------------------------------------------------- */

namespace akantu {

template <class ConstitutiveLawType, class ModelType>
void register_constitutive_law_handler(py::module & mod,
                                       const std::string & name) {
  py::class_<ConstitutiveLawsHandler<ConstitutiveLawType, ModelType>,
             ModelType>(mod, name.c_str(), py::multiple_inheritance())
      .def(
          "getConstitutiveLaw",
          [](ConstitutiveLawsHandler<ConstitutiveLawType, ModelType> & self,
             UInt constitutive_law_id) -> decltype(auto) {
            return self.getConstitutiveLaw(constitutive_law_id);
          },
          py::arg("constitutive_law_id"), py::return_value_policy::reference)
      .def(
          "getConstitutiveLaw",
          [](ConstitutiveLawsHandler<ConstitutiveLawType, ModelType> & self,
             const ID & constitutive_law_name) -> decltype(auto) {
            return self.getConstitutiveLaw(constitutive_law_name);
          },
          py::arg("constitutive_law_name"), py::return_value_policy::reference)
      .def(
          "getConstitutiveLaw",
          [](ConstitutiveLawsHandler<ConstitutiveLawType, ModelType> & self,
             const Element & element) -> decltype(auto) {
            return self.getConstitutiveLaw(element);
          },
          py::arg("element"), py::return_value_policy::reference)

      .def("getNbConstitutiveLaws",
           &ConstitutiveLawsHandler<ConstitutiveLawType,
                                    ModelType>::getNbConstitutiveLaws)
      .def("getConstitutiveLawIndex",
           &ConstitutiveLawsHandler<ConstitutiveLawType,
                                    ModelType>::getConstitutiveLawIndex)
      .def("setConstitutiveLawSelector",
           [](ConstitutiveLawsHandler<ConstitutiveLawType, ModelType> & self,
              std::shared_ptr<ConstitutiveLawSelector>
                  constitutive_law_selector) {
             self.setConstitutiveLawSelector(constitutive_law_selector);
           })
      .def("getConstitutiveLawSelector",
           &ConstitutiveLawsHandler<ConstitutiveLawType,
                                    ModelType>::getConstitutiveLawSelector)
      .def(
          "getConstitutiveLawByElement",
          [](const ConstitutiveLawsHandler<ConstitutiveLawType, ModelType> &
                 self) -> decltype(auto) {
            return self.getConstitutiveLawByElement();
          },
          py::return_value_policy::reference, py::keep_alive<0, 1>())
      .def("reassignConstitutiveLaw",
           &ConstitutiveLawsHandler<ConstitutiveLawType,
                                    ModelType>::reassignConstitutiveLaw)
      .def(
          "registerNewConstitutiveLaw",
          [](ConstitutiveLawsHandler<ConstitutiveLawType, ModelType> & self,
             const ID & cl_name, const ID & cl_type,
             const ID & opt_param) -> decltype(auto) {
            return self.registerNewConstitutiveLaw(cl_name, cl_type, opt_param);
          },
          py::arg("constitutive_law_name"), py::arg("constitutive_law_type"),
          py::arg("option") = "", py::return_value_policy::reference)
      .def("initConstitutiveLaws",
           &ConstitutiveLawsHandler<ConstitutiveLawType,
                                    ModelType>::initConstitutiveLaws)
      .def("flattenInternal",
           &ConstitutiveLawsHandler<ConstitutiveLawType,
                                    ModelType>::flattenInternal,
           py::return_value_policy::reference)
      .def("inflateInternal",
           &ConstitutiveLawsHandler<ConstitutiveLawType,
                                    ModelType>::inflateInternal,
           py::return_value_policy::reference);
}

/* -------------------------------------------------------------------------- */
void register_model(py::module & mod) {
  py::class_<ModelSolver, SolverCallback, Parsable>(mod, "ModelSolver",
                                                    py::multiple_inheritance())
      .def(
          "getNonLinearSolver",
          [](ModelSolver & self, const ID & solver_id) -> NonLinearSolver & {
            return self.getNonLinearSolver(solver_id);
          },
          py::arg("solver_id") = "", py::return_value_policy::reference)
      .def(
          "getTimeStepSolver",
          [](ModelSolver & self, const ID & solver_id) -> TimeStepSolver & {
            return self.getTimeStepSolver(solver_id);
          },
          py::arg("solver_id") = "", py::return_value_policy::reference)
      .def(
          "solveStep",
          [](ModelSolver & self, const ID & solver_id) {
            self.solveStep(solver_id);
          },
          py::arg("solver_id") = "")
      .def(
          "solveStep",
          [](ModelSolver & self, SolverCallback & callback,
             const ID & solver_id) { self.solveStep(callback, solver_id); },
          py::arg("callback"), py::arg("solver_id") = "");

  py::class_<Model, ModelSolver>(mod, "Model", py::multiple_inheritance())
      .def("setBaseName", &Model::setBaseName)
      .def("setDirectory", &Model::setDirectory)
      .def("getFEEngine", &Model::getFEEngine, py::arg("name") = "",
           py::return_value_policy::reference)
      .def("getFEEngineBoundary", &Model::getFEEngine, py::arg("name") = "",
           py::return_value_policy::reference)
      .def("addDumpFieldVector", &Model::addDumpFieldVector)
      .def("addDumpField", &Model::addDumpField)
      .def("setBaseNameToDumper", &Model::setBaseNameToDumper)
      .def("addDumpFieldVectorToDumper", &Model::addDumpFieldVectorToDumper)
      .def("addDumpFieldToDumper", &Model::addDumpFieldToDumper)
      .def("dump", [](Model & self) { self.dump(); })
      .def(
          "dump", [](Model & self, UInt step) { self.dump(step); },
          py::arg("step"))
      .def(
          "dump",
          [](Model & self, Real time, UInt step) { self.dump(time, step); },
          py::arg("time"), py::arg("step"))
      .def(
          "dump",
          [](Model & self, const std::string & dumper) { self.dump(dumper); },
          py::arg("dumper_name"))
      .def(
          "dump",
          [](Model & self, const std::string & dumper, UInt step) {
            self.dump(dumper, step);
          },
          py::arg("dumper_name"), py::arg("step"))
      .def(
          "dump",
          [](Model & self, const std::string & dumper, Real time, UInt step) {
            self.dump(dumper, time, step);
          },
          py::arg("dumper_name"), py::arg("time"), py::arg("step"))

      .def(
          "dumpGroup",
          [](Model & self, const std::string & group_name) {
            self.dumpGroup(group_name);
          },
          py::arg("group_name"))
      .def(
          "setGroupDirectory",
          [](Model & self, const std::string & directory,
             const std::string & group_name) {
            self.setGroupDirectory(directory, group_name);
          },
          py::arg("directory"), py::arg("group_name"))
      .def(
          "setGroupBaseName",
          [](Model & self, const std::string & basename,
             const std::string & group_name) {
            self.setGroupBaseName(basename, group_name);
          },
          py::arg("basename"), py::arg("group_name"))
      .def(
          "addDumpGroupField",
          [](Model & self, const std::string & field_id,
             const std::string & group_name) {
            self.addDumpGroupField(field_id, group_name);
          },
          py::arg("field_id"), py::arg("group_name"))

      .def(
          "addDumpGroupFieldVector",
          [](Model & self, const std::string & field_id,
             const std::string & group_name) {
            self.addDumpGroupFieldVector(field_id, group_name);
          },
          py::arg("field_id"), py::arg("group_name"))
      .def("initNewSolver", &Model::initNewSolver)
      .def(
          "getNewSolver",
          [](Model & self, const std::string id,
             ModelSolverOptions & solver_options) {
            self.getNewSolver(id, solver_options);
          },
          py::return_value_policy::reference)
      .def(
          "setIntegrationScheme",
          [](Model & self, const std::string id, const std::string primal,
             const IntegrationSchemeType & scheme_type,
             IntegrationScheme::SolutionType solution_type) {
            self.setIntegrationScheme(id, primal, scheme_type, solution_type);
          },
          py::arg("id"), py::arg("primal"), py::arg("scheme_type"),
          py::arg("solution_type") =
              IntegrationScheme::SolutionType::_not_defined)
      .def("getDOFManager", &Model::getDOFManager,
           py::return_value_policy::reference)
      .def("assembleMatrix", &Model::assembleMatrix);

  py::class_<ModelSolverOptions>(mod, "ModelSolverOptions").def(py::init());
  py::class_<ModelOptions>(mod, "ModelOptions")
      .def(py::init<AnalysisMethod, ModelSolverOptions>(),
           py::arg("analysis_method") = _explicit_lumped_mass,
           py::arg("solver_options") = detail::_default_solver_options);
}

} // namespace akantu
