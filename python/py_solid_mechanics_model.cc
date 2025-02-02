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
#include "py_constitutive_laws_handler.hh"
/* -------------------------------------------------------------------------- */
#include <constitutive_law.hh>
#include <constitutive_laws_handler.hh>
#include <material.hh>
#include <non_linear_solver.hh>
#include <solid_mechanics_model.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/pybind11.h>
/* -------------------------------------------------------------------------- */
namespace py = pybind11;
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
#define def_function_nocopy(func_name)                                         \
  def(                                                                         \
      #func_name,                                                              \
      [](SolidMechanicsModel & self) -> decltype(auto) {                       \
        return self.func_name();                                               \
      },                                                                       \
      py::return_value_policy::reference)

#define def_function(func_name)                                                \
  def(#func_name, [](SolidMechanicsModel & self) -> decltype(auto) {           \
    return self.func_name();                                                   \
  })

/* -------------------------------------------------------------------------- */
void register_solid_mechanics_model(py::module & mod) {
  register_constitutive_laws_handler<Material, Model>(mod);

  py::class_<SolidMechanicsModelOptions>(mod, "SolidMechanicsModelOptions")
      .def(py::init<AnalysisMethod>(),
           py::arg("_analysis_method") = _explicit_lumped_mass);

  py::class_<SolidMechanicsModel, ConstitutiveLawsHandler<Material, Model>>(
      mod, "SolidMechanicsModel", py::multiple_inheritance())
      .def(py::init<Mesh &, Int, const ID &, std::shared_ptr<DOFManager>,
                    const ModelType>(),
           py::arg("mesh"), py::arg("spatial_dimension") = _all_dimensions,
           py::arg("id") = "solid_mechanics_model",
           py::arg("dof_manager") = nullptr,
           py::arg("model_type") = ModelType::_solid_mechanics_model)
      .def(
          "initFull",
          [](SolidMechanicsModel & self,
             const SolidMechanicsModelOptions & options) {
            self.initFull(options);
          },
          py::arg("option") = SolidMechanicsModelOptions())
      .def(
          "initFull",
          [](SolidMechanicsModel & self,
             const AnalysisMethod & analysis_method) {
            self.initFull(_analysis_method = analysis_method);
          },
          py::arg("_analysis_method"))
      .def("applyDirichletBC",
           [](SolidMechanicsModel & /*self*/) {
             PyErr_WarnEx(
                 PyExc_DeprecationWarning,
                 "applyDirichletBC() is deprecated, use applyBC instead", 1);
           })
      .def("applyBC",
           [](SolidMechanicsModel & self,
              BC::Dirichlet::DirichletFunctor & func,
              const std::string & element_group) {
             self.applyBC(func, element_group);
           })
      .def("applyBC",
           [](SolidMechanicsModel & self, BC::Neumann::NeumannFunctor & func,
              const std::string & element_group) {
             self.applyBC(func, element_group);
           })
      .def("setTimeStep", &SolidMechanicsModel::setTimeStep,
           py::arg("time_step"), py::arg("solver_id") = "")
      .def(
          "getEnergy",
          [](SolidMechanicsModel & self, const std::string & energy_id) {
            return self.getEnergy(energy_id);
          },
          py::arg("energy_id"))
      .def(
          "getEnergy",
          [](SolidMechanicsModel & self, const std::string & energy_id,
             const std::string & group_id) {
            return self.getEnergy(energy_id, group_id);
          },
          py::arg("energy_id"), py::arg("group_id"))

      .def_function(assembleStiffnessMatrix)
      .def_function(assembleInternalForces)
      .def_function(assembleMass)
      .def_function(assembleMassLumped)
      .def_function(getStableTimeStep)
      .def_function_nocopy(getExternalForce)
      .def_function_nocopy(getDisplacement)
      .def_function_nocopy(getPreviousDisplacement)
      .def_function_nocopy(getCurrentPosition)
      .def_function_nocopy(getIncrement)
      .def_function_nocopy(getInternalForce)
      .def_function_nocopy(getMass)
      .def_function_nocopy(getVelocity)
      .def_function_nocopy(getAcceleration)
      .def_function_nocopy(getInternalForce)
      .def_function_nocopy(getBlockedDOFs)
      .def_function_nocopy(getMesh)
      .def(
          "getMaterial",
          [](SolidMechanicsModel & self, UInt material_id) -> decltype(auto) {
            return self.getMaterial(material_id);
          },
          py::arg("material_id"), py::return_value_policy::reference)
      .def(
          "getMaterial",
          [](SolidMechanicsModel & self, const ID & material_name)
              -> decltype(auto) { return self.getMaterial(material_name); },
          py::arg("material_name"), py::return_value_policy::reference)
      .def(
          "getMaterial",
          [](SolidMechanicsModel & self, const Element & element)
              -> decltype(auto) { return self.getMaterial(element); },
          py::arg("element"), py::return_value_policy::reference)

      .def("getNbMaterials", &SolidMechanicsModel::getNbMaterials)
      .def("getMaterialIndex", &SolidMechanicsModel::getMaterialIndex)
      .def("setMaterialSelector",
           [](SolidMechanicsModel & self,
              std::shared_ptr<ConstitutiveLawSelector> material_selector) {
             self.setMaterialSelector(material_selector);
           })
      .def("getMaterialSelector", &SolidMechanicsModel::getMaterialSelector)
      .def(
          "getMaterialByElement",
          [](const SolidMechanicsModel & self) -> decltype(auto) {
            return self.getMaterialByElement();
          },
          py::return_value_policy::reference, py::keep_alive<0, 1>())
      .def("reassignMaterial", &SolidMechanicsModel::reassignMaterial)
      .def(
          "registerNewMaterial",
          [](SolidMechanicsModel & self, const ID & mat_name,
             const ID & mat_type, const ID & opt_param) -> decltype(auto) {
            return self.registerNewMaterial(mat_name, mat_type, opt_param);
          },
          py::arg("material_name"), py::arg("material_type"),
          py::arg("option") = "", py::return_value_policy::reference);
  //      .def("initMaterials", &SolidMechanicsModel::initMaterials)
}

} // namespace akantu
