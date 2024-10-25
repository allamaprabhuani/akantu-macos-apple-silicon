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
#include <contact_detector.hh>
#include <contact_element.hh>
#include <contact_mechanics_model.hh>
#include <geometry_utils.hh>
#include <mesh_events.hh>
#include <parsable.hh>
#include <surface_selector.hh>
/* -------------------------------------------------------------------------- */
#include <algorithm>
/* -------------------------------------------------------------------------- */
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
namespace py = pybind11;
/* -------------------------------------------------------------------------- */

namespace akantu {
/* -------------------------------------------------------------------------- */
#define def_function_nocopy(func_name)                                         \
  def(                                                                         \
      #func_name,                                                              \
      [](ContactMechanicsModel & self) -> decltype(auto) {                     \
        return self.func_name();                                               \
      },                                                                       \
      py::return_value_policy::reference)

#define def_function(func_name)                                                \
  def(#func_name, [](ContactMechanicsModel & self) -> decltype(auto) {         \
    return self.func_name();                                                   \
  })

namespace {
  class ContactElementsView {
  public:
    ContactElementsView(const Array<ContactElement> & contact_elements)
        : contact_elements(contact_elements) {}

    auto begin() const { return contact_elements.begin(); }
    auto end() const { return contact_elements.end(); }

    auto size() const { return contact_elements.size(); }
    auto operator[](size_t i) const { return contact_elements(i); }

    auto contains(const ContactElement & contact_element) const {
      return std::find(contact_elements.begin(), contact_elements.end(),
                       contact_element) != contact_elements.end();
    }

  private:
    const Array<ContactElement> & contact_elements;
  };
} // namespace

/* -------------------------------------------------------------------------- */
void register_contact_mechanics_model(py::module & mod) {
  py::class_<ContactDetector>(mod, "ContactDetector",
                              py::multiple_inheritance())
      .def(py::init<Mesh &, const ID &>(), py::arg("mesh"),
           py::arg("id") = "contact_detector")
      .def(py::init<Mesh &, Array<Real>, const ID &>(), py::arg("mesh"),
           py::arg("positions"), py::arg("id") = "contact_detector")
      .def("setSurfaceSelector", &ContactDetector::setSurfaceSelector);

  py::class_<SurfaceSelector, std::shared_ptr<SurfaceSelector>>(
      mod, "SurfaceSelector", py::multiple_inheritance())
      .def(py::init<Mesh &>(), py::arg("mesh"));

  py::class_<PhysicalSurfaceSelector, SurfaceSelector,
             std::shared_ptr<PhysicalSurfaceSelector>>(
      mod, "PhysicalSurfaceSelector")
      .def(py::init<Mesh &>(), py::arg("mesh"));

  py::class_<CohesiveSurfaceSelector, SurfaceSelector,
             std::shared_ptr<CohesiveSurfaceSelector>>(
      mod, "CohesiveSurfaceSelector")
      .def(py::init<Mesh &>(), py::arg("mesh"));

  py::class_<AllSurfaceSelector, SurfaceSelector,
             std::shared_ptr<AllSurfaceSelector>>(mod, "AllSurfaceSelector")
      .def(py::init<Mesh &>(), py::arg("mesh"));

  mod.attr("ContactMechanicsModelOptions") = mod.attr("ModelOptions");

  /* ------------------------------------------------------------------------ */
  py::class_<ContactElementsView>(mod, "ContactElementsView")
      .def("__iter__",
           [](const ContactElementsView & self) {
             return py::make_iterator(self.begin(), self.end());
           })
      .def("__size__",
           [](const ContactElementsView & self) { return self.size(); })
      .def(
          "__contains__",
          [](const ContactElementsView & self, const ContactElement & element) {
            return self.contains(element);
          })
      .def("__getitem__",
           [](const ContactElementsView & self, size_t i) { return self[i]; });

  /* ------------------------------------------------------------------------ */
  py::class_<ContactMechanicsModel, Model>(mod, "ContactMechanicsModel",
                                           py::multiple_inheritance())
      .def(py::init<Mesh &, UInt, const ID &, std::shared_ptr<DOFManager>,
                    const ModelType>(),
           py::arg("mesh"), py::arg("spatial_dimension") = _all_dimensions,
           py::arg("id") = "contact_mechanics_model",
           py::arg("dof_manager") = nullptr,
           py::arg("model_type") = ModelType::_contact_mechanics_model)
      .def(
          "initFull",
          [](ContactMechanicsModel & self,
             const ContactMechanicsModelOptions & options) {
            self.initFull(options);
          },
          py::arg("options") = ContactMechanicsModelOptions())
      .def(
          "initFull",
          [](ContactMechanicsModel & self,
             const AnalysisMethod & analysis_method) {
            self.initFull(_analysis_method = analysis_method);
          },
          py::arg("_analysis_method"))
      .def_function(search)
      .def_function(assembleStiffnessMatrix)
      .def_function(assembleInternalForces)
      .def_function_nocopy(getExternalForce)
      .def_function_nocopy(getNormalForce)
      .def_function_nocopy(getTangentialForce)
      .def_function_nocopy(getInternalForce)
      .def_function_nocopy(getGaps)
      .def_function_nocopy(getNormals)
      .def_function_nocopy(getNodalArea)
      .def("getContactDetector", &ContactMechanicsModel::getContactDetector,
           py::return_value_policy::reference)
      .def("getContactElements", [](ContactMechanicsModel & self) {
        return ContactElementsView(self.getContactElements());
      });

  py::class_<ContactElement>(mod, "ContactElement")
      .def(py::init<>())
      .def_readwrite("master", &ContactElement::master)
      .def_readwrite("slave", &ContactElement::slave)
      .def("__repr__", [](ContactElement & self) {
        return "{master: " + std::to_string(self.master) +
               ", slave: " + std::to_string(self.slave) + "}";
      });

  py::class_<GeometryUtils>(mod, "GeometryUtils")
      .def_static(
          "normal",
          [](const Mesh & mesh, const Array<Real> & positions,
             const Element & element, bool outward) {
            auto && coords =
                mesh.extractNodalValuesFromElement(positions, element);
            return GeometryUtils::normal(mesh, coords, element, outward);
          },
          py::arg("mesh"), py::arg("positions"), py::arg("element"),
          py::arg("outward") = true)
      .def_static(
          "covariantBasis",
          [](const Mesh & mesh, const Array<Real> & positions,
             const Element & element, const Vector<Real> & normal,
             Vector<Real> & natural_coord) {
            auto && coords =
                mesh.extractNodalValuesFromElement(positions, element);
            return GeometryUtils::covariantBasis(coords, element, normal,
                                                 natural_coord);
          },
          py::arg("mesh"), py::arg("positions"), py::arg("element"),
          py::arg("normal"), py::arg("natural_projection"))
      .def_static(
          "curvature",
          [](const Mesh & mesh, const Array<Real> & positions,
             const Element & element, const Vector<Real> & natural_coord) {
            auto && coords =
                mesh.extractNodalValuesFromElement(positions, element);
            return GeometryUtils::curvature(coords, element, natural_coord);
          })
      .def_static(
          "contravariantBasis",
          [](const Vector<Real> & covariant) {
            return GeometryUtils::contravariantBasis(covariant);
          },
          py::arg("covariant_basis"))
      .def_static(
          "realProjection",
          [](const Mesh & mesh, const Array<Real> & positions,
             const Vector<Real> & slave, const Element & element,
             const Vector<Real> & normal) {
            auto && coords =
                mesh.extractNodalValuesFromElement(positions, element);
            return GeometryUtils::realProjection(coords, slave, normal);
          },
          py::arg("mesh"), py::arg("positions"), py::arg("slave"),
          py::arg("element"), py::arg("normal"))
      .def_static("isBoundaryElement", &GeometryUtils::isBoundaryElement);
}

} // namespace akantu
