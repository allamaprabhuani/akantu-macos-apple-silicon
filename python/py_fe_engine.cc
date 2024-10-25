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
#include "py_aka_common.hh"
/* -------------------------------------------------------------------------- */
#include <element.hh>
#include <fe_engine.hh>
#include <integration_point.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
#include <memory>
/* -------------------------------------------------------------------------- */
namespace py = pybind11;
/* -------------------------------------------------------------------------- */

namespace akantu {

void register_fe_engine(py::module & mod) {
  py::class_<Element>(mod, "Element")
      .def(py::init([](ElementType type, Int id, GhostType ghost_type) {
             return std::make_unique<Element>(Element{type, id, ghost_type});
           }),
           py::arg("type"), py::arg("ghost_type"),
           py::arg("ghost_type") = _not_ghost)
      .def("__lt__",
           [](Element & self, const Element & other) { return (self < other); })
      .def("__repr__", [](Element & self) { return std::to_string(self); });

  mod.attr("ElementNull") = ElementNull;

  py::class_<FEEngine>(mod, "FEEngine")
      .def(
          "getNbIntegrationPoints",
          [](FEEngine & fem, ElementType type, GhostType ghost_type) {
            return fem.getNbIntegrationPoints(type, ghost_type);
          },
          py::arg("type"), py::arg("ghost_type") = _not_ghost)
      .def("initShapeFunctions", &FEEngine::initShapeFunctions,
           py::arg("ghost_type") = _not_ghost)
      .def(
          "integrate",
          [](FEEngine & self,
             const std::shared_ptr<ElementTypeMapArray<Real>> f,
             std::shared_ptr<ElementTypeMapArray<Real>> intf,
             std::shared_ptr<const ElementTypeMapArray<Idx>> filter_elements) {
            self.integrate(*f, *intf, filter_elements.get());
          },
          py::arg("f"), py::arg("intf"), py::arg("filter_elements") = nullptr)

      .def(
          "integrate",
          [](FEEngine & self, const Array<Real> & f, Array<Real> & intf,
             Int nb_degree_of_freedom, ElementType type,
             GhostType ghost_type = _not_ghost,
             const Array<Idx> & filter_elements = empty_filter) {
            // Check to overcome the empty_filter that is maped as a numpy array
            // without an id
            if (filter_elements.data() == empty_filter.data()) {
              self.integrate(f, intf, nb_degree_of_freedom, type, ghost_type,
                             empty_filter);
            } else {
              self.integrate(f, intf, nb_degree_of_freedom, type, ghost_type,
                             filter_elements);
            }
          },
          py::arg("f"), py::arg("intf"), py::arg("nb_degree_of_freedom"),
          py::arg("type"), py::arg("ghost_type") = _not_ghost,
          py::arg("filter_elements") = empty_filter)
      .def(
          "integrate",
          [](FEEngine & self, const Array<Real> & f, ElementType type,
             GhostType ghost_type = _not_ghost,
             const Array<Idx> & filter_elements = empty_filter) -> Real {
            // Check to overcome the empty_filter that is maped as a numpy array
            // without an id
            if (filter_elements.data() == empty_filter.data()) {
              return self.integrate(f, type, ghost_type, empty_filter);
            } else {
              return self.integrate(f, type, ghost_type, filter_elements);
            }
          },
          py::arg("f"), py::arg("type"), py::arg("ghost_type") = _not_ghost,
          py::arg("filter_elements") = empty_filter)
      .def(
          "gradientOnIntegrationPoints",
          [](FEEngine & fem, const Array<Real> & u, Array<Real> & nablauq,
             const Int nb_degree_of_freedom, ElementType type,
             GhostType ghost_type, const Array<Idx> & filter_elements) {
            // Check to overcome the empty_filter that is maped as a numpy array
            // without an id
            if (filter_elements.data() == empty_filter.data()) {
              fem.gradientOnIntegrationPoints(u, nablauq, nb_degree_of_freedom,
                                              type, ghost_type, empty_filter);
            } else {
              fem.gradientOnIntegrationPoints(u, nablauq, nb_degree_of_freedom,
                                              type, ghost_type,
                                              filter_elements);
            }
          },
          py::arg("u"), py::arg("nablauq"), py::arg("nb_degree_of_freedom"),
          py::arg("type"), py::arg("ghost_type") = _not_ghost,
          py::arg("filter_elements") = empty_filter)
      .def(
          "interpolateOnIntegrationPoints",
          [](FEEngine & self, const Array<Real> & u, Array<Real> & uq,
             Int nb_degree_of_freedom, ElementType type, GhostType ghost_type,
             const Array<Idx> & filter_elements) {
            // Check to overcome the empty_filter that is maped as a numpy array
            // without an id
            if (filter_elements.data() == empty_filter.data()) {
              self.interpolateOnIntegrationPoints(
                  u, uq, nb_degree_of_freedom, type, ghost_type, empty_filter);
            } else {
              self.interpolateOnIntegrationPoints(u, uq, nb_degree_of_freedom,
                                                  type, ghost_type,
                                                  filter_elements);
            }
          },
          py::arg("u"), py::arg("uq"), py::arg("nb_degree_of_freedom"),
          py::arg("type"), py::arg("ghost_type") = _not_ghost,
          py::arg("filter_elements") = empty_filter)
      .def(
          "interpolateOnIntegrationPoints",
          [](FEEngine & self, const Array<Real> & u,
             std::shared_ptr<ElementTypeMapArray<Real>> uq,
             std::shared_ptr<const ElementTypeMapArray<Idx>> filter_elements) {
            self.interpolateOnIntegrationPoints(u, *uq, filter_elements.get());
          },
          py::arg("u"), py::arg("uq"), py::arg("filter_elements") = nullptr)

      .def(
          "computeBtD",
          [](FEEngine & self, const Array<Real> & Ds, Array<Real> & BtDs,
             ElementType type, GhostType ghost_type = _not_ghost,
             const Array<Idx> & filter_elements = empty_filter) {
            // Check to overcome the empty_filter that is maped as a numpy array
            // without an id
            if (filter_elements.data() == empty_filter.data()) {
              self.computeBtD(Ds, BtDs, type, ghost_type, empty_filter);
            } else {
              self.computeBtD(Ds, BtDs, type, ghost_type, filter_elements);
            }
          },
          py::arg("Ds"), py::arg("BtDs"), py::arg("type"),
          py::arg("ghost_type") = _not_ghost,
          py::arg("filter_elements") = empty_filter)

      .def(
          "computeBtDB",
          [](FEEngine & self, const Array<Real> & Ds, Array<Real> & BtDBs,
             Int order_d, ElementType type, GhostType ghost_type = _not_ghost,
             const Array<Idx> & filter_elements = empty_filter) {
            // Check to overcome the empty_filter that is maped as a numpy array
            // without an id
            if (filter_elements.data() == empty_filter.data()) {
              self.computeBtDB(Ds, BtDBs, order_d, type, ghost_type,
                               empty_filter);
            } else {
              self.computeBtDB(Ds, BtDBs, order_d, type, ghost_type,
                               filter_elements);
            }
          },
          py::arg("Ds"), py::arg("BtDBs"), py::arg("order_d"), py::arg("type"),
          py::arg("ghost_type") = _not_ghost,
          py::arg("filter_elements") = empty_filter)

      .def(
          "computeNtb",
          [](FEEngine & self, const Array<Real> & bs, Array<Real> & Ntbs,
             ElementType type, GhostType ghost_type = _not_ghost,
             const Array<Idx> & filter_elements = empty_filter) {
            // Check to overcome the empty_filter that is maped as a numpy array
            // without an id
            if (filter_elements.data() == empty_filter.data()) {
              self.computeNtb(bs, Ntbs, type, ghost_type, empty_filter);
            } else {
              self.computeNtb(bs, Ntbs, type, ghost_type, filter_elements);
            }
          },
          py::arg("bs"), py::arg("Ntbs"), py::arg("type"),
          py::arg("ghost_type") = _not_ghost,
          py::arg("filter_elements") = empty_filter)
      .def(
          "computeNtbN",
          [](FEEngine & self, const Array<Real> & bs, Array<Real> & NtbNs,
             ElementType type, GhostType ghost_type = _not_ghost,
             const Array<Idx> & filter_elements = empty_filter) {
            // Check to overcome the empty_filter that is maped as a numpy array
            // without an id
            if (filter_elements.data() == empty_filter.data()) {
              self.computeNtbN(bs, NtbNs, type, ghost_type, empty_filter);
            } else {
              self.computeNtbN(bs, NtbNs, type, ghost_type, filter_elements);
            }
          },
          py::arg("bs"), py::arg("NtbNs"), py::arg("type"),
          py::arg("ghost_type") = _not_ghost,
          py::arg("filter_elements") = empty_filter)
      .def(
          "computeIntegrationPointsCoordinates",
          [](FEEngine & self,
             std::shared_ptr<ElementTypeMapArray<Real>> coordinates,
             std::shared_ptr<const ElementTypeMapArray<Idx>> filter_elements)
              -> decltype(auto) {
            return self.computeIntegrationPointsCoordinates(
                *coordinates, filter_elements.get());
          },
          py::arg("coordinates"), py::arg("filter_elements") = nullptr)
      .def(
          "assembleFieldLumped",
          [](FEEngine & fem,
             const std::function<void(Matrix<Real> &, const Element &)> &
                 field_funct,
             const ID & matrix_id, const ID & dof_id, DOFManager & dof_manager,
             ElementType type, GhostType ghost_type) {
            fem.assembleFieldLumped(field_funct, matrix_id, dof_id, dof_manager,
                                    type, ghost_type);
          },
          py::arg("field_funct"), py::arg("matrix_id"), py::arg("dof_id"),
          py::arg("dof_manager"), py::arg("type"),
          py::arg("ghost_type") = _not_ghost)
      .def(
          "assembleFieldMatrix",
          [](FEEngine & fem,
             const std::function<void(Matrix<Real> &, const Element &)> &
                 field_funct,
             const ID & matrix_id, const ID & dof_id, DOFManager & dof_manager,
             ElementType type, GhostType ghost_type = _not_ghost) {
            fem.assembleFieldMatrix(field_funct, matrix_id, dof_id, dof_manager,
                                    type, ghost_type);
          },
          py::arg("field_funct"), py::arg("matrix_id"), py::arg("dof_id"),
          py::arg("dof_manager"), py::arg("type"),
          py::arg("ghost_type") = _not_ghost)
      .def("getElementInradius",
           [](FEEngine & self, const Element & element) {
             return self.getElementInradius(element);
           })
      .def("getNormalsOnIntegrationPoints",
           &FEEngine::getNormalsOnIntegrationPoints, py::arg("type"),
           py::arg("ghost_type") = _not_ghost,
           py::return_value_policy::reference)
      .def("getShapes", &FEEngine::getShapes, py::arg("type"),
           py::arg("ghost_type") = _not_ghost, py::arg("id") = 0,
           py::return_value_policy::reference)
      .def("getShapesDerivatives", &FEEngine::getShapesDerivatives,
           py::arg("type"), py::arg("ghost_type") = _not_ghost,
           py::arg("id") = 0, py::return_value_policy::reference)
      .def(
          "getMesh", [](FEEngine & fem) -> Mesh & { return fem.getMesh(); },
          py::return_value_policy::reference);

  py::class_<IntegrationPoint>(mod, "IntegrationPoint");
}
} // namespace akantu
