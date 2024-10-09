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
#include <dumper_iohelper_paraview.hh>
#include <element_group.hh>
#include <mesh.hh>
#include <node_group.hh>
/* -------------------------------------------------------------------------- */
#include <dumpable_inline_impl.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
namespace py = pybind11;
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
void register_group_manager(py::module & mod) {
  /* ------------------------------------------------------------------------ */
  py::class_<NodeGroup>(mod, "NodeGroup")
      .def(
          "getNodes",
          [](NodeGroup & self) -> decltype(auto) { return self.getNodes(); },
          py::return_value_policy::reference)
      .def("__len__", &NodeGroup::size)
      .def(
          "__iter__",
          [](const NodeGroup & self) {
            return py::make_iterator(self.begin(), self.end());
          },
          py::keep_alive<0, 1>())
      .def("__contains__", [](const NodeGroup & self,
                              UInt node) { return self.find(node) != -1; })
      .def("getName", &NodeGroup::getName)
      .def("clear", &NodeGroup::clear)
      .def("empty", &NodeGroup::empty)
      .def("append", &NodeGroup::append)
      .def(
          "add",
          [](NodeGroup & self, UInt node, bool check_for_duplicate) {
            auto && it = self.add(node, check_for_duplicate);
            return *it;
          },
          py::arg("node"), py::arg("check_for_duplicate") = true)
      .def("remove", &NodeGroup::add);

  /* ------------------------------------------------------------------------ */
  py::class_<ElementGroup, Dumpable>(mod, "ElementGroup")
      .def(
          "getNodeGroup",
          [](ElementGroup & self) -> decltype(auto) {
            return self.getNodeGroup();
          },
          py::return_value_policy::reference)
      .def("getName", &ElementGroup::getName)
      .def(
          "getElements",
          [](ElementGroup & self) -> decltype(auto) {
            return self.getElements();
          },
          py::return_value_policy::reference)
      .def(
          "getNodeGroup",
          [](ElementGroup & self) -> decltype(auto) {
            return self.getNodeGroup();
          },
          py::return_value_policy::reference)
      .def("__len__", [](const ElementGroup & self) { return self.size(); })
      .def("clear", [](ElementGroup & self) { self.clear(); })
      .def("empty", &ElementGroup::empty)
      .def("append", &ElementGroup::append)
      .def("optimize", &ElementGroup::optimize)
      .def(
          "add",
          [](ElementGroup & self, const Element & element, bool add_nodes,
             bool check_for_duplicate) {
            self.add(element, add_nodes, check_for_duplicate);
          },
          py::arg("element"), py::arg("add_nodes") = false,
          py::arg("check_for_duplicate") = true)
      .def("fillFromNodeGroup", &ElementGroup::fillFromNodeGroup)
      .def("addDimension", &ElementGroup::addDimension);

  /* ------------------------------------------------------------------------ */
  py::class_<GroupManager>(mod, "GroupManager")
      .def(
          "getElementGroup",
          [](GroupManager & self, const std::string & name) -> decltype(auto) {
            return self.getElementGroup(name);
          },
          py::return_value_policy::reference)
      .def("iterateElementGroups",
           [](GroupManager & self) -> decltype(auto) {
             std::vector<std::reference_wrapper<ElementGroup>> groups;
             for (auto & group : self.iterateElementGroups()) {
               groups.emplace_back(group);
             }
             return groups;
           })
      .def("iterateNodeGroups",
           [](GroupManager & self) -> decltype(auto) {
             std::vector<std::reference_wrapper<NodeGroup>> groups;
             for (auto & group : self.iterateNodeGroups()) {
               groups.emplace_back(group);
             }
             return groups;
           })
      .def("createNodeGroup", &GroupManager::createNodeGroup,
           py::return_value_policy::reference)
      .def(
          "createElementGroup",
          [](GroupManager & self, const std::string & id, Int spatial_dimension,
             bool replace_group) -> decltype(auto) {
            return self.createElementGroup(id, spatial_dimension,
                                           replace_group);
          },
          py::arg("id"), py::arg("spatial_dimension"),
          py::arg("replace_group") = false, py::return_value_policy::reference)
      .def("createGroupsFromMeshDataUInt",
           &GroupManager::createGroupsFromMeshData<UInt>)
      .def("createElementGroupFromNodeGroup",
           &GroupManager::createElementGroupFromNodeGroup, py::arg("name"),
           py::arg("node_group"), py::arg("dimension") = _all_dimensions)
      .def(
          "getNodeGroup",
          [](GroupManager & self, const std::string & name) -> decltype(auto) {
            return self.getNodeGroup(name);
          },
          py::return_value_policy::reference)
      .def(
          "nodeGroups",
          [](GroupManager & self) {
            std::vector<NodeGroup *> groups;
            for (auto & g : self.iterateNodeGroups()) {
              groups.push_back(&g);
            }
            return groups;
          },
          py::return_value_policy::reference)
      .def(
          "elementGroups",
          [](GroupManager & self) {
            std::vector<ElementGroup *> groups;
            for (auto & g : self.iterateElementGroups()) {
              groups.push_back(&g);
            }
            return groups;
          },
          py::return_value_policy::reference)
      .def("createBoundaryGroupFromGeometry",
           &GroupManager::createBoundaryGroupFromGeometry)
      .def("createClusters", [](GroupManager & self, Int element_dimension,
                                std::string cluster_name_prefix,
                                const GroupManager::ClusteringFilter & filter) {
        return self.createClusters(element_dimension, cluster_name_prefix,
                                   filter);
      });

  /* ------------------------------------------------------------------------ */
  py::class_<GroupManager::ClusteringFilter>(mod, "ClusteringFilter")
      .def(py::init<>());
}

} // namespace akantu
