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
#include "py_solver.hh"
#include "py_aka_array.hh"
/* -------------------------------------------------------------------------- */
#include <model.hh>
#include <non_linear_solver.hh>
#include <sparse_matrix_aij.hh>
#include <sparse_solver.hh>
#include <terms_to_assemble.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
namespace py = pybind11;
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
void register_solvers(py::module & mod) {
  py::class_<SparseMatrix>(mod, "SparseMatrix")
      .def("getMatrixType", &SparseMatrix::getMatrixType)
      .def("size", &SparseMatrix::size)
      .def("zero", &SparseMatrix::zero)
      .def("saveProfile", &SparseMatrix::saveProfile)
      .def("saveMatrix", &SparseMatrix::saveMatrix)
      .def(
          "add", [](SparseMatrix & self, Int i, Int j) { self.add(i, j); },
          "Add entry in the profile")
      .def(
          "add",
          [](SparseMatrix & self, Int i, Int j, Real value) {
            self.add(i, j, value);
          },
          "Add the value to the matrix")
      .def(
          "add",
          [](SparseMatrix & self, SparseMatrix & A, Real alpha) {
            self.add(A, alpha);
          },
          "Add a matrix to the matrix", py::arg("A"), py::arg("alpha") = 1.)

      .def("isFinite", &SparseMatrix::isFinite)

      .def("getRelease",
           [](const SparseMatrix & self) -> Int { return self.getRelease(); })
      .def("__call__",
           [](const SparseMatrix & self, Int i, Int j) { return self(i, j); })
      .def("getRelease", &SparseMatrix::getRelease);

  py::class_<SparseMatrixAIJ, SparseMatrix>(mod, "SparseMatrixAIJ")
      .def("getIRN", &SparseMatrixAIJ::getIRN)
      .def("getJCN", &SparseMatrixAIJ::getJCN)
      .def("getA", &SparseMatrixAIJ::getA);

  py::class_<SolverVector>(mod, "SolverVector")
      .def(
          "getValues",
          [](SolverVector & self) -> decltype(auto) {
            return static_cast<const Array<Real> &>(self);
          },
          py::return_value_policy::reference_internal,
          "Transform this into a vector, Is not copied.")
      .def("isDistributed",
           [](const SolverVector & self) { return self.isDistributed(); });

  py::class_<SparseSolver, Parsable>(mod, "Solver");

  py::class_<TermsToAssemble::TermToAssemble>(mod, "TermToAssemble")
      .def(py::init<Int, Int>())
      .def(py::self += Real())
      .def_property_readonly("i", &TermsToAssemble::TermToAssemble::i)
      .def_property_readonly("j", &TermsToAssemble::TermToAssemble::j);

  py::class_<TermsToAssemble>(mod, "TermsToAssemble")
      .def(py::init<const ID &, const ID &>())
      .def("getDOFIdM", &TermsToAssemble::getDOFIdM)
      .def("getDOFIdN", &TermsToAssemble::getDOFIdN)
      .def(
          "__call__",
          [](TermsToAssemble & self, Int i, Int j, Real val) {
            auto & term = self(i, j);
            term = val;
            return term;
          },
          py::arg("i"), py::arg("j"), py::arg("val") = 0.,
          py::return_value_policy::reference);
}

} // namespace akantu
