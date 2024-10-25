/**
 * Copyright (©) 2021-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "py_dof_manager.hh"
#include "py_aka_array.hh"
#include "py_akantu_pybind11_compatibility.hh"
/* -------------------------------------------------------------------------- */
#include <dof_manager.hh>
#include <non_linear_solver.hh>
#include <non_linear_solver_lumped.hh>
#include <non_linear_solver_newton_raphson.hh>
#include <solver_callback.hh>
#include <sparse_solver.hh>
#include <time_step_solver.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
namespace py = pybind11;
/* -------------------------------------------------------------------------- */

namespace akantu {

namespace {
  class PySolverCallback : public SolverCallback {
  public:
    using SolverCallback::SolverCallback;

    /// get the type of matrix needed
    MatrixType getMatrixType(const ID & matrix_id) const override {
      // NOLINTNEXTLINE
      PYBIND11_OVERRIDE_PURE(MatrixType, SolverCallback, getMatrixType,
                             matrix_id);
    }

    /// callback to assemble a Matrix
    void assembleMatrix(const ID & matrix_id) override {
      // NOLINTNEXTLINE
      PYBIND11_OVERRIDE_PURE(void, SolverCallback, assembleMatrix, matrix_id);
    }

    /// callback to assemble a lumped Matrix
    void assembleLumpedMatrix(const ID & matrix_id) override {
      // NOLINTNEXTLINE
      PYBIND11_OVERRIDE_PURE(void, SolverCallback, assembleLumpedMatrix,
                             matrix_id);
    }

    /// callback to assemble the residual (rhs)
    void assembleResidual() override {
      // NOLINTNEXTLINE
      PYBIND11_OVERRIDE_PURE(void, SolverCallback, assembleResidual, );
    }

    /// callback for the predictor (in case of dynamic simulation)
    void predictor() override {
      // NOLINTNEXTLINE
      PYBIND11_OVERRIDE(void, SolverCallback, predictor, );
    }

    /// callback for the corrector (in case of dynamic simulation)
    void corrector() override {
      // NOLINTNEXTLINE
      PYBIND11_OVERRIDE(void, SolverCallback, corrector, );
    }

    void beforeSolveStep() override {
      // NOLINTNEXTLINE
      PYBIND11_OVERRIDE(void, SolverCallback, beforeSolveStep, );
    }

    void afterSolveStep(bool converged) override {
      // NOLINTNEXTLINE
      PYBIND11_OVERRIDE(void, SolverCallback, afterSolveStep, converged);
    }
  };

  class PyInterceptSolverCallback : public InterceptSolverCallback {
  public:
    using InterceptSolverCallback::InterceptSolverCallback;

    MatrixType getMatrixType(const ID & matrix_id) const override {
      // NOLINTNEXTLINE
      PYBIND11_OVERRIDE(MatrixType, InterceptSolverCallback, getMatrixType,
                        matrix_id);
    }

    void assembleMatrix(const ID & matrix_id) override {
      // NOLINTNEXTLINE
      PYBIND11_OVERRIDE(void, InterceptSolverCallback, assembleMatrix,
                        matrix_id);
    }

    /// callback to assemble a lumped Matrix
    void assembleLumpedMatrix(const ID & matrix_id) override {
      // NOLINTNEXTLINE
      PYBIND11_OVERRIDE(void, InterceptSolverCallback, assembleLumpedMatrix,
                        matrix_id);
    }

    void assembleResidual() override {
      // NOLINTNEXTLINE
      PYBIND11_OVERRIDE(void, InterceptSolverCallback, assembleResidual, );
    }

    void predictor() override {
      // NOLINTNEXTLINE
      PYBIND11_OVERRIDE(void, InterceptSolverCallback, predictor, );
    }

    void corrector() override {
      // NOLINTNEXTLINE
      PYBIND11_OVERRIDE(void, InterceptSolverCallback, corrector, );
    }

    void beforeSolveStep() override {
      // NOLINTNEXTLINE
      PYBIND11_OVERRIDE(void, InterceptSolverCallback, beforeSolveStep, );
    }

    void afterSolveStep(bool converged) override {
      // NOLINTNEXTLINE
      PYBIND11_OVERRIDE(void, InterceptSolverCallback, afterSolveStep,
                        converged);
    }
  };

} // namespace

/* -------------------------------------------------------------------------- */
void register_dof_manager(py::module & mod) {
  py::class_<DOFManager, std::shared_ptr<DOFManager>>(mod, "DOFManager")
      .def("getMatrix", &DOFManager::getMatrix,
           py::return_value_policy::reference)
      .def(
          "getNewMatrix",
          [](DOFManager & self, const std::string & name,
             const std::string & matrix_to_copy_id) -> decltype(auto) {
            return self.getNewMatrix(name, matrix_to_copy_id);
          },
          py::return_value_policy::reference)
      .def(
          "getResidual",
          [](DOFManager & self) -> decltype(auto) {
            return self.getResidual();
          },
          py::return_value_policy::reference)
      .def("getArrayPerDOFs", &DOFManager::getArrayPerDOFs)
      .def(
          "hasMatrix",
          [](DOFManager & self, const ID & name) -> bool {
            return self.hasMatrix(name);
          },
          py::arg("name"))
      .def("assembleToResidual", &DOFManager::assembleToResidual,
           py::arg("dof_id"), py::arg("array_to_assemble"),
           py::arg("scale_factor") = 1.)
      .def("assembleToLumpedMatrix", &DOFManager::assembleToLumpedMatrix,
           py::arg("dof_id"), py::arg("array_to_assemble"),
           py::arg("lumped_mtx"), py::arg("scale_factor") = 1.)
      .def("assemblePreassembledMatrix",
           &DOFManager::assemblePreassembledMatrix, py::arg("matrix_id"),
           py::arg("terms"))
      .def("zeroResidual", &DOFManager::zeroResidual);

  py::class_<NonLinearSolver, Parsable>(mod, "NonLinearSolver")
      .def(
          "set",
          [](NonLinearSolver & self, const std::string & id, const Real & val) {
            if (id == "max_iterations") {
              self.set(id, int(val));
            } else {
              self.set(id, val);
            }
          })
      .def("set",
           [](NonLinearSolver & self, const std::string & id,
              const SolveConvergenceCriteria & val) { self.set(id, val); });

  py::class_<NonLinearSolverNewtonRaphson, NonLinearSolver>(
      mod, "NonLinearSolverNewtonRaphson")
      .def(
          "getSparseSolver",
          [](NonLinearSolverNewtonRaphson & self) -> SparseSolver & {
            return self.getSparseSolver();
          },
          py::return_value_policy::reference);

  py::class_<NonLinearSolverLumped, NonLinearSolver>(mod,
                                                     "NonLinearSolverLumped");

  py::class_<TimeStepSolver>(mod, "TimeStepSolver")
      .def("getIntegrationScheme", &TimeStepSolver::getIntegrationScheme);

  py::class_<SolverCallback, PySolverCallback>(mod, "SolverCallback")
      .def(py::init_alias<DOFManager &>())
      .def("getMatrixType", &SolverCallback::getMatrixType)
      .def("assembleMatrix", &SolverCallback::assembleMatrix)
      .def("assembleLumpedMatrix", &SolverCallback::assembleLumpedMatrix)
      .def("assembleResidual",
           [](SolverCallback & self) { self.assembleResidual(); })
      .def("predictor", &SolverCallback::predictor)
      .def("corrector", &SolverCallback::corrector)
      .def("beforeSolveStep", &SolverCallback::beforeSolveStep)
      .def("afterSolveStep", &SolverCallback::afterSolveStep)
      .def_property_readonly("dof_manager", &SolverCallback::getSCDOFManager,
                             py::return_value_policy::reference);

  py::class_<InterceptSolverCallback, SolverCallback,
             PyInterceptSolverCallback>(mod, "InterceptSolverCallback")
      .def(py::init_alias<SolverCallback &>());
}

} // namespace akantu
