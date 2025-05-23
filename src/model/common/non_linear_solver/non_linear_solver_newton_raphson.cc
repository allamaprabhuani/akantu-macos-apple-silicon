/**
 * Copyright (©) 2015-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "non_linear_solver_newton_raphson.hh"
#include "communicator.hh"
#include "dof_manager_default.hh"
#include "solver_callback.hh"
#include "solver_vector.hh"
#include "sparse_solver_eigen.hh"
#if defined(AKANTU_USE_MUMPS)
#include "sparse_solver_mumps.hh"
#endif
#if defined(AKANTU_USE_PETSC)
#include "dof_manager_petsc.hh"
#include "sparse_solver_petsc.hh"
#endif
/* -------------------------------------------------------------------------- */

// #if !defined(AKANTU_USE_MUMPS) && !defined(AKANTU_USE_PETSC)
// #include "sparse_solver_eigen.hh"
// namespace akantu {
// using SparseSolverType = SparseSolverEigen;
// }
// #else
// #include "sparse_solver_mumps.hh"
// namespace akantu {
// using SparseSolverType = SparseSolverMumps;
// }
// #endif

namespace akantu {

/* -------------------------------------------------------------------------- */

// NonLinearSolverLinear::NonLinearSolverLinear(
//     DOFManagerDefault & dof_manager,
//     const NonLinearSolverType & non_linear_solver_type, const ID & id)
//     : NonLinearSolver(dof_manager, non_linear_solver_type, id),
//       solver(dof_manager, "J", id + ":sparse_solver") {

//   this->supported_type.insert(NonLinearSolverType::_linear);
//   this->checkIfTypeIsSupported();
// }

/* -------------------------------------------------------------------------- */

NonLinearSolverNewtonRaphson::NonLinearSolverNewtonRaphson(
    DOFManager & dof_manager, const ModelSolverOptions & solver_options,
    const ID & id)
    : NonLinearSolver(dof_manager, solver_options, id) {

  switch (solver_options.sparse_solver_type) {
#if defined(AKANTU_USE_MUMPS)
  case SparseSolverType::_mumps:
    if (aka::is_of_type<DOFManagerDefault>(dof_manager)) {
      sparse_solver = std::make_unique<SparseSolverMumps>(
          aka::as_type<DOFManagerDefault>(dof_manager), "J",
          id + ":sparse_solver");
    } else {
      AKANTU_EXCEPTION(
          "Cannot use MUMPS sparse solver without a DOFManagerDefault");
    }
    break;
#endif
  case SparseSolverType::_eigen:
    if (aka::is_of_type<DOFManagerDefault>(dof_manager)) {
      sparse_solver = std::make_unique<SparseSolverEigen>(
          aka::as_type<DOFManagerDefault>(dof_manager), "J",
          id + ":sparse_solver");
    } else {
      AKANTU_EXCEPTION(
          "Cannot use EIGEN sparse solver without a DOFManagerDefault");
    }
    break;
#if defined(AKANTU_USE_PETSC)
  case SparseSolverType::_petsc:
    if (aka::is_of_type<DOFManagerPETSc>(dof_manager)) {
      sparse_solver = std::make_unique<SparseSolverPETSc>(
          aka::as_type<DOFManagerPETSc>(dof_manager), "J",
          id + ":sparse_solver");
    } else {
      AKANTU_EXCEPTION(
          "Cannot use PETSC sparse solver without a DOFManagerPETSc");
    }
    break;
#endif
  case SparseSolverType::_auto:
    AKANTU_TO_IMPLEMENT();
    break;
  default:
    AKANTU_EXCEPTION(solver_options.sparse_solver_type
                     << " compilation was not activated");
  }
  this->supported_type.insert(NonLinearSolverType::_newton_raphson_modified);
  this->supported_type.insert(NonLinearSolverType::_newton_raphson_contact);
  this->supported_type.insert(NonLinearSolverType::_newton_raphson);
  this->supported_type.insert(NonLinearSolverType::_linear);

  this->checkIfTypeIsSupported();

  this->registerParam("threshold", convergence_criteria, 1e-10, _pat_parsmod,
                      "Threshold to consider results as converged");
  this->registerParam("convergence_type", convergence_criteria_type,
                      SolveConvergenceCriteria::_solution, _pat_parsmod,
                      "Type of convergence criteria");
  this->registerParam("max_iterations", max_iterations, 10, _pat_parsmod,
                      "Max number of iterations");
  this->registerParam("error", error, _pat_readable, "Last reached error");
  this->registerParam("nb_iterations", n_iter, _pat_readable,
                      "Last reached number of iterations");
  this->registerParam("converged", converged, _pat_readable,
                      "Did last solve converged");
  this->registerParam("force_linear_recompute", force_linear_recompute, true,
                      _pat_modifiable,
                      "Force reassembly of the jacobian matrix");
  this->registerParam("verbose_iteration", verbose_iteration, false,
                      _pat_modifiable, "Make internal iterations verbose");
}

/* ------------------------------------------------------------------------ */
NonLinearSolverNewtonRaphson::~NonLinearSolverNewtonRaphson() = default;

/* ------------------------------------------------------------------------ */
void NonLinearSolverNewtonRaphson::solve(SolverCallback & solver_callback) {
  if (this->linear) {
    NonLinearSolverNewtonRaphson::solve_linear(solver_callback);
    return;
  }
  solver_callback.beforeSolveStep();
  this->dof_manager.updateGlobalBlockedDofs();

  solver_callback.predictor();

  if (non_linear_solver_type == NonLinearSolverType::_linear and
      solver_callback.canSplitResidual()) {
    solver_callback.assembleMatrix("K");
  }

  this->assembleResidual(solver_callback);

  if (this->non_linear_solver_type ==
          NonLinearSolverType::_newton_raphson_modified ||
      (this->non_linear_solver_type == NonLinearSolverType::_linear &&
       this->force_linear_recompute)) {
    solver_callback.assembleMatrix("J");
    this->force_linear_recompute = false;
  }

  this->n_iter = 0;
  this->converged = false;

  this->convergence_criteria_normalized = this->convergence_criteria;

  if (this->convergence_criteria_type == SolveConvergenceCriteria::_residual) {
    this->converged = this->testConvergence(this->dof_manager.getResidual());

    if (this->converged) {
      return;
    }

    this->convergence_criteria_normalized =
        this->error * this->convergence_criteria;
  }

  do {
    if (this->non_linear_solver_type == NonLinearSolverType::_newton_raphson or
        this->non_linear_solver_type ==
            NonLinearSolverType::_newton_raphson_contact) {
      solver_callback.assembleMatrix("J");
    }

    // static bool first{true};
    // if (first) {
    //   this->dof_manager.getMatrix("J").saveMatrix("J_default.mtx");
    // }
    // first = false;
    this->sparse_solver->solve();
    solver_callback.corrector();

    // EventManager::sendEvent(NonLinearSolver::AfterSparseSolve(method));

    if (this->convergence_criteria_type ==
        SolveConvergenceCriteria::_residual) {
      this->assembleResidual(solver_callback);
      this->converged = this->testConvergence(this->dof_manager.getResidual());
    } else {
      this->converged = this->testConvergence(this->dof_manager.getSolution());
    }

    if (this->convergence_criteria_type ==
            SolveConvergenceCriteria::_solution and
        not this->converged) {
      this->assembleResidual(solver_callback);
    }

    this->n_iter++;
    AKANTU_DEBUG_INFO(
        "[" << this->convergence_criteria_type << "] Convergence iteration "
            << std::setw(std::log10(this->max_iterations)) << this->n_iter
            << ": error " << this->error << (this->converged ? " < " : " > ")
            << this->convergence_criteria_normalized);

    if (verbose_iteration and dof_manager.getCommunicator().whoAmI() == 0) {
      std::cout << "[" << this->convergence_criteria_type
                << "] Convergence iteration "
                << std::setw(std::log10(this->max_iterations)) << this->n_iter
                << ": error " << std::scientific << this->error
                << (this->converged ? " < " : " > ")
                << this->convergence_criteria_normalized << " ["
                << this->convergence_criteria << "/"
                << (this->convergence_criteria_normalized /
                    this->convergence_criteria)
                << "]\n";
    }

  } while (not this->converged and this->n_iter < this->max_iterations);

  // this makes sure that you have correct strains and stresses after the
  // solveStep function (e.g., for dumping)
  if (this->convergence_criteria_type == SolveConvergenceCriteria::_solution) {
    this->assembleResidual(solver_callback);
  }

  this->converged =
      this->converged and not(this->n_iter > this->max_iterations);

  solver_callback.afterSolveStep(this->converged);

  if (not this->converged) {
    AKANTU_CUSTOM_EXCEPTION(debug::NLSNewtonRaphsonNotConvergedException(
        this->convergence_criteria_normalized, this->n_iter, this->error,
        this->convergence_criteria_type));

    AKANTU_DEBUG_WARNING("[" << this->convergence_criteria_type
                             << "] Convergence not reached after "
                             << std::setw(std::log10(this->max_iterations))
                             << this->n_iter << " iteration"
                             << (this->n_iter == 1 ? "" : "s") << "!");
  }
}

/* ------------------------------------------------------------------------ */
bool NonLinearSolverNewtonRaphson::testConvergence(
    const SolverVector & solver_vector) {
  AKANTU_DEBUG_IN();

  const auto & blocked_dofs = this->dof_manager.getGlobalBlockedDOFs();

  const Array<Real> & array(solver_vector);
  Int nb_degree_of_freedoms = array.size();

  auto arr_it = array.begin();
  auto bld_it = blocked_dofs.begin();

  Real norm = 0.;
  for (Int n = 0; n < nb_degree_of_freedoms; ++n, ++arr_it, ++bld_it) {
    bool is_local_node = this->dof_manager.isLocalOrMasterDOF(n);
    if ((!*bld_it) && is_local_node) {
      norm += *arr_it * *arr_it;
    }
  }

  dof_manager.getCommunicator().allReduce(norm, SynchronizerOperation::_sum);

  norm = std::sqrt(norm);

  AKANTU_DEBUG_ASSERT(!Math::isnan(norm),
                      "Something went wrong in the solve phase");

  this->error = norm;

  return (error < this->convergence_criteria_normalized);
}

/* ------------------------------------------------------------------------ */

void NonLinearSolverNewtonRaphson::solve_linear(
    SolverCallback & solver_callback) {
  solver_callback.beforeSolveStep();
  this->dof_manager.updateGlobalBlockedDofs();

  solver_callback.predictor();

  solver_callback.assembleMatrix("J");

  // Residual computed after J to allow the model to use K to compute the
  // residual
  this->assembleResidual(solver_callback);

  this->sparse_solver->solve();

  solver_callback.corrector();

  if (solver_callback.canSplitResidual()) {
    solver_callback.assembleResidual("internal");
  } else {
    this->assembleResidual(solver_callback);
  }

  solver_callback.afterSolveStep(true);
}

/* ------------------------------------------------------------------------- */

} // namespace akantu
