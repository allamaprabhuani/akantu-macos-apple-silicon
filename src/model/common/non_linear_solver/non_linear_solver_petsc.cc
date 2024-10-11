/**
 * Copyright (©) 2018-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "non_linear_solver_petsc.hh"
#include "dof_manager_petsc.hh"
#include "mpi_communicator_data.hh"
#include "solver_callback.hh"
#include "solver_vector_petsc.hh"
#include "sparse_matrix_petsc.hh"
/* -------------------------------------------------------------------------- */
#include <petscoptions.h>
/* -------------------------------------------------------------------------- */

namespace akantu {

NonLinearSolverPETSc::NonLinearSolverPETSc(
    DOFManagerPETSc & dof_manager,
    const NonLinearSolverType & non_linear_solver_type,
    const SparseSolverType & sparse_solver_type, const ID & id)
    : NonLinearSolver(dof_manager, non_linear_solver_type, id) {

  if (sparse_solver_type != SparseSolverType::_petsc)
    AKANTU_EXCEPTION(
        "petsc non linear solver works only with petsc sparse solver");

  std::unordered_map<NonLinearSolverType, SNESType>
      petsc_non_linear_solver_types{
          {NonLinearSolverType::_newton_raphson, SNESNEWTONLS},
          {NonLinearSolverType::_linear, SNESKSPONLY},
          {NonLinearSolverType::_gmres, SNESNGMRES},
          {NonLinearSolverType::_bfgs, SNESQN},
          {NonLinearSolverType::_cg, SNESNCG}};

  this->has_internal_set_param = true;

  for (const auto & pair : petsc_non_linear_solver_types) {
    supported_type.insert(pair.first);
  }

  this->checkIfTypeIsSupported();

  auto && mpi_comm = dof_manager.getMPIComm();

  SNESCreate(mpi_comm, &snes);

  auto it = petsc_non_linear_solver_types.find(non_linear_solver_type);
  if (it != petsc_non_linear_solver_types.end()) {
    SNESSetType(snes, it->second);
  }

  SNESSetFromOptions(snes);
}

/* -------------------------------------------------------------------------- */
NonLinearSolverPETSc::~NonLinearSolverPETSc() { SNESDestroy(&snes); }

/* -------------------------------------------------------------------------- */
class NonLinearSolverPETScCallback {
public:
  NonLinearSolverPETScCallback(DOFManagerPETSc & dof_manager, SNES snes,
                               SolverVectorPETSc & x)
      : dof_manager(dof_manager), snes(snes), x_prev(x, "previous solution") {}

  void corrector(Vec & x) {
    PetscInt iteration;
    SNESGetIterationNumber(snes, &iteration);

    if (prev_iteration == iteration) {
      return;
    }

    prev_iteration = iteration;

    auto & dx = dof_manager._getSolution();
    VecWAXPY(dx, -1., x_prev, x);
    dof_manager.splitSolutionPerDOFs();

    callback->corrector();

    VecCopy(x, x_prev);
  }

  void assembleResidual(Vec x) {
    corrector(x);
    callback->assembleResidual();
  }

  void assembleJacobian(Vec x) {
    corrector(x);
    callback->assembleMatrix("J");
  }

  void reset() { prev_iteration = -1; }

  void setInitialSolution(SolverVectorPETSc & x) { VecCopy(x, x_prev); }

  void setCallback(SolverCallback & callback) { this->callback = &callback; }

private:
  SolverCallback * callback;
  DOFManagerPETSc & dof_manager;
  SNES snes;

  // SolverVectorPETSc & x;
  SolverVectorPETSc x_prev;
  PetscInt prev_iteration{-1};
}; // namespace akantu

/* -------------------------------------------------------------------------- */
PetscErrorCode NonLinearSolverPETSc::FormFunction(SNES /*snes*/, Vec x,
                                                  Vec /*f*/, void * ctx) {
  auto * _this = reinterpret_cast<NonLinearSolverPETScCallback *>(ctx);
  _this->assembleResidual(x);
  return 0;
}

/* -------------------------------------------------------------------------- */
PetscErrorCode NonLinearSolverPETSc::FormJacobian(SNES /*snes*/, Vec x,
                                                  Mat /*J*/, Mat /*P*/,
                                                  void * ctx) {
  auto * _this = reinterpret_cast<NonLinearSolverPETScCallback *>(ctx);
  _this->assembleJacobian(x);
  return 0;
}

/* -------------------------------------------------------------------------- */
void NonLinearSolverPETSc::solve(SolverCallback & callback) {
  callback.beforeSolveStep();
  this->dof_manager.updateGlobalBlockedDofs();

  callback.assembleMatrix("J");
  auto & global_x =
      dynamic_cast<SolverVectorPETSc &>(dof_manager.getSolution());
  global_x.zero();

  if (not x) {
    x = std::make_unique<SolverVectorPETSc>(global_x, "temporary_solution");
  }

  *x = global_x;

  if (not ctx) {
    ctx = std::make_unique<NonLinearSolverPETScCallback>(
        dynamic_cast<DOFManagerPETSc &>(dof_manager), snes, *x);
  } else {
    ctx->reset();
  }

  ctx->setCallback(callback);
  ctx->setInitialSolution(global_x);

  auto & rhs = dynamic_cast<SolverVectorPETSc &>(dof_manager.getResidual());
  auto & J = dynamic_cast<SparseMatrixPETSc &>(dof_manager.getMatrix("J"));

  SNESSetFunction(snes, rhs, NonLinearSolverPETSc::FormFunction, ctx.get());
  SNESSetJacobian(snes, J, J, NonLinearSolverPETSc::FormJacobian, ctx.get());

  rhs.zero();

  callback.predictor();
  //  callback.assembleResidual();

  SNESSolve(snes, nullptr, *x);
  SNESGetConvergedReason(snes, &reason);
  SNESGetIterationNumber(snes, &n_iter);

  VecAXPY(global_x, -1.0, *x);
  dof_manager.splitSolutionPerDOFs();
  callback.corrector();

  bool converged = reason >= 0;
  callback.afterSolveStep(converged);

  if (not converged) {
    PetscReal atol;
    PetscReal rtol;
    PetscReal stol;
    PetscInt maxit;
    PetscInt maxf;

    SNESGetTolerances(snes, &atol, &rtol, &stol, &maxit, &maxf);
    AKANTU_CUSTOM_EXCEPTION(debug::SNESNotConvergedException(
        this->reason, this->n_iter, stol, atol, rtol, maxit));
  }
}

/* -------------------------------------------------------------------------- */
void NonLinearSolverPETSc::updateInternalParameters() {

  std::map<ID, ID> akantu_to_petsc_option = {{"max_iterations", "snes_max_it"},
                                             {"threshold", "snes_stol"}};

  for (auto && [param, param_akantu] : akantu_to_petsc_option) {
    Real value = this->get(param);
    PetscOptionsSetValue(nullptr, param_akantu.c_str(),
                         std::to_string(value).c_str());
  }
  SNESSetFromOptions(snes);
  PetscOptionsClear(nullptr);
}
/* --------------------------------------------------------------------------
 */
void NonLinearSolverPETSc::parseSection(const ParserSection & section) {
  auto parameters = section.getParameters();
  for (auto && param : range(parameters.first, parameters.second)) {
    PetscOptionsSetValue(nullptr, param.getName().c_str(),
                         param.getValue().c_str());
  }
  SNESSetFromOptions(snes);
  PetscOptionsClear(nullptr);
}

} // namespace akantu
