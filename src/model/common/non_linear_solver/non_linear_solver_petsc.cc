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
#include "petscsnes.h"
/* -------------------------------------------------------------------------- */

namespace akantu {

NonLinearSolverPETSc::NonLinearSolverPETSc(
    DOFManagerPETSc & dof_manager, const ModelSolverOptions & solver_options,
    const ID & id)
    : NonLinearSolver(dof_manager, solver_options, id) {

  if (solver_options.sparse_solver_type != SparseSolverType::_petsc)
    AKANTU_EXCEPTION(
        "petsc non linear solver works only with petsc sparse solver");

  this->has_internal_set_param = true;

  supported_type.insert(NonLinearSolverType::_petsc_snes);

  this->checkIfTypeIsSupported();

  auto && mpi_comm = dof_manager.getMPIComm();

  SNESCreate(mpi_comm, &snes);

  SNESSetType(snes, SNESNEWTONLS);
  SNESSetFromOptions(snes);

  this->registerParam("max_iterations", max_iterations, 10, _pat_parsmod,
                      "Max number of iterations");

  this->registerParam("threshold", convergence_criteria, 1e-10, _pat_parsmod,
                      "Threshold to consider results as converged");

  this->registerParam("convergence_type", convergence_criteria_type,
                      SolveConvergenceCriteria::_solution, _pat_parsmod,
                      "Type of convergence criteria");
}

/* -------------------------------------------------------------------------- */
NonLinearSolverPETSc::~NonLinearSolverPETSc() { SNESDestroy(&snes); }

/* -------------------------------------------------------------------------- */

void NonLinearSolverPETSc::saveSolution() { AKANTU_TO_IMPLEMENT(); }

/* -------------------------------------------------------------------------- */

void NonLinearSolverPETSc::restoreSolution() { AKANTU_TO_IMPLEMENT(); }

/* -------------------------------------------------------------------------- */

void NonLinearSolverPETSc::corrector(Vec x) {

  auto & solution = aka::as_type<SolverVectorPETSc>(dof_manager.getSolution());
  if (x != solution.getVec())
    VecCopy(x, solution);

  dof_manager.splitSolutionPerDOFs();
  callback->restoreLastConvergedStep();
  callback->corrector();
}

/* -------------------------------------------------------------------------- */

void NonLinearSolverPETSc::assembleResidual(Vec x, Vec f) {
  corrector(x);
  auto & residual =
      dynamic_cast<SolverVectorPETSc &>(dof_manager.getResidual());

  callback->assembleResidual();

  const auto & blocked_dofs = this->dof_manager.getGlobalBlockedDOFsIndexes();
  std::vector<Real> zeros_to_set(blocked_dofs.size());
  VecSetValuesLocal(residual, blocked_dofs.size(), blocked_dofs.data(),
                    zeros_to_set.data(), INSERT_VALUES);

  // for PETSc F is -akantu::residual
  VecScale(residual, -1);

  if (residual.getVec() != f) {
    VecCopy(residual, f);
  }
}

/* -------------------------------------------------------------------------- */

void NonLinearSolverPETSc::assembleJacobian(Vec x, Mat J) {
  corrector(x);
  callback->assembleMatrix("J");
  auto & _J = aka::as_type<SparseMatrixPETSc>(dof_manager.getMatrix("J"));
  if (_J.getMat() != J) {
    MatCopy(_J, J, SAME_NONZERO_PATTERN);
  }
}

/* -------------------------------------------------------------------------- */

PetscErrorCode NonLinearSolverPETSc::FormFunction(SNES /*snes*/, Vec x, Vec f,
                                                  void * ctx) {
  reinterpret_cast<NonLinearSolverPETSc *>(ctx)->assembleResidual(x, f);
  return 0;
}

/* -------------------------------------------------------------------------- */
PetscErrorCode NonLinearSolverPETSc::FormJacobian(SNES /*snes*/, Vec x, Mat J,
                                                  Mat /*P*/, void * ctx) {
  reinterpret_cast<NonLinearSolverPETSc *>(ctx)->assembleJacobian(x, J);
  return 0;
}

/* -------------------------------------------------------------------------- */
void NonLinearSolverPETSc::solve(SolverCallback & callback) {
  callback.beforeSolveStep();
  this->dof_manager.updateGlobalBlockedDofs();

  callback.assembleMatrix("J");
  auto & x = dynamic_cast<SolverVectorPETSc &>(dof_manager.getSolution());
  x.zero();

  this->callback = &callback;

  auto & rhs = aka::as_type<SolverVectorPETSc>(dof_manager.getResidual());
  rhs.zero();

  auto & J = aka::as_type<SparseMatrixPETSc>(dof_manager.getMatrix("J"));

  SNESSetFunction(snes, rhs, NonLinearSolverPETSc::FormFunction, this);
  SNESSetJacobian(snes, J, J, NonLinearSolverPETSc::FormJacobian, this);

  callback.predictor();

  // SNESView(snes, PETSC_VIEWER_STDOUT_WORLD);
  SNESSolve(snes, nullptr, x);
  SNESGetConvergedReason(snes, &reason);
  SNESGetIterationNumber(snes, &n_iter);

  // access the model solution counter part: only for debug
  // auto & model_x = this->dof_manager.getDOFs("displacement");
  dof_manager.splitSolutionPerDOFs();
  callback.restoreLastConvergedStep();
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
    auto & value = this->get(param);
    PetscOptionsSetValue(nullptr, ("-" + param_akantu).c_str(),
                         value.to_string().c_str());
  }
  SNESSetFromOptions(snes);
  PetscOptionsClear(nullptr);
}
/* -------------------------------------------------------------------------- */
void NonLinearSolverPETSc::parseSection(const ParserSection & section) {
  auto parameters = section.getParameters();
  for (auto && param : range(parameters.first, parameters.second)) {
    PetscOptionsSetValue(nullptr, param.getName().c_str(),
                         param.getValue().c_str());
  }
  SNESSetFromOptions(snes);
  PetscOptionsClear(nullptr);
}
/* -------------------------------------------------------------------------- */

} // namespace akantu
