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
#include "non_linear_solver_tao.hh"
#include "aka_common.hh"
#include "dof_manager_petsc.hh"
#include "mpi_communicator_data.hh"
#include "solver_callback.hh"
#include "solver_vector_petsc.hh"
#include "sparse_matrix_petsc.hh"
/* -------------------------------------------------------------------------- */
#include <cstdlib>
#include <petscsnes.h>
#include <petsctao.h>
#include <petscvec.h>
/* -------------------------------------------------------------------------- */

namespace akantu {

NonLinearSolverTAO::NonLinearSolverTAO(
    DOFManagerPETSc & dof_manager, const ModelSolverOptions & solver_options,
    const ID & id)
    : NonLinearSolver(dof_manager, solver_options, id) {

  if (solver_options.sparse_solver_type != SparseSolverType::_petsc)
    AKANTU_EXCEPTION(
        "petsc non linear solver works only with petsc sparse solver");

  this->has_internal_set_param = true;

  supported_type.insert(NonLinearSolverType::_petsc_tao);

  this->checkIfTypeIsSupported();

  auto && mpi_comm = dof_manager.getMPIComm();

  TaoCreate(mpi_comm, &tao);

  TaoSetType(tao, TAOGPCG);
  TaoSetFromOptions(tao);

  this->registerParam("max_iterations", max_iterations, 10, _pat_parsmod,
                      "Max number of iterations");

  this->registerParam("threshold", convergence_criteria, 1e-10, _pat_parsmod,
                      "Threshold to consider results as converged");

  this->registerParam("convergence_type", convergence_criteria_type,
                      SolveConvergenceCriteria::_solution, _pat_parsmod,
                      "Type of convergence criteria");
}

/* -------------------------------------------------------------------------- */
NonLinearSolverTAO::~NonLinearSolverTAO() { TaoDestroy(&tao); }

/* -------------------------------------------------------------------------- */
void NonLinearSolverTAO::setTAOType(const ID & type) {
  PetscOptionsSetValue(NULL, "-tao_type", type.c_str());
}

/* -------------------------------------------------------------------------- */

void NonLinearSolverTAO::saveSolution() { AKANTU_TO_IMPLEMENT(); }

/* -------------------------------------------------------------------------- */

void NonLinearSolverTAO::restoreSolution() { AKANTU_TO_IMPLEMENT(); }

/* -------------------------------------------------------------------------- */

void NonLinearSolverTAO::corrector(Vec x) {

  auto & solution = aka::as_type<SolverVectorPETSc>(dof_manager.getSolution());
  if (x != solution.getVec())
    VecCopy(x, solution);

  dof_manager.splitSolutionPerDOFs();
  callback->restoreLastConvergedStep();
  callback->corrector();
}

/* -------------------------------------------------------------------------- */

void NonLinearSolverTAO::assembleResidual(Vec x, Vec f) {
  corrector(x);
  auto & residual =
      dynamic_cast<SolverVectorPETSc &>(dof_manager.getResidual());

  callback->assembleResidual();

  const auto & blocked_dofs = this->dof_manager.getGlobalBlockedDOFsIndexes();
  std::vector<Real> zeros_to_set(blocked_dofs.size());
  VecSetValuesLocal(residual, blocked_dofs.size(), blocked_dofs.data(),
                    zeros_to_set.data(), INSERT_VALUES);
  VecAssemblyBegin(residual);
  VecAssemblyEnd(residual);

  // for PETSc F is -akantu::residual
  VecScale(residual, -1);

  if (residual.getVec() != f) {
    VecCopy(residual, f);
  }
}

/* -------------------------------------------------------------------------- */

void NonLinearSolverTAO::assembleJacobian(Vec x, Mat J) {
  corrector(x);
  callback->assembleMatrix("J");
  auto & _J = aka::as_type<SparseMatrixPETSc>(dof_manager.getMatrix("J"));
  if (_J.getMat() != J) {
    // PetscPrint(_J);
    // PetscPrint(J);
    MatCopy(_J, J, SAME_NONZERO_PATTERN);
  }
}

/* -------------------------------------------------------------------------- */
void NonLinearSolverTAO::setBounds(const Vec lower_bound,
                                   const Vec upper_bound) {
  TaoSetVariableBounds(tao, lower_bound, upper_bound);
}

/* -------------------------------------------------------------------------- */
void NonLinearSolverTAO::computeObjectiveGradient(Vec x, PetscReal * obj,
                                                  Vec grad) {
  auto & J = aka::as_type<SparseMatrixPETSc>(dof_manager.getMatrix("J"));
  auto & rhs = aka::as_type<SolverVectorPETSc>(dof_manager.getResidual());

  assembleResidual(x, rhs);
  assembleJacobian(x, J);
  SolverVectorPETSc Jx(x, aka::as_type<DOFManagerPETSc>(this->dof_manager),
                       this->id + ":Kx");
  Real fx{};
  Real xJx{};

  MatMult(J, x, Jx);
  VecWAXPY(grad, 1, Jx, rhs);
  VecDot(rhs, x, &fx);
  VecDot(x, Jx, &xJx);
  *obj = 0.5 * xJx + fx;
}

/* -------------------------------------------------------------------------- */
PetscErrorCode NonLinearSolverTAO::FormFunctionGradient(Tao /*tao*/, Vec x,
                                                        PetscReal * obj,
                                                        Vec grad, void * ctx) {
  reinterpret_cast<NonLinearSolverTAO *>(ctx)->computeObjectiveGradient(x, obj,
                                                                        grad);

  return 0;
}

/* -------------------------------------------------------------------------- */
PetscErrorCode NonLinearSolverTAO::FormHessian(Tao /*tao*/, Vec x, Mat hess,
                                               Mat /*pre*/, void * ctx) {
  reinterpret_cast<NonLinearSolverTAO *>(ctx)->assembleJacobian(x, hess);
  return 0;
}

/* -------------------------------------------------------------------------- */
void NonLinearSolverTAO::solve(SolverCallback & callback) {
  callback.beforeSolveStep();
  this->dof_manager.updateGlobalBlockedDofs();

  callback.assembleMatrix("J");
  auto & x = aka::as_type<SolverVectorPETSc>(dof_manager.getSolution());
  x.zero();

  this->callback = &callback;

  auto & rhs = aka::as_type<SolverVectorPETSc>(dof_manager.getResidual());
  rhs.zero();

  auto & J = aka::as_type<SparseMatrixPETSc>(dof_manager.getMatrix("J"));

#if PETSC_VERSION_GE(3, 17, 0)
  TaoSetSolution(tao, x);
  TaoSetObjectiveAndGradient(tao, PETSC_NULLPTR,
                             NonLinearSolverTAO::FormFunctionGradient, this);
  TaoSetHessian(tao, J, J, NonLinearSolverTAO::FormHessian, this);
#else
  TaoSetInitialVector(tao, x);

  TaoSetObjectiveAndGradientRoutine(
      tao, NonLinearSolverTAO::FormFunctionGradient, this);
  TaoSetHessianRoutine(tao, J, J, NonLinearSolverTAO::FormHessian, this);

#endif

  callback.predictor();

  TaoSolve(tao);
  TaoGetConvergedReason(tao, &reason);
  TaoGetIterationNumber(tao, &n_iter);

  // access the model solution counter part: only for debug
  // auto & model_x = this->dof_manager.getDOFs("displacement");
  dof_manager.splitSolutionPerDOFs();
  callback.restoreLastConvergedStep();

  // \TODO: not efficient, TAO returns full solution, pseudo_time add it to
  // previous one, need to substract in model or script
  auto & us = this->dof_manager.getDOFs("damage");
  const auto & blocked_dofs = this->dof_manager.getBlockedDOFs("damage");
  for (auto && [u, bld] : zip(make_view(us), make_view(blocked_dofs))) {
    if (not bld) {
      u = 0;
    }
  }

  callback.corrector();
  bool converged = reason >= 0;
  callback.afterSolveStep(converged);

  if (not converged) {
    PetscReal atol;
    PetscReal rtol;
    PetscReal ttol;
    PetscInt maxit;

    TaoGetTolerances(tao, &atol, &rtol, &ttol);
    TaoGetMaximumIterations(tao, &maxit);
    AKANTU_CUSTOM_EXCEPTION(debug::TAONotConvergedException(
        this->reason, this->n_iter, ttol, atol, rtol, maxit));
  }
}

/* -------------------------------------------------------------------------- */
void NonLinearSolverTAO::updateInternalParameters() {

  std::map<ID, ID> akantu_to_petsc_option = {
      {"solver_type", "tao_type"},
      {"max_iterations", "tao_max_it"},
      {"absolute_threshold", "tao_gatol"},
      {"relative_threshold", "tao_grtol"}};

  for (auto && [param, param_akantu] : akantu_to_petsc_option) {
    auto & value = this->get(param);
    PetscOptionsSetValue(nullptr, ("-" + param_akantu).c_str(),
                         value.to_string().c_str());
  }
  TaoSetFromOptions(tao);
  PetscOptionsClear(nullptr);
}
/* -------------------------------------------------------------------------- */
void NonLinearSolverTAO::parseSection(const ParserSection & section) {
  auto parameters = section.getParameters();
  for (auto && param : range(parameters.first, parameters.second)) {
    PetscOptionsSetValue(nullptr, param.getName().c_str(),
                         param.getValue().c_str());
  }
  TaoSetFromOptions(tao);
  PetscOptionsClear(nullptr);
}
/* -------------------------------------------------------------------------- */

} // namespace akantu
