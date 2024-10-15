/**
 * Copyright (©) 2016-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "non_linear_solver_lumped.hh"
#include "communicator.hh"
#include "dof_manager_default.hh"
#include "solver_callback.hh"
#include "solver_vector_default.hh"
#if defined(AKANTU_USE_PETSC)
#include "solver_vector_petsc.hh"
#endif
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
NonLinearSolverLumped::NonLinearSolverLumped(
    DOFManager & dof_manager, const ModelSolverOptions & solver_options,
    const ID & id)
    : NonLinearSolver(dof_manager, solver_options, id) {
  this->supported_type.insert(NonLinearSolverType::_lumped);
  this->checkIfTypeIsSupported();

  this->registerParam("b_a2x", this->alpha, 1., _pat_parsmod,
                      "Conversion coefficient between x and A^{-1} b");
}

/* -------------------------------------------------------------------------- */
NonLinearSolverLumped::~NonLinearSolverLumped() = default;

/* ------------------------------------------------------------------------ */
void NonLinearSolverLumped::solve(SolverCallback & solver_callback) {
  solver_callback.beforeSolveStep();
  this->dof_manager.updateGlobalBlockedDofs();
  solver_callback.predictor();

  solver_callback.assembleResidual();

  SolverVector & x =
      aka::as_type<SolverVector>(this->dof_manager.getSolution());
  const auto & b = this->dof_manager.getResidual();

  x.resize();

  // const auto & blocked_dofs = this->dof_manager.getGlobalBlockedDOFs();
  const auto & A = this->dof_manager.getLumpedMatrix("M");

  // alpha is the conversion factor from from force/mass to acceleration needed
  // in model coupled with atomistic \todo find a way to define alpha per dof
  // type
  x.zero();

  if (1 == 2) {
#if defined(AKANTU_USE_PETSC)
  } else if (aka::is_of_type<SolverVectorPETSc>(x)) {
    auto & _x = aka::as_type<SolverVectorPETSc>(x);
    auto & _A = aka::as_type<SolverVectorPETSc>(A);
    auto & _b = aka::as_type<SolverVectorPETSc>(b);
    // VecView(aka::as_type<SolverVectorPETSc>(b).getVec(),
    //         PETSC_VIEWER_STDOUT_WORLD);
    NonLinearSolverLumped::solveLumped(_A, _x, _b, alpha);
#endif
  } else {
    auto & _x = aka::as_type<SolverVectorDefault>(x);
    NonLinearSolverLumped::solveLumped(A, _x, b, alpha);
  }

  this->dof_manager.splitSolutionPerDOFs();

  solver_callback.corrector();
  solver_callback.afterSolveStep(true);
}

/* -------------------------------------------------------------------------- */
void NonLinearSolverLumped::solveLumped(const Array<Real> & As,
                                        Array<Real> & xs,
                                        const Array<Real> & bs, Real alpha) {

  for (auto && [A, x, b] : zip(make_view(As), make_view(xs), make_view(bs))) {
    x = alpha * (b / A);
  }
}

/* -------------------------------------------------------------------------- */
#if defined(AKANTU_USE_PETSC)

void NonLinearSolverLumped::solveLumped(const SolverVectorPETSc & As,
                                        SolverVectorPETSc & xs,
                                        const SolverVectorPETSc & bs,
                                        Real alpha // ,
                                        // const Array<bool> & blocked_dofs
) {

  // Array<Real> _xs(As.size(), As.getNbComponent());
  // NonLinearSolverLumped::solveLumped(As, _xs, bs, alpha, blocked_dofs);
  // VecCopy(internal::make_petsc_wraped_vector(_xs), xs);

  VecPointwiseDivide(xs, bs, As);
  VecScale(xs, alpha);
  // VecView(xs.getVec(), PETSC_VIEWER_STDOUT_WORLD);
}
#endif
/* -------------------------------------------------------------------------- */

} // namespace akantu
