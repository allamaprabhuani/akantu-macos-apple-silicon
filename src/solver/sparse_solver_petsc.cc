/**
 * Copyright (©) 2014-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "sparse_solver_petsc.hh"
#include "dof_manager_petsc.hh"
#include "mpi_communicator_data.hh"
#include "solver_vector_petsc.hh"
#include "sparse_matrix_petsc.hh"
/* -------------------------------------------------------------------------- */
#include <petscksp.h>
// #include <petscsys.h>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
SparseSolverPETSc::SparseSolverPETSc(DOFManager & dof_manager,
                                     const ID & matrix_id, const ID & id)
    : ParsablePETSc<SparseSolver, KSP>(ksp, KSPSetFromOptions, dof_manager,
                                       matrix_id, id) {
  auto && mpi_comm = getDOFManager().getMPIComm();

  /// create a solver context
  KSPCreate(mpi_comm, &this->ksp);
}
/* -------------------------------------------------------------------------- */
void SparseSolverPETSc::initialize() {}

/* -------------------------------------------------------------------------- */
SparseSolverPETSc::~SparseSolverPETSc() { KSPDestroy(&ksp); }

/* -------------------------------------------------------------------------- */
void SparseSolverPETSc::setOperators() {
  auto && matrix = getDOFManager().getMatrix(matrix_id);

  // set the matrix that defines the linear system and the matrix for
  // preconditioning (here they are the same)
#if PETSC_VERSION_GE(3, 5, 0)
  KSPSetOperators(ksp, matrix.getMat(), matrix.getMat());
#else
  KSPSetOperators(ksp, matrix.getMat(), matrix.getMat(), SAME_NONZERO_PATTERN);
#endif

  // If this is not called the solution vector is zeroed in the call to
  // KSPSolve().
  // @Nicolas: when this is done, the increment guess is completely wrong
  // (as it decays fast to zero in principle)
  // For some reason (I do not understand) a second solve in the phase field
  // problem is not made correctly with this setup=>to investigate
  // => uncommenting this line breaks test_phase_solid_coupling.cc when using a
  // sparse_petsc_solver
  // KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
  try {
    KSPSetUp(ksp);
  } catch (std::exception & e) {
    AKANTU_EXCEPTION("PETSc KSP sparse solver initialization failed: is your "
                     "matrix singular ?\n"
                     << e.what());
  }
}

/* -------------------------------------------------------------------------- */
void SparseSolverPETSc::solve() {
  Vec & rhs(getDOFManager()._getResidual());
  Vec & solution(getDOFManager()._getSolution());

  this->setOperators();

  KSPSolve(ksp, rhs, solution);

  this->dof_manager.splitSolutionPerDOFs();
}

/* -------------------------------------------------------------------------- */
DOFManagerPETSc & SparseSolverPETSc::getDOFManager() {
  return aka::as_type<DOFManagerPETSc &>(this->dof_manager);
}

} // namespace akantu
