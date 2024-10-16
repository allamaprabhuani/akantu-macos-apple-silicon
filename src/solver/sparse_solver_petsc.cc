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
    : SparseSolver(dof_manager, matrix_id, id) {
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
#if PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 5
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
  KSPSetFromOptions(ksp);
  try {
    KSPSetUp(ksp);
  } catch (std::exception & e) {
    AKANTU_EXCEPTION("KSP(PETSc sparse solver) initialization failed: is your "
                     "matrix singular ?"
                     << std::endl
                     << e.what());
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SparseSolverPETSc::solve() {
  Vec & rhs(getDOFManager()._getResidual());
  Vec & solution(getDOFManager()._getSolution());

  // auto && matrix = getDOFManager().getMatrix(matrix_id);

  this->setOperators();
  // MatView(matrix.getMat(), PETSC_VIEWER_STDOUT_WORLD);
  // VecView(rhs, PETSC_VIEWER_STDOUT_WORLD);

  KSPSolve(ksp, rhs, solution);
  // VecView(solution, PETSC_VIEWER_STDOUT_WORLD);

  this->dof_manager.splitSolutionPerDOFs();
}

/* -------------------------------------------------------------------------- */

DOFManagerPETSc & SparseSolverPETSc::getDOFManager() {
  return aka::as_type<DOFManagerPETSc &>(this->dof_manager);
}

/* -------------------------------------------------------------------------- */
void SparseSolverPETSc::updateInternalParameters() {

  PetscOptionsInsertString(nullptr, petsc_options.c_str());
  KSPSetFromOptions(ksp);
  PetscOptionsClear(nullptr);
}

/* -------------------------------------------------------------------------- */
void SparseSolverPETSc::parseSection(const ParserSection & section) {
  auto parameters = section.getParameters();
  for (auto && param : range(parameters.first, parameters.second)) {
    PetscOptionsSetValue(nullptr, param.getName().c_str(),
                         param.getValue().c_str());
  }
  KSPSetFromOptions(ksp);
  PetscOptionsClear(nullptr);
}
/* -------------------------------------------------------------------------- */

void SparseSolverPETSc::set(const std::string & name, std::any value) {
  if (this->hasParameter(name)) {
    SparseSolver::set(name, value);
  } else {
    try {
      // PetscOptionsView(nullptr, PETSC_VIEWER_STDOUT_WORLD);
      std::string option = std::any_cast<const char *>(value);
      // std::cout << name << ": " << option << std::endl;
      PetscOptionsSetValue(nullptr, ("-" + name).c_str(), option.c_str());
      KSPSetFromOptions(ksp);
      // PetscOptionsClear(nullptr);
    } catch (std::bad_any_cast & c) {
      std::cout << c.what() << std::endl;
      std::cout << value.type().name() << std::endl;
    }
  }
}

} // namespace akantu
