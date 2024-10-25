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
#include "sparse_solver.hh"
/* -------------------------------------------------------------------------- */
#include <petscksp.h>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_SOLVER_PETSC_HH_
#define AKANTU_SOLVER_PETSC_HH_

namespace akantu {
class SparseMatrixPETSc;
class DOFManagerPETSc;
} // namespace akantu

namespace akantu {

class SparseSolverPETSc : public SparseSolver {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  SparseSolverPETSc(DOFManager & dof_manager, const ID & matrix_id,
                    const ID & id = "solver_petsc");

  ~SparseSolverPETSc() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  DOFManagerPETSc & getDOFManager();
  /// create the solver context and set the matrices
  virtual void setOperators();
  void solve() override;

  void set(const std::string & name, std::any value);
  void parseSection(const ParserSection & section) override;
  void updateInternalParameters() override;

  void initialize() override;

private:
  /// PETSc linear solver
  KSP ksp;

  /// options to pass to petsc
  std::string petsc_options;
};

} // namespace akantu

#endif /* AKANTU_SOLVER_PETSC_HH_ */
