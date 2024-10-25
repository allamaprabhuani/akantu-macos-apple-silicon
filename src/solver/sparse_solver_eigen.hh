/**
 * Copyright (©) 2010-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_SOLVER_EIGEN_HH_
#define AKANTU_SOLVER_EIGEN_HH_

namespace akantu {
class DOFManagerDefault;
class SparseMatrixAIJ;
} // namespace akantu

namespace akantu {

class SparseSolverEigen : public SparseSolver {
  /* ------------------------------------------------------------------------
   */
  /* Constructors/Destructors */
  /* ------------------------------------------------------------------------
   */
public:
  SparseSolverEigen(DOFManager & dof_manager, const ID & matrix_id,
                    const ID & id = "sparse_solver_eigen");

  /* ------------------------------------------------------------------------
   */
  /* Methods */
  /* ------------------------------------------------------------------------
   */
public:
  void initialize() override{};

  /// solve using residual and solution from the dof_manager
  void solve() override;

  DOFManagerDefault & getDOFManager();

private:
  void setA();

  /* ------------------------------------------------------------------------
   */
  /* Class Members */
  /* ------------------------------------------------------------------------
   */
private:
  /// matrix release at last solve
  Int last_profile_release{-1};

  /// matrix release at last solve
  Int last_value_release{-1};

  Eigen::SparseMatrix<Real> A;
  Eigen::SparseLU<Eigen::SparseMatrix<Real>> solver;
};

} // namespace akantu

#endif /* AKANTU_SOLVER_EIGEN_HH_ */
