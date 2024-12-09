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
#include "non_linear_solver.hh"
#include "solver_vector_petsc.hh"
/* -------------------------------------------------------------------------- */
#include <petsctao.h>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_NON_LINEAR_SOLVER_TAO_HH_
#define AKANTU_NON_LINEAR_SOLVER_TAO_HH_

namespace akantu {
class DOFManagerPETSc;
class SolverVectorPETSc;
} // namespace akantu

namespace akantu {

class NonLinearSolverTAO : public NonLinearSolver {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  NonLinearSolverTAO(DOFManagerPETSc & dof_manager,
                       const ModelSolverOptions & solver_options,
                       const ID & id = "non_linear_solver_tao");

  ~NonLinearSolverTAO() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// solve the system described by the jacobian matrix, and rhs contained in
  /// the dof manager
  void solve(SolverCallback & callback) override;

  /// parse the arguments from the input file
  void parseSection(const ParserSection & section) override;

  /// set solution bounds
  void setBounds(const Vec lower_bound, const Vec upper_bound);
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:

  static PetscErrorCode FormFunctionGradient(Tao tao, Vec x, PetscReal* obj,
                                             Vec grad, void* ctx);

  static PetscErrorCode FormHessian(Tao tao, Vec x, Mat hess, Mat pre,
                                    void* ctx);

  void computeObjectiveGradient(Vec x, PetscReal* obj, Vec grad);

  void corrector(Vec x);
  void assembleResidual(Vec x, Vec f);
  void assembleJacobian(Vec x, Mat J);
  void updateInternalParameters() override;

  void saveSolution();
  void restoreSolution();

  /// PETSc non linear solver
  Tao tao;
  TaoConvergedReason reason;

  SolverCallback * callback{nullptr};

  std::unique_ptr<SolverVectorPETSc> x;

  Int n_iter{0};
  Int max_iterations;
  /// Type of convergence criteria
  SolveConvergenceCriteria convergence_criteria_type;
  /// convergence threshold
  Real convergence_criteria;
};

namespace debug {
  class TAONotConvergedException : public NLSNotConvergedException {
  public:
    TAONotConvergedException(TaoConvergedReason reason, Int niter, Real error,
                              Real absolute_tolerance, Real relative_tolerance,
                              Int max_iterations)
        : NLSNotConvergedException(relative_tolerance, niter, error),
          reason(reason), absolute_tolerance(absolute_tolerance),
          max_iterations(max_iterations) {}
    TaoConvergedReason reason;
    Real absolute_tolerance;
    Int max_iterations;
  };
} // namespace debug

} // namespace akantu

#endif /* AKANTU_NON_LINEAR_SOLVER_PETSC_HH_ */
