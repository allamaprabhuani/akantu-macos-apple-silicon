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
#include "non_linear_solver.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_NON_LINEAR_SOLVER_NEWTON_RAPHSON_HH_
#define AKANTU_NON_LINEAR_SOLVER_NEWTON_RAPHSON_HH_

namespace akantu {
class DOFManagerDefault;
class SparseSolver;
class SolverVector;
} // namespace akantu

namespace akantu {

class NonLinearSolverNewtonRaphson : public NonLinearSolver {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  NonLinearSolverNewtonRaphson(
      DOFManager & dof_manager, const ModelSolverOptions & solver_otpions,
      const ID & id = "non_linear_solver_newton_raphson");
  ~NonLinearSolverNewtonRaphson() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// Function that solve the non linear system described by the dof manager and
  /// the solver callback functions
  void solve(SolverCallback & solver_callback) override;

  AKANTU_GET_MACRO_NOT_CONST(SparseSolver, *sparse_solver, SparseSolver &);
  AKANTU_GET_MACRO(SparseSolver, *sparse_solver, const SparseSolver &);

protected:
  //! function to call when a single sparse solve is demanded (1 iteration)
  void solve_linear(SolverCallback & solver_callback);

  /// test the convergence compare norm of array to convergence_criteria
  bool testConvergence(const SolverVector & solver_vector);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// Sparse solver used for the linear solves
  std::unique_ptr<SparseSolver> sparse_solver;

  /// Type of convergence criteria
  SolveConvergenceCriteria convergence_criteria_type;

  /// convergence threshold
  Real convergence_criteria;

  /// convergence threshold
  Real convergence_criteria_normalized;

  /// Max number of iterations
  int max_iterations;

  /// Number of iterations at last solve call
  int n_iter{0};

  /// Convergence error at last solve call
  Real error{0.};

  /// Did the last call to solve reached convergence
  bool converged{false};

  /// Force a re-computation of the jacobian matrix
  bool force_linear_recompute{true};

protected:
  /// flag do decide if one iteration only to be done
  bool linear{false};
};

/* -------------------------------------------------------------------------- */
class NonLinearSolverLinear : public NonLinearSolverNewtonRaphson {
public:
  NonLinearSolverLinear(DOFManager & dof_manager,
                        const ModelSolverOptions & solver_options,
                        const ID & id = "non_linear_solver_linear")
      : NonLinearSolverNewtonRaphson(dof_manager, solver_options, id) {
    this->linear = true;
  }
};

namespace debug {
class NLSNewtonRaphsonNotConvergedException : public NLSNotConvergedException {
public:
  NLSNewtonRaphsonNotConvergedException(
      Real threshold, Int niter, Real error,
      SolveConvergenceCriteria convergence_type)
      : NLSNotConvergedException(threshold, niter, error),
        convergence_type(convergence_type) {
    _info = _info + " - convergence type: " + std::to_string(convergence_type);
  }
  SolveConvergenceCriteria convergence_type;
};
} // namespace debug

} // namespace akantu

#endif /* AKANTU_NON_LINEAR_SOLVER_NEWTON_RAPHSON_HH_ */
