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
#include "aka_common.hh"
#include "parsable.hh"
/* -------------------------------------------------------------------------- */
#include <set>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_NON_LINEAR_SOLVER_HH_
#define AKANTU_NON_LINEAR_SOLVER_HH_

namespace akantu {
class DOFManager;
class SolverCallback;
struct ModelSolverOptions;
} // namespace akantu

namespace akantu {

class NonLinearSolver : public Parsable {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  NonLinearSolver(DOFManager & dof_manager,
                  const ModelSolverOptions & solver_options,
                  const ID & id = "non_linear_solver");
  ~NonLinearSolver() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// solve the system described by the jacobian matrix, and rhs contained in
  /// the dof manager
  virtual void solve(SolverCallback & callback) = 0;

protected:
  void checkIfTypeIsSupported();

  void assembleResidual(SolverCallback & callback);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  ID id;

  DOFManager & dof_manager;

  /// type of non linear solver
  NonLinearSolverType non_linear_solver_type;

  /// list of supported non linear solver types
  std::set<NonLinearSolverType> supported_type;

  /// specifies if the set param should be redirected
  bool has_internal_set_param{false};
};

namespace debug {
class NLSNotConvergedException : public Exception {
public:
  NLSNotConvergedException(Real threshold, Int niter, Real error)
      : Exception(""), threshold(threshold), niter(niter), error(error) {
    std::stringstream sstr;
    sstr << "The non linear solver did not converge after " << niter
         << " iterations (error " << error << " > " << threshold << ")";
    _info = sstr.str();
  }
  Real threshold;
  Int niter;
  Real error;
};
} // namespace debug

} // namespace akantu

#endif /* AKANTU_NON_LINEAR_SOLVER_HH_ */
