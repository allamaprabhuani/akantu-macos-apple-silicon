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
#include "parsable_petsc.hh"
#include "solver_vector_petsc.hh"
/* -------------------------------------------------------------------------- */
#include <petscsnes.h>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_NON_LINEAR_SOLVER_PETSC_HH_
#define AKANTU_NON_LINEAR_SOLVER_PETSC_HH_

namespace akantu {
class DOFManagerPETSc;
class SolverVectorPETSc;
} // namespace akantu

namespace akantu {

class NonLinearSolverPETSc : public ParsablePETSc<NonLinearSolver, SNES> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  NonLinearSolverPETSc(DOFManagerPETSc & dof_manager,
                       const ModelSolverOptions & solver_options,
                       const ID & id = "non_linear_solver_petsc");

  ~NonLinearSolverPETSc() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// solve the system described by the jacobian matrix, and rhs contained in
  /// the dof manager
  void solve(SolverCallback & callback) override;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  static PetscErrorCode FormFunction(SNES snes, Vec dx, Vec f, void * ctx);
  static PetscErrorCode FormJacobian(SNES snes, Vec dx, Mat J, Mat P,
                                     void * ctx);

  void corrector(Vec x);
  void assembleResidual(Vec x, Vec f);
  void assembleJacobian(Vec x, Mat J);

  void saveSolution();
  void restoreSolution();

private:
  /// PETSc non linear solver
  SNES snes{};

  SolverCallback * callback{nullptr};

  std::unique_ptr<SolverVectorPETSc> x;

  SNESConvergedReason reason{};
  std::string reason_str{};
  Int n_iter{0};
  PetscReal norm{}, xnorm{}, ynorm{};
  PetscInt petsc_na{1000};
  Array<PetscInt> petsc_its;
  Array<PetscReal> petsc_a;
  Real error{};
};

namespace debug {
class SNESNotConvergedException : public NLSNotConvergedException {
public:
  SNESNotConvergedException(std::string reason, Int niter, Int max_iterations,
                            Real absolute_tolerance, Real relative_tolerance,
                            Real solution_tolerance, Real absolute_norm,
                            Real relative_norm, Real solution_norm)
      : NLSNotConvergedException(relative_tolerance, niter, relative_norm),
        reason(reason), max_iterations(max_iterations),
        absolute_tolerance(absolute_tolerance),
        relative_tolerance(relative_tolerance),
        solution_tolerance(solution_tolerance), absolute_norm(absolute_norm),
        relative_norm(relative_norm), solution_norm(solution_norm) {
    std::stringstream sstr;
    sstr << "The PETSc solver did not converge for the reason " << reason
         << "\nLast norm:\n - atol " << absolute_tolerance << " <> "
         << absolute_norm << "\n - rtol " << relative_tolerance << " <> "
         << relative_norm << "\n - stol " << solution_tolerance << " <> "
         << solution_norm;
    this->_info = sstr.str();
  }
  std::string reason{};
  Int max_iterations{};
  Real absolute_tolerance{};
  Real relative_tolerance{};
  Real solution_tolerance{};
  Real absolute_norm{};
  Real relative_norm{};
  Real solution_norm{};
};
} // namespace debug

} // namespace akantu

#endif /* AKANTU_NON_LINEAR_SOLVER_PETSC_HH_ */
