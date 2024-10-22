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
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_SOLVER_CALLBACK_HH_
#define AKANTU_SOLVER_CALLBACK_HH_

namespace akantu {
class DOFManager;
}

namespace akantu {

class SolverCallback {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  explicit SolverCallback(DOFManager & dof_manager);
  explicit SolverCallback() = default;
  /* ------------------------------------------------------------------------ */
  SolverCallback(const SolverCallback & other) = default;
  SolverCallback(SolverCallback && other) = default;
  SolverCallback & operator=(const SolverCallback & other) = default;
  SolverCallback & operator=(SolverCallback && other) = default;
  /* ------------------------------------------------------------------------ */
  virtual ~SolverCallback() = default;

protected:
  void setDOFManager(DOFManager & dof_manager);
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// get the type of matrix needed
  [[nodiscard]] virtual MatrixType
  getMatrixType(const ID & /*matrix_id*/) const {
    return _mt_not_defined;
  }

  /// callback to assemble a Matrix
  virtual void assembleMatrix(const ID & /*matrix_id*/) {
    AKANTU_TO_IMPLEMENT();
  }

  /// callback to assemble a lumped Matrix
  virtual void assembleLumpedMatrix(const ID & /*matrix_id*/) {
    AKANTU_TO_IMPLEMENT();
  }

  /// callback to assemble the residual (rhs)
  virtual void assembleResidual() { AKANTU_TO_IMPLEMENT(); }

  /// callback to assemble the rhs parts, (e.g. internal_forces +
  /// external_forces)
  virtual void assembleResidual(const ID & /*residual_part*/) {}

  /* ------------------------------------------------------------------------ */
  /* Dynamic simulations part                                                 */
  /* ------------------------------------------------------------------------ */
  /// callback for the predictor (in case of dynamic simulation)
  virtual void predictor() {}

  /// callback for the corrector (in case of dynamic simulation)
  virtual void corrector() {}

  /// tells if the residual can be computed in separated parts
  [[nodiscard]] virtual bool canSplitResidual() const { return false; }

  /* ------------------------------------------------------------------------ */
  /* management callbacks                                                     */
  /* ------------------------------------------------------------------------ */
  virtual void beforeSolveStep() {}
  virtual void afterSolveStep(bool /*converged*/ = true) {}
  virtual void restoreLastConvergedStep() {}

  [[nodiscard]] DOFManager & getSCDOFManager() { return *sc_dof_manager; }

protected:
  /// DOFManager prefixed to avoid collision in multiple inheritance cases
  DOFManager * sc_dof_manager{nullptr};
};

namespace debug {
  class SolverCallbackResidualPartUnknown : public Exception {
  public:
    SolverCallbackResidualPartUnknown(const ID & residual_part)
        : Exception(residual_part + " is not known here.") {}
  };
} // namespace debug

/* -------------------------------------------------------------------------- */
class InterceptSolverCallback : public SolverCallback {
public:
  InterceptSolverCallback(SolverCallback & solver_callback)
      : solver_callback(solver_callback) {}

  /// get the type of matrix needed
  [[nodiscard]] MatrixType getMatrixType(const ID & matrix_id) const override {
    return solver_callback.getMatrixType(matrix_id);
  }

  /// callback to assemble a Matrix
  void assembleMatrix(const ID & matrix_id) override {
    solver_callback.assembleMatrix(matrix_id);
  }

  /// callback to assemble a lumped Matrix
  void assembleLumpedMatrix(const ID & matrix_id) override {
    solver_callback.assembleLumpedMatrix(matrix_id);
  }

  /// callback to assemble the residual (rhs)
  void assembleResidual() override { solver_callback.assembleResidual(); }

  /// callback to assemble the rhs parts, (e.g. internal_forces +
  /// external_forces)
  void assembleResidual(const ID & residual_part) override {
    solver_callback.assembleResidual(residual_part);
  }

  /* ------------------------------------------------------------------------ */
  /* Dynamic simulations part                                                 */
  /* ------------------------------------------------------------------------ */
  /// callback for the predictor (in case of dynamic simulation)
  void predictor() override { solver_callback.predictor(); }

  /// callback for the corrector (in case of dynamic simulation)
  void corrector() override { solver_callback.corrector(); }

  /// tells if the residual can be computed in separated parts
  [[nodiscard]] bool canSplitResidual() const override {
    return solver_callback.canSplitResidual();
  }

  /* ------------------------------------------------------------------------ */
  /* management callbacks                                                     */
  /* ------------------------------------------------------------------------ */
  void beforeSolveStep() override { solver_callback.beforeSolveStep(); }
  void afterSolveStep(bool converged = true) override {
    solver_callback.afterSolveStep(converged);
  }

protected:
  SolverCallback & solver_callback;
};

} // namespace akantu

#endif /* AKANTU_SOLVER_CALLBACK_HH_ */
