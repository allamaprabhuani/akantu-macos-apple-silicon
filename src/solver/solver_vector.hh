/**
 * Copyright (©) 2013-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_array.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_SOLVER_VECTOR_HH_
#define AKANTU_SOLVER_VECTOR_HH_

namespace akantu {
class DOFManager;
}

namespace akantu {

class SolverVector {
public:
  SolverVector(DOFManager & dof_manager, const ID & id = "solver_vector")
      : id(id), dof_manager(dof_manager) {}

  SolverVector(const SolverVector & vector, const ID & id = "solver_vector")
      : id(id), dof_manager(vector.dof_manager) {}

  virtual ~SolverVector() = default;

  // resize the vector to the size of the problem
  virtual void resize() = 0;

  // clear the vector
  virtual void set(Real val) = 0;
  void zero() { this->set({}); }

  virtual operator const Array<Real> &() const = 0;

  virtual Int size() const = 0;
  virtual Int localSize() const = 0;

  virtual SolverVector & operator+(const SolverVector & y) = 0;
  virtual SolverVector & copy(const SolverVector & y) = 0;

  Int & release() { return release_; }
  Int release() const { return release_; }

  virtual void printself(std::ostream & stream, int indent = 0) const = 0;

  virtual bool isFinite() const = 0;

  /// Returns `true` if `*this` is distributed or not.
  virtual bool isDistributed() const { return false; }

protected:
  ID id;

  /// Underlying dof manager
  DOFManager & dof_manager;

  Int release_{0};
};

inline std::ostream & operator<<(std::ostream & stream, SolverVector & _this) {
  _this.printself(stream);
  return stream;
}

} // namespace akantu

#endif /* AKANTU_SOLVER_VECTOR_HH_ */
