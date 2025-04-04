/**
 * Copyright (©) 2019-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "dof_manager_default.hh"
#include "solver_vector_default.hh"
/* -------------------------------------------------------------------------- */
#include <filesystem>
#include <fstream>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_SOLVER_VECTOR_DEFAULT_TMPL_HH_
#define AKANTU_SOLVER_VECTOR_DEFAULT_TMPL_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
inline SolverVectorArray::SolverVectorArray(DOFManagerDefault & dof_manager,
                                            const ID & id)
    : SolverVector(dof_manager, id) {}

/* -------------------------------------------------------------------------- */
inline SolverVectorArray::SolverVectorArray(const SolverVectorArray & vector,
                                            const ID & id)
    : SolverVector(vector, id) {}

/* -------------------------------------------------------------------------- */
template <class Array_>
SolverVector &
SolverVectorArrayTmpl<Array_>::operator+(const SolverVector & y) {
  const auto & y_ = aka::as_type<SolverVectorArray>(y);
  this->vector += y_.getVector();

  ++this->release_;
  return *this;
}

/* -------------------------------------------------------------------------- */
template <class Array_>
SolverVector & SolverVectorArrayTmpl<Array_>::copy(const SolverVector & y) {
  const auto & y_ = aka::as_type<SolverVectorArray>(y);
  this->vector.copy(y_.getVector());

  this->release_ = y.release();
  return *this;
}

/* -------------------------------------------------------------------------- */
template <class Array_> inline Int SolverVectorArrayTmpl<Array_>::size() const {
  return this->dof_manager.getSystemSize();
}

/* -------------------------------------------------------------------------- */
template <class Array_>
inline Int SolverVectorArrayTmpl<Array_>::localSize() const {
  return dof_manager.getLocalSystemSize();
}

/* -------------------------------------------------------------------------- */
template <class Array_>
inline void
SolverVectorArrayTmpl<Array_>::saveVector(const std::string & filename) const {
  std::filesystem::path file = filename;
  if (not file.has_extension()) {
    file.replace_extension(".mtx");
  }
  // open and set the properties of the stream
  std::ofstream outfile;
  auto range = arange(this->vector.size());
  auto size = std::count_if(range.begin(), range.end(), [&](auto n) {
    return dof_manager.isLocalOrMasterDOF(n);
  });

  outfile.open(file);
  outfile.precision(std::numeric_limits<Real>::digits10);
  outfile << "%%MatrixMarket matrix coordinate real general\n"
          << this->size() << " 1 " << size << "\n";
  for (auto && [i, a] : enumerate(this->vector)) {
    if (dof_manager.isLocalOrMasterDOF(i)) {
      outfile << (dof_manager.localToGlobalEquationNumber(i) + 1) << " 1 " << a
              << "\n";
    }
  }
  outfile.close();
}
} // namespace akantu

#endif /* AKANTU_SOLVER_VECTOR_DEFAULT_TMPL_HH_ */
