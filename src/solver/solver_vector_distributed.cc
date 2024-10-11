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
#include "solver_vector_distributed.hh"
#include "dof_manager_default.hh"
#include "dof_synchronizer.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
SolverVectorDistributed::SolverVectorDistributed(
    DOFManagerDefault & dof_manager, const ID & id)
    : SolverVectorDefault(dof_manager, id) {}

/* -------------------------------------------------------------------------- */
SolverVectorDistributed::SolverVectorDistributed(
    const SolverVectorDefault & vector, const ID & id)
    : SolverVectorDefault(vector, id) {}

/* -------------------------------------------------------------------------- */
Array<Real> & SolverVectorDistributed::getGlobalVector() {
  auto & synchronizer =
      dynamic_cast<DOFManagerDefault &>(dof_manager).getSynchronizer();

  if (not this->global_vector) {
    this->global_vector =
        std::make_unique<Array<Real>>(0, 1, "global_residual");
  }

  if (synchronizer.getCommunicator().whoAmI() == 0) {
    this->global_vector->resize(dof_manager.getSystemSize());
    synchronizer.gather(this->vector, *this->global_vector);
  } else {
    synchronizer.gather(this->vector);
  }

  return *this->global_vector;
}

/* -------------------------------------------------------------------------- */
void SolverVectorDistributed::setGlobalVector(const Array<Real> & solution) {
  auto & synchronizer =
      dynamic_cast<DOFManagerDefault &>(dof_manager).getSynchronizer();
  if (synchronizer.getCommunicator().whoAmI() == 0) {
    synchronizer.scatter(this->vector, solution);
  } else {
    synchronizer.scatter(this->vector);
  }
}

/* -------------------------------------------------------------------------- */
bool SolverVectorDistributed::isFinite() const {
  auto & synchronizer =
      dynamic_cast<DOFManagerDefault &>(dof_manager).getSynchronizer();
  bool is_finite = this->vector.isFinite();
  synchronizer.getCommunicator().allReduce(is_finite,
                                           SynchronizerOperation::_land);
  return is_finite;
}

/* -------------------------------------------------------------------------- */

} // namespace akantu
