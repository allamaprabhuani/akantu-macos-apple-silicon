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
#include "communicator.hh"
#include "dof_manager.hh"
#include "mesh.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
SparseSolver::SparseSolver(DOFManager & dof_manager, const ID & matrix_id,
                           const ID & id)
    : Parsable(ParserType::_solver, id), dof_manager(dof_manager),
      matrix_id(matrix_id), communicator(dof_manager.getCommunicator()) {
  AKANTU_DEBUG_IN();

  // OK this is fishy...
  this->communicator.registerEventHandler(*this);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
SparseSolver::~SparseSolver() {
  AKANTU_DEBUG_IN();

  // this->destroyInternalData();
  this->communicator.unregisterEventHandler(*this);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SparseSolver::beforeStaticSolverDestroy() {
  AKANTU_DEBUG_IN();

  try {
    this->destroyInternalData();
  } catch (...) {
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SparseSolver::createSynchronizerRegistry() {
  // this->synch_registry = new SynchronizerRegistry(this);
}

void SparseSolver::onCommunicatorFinalize() { this->destroyInternalData(); }

} // namespace akantu
