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
#include "solver_vector_petsc.hh"
#include "dof_manager_petsc.hh"
#include "mpi_communicator_data.hh"
/* -------------------------------------------------------------------------- */
#include <numeric>
#include <petscvec.h>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
SparseSolverVectorPETSc::SparseSolverVectorPETSc(DOFManagerPETSc & dof_manager,
                                                 const ID & id)
    : SparseSolverVector(dof_manager, id) {
  auto && mpi_comm = dof_manager.getMPIComm();
  PETSc_call(VecCreate, mpi_comm, &x);
  detail::PETScSetName(x, id);

  PETSc_call(VecSetFromOptions, x);

  auto local_system_size = dof_manager.getLocalSystemSize();
  auto nb_local_dofs = dof_manager.getPureLocalSystemSize();
  PETSc_call(VecSetSizes, x, nb_local_dofs, PETSC_DECIDE);

  VecType vec_type;
  PETSc_call(VecGetType, x, &vec_type);
  if (std::string(vec_type) == std::string(VECMPI)) {
    PetscInt lowest_gidx;
    PetscInt highest_gidx;
    PETSc_call(VecGetOwnershipRange, x, &lowest_gidx, &highest_gidx);

    std::vector<PetscInt> ghost_idx;
    for (auto && d : arange(local_system_size)) {
      int gidx = dof_manager.localToGlobalEquationNumber(d);
      if (gidx != -1) {
        if ((gidx < lowest_gidx) or (gidx >= highest_gidx)) {
          ghost_idx.push_back(gidx);
        }
      }
    }

    PETSc_call(VecMPISetGhost, x, ghost_idx.size(), ghost_idx.data());
  } else {
    std::vector<int> idx(nb_local_dofs);
    std::iota(idx.begin(), idx.end(), 0);
    ISLocalToGlobalMapping is;
    PETSc_call(ISLocalToGlobalMappingCreate, PETSC_COMM_SELF, 1, idx.size(),
               idx.data(), PETSC_COPY_VALUES, &is);
    PETSc_call(VecSetLocalToGlobalMapping, x, is);
    PETSc_call(ISLocalToGlobalMappingDestroy, &is);
  }
}

/* -------------------------------------------------------------------------- */
SparseSolverVectorPETSc::
    SparseSolverVectorPETSc( // NOLINT(bugprone-copy-constructor-init)
        const SparseSolverVectorPETSc & vector, const ID & id)
    : SparseSolverVector(vector, id) {
  if (vector.x != nullptr) {
    PETSc_call(VecDuplicate, vector.x, &x);
    PETSc_call(VecCopy, vector.x, x);
    detail::PETScSetName(x, id);
  }
}

/* -------------------------------------------------------------------------- */
void SparseSolverVectorPETSc::printself(std::ostream & stream,
                                        int indent) const {
  std::string space(indent, AKANTU_INDENT);
  stream << space << "SolverVectorPETSc [" << std::endl;
  stream << space << " + id: " << id << std::endl;
  PETSc_call(PetscViewerPushFormat, PETSC_VIEWER_STDOUT_WORLD,
             PETSC_VIEWER_ASCII_INDEX);
  PETSc_call(VecView, x, PETSC_VIEWER_STDOUT_WORLD);
  PETSc_call(PetscViewerPopFormat, PETSC_VIEWER_STDOUT_WORLD);
  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */
SparseSolverVectorPETSc::SparseSolverVectorPETSc(Vec x,
                                                 DOFManagerPETSc & dof_manager,
                                                 const ID & id)
    : SparseSolverVector(dof_manager, id) {
  PETSc_call(VecDuplicate, x, &this->x);

  PETSc_call(VecCopy, x, this->x);
  detail::PETScSetName(x, id);
}

/* -------------------------------------------------------------------------- */
SparseSolverVectorPETSc::~SparseSolverVectorPETSc() {
  if (x != nullptr) {
    PETSc_call(VecDestroy, &x);
  }
}

/* -------------------------------------------------------------------------- */
void SparseSolverVectorPETSc::resize() {
  // the arrays are destroyed and recreated in the dof manager
  // resize is so not implemented
}

/* -------------------------------------------------------------------------- */
void SparseSolverVectorPETSc::set(Real val) {
  PETSc_call(VecSet, x, val);
  applyModifications();
}

/* -------------------------------------------------------------------------- */
void SparseSolverVectorPETSc::applyModifications() {
  PETSc_call(VecAssemblyBegin, x);
  PETSc_call(VecAssemblyEnd, x);
  updateGhost();
}

/* -------------------------------------------------------------------------- */
void SparseSolverVectorPETSc::updateGhost() {
  Vec x_ghosted{nullptr};
  PETSc_call(VecGhostGetLocalForm, x, &x_ghosted);
  if (x_ghosted != nullptr) {
    PETSc_call(VecGhostUpdateBegin, x, INSERT_VALUES, SCATTER_FORWARD);
    PETSc_call(VecGhostUpdateEnd, x, INSERT_VALUES, SCATTER_FORWARD);
  }
  PETSc_call(VecGhostRestoreLocalForm, x, &x_ghosted);
}

/* -------------------------------------------------------------------------- */
void SparseSolverVectorPETSc::getValues(const Array<Int> & idx,
                                        Array<Real> & values) const {
  if (idx.empty()) {
    return;
  }

  ISLocalToGlobalMapping is_ltog_map;
  PETSc_call(VecGetLocalToGlobalMapping, x, &is_ltog_map);

  PetscInt n;
  Array<PetscInt> lidx(idx.size());
  PETSc_call(ISGlobalToLocalMappingApply, is_ltog_map, IS_GTOLM_MASK,
             idx.size(), idx.data(), &n, lidx.data());

  getValuesLocal(lidx, values);
}
/* -------------------------------------------------------------------------- */
void SparseSolverVectorPETSc::getValuesLocal(const Array<Int> & idx,
                                             Array<Real> & values) const {
  if (idx.empty()) {
    return;
  }

  Vec x_ghosted{nullptr};
  PETSc_call(VecGhostGetLocalForm, x, &x_ghosted);
  // VecScatterBegin(scatter, x, x_local, INSERT_VALUES, SCATTER_FORWARD);
  // VecScatterEnd(scatter, x, x_local, INSERT_VALUES, SCATTER_FORWARD);

  if (x_ghosted == nullptr) {
    const PetscScalar * array;
    PETSc_call(VecGetArrayRead, x, &array);

    for (auto && data : zip(idx, make_view(values))) {
      auto i = std::get<0>(data);
      if (i != -1) {
        std::get<1>(data) = array[i];
      }
    }

    PETSc_call(VecRestoreArrayRead, x, &array);
    return;
  }

  PETSc_call(VecSetOption, x_ghosted, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
  PETSc_call(VecGetValues, x_ghosted, idx.size(), idx.data(), values.data());
  PETSc_call(VecGhostRestoreLocalForm, x, &x_ghosted);
}

/* -------------------------------------------------------------------------- */
void SparseSolverVectorPETSc::addValues(const Array<Int> & gidx,
                                        const Array<Real> & values,
                                        Real scale_factor) {
  auto to_add = values.data();
  Array<Real> scaled_array(0, values.getNbComponent());
  if (scale_factor != 1.) {
    scaled_array.copy(values, false);
    scaled_array *= scale_factor;
    to_add = scaled_array.data();
  }

  PETSc_call(VecSetOption, x, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
  PETSc_call(VecSetValues, x, gidx.size(), gidx.data(), to_add, ADD_VALUES);

  applyModifications();
}

/* -------------------------------------------------------------------------- */
void SparseSolverVectorPETSc::addValuesLocal(const Array<Int> & lidx,
                                             const Array<Real> & values,
                                             Real scale_factor) {
  Vec x_ghosted{nullptr};
  PETSc_call(VecGhostGetLocalForm, x, &x_ghosted);

  if (x_ghosted == nullptr) {
    auto to_add = values.data();
    Array<Real> scaled_array;
    if (scale_factor != 1.) {
      scaled_array.copy(values, false);
      scaled_array *= scale_factor;
      to_add = scaled_array.data();
    }

    PETSc_call(VecSetOption, x, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
    PETSc_call(VecSetValuesLocal, x, lidx.size(), lidx.data(), to_add,
               ADD_VALUES);
    return;
  }

  PETSc_call(VecGhostRestoreLocalForm, x, &x_ghosted);

  ISLocalToGlobalMapping is_ltog_map;
  PETSc_call(VecGetLocalToGlobalMapping, x, &is_ltog_map);

  Array<Int> gidx(lidx.size());
  PETSc_call(ISLocalToGlobalMappingApply, is_ltog_map, lidx.size(), lidx.data(),
             gidx.data());
  addValues(gidx, values, scale_factor);
}

/* -------------------------------------------------------------------------- */
SparseSolverVectorPETSc::operator const Array<Real> &() const {
  const_cast<Array<Real> &>(cache).resize(local_size());

  auto xl = internal::make_petsc_local_vector(x);
  auto cachep = internal::make_petsc_wraped_vector(this->cache);

  PETSc_call(VecCopy, xl, cachep);
  return cache;
}

/* -------------------------------------------------------------------------- */
SparseSolverVectorPETSc &
SparseSolverVectorPETSc::operator=(const SparseSolverVectorPETSc & y) {
  if (size() != y.size()) {
    PETSc_call(VecDuplicate, y, &x);
  }

  PETSc_call(VecCopy, y.x, x);
  release_ = y.release_;
  return *this;
}

/* -------------------------------------------------------------------------- */
SparseSolverVector &
SparseSolverVectorPETSc::copy(const SparseSolverVector & y) {
  const auto & y_ = aka::as_type<SparseSolverVectorPETSc>(y);
  return operator=(y_);
}

/* -------------------------------------------------------------------------- */
SparseSolverVector &
SparseSolverVectorPETSc::operator+(const SparseSolverVector & y) {
  const auto & y_ = aka::as_type<SparseSolverVectorPETSc>(y);
  PETSc_call(VecAXPY, x, 1., y_.x);
  release_ = y_.release_;
  return *this;
}

bool SparseSolverVectorPETSc::isFinite() const {
  Real max, min;
  PETSc_call(VecMax, x, PETSC_NULLPTR, &max);
  PETSc_call(VecMin, x, PETSC_NULLPTR, &min);
  return std::isfinite(min) and std::isfinite(max);
}

} // namespace akantu
