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
SolverVectorPETSc::SolverVectorPETSc(DOFManagerPETSc & dof_manager,
                                     const ID & id)
    : SolverVector(dof_manager, id) {
  auto && mpi_comm = dof_manager.getMPIComm();
  VecCreate(mpi_comm, &x);
  detail::PETScSetName(x, id);

  resize();
  VecSetFromOptions(x);
}

/* -------------------------------------------------------------------------- */
SolverVectorPETSc::SolverVectorPETSc( // NOLINT(bugprone-copy-constructor-init)
    const SolverVectorPETSc & vector, const ID & id)
    : SolverVector(vector, id) {
  if (vector.x != nullptr) {
    VecDuplicate(vector.x, &x);
    VecCopy(vector.x, x);
    detail::PETScSetName(x, id);
  }
}

/* -------------------------------------------------------------------------- */
void SolverVectorPETSc::printself(std::ostream & stream, int indent) const {
  std::string space(indent, AKANTU_INDENT);
  stream << space << "SolverVectorPETSc [" << std::endl;
  stream << space << " + id: " << id << std::endl;
  PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_INDEX);
  VecView(x, PETSC_VIEWER_STDOUT_WORLD);
  PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD);
  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */
SolverVectorPETSc::SolverVectorPETSc(Vec x, DOFManagerPETSc & dof_manager,
                                     const ID & id)
    : SolverVector(dof_manager, id) {
  VecDuplicate(x, &this->x);

  VecCopy(x, this->x);
  detail::PETScSetName(x, id);
}

/* -------------------------------------------------------------------------- */
SolverVectorPETSc::~SolverVectorPETSc() {
  if (x != nullptr) {
    VecDestroy(&x);
  }
}

/* -------------------------------------------------------------------------- */
void SolverVectorPETSc::resize() {
  if (x != nullptr) {
    VecDestroy(&x);
    auto && mpi_comm = aka::as_type<DOFManagerPETSc>(dof_manager).getMPIComm();
    VecCreate(mpi_comm, &x);
    detail::PETScSetName(x, id);
    VecSetFromOptions(x);
  }

  auto nb_local_dofs = dof_manager.getPureLocalSystemSize();
  VecSetSizes(x, nb_local_dofs, dof_manager.getSystemSize());

  auto & is_ltog_mapping =
      aka::as_type<DOFManagerPETSc>(dof_manager).getISLocalToGlobalMapping();
  VecSetLocalToGlobalMapping(x, is_ltog_mapping);
}

/* -------------------------------------------------------------------------- */
void SolverVectorPETSc::set(Real val) {
  VecSet(x, val);
  applyModifications();
}

/* -------------------------------------------------------------------------- */
void SolverVectorPETSc::applyModifications() {
  VecAssemblyBegin(x);
  VecAssemblyEnd(x);
  updateGhost();
}

/* -------------------------------------------------------------------------- */
void SolverVectorPETSc::updateGhost() {
  Vec x_ghosted{nullptr};
  VecGhostGetLocalForm(x, &x_ghosted);
  if (x_ghosted != nullptr) {
    VecGhostUpdateBegin(x, INSERT_VALUES, SCATTER_FORWARD);
    VecGhostUpdateEnd(x, INSERT_VALUES, SCATTER_FORWARD);
  }
  VecGhostRestoreLocalForm(x, &x_ghosted);
}

/* -------------------------------------------------------------------------- */
void SolverVectorPETSc::getValues(const Array<Int> & idx,
                                  Array<Real> & values) const {
  if (idx.empty()) {
    return;
  }

  ISLocalToGlobalMapping is_ltog_map;
  VecGetLocalToGlobalMapping(x, &is_ltog_map);

  PetscInt n;
  Array<PetscInt> lidx(idx.size());
  ISGlobalToLocalMappingApply(is_ltog_map, IS_GTOLM_MASK, idx.size(),
                              idx.data(), &n, lidx.data());

  getValuesLocal(lidx, values);
}
/* -------------------------------------------------------------------------- */
void SolverVectorPETSc::getValuesLocal(const Array<Int> & idx,
                                       Array<Real> & values) const {
  if (idx.empty()) {
    return;
  }

  Vec x_ghosted{nullptr};
  VecGhostGetLocalForm(x, &x_ghosted);
  // VecScatterBegin(scatter, x, x_local, INSERT_VALUES, SCATTER_FORWARD);
  // VecScatterEnd(scatter, x, x_local, INSERT_VALUES, SCATTER_FORWARD);

  // VecView(x, PETSC_VIEWER_STDOUT_WORLD);

  if (x_ghosted == nullptr) {
    const PetscScalar * array;
    VecGetArrayRead(x, &array);

    for (auto && data : zip(idx, make_view(values))) {
      auto i = std::get<0>(data);
      if (i != -1) {
        std::get<1>(data) = array[i];
      }
    }

    VecRestoreArrayRead(x, &array);
    return;
  }

  VecSetOption(x_ghosted, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
  VecGetValues(x_ghosted, idx.size(), idx.data(), values.data());
  VecGhostRestoreLocalForm(x, &x_ghosted);
}

/* -------------------------------------------------------------------------- */
void SolverVectorPETSc::addValues(const Array<Int> & gidx,
                                  const Array<Real> & values,
                                  Real scale_factor) {
  auto to_add = values.data();
  Array<Real> scaled_array(0, values.getNbComponent());
  if (scale_factor != 1.) {
    scaled_array.copy(values, false);
    scaled_array *= scale_factor;
    to_add = scaled_array.data();
  }

  VecSetOption(x, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
  VecSetValues(x, gidx.size(), gidx.data(), to_add, ADD_VALUES);
  // VecView(x, PETSC_VIEWER_STDOUT_WORLD);
  applyModifications();
  // VecView(x, PETSC_VIEWER_STDOUT_WORLD);
}

/* -------------------------------------------------------------------------- */
void SolverVectorPETSc::addValuesLocal(const Array<Int> & lidx,
                                       const Array<Real> & values,
                                       Real scale_factor) {
  Vec x_ghosted{nullptr};
  VecGhostGetLocalForm(x, &x_ghosted);

  if (x_ghosted == nullptr) {
    auto to_add = values.data();
    Array<Real> scaled_array;
    if (scale_factor != 1.) {
      scaled_array.copy(values, false);
      scaled_array *= scale_factor;
      to_add = scaled_array.data();
    }

    VecSetOption(x, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
    VecSetValuesLocal(x, lidx.size(), lidx.data(), to_add, ADD_VALUES);
    return;
  }

  VecGhostRestoreLocalForm(x, &x_ghosted);

  ISLocalToGlobalMapping is_ltog_map;
  VecGetLocalToGlobalMapping(x, &is_ltog_map);

  Array<Int> gidx(lidx.size());
  ISLocalToGlobalMappingApply(is_ltog_map, lidx.size(), lidx.data(),
                              gidx.data());
  addValues(gidx, values, scale_factor);
}

/* -------------------------------------------------------------------------- */
SolverVectorPETSc::operator const Array<Real> &() const {
  const_cast<Array<Real> &>(cache).resize(local_size());

  auto xl = internal::make_petsc_local_vector(x);
  auto cachep = internal::make_petsc_wraped_vector(this->cache);

  VecCopy(xl, cachep);
  return cache;
}

/* -------------------------------------------------------------------------- */
SolverVectorPETSc & SolverVectorPETSc::operator=(const SolverVectorPETSc & y) {
  if (size() != y.size()) {
    VecDuplicate(y, &x);
  }

  VecCopy(y.x, x);
  release_ = y.release_;
  return *this;
}

/* -------------------------------------------------------------------------- */
SolverVector & SolverVectorPETSc::copy(const SolverVector & y) {
  const auto & y_ = aka::as_type<SolverVectorPETSc>(y);
  return operator=(y_);
}

/* -------------------------------------------------------------------------- */
SolverVector & SolverVectorPETSc::operator+(const SolverVector & y) {
  const auto & y_ = aka::as_type<SolverVectorPETSc>(y);
  VecAXPY(x, 1., y_.x);
  release_ = y_.release_;
  return *this;
}
/* -------------------------------------------------------------------------- */
bool SolverVectorPETSc::isFinite() const {
  Real max, min;
  VecMax(x, PETSC_NULLPTR, &max);
  VecMin(x, PETSC_NULLPTR, &min);
  return std::isfinite(min) and std::isfinite(max);
}
/* -------------------------------------------------------------------------- */
void PetscPrint(Vec x) { VecView(x, PETSC_VIEWER_STDOUT_WORLD); }

} // namespace akantu
