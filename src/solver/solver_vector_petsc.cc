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
#include "aka_array_printer.hh"
#include "aka_common.hh"
#include "dof_manager_petsc.hh"
#include "mesh.hh"
#include "mpi_communicator_data.hh"
/* -------------------------------------------------------------------------- */
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <ostream>
#include <petscis.h>
#include <petscistypes.h>
#include <petscsystypes.h>
#include <petscvec.h>
#include <petscviewer.h>
#include <string>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
SolverVectorPETSc::SolverVectorPETSc(DOFManagerPETSc & dof_manager,
                                     const ID & id)
    : SolverVector(dof_manager, id) {
  auto && mpi_comm = dof_manager.getMPIComm();
  VecCreate(mpi_comm, &x);
  detail::PETScSetName(x, id);

  VecSetFromOptions(x);
  dof_manager.setSolverVectorDataForParallelism(*this);
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
  auto & dof_manager_petsc = aka::as_type<DOFManagerPETSc>(dof_manager);
  if (x != nullptr) {
    VecDestroy(&x);
    auto && mpi_comm = dof_manager_petsc.getMPIComm();
    VecCreate(mpi_comm, &x);
    detail::PETScSetName(x, id);
    VecSetFromOptions(x);
  }

  dof_manager_petsc.setSolverVectorDataForParallelism(*this);
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
  // updateGhost();
}

/* -------------------------------------------------------------------------- */
void SolverVectorPETSc::updateGhost() const {
  // Vec x_ghosted{nullptr};
  // VecGhostGetLocalForm(x, &x_ghosted);
  // if (x_ghosted != nullptr) {
  VecGhostUpdateBegin(x, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd(x, INSERT_VALUES, SCATTER_FORWARD);
  //}
  // VecGhostRestoreLocalForm(x, &x_ghosted);
}

/* -------------------------------------------------------------------------- */
void SolverVectorPETSc::getValues(const Array<Int> & gidx,
                                  Array<Real> & values) const {
  if (gidx.empty()) {
    return;
  }

  ISLocalToGlobalMapping is_ltog_map;
  VecGetLocalToGlobalMapping(x, &is_ltog_map);

  PetscInt n;
  Array<PetscInt> lidx(gidx.size());
  ISGlobalToLocalMappingApply(is_ltog_map, IS_GTOLM_MASK, gidx.size(),
                              gidx.data(), &n, lidx.data());

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

  if (x_ghosted == nullptr) {
    Vec x_local{};

#if PETSC_VERSION_GE(3, 18, 0)
    VecCreateLocalVector(x, &x_local);
#else
    VecType roottype{};
    PetscInt n{};

    VecCreate(PETSC_COMM_SELF, &x_local);
    VecGetLocalSize(x, &n);
    VecSetSizes(x_local, n, n);
    VecGetBlockSize(x, &n);
    VecSetBlockSize(x_local, n);
    VecGetType(x, &roottype);
    VecSetType(x, roottype);
#endif

    VecGetLocalVector(x, x_local);

    VecSetOption(x_local, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
    VecGetValues(x_local, idx.size(), idx.data(), values.data());

    VecRestoreLocalVector(x, x_local);
    VecDestroy(&x_local);
    return;
  }

  updateGhost();

  VecSetOption(x_ghosted, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
  VecGetValues(x_ghosted, idx.size(), idx.data(), values.data());
  VecGhostRestoreLocalForm(x, &x_ghosted);
}

/* -------------------------------------------------------------------------- */
void SolverVectorPETSc::addValues(const Array<Int> & gidx,
                                  const Array<Real> & values,
                                  Real scale_factor) {
  const auto * to_add = values.data();
  Array<Real> scaled_array(0, values.getNbComponent());
  if (scale_factor != 1.) {
    scaled_array.copy(values, false);
    scaled_array *= scale_factor;
    to_add = scaled_array.data();
  }

  VecSetOption(x, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
  VecSetValues(x, gidx.size(), gidx.data(), to_add, ADD_VALUES);
}

/* -------------------------------------------------------------------------- */
void SolverVectorPETSc::addValuesLocal(const Array<Int> & lidx,
                                       const Array<Real> & values,
                                       Real scale_factor) {
  const auto * to_add = values.data();
  Array<Real> scaled_array;
  if (scale_factor != 1.) {
    scaled_array.copy(values, true);
    scaled_array *= scale_factor;
    to_add = scaled_array.data();
  }

  Vec x_ghosted{nullptr};
  VecGhostGetLocalForm(x, &x_ghosted);

  if (x_ghosted == nullptr) {
    VecSetOption(x, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
    VecSetValuesLocal(x, lidx.size(), lidx.data(), to_add, ADD_VALUES);
    return;
  }

  VecSetOption(x_ghosted, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
  VecSetValuesLocal(x_ghosted, lidx.size(), lidx.data(), to_add, ADD_VALUES);

  VecGhostRestoreLocalForm(x, &x_ghosted);
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
  Real max{};
  Real min{};
  VecMax(x, nullptr, &max);
  VecMin(x, nullptr, &min);
  return std::isfinite(min) and std::isfinite(max);
}
/* -------------------------------------------------------------------------- */
void PetscPrint(Vec x) { VecView(x, PETSC_VIEWER_STDOUT_WORLD); }

/* -------------------------------------------------------------------------- */
void SolverVectorPETSc::saveVector(const std::string & filename) const {

  AKANTU_DEBUG_IN();

  auto & comm = dof_manager.getCommunicator();

  std::filesystem::path file = filename;
  if (comm.getNbProc() > 1) {
    file.replace_extension("");
    file += "_rank-" + std::to_string(comm.whoAmI()) + ".mtx";
  }

  if (not file.has_extension()) {
    file.replace_extension(".mtx");
  }

  // open and set the properties of the stream
  std::ofstream outfile;
  // auto range = arange(this->localSize());
  // auto size = std::count_if(range.begin(), range.end(), [&](auto n) {
  //   return dof_manager.isLocalOrMasterDOF(n);
  // });

  outfile.open(file);
  outfile.precision(std::numeric_limits<Real>::digits10);
  outfile << "%%MatrixMarket matrix coordinate real general\n"
          << this->size() << " 1 " << this->localSize() << "\n";
  // VecView(x, PETSC_VIEWER_STDOUT_WORLD);
  Int start{};
  const Real * array{};
  VecGetOwnershipRange(x, &start, nullptr);
  VecGetArrayRead(x, &array);
  // const Array<Real> & vector = *this;
  // ArrayPrinter<Array<Real>> printer(vector);
  // printer.printself(std::cout);
  for (Int i = 0; i < this->localSize(); ++i) {
    outfile << (start + i + 1) << " 1 " << array[i] << "\n";
    //   if (dof_manager.isLocalOrMasterDOF(i)) {
    //     outfile << (dof_manager.localToGlobalEquationNumber(i) + 1) << " 1 "
    //     << a
    //             << "\n";
    //   }
  }
  outfile.close();
  AKANTU_DEBUG_OUT();
}

} // namespace akantu
