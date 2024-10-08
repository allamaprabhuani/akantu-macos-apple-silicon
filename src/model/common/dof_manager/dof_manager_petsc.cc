/**
 * Copyright (©) 2015-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "dof_manager_petsc.hh"
#include "aka_iterators.hh"
#include "communicator.hh"
#include "cppargparse.hh"
#include "non_linear_solver_default.hh"
#include "non_linear_solver_petsc.hh"
#include "solver_vector_petsc.hh"
#include "sparse_matrix_petsc.hh"
#include "time_step_solver_default.hh"
#if defined(AKANTU_USE_MPI)
#include "mpi_communicator_data.hh"
#endif
/* -------------------------------------------------------------------------- */
#include <petscis.h>
#include <petscsys.h>
/* -------------------------------------------------------------------------- */

namespace akantu {

class PETScSingleton {
private:
  PETScSingleton() {
    _PETSc_call(PetscInitialized, &is_initialized);

    if (is_initialized == 0U) {
      cppargparse::ArgumentParser & argparser = getStaticArgumentParser();
      int & argc = argparser.getArgC();
      char **& argv = argparser.getArgV();
      _PETSc_call(PetscInitialize, &argc, &argv, nullptr, nullptr);
      _PETSc_call(
          PetscPopErrorHandler, ); // remove the default PETSc signal handler
      _PETSc_call(PetscPushErrorHandler, petscErrorHandler, nullptr);
    }
  }

public:
  PETScSingleton(const PETScSingleton &) = delete;
  PETScSingleton & operator=(const PETScSingleton &) = delete;

  ~PETScSingleton() {
    if (is_initialized == 0U) {
      int mpi_finalized;
      MPI_Finalized(&mpi_finalized);
      if (not mpi_finalized)
        PetscFinalize();
    }
  }

  static PETScSingleton & getInstance() {
    static PETScSingleton instance;
    return instance;
  }

private:
  PetscBool is_initialized;
};

/* -------------------------------------------------------------------------- */
DOFManagerPETSc::DOFDataPETSc::DOFDataPETSc(const ID & dof_id)
    : DOFData(dof_id) {}

/* -------------------------------------------------------------------------- */
DOFManagerPETSc::DOFManagerPETSc(const ID & id) : DOFManager(id) { init(); }

/* -------------------------------------------------------------------------- */
DOFManagerPETSc::DOFManagerPETSc(Mesh & mesh, const ID & id)
    : DOFManager(mesh, id) {
  init();
}

/* -------------------------------------------------------------------------- */
void DOFManagerPETSc::init() {
  // check if the akantu types and PETSc one are consistant
  static_assert(sizeof(Int) == sizeof(PetscInt),
                "The integer type of Akantu does not match the one from PETSc");
  static_assert(sizeof(Real) == sizeof(PetscReal),
                "The integer type of Akantu does not match the one from PETSc");

#if defined(AKANTU_USE_MPI)

  const auto & mpi_data =
      aka::as_type<MPICommunicatorData>(communicator.getCommunicatorData());
  MPI_Comm mpi_comm = mpi_data.getMPICommunicator();
  this->mpi_communicator = mpi_comm;
#else
  this->mpi_communicator = PETSC_COMM_SELF;
#endif

  PETScSingleton & instance [[gnu::unused]] = PETScSingleton::getInstance();
}

/* -------------------------------------------------------------------------- */
auto DOFManagerPETSc::getNewDOFData(const ID & dof_id)
    -> std::unique_ptr<DOFData> {
  return std::make_unique<DOFDataPETSc>(dof_id);
}

/* -------------------------------------------------------------------------- */
std::tuple<Int, Int, Int>
DOFManagerPETSc::registerDOFsInternal(const ID & dof_id,
                                      Array<Real> & dofs_array) {
  dofs_ids.push_back(dof_id);
  auto ret = DOFManager::registerDOFsInternal(dof_id, dofs_array);
  UInt nb_dofs;
  UInt nb_pure_local_dofs;
  std::tie(nb_dofs, nb_pure_local_dofs, std::ignore) = ret;

  auto && vector =
      std::make_unique<SparseSolverVectorPETSc>(*this, id + ":solution");
  auto * x = vector->getVec();
  VecGetLocalToGlobalMapping(x, &is_ltog_map);

  // redoing the indexes based on the petsc numbering
  for (auto & dof_id : dofs_ids) {
    auto & dof_data = this->getDOFDataTyped<DOFDataPETSc>(dof_id);

    Array<PetscInt> gidx(dof_data.local_equation_number.size());
    for (auto && data : zip(dof_data.local_equation_number, gidx)) {
      std::get<1>(data) = localToGlobalEquationNumber(std::get<0>(data));
    }

    auto & lidx = dof_data.local_equation_number_petsc;
    if (is_ltog_map != nullptr) {
      lidx.resize(gidx.size());

      PetscInt n;
      ISGlobalToLocalMappingApply(is_ltog_map, IS_GTOLM_MASK, gidx.size(),
                                  gidx.data(), &n, lidx.data());
    }
  }

  residual =
      std::make_unique<SparseSolverVectorPETSc>(*vector, id + ":residual");
  data_cache =
      std::make_unique<SparseSolverVectorPETSc>(*vector, id + ":data_cache");
  solution = std::move(vector);

  for (auto & mat : matrices) {
    auto & A = this->getMatrix(mat.first);
    A.resize();
  }

  return ret;
}

/* -------------------------------------------------------------------------- */
void DOFManagerPETSc::assembleToGlobalArray(
    const ID & dof_id, const Array<Real> & array_to_assemble,
    SparseSolverVector & global_array, Real scale_factor) {
  const auto & dof_data = getDOFDataTyped<DOFDataPETSc>(dof_id);
  auto & g = aka::as_type<SparseSolverVectorPETSc>(global_array);

  AKANTU_DEBUG_ASSERT(array_to_assemble.size() *
                              array_to_assemble.getNbComponent() ==
                          dof_data.local_nb_dofs,
                      "The array to assemble does not have the proper size");

  g.addValuesLocal(dof_data.local_equation_number_petsc, array_to_assemble,
                   scale_factor);
}

/* -------------------------------------------------------------------------- */
void DOFManagerPETSc::getArrayPerDOFs(const ID & dof_id,
                                      const SparseSolverVector & global_array,
                                      Array<Real> & local) {
  const auto & dof_data = getDOFDataTyped<DOFDataPETSc>(dof_id);
  const auto & petsc_vector =
      aka::as_type<SparseSolverVectorPETSc>(global_array);

  AKANTU_DEBUG_ASSERT(
      local.size() * local.getNbComponent() == dof_data.local_nb_dofs,
      "The array to get the values does not have the proper size");

  petsc_vector.getValuesLocal(dof_data.local_equation_number_petsc, local);
}

/* -------------------------------------------------------------------------- */
void DOFManagerPETSc::assembleElementalMatricesToMatrix(
    const ID & matrix_id, const ID & dof_id, const Array<Real> & elementary_mat,
    ElementType type, GhostType ghost_type,
    const MatrixType & elemental_matrix_type,
    const Array<Int> & filter_elements) {
  auto & A = getMatrix(matrix_id);
  DOFManager::assembleElementalMatricesToMatrix_(
      A, dof_id, elementary_mat, type, ghost_type, elemental_matrix_type,
      filter_elements);

  A.applyModifications();
  // MatView(A.getMat(), PETSC_VIEWER_STDOUT_WORLD);
}

/* -------------------------------------------------------------------------- */
void DOFManagerPETSc::assemblePreassembledMatrix(
    const ID & matrix_id, const TermsToAssemble & terms) {
  auto & A = getMatrix(matrix_id);
  DOFManager::assemblePreassembledMatrix_(A, terms);

  A.applyModifications();
}

/* -------------------------------------------------------------------------- */
void DOFManagerPETSc::assembleMatMulVectToArray(const ID & dof_id,
                                                const ID & A_id,
                                                const Array<Real> & x,
                                                Array<Real> & array,
                                                Real scale_factor) {
  DOFManager::assembleMatMulVectToArray_<SparseSolverVectorPETSc>(
      dof_id, A_id, x, array, scale_factor);
}

/* -------------------------------------------------------------------------- */
void DOFManagerPETSc::makeConsistentForPeriodicity(
    const ID & /*dof_id*/, SparseSolverVector & /*array*/) {}

/* -------------------------------------------------------------------------- */
NonLinearSolver &
DOFManagerPETSc::getNewNonLinearSolver(const ID & id,
                                       const NonLinearSolverType & type,
                                       const SparseSolverType & s_type) {
  switch (type) {
  case NonLinearSolverType::_newton_raphson:
    /* FALLTHRU */
    /* [[fallthrough]]; un-comment when compiler will get it */
  case NonLinearSolverType::_newton_raphson_contact:
  case NonLinearSolverType::_newton_raphson_modified: {
    return this->registerNonLinearSolver<NonLinearSolverNewtonRaphson>(
        *this, id, type, s_type);
  }
  case NonLinearSolverType::_linear: {
    return this->registerNonLinearSolver<NonLinearSolverLinear>(*this, id, type,
                                                                s_type);
  }
  case NonLinearSolverType::_lumped: {
    return this->registerNonLinearSolver<NonLinearSolverLumped>(*this, id, type,
                                                                s_type);
  }
  default:
    return this->registerNonLinearSolver<NonLinearSolverPETSc>(*this, id, type,
                                                               s_type);
  }
}

/* -------------------------------------------------------------------------- */
TimeStepSolver & DOFManagerPETSc::getNewTimeStepSolver(
    const ID & id, const TimeStepSolverType & type,
    NonLinearSolver & non_linear_solver, SolverCallback & callback) {
  return this->registerTimeStepSolver<TimeStepSolverDefault>(
      *this, id, type, non_linear_solver, callback);
}

/* -------------------------------------------------------------------------- */
SparseMatrix & DOFManagerPETSc::getNewMatrix(const ID & id,
                                             const MatrixType & matrix_type) {
  return this->registerSparseMatrix<SparseMatrixPETSc>(*this, id, matrix_type);
}

/* -------------------------------------------------------------------------- */
SparseMatrix & DOFManagerPETSc::getNewMatrix(const ID & id,
                                             const ID & matrix_to_copy_id) {
  return this->registerSparseMatrix<SparseMatrixPETSc>(id, matrix_to_copy_id);
}

/* -------------------------------------------------------------------------- */
SparseMatrixPETSc & DOFManagerPETSc::getMatrix(const ID & id) {
  auto & matrix = DOFManager::getMatrix(id);
  return aka::as_type<SparseMatrixPETSc>(matrix);
}

/* -------------------------------------------------------------------------- */
SparseSolverVector & DOFManagerPETSc::getNewLumpedMatrix(const ID & id) {
  return this->registerLumpedMatrix<SparseSolverVectorPETSc>(*this, id);
}

/* -------------------------------------------------------------------------- */
SparseSolverVectorPETSc & DOFManagerPETSc::_getSolution() {
  return aka::as_type<SparseSolverVectorPETSc>(*this->solution);
}

const SparseSolverVectorPETSc & DOFManagerPETSc::_getSolution() const {
  return aka::as_type<SparseSolverVectorPETSc>(*this->solution);
}

SparseSolverVectorPETSc & DOFManagerPETSc::_getResidual() {
  return aka::as_type<SparseSolverVectorPETSc>(*this->residual);
}

const SparseSolverVectorPETSc & DOFManagerPETSc::_getResidual() const {
  return aka::as_type<SparseSolverVectorPETSc>(*this->residual);
}

/* -------------------------------------------------------------------------- */
void DOFManagerPETSc::assembleLumpedMatMulVectToResidual(const ID & dof_id,
                                                         const ID & A_id,
                                                         const Array<Real> & x,
                                                         Real scale_factor) {
  const auto & A =
      aka::as_type<SparseSolverVectorPETSc>(this->getLumpedMatrix(A_id))
          .getVec();
  auto & cache = aka::as_type<SparseSolverVectorPETSc>(*this->data_cache);

  // int sz;
  // VecGetSize(A, &sz);
  // std::cout << "AAAAAA: A " << sz << std::endl;

  // VecGetSize(cache, &sz);
  // std::cout << "AAAAAA: cache " << sz << std::endl;

  cache.zero();
  this->assembleToGlobalArray(dof_id, x, cache, scale_factor);

  auto & r = aka::as_type<SparseSolverVectorPETSc>(this->getResidual());

  // VecGetSize(r, &sz);
  // std::cout << "AAAAAA: r " << sz << std::endl;

  VecPointwiseMult(cache, A, cache);
  VecAXPY(r.getVec(), 1., cache.getVec());
}
/* -------------------------------------------------------------------------- */
static bool dof_manager_is_registered =
    DOFManagerFactory::getInstance().registerAllocator(
        "petsc", [](Mesh & mesh, const ID & id) -> std::unique_ptr<DOFManager> {
          return std::make_unique<DOFManagerPETSc>(mesh, id);
        });

} // namespace akantu
