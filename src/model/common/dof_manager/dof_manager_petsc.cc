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
#include "aka_common.hh"
#include "aka_iterators.hh"
#include "communicator.hh"
#include "cppargparse.hh"
#include "mesh.hh"
#include "non_linear_solver_default.hh"
#include "non_linear_solver_petsc.hh"
#include "non_linear_solver_tao.hh"
#include "solver_vector_petsc.hh"
#include "sparse_matrix_petsc.hh"
#include "time_step_solver_default.hh"
#if defined(AKANTU_USE_MPI)
#include "mpi_communicator_data.hh"
#endif
/* -------------------------------------------------------------------------- */
#include <memory>
#include <mpi.h>
#include <numeric>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
/* -------------------------------------------------------------------------- */
#include <petscao.h>
#include <petscerror.h>
#include <petscis.h>
#include <petscistypes.h>
#include <petscsys.h>
#include <petscsystypes.h>
#include <petscvec.h>
/* -------------------------------------------------------------------------- */

namespace akantu {

template <class Func, class... Args>
auto PETSc_call(Func && func, Args... args) -> decltype(auto) {
  auto ierr = std::forward<Func>(func)(std::forward<Args>(args)...);
  if (PetscUnlikely(ierr != 0)) {
    const char * desc{nullptr};
    PetscErrorMessage(ierr, &desc, nullptr);
    AKANTU_EXCEPTION("Error in PETSc call: " << std::string(desc));
  }
  return ierr;
}

inline auto petscErrorHandler(MPI_Comm /* comm */, int line, const char * fun,
                              const char * file, PetscErrorCode n,
                              PetscErrorType p, const char * mess,
                              void * /* ctx */) -> PetscErrorCode {
  if (PetscUnlikely(n != 0)) {
    const char * desc{nullptr};
    PetscErrorMessage(n, &desc, nullptr);
    AKANTU_EXCEPTION(file << ":" << line << ": Error(" << p
                          << ") in PETSc call to \'" << fun << "\': " << mess);
  }
  return n;
}

/* -------------------------------------------------------------------------- */
class PETScSingleton {
private:
  PETScSingleton() {
    // Using PETSc_call since the Error handler is not set yet

    PETSc_call(PetscInitialized, &is_initialized);

    if (is_initialized == PETSC_FALSE) {
      cppargparse::ArgumentParser & argparser = getStaticArgumentParser();
      int & argc = argparser.getArgC();
      char **& argv = argparser.getArgV();
      PETSc_call(PetscInitialize, &argc, &argv, nullptr, nullptr);

      // remove the default PETSc signal handler
      PETSc_call(PetscPopErrorHandler);
      PETSc_call(PetscPushErrorHandler, petscErrorHandler, nullptr);
    }
  }

public:
  PETScSingleton(PETScSingleton &&) = delete;
  auto operator=(PETScSingleton &&) -> PETScSingleton & = delete;
  PETScSingleton(const PETScSingleton &) = delete;
  auto operator=(const PETScSingleton &) -> PETScSingleton & = delete;

  ~PETScSingleton() {
    if (is_initialized == 0U) {
      int mpi_finalized{0};
      MPI_Finalized(&mpi_finalized);
      if (mpi_finalized == 0) {
        PetscFinalize();
      }
    }
  }

  static auto getInstance() -> PETScSingleton & {
    static PETScSingleton instance;
    return instance;
  }

private:
  PetscBool is_initialized{PETSC_FALSE};
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
void DOFManagerPETSc::updateLocalEquationNumber(const ID & dof_id) {
  auto & dof_data = this->getDOFDataTyped<DOFDataPETSc>(dof_id);

  Array<PetscInt> gidx(dof_data.local_equation_number.size());
  for (auto && [local, global] : zip(dof_data.local_equation_number, gidx)) {
    global = localToGlobalEquationNumber(local);
  }

  AOApplicationToPetsc(ao, PetscInt(gidx.size()), gidx.data());

  auto & lidx = dof_data.local_equation_number_petsc;
  if (is_ltog_map != nullptr) {
    lidx.resize(gidx.size());

    PetscInt n;
    ISGlobalToLocalMappingApply(is_ltog_map, IS_GTOLM_MASK, gidx.size(),
                                gidx.data(), &n, lidx.data());
  }
}

/* -------------------------------------------------------------------------- */
auto DOFManagerPETSc::updateNodalDOFs(const ID & dof_id,
                                      const Array<Idx> & nodes_list)
    -> std::pair<Int, Int> {
  auto && ret = DOFManager::updateNodalDOFs(dof_id, nodes_list);
  this->setISLocalToGlobalMapping();
  this->updateLocalEquationNumber(dof_id);
  return ret;
}

/* -------------------------------------------------------------------------- */
void DOFManagerPETSc::setSolverVectorDataForParallelism(
    SolverVectorPETSc & vector) const {
  Vec & x = vector.getVec();

  auto nb_local_dofs = this->getPureLocalSystemSize();
  auto system_size = this->getSystemSize();
  VecSetSizes(x, nb_local_dofs, system_size);

  VecSetLocalToGlobalMapping(x, is_ltog_map);
}

/* -------------------------------------------------------------------------- */
void DOFManagerPETSc::setISLocalToGlobalMapping() {
  auto local_system_size = this->getLocalSystemSize();
  auto nb_local_dofs = this->getPureLocalSystemSize();

  if (ao != nullptr) {
    AODestroy(&ao);
  }

  Array<PetscInt> app_indexes(nb_local_dofs);
  for (auto && [lidx, app_idx] :
       zip(filter_if(arange(local_system_size),
                     [&](auto lidx) {
                       auto is_local = this->isLocalOrMasterDOF(lidx);
                       return is_local;
                     }),
           app_indexes)) {
    app_idx = this->localToGlobalEquationNumber(lidx);
  }

  AOCreateBasic(mpi_communicator, nb_local_dofs, app_indexes.data(), nullptr,
                &ao);

  Vec x{nullptr};
  VecCreate(this->getMPIComm(), &x);
  VecSetFromOptions(x);
  VecSetSizes(x, nb_local_dofs, PETSC_DECIDE);

  VecType vec_type{};
  VecGetType(x, &vec_type);
  if (std::string(vec_type) == std::string(VECMPI)) {
    std::vector<Int> app_ghosts;
    for (auto lidx : filter_if(arange(local_system_size), [&](auto lidx) {
           return this->isSlaveDOF(lidx);
         })) {
      app_ghosts.push_back(this->localToGlobalEquationNumber(lidx));
    }

    auto nghosts = PetscInt(app_ghosts.size());
    ghost_idx.resize(nghosts);
    AOApplicationToPetsc(ao, nghosts, app_ghosts.data());

    VecMPISetGhost(x, nghosts, ghost_idx.data());
  } else {
    std::vector<int> idx(nb_local_dofs);
    std::iota(idx.begin(), idx.end(), 0);

    ISLocalToGlobalMapping is{};
    ISLocalToGlobalMappingCreate(PETSC_COMM_SELF, 1, PetscInt(idx.size()),
                                 idx.data(), PETSC_COPY_VALUES, &is);
    VecSetLocalToGlobalMapping(x, is);
    ISLocalToGlobalMappingDestroy(&is);
  }

  VecGetLocalToGlobalMapping(x, &is_ltog_map);
}

/* -------------------------------------------------------------------------- */
auto DOFManagerPETSc::registerDOFsInternal(const ID & dof_id,
                                           Array<Real> & dofs_array)
    -> std::tuple<Int, Int, Int> {
  dofs_ids.push_back(dof_id);

  auto && [nb_dofs, nb_pure_local_dofs, nb_total_pure_local_dofs] =
      DOFManager::registerDOFsInternal(dof_id, dofs_array);

  this->setISLocalToGlobalMapping();

  // redoing the indexes based on the petsc numbering
  for (auto & dof_id : this->dofs_ids) {
    this->updateLocalEquationNumber(dof_id);
  }

  solution = std::make_unique<SolverVectorPETSc>(*this, id + ":solution");
  residual = std::make_unique<SolverVectorPETSc>(*this, id + ":residual");
  data_cache = std::make_unique<SolverVectorPETSc>(*this, id + ":data_cache");

  for (auto & mat : matrices) {
    auto & A = this->getMatrix(mat.first);
    A.resize();
  }

  return {nb_dofs, nb_pure_local_dofs, nb_total_pure_local_dofs};
}

/* -------------------------------------------------------------------------- */
void DOFManagerPETSc::assembleToGlobalArray(
    const ID & dof_id, const Array<Real> & array_to_assemble,
    SolverVector & global_array, Real scale_factor) {
  const auto & dof_data = getDOFDataTyped<DOFDataPETSc>(dof_id);
  auto & g = aka::as_type<SolverVectorPETSc>(global_array);

  AKANTU_DEBUG_ASSERT(dof_data.local_equation_number.size() ==
                          array_to_assemble.size() *
                              array_to_assemble.getNbComponent(),
                      "The array to assemble does not have a correct size."
                          << " (" << array_to_assemble.getID() << ")");

  AKANTU_DEBUG_ASSERT(dof_data.local_equation_number_petsc.size() ==
                          array_to_assemble.size() *
                              array_to_assemble.getNbComponent(),
                      "The array to assemble does not have a correct size."
                          << " (" << array_to_assemble.getID() << ")");

  g.addValuesLocal(dof_data.local_equation_number_petsc, array_to_assemble,
                   scale_factor);
}

/* -------------------------------------------------------------------------- */
void DOFManagerPETSc::getArrayPerDOFs(const ID & dof_id,
                                      const SolverVector & global_array,
                                      Array<Real> & local) {
  const auto & dof_data = getDOFDataTyped<DOFDataPETSc>(dof_id);
  const auto & petsc_vector = aka::as_type<SolverVectorPETSc>(global_array);

  AKANTU_DEBUG_ASSERT(
      local.size() * local.getNbComponent() ==
          dof_data.local_equation_number_petsc.size(),
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
  DOFManager::assembleMatMulVectToArray_<SolverVectorPETSc>(
      dof_id, A_id, x, array, scale_factor);
}

/* -------------------------------------------------------------------------- */
void DOFManagerPETSc::makeConsistentForPeriodicity(const ID & /*dof_id*/,
                                                   SolverVector & /*array*/) {}

/* -------------------------------------------------------------------------- */
NonLinearSolver & DOFManagerPETSc::getNewNonLinearSolver(
    const ID & id, const ModelSolverOptions & solver_options) {
  switch (solver_options.non_linear_solver_type) {
  case NonLinearSolverType::_petsc_snes:
    return this->registerNonLinearSolver<NonLinearSolverPETSc>(*this, id,
                                                               solver_options);
  case NonLinearSolverType::_petsc_tao:
    return this->registerNonLinearSolver<NonLinearSolverTAO>(*this, id,
                                                             solver_options);
  case NonLinearSolverType::_newton_raphson:
    /* FALLTHRU */
    /* [[fallthrough]]; un-comment when compiler will get it */
  case NonLinearSolverType::_newton_raphson_contact:
  case NonLinearSolverType::_newton_raphson_modified: {
    return this->registerNonLinearSolver<NonLinearSolverNewtonRaphson>(
        *this, id, solver_options);
  }
  case NonLinearSolverType::_linear: {
    return this->registerNonLinearSolver<NonLinearSolverLinear>(*this, id,
                                                                solver_options);
  }
  case NonLinearSolverType::_lumped: {
    return this->registerNonLinearSolver<NonLinearSolverLumped>(*this, id,
                                                                solver_options);
  }
  default:
    return this->registerNonLinearSolver<NonLinearSolverPETSc>(*this, id,
                                                               solver_options);
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
SolverVector & DOFManagerPETSc::getNewLumpedMatrix(const ID & id) {
  return this->registerLumpedMatrix<SolverVectorPETSc>(*this, id);
}

/* -------------------------------------------------------------------------- */
SolverVectorPETSc & DOFManagerPETSc::_getSolution() {
  return aka::as_type<SolverVectorPETSc>(*this->solution);
}

const SolverVectorPETSc & DOFManagerPETSc::_getSolution() const {
  return aka::as_type<SolverVectorPETSc>(*this->solution);
}

SolverVectorPETSc & DOFManagerPETSc::_getResidual() {
  return aka::as_type<SolverVectorPETSc>(*this->residual);
}

const SolverVectorPETSc & DOFManagerPETSc::_getResidual() const {
  return aka::as_type<SolverVectorPETSc>(*this->residual);
}

/* -------------------------------------------------------------------------- */
void DOFManagerPETSc::assembleLumpedMatMulVectToResidual(const ID & dof_id,
                                                         const ID & A_id,
                                                         const Array<Real> & x,
                                                         Real scale_factor) {
  const auto & A = aka::as_type<SolverVectorPETSc>(this->getLumpedMatrix(A_id));
  auto & cache = aka::as_type<SolverVectorPETSc>(*this->data_cache);

  // int sz;
  // VecGetSize(A, &sz);
  // std::cout << "AAAAAA: A " << sz << std::endl;

  // VecGetSize(cache, &sz);
  // std::cout << "AAAAAA: cache " << sz << std::endl;

  cache.zero();
  this->assembleToGlobalArray(dof_id, x, cache, scale_factor);

  auto & r = aka::as_type<SolverVectorPETSc>(this->getResidual());

  // VecGetSize(r, &sz);
  // std::cout << "AAAAAA: r " << sz << std::endl;

  // VecView(cache, PETSC_VIEWER_STDOUT_WORLD);
  // VecView(A, PETSC_VIEWER_STDOUT_WORLD);

  VecPointwiseMult(cache, A, cache);
  VecAXPY(r, 1., cache);
}
/* -------------------------------------------------------------------------- */
static bool dof_manager_is_registered =
    DOFManagerFactory::getInstance().registerAllocator(
        "petsc", [](Mesh & mesh, const ID & id) -> std::unique_ptr<DOFManager> {
          return std::make_unique<DOFManagerPETSc>(mesh, id);
        });

} // namespace akantu
