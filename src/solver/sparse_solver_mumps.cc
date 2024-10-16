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
#include "dof_manager_default.hh"
#include "dof_synchronizer.hh"
#include "solver_vector_default.hh"
#include "sparse_matrix_aij.hh"
#if defined(AKANTU_USE_MPI)
#include "mpi_communicator_data.hh"
#endif
#include "sparse_solver_mumps.hh"
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
// static std::ostream & operator <<(std::ostream & stream, const DMUMPS_STRUC_C
// & _this) {
//   stream << "DMUMPS Data [" << std::endl;
//   stream << " + job          : " << _this.job          << std::endl;
//   stream << " + par          : " << _this.par          << std::endl;
//   stream << " + sym          : " << _this.sym          << std::endl;
//   stream << " + comm_fortran : " << _this.comm_fortran << std::endl;
//   stream << " + nz           : " << _this.nz           << std::endl;
//   stream << " + irn          : " << _this.irn          << std::endl;
//   stream << " + jcn          : " << _this.jcn          << std::endl;
//   stream << " + nz_loc       : " << _this.nz_loc       << std::endl;
//   stream << " + irn_loc      : " << _this.irn_loc      << std::endl;
//   stream << " + jcn_loc      : " << _this.jcn_loc      << std::endl;
//   stream << "]";
//   return stream;
// }

namespace akantu {

/* -------------------------------------------------------------------------- */
SparseSolverMumps::SparseSolverMumps(DOFManagerDefault & dof_manager,
                                     const ID & matrix_id, const ID & id)
    : SparseSolver(dof_manager, matrix_id, id), master_rhs_solution(0, 1) {
  AKANTU_DEBUG_IN();

  this->prank = communicator.whoAmI();

#ifdef AKANTU_USE_MPI
  this->parallel_method = _fully_distributed;
#else  // AKANTU_USE_MPI
  this->parallel_method = _not_parallel;
#endif // AKANTU_USE_MPI

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
SparseSolverMumps::~SparseSolverMumps() {
  AKANTU_DEBUG_IN();

  mumpsDataDestroy();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SparseSolverMumps::mumpsDataDestroy() {
#ifdef AKANTU_USE_MPI
  int finalized = 0;
  MPI_Finalized(&finalized);
  if (finalized != 0) { // Da fuck !?
    return;
  }
#endif

  if (this->is_initialized) {
    this->mumps_data.job = _smj_destroy; // destroy
    dmumps_c(&this->mumps_data);
    this->is_initialized = false;
  }
}

/* -------------------------------------------------------------------------- */
void SparseSolverMumps::destroyInternalData() { mumpsDataDestroy(); }

/* -------------------------------------------------------------------------- */
void SparseSolverMumps::checkInitialized() {
  if (this->is_initialized) {
    return;
  }

  this->initialize();
}

/* -------------------------------------------------------------------------- */
void SparseSolverMumps::setOutputLevel() {
  // Output setup
  icntl(1) = 0; // error output
  icntl(2) = 0; // diagnostics output
  icntl(3) = 0; // information
  icntl(4) = 0;

#if !defined(AKANTU_NDEBUG)
  DebugLevel dbg_lvl = debug::debugger.getDebugLevel();

  if (AKANTU_DEBUG_TEST(dblDump)) {
    strcpy(this->mumps_data.write_problem, "mumps_matrix.mtx");
  }

  // clang-format off
  icntl(1) = (dbg_lvl >= dblWarning) ? 6 : 0;
  icntl(3) = (dbg_lvl >= dblInfo)    ? 6 : 0;
  icntl(2) = (dbg_lvl >= dblTrace)   ? 6 : 0;

  icntl(4) =
    dbg_lvl >= dblDump    ? 4 :
    dbg_lvl >= dblTrace   ? 3 :
    dbg_lvl >= dblInfo    ? 2 :
    dbg_lvl >= dblWarning ? 1 :
                            0;

// clang-format on
#endif
}

/* -------------------------------------------------------------------------- */
void SparseSolverMumps::initMumpsData() {
  auto & A = this->getDOFManager().getMatrix(matrix_id);

  // Default Scaling
  icntl(8) = 77;

  // Assembled matrix
  icntl(5) = 0;

  /// Default centralized dense second member
  icntl(20) = 0;
  icntl(21) = 0;

  // automatic choice for analysis
  icntl(28) = 0;

  auto size = A.size();

  if (prank == 0) {
    this->master_rhs_solution.resize(size);
  }

  this->mumps_data.nz_alloc = 0;
  this->mumps_data.n = size;

  switch (this->parallel_method) {
  case _fully_distributed:
    icntl(18) = 3; // fully distributed

    this->mumps_data.nz_loc = A.getNbNonZero();
    this->mumps_data.irn_loc = A.irn.data();
    this->mumps_data.jcn_loc = A.jcn.data();

    break;
  case _not_parallel:
  case _master_slave_distributed:
    icntl(18) = 0; // centralized

    if (prank == 0) {
      this->mumps_data.nz = A.getNbNonZero();
      this->mumps_data.irn = A.irn.data();
      this->mumps_data.jcn = A.jcn.data();
    } else {
      this->mumps_data.nz = 0;
      this->mumps_data.irn = nullptr;
      this->mumps_data.jcn = nullptr;
    }
    break;
  default:
    AKANTU_ERROR("This case should not happen!!");
  }
}

/* -------------------------------------------------------------------------- */
void SparseSolverMumps::initialize() {
  AKANTU_DEBUG_IN();

  this->mumps_data.par = 1; // The host is part of computations

  switch (this->parallel_method) {
  case _not_parallel:
    break;
  case _master_slave_distributed:
    this->mumps_data.par = 0; // The host is not part of the computations
    /* FALLTHRU */
    /* [[fallthrough]]; un-comment when compiler will get it */
  case _fully_distributed:
#ifdef AKANTU_USE_MPI
    const auto & mpi_data =
        aka::as_type<MPICommunicatorData>(communicator.getCommunicatorData());
    MPI_Comm mpi_comm = mpi_data.getMPICommunicator();
    this->mumps_data.comm_fortran = MPI_Comm_c2f(mpi_comm);
#else
    AKANTU_ERROR(
        "You cannot use parallel method to solve without activating MPI");
#endif
    break;
  }

  const auto & A = dof_manager.getMatrix(matrix_id);

  this->mumps_data.sym = 2 * static_cast<int>(A.getMatrixType() == _symmetric);
  this->prank = communicator.whoAmI();

  this->setOutputLevel();

  this->mumps_data.job = _smj_initialize; // initialize
  dmumps_c(&this->mumps_data);

  this->setOutputLevel();

  this->is_initialized = true;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SparseSolverMumps::analysis() {
  AKANTU_DEBUG_IN();

  initMumpsData();

  this->mumps_data.job = _smj_analyze; // analyze
  dmumps_c(&this->mumps_data);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SparseSolverMumps::factorize() {
  AKANTU_DEBUG_IN();

  auto & A = this->getDOFManager().getMatrix(matrix_id);

  if (parallel_method == _fully_distributed) {
    this->mumps_data.a_loc = A.a.data();
  } else {
    if (prank == 0) {
      this->mumps_data.a = A.a.data();
    }
  }

  this->mumps_data.job = _smj_factorize; // factorize
  dmumps_c(&this->mumps_data);

  this->printError();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SparseSolverMumps::solve() {
  this->master_rhs_solution.copy(
      aka::as_type<SolverVectorDefault>(dof_manager.getResidual())
          .getGlobalVector());

  // this->dof_manager.getMatrix(matrix_id).saveMatrix("toto.mtx");
  this->solveInternal();

  aka::as_type<SolverVectorDefault>(dof_manager.getSolution())
      .setGlobalVector(this->master_rhs_solution);

  this->dof_manager.splitSolutionPerDOFs();
}

/* -------------------------------------------------------------------------- */
void SparseSolverMumps::solveInternal() {
  AKANTU_DEBUG_IN();

  this->checkInitialized();

  const auto & A = this->getDOFManager().getMatrix(matrix_id);

  this->setOutputLevel();

  if (this->last_profile_release != A.getProfileRelease()) {
    this->analysis();
    this->last_profile_release = A.getProfileRelease();
  }

  if (AKANTU_DEBUG_TEST(dblDump)) {
    A.saveMatrix("solver_mumps" + std::to_string(prank) + ".mtx");
  }

  if (this->last_value_release != A.getValueRelease()) {
    this->factorize();
    this->last_value_release = A.getValueRelease();
  }

  if (prank == 0) {
    this->mumps_data.rhs = this->master_rhs_solution.data();
  }

  this->mumps_data.job = _smj_solve; // solve
  dmumps_c(&this->mumps_data);

  this->printError();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SparseSolverMumps::printError() {
  Vector<Int> _info_v(2);
  _info_v[0] = info(1);  // to get errors
  _info_v[1] = -info(1); // to get warnings
  dof_manager.getCommunicator().allReduce(_info_v, SynchronizerOperation::_min);
  _info_v[1] = -_info_v[1];

  if (_info_v[0] < 0) { // < 0 is an error
    switch (_info_v[0]) {
    case -10: {
      AKANTU_CUSTOM_EXCEPTION(
          debug::SingularMatrixException(dof_manager.getMatrix(matrix_id)));
      break;
    }
    case -9: {
      icntl(14) += 10;
      if (icntl(14) != 90) {
        // std::cout << "Dynamic memory increase of 10%" << std::endl;
        AKANTU_DEBUG_WARNING("MUMPS dynamic memory is insufficient it will be "
                             "increased allowed to use 10% more");

        // change releases to force a recompute
        this->last_value_release--;
        this->last_profile_release--;

        this->solve();
      } else {
        AKANTU_ERROR("The MUMPS workarray is too small INFO(2)="
                     << info(2) << "No further increase possible");
      }
      break;
    }
    default:
      AKANTU_ERROR("Error in mumps during solve process, check mumps "
                   "user guide INFO(1) = "
                   << _info_v[1]);
    }
  } else if (_info_v[1] > 0) {
    AKANTU_DEBUG_WARNING("Warning in mumps during solve process, check mumps "
                         "user guide INFO(1) = "
                         << _info_v[1]);
  }
}
/* -------------------------------------------------------------------------- */

DOFManagerDefault & SparseSolverMumps::getDOFManager() {
  return aka::as_type<DOFManagerDefault &>(this->dof_manager);
}

} // namespace akantu
