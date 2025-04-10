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
#include "dof_manager.hh"
/* -------------------------------------------------------------------------- */
#include <memory>
#include <mpi.h>
#include <petscao.h>
#include <petscerror.h>
#include <petscis.h>
#include <petscistypes.h>
#include <petscmacros.h>
#include <petscsys.h>
#include <petscsystypes.h>
#include <tuple>
#include <utility>
#include <vector>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_DOF_MANAGER_PETSC_HH
#define AKANTU_DOF_MANAGER_PETSC_HH

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

namespace akantu::detail {
template <typename T> void PETScSetName(T t, const ID & id) {
  // NOLINT(cppcoregui)
  PetscObjectSetName(reinterpret_cast<PetscObject>(t), id.c_str());
}
} // namespace akantu::detail

namespace akantu {
class SparseMatrixPETSc;
class SolverVectorPETSc;
} // namespace akantu

namespace akantu {

class DOFManagerPETSc : public DOFManager {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors */
  /* ------------------------------------------------------------------------ */
public:
  DOFManagerPETSc(const ID & id = "dof_manager_petsc");
  DOFManagerPETSc(Mesh & mesh, const ID & id = "dof_manager_petsc");

  ~DOFManagerPETSc() override = default;

protected:
  void init();

  struct DOFDataPETSc : public DOFData {
    explicit DOFDataPETSc(const ID & dof_id);

    /**
       This is the petsc local numbering which differs from local_numbering from
       parent class. In the petsc linguo parent::local_numbering is in
       "Application Ordering"
     */
    Array<PetscInt> local_equation_number_petsc;

    auto getLocalEquationsNumbers() -> Array<Int> & override {
      return local_equation_number_petsc;
    }
  };

  void updateLocalEquationNumber(const ID & dof_id);

  void setSolverVectorDataForParallelism(SolverVectorPETSc & vector) const;

  friend class SolverVectorPETSc;

  /* ------------------------------------------------------------------------ */
  /* Methods */
  /* ------------------------------------------------------------------------ */
public:
  void assembleElementalMatricesToMatrix(
      const ID & matrix_id, const ID & dof_id,
      const Array<Real> & elementary_mat, ElementType type,
      GhostType ghost_type, const MatrixType & elemental_matrix_type,
      const Array<Idx> & filter_elements) override;

  void assembleMatMulVectToArray(const ID & dof_id, const ID & A_id,
                                 const Array<Real> & x, Array<Real> & array,
                                 Real scale_factor = 1.) override;

  void assembleLumpedMatMulVectToResidual(const ID & dof_id, const ID & A_id,
                                          const Array<Real> & x,
                                          Real scale_factor = 1) override;

  void assemblePreassembledMatrix(const ID & matrix_id,
                                  const TermsToAssemble & /*terms*/) override;

  void assembleToGlobalArray(const ID & dof_id,
                             const Array<Real> & array_to_assemble,
                             SolverVector & global_array,
                             Real scale_factor) override;
  void getArrayPerDOFs(const ID & dof_id, const SolverVector & global,
                       Array<Real> & local) override;

protected:
  void makeConsistentForPeriodicity(const ID & dof_id,
                                    SolverVector & array) override;

  auto getNewDOFData(const ID & dof_id) -> std::unique_ptr<DOFData> override;

  auto registerDOFsInternal(const ID & dof_id, Array<Real> & dofs_array)
      -> std::tuple<Int, Int, Int> override;

  auto updateNodalDOFs(const ID & dof_id, const Array<Idx> & nodes_list)
      -> std::pair<Int, Int> override;

  void setISLocalToGlobalMapping();

  auto getNewNonLinearSolver(const ID & nls_solver_id,
                             const ModelSolverOptions & solver_options)
      -> NonLinearSolver & override;

  auto getNewTimeStepSolver(const ID & id, const TimeStepSolverType & type,
                            NonLinearSolver & non_linear_solver,
                            SolverCallback & solver_callback)
      -> TimeStepSolver & override;

  /* ------------------------------------------------------------------------ */
  /* Accessors */
  /* ------------------------------------------------------------------------ */
public:
  /// Get an instance of a new SparseMatrix
  auto getNewMatrix(const ID & matrix_id, const MatrixType & matrix_type)
      -> SparseMatrix & override;

  /// Get an instance of a new SparseMatrix as a copy of the SparseMatrix
  /// matrix_to_copy_id
  auto getNewMatrix(const ID & matrix_id, const ID & matrix_to_copy_id)
      -> SparseMatrix & override;

  /// Get the reference of an existing matrix
  auto getMatrix(const ID & matrix_id) -> SparseMatrixPETSc &;

  /// Get an instance of a new lumped matrix
  auto getNewLumpedMatrix(const ID & matrix_id) -> SolverVector & override;

  /// Get the blocked dofs array
  AKANTU_GET_MACRO(MPIComm, mpi_communicator, MPI_Comm);

  AKANTU_GET_MACRO_NOT_CONST(ISLocalToGlobalMapping, is_ltog_map,
                             ISLocalToGlobalMapping &);

  SolverVectorPETSc & _getSolution();
  const SolverVectorPETSc & _getSolution() const;

  SolverVectorPETSc & _getResidual();
  const SolverVectorPETSc & _getResidual() const;

  /* ------------------------------------------------------------------------ */
  /* Class Members */
  /* ------------------------------------------------------------------------ */
private:
  using PETScMatrixMap = std::map<ID, SparseMatrixPETSc *>;
  using PETScLumpedMatrixMap = std::map<ID, SolverVectorPETSc *>;

  /// list of matrices registered to the dof manager
  PETScMatrixMap petsc_matrices;

  /// list of lumped matrices registered
  PETScLumpedMatrixMap petsc_lumped_matrices;

  /// PETSc local to global mapping of dofs
  ISLocalToGlobalMapping is_ltog_map{nullptr};

  /// Mapping of akantu global numbering to petsc global numbering
  AO ao{nullptr};

  /// List of ghost dofs in petsc global indexes
  std::vector<PetscInt> ghost_idx;

  /// Communicator associated to PETSc
  MPI_Comm mpi_communicator{MPI_COMM_SELF};

  /// list of the dof ids to be able to always iterate in the same order
  std::vector<ID> dofs_ids;
};

/* -------------------------------------------------------------------------- */

} // namespace akantu

#endif /* AKANTU_DOF_MANAGER_PETSC_HH */
