#include "sparse_solver_eigen.hh"
#include "dof_manager_default.hh"
#include "solver_vector_default.hh"
#include "sparse_matrix_aij.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
SparseSolverEigen::SparseSolverEigen(DOFManager & dof_manager,
                                     const ID & matrix_id, const ID & id)
    : SparseSolver(dof_manager, matrix_id, id) {
  AKANTU_DEBUG_ASSERT(communicator.getNbProc() == 1,
                      "This solver does not work in parallel");
}

/* -------------------------------------------------------------------------- */
void SparseSolverEigen::setA() {
  const auto & Aij =
      aka::as_type<SparseMatrixAIJ>(dof_manager.getMatrix(matrix_id));

  if (this->last_value_release == Aij.getRelease() and
      this->last_value_release != -1) {
    return;
  }

  using Triplet = Eigen::Triplet<Real>;
  std::vector<Triplet> triplets;
  triplets.reserve(Aij.getNbNonZero());
  for (auto && [i, j, a] : zip(Aij.getIRN(), Aij.getJCN(), Aij.getA())) {
    triplets.emplace_back(i - 1, j - 1, a);
    if (Aij.getMatrixType() == _symmetric and i != j) {
      triplets.emplace_back(j - 1, i - 1, a);
    }
  }

  A.resize(Aij.size(), Aij.size());
  A.setZero();
  A.setFromTriplets(triplets.begin(), triplets.end());
}

/* -------------------------------------------------------------------------- */
void SparseSolverEigen::solve() {
  const auto & Aij =
      aka::as_type<SparseMatrixAIJ>(dof_manager.getMatrix(matrix_id));

  auto & b_ = aka::as_type<SolverVectorDefault>(this->dof_manager.getResidual())
                  .getGlobalVector();
  Array<Real> x_(b_.size() * b_.getNbComponent());

  VectorProxy<Real> x(x_.data(), x_.size() * x_.getNbComponent());
  VectorProxy<Real> b(b_.data(), b_.size() * b_.getNbComponent());

  setA();

  if (this->last_profile_release != Aij.getProfileRelease()) {
    solver.analyzePattern(A);
    this->last_profile_release = Aij.getProfileRelease();
  }

  if (this->last_value_release != Aij.getValueRelease()) {
    solver.factorize(A);
    this->last_value_release = Aij.getValueRelease();
  }

  x = solver.solve(b);

  aka::as_type<SolverVectorDefault>(this->dof_manager.getSolution())
      .setGlobalVector(x_);

  this->dof_manager.splitSolutionPerDOFs();
}

/* -------------------------------------------------------------------------- */
DOFManagerDefault & SparseSolverEigen::getDOFManager() {
  return aka::as_type<DOFManagerDefault &>(this->dof_manager);
}

} // namespace akantu
