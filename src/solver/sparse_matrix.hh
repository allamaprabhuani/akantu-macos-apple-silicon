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
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_SPARSE_MATRIX_HH_
#define AKANTU_SPARSE_MATRIX_HH_

/* -------------------------------------------------------------------------- */
namespace akantu {
class DOFManager;
class TermsToAssemble;
class SparseSolverVector;
} // namespace akantu

namespace akantu {

class SparseMatrix {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  SparseMatrix(DOFManager & dof_manager, const MatrixType & matrix_type,
               const ID & id = "sparse_matrix");

  SparseMatrix(const SparseMatrix & matrix, const ID & id = "sparse_matrix");

  virtual ~SparseMatrix();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// remove the existing profile
  virtual void clearProfile();

  /// set the matrix to 0
  virtual void set(Real val) = 0;

  virtual void zero() { this->set(0); }

  /// add a non-zero element to the profile
  virtual Idx add(Idx i, Idx j) = 0;

  /// assemble a local matrix in the sparse one
  virtual void add(Idx i, Idx j, Real value) = 0;

  /// save the profil in a file using the MatrixMarket file format
  virtual void saveProfile(const std::string & /* filename */) const {
    AKANTU_TO_IMPLEMENT();
  }

  /// save the matrix in a file using the MatrixMarket file format
  virtual void saveMatrix(const std::string & /* filename */) const {
    AKANTU_TO_IMPLEMENT();
  };

  /// multiply the matrix by a coefficient
  virtual void mul(Real alpha) = 0;

  /// add matrices
  virtual void add(const SparseMatrix & B, Real alpha = 1.);

  /// Equivalent of *gemv in blas
  virtual void matVecMul(const SparseSolverVector & x, SparseSolverVector & y,
                         Real alpha = 1., Real beta = 0.) const = 0;

  /// modify the matrix to "remove" the blocked dof
  virtual void applyBoundary(Real block_val = 1.) = 0;

  /// copy the profile of another matrix
  virtual void copyProfile(const SparseMatrix & other) = 0;

  /// operator *=
  SparseMatrix & operator*=(Real alpha) {
    this->mul(alpha);
    return *this;
  }

  /// Check if all entries are finite. The default implementation throws.
  virtual bool isFinite() const { AKANTU_TO_IMPLEMENT(); }

protected:
  /// This is the revert of add \f[B += \alpha * *this\f];
  virtual void addMeTo(SparseMatrix & B, Real alpha) const = 0;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// return the values at potition i, j
  virtual inline Real operator()(Idx /*i*/, Idx /*j*/) const {
    AKANTU_TO_IMPLEMENT();
  }
  /// return the values at potition i, j
  virtual inline Real & operator()(Idx /*i*/, Idx /*j*/) {
    AKANTU_TO_IMPLEMENT();
  }

  AKANTU_GET_MACRO_AUTO(NbNonZero, nb_non_zero);
  Int size() const { return size_; }
  AKANTU_GET_MACRO_AUTO(MatrixType, matrix_type);

  virtual Int getRelease() const = 0;

  virtual Real min() = 0;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  ID id;

  /// Underlying dof manager
  DOFManager & dof_manager;

  /// sparce matrix type
  MatrixType matrix_type;

  /// Size of the matrix
  Int size_;

  /// number of processors
  Int nb_proc;

  /// number of non zero element
  Int nb_non_zero;
};

// Array<Real> & operator*=(Array<Real> & vect, const SparseMatrix & mat);

} // namespace akantu

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
#include "sparse_matrix_inline_impl.hh"

#endif /* AKANTU_SPARSE_MATRIX_HH_ */
