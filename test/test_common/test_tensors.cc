/**
 * Copyright (©) 2017-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_array.hh"
#include "aka_iterators.hh"
#include "aka_types.hh"
/* -------------------------------------------------------------------------- */
#include <cstdlib>
#include <gtest/gtest.h>
/* -------------------------------------------------------------------------- */

using namespace akantu;

namespace {

/* -------------------------------------------------------------------------- */
class TensorConstructorFixture : public ::testing::Test {
public:
  void SetUp() override {
    for (auto & r : reference) {
      r = rand(); // google-test seeds srand()
    }
  }
  void TearDown() override {}

  template <typename V> void compareToRef(const V & v) {
    for (int i = 0; i < size_; ++i) {
      EXPECT_DOUBLE_EQ(reference[i], v.data()[i]);
    }
  }

protected:
  const int size_{24};
  const std::array<int, 2> mat_size{{4, 6}};
  // const std::array<int, 3> tens3_size{{4, 2, 3}};
  std::array<double, 24> reference;
};

/* -------------------------------------------------------------------------- */
class TensorFixture : public TensorConstructorFixture {
public:
  TensorFixture()
      : vref(reference.data(), size_),
        mref(reference.data(), mat_size[0], mat_size[1]) {}

protected:
  VectorProxy<double> vref;
  MatrixProxy<double> mref;
};

/* -------------------------------------------------------------------------- */
// Vector ----------------------------------------------------------------------
TEST_F(TensorConstructorFixture, VectorDefaultConstruct) {
  Vector<double> v;
  EXPECT_EQ(0, v.size());
  EXPECT_EQ(nullptr, v.data());
}

TEST_F(TensorConstructorFixture, VectorConstruct1) {
  Vector<double> v(size_);
  EXPECT_EQ(size_, v.size());
}

TEST_F(TensorConstructorFixture, VectorConstructWrapped) {
  VectorProxy<double> v(reference.data(), size_);
  EXPECT_EQ(size_, v.size());
  EXPECT_EQ(v.data(), reference.data());

  for (int i = 0; i < size_; ++i) {
    EXPECT_DOUBLE_EQ(reference[i], v(i));
    EXPECT_DOUBLE_EQ(reference[i], v[i]);
  }
}

TEST_F(TensorConstructorFixture, VectorConstructInitializer) {
  Vector<double> v{0., 1., 2., 3., 4., 5.};
  EXPECT_EQ(6, v.size());
  for (int i = 0; i < 6; ++i) {
    EXPECT_DOUBLE_EQ(i, v(i));
  }
}

TEST_F(TensorConstructorFixture, VectorConstructCopy1) {
  VectorProxy<double> vref(reference.data(), reference.size());
  Vector<double> v(vref);
  EXPECT_EQ(size_, v.size());
  compareToRef(v);
}

TEST_F(TensorConstructorFixture, VectorConstructProxy1) {
  VectorProxy<double> vref(reference.data(), reference.size());
  EXPECT_EQ(size_, vref.size());
  compareToRef(vref);

  Vector<double> v(vref);
  EXPECT_EQ(size_, v.size());
  compareToRef(v);
}

/* -------------------------------------------------------------------------- */
TEST_F(TensorFixture, VectorEqual) {
  Vector<double> v;

  v = vref;
  compareToRef(v);

  EXPECT_EQ(size_, v.size());
}

TEST_F(TensorFixture, VectorEqualProxy) {
  VectorProxy<double> vref_proxy(vref);
  Vector<double> v;

  v = vref;
  compareToRef(v);

  EXPECT_EQ(size_, v.size());
}

TEST_F(TensorFixture, VectorEqualProxy2) {
  Vector<double> v_store(size_);
  v_store.zero();
  VectorProxy<double> v(v_store.data(), size_);

  v = vref;
  compareToRef(v);
  compareToRef(v_store);
}

/* -------------------------------------------------------------------------- */
TEST_F(TensorFixture, VectorSet) {
  Vector<double> v(vref);
  compareToRef(v);

  double r = rand();
  v.set(r);

  for (int i = 0; i < size_; ++i)
    EXPECT_DOUBLE_EQ(r, v[i]);
}

TEST_F(TensorFixture, VectorClear) {
  Vector<double> v(vref);
  compareToRef(v);

  v.zero();

  for (int i = 0; i < size_; ++i)
    EXPECT_DOUBLE_EQ(0, v[i]);
}

/* -------------------------------------------------------------------------- */
TEST_F(TensorFixture, VectorDivide) {
  Vector<double> v;
  double r = rand();
  v = vref / r;

  for (int i = 0; i < size_; ++i) {
    EXPECT_DOUBLE_EQ(reference[i] / r, v[i]);
  }
}

TEST_F(TensorFixture, VectorMultiply1) {
  Vector<double> v;
  double r = rand();
  v = vref * r;

  for (int i = 0; i < size_; ++i) {
    EXPECT_DOUBLE_EQ(reference[i] * r, v[i]);
  }
}

TEST_F(TensorFixture, VectorMultiply2) {
  Vector<double> v;
  double r = rand();
  v = r * vref;

  for (int i = 0; i < size_; ++i) {
    EXPECT_DOUBLE_EQ(reference[i] * r, v[i]);
  }
}

TEST_F(TensorFixture, VectorAddition) {
  Vector<double> v;
  v = vref + vref;

  for (int i = 0; i < size_; ++i) {
    EXPECT_DOUBLE_EQ(reference[i] * 2., v[i]);
  }
}

TEST_F(TensorFixture, VectorSubstract) {
  Vector<double> v;
  v = vref - vref;

  for (int i = 0; i < size_; ++i) {
    EXPECT_DOUBLE_EQ(0., v[i]);
  }
}

TEST_F(TensorFixture, VectorDivideEqual) {
  Vector<double> v(vref);
  double r = rand();
  v /= r;

  for (int i = 0; i < size_; ++i) {
    EXPECT_DOUBLE_EQ(reference[i] / r, v[i]);
  }
}

TEST_F(TensorFixture, VectorMultiplyEqual1) {
  Vector<double> v(vref);
  double r = rand();
  v *= r;

  for (int i = 0; i < size_; ++i) {
    EXPECT_DOUBLE_EQ(reference[i] * r, v[i]);
  }
}

TEST_F(TensorFixture, VectorMultiplyEqual2) {
  Vector<double> v(vref);
  v.array() *= v.array();

  for (int i = 0; i < size_; ++i) {
    EXPECT_DOUBLE_EQ(reference[i] * reference[i], v[i]);
  }
}

TEST_F(TensorFixture, VectorAdditionEqual) {
  Vector<double> v(vref);
  v += vref;

  for (int i = 0; i < size_; ++i) {
    EXPECT_DOUBLE_EQ(reference[i] * 2., v[i]);
  }
}

TEST_F(TensorFixture, VectorSubstractEqual) {
  Vector<double> v(vref);
  v -= vref;

  for (int i = 0; i < size_; ++i) {
    EXPECT_DOUBLE_EQ(0., v[i]);
  }
}

/* -------------------------------------------------------------------------- */
// Matrix ----------------------------------------------------------------------

TEST_F(TensorConstructorFixture, MatrixDefaultConstruct) {
  Matrix<double> m;
  EXPECT_EQ(0, m.size());
  EXPECT_EQ(0, m.rows());
  EXPECT_EQ(0, m.cols());
  EXPECT_EQ(nullptr, m.data());
}

TEST_F(TensorConstructorFixture, MatrixConstruct1) {
  Matrix<double> m(mat_size[0], mat_size[1]);
  EXPECT_EQ(size_, m.size());
  EXPECT_EQ(mat_size[0], m.rows());
  EXPECT_EQ(mat_size[1], m.cols());
}

TEST_F(TensorConstructorFixture, MatrixConstructWrapped) {
  MatrixProxy<double> m(reference.data(), mat_size[0], mat_size[1]);
  EXPECT_EQ(size_, m.size());
  EXPECT_EQ(mat_size[0], m.rows());
  EXPECT_EQ(mat_size[1], m.cols());

  for (int i = 0; i < mat_size[0]; ++i) {
    for (int j = 0; j < mat_size[1]; ++j) {
      EXPECT_DOUBLE_EQ(reference[i + j * mat_size[0]], m(i, j));
    }
  }
  compareToRef(m);
}

TEST_F(TensorConstructorFixture, MatrixConstructInitializer) {
  Matrix<double> m{{0., 1., 2.}, {3., 4., 5.}};
  EXPECT_EQ(6, m.size());
  EXPECT_EQ(2, m.rows());
  EXPECT_EQ(3, m.cols());

  int c = 0;
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 3; ++j, ++c) {
      EXPECT_DOUBLE_EQ(c, m(i, j));
    }
  }
}

TEST_F(TensorConstructorFixture, MatrixConstructCopy1) {
  MatrixProxy<double> mref(reference.data(), mat_size[0], mat_size[1]);
  Matrix<double> m(mref);
  EXPECT_EQ(size_, m.size());
  EXPECT_EQ(mat_size[0], m.rows());
  EXPECT_EQ(mat_size[1], m.cols());
  compareToRef(m);
}

TEST_F(TensorConstructorFixture, MatrixConstructProxy1) {
  MatrixProxy<double> mref(reference.data(), mat_size[0], mat_size[1]);
  EXPECT_EQ(size_, mref.size());
  EXPECT_EQ(mat_size[0], mref.size(0));
  EXPECT_EQ(mat_size[1], mref.size(1));
  compareToRef(mref);

  Matrix<double> m(mref);
  EXPECT_EQ(size_, m.size());
  EXPECT_EQ(mat_size[0], m.rows());
  EXPECT_EQ(mat_size[1], m.cols());
  compareToRef(m);
}

/* -------------------------------------------------------------------------- */
TEST_F(TensorFixture, MatrixEqual) {
  Matrix<double> m;

  m = mref;
  compareToRef(m);

  EXPECT_EQ(size_, m.size());
  EXPECT_EQ(mat_size[0], m.rows());
  EXPECT_EQ(mat_size[1], m.cols());
}

TEST_F(TensorFixture, MatrixEqualProxy1) {
  MatrixProxy<double> mref_proxy(mref.data(), mref.rows(), mref.cols());
  Matrix<double> m;

  m = mref;
  compareToRef(m);

  EXPECT_EQ(size_, m.size());
  EXPECT_EQ(mat_size[0], m.rows());
  EXPECT_EQ(mat_size[1], m.cols());
}

TEST_F(TensorFixture, MatrixEqualProxy2) {
  Matrix<double> m_store(mat_size[0], mat_size[1]);
  m_store.zero();
  MatrixProxy<double> m(m_store.data(), mat_size[0], mat_size[1]);

  m = mref;
  compareToRef(m);
  compareToRef(m_store);
}

TEST_F(TensorFixture, MatrixEqualSlice) {
  Matrix<double> m(mat_size[0], mat_size[1]);
  m.zero();

  for (Int i = 0; i < m.cols(); ++i) {
    m(i) = mref(i);
  }

  compareToRef(m);
}

/* -------------------------------------------------------------------------- */
TEST_F(TensorFixture, MatrixSet) {
  Matrix<double> m(mref);
  compareToRef(m);

  double r = rand();
  m.set(r);

  for (int i = 0; i < size_; ++i) {
    EXPECT_DOUBLE_EQ(r, m.array()(i));
  }
}

TEST_F(TensorFixture, MatrixClear) {
  Matrix<double> m(mref);
  compareToRef(m);

  m.zero();

  for (int i = 0; i < size_; ++i) {
    EXPECT_DOUBLE_EQ(0, m.array()(i));
  }
}

/* -------------------------------------------------------------------------- */
TEST_F(TensorFixture, MatrixDivide) {
  Matrix<double> m;
  double r = rand();
  m = mref / r;

  for (int i = 0; i < size_; ++i) {
    EXPECT_DOUBLE_EQ(reference[i] / r, m.array()(i));
  }
}

TEST_F(TensorFixture, MatrixMultiply1) {
  Matrix<double> m;
  double r = rand();
  m = mref * r;

  for (int i = 0; i < size_; ++i) {
    EXPECT_DOUBLE_EQ(reference[i] * r, m.array()(i));
  }
}

TEST_F(TensorFixture, MatrixMultiply2) {
  Matrix<double> m;
  double r = rand();
  m = r * mref;

  for (int i = 0; i < size_; ++i) {
    EXPECT_DOUBLE_EQ(reference[i] * r, m.array()(i));
  }
}

TEST_F(TensorFixture, MatrixAddition) {
  Matrix<double> m;
  m = mref + mref;

  for (int i = 0; i < size_; ++i) {
    EXPECT_DOUBLE_EQ(reference[i] * 2., m.array()(i));
  }
}

TEST_F(TensorFixture, MatrixSubstract) {
  Matrix<double> m;
  m = mref - mref;

  for (int i = 0; i < size_; ++i) {
    EXPECT_DOUBLE_EQ(0., m.array()(i));
  }
}

TEST_F(TensorFixture, MatrixDivideEqual) {
  Matrix<double> m(mref);
  double r = rand();
  m /= r;

  for (int i = 0; i < size_; ++i) {
    EXPECT_DOUBLE_EQ(reference[i] / r, m.array()(i));
  }
}

TEST_F(TensorFixture, MatrixMultiplyEqual1) {
  Matrix<double> m(mref);
  double r = rand();
  m *= r;

  for (int i = 0; i < size_; ++i) {
    EXPECT_DOUBLE_EQ(reference[i] * r, m.array()(i));
  }
}

TEST_F(TensorFixture, MatrixAdditionEqual) {
  Matrix<double> m(mref);
  m += mref;

  for (int i = 0; i < size_; ++i) {
    EXPECT_DOUBLE_EQ(reference[i] * 2., m.array()(i));
  }
}

TEST_F(TensorFixture, MatrixSubstractEqual) {
  Matrix<double> m(mref);
  m -= mref;

  for (int i = 0; i < size_; ++i) {
    EXPECT_DOUBLE_EQ(0., m.array()(i));
  }
}

TEST_F(TensorFixture, MatrixIterator) {
  Matrix<double> m(mref);

  UInt col_count = 0;
  for (auto && col : m) {
    VectorProxy<Real> col_hand(m.data() + col_count * m.rows(), m.rows());
    Vector<Real> col_wrap(col);

    auto comp = (col_wrap - col_hand).lpNorm<Eigen::Infinity>();
    EXPECT_DOUBLE_EQ(0., comp);
    ++col_count;
  }
}

TEST_F(TensorFixture, MatrixIteratorZip) {
  Matrix<double> m1(mref);
  Matrix<double> m2(mref);

  for (auto && [col1, col2] : zip(m1, m2)) {
    auto comp = (col1 - col2).lpNorm<Eigen::Infinity>();
    EXPECT_DOUBLE_EQ(0., comp);
  }
}

TEST_F(TensorFixture, VectorViseIteratorZip) {
  Matrix<double> m1(mref);
  Matrix<double> m2(mref);

  for (auto && [col1, col2] : zip(m1.colwise(), m2.colwise())) {
    auto comp = (col1 - col2).lpNorm<Eigen::Infinity>();
    EXPECT_DOUBLE_EQ(0., comp);
  }
}

TEST_F(TensorFixture, VectorIteratorZip) {
  Matrix<double> m1(mref);
  Matrix<double> m2(mref);

  for (auto && [col1, col2] : zip(m1.col(0), m2.col(0))) {
    auto comp = (col1 - col2);
    EXPECT_DOUBLE_EQ(0., comp);
  }
}

TEST_F(TensorFixture, MatrixEigs) {
  Matrix<double, 4, 4> A{
      {0, 1., 0, 0}, {1., 0, 0, 0}, {0, 1., 0, 1.}, {0, 0, 4., 0}};

  Matrix<double, 4, 4> v;
  Vector<double, 4> lambda;
  lambda.zero();
  v.zero();
  A.eig(lambda, v);

  Vector<double> eigs_ref{2, 1., -1., -2};

  auto Av = (A * v).eval();

  for (int i = 0; i < 4; ++i) {
    EXPECT_NEAR(eigs_ref(i), lambda(i), 1e-14);
  }

  for (int i = 0; i < 4; ++i) {
    auto lambda_v_minus_a_v =
        (lambda(i) * v(i) - Av(i)).template lpNorm<Eigen::Infinity>();

    EXPECT_NEAR(lambda_v_minus_a_v, 0., 1e-14);
  }
}

/* -------------------------------------------------------------------------- */

} // namespace
