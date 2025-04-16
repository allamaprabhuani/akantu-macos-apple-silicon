/**
 * Copyright (©) 2011-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_compatibilty_with_cpp_standard.hh"
#include "aka_error.hh"
/* -------------------------------------------------------------------------- */
#include <algorithm>
#include <array>
#include <initializer_list>
#include <numeric>
#include <type_traits>

#ifndef AKANTU_AKA_TYPES_HH
#define AKANTU_AKA_TYPES_HH

/* -------------------------------------------------------------------------- */
#define EIGEN_DEFAULT_DENSE_INDEX_TYPE akantu::Idx
#define EIGEN_DEFAULT_MATRIX_STORAGE_ORDER_OPTION Eigen::ColMajor
#define EIGEN_DEFAULT_IO_FORMAT                                                \
  Eigen::IOFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ",    \
                  "[", "]", "[", "]")
/* -------------------------------------------------------------------------- */
#define EIGEN_MATRIXBASE_PLUGIN "aka_types_eigen_matrix_base_plugin.hh"
#define EIGEN_MATRIX_PLUGIN "aka_types_eigen_matrix_plugin.hh"
#define EIGEN_PLAINOBJECTBASE_PLUGIN                                           \
  "aka_types_eigen_plain_object_base_plugin.hh"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
/* -------------------------------------------------------------------------- */

namespace aka {

template <typename T> struct is_eigen_map : public std::false_type {};

template <typename PlainObjectType, int MapOptions, typename StrideType>
struct is_eigen_map<Eigen::Map<PlainObjectType, MapOptions, StrideType>>
    : public std::true_type {};

/* -------------------------------------------------------------------------- */
template <typename T> struct is_eigen_matrix : public std::false_type {};

template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows,
          int _MaxCols>
struct is_eigen_matrix<
    Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>>
    : public std::true_type {};

/* -------------------------------------------------------------------------- */
template <typename T> struct is_eigen_matrix_base : public std::false_type {};

template <typename Derived>
struct is_eigen_matrix_base<Eigen::MatrixBase<Derived>>
    : public std::true_type {};
} // namespace aka

namespace akantu {

using Eigen::Ref;

template <typename T, Eigen::Index n = Eigen::Dynamic>
using Vector = Eigen::Matrix<T, n, 1>;

template <typename T, Eigen::Index m = Eigen::Dynamic,
          Eigen::Index n = Eigen::Dynamic>
using Matrix = Eigen::Matrix<T, m, n>;

template <typename T, Eigen::Index n = Eigen::Dynamic>
using VectorProxy =
    Eigen::Map<std::conditional_t<std::is_const<T>::value,
                                  const Vector<std::remove_const_t<T>, n>,
                                  Vector<std::remove_const_t<T>, n>>>;

template <typename T, Eigen::Index m = Eigen::Dynamic,
          Eigen::Index n = Eigen::Dynamic>
using MatrixProxy =
    Eigen::Map<std::conditional_t<std::is_const<T>::value,
                                  const Matrix<std::remove_const_t<T>, m, n>,
                                  Matrix<std::remove_const_t<T>, m, n>>>;

using VectorXr = Vector<Real>;
using MatrixXr = Matrix<Real>;

enum NormType : int8_t { L_1 = 1, L_2 = 2, L_inf = -1 };

struct TensorTraitBase {};

template <size_t n> struct TensorTrait : public TensorTraitBase {};

} // namespace akantu

namespace aka {
template <typename Derived>
using is_vector = aka::bool_constant<
    std::remove_reference_t<std::decay_t<Derived>>::IsVectorAtCompileTime>;
template <class V> inline constexpr bool is_vector_v = is_vector<V>::value;

template <typename Derived> using is_matrix = aka::negation<is_vector<Derived>>;
template <class M> inline constexpr bool is_matrix_v = is_matrix<M>::value;

template <typename... Ds>
using are_vectors = aka::conjunction<is_vector<Ds>...>;

template <class... Vs>
inline constexpr bool are_vectors_v = are_vectors<Vs...>::value;

template <typename... Ds>
using are_matrices = aka::conjunction<is_matrix<Ds>...>;

template <class... Ms>
inline constexpr bool are_matrices_v = are_matrices<Ms...>::value;

template <typename... Ds>
using enable_if_matrices_t = std::enable_if_t<are_matrices<Ds...>::value>;

template <typename... Ds>
using enable_if_vectors_t = std::enable_if_t<are_vectors<Ds...>::value>;

/* -------------------------------------------------------------------------- */

template <typename T>
struct is_tensor : public std::is_base_of<akantu::TensorTraitBase, T> {};

template <typename PlainObjectType, int MapOptions, typename StrideType>
struct is_tensor<Eigen::Map<PlainObjectType, MapOptions, StrideType>>
    : public std::true_type {};

template <typename XprType, int BlockRows, int BlockCols, bool InnerPanel>
struct is_tensor<Eigen::Block<XprType, BlockRows, BlockCols, InnerPanel>>
    : public std::true_type {};

template <typename T, Eigen::Index m, Eigen::Index n>
struct is_tensor<Eigen::Matrix<T, m, n>> : public std::true_type {};

template <typename T, size_t n>
using is_tensor_n = std::is_base_of<akantu::TensorTrait<n>, T>;

template <class T> inline constexpr bool is_tensor_v = is_tensor<T>::value;

template <class T, size_t n>
inline constexpr bool is_tensor_n_v = is_tensor_n<T, n>::value;

template <std::size_t n, typename T = void, typename... Ts>
using enable_if_tensors_n = std::enable_if<
    aka::conjunction<
        is_tensor_n<Ts, n>...,
        std::is_same<
            std::common_type_t<std::decay_t<typename Ts::value_type>...>,
            std::decay_t<typename Ts::value_type>>...>::value,
    T>;

template <std::size_t n, typename T = void, typename... Ts>
using enable_if_tensors_n_t = typename enable_if_tensors_n<n, T, Ts...>::type;

template <typename T> struct size_type<akantu::Vector<T>> {
  using type = ::akantu::Int;
};

template <typename T> struct size_type<Eigen::Ref<akantu::Vector<T>>> {
  using type = ::akantu::Int;
};

template <typename ExpressionType, int Direction>
struct size_type<Eigen::VectorwiseOp<ExpressionType, Direction>> {
  using type = ::akantu::Int;
};

} // namespace aka

namespace akantu { // fwd declaration
template <typename T, Int ndim> class TensorBase;
template <typename T, Int ndim> class TensorProxy;
template <typename T, Int ndim> class Tensor;
} // namespace akantu

/* -------------------------------------------------------------------------- */
#include "aka_view_iterators.hh"
/* -------------------------------------------------------------------------- */
#include "aka_tensor.hh" // NOLINT(unused-includes)
/* -------------------------------------------------------------------------- */

namespace akantu {

class ArrayBase;

/* -------------------------------------------------------------------------- */
namespace details {
template <typename T> struct MapPlainObjectType {
  using type = T;
};

template <typename PlainObjectType, int MapOptions, typename StrideType>
struct MapPlainObjectType<Eigen::Map<PlainObjectType, MapOptions, StrideType>> {
  using type = PlainObjectType;
};

template <typename T>
using MapPlainObjectType_t = typename MapPlainObjectType<T>::type;

template <typename Scalar, Idx...> struct EigenMatrixViewHelper {};

template <typename Scalar, Idx RowsAtCompileTime>
struct EigenMatrixViewHelper<Scalar, RowsAtCompileTime> {
  using type = Eigen::Matrix<Scalar, RowsAtCompileTime, 1>;
};

template <typename Scalar, Idx RowsAtCompileTime, Idx ColsAtCompileTime>
struct EigenMatrixViewHelper<Scalar, RowsAtCompileTime, ColsAtCompileTime> {
  using type = Eigen::Matrix<Scalar, RowsAtCompileTime, ColsAtCompileTime>;
};

template <typename Scalar, Idx... sizes>
using EigenMatrixViewHelper_t =
    typename EigenMatrixViewHelper<Scalar, sizes...>::type;

template <typename Array, Idx... sizes> class EigenView {
  static_assert(sizeof...(sizes) == 1 or sizeof...(sizes) == 2,
                "Eigen only supports Vector and Matrices");

private:
  template <class A = Array,
            std::enable_if_t<aka::is_array<std::decay_t<A>>::value> * = nullptr>
  auto array_size() const {
    return array.get().size() * array.get().getNbComponent();
  }

  template <
      class A = Array,
      std::enable_if_t<not aka::is_array<std::decay_t<A>>::value> * = nullptr>
  auto array_size() const {
    return array.get().size();
  }

  using ArrayRef_t = decltype(std::ref(std::declval<Array>()));

public:
  using size_type = typename std::decay_t<Array>::size_type;
  using value_type = typename std::decay_t<Array>::value_type;

  EigenView(Array && array, decltype(sizes)... sizes_) // NOLINT
      : array(std::ref(array)), sizes_(sizes_...) {}

  EigenView(Array && array)
      : array(std::ref(array)), sizes_(sizes...) {} // NOLINT

  EigenView(const EigenView & other) = default;
  EigenView(EigenView && other) noexcept = default;
  ~EigenView() = default;

  auto operator=(const EigenView & other) -> EigenView & = default;
  auto operator=(EigenView && other) noexcept -> EigenView & = default;

  template <typename A = Array, std::enable_if_t<not std::is_const_v<
                                    std::remove_reference_t<A>>> * = nullptr>
  decltype(auto) begin() {
    return aka::make_from_tuple<::akantu::view_iterator<
        Eigen::Map<EigenMatrixViewHelper_t<value_type, sizes...>>>>(
        std::tuple_cat(std::make_tuple(array.get().data()), sizes_));
  }

  template <typename A = Array, std::enable_if_t<not std::is_const_v<
                                    std::remove_reference_t<A>>> * = nullptr>
  decltype(auto) end() {
    return aka::make_from_tuple<::akantu::view_iterator<
        Eigen::Map<EigenMatrixViewHelper_t<value_type, sizes...>>>>(
        std::tuple_cat(std::make_tuple(array.get().data() + array_size()),
                       sizes_));
  }

  decltype(auto) begin() const {
    return aka::make_from_tuple<::akantu::view_iterator<
        Eigen::Map<const EigenMatrixViewHelper_t<value_type, sizes...>>>>(
        std::tuple_cat(std::make_tuple(array.get().data()), sizes_));
  }
  decltype(auto) end() const {
    return aka::make_from_tuple<::akantu::view_iterator<
        Eigen::Map<const EigenMatrixViewHelper_t<value_type, sizes...>>>>(
        std::tuple_cat(std::make_tuple(array.get().data() + array_size()),
                       sizes_));
  }

private:
  ArrayRef_t array;
  std::tuple<decltype(sizes)...> sizes_;
};

} // namespace details

template <Idx RowsAtCompileTime, typename Array>
decltype(auto) make_view(Array && array, Idx rows = RowsAtCompileTime) {
  return details::EigenView<Array, RowsAtCompileTime>(
      std::forward<Array>(array), rows);
}

template <Idx RowsAtCompileTime, Idx ColsAtCompileTime, typename Array>
decltype(auto) make_view(Array && array, Idx rows = RowsAtCompileTime,
                         Idx cols = ColsAtCompileTime) {
  return details::EigenView<Array, RowsAtCompileTime, ColsAtCompileTime>(
      std::forward<Array>(array), rows, cols);
}

template <Idx RowsAtCompileTime, typename Array>
decltype(auto) make_const_view(const Array & array,
                               Idx rows = RowsAtCompileTime) {
  return make_view<RowsAtCompileTime>(array, rows);
}

template <Idx RowsAtCompileTime, Idx ColsAtCompileTime, typename Array>
decltype(auto) make_const_view(const Array & array,
                               Idx rows = RowsAtCompileTime,
                               Idx cols = ColsAtCompileTime) {
  return make_view<RowsAtCompileTime, ColsAtCompileTime>(array, rows, cols);
}

} // namespace akantu

namespace Eigen {

template <typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void MatrixBase<Derived>::zero() {
  return this->fill(0.);
}

/* -------------------------------------------------------------------------- */
/* Vector                                                                     */
/* -------------------------------------------------------------------------- */
template <typename Derived>
template <typename ED, typename T,
          std::enable_if_t<not std::is_const<T>::value and
                           ED::IsVectorAtCompileTime> *>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE decltype(auto)
MatrixBase<Derived>::begin() {
  return ::akantu::view_iterator<typename Derived::Scalar>(
      this->derived().data());
}

template <typename Derived>
template <typename ED, typename T,
          std::enable_if_t<not std::is_const<T>::value and
                           ED::IsVectorAtCompileTime> *>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE decltype(auto)
MatrixBase<Derived>::end() {
  return ::akantu::view_iterator<typename Derived::Scalar>(
      this->derived().data() + this->size());
}

template <typename Derived>
template <typename ED, std::enable_if_t<ED::IsVectorAtCompileTime> *>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE decltype(auto)
MatrixBase<Derived>::begin() const {
  using Scalar = typename Derived::Scalar;
  return ::akantu::const_view_iterator<Scalar>(this->derived().data());
}

template <typename Derived>
template <typename ED, std::enable_if_t<ED::IsVectorAtCompileTime> *>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE decltype(auto)
MatrixBase<Derived>::end() const {
  using Scalar = typename Derived::Scalar;
  return ::akantu::const_view_iterator<Scalar>(this->derived().data() +
                                               this->size());
}

/* -------------------------------------------------------------------------- */
/* Matrix                                                                     */
/* -------------------------------------------------------------------------- */
template <typename Derived>
template <typename ED, typename T,
          std::enable_if_t<not std::is_const<T>::value and
                           not ED::IsVectorAtCompileTime> *>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE decltype(auto)
MatrixBase<Derived>::begin() {
  return ::akantu::view_iterator<
      Map<Matrix<typename Derived::Scalar, Derived::RowsAtCompileTime, 1>>>(
      this->derived().data(), this->rows());
}

template <typename Derived>
template <typename ED, typename T,
          std::enable_if_t<not std::is_const<T>::value and
                           not ED::IsVectorAtCompileTime> *>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE decltype(auto)
MatrixBase<Derived>::end() {
  return ::akantu::view_iterator<
      Map<Matrix<typename Derived::Scalar, Derived::RowsAtCompileTime, 1>>>(
      this->derived().data() + this->size(), this->rows());
}

template <typename Derived>
template <typename ED, std::enable_if_t<not ED::IsVectorAtCompileTime> *>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE decltype(auto)
MatrixBase<Derived>::begin() const {
  using Scalar = typename Derived::Scalar;
  return ::akantu::const_view_iterator<
      Map<const Matrix<Scalar, Derived::RowsAtCompileTime, 1>>>(
      const_cast<Scalar *>(this->derived().data()), this->rows()); // NOLINT
}

template <typename Derived>
template <typename ED, std::enable_if_t<not ED::IsVectorAtCompileTime> *>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE decltype(auto)
MatrixBase<Derived>::end() const {
  using Scalar = typename Derived::Scalar;
  return ::akantu::const_view_iterator<
      Map<const Matrix<Scalar, Derived::RowsAtCompileTime, 1>>>(
      const_cast<Scalar *>(this->derived().data()) + this->size(), // NOLINT
      this->rows());
}

template <typename Derived>
template <typename OtherScalar>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE auto MatrixBase<Derived>::eig() const {
  EigenSolver<akantu::details::MapPlainObjectType_t<std::decay_t<Derived>>>
      solver(*this, false);
  if constexpr (std::is_floating_point<OtherScalar>{}) {
    return solver.eigenvalues().real();
  } else {
    return solver.eigenvalues();
  }
}

template <typename Derived>
template <typename OtherDerived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void
MatrixBase<Derived>::eig(MatrixBase<OtherDerived> & values) const {
  EigenSolver<akantu::details::MapPlainObjectType_t<std::decay_t<Derived>>>
      solver(*this, false);
  using OtherScalar = typename OtherDerived::Scalar;
  // as advised by the Eigen developers even though this is a UB
  // auto & values = const_cast<MatrixBase<OtherDerived> &>(values_);
  if constexpr (std::is_floating_point<OtherScalar>{}) {
    values = solver.eigenvalues().real();
  } else {
    values = solver.eigenvalues();
  }
}

template <typename Derived>
template <typename D1, typename D2>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void
MatrixBase<Derived>::eig(MatrixBase<D1> & values, MatrixBase<D2> & vectors,
                         bool sort) const {
  EigenSolver<akantu::details::MapPlainObjectType_t<std::decay_t<Derived>>>
      solver(*this, true);

  // as advised by the Eigen developers even though this is a UB
  // auto & values = const_cast<MatrixBase<D1> &>(values_);
  // auto & vectors = const_cast<MatrixBase<D2> &>(vectors_);

  auto norm = this->norm();

  using OtherScalar = typename D1::Scalar;

  if ((solver.eigenvectors().imag().template lpNorm<Infinity>() >
       1e-15 * norm) and
      std::is_floating_point<OtherScalar>::value) {
    AKANTU_EXCEPTION("This matrix has complex eigenvectors()");
  }

  if (not sort) {
    if constexpr (std::is_floating_point<OtherScalar>{}) {
      values = solver.eigenvalues().real();
      vectors = solver.eigenvectors().real();
    } else {
      values = solver.eigenvalues();
      vectors = solver.eigenvectors();
    }
    return;
  }

  if (not std::is_floating_point<OtherScalar>::value) {
    AKANTU_EXCEPTION("Cannot sort complex eigen values");
  }

  values = solver.eigenvalues().real();

  PermutationMatrix<Dynamic> P(values.size());
  P.setIdentity();

  std::sort(P.indices().begin(), P.indices().end(),
            [&values](const Index & a, const Index & b) {
              return (values(a) - values(b)) > 0;
            });

  if constexpr (std::is_floating_point<OtherScalar>{}) {
    values = P.transpose() * values;
    vectors = solver.eigenvectors().real() * P;
  }
}

template <typename Derived>
template <typename OtherDerived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void
MatrixBase<Derived>::eigh(const MatrixBase<OtherDerived> & values_) const {
  SelfAdjointEigenSolver<
      akantu::details::MapPlainObjectType_t<std::decay_t<Derived>>>
      solver(*this, EigenvaluesOnly);
  // as advised by the Eigen developers even though this is a UB
  auto & values = const_cast<MatrixBase<OtherDerived> &>(values_); // NOLINT
  values = solver.eigenvalues();
}

template <typename Derived>
template <typename D1, typename D2>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void
MatrixBase<Derived>::eigh(const MatrixBase<D1> & values_,
                          const MatrixBase<D2> & vectors_, bool sort) const {
  SelfAdjointEigenSolver<
      akantu::details::MapPlainObjectType_t<std::decay_t<Derived>>>
      solver(*this, ComputeEigenvectors);

  // as advised by the Eigen developers, even though this is a UB
  auto & values = const_cast<MatrixBase<D1> &>(values_);   // NOLINT
  auto & vectors = const_cast<MatrixBase<D2> &>(vectors_); // NOLINT

  if (not sort) {
    values = solver.eigenvalues();
    vectors = solver.eigenvectors();
    return;
  }

  values = solver.eigenvalues();
  PermutationMatrix<Dynamic> P(values.size());
  P.setIdentity();

  std::sort(P.indices().begin(), P.indices().end(),
            [&values](const Index & a, const Index & b) {
              return (values(a) - values(b)) > 0;
            });

  values = P.transpose() * values;
  vectors = solver.eigenvectors() * P; // permutes the columns (eigen vectors)
}

} // namespace Eigen

namespace std {
template <typename POD1, typename POD2, int MapOptions, typename StrideType>
struct is_convertible<Eigen::Map<POD1, MapOptions, StrideType>,
                      Eigen::Map<POD2, MapOptions, StrideType>>
    : aka::bool_constant<is_convertible<POD1, POD2>::value> {};
} // namespace std

#endif /* AKANTU_AKA_TYPES_HH */
