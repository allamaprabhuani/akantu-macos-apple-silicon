/**
 * Copyright (©) 2018-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include <memory>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_AKA_VIEW_ITERATORS_HH__
#define __AKANTU_AKA_VIEW_ITERATORS_HH__

namespace akantu {
template <typename T, Int ndim> class TensorBase;
}

namespace akantu {
/* -------------------------------------------------------------------------- */
/* Iterators                                                                  */
/* -------------------------------------------------------------------------- */
namespace detail {
  template <typename... V> constexpr auto product_all(V &&... v) {
    std::common_type_t<int, V...> result = 1;
    (void)std::initializer_list<int>{(result *= v, 0)...};
    return result;
  }

  template <class R> struct IteratorHelper { static constexpr Int dim = 0; };

  template <class Derived> struct IteratorHelper<Eigen::MatrixBase<Derived>> {
  private:
    using T = typename Derived::Scalar;
    static constexpr Int m = Derived::RowsAtCompileTime;
    static constexpr Int n = Derived::ColsAtCompileTime;

  public:
    static constexpr Int dim =
        Eigen::MatrixBase<Derived>::IsVectorAtCompileTime ? 1 : 2;
    using pointer = T *;
    using proxy = Eigen::Map<Eigen::Matrix<T, m, n>>;
    using const_proxy = Eigen::Map<const Eigen::Matrix<T, m, n>>;
  };

  template <class Derived> struct IteratorHelper<Eigen::Map<Derived>> {
  private:
    using T = typename Derived::Scalar;
    static constexpr Int m = Derived::RowsAtCompileTime;
    static constexpr Int n = Derived::ColsAtCompileTime;

  public:
    static constexpr Int dim =
        Derived::IsVectorAtCompileTime and m != 1 ? 1 : 2;
    using pointer =
        std::conditional_t<std::is_const_v<Derived>, const T *, T *>;
    using proxy = Eigen::Map<Derived>;
    using const_proxy = Eigen::Map<const Derived>;
  };

  template <typename T, Int _dim> struct IteratorHelper<TensorBase<T, _dim>> {
    static constexpr Int dim = _dim;
    using pointer = T *;
    using proxy = TensorProxy<T, _dim>;
    using const_proxy = TensorProxy<const T, _dim>;
  };

  template <typename T, Int _dim> struct IteratorHelper<TensorProxy<T, _dim>> {
    static constexpr Int dim = _dim;
    using pointer = T *;
    using proxy = TensorProxy<T, _dim>;
    using const_proxy = TensorProxy<const T, _dim>;
  };

  /* ------------------------------------------------------------------------ */
  template <class R, class daughter, class IR = R,
            Int dim = IteratorHelper<std::decay_t<R>>::dim>
  class internal_view_iterator {
  protected:
    using helper = IteratorHelper<std::decay_t<R>>;
    using internal_value_type = IR;
    using internal_pointer = IR *;
    using scalar_pointer = typename helper::pointer;
    using proxy_t = typename helper::proxy;
    static constexpr int dim_ = dim;

  public:
    using pointer = proxy_t *;
    using value_type = proxy_t;
    using reference = proxy_t &;
    using difference_type = Int;
    using iterator_category = std::random_access_iterator_tag;

  private:
    template <class ProxyType, std::size_t... I>
    constexpr auto get_new_proxy(scalar_pointer data,
                                 std::index_sequence<I...>) const {
      return ProxyType(data, dims[I]...);
    }

    constexpr auto get_new_proxy(scalar_pointer data) {
      return this->template get_new_proxy<proxy_t>(
          data, std::make_index_sequence<dim>());
    }

    constexpr auto get_new_proxy(scalar_pointer data) const {
      return this->template get_new_proxy<proxy_t>(
          data, std::make_index_sequence<dim>());
    }

    template <std::size_t... I>
    constexpr void reset_proxy(std::index_sequence<I...>) {
      new (&proxy) proxy_t(ret_ptr, dims[I]...);
    }

    constexpr auto reset_proxy() {
      return reset_proxy(std::make_index_sequence<dim>());
    }

  protected:
    template <typename OR, typename OD, typename OIR,
              std::enable_if_t<std::is_convertible<
                  decltype(std::declval<OIR>().data()),
                  decltype(std::declval<IR>().data())>::value> * = nullptr>
    explicit internal_view_iterator(
        internal_view_iterator<OR, OD, OIR, dim> & it)
        : dims(it.dims), _offset(it._offset), initial(it.initial),
          ret_ptr(it.ret_ptr), proxy(get_new_proxy(it.ret_ptr)) {}

    template <typename OR, typename OD, typename OIR, Int _dim>
    friend class internal_view_iterator;

    template <typename... Args>
    using valid_args_t = std::enable_if_t<
        aka::conjunction<aka::disjunction<std::is_integral<Args>,
                                          std::is_enum<Args>>...>::value and
            dim == sizeof...(Args),
        int>;

  public:
    /// Generic constructor dor any tensor dimension
    template <typename... Ns, valid_args_t<Ns...> = 0>
    internal_view_iterator(scalar_pointer data, Ns... ns)
        : dims({Int(ns)...}),
          _offset(detail::product_all(std::forward<Ns>(ns)...)), initial(data),
          ret_ptr(data), proxy(data, ns...) {}

    // Specific constructor for Vector of static size 1
    template <typename RD = std::decay_t<R>,
              std::enable_if_t<aka::is_eigen_map<RD>::value and
                               RD::RowsAtCompileTime == 1 and
                               RD::ColsAtCompileTime == 1> * = nullptr>
    constexpr internal_view_iterator(scalar_pointer data, Idx rows)
        : dims({rows, 1}), _offset(rows), initial(data), ret_ptr(data),
          proxy(data, rows, 1) {
      assert(rows == 1 && "1x1 Matrix");
    }

    /// Specific constructor for Eigen::Map<Matrix> that look like
    /// Eigen::Map<Vector>
    template <typename RD = std::decay_t<R>,
              std::enable_if_t<(RD::RowsAtCompileTime != 1) and
                               RD::ColsAtCompileTime == 1> * = nullptr>
    constexpr internal_view_iterator(scalar_pointer data, Idx rows,
                                     [[gnu::unused]] Idx cols)
        : dims({rows}), _offset(rows), initial(data), ret_ptr(data),
          proxy(data, rows, 1) {
      assert(cols == 1 && "nx1 Matrix");
    }

    /// Default constructor for Eigen::Map<Vector>
    template <Int _dim = dim, typename RD = std::decay_t<R>,
              std::enable_if_t<_dim == 1> * = nullptr>
    internal_view_iterator()
        : proxy(reinterpret_cast<scalar_pointer>(0xdeadbeaf),
                RD::RowsAtCompileTime == Eigen::Dynamic
                    ? 1
                    : RD::RowsAtCompileTime) {
      // initialized to a fake pointer to pass the static_assert in Eigen
      // this proxy should not be returned
    }

    /// Default constructor for Eigen::Map<Matrix>
    template <Int _dim = dim, typename RD = std::decay_t<R>,
              std::enable_if_t<_dim == 2> * = nullptr>
    internal_view_iterator()
        : proxy(reinterpret_cast<scalar_pointer>(0xdeadbeaf),
                RD::RowsAtCompileTime == Eigen::Dynamic ? 1
                                                        : RD::RowsAtCompileTime,
                RD::ColsAtCompileTime == Eigen::Dynamic
                    ? 1
                    : RD::ColsAtCompileTime) {
      // initialized to a fake pointer to pass the `static_assert` in `Eigen
      // this proxy should not be returned
    }

    template <Int _dim = dim, typename RD = std::decay_t<R>,
              std::enable_if_t<(_dim > 2)> * = nullptr>
    internal_view_iterator() {}

    internal_view_iterator(const internal_view_iterator & it)
        : dims(it.dims), _offset(it._offset), initial(it.initial),
          ret_ptr(it.ret_ptr), proxy(get_new_proxy(it.ret_ptr)) {}

    internal_view_iterator &
    operator=(internal_view_iterator && it) noexcept = default;
    internal_view_iterator(internal_view_iterator && it) noexcept = default;

    virtual ~internal_view_iterator() = default;

    template <typename OR, typename OD, typename OIR,
              std::enable_if_t<std::is_convertible<
                  decltype(std::declval<OIR>().data()),
                  decltype(std::declval<IR>().data())>::value> * = nullptr>
    inline internal_view_iterator &
    operator=(const internal_view_iterator<OR, OD, OIR, dim> & it) {
      this->dims = it.dims;
      this->_offset = it._offset;
      this->initial = it.initial;
      this->ret_ptr = it.ret_ptr;
      reset_proxy();
      return *this;
    }

    inline internal_view_iterator &
    operator=(const internal_view_iterator & it) {
      if (this != &it) {
        this->dims = it.dims;
        this->_offset = it._offset;
        this->initial = it.initial;
        this->ret_ptr = it.ret_ptr;
        reset_proxy();
      }
      return *this;
    }

  public:
    Idx getCurrentIndex() {
      return (this->ret_ptr - this->initial) / this->_offset;
    }

    inline reference operator*() {
      this->reset_proxy();
      return proxy;
    }

    inline reference operator*() const {
      this->reset_proxy();
      return proxy;
    }


    inline pointer operator->() {
      reset_proxy();
      return &proxy;
    }

    inline daughter & operator++() {
      ret_ptr += _offset;
      return static_cast<daughter &>(*this);
    }

    inline daughter & operator--() {
      ret_ptr -= _offset;
      return static_cast<daughter &>(*this);
    }

    inline daughter & operator+=(Idx n) {
      ret_ptr += _offset * n;
      return static_cast<daughter &>(*this);
    }

    inline daughter & operator-=(Idx n) {
      ret_ptr -= _offset * n;
      return static_cast<daughter &>(*this);
    }

    inline auto operator[](Idx n) {
      return get_new_proxy(ret_ptr + n * _offset);
    }

    inline auto operator[](Idx n) const {
      return get_new_proxy(ret_ptr + n * _offset);
    }

    inline bool operator==(const internal_view_iterator & other) const {
      return this->ret_ptr == other.ret_ptr;
    }

    inline bool operator!=(const internal_view_iterator & other) const {
      return this->ret_ptr != other.ret_ptr;
    }

    inline bool operator<(const internal_view_iterator & other) const {
      return this->ret_ptr < other.ret_ptr;
    }

    inline bool operator<=(const internal_view_iterator & other) const {
      return this->ret_ptr <= other.ret_ptr;
    }

    inline bool operator>(const internal_view_iterator & other) const {
      return this->ret_ptr > other.ret_ptr;
    }

    inline bool operator>=(const internal_view_iterator & other) const {
      return this->ret_ptr >= other.ret_ptr;
    }

    inline auto operator+(difference_type n) const {
      daughter tmp(static_cast<const daughter &>(*this));
      tmp += n;
      return tmp;
    }

    inline auto operator-(difference_type n) const {
      daughter tmp(static_cast<const daughter &>(*this));
      tmp -= n;
      return tmp;
    }

    inline difference_type operator-(const internal_view_iterator & b) const {
      return (this->ret_ptr - b.ret_ptr) / _offset;
    }

    inline scalar_pointer data() const { return ret_ptr; }
    inline difference_type offset() const { return _offset; }
    inline decltype(auto) getDims() const { return dims; }

  protected:
    std::array<Int, dim> dims;
    difference_type _offset{0};
    scalar_pointer initial{nullptr};
    scalar_pointer ret_ptr{nullptr};
    proxy_t proxy;
  };

  /* ------------------------------------------------------------------------ */
  /**
   * Specialization for scalar types
   */
  template <class R, class daughter, class IR>
  class internal_view_iterator<R, daughter, IR, 0> {
  public:
    using value_type = R;
    using pointer = R *;
    using reference = R &;
    using const_reference = const R &;
    using difference_type = Idx; // std::ptrdiff_t;
    using iterator_category = std::random_access_iterator_tag;
    static constexpr int dim_ = 0;

  protected:
    using internal_value_type = IR;
    using internal_pointer = IR *;

  public:
    internal_view_iterator(pointer data = nullptr) : ret(data), initial(data){};
    internal_view_iterator(const internal_view_iterator & it) = default;
    internal_view_iterator(internal_view_iterator && it) = default;

    virtual ~internal_view_iterator() = default;

    inline internal_view_iterator &
    operator=(const internal_view_iterator & it) = default;

    Idx getCurrentIndex() { return (this->ret - this->initial); };

    inline reference operator*() { return *ret; }
    inline reference operator*() const { return *ret; }

    inline pointer operator->() { return ret; };

    inline daughter & operator++() {
      ++ret;
      return static_cast<daughter &>(*this);
    }
    inline daughter & operator--() {
      --ret;
      return static_cast<daughter &>(*this);
    }

    inline daughter & operator+=(const Idx n) {
      ret += n;
      return static_cast<daughter &>(*this);
    }
    inline daughter & operator-=(const Idx n) {
      ret -= n;
      return static_cast<daughter &>(*this);
    }

    inline reference operator[](const Idx n) { return ret[n]; }

    inline bool operator==(const internal_view_iterator & other) const {
      return ret == other.ret;
    }
    inline bool operator!=(const internal_view_iterator & other) const {
      return ret != other.ret;
    }
    inline bool operator<(const internal_view_iterator & other) const {
      return ret < other.ret;
    }
    inline bool operator<=(const internal_view_iterator & other) const {
      return ret <= other.ret;
    }
    inline bool operator>(const internal_view_iterator & other) const {
      return ret > other.ret;
    }
    inline bool operator>=(const internal_view_iterator & other) const {
      return ret >= other.ret;
    }

    inline daughter operator-(difference_type n) const {
      return daughter(ret - n);
    }
    inline daughter operator+(difference_type n) const {
      return daughter(ret + n);
    }

    inline difference_type operator-(const internal_view_iterator & b) const {
      return ret - b.ret;
    }

    inline pointer data() const { return ret; }
    inline decltype(auto) getDims() const { return dims; }

  protected:
    std::array<int, 0> dims;
    pointer ret{nullptr};
    pointer initial{nullptr};
  };
} // namespace detail

/* -------------------------------------------------------------------------- */
template <typename R> class view_iterator;

template <typename R>
class const_view_iterator
    : public detail::internal_view_iterator<const R, const_view_iterator<R>,
                                            R> {
public:
  using parent =
      detail::internal_view_iterator<const R, const_view_iterator, R>;
  using value_type = typename parent::value_type;
  using pointer = typename parent::pointer;
  using reference = typename parent::reference;
  using difference_type = typename parent::difference_type;
  using iterator_category = typename parent::iterator_category;

protected:
  template <typename Iterator, std::size_t... I>
  static inline auto convert_helper(const Iterator & it,
                                    std::index_sequence<I...>) {
    return const_view_iterator(it.data(), it.getDims()[I]...);
  }
  template <typename Iterator> static inline auto convert(const Iterator & it) {
    return convert_helper(it, std::make_index_sequence<parent::dim_>());
  }

public:
  const_view_iterator() : parent(){};
  const_view_iterator(const const_view_iterator & it) = default;
  const_view_iterator(const_view_iterator && it) noexcept = default;

  template <typename P, typename... Ns>
  const_view_iterator(P * data, Ns... ns) : parent(data, ns...) {}

  const_view_iterator & operator=(const const_view_iterator & it) = default;

  template <typename OR,
            std::enable_if_t<not std::is_same<R, OR>::value> * = nullptr>
  const_view_iterator(const const_view_iterator<OR> & it)
      : parent(convert(it)) {}

  template <typename OR,
            std::enable_if_t<std::is_convertible<R, OR>::value> * = nullptr>
  const_view_iterator(const view_iterator<OR> & it) : parent(convert(it)) {}

  template <typename OR,
            std::enable_if_t<not std::is_same<R, OR>::value and
                             std::is_convertible<R, OR>::value> * = nullptr>
  const_view_iterator & operator=(const const_view_iterator<OR> & it) {
    return dynamic_cast<const_view_iterator &>(parent::operator=(it));
  }

  template <typename OR,
            std::enable_if_t<std::is_convertible<R, OR>::value> * = nullptr>
  const_view_iterator operator=(const view_iterator<OR> & it) {
    return convert(it);
  }
};

template <class R, bool is_tensor_ = aka::is_tensor<R>::value>
struct ConstConverterIteratorHelper {
protected:
  template <std::size_t... I>
  static inline auto convert_helper(const view_iterator<R> & it,
                                    std::index_sequence<I...>) {
    return const_view_iterator<R>(it.data(), it.getDims()[I]...);
  }

public:
  static inline auto convert(const view_iterator<R> & it) {
    return convert_helper(
        it, std::make_index_sequence<
                std::tuple_size<decltype(it.getDims())>::value>());
  }
};

template <class R> struct ConstConverterIteratorHelper<R, false> {
  static inline auto convert(const view_iterator<R> & it) {
    return const_view_iterator<R>(it.data());
  }
};

template <typename R>
class view_iterator
    : public detail::internal_view_iterator<R, view_iterator<R>> {
public:
  using parent = detail::internal_view_iterator<R, view_iterator>;
  using value_type = typename parent::value_type;
  using pointer = typename parent::pointer;
  using reference = typename parent::reference;
  using difference_type = typename parent::difference_type;
  using iterator_category = typename parent::iterator_category;

public:
  view_iterator() : parent(){};
  view_iterator(const view_iterator & it) = default;
  view_iterator(view_iterator && it) = default;

  template <typename P, typename... Ns>
  view_iterator(P * data, Ns... ns) : parent(data, ns...) {}

  view_iterator & operator=(const view_iterator & it) = default;

  operator const_view_iterator<R>() {
    return ConstConverterIteratorHelper<R>::convert(*this);
  }
};

namespace {
  template <std::size_t dim, typename T> struct ViewIteratorHelper {
    using type = TensorProxy<T, dim>;
  };

  template <typename T> struct ViewIteratorHelper<0, T> { using type = T; };
  template <typename T> struct ViewIteratorHelper<1, T> {
    using type = Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, 1>>;
  };
  template <typename T> struct ViewIteratorHelper<1, const T> {
    using type = Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, 1>>;
  };

  template <typename T> struct ViewIteratorHelper<2, T> {
    using type = Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>;
  };

  template <typename T> struct ViewIteratorHelper<2, const T> {
    using type =
        Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>;
  };

  template <std::size_t dim, typename T>
  using ViewIteratorHelper_t = typename ViewIteratorHelper<dim, T>::type;
} // namespace

} // namespace akantu

#include <iterator>

namespace std {

template <typename R> struct iterator_traits<::akantu::const_view_iterator<R>> {
protected:
  using iterator = ::akantu::const_view_iterator<R>;

public:
  using iterator_category = typename iterator::iterator_category;
  using value_type = typename iterator::value_type;
  using difference_type = typename iterator::difference_type;
  using pointer = typename iterator::pointer;
  using reference = typename iterator::reference;
};

template <typename R> struct iterator_traits<::akantu::view_iterator<R>> {
protected:
  using iterator = ::akantu::view_iterator<R>;

public:
  using iterator_category = typename iterator::iterator_category;
  using value_type = typename iterator::value_type;
  using difference_type = typename iterator::difference_type;
  using pointer = typename iterator::pointer;
  using reference = typename iterator::reference;
};

} // namespace std

#endif /* !__AKANTU_AKA_VIEW_ITERATORS_HH__ */
