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
#include "aka_types.hh"
#include "aka_view_iterators.hh"
/* -------------------------------------------------------------------------- */
#include <typeinfo>
#include <vector>
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
#ifndef AKANTU_ARRAY_HH_
#define AKANTU_ARRAY_HH_

namespace akantu {

/// class that afford to store vectors in static memory
// NOLINTNEXTLINE(cppcoreguidelines-special-member-functions)
class ArrayBase {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  using size_type = Int;

  explicit ArrayBase(const ID & id = "") : id(id) {}
  ArrayBase(const ArrayBase & other, const ID & id = "") {
    this->id = (id.empty()) ? other.id : id;
  }

  ArrayBase(ArrayBase && other) = default;
  ArrayBase & operator=(const ArrayBase & other) = default;
  ArrayBase & operator=(ArrayBase && other) noexcept = default;

  virtual ~ArrayBase() = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// get the amount of space allocated in bytes
  virtual Int getMemorySize() const = 0;

  // changed empty to match std::vector empty
  [[nodiscard]] inline bool empty() const { return size_ == 0; }

  /// function to print the content of the class
  virtual void printself(std::ostream & stream, int indent = 0) const = 0;

  /* ------------------------------------------------------------------------ */
  /* Accessors */
  /* ------------------------------------------------------------------------ */
public:
  /// Get the Size of the Array
  [[nodiscard]] decltype(auto) size() const { return size_; }
  /// Get the number of components
  [[nodiscard]] decltype(auto) getNbComponent() const { return nb_component; }
  /// Get the name of the array
  AKANTU_GET_MACRO_AUTO(ID, id);
  /// Set the name of th array
  AKANTU_SET_MACRO(ID, id, const ID &);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// id of the vector
  ID id;

  /// the size used
  Int size_{0};

  /// number of components
  Int nb_component{1};
};

/* -------------------------------------------------------------------------- */
/* Memory handling layer                                                      */
/* -------------------------------------------------------------------------- */
enum class ArrayAllocationType {
  _default,
  _pod,
};

template <typename T>
struct ArrayAllocationTrait
    : public std::conditional_t<
          aka::is_scalar<T>::value,
          std::integral_constant<ArrayAllocationType,
                                 ArrayAllocationType::_pod>,
          std::integral_constant<ArrayAllocationType,
                                 ArrayAllocationType::_default>> {};

/* -------------------------------------------------------------------------- */
template <typename T,
          ArrayAllocationType allocation_trait = ArrayAllocationTrait<T>::value>
class ArrayDataLayer : public ArrayBase {
public:
  using value_type = T;
  using size_type = typename ArrayBase::size_type;
  using reference = value_type &;
  using pointer_type = value_type *;
  using const_reference = const value_type &;

public:
  ~ArrayDataLayer() override = default;

  /// Allocation of a new vector
  explicit ArrayDataLayer(Int size = 0, Int nb_component = 1,
                          const ID & id = "");

  /// Allocation of a new vector with a default value
  ArrayDataLayer(Int size, Int nb_component, const_reference value,
                 const ID & id = "");

  /// Copy constructor (deep copy)
  ArrayDataLayer(const ArrayDataLayer & vect, const ID & id = "");

  /// Copy constructor (deep copy)
  explicit ArrayDataLayer(const std::vector<value_type> & vect);

  // copy operator
  ArrayDataLayer & operator=(const ArrayDataLayer & other);

  // move constructor
  ArrayDataLayer(ArrayDataLayer && other) noexcept = default;

  // move assign
  ArrayDataLayer & operator=(ArrayDataLayer && other) noexcept = default;

protected:
  // deallocate the memory
  virtual void deallocate() {}

  // allocate the memory
  virtual void allocate(Int size, Int nb_component);

  // allocate and initialize the memory
  virtual void allocate(Int size, Int nb_component, const T & value);

public:
  /// append a tuple of size nb_component containing value
  inline void push_back(const_reference value);
  /// append a vector
  // inline void push_back(const value_type new_elem[]);

  /// append a Vector or a Matrix
  template <typename Derived>
  inline void push_back(const Eigen::MatrixBase<Derived> & new_elem);

  /// changes the allocated size but not the size, if new_size = 0, the size is
  /// set to min(current_size and reserve size)
  virtual void reserve(Int size, Int new_size = Int(-1));

  /// change the size of the Array
  virtual void resize(Int size);

  /// change the size of the Array and initialize the values
  virtual void resize(Int size, const T & val);

  /// get the amount of space allocated in bytes
  inline Int getMemorySize() const override;

  /// Get the real size allocated in memory
  inline Int getAllocatedSize() const;

  /// give the address of the memory allocated for this vector
  [[deprecated("use data instead to be stl compatible")]] T * storage() const {
    return values;
  };

  const T * data() const { return values; };
  T * data() { return values; };

protected:
  /// allocation type agnostic  data access
  T * values{nullptr};

  /// data storage
  std::vector<T> data_storage;
};

/* -------------------------------------------------------------------------- */
/* Actual Array                                                               */
/* -------------------------------------------------------------------------- */
template <typename T, bool is_scal> class Array : public ArrayDataLayer<T> {
private:
  using parent = ArrayDataLayer<T>;
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  using value_type = typename parent::value_type;
  using size_type = typename parent::size_type;
  using reference = typename parent::reference;
  using pointer_type = typename parent::pointer_type;
  using const_reference = typename parent::const_reference;
  using array_type = Array<T>;

  ~Array() override;

  Array() : Array(0) {};

  /// Allocation of a new vector
  explicit Array(Int size, Int nb_component = 1, const ID & id = "");

  /// Allocation of a new vector with a default value
  explicit Array(Int size, Int nb_component, const_reference value,
                 const ID & id = "");

  /// Copy constructor
  Array(const Array & vect, const ID & id = "");

  /// Copy constructor (deep copy)
  explicit Array(const std::vector<T> & vect);

  // copy operator
  Array & operator=(const Array & other);

  // move constructor
  Array(Array && other) noexcept = default;

  // move assign
  Array & operator=(Array && other) noexcept = default;

  /* ------------------------------------------------------------------------ */
  /* Iterator                                                                 */
  /* ------------------------------------------------------------------------ */
  /// iterator for Array of nb_component = 1
  using scalar_iterator = view_iterator<T>;
  /// const_iterator for Array of nb_component = 1
  using const_scalar_iterator = const_view_iterator<const T>;

  /// iterator returning Vectors of size n  on entries of Array with
  /// nb_component = n
  using vector_iterator = view_iterator<VectorProxy<T>>;
  /// const_iterator returning Vectors of n size on entries of Array with
  /// nb_component = n
  using const_vector_iterator = const_view_iterator<VectorProxy<const T>>;

  /// iterator returning Matrices of size (m, n) on entries of Array with
  /// nb_component = m*n
  using matrix_iterator = view_iterator<MatrixProxy<T>>;
  /// const iterator returning Matrices of size (m, n) on entries of Array with
  /// nb_component = m*n
  using const_matrix_iterator = const_view_iterator<MatrixProxy<const T>>;

  /* ------------------------------------------------------------------------ */
  template <typename... Ns> inline auto begin(Ns &&... n);
  template <typename... Ns> inline auto end(Ns &&... n);
  template <typename... Ns> inline auto begin(Ns &&... n) const;
  template <typename... Ns> inline auto end(Ns &&... n) const;

  template <typename... Ns> inline auto cbegin(Ns &&... n) const;
  template <typename... Ns> inline auto cend(Ns &&... n) const;

  template <typename... Ns>
  [[deprecated("use make_view instead")]] inline auto
  begin_reinterpret(Ns &&... n);
  template <typename... Ns>
  [[deprecated("use make_view instead")]] inline auto
  end_reinterpret(Ns &&... n);
  template <typename... Ns>
  [[deprecated("use make_view instead")]] inline auto
  begin_reinterpret(Ns &&... n) const;
  template <typename... Ns>
  [[deprecated("use make_view instead")]] inline auto
  end_reinterpret(Ns &&... n) const;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// search elem in the vector, return  the position of the first occurrence or
  /// -1 if not found
  Idx find(const_reference elem) const;

  /// @see Array::find(const_reference elem) const
  //  Int find(T elem[]) const;

  /// append a value to the end of the Array
  inline void push_back(const_reference value) { parent::push_back(value); }

  /// append a Vector or a Matrix
  template <typename Derived>
  inline void push_back(const Eigen::MatrixBase<Derived> & new_elem) {
    parent::push_back(new_elem);
  }

  template <typename Ret>
  inline void push_back(const const_view_iterator<Ret> & it) {
    push_back(*it);
  }

  template <typename Ret> inline void push_back(const view_iterator<Ret> & it) {
    push_back(*it);
  }

  /// erase the value at position i
  inline void erase(Idx i);

  /// erase the entry corresponding to the iterator
  template <typename R> inline auto erase(const view_iterator<R> & it);

  /// @see Array::find(const_reference elem) const
  template <typename C, std::enable_if_t<aka::is_tensor<C>::value> * = nullptr>
  inline Idx find(const C & elem);

  /// set all entries of the array to the value t
  /// @param t value to fill the array with
  inline void set(T t) {
    std::fill_n(this->values, this->size_ * this->nb_component, t);
  }

  /// set all tuples of the array to a given vector or matrix
  /// @param vm Matrix or Vector to fill the array with
  template <typename C, std::enable_if_t<aka::is_tensor<C>::value> * = nullptr>
  inline void set(const C & vm);

  /// set the array to T{}
  inline void zero() { this->set({}); }

  /// resize the array to 0
  inline void clear() { this->resize(0); }

  /// Append the content of the other array to the current one
  void append(const Array & other);

  /// copy another Array in the current Array, the no_sanity_check allows you to
  /// force the copy in cases where you know what you do with two non matching
  /// Arrays in terms of n
  void copy(const Array & other, bool no_sanity_check = false);

  /// function to print the containt of the class
  void printself(std::ostream & stream, int indent = 0) const override;

  /// Tests if all elements are finite.
  template <typename OT = T,
            std::enable_if_t<std::is_arithmetic<OT>::value> * = nullptr>
  bool isFinite() const noexcept;

  /* ------------------------------------------------------------------------ */
  /* Operators                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// substraction entry-wise
  Array & operator-=(const Array & vect);
  /// addition entry-wise
  Array & operator+=(const Array & vect);
  /// multiply evry entry by alpha
  Array & operator*=(const T & alpha);

  /// check if the array are identical entry-wise
  bool operator==(const Array<T, is_scal> & other) const;
  /// @see Array::operator==(const Array<T, is_scal> & other) const
  bool operator!=(const Array<T, is_scal> & other) const;

  /// return a reference to the j-th entry of the i-th tuple
  inline reference operator()(Idx i, Idx j = 0);
  /// return a const reference to the j-th entry of the i-th tuple
  inline const_reference operator()(Idx i, Idx j = 0) const;

  /// return a reference to the ith component of the 1D array
  inline reference operator[](Idx i);
  /// return a const reference to the ith component of the 1D array
  inline const_reference operator[](Idx i) const;
};

/* -------------------------------------------------------------------------- */
/* Inline Functions Array<T, is_scal>                                         */
/* -------------------------------------------------------------------------- */
template <typename T, bool is_scal>
inline std::ostream & operator<<(std::ostream & stream,
                                 const Array<T, is_scal> & _this) {
  _this.printself(stream);
  return stream;
}

/* -------------------------------------------------------------------------- */
/* Inline Functions ArrayBase                                                 */
/* -------------------------------------------------------------------------- */
inline std::ostream & operator<<(std::ostream & stream,
                                 const ArrayBase & _this) {
  _this.printself(stream);
  return stream;
}

} // namespace akantu

#include "aka_array_tmpl.hh"
// IWYU pragma: begin_exports, begin_keep
#include "aka_array_filter.hh"
#include "aka_array_view.hh"
// IWYU pragma: end_exports, end_keep

#endif /* AKANTU_ARRAY_HH_ */
