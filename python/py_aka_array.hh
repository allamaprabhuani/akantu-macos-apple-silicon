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
#include <aka_array.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_PY_AKA_ARRAY_HH_
#define AKANTU_PY_AKA_ARRAY_HH_

namespace py = pybind11;
namespace _aka = akantu;

namespace akantu {
namespace detail {
  template <class T> struct is_array_type : public std::false_type {};
  // template <class T> struct is_array_type<Vector<T>> : public std::true_type
  // {}; template <class T> struct is_array_type<Matrix<T>> : public
  // std::true_type {};
  template <class T> struct is_array_type<Array<T>> : public std::true_type {};

  /* ------------------------------------------------------------------------ */
  template <typename T> class ArrayProxy : public Array<T> {
  protected:
    // deallocate the memory
    void deallocate() final {}

    // allocate the memory
    void allocate(Int /*size*/, Int /*nb_component*/) final {}

    // allocate and initialize the memory
    void allocate(Int /*size*/, Int /*nb_component*/,
                  const T & /*value*/) final {}

  public:
    ArrayProxy(T * data, Int size, Int nb_component) {
      this->values = data;
      this->size_ = size;
      this->nb_component = nb_component;
    }

    ArrayProxy(const Array<T> & src) {
      this->values = src.data();
      this->size_ = src.size();
      this->nb_component = src.getNbComponent();
    }

    ~ArrayProxy() override { this->values = nullptr; }

    void resize(Int size, const T & /*val*/) final {
      if (size != this->size()) {
        AKANTU_EXCEPTION("cannot resize a temporary array");
      }
      // std::fill(this->begin(), this->end(), val);
    }

    void resize(Int new_size) final {
      if (new_size != this->size()) {
        AKANTU_EXCEPTION("cannot resize a temporary array "
                         << new_size << " != " << this->size());
      }
    }

    void reserve(Int /*size*/, Int /*new_size*/) final {
      AKANTU_EXCEPTION("cannot resize a temporary array");
    }
  };

  /* ------------------------------------------------------------------------ */
  template <typename T> struct ProxyType {};
  // template <typename T> struct ProxyType<Vector<T>> { using type = Vector<T>;
  // }; template <typename T> struct ProxyType<Matrix<T>> { using type =
  // Matrix<T>; };
  template <typename T> struct ProxyType<Array<T>> {
    using type = ArrayProxy<T>;
  };
  template <typename array> using ProxyType_t = typename ProxyType<array>::type;
} // namespace detail
} // namespace akantu

namespace pybind11 {
namespace detail {

  template <typename T> struct AkaArrayType {
    using type =
        array_t<typename T::value_type, array::c_style | array::forcecast>;
  };

  template <typename U> using array_type_t = typename AkaArrayType<U>::type;

  template <typename T>
  decltype(auto) create_proxy(array_type_t<_aka::Array<T>> & ref,
                              const _aka::Array<T> * /*unused*/) {
    return std::make_unique<_aka::detail::ProxyType_t<_aka::Array<T>>>(
        ref.mutable_data(), ref.shape(0), ref.shape(1));
  }

  /* ------------------------------------------------------------------------ */
  template <typename T>
  py::handle aka_array_cast(const _aka::Array<T> & src,
                            py::handle base = handle(), bool writeable = true) {
    array a;
    a = array_type_t<_aka::Array<T>>({src.size(), src.getNbComponent()},
                                     src.data(), base);

    if (not writeable) {
      array_proxy(a.ptr())->flags &= ~detail::npy_api::NPY_ARRAY_WRITEABLE_;
    }

    return a.release();
  }

  /* ------------------------------------------------------------------------ */
  template <typename AkaArrayType>
  class type_caster<
      AkaArrayType,
      std::enable_if_t<_aka::detail::is_array_type<AkaArrayType>::value>> {
  protected:
    using T = typename AkaArrayType::value_type;
    using type = AkaArrayType;
    using proxy_type = _aka::detail::ProxyType_t<AkaArrayType>;
    using array_type = array_type_t<AkaArrayType>;

    std::unique_ptr<proxy_type> array_proxy;
    array_type_t<AkaArrayType> copy_or_ref;

  public:
#if PYBIND11_VERSION_MAJOR >= 2 && PYBIND11_VERSION_MINOR >= 3
    static constexpr auto name = _("AkaArray");
    operator type &&() && { return std::move(*array_proxy); }
    template <typename T_>
    using cast_op_type = pybind11::detail::movable_cast_op_type<T_>;
#else
    static PYBIND11_DESCR name() { return type_descr(_("AkaArray")); };
    template <typename _T>
    using cast_op_type = pybind11::detail::cast_op_type<_T>;
#endif
    operator type *() { return array_proxy.get(); }
    operator type &() { return *array_proxy; }

    /**
     * Conversion part 1 (Python->C++)
     */
    bool load(handle src, bool convert) {
      bool need_copy = not isinstance<array_type>(src);

      auto && fits = [&](auto && aref) {
        auto && dims = aref.ndim();
        if (dims < 1 || dims > 2) {
          return false;
        }

        return true;
      };

      if (not need_copy) {
        // We don't need a converting copy, but we also need to check whether
        // the strides are compatible with the Ref's stride requirements
        auto aref = py::cast<array_type>(src);

        if (not fits(aref)) {
          return false;
        }
        copy_or_ref = std::move(aref);
      } else {
        if (not convert) {
          return false;
        }

        auto copy = array_type::ensure(src);
        if (not copy) {
          return false;
        }

        if (not fits(copy)) {
          return false;
        }
        copy_or_ref = std::move(array_type::ensure(src));
        loader_life_support::add_patient(copy_or_ref);
      }

      AkaArrayType * dispatch = nullptr; // cannot detect T from the expression
      array_proxy = create_proxy(copy_or_ref, dispatch);
      return true;
    }

    /**
     * Conversion part 2 (C++ -> Python)
     */
    static handle cast(const type & src, return_value_policy policy,
                       handle parent) {
      switch (policy) {
      case return_value_policy::copy:
        return aka_array_cast<T>(src);
      case return_value_policy::reference_internal:
        return aka_array_cast<T>(src, parent);
      case return_value_policy::reference:
      case return_value_policy::automatic:
      case return_value_policy::automatic_reference:
        return aka_array_cast<T>(src, none());
      default:
        pybind11_fail("Invalid return_value_policy for ArrayProxy type");
      }
    }
  };
} // namespace detail
} // namespace pybind11

#endif /* AKANTU_PY_AKA_ARRAY_HH_ */
