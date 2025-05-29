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
// IWYU pragma: private, include "aka_array.hh"

/* -------------------------------------------------------------------------- */
#include "aka_array.hh"
#include "aka_common.hh"
/* -------------------------------------------------------------------------- */
#include <cstddef>
#include <functional>
#include <tuple>
#include <type_traits>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_ARRAY_VIEW_HH
#define AKANTU_ARRAY_VIEW_HH

namespace akantu {

/* ------------------------------------------------------------------------ */
/* Views                                                                    */
/* ------------------------------------------------------------------------ */
namespace detail {
template <typename Array, typename... Ns> class ArrayView {
  using tuple = std::tuple<Ns...>;

public:
  using size_type = Idx;
  using pointer = decltype(std::declval<Array>().data());

  ~ArrayView() = default;
  constexpr ArrayView(Array && array, Ns... ns) noexcept
      : array(array), sizes(std::move(ns)...) {}

  constexpr ArrayView(const ArrayView & array_view) = default;
  constexpr ArrayView(ArrayView && array_view) noexcept = default;

  constexpr ArrayView & operator=(const ArrayView & array_view) = default;
  constexpr ArrayView & operator=(ArrayView && array_view) noexcept = default;

  auto begin() {
    return std::apply(
        [&](auto &&... ns) {
          return detail::get_iterator(array.get(), array.get().data(),
                                      std::forward<decltype(ns)>(ns)...);
        },
        sizes);
  }

  auto begin() const {
    return std::apply(
        [&](auto &&... ns) {
          return detail::get_iterator(array.get(), array.get().data(),
                                      std::forward<decltype(ns)>(ns)...);
        },
        sizes);
  }

  auto end() {
    return std::apply(
        [&](auto &&... ns) {
          return detail::get_iterator(
              array.get(),
              array.get().data() +
                  detail::product_all(std::forward<decltype(ns)>(ns)...),
              std::forward<decltype(ns)>(ns)...);
        },
        sizes);
  }

  auto end() const {
    return std::apply(
        [&](auto &&... ns) {
          return detail::get_iterator(
              array.get(),
              array.get().data() +
                  detail::product_all(std::forward<decltype(ns)>(ns)...),
              std::forward<decltype(ns)>(ns)...);
        },
        sizes);
  }

  auto cbegin() const { return this->begin(); }
  auto cend() const { return this->end(); }

  constexpr auto size() const {
    return std::get<std::tuple_size<tuple>::value - 1>(sizes);
  }

  constexpr auto dims() const { return std::tuple_size<tuple>::value - 1; }

private:
  std::reference_wrapper<std::remove_reference_t<Array>> array;
  tuple sizes;
};

/* ---------------------------------------------------------------------- */
template <typename T, typename... Ns>
class ArrayView<const ArrayFilter<T> &, Ns...> {
  using tuple = std::tuple<Ns...>;

public:
  constexpr ArrayView(const ArrayFilter<T> & array, Ns... ns)
      : array(array), sizes(std::move(ns)...) {}

  constexpr ArrayView(const ArrayView & array_view) = default;
  constexpr ArrayView(ArrayView && array_view) = default;

  constexpr ArrayView & operator=(const ArrayView & array_view) = default;
  constexpr ArrayView & operator=(ArrayView && array_view) = default;

  auto begin() const {
    return std::apply(
        [&](auto &&... ns) {
          return array.get().begin_reinterpret(
              std::forward<decltype(ns)>(ns)...);
        },
        sizes);
  }

  auto end() const {
    return std::apply(
        [&](auto &&... ns) {
          return array.get().end_reinterpret(std::forward<decltype(ns)>(ns)...);
        },
        sizes);
  }

  auto cbegin() const { return this->begin(); }
  auto cend() const { return this->end(); }

  constexpr auto size() const {
    return std::get<std::tuple_size<tuple>::value - 1>(sizes);
  }

  constexpr auto dims() const { return std::tuple_size<tuple>::value - 1; }

private:
  std::reference_wrapper<const ArrayFilter<T>> array;
  tuple sizes;
};
} // namespace detail

/* ------------------------------------------------------------------------ */
template <typename Array, typename... Ns,
          std::enable_if_t<aka::conjunction<
              std::is_integral<std::decay_t<Ns>>...>::value> * = nullptr>
auto make_view(Array && array, const Ns... ns) -> decltype(auto) {
  AKANTU_DEBUG_ASSERT((detail::product_all(ns...) != 0),
                      "You must specify non zero dimensions");
  // auto size = std::forward<decltype(array)>(array).size() *
  //             std::forward<decltype(array)>(array).getNbComponent() /
  //             detail::product_all(ns...);
  auto array_size = std::forward<decltype(array)>(array).size();
  auto nb_component = std::forward<decltype(array)>(array).getNbComponent();
  // detail::GetNbComponent<std::decay_t<Array>>::getNbComponent(
  //     std::forward<decltype(array)>(array));
  auto product_all = detail::product_all(ns...);
  auto size = array_size * nb_component / product_all;

  return detail::ArrayView<Array, std::common_type_t<size_t, Ns>...,
                           std::common_type_t<size_t, decltype(size)>>(
      std::forward<Array>(array), std::move(ns)..., size);
}

template <typename Array, typename... Ns>
auto make_const_view(const Array & array, const Ns... ns) -> decltype(auto) {
  return make_view(array, std::move(ns)...);
}

} // namespace akantu

#endif // AKANTU_ARRAY_VIEW_HH
