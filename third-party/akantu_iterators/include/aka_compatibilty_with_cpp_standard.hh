/**
 * @file   aka_compatibilty_with_cpp_standard.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Wed Jan 10 2018
 *
 * @brief  The content of this file is taken from the possible implementations
 * on
 * http://en.cppreference.com
 *
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * akantu-iterators is free  software: you can redistribute it and/or  modify it
 * under the terms  of the  GNU Lesser  General Public  License as published by
 * the Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * akantu-iterators is  distributed in the  hope that it  will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public
 * License  for more details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with akantu-iterators. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#if __cplusplus < 201703L
#include "aka_static_if.hh"
#endif
/* -------------------------------------------------------------------------- */
#include <iterator>
#include <tuple>
#include <type_traits>
#include <utility>

#if __cplusplus >= 201703L
#include <functional>
#include <algorithm>
#endif
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_AKA_COMPATIBILTY_WITH_CPP_STANDARD_HH
#define AKANTU_AKA_COMPATIBILTY_WITH_CPP_STANDARD_HH

namespace aka {

/* -------------------------------------------------------------------------- */
// Part taken from C++14
#if __cplusplus < 201402L
template <bool B, class T = void>
using enable_if_t = typename enable_if<B, T>::type;
#else
template <bool B, class T = void> using enable_if_t = std::enable_if_t<B, T>;
#endif

/* -------------------------------------------------------------------------- */
// Part taken from C++17
#if __cplusplus < 201703L

/* -------------------------------------------------------------------------- */
// bool_constant
template <bool B> using bool_constant = std::integral_constant<bool, B>;
namespace {
  template <bool B> constexpr bool bool_constant_v = bool_constant<B>::value;
}

/* -------------------------------------------------------------------------- */
// conjunction
template <class...> struct conjunction : std::true_type {};
template <class B1> struct conjunction<B1> : B1 {};
template <class B1, class... Bn>
struct conjunction<B1, Bn...>
    : std::conditional_t<bool(B1::value), conjunction<Bn...>, B1> {};

/* -------------------------------------------------------------------------- */
// disjunction
template <class...> struct disjunction : std::false_type {};
template <class B1> struct disjunction<B1> : B1 {};
template <class B1, class... Bn>
struct disjunction<B1, Bn...>
    : std::conditional_t<bool(B1::value), B1, disjunction<Bn...>> {};

/* -------------------------------------------------------------------------- */
// negations
template <class B> struct negation : bool_constant<!bool(B::value)> {};

/* -------------------------------------------------------------------------- */
// invoke
namespace detail {
  template <class T> struct is_reference_wrapper : std::false_type {};
  template <class U>
  struct is_reference_wrapper<std::reference_wrapper<U>> : std::true_type {};

  template <class T, class Type, class T1, class... Args>
  decltype(auto) INVOKE(Type T::*f, T1 && t1, Args &&... args) {
    static_if(std::is_member_function_pointer<decltype(f)>{})
        .then_([&](auto && f) {
          static_if(std::is_base_of<T, std::decay_t<T1>>{})
              .then_([&](auto && f) {
                return (std::forward<T1>(t1).*f)(std::forward<Args>(args)...);
              })
              .elseif(is_reference_wrapper<std::decay_t<T1>>{})([&](auto && f) {
                return (t1.get().*f)(std::forward<Args>(args)...);
              })
              .else_([&](auto && f) {
                return ((*std::forward<T1>(t1)).*
                        f)(std::forward<Args>(args)...);
              })(std::forward<decltype(f)>(f));
        })
        .else_([&](auto && f) {
          static_assert(std::is_member_object_pointer<decltype(f)>::value,
                        "f is not a member object");
          static_assert(sizeof...(args) == 0, "f takes arguments");
          static_if(std::is_base_of<T, std::decay_t<T1>>{})
              .then_([&](auto && f) { return std::forward<T1>(t1).*f; })
              .elseif(std::is_base_of<T, std::decay_t<T1>>{})(
                  [&](auto && f) { return t1.get().*f; })
              .else_([&](auto && f) { return (*std::forward<T1>(t1)).*f; })(
                  std::forward<decltype(f)>(f));
        })(std::forward<decltype(f)>(f));
  }

  template <class F, class... Args>
  decltype(auto) INVOKE(F && f, Args &&... args) {
    return std::forward<F>(f)(std::forward<Args>(args)...);
  }
} // namespace detail

template <class F, class... Args>
decltype(auto) invoke(F && f, Args &&... args) {
  return detail::INVOKE(std::forward<F>(f), std::forward<Args>(args)...);
}

/* -------------------------------------------------------------------------- */
// apply
namespace detail {
  template <class F, class Tuple, std::size_t... Is>
  constexpr decltype(auto) apply_impl(F && f, Tuple && t,
                                      std::index_sequence<Is...> /*unused*/) {
    return invoke(std::forward<F>(f), std::get<Is>(std::forward<Tuple>(t))...);
  }
} // namespace detail

/* -------------------------------------------------------------------------- */
template <class F, class Tuple>
constexpr decltype(auto) apply(F && f, Tuple && t) {
  return detail::apply_impl(
      std::forward<F>(f), std::forward<Tuple>(t),
      std::make_index_sequence<std::tuple_size<std::decay_t<Tuple>>::value>{});
}

/* -------------------------------------------------------------------------- */
// count_if
template <class InputIt, class UnaryPredicate>
auto count_if(InputIt first, InputIt last, UnaryPredicate p) ->
    typename std::iterator_traits<InputIt>::difference_type {
  typename std::iterator_traits<InputIt>::difference_type ret = 0;
  for (; first != last; ++first) {
    if (p(*first)) {
      ret++;
    }
  }
  return ret;
}

namespace detail {
  template <class T, class Tuple, std::size_t... I>
  constexpr T make_from_tuple_impl(Tuple && t, std::index_sequence<I...>) {
    static_assert(
        std::is_constructible<T, decltype(std::get<I>(
                                     std::declval<Tuple>()))...>::value,
        "make_from_tuple needs T to be constructible");
    return T(std::get<I>(std::forward<Tuple>(t))...);
  }
} // namespace detail

template <class T, class Tuple> constexpr T make_from_tuple(Tuple && t) {
  return detail::make_from_tuple_impl<T>(
      std::forward<Tuple>(t),
      std::make_index_sequence<
          std::tuple_size<std::remove_reference_t<Tuple>>::value>{});
}

#else
template <bool B> using bool_constant = std::bool_constant<B>;
template <bool B> constexpr bool bool_constant_v = std::bool_constant<B>::value;

template <class... Args> using conjunction = std::conjunction<Args...>;
template <class... Args> using disjunction = std::disjunction<Args...>;
template <class B> using negation = std::negation<B>;

template <class F, class Tuple>
constexpr decltype(auto) apply(F && f, Tuple && t) {
  return std::apply(std::forward<F>(f), std::forward<Tuple>(t));
}

template <class InputIt, class UnaryPredicate>
decltype(auto) count_if(InputIt first, InputIt last, UnaryPredicate p) {
  return std::count_if(std::forward<InputIt>(first),
                       std::forward<InputIt>(last),
                       std::forward<UnaryPredicate>(p));
}

template <class T, class Tuple>
constexpr auto make_from_tuple(Tuple && t) -> T {
  return std::make_from_tuple<T>(std::forward<Tuple>(t));
}
#endif

template <typename cat1, typename cat2>
using is_iterator_category_at_least =
    std::is_same<std::common_type_t<cat1, cat2>, cat2>;

template <typename T> struct size_type {
  using type = typename std::decay_t<T>::size_type;
};

template <typename T> using size_type_t = typename size_type<T>::type;

} // namespace aka

#endif /* AKANTU_AKA_COMPATIBILTY_WITH_CPP_STANDARD_HH */
