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
#include "aka_config.hh"
#include "aka_tuple_tools.hh"
#if defined(AKANTU_CAN_COMPILE_CONSTEXPR_MAP)
#include "aka_constexpr_map.hh"
#endif
/* -------------------------------------------------------------------------- */
#include <tuple>
#include <unordered_map>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_AKA_ELEMENT_CLASSES_INFO_INLINE_IMPL_HH_
#define AKANTU_AKA_ELEMENT_CLASSES_INFO_INLINE_IMPL_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
// BOOST PART: TOUCH ONLY IF YOU KNOW WHAT YOU ARE DOING
#define AKANTU_BOOST_CASE_MACRO(r, macro, _type)                               \
  case _type: {                                                                \
    macro(_type);                                                              \
    break;                                                                     \
  }

#define AKANTU_BOOST_LIST_SWITCH(macro1, list1, var)                           \
  do {                                                                         \
    switch (var) {                                                             \
      BOOST_PP_SEQ_FOR_EACH(AKANTU_BOOST_CASE_MACRO, macro1, list1)            \
    default: {                                                                 \
      AKANTU_ERROR("Type (" << var << ") not handled by this function");       \
    }                                                                          \
    }                                                                          \
  } while (0)

#define AKANTU_BOOST_LIST_SWITCH_NO_DEFAULT(macro1, list1, var)                \
  do {                                                                         \
    switch (var) {                                                             \
      BOOST_PP_SEQ_FOR_EACH(AKANTU_BOOST_CASE_MACRO, macro1, list1)            \
    case _not_defined:                                                         \
      break;                                                                   \
    case _max_element_type:                                                    \
      break;                                                                   \
    }                                                                          \
  } while (0)

#define AKANTU_BOOST_LIST_SWITCH_CONSTEXPR(macro1, list1, var)                 \
  do {                                                                         \
    switch (var) {                                                             \
      BOOST_PP_SEQ_FOR_EACH(AKANTU_BOOST_CASE_MACRO, macro1, list1)            \
    default: {                                                                 \
      macro1(_not_defined);                                                    \
    }                                                                          \
    }                                                                          \
  } while (0)

#define AKANTU_BOOST_LIST_SWITCH_WITH_DEFAULT_MACRO(macro1, def_macro, list1,  \
                                                    var)                       \
  do {                                                                         \
    switch (var) {                                                             \
      BOOST_PP_SEQ_FOR_EACH(AKANTU_BOOST_CASE_MACRO, macro1, list1)            \
    default: {                                                                 \
      def_macro(var);                                                          \
    }                                                                          \
    }                                                                          \
  } while (0)

#define AKANTU_BOOST_ELEMENT_SWITCH(macro1, list1)                             \
  AKANTU_BOOST_LIST_SWITCH(macro1, list1, type)

#define AKANTU_BOOST_ELEMENT_SWITCH_NO_DEFAULT(macro1, list1)                  \
  AKANTU_BOOST_LIST_SWITCH_NO_DEFAULT(macro1, list1, type)

#define AKANTU_BOOST_ELEMENT_SWITCH_CONSTEXPR(macro1, list1)                   \
  AKANTU_BOOST_LIST_SWITCH_CONSTEXPR(macro1, list1, type)

#define AKANTU_BOOST_ELEMENT_SWITCH_WITH_DEFAULT_MACRO(macro1, def_macro,      \
                                                       list1)                  \
  AKANTU_BOOST_LIST_SWITCH_WITH_DEFAULT_MACRO(macro1, def_macro, list1, type)

#define AKANTU_BOOST_ALL_ELEMENT_SWITCH(macro)                                 \
  AKANTU_BOOST_ELEMENT_SWITCH(macro, AKANTU_ALL_ELEMENT_TYPE)

#define AKANTU_BOOST_ALL_ELEMENT_SWITCH_NO_DEFAULT(macro)                      \
  AKANTU_BOOST_ELEMENT_SWITCH_NO_DEFAULT(macro, AKANTU_ALL_ELEMENT_TYPE)

#define AKANTU_BOOST_ALL_ELEMENT_SWITCH_CONSTEXPR(macro)                       \
  AKANTU_BOOST_ELEMENT_SWITCH_CONSTEXPR(macro, AKANTU_ALL_ELEMENT_TYPE)

#define AKANTU_BOOST_LIST_MACRO(r, macro, type) macro(type)

#define AKANTU_BOOST_APPLY_ON_LIST(macro, list)                                \
  BOOST_PP_SEQ_FOR_EACH(AKANTU_BOOST_LIST_MACRO, macro, list)

#define AKANTU_BOOST_ALL_ELEMENT_LIST(macro)                                   \
  AKANTU_BOOST_APPLY_ON_LIST(macro, AKANTU_ALL_ELEMENT_TYPE)

#define AKANTU_GET_ELEMENT_LIST(kind) AKANTU##kind##_ELEMENT_TYPE

#define AKANTU_BOOST_KIND_ELEMENT_SWITCH(macro, kind)                          \
  AKANTU_BOOST_ELEMENT_SWITCH(macro, AKANTU_GET_ELEMENT_LIST(kind))

// BOOST_PP_SEQ_TO_LIST does not exists in Boost < 1.49
#define AKANTU_GENERATE_KIND_LIST(seq)                                         \
  BOOST_PP_TUPLE_TO_LIST(BOOST_PP_SEQ_SIZE(seq), BOOST_PP_SEQ_TO_TUPLE(seq))

#define AKANTU_ELEMENT_KIND_BOOST_LIST                                         \
  AKANTU_GENERATE_KIND_LIST(AKANTU_ELEMENT_KIND)

#define AKANTU_BOOST_ALL_KIND_LIST(macro, list)                                \
  BOOST_PP_LIST_FOR_EACH(AKANTU_BOOST_LIST_MACRO, macro, list)

#define AKANTU_BOOST_ALL_KIND(macro)                                           \
  AKANTU_BOOST_ALL_KIND_LIST(macro, AKANTU_ELEMENT_KIND_BOOST_LIST)

#define AKANTU_BOOST_ALL_KIND_SWITCH(macro)                                    \
  AKANTU_BOOST_LIST_SWITCH(macro, AKANTU_ELEMENT_KIND, kind)

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

AKANTU_ENUM_OUTPUT_STREAM(
    ElementType, AKANTU_ALL_ELEMENT_TYPE(_not_defined)(_max_element_type))
AKANTU_ENUM_INPUT_STREAM(ElementType, AKANTU_ALL_ELEMENT_TYPE)

AKANTU_ENUM_OUTPUT_STREAM(InterpolationType, AKANTU_INTERPOLATION_TYPES)
AKANTU_ENUM_INPUT_STREAM(InterpolationType, AKANTU_INTERPOLATION_TYPES)

AKANTU_ENUM_OUTPUT_STREAM(ElementKind, AKANTU_ELEMENT_KIND)
AKANTU_ENUM_INPUT_STREAM(ElementKind, AKANTU_ELEMENT_KIND)

template <::akantu::ElementType t>
using element_type_t = std::integral_constant<::akantu::ElementType, t>;

#define OP_CAT(s, data, elem) BOOST_PP_CAT(_element_type, elem)

// creating a type instead of a using helps to debug
#define AKANTU_DECLARE_ELEMENT_TYPE_STRUCT(r, data, elem)                      \
  struct BOOST_PP_CAT(_element_type, elem)                                     \
      : public element_type_t<::akantu::elem> {};

BOOST_PP_SEQ_FOR_EACH(AKANTU_DECLARE_ELEMENT_TYPE_STRUCT, _,
                      AKANTU_ALL_ELEMENT_TYPE)

#undef AKANTU_DECLARE_ELEMENT_TYPE_STRUCT

template <ElementKind kind> struct ElementTypes;

template <> struct ElementTypes<_ek_regular> {
  using type = std::tuple<BOOST_PP_SEQ_ENUM(
      BOOST_PP_SEQ_TRANSFORM(OP_CAT, _, AKANTU_ek_regular_ELEMENT_TYPE))>;
};

#if defined(AKANTU_COHESIVE_ELEMENT)
template <> struct ElementTypes<_ek_cohesive> {
  using type = std::tuple<BOOST_PP_SEQ_ENUM(
      BOOST_PP_SEQ_TRANSFORM(OP_CAT, _, AKANTU_ek_cohesive_ELEMENT_TYPE))>;
};
#else
template <> struct ElementTypes<_ek_cohesive> { using type = std::tuple<>; };
#endif

#if defined(AKANTU_STRUCTURAL_MECHANICS)
template <> struct ElementTypes<_ek_structural> {
  using type = std::tuple<BOOST_PP_SEQ_ENUM(
      BOOST_PP_SEQ_TRANSFORM(OP_CAT, _, AKANTU_ek_structural_ELEMENT_TYPE))>;
};
#else
template <> struct ElementTypes<_ek_structural> { using type = std::tuple<>; };
#endif

#undef OP_CAT

template <ElementKind kind>
using ElementTypes_t = typename ElementTypes<kind>::type;

#define OP_CAT(s, data, elem) ElementTypes_t<elem>
using AllElementTypes = tuple::cat_t<BOOST_PP_SEQ_ENUM(
    BOOST_PP_SEQ_TRANSFORM(OP_CAT, _, AKANTU_ELEMENT_KIND))>;
#undef OP_CAT

namespace details {

#if !defined(AKANTU_CAN_COMPILE_CONSTEXPR_MAP)
  template <class Function, class DynamicType, class DefaultFunction,
            std::size_t... Is>
  constexpr decltype(auto)
  static_switch_dispatch(const AllSpatialDimensions &, Function && function,
                         const DynamicType & type,
                         DefaultFunction && default_function,
                         std::index_sequence<Is...> /*is*/) {
    switch (type) {
    case 1: {
      return function(dim_1_t{});
    }
    case 2: {
      return function(dim_2_t{});
    }
    case 3: {
      return function(dim_3_t{});
    }
    default:
      return default_function(type);
    }
  }

#define AKANTU_STATIC_SWITCH_DISPATCH_FUNCTION(type)                           \
  return function(element_type_t<type>{})
#define AKANTU_STATIC_SWITCH_DISPATCH_DEFAULT(type)                            \
  return default_function(type)
  template <class Function, class DynamicType, class DefaultFunction,
            std::size_t... Is>
  constexpr decltype(auto) static_switch_dispatch(
      const AllElementTypes &, Function && function, const DynamicType & type,
      DefaultFunction && default_function, std::index_sequence<Is...> /*is*/) {

    AKANTU_BOOST_ELEMENT_SWITCH_WITH_DEFAULT_MACRO(
        AKANTU_STATIC_SWITCH_DISPATCH_FUNCTION,
        AKANTU_STATIC_SWITCH_DISPATCH_DEFAULT, AKANTU_ALL_ELEMENT_TYPE);
  }

  template <class Function, class DynamicType, class DefaultFunction,
            std::size_t... Is>
  constexpr decltype(auto)
  static_switch_dispatch(const ElementTypes_t<_ek_regular> &,
                         Function && function, const DynamicType & type,
                         DefaultFunction && default_function,
                         std::index_sequence<Is...> /*is*/) {
    AKANTU_BOOST_ELEMENT_SWITCH_WITH_DEFAULT_MACRO(
        AKANTU_STATIC_SWITCH_DISPATCH_FUNCTION,
        AKANTU_STATIC_SWITCH_DISPATCH_DEFAULT, AKANTU_ek_regular_ELEMENT_TYPE);
  }

#if defined(AKANTU_STRUCTURAL_MECHANICS)
  template <class Function, class DynamicType, class DefaultFunction,
            std::size_t... Is>
  constexpr decltype(auto)
  static_switch_dispatch(const ElementTypes_t<_ek_structural> &,
                         Function && function, const DynamicType & type,
                         DefaultFunction && default_function,
                         std::index_sequence<Is...> /*is*/) {
    AKANTU_BOOST_ELEMENT_SWITCH_WITH_DEFAULT_MACRO(
        AKANTU_STATIC_SWITCH_DISPATCH_FUNCTION,
        AKANTU_STATIC_SWITCH_DISPATCH_DEFAULT,
        AKANTU_ek_structural_ELEMENT_TYPE);
  }
#endif

#if defined(AKANTU_COHESIVE_ELEMENT)
  template <class Function, class DynamicType, class DefaultFunction,
            std::size_t... Is>
  constexpr decltype(auto)
  static_switch_dispatch(const ElementTypes_t<_ek_cohesive> &,
                         Function && function, const DynamicType & type,
                         DefaultFunction && default_function,
                         std::index_sequence<Is...> /*is*/) {
    AKANTU_BOOST_ELEMENT_SWITCH_WITH_DEFAULT_MACRO(
        AKANTU_STATIC_SWITCH_DISPATCH_FUNCTION,
        AKANTU_STATIC_SWITCH_DISPATCH_DEFAULT, AKANTU_ek_cohesive_ELEMENT_TYPE);
  }

  template <class Function, class DynamicType, class DefaultFunction,
            std::size_t... Is>
  constexpr decltype(auto)
  static_switch_dispatch(const tuple::cat_t<ElementTypes_t<_ek_regular>,
                                            ElementTypes_t<_ek_cohesive>> &,
                         Function && function, const DynamicType & type,
                         DefaultFunction && default_function,
                         std::index_sequence<Is...> /*is*/) {
#define AKANTU_REGULAR_AND_COHESIVE_ELEMENT_TYPE                               \
  AKANTU_ek_regular_ELEMENT_TYPE AKANTU_ek_cohesive_ELEMENT_TYPE

    AKANTU_BOOST_ELEMENT_SWITCH_WITH_DEFAULT_MACRO(
        AKANTU_STATIC_SWITCH_DISPATCH_FUNCTION,
        AKANTU_STATIC_SWITCH_DISPATCH_DEFAULT,
        AKANTU_REGULAR_AND_COHESIVE_ELEMENT_TYPE);

#undef AKANTU_REGULAR_AND_COHESIVE_ELEMENT_TYPE
  }

#endif

#undef AKANTU_STATIC_SWITCH_DISPATCH_FUNCTION
#undef AKANTU_STATIC_SWITCH_DISPATCH_DEFAULT

#endif
} // namespace details

template <class Tuple, class Function, class DynamicType, class DefaultFunction>
[[gnu::visibility("hidden")]] constexpr decltype(auto)
tuple_dispatch_with_default(Function && function,
                            DefaultFunction && default_function,
                            const DynamicType & type) {
  return details::static_switch_dispatch(
      Tuple{}, std::forward<Function>(function), type,
      std::forward<DefaultFunction>(default_function),
      std::make_index_sequence<std::tuple_size<Tuple>::value>{});
}

template <class Tuple, class Function, class DynamicType>
[[gnu::visibility("hidden")]] constexpr decltype(auto)
tuple_dispatch(Function && function, const DynamicType & type) {
  return details::static_switch_dispatch(
      Tuple{}, std::forward<Function>(function), type,
      [](auto && type) -> decltype(function(std::tuple_element_t<0, Tuple>{})) {
        throw std::range_error(
            "Unknown " + debug::demangle(typeid(DynamicType).name()) +
            " type in dispatch: " + std::to_string(Int(type)));
      },
      std::make_index_sequence<std::tuple_size<Tuple>::value>{});
}

} // namespace akantu

#endif /* AKANTU_AKA_ELEMENT_CLASSES_INFO_INLINE_IMPL_HH_ */
