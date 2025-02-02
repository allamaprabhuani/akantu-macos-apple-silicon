/**
 * Copyright (©) 2013-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "algebraic_parser_x3.hh"
/* -------------------------------------------------------------------------- */
// #define BOOST_SPIRIT_X3_DEBUG
#include <boost/config/warning_disable.hpp>
#include <boost/fusion/include/at_c.hpp>
#include <boost/spirit/home/x3.hpp>
#if __cplusplus >= 202000L
#include <numbers>
#endif
/* -------------------------------------------------------------------------- */

namespace akantu {
namespace parser {
namespace algebraic {

namespace x3 = boost::spirit::x3;
namespace fusion = boost::fusion;

const x3::rule<struct start_class, Real> start{"start"};
const x3::rule<struct expr_class, Real> expr{"expr"};
const x3::rule<struct term_class, Real> term{"term"};
const x3::rule<struct factor_class, Real> factor{"factor"};
const x3::rule<struct number_class, Real> number{"number"};
const x3::rule<struct variable_class, Real> variable{"variable"};
const x3::rule<struct func_class, Real> function{"function"};
const x3::rule<struct key_class, std::string> key{"key"};

struct binary_function_ : x3::symbols<std::function<Real(Real, Real)>> {
  binary_function_() {
    add("pow", [](Real x, Real y) -> Real { return std::pow(x, y); })(
        "min", [](Real x, Real y) -> Real { return std::min(x, y); })(
        "max", [](Real x, Real y) -> Real { return std::max(x, y); })(
        "atan2", [](Real x, Real y) -> Real { return std::atan2(x, y); })(
        "fmod", [](Real x, Real y) -> Real { return std::fmod(x, y); })(
        "hypot", [](Real x, Real y) -> Real { return std::hypot(x, y); });
  }
} binary_function;

struct unary_function_ : x3::symbols<std::function<Real(Real)>> {
  unary_function_() {
    add("abs", [](Real x) -> Real { return std::abs(x); })(
        "acos", [](Real x) -> Real { return std::acos(x); })(
        "asin", [](Real x) -> Real { return std::asin(x); })(
        "atan", [](Real x) -> Real { return std::atan(x); })(
        "ceil", [](Real x) -> Real { return std::ceil(x); })(
        "cos", [](Real x) -> Real { return std::cos(x); })(
        "cosh", [](Real x) -> Real { return std::cosh(x); })(
        "exp", [](Real x) -> Real { return std::exp(x); })(
        "floor", [](Real x) -> Real { return std::floor(x); })(
        "log10", [](Real x) -> Real { return std::log10(x); })(
        "log", [](Real x) -> Real { return std::log(x); })(
        "sin", [](Real x) -> Real { return std::sin(x); })(
        "sinh", [](Real x) -> Real { return std::sinh(x); })(
        "sqrt", [](Real x) -> Real { return std::sqrt(x); })(
        "tan", [](Real x) -> Real { return std::tan(x); })(
        "tanh", [](Real x) -> Real { return std::tanh(x); })(
        "acosh", [](Real x) -> Real { return std::acosh(x); })(
        "asinh", [](Real x) -> Real { return std::asinh(x); })(
        "atanh", [](Real x) -> Real { return std::atanh(x); })(
        "exp2", [](Real x) -> Real { return std::exp2(x); })(
        "expm1", [](Real x) -> Real { return std::expm1(x); })(
        "log1p", [](Real x) -> Real { return std::log1p(x); })(
        "log2", [](Real x) -> Real { return std::log2(x); })(
        "erf", [](Real x) -> Real { return std::erf(x); })(
        "erfc", [](Real x) -> Real { return std::erfc(x); })(
        "lgamma", [](Real x) -> Real { return std::lgamma(x); })(
        "tgamma", [](Real x) -> Real { return std::tgamma(x); })(
        "trunc", [](Real x) -> Real { return std::trunc(x); })(
        "round", [](Real x) -> Real { return std::round(x); })
        //      ("crbt"  , &std::crbt  )
        ;
  }
} unary_function;

struct constant_ : x3::symbols<Real> {
  constant_() {
    add
#if __cplusplus >= 202000L
        ("pi", std::numbers::pi_v<Real>)("e", std::numbers::e_v<Real>)
#else
        ("pi", 3.14159265358979323846)("e", 2.7182818284590452354)
#endif
            ;
  }
} constant;

auto assign = [](auto & ctx) { _val(ctx) = _attr(ctx); };
auto negate = [](auto & ctx) { _val(ctx) = -_attr(ctx); };
auto add = [](auto & ctx) { _val(ctx) = _val(ctx) + _attr(ctx); };
auto sub = [](auto & ctx) { _val(ctx) = _val(ctx) - _attr(ctx); };
auto mul = [](auto & ctx) { _val(ctx) = _val(ctx) * _attr(ctx); };
auto div = [](auto & ctx) { _val(ctx) = _val(ctx) / _attr(ctx); };
auto exponatiate = [](auto & ctx) {
  _val(ctx) = std::pow(_val(ctx), _attr(ctx));
};
auto call1 = [](auto & ctx) {
  auto & attr = _attr(ctx);
  auto & op = fusion::at_c<0>(attr);
  auto & x = fusion::at_c<1>(attr);
  _val(ctx) = op(x);
};
auto call2 = [](auto & ctx) {
  auto & attr = _attr(ctx);
  auto & op = fusion::at_c<0>(attr);
  auto & x = fusion::at_c<1>(attr);
  auto & y = fusion::at_c<2>(attr);
  _val(ctx) = op(x, y);
};

auto eval = [](auto & ctx) {
  auto && section = x3::get<const ParserSection &>(ctx);
  _val(ctx) = section.getParameter(_attr(ctx), _ppsc_current_and_parent_scope);
};

/* clang-format off */
const auto function_def
  =   (binary_function
           > '(' >> expr
           > ',' >> expr
           > ')')        [call2]
      | (unary_function
           > '('
           > expr
           > ')')        [call1]
      ;

const auto variable_def
  =   key [eval]
      ;

const auto key_def
  =   x3::no_skip[x3::char_("a-zA-Z_") >> *x3::char_("a-zA-Z_0-9")]
      ;

const auto number_def
  =   x3::double_            [assign]
      |   ('-' > number      [negate])
      |   ('+' > number      [assign])
      |   function           [assign]
      |   constant           [assign]
      |   variable           [assign]
      |   ('(' > expr > ')') [assign]
      ;

const auto factor_def
  =   number             [assign]
      >> *("**" > number [exponatiate])
      ;

const auto term_def
  =   factor                [assign]
      >> *( ('*' > factor   [mul])
            | ('/' > factor [div])
            )
      ;

const auto expr_def
  =   term                [assign]
      >> *( ('+' > term   [add])
            | ('-' > term [sub])
            )
      ;

const auto start_def
  =   expr
      ;
/* clang-format on */

BOOST_SPIRIT_DEFINE(start, expr, term, factor, number, variable, function, key);

Real parse_real(const std::string & value, const ParserSection & section) {
  //  using boost::ascii::space;
  Real res;
  auto && it = value.begin();
  auto ret = x3::phrase_parse(it, value.end(),
                              x3::with<const ParserSection &>(section)[start],
                              x3::ascii::space, res);

  if (not ret or it != value.end()) {
    AKANTU_EXCEPTION("Could not parse the expression " << value);
  }

  return res;
}

} // namespace algebraic
} // namespace parser
} // namespace akantu
