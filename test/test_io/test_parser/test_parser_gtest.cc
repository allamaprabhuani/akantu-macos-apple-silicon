#include "algebraic_parser_x3.hh"
#include "parser.hh"
#include <cmath>
#include <gtest/gtest.h>
#include <limits>
#include <memory>

using namespace akantu;

constexpr Real eps = std::numeric_limits<Real>::epsilon();

class AlgebraicParser : public testing::Test {
protected:
  void SetUp() override { parser = std::make_unique<Parser>(); }

  void TearDown() override { parser.reset(nullptr); }

  std::unique_ptr<Parser> parser;
};

/* -------------------------------------------------------------------------- */
TEST_F(AlgebraicParser, Number) {
  auto number = parser::algebraic::parse_real("3.", *parser);
  EXPECT_NEAR(number, 3., eps);

  number = parser::algebraic::parse_real("10", *parser);
  EXPECT_NEAR(number, 10., eps);

  number = parser::algebraic::parse_real("-100.0", *parser);
  EXPECT_NEAR(number, -100., eps);

  number = parser::algebraic::parse_real("+1e10", *parser);
  EXPECT_NEAR(number, 1e10, eps);
}

/* -------------------------------------------------------------------------- */
TEST_F(AlgebraicParser, Expr) {
  auto number = parser::algebraic::parse_real("3. + 10.", *parser);
  EXPECT_NEAR(number, 13., eps);

  number = parser::algebraic::parse_real("10 - 100", *parser);
  EXPECT_NEAR(number, -90., eps);

  number = parser::algebraic::parse_real("-1 - -1", *parser);
  EXPECT_NEAR(number, 0., eps);

  number = parser::algebraic::parse_real("1e10 + 4e5 - 1e9 + 1e10", *parser);
  EXPECT_NEAR(number, (1e10 + 4e5 - 1e9 + 1e10), eps);
}

/* -------------------------------------------------------------------------- */
TEST_F(AlgebraicParser, Term) {
  auto number = parser::algebraic::parse_real("3. * 10.", *parser);
  EXPECT_NEAR(number, 30., eps);

  number = parser::algebraic::parse_real("10 * 100", *parser);
  EXPECT_NEAR(number, 1000, eps);

  number = parser::algebraic::parse_real("-1 / -1", *parser);
  EXPECT_NEAR(number, 1., eps);

  number = parser::algebraic::parse_real("1e10 / 4e5", *parser);
  EXPECT_NEAR(number, (1e10 / 4e5), eps);
}

/* -------------------------------------------------------------------------- */
TEST_F(AlgebraicParser, Factor) {
  auto number = parser::algebraic::parse_real("3. ** 10.", *parser);
  EXPECT_NEAR(number, std::pow(3, 10), eps);
}

/* -------------------------------------------------------------------------- */
TEST_F(AlgebraicParser, Func) {
  /* ------------------------------------------------------------------------ */
  // Binary Functrions
  auto number = parser::algebraic::parse_real("pow(3., 10.)", *parser);
  EXPECT_NEAR(number, std::pow(3., 10.), eps);

  number = parser::algebraic::parse_real("min(3., 10.)", *parser);
  EXPECT_NEAR(number, std::min(3., 10.), eps);

  number = parser::algebraic::parse_real("max(3., 10.)", *parser);
  EXPECT_NEAR(number, std::max(3., 10.), eps);

  number = parser::algebraic::parse_real("atan2(3.1415, 2.)", *parser);
  EXPECT_NEAR(number, std::atan2(3.1415, 2.), eps);

  number = parser::algebraic::parse_real("fmod(10., 3.)", *parser);
  EXPECT_NEAR(number, std::fmod(10., 3.), eps);

  number = parser::algebraic::parse_real("hypot(3., 10.)", *parser);
  EXPECT_NEAR(number, std::hypot(3., 10.), eps);

  /* ------------------------------------------------------------------------ */
  // Unary Functrions
  number = parser::algebraic::parse_real("abs(3.)", *parser);
  EXPECT_NEAR(number, std::abs(3.), eps);

  number = parser::algebraic::parse_real("acos(.5)", *parser);
  EXPECT_NEAR(number, std::acos(0.5), eps);

  number = parser::algebraic::parse_real("asin(.5)", *parser);
  EXPECT_NEAR(number, std::asin(.5), eps);

  number = parser::algebraic::parse_real("atan(.5)", *parser);
  EXPECT_NEAR(number, std::atan(.5), eps);

  number = parser::algebraic::parse_real("ceil(3.)", *parser);
  EXPECT_NEAR(number, std::ceil(3.), eps);

  number = parser::algebraic::parse_real("cos(3.)", *parser);
  EXPECT_NEAR(number, std::cos(3.), eps);

  number = parser::algebraic::parse_real("cosh(3.)", *parser);
  EXPECT_NEAR(number, std::cosh(3.), eps);

  number = parser::algebraic::parse_real("exp(3.)", *parser);
  EXPECT_NEAR(number, std::exp(3.), eps);

  number = parser::algebraic::parse_real("floor(3.)", *parser);
  EXPECT_NEAR(number, std::floor(3.), eps);

  number = parser::algebraic::parse_real("log10(3.)", *parser);
  EXPECT_NEAR(number, std::log10(3.), eps);

  number = parser::algebraic::parse_real("log(3.)", *parser);
  EXPECT_NEAR(number, std::log(3.), eps);

  number = parser::algebraic::parse_real("sin(3.)", *parser);
  EXPECT_NEAR(number, std::sin(3.), eps);

  number = parser::algebraic::parse_real("sinh(3.)", *parser);
  EXPECT_NEAR(number, std::sinh(3.), eps);

  number = parser::algebraic::parse_real("sqrt(3.)", *parser);
  EXPECT_NEAR(number, std::sqrt(3.), eps);

  number = parser::algebraic::parse_real("tan(3.)", *parser);
  EXPECT_NEAR(number, std::tan(3.), eps);

  number = parser::algebraic::parse_real("tanh(3.)", *parser);
  EXPECT_NEAR(number, std::tanh(3.), eps);

  number = parser::algebraic::parse_real("acosh(3.)", *parser);
  EXPECT_NEAR(number, std::acosh(3.), eps);

  number = parser::algebraic::parse_real("asinh(3.)", *parser);
  EXPECT_NEAR(number, std::asinh(3.), eps);

  number = parser::algebraic::parse_real("atanh(.5)", *parser);
  EXPECT_NEAR(number, std::atanh(.5), eps);

  number = parser::algebraic::parse_real("exp2(3.)", *parser);
  EXPECT_NEAR(number, std::exp2(3.), eps);

  number = parser::algebraic::parse_real("expm1(3.)", *parser);
  EXPECT_NEAR(number, std::expm1(3.), eps);

  number = parser::algebraic::parse_real("log1p(3.)", *parser);
  EXPECT_NEAR(number, std::log1p(3.), eps);

  number = parser::algebraic::parse_real("log2(3.)", *parser);
  EXPECT_NEAR(number, std::log2(3.), eps);

  number = parser::algebraic::parse_real("erf(3.)", *parser);
  EXPECT_NEAR(number, std::erf(3.), eps);

  number = parser::algebraic::parse_real("erfc(3.)", *parser);
  EXPECT_NEAR(number, std::erfc(3.), eps);

  number = parser::algebraic::parse_real("lgamma(3.)", *parser);
  EXPECT_NEAR(number, std::lgamma(3.), eps);

  number = parser::algebraic::parse_real("tgamma(3.)", *parser);
  EXPECT_NEAR(number, std::tgamma(3.), eps);

  number = parser::algebraic::parse_real("trunc(3.14)", *parser);
  EXPECT_NEAR(number, std::trunc(3.), eps);

  number = parser::algebraic::parse_real("round(3.14)", *parser);
  EXPECT_NEAR(number, std::round(3.14), eps);
}

/* -------------------------------------------------------------------------- */
TEST_F(AlgebraicParser, Constant) {
  auto number = parser::algebraic::parse_real("e", *parser);
  EXPECT_NEAR(number, 2.7182818284590452354, eps);

  number = parser::algebraic::parse_real("pi", *parser);
  EXPECT_NEAR(number, 3.14159265358979323846, eps);
}

/* -------------------------------------------------------------------------- */
TEST_F(AlgebraicParser, Complex) {
  constexpr auto pi = 3.14159265358979323846;
  auto number = parser::algebraic::parse_real(
      "10. * (cos(3 * pi) + sin(5 * pi / 2)) ** 2", *parser);
  EXPECT_NEAR(number,
              10. * std::pow(std::cos(3 * pi) + std::sin(5 * pi / 2), 2), eps);
}
