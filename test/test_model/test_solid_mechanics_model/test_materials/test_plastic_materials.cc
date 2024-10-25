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
#include "material_linear_isotropic_hardening.hh"
#include "solid_mechanics_model.hh"
#include "test_material_fixtures.hh"
#include <gtest/gtest.h>
#include <type_traits>
/* -------------------------------------------------------------------------- */

using namespace akantu;

using mat_types = ::testing::Types<
    // Traits<MaterialLinearIsotropicHardening, 1>,
    // Traits<MaterialLinearIsotropicHardening, 2>,
    Traits<MaterialLinearIsotropicHardening, 3>>;

/* -------------------------------------------------------------------------- */

template <>
void FriendMaterial<MaterialLinearIsotropicHardening<3>>::testComputeStress() {

  Real E = 1.;
  // Real nu = .3;
  Real nu = 0.;
  Real rho = 1.;
  Real sigma_0 = 1.;
  Real h = 0.;
  Real bulk_modulus_K = E / 3. / (1 - 2. * nu);
  Real shear_modulus_mu = 0.5 * E / (1 + nu);

  set("E", E);
  set("nu", nu);
  set("rho", rho);
  set("sigma_y", sigma_0);
  set("h", h);

  auto rotation_matrix = getRandomRotation();

  Real max_strain = 10.;
  Real strain_steps = 100;
  Real dt = max_strain / strain_steps;
  std::vector<double> steps(strain_steps);
  std::iota(steps.begin(), steps.end(), 0.);

  Matrix<Real, 3, 3> previous_grad_u_rot = Matrix<Real, 3, 3>::Zero();
  Matrix<Real, 3, 3> previous_sigma = Matrix<Real, 3, 3>::Zero();
  Matrix<Real, 3, 3> previous_sigma_rot = Matrix<Real, 3, 3>::Zero();
  Matrix<Real, 3, 3> inelastic_strain_rot = Matrix<Real, 3, 3>::Zero();
  // Matrix<Real, 3, 3> inelastic_strain = Matrix<Real, 3, 3>::Zero();
  Matrix<Real, 3, 3> previous_inelastic_strain = Matrix<Real, 3, 3>::Zero();
  Matrix<Real, 3, 3> previous_inelastic_strain_rot = Matrix<Real, 3, 3>::Zero();
  Matrix<Real, 3, 3> sigma_rot = Matrix<Real, 3, 3>::Zero();

  Real iso_hardening = 0.;
  Real previous_iso_hardening = 0.;

  // hydrostatic loading (should not plastify)
  for (auto && i : steps) {
    auto t = i * dt;

    Matrix<Real, 3, 3> grad_u = this->getHydrostaticStrain(t);
    Matrix<Real, 3, 3> grad_u_rot =
        this->applyRotation(grad_u, rotation_matrix);

    this->computeStressOnQuad(make_named_tuple(
        "grad_u"_n = grad_u_rot, "sigma"_n = sigma_rot,
        "previous_sigma"_n = previous_sigma_rot,
        "previous_grad_u"_n = previous_grad_u_rot,
        "inelastic_strain"_n = inelastic_strain_rot,
        "previous_inelastic_strain"_n = previous_inelastic_strain_rot,
        "iso_hardening"_n = iso_hardening,
        "previous_iso_hardening"_n = previous_iso_hardening));

    Matrix<Real, 3, 3> sigma =
        this->reverseRotation(sigma_rot, rotation_matrix);

    Matrix<Real, 3, 3> sigma_expected =
        t * 3. * bulk_modulus_K * Matrix<Real, 3, 3>::Identity();

    Real stress_error = (sigma - sigma_expected).lpNorm<Eigen::Infinity>();

    ASSERT_NEAR(stress_error, 0., 1e-13);
    ASSERT_NEAR(inelastic_strain_rot.lpNorm<Eigen::Infinity>(), 0., 1e-13);

    previous_grad_u_rot = grad_u_rot;
    previous_sigma_rot = sigma_rot;
    previous_inelastic_strain_rot = inelastic_strain_rot;
    previous_iso_hardening = iso_hardening;
  }

  // deviatoric loading (should plastify)
  // stress at onset of plastication
  Real beta = sqrt(42);
  Real t_P = sigma_0 / 2. / shear_modulus_mu / beta;
  // Matrix<Real, 3, 3> sigma_P = sigma_0 / beta *
  // this->getDeviatoricStrain(1.);

  for (auto && i : steps) {

    auto t = i * dt;
    auto grad_u = this->getDeviatoricStrain(t);
    auto grad_u_rot = this->applyRotation(grad_u, rotation_matrix);
    Real iso_hardening{0.};
    Real previous_iso_hardening{0.};

    this->computeStressOnQuad(make_named_tuple(
        "grad_u"_n = grad_u_rot, "sigma"_n = sigma_rot,
        "previous_grad_u"_n = previous_grad_u_rot,
        "previous_sigma"_n = previous_sigma_rot,
        "inelastic_strain"_n = inelastic_strain_rot,
        "previous_inelastic_strain"_n = previous_inelastic_strain_rot,
        "iso_hardening"_n = iso_hardening,
        "previous_iso_hardening"_n = previous_iso_hardening));

    auto sigma = this->reverseRotation(sigma_rot, rotation_matrix);
    auto inelastic_strain =
        this->reverseRotation(inelastic_strain_rot, rotation_matrix);

    if (t < t_P) {
      Matrix<Real, 3, 3> sigma_expected =
          shear_modulus_mu * (grad_u + grad_u.transpose());

      Real stress_error = (sigma - sigma_expected).lpNorm<Eigen::Infinity>();
      ASSERT_NEAR(stress_error, 0., 1e-13);
      ASSERT_NEAR(inelastic_strain_rot.lpNorm<Eigen::Infinity>(), 0., 1e-13);
    } else if (t > t_P + dt) {
      // skip the transition from non plastic to plastic

      auto delta_lambda_expected =
          dt / t * previous_sigma.doubleDot(grad_u + grad_u.transpose()) / 2.;
      Matrix<Real, 3, 3> delta_inelastic_strain_expected =
          delta_lambda_expected * 3. / 2. / sigma_0 * previous_sigma;
      Matrix<Real, 3, 3> inelastic_strain_expected =
          delta_inelastic_strain_expected + previous_inelastic_strain;
      ASSERT_NEAR((inelastic_strain - inelastic_strain_expected)
                      .lpNorm<Eigen::Infinity>(),
                  0., 1e-13);
      Matrix<Real, 3, 3> delta_sigma_expected =
          2. * shear_modulus_mu *
          (0.5 * dt / t * (grad_u + grad_u.transpose()) -
           delta_inelastic_strain_expected);

      Matrix<Real, 3, 3> delta_sigma = sigma - previous_sigma;
      ASSERT_NEAR(
          (delta_sigma_expected - delta_sigma).lpNorm<Eigen::Infinity>(), 0.,
          1e-13);
    }
    previous_sigma = sigma;
    previous_sigma_rot = sigma_rot;
    previous_grad_u_rot = grad_u_rot;
    previous_inelastic_strain = inelastic_strain;
    previous_inelastic_strain_rot = inelastic_strain_rot;
  }
}

namespace {

template <typename T>
class TestPlasticMaterialFixture : public ::TestMaterialFixture<T> {};

TYPED_TEST_SUITE(TestPlasticMaterialFixture, mat_types, );

TYPED_TEST(TestPlasticMaterialFixture, ComputeStress) {
  this->material->testComputeStress();
}
TYPED_TEST(TestPlasticMaterialFixture, DISABLED_EnergyDensity) {
  this->material->testEnergyDensity();
}
TYPED_TEST(TestPlasticMaterialFixture, DISABLED_ComputeTangentModuli) {
  this->material->testComputeTangentModuli();
}
TYPED_TEST(TestPlasticMaterialFixture, DISABLED_ComputeCelerity) {
  this->material->testCelerity();
}
} // namespace

/*****************************************************************/
