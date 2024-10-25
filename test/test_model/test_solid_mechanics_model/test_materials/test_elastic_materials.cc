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
#include "test_gtest_utils.hh"
#include "test_material_fixtures.hh"
/* -------------------------------------------------------------------------- */
#include <aka_voigthelper.hh>
#include <material_elastic.hh>
#include <material_elastic_orthotropic.hh>
#include <solid_mechanics_model.hh>
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>
#include <type_traits>
/* -------------------------------------------------------------------------- */

using namespace akantu;

using mat_types =
    ::testing::Types<Traits<MaterialElastic, 1>, Traits<MaterialElastic, 2>,
                     Traits<MaterialElastic, 3>,

                     Traits<MaterialElasticOrthotropic, 2>,
                     Traits<MaterialElasticOrthotropic, 3>,

                     Traits<MaterialElasticLinearAnisotropic, 2>,
                     Traits<MaterialElasticLinearAnisotropic, 3>>;

/* -------------------------------------------------------------------------- */
template <> void FriendMaterial<MaterialElastic<1>>::setParams() {
  Real E = 3.;
  Real rho = 2;
  set("E", E);
  set("rho", rho);
}

/* -------------------------------------------------------------------------- */
template <> void FriendMaterial<MaterialElastic<1>>::testComputeStress() {
  Matrix<Real, 1, 1> eps{{2}};
  Matrix<Real, 1, 1> sigma;
  Real sigma_th = 2.;
  this->computeStressOnQuad(make_named_tuple(
      "grad_u"_n = eps, "sigma"_n = sigma, "sigma_th"_n = sigma_th));

  auto solution = E * eps(0, 0) + sigma_th;
  EXPECT_NEAR(sigma(0, 0), solution, 1e-14);
}

/* -------------------------------------------------------------------------- */
template <> void FriendMaterial<MaterialElastic<1>>::testEnergyDensity() {
  Real eps = 2, sigma = 2;
  Real epot = 0;
  this->computePotentialEnergyOnQuad(
      make_named_tuple("grad_u"_n = Matrix<Real>{{eps}},
                       "sigma"_n = Matrix<Real>{{sigma}}),
      epot);
  Real solution = 2;
  EXPECT_NEAR(epot, solution, 1e-14);
}

/* -------------------------------------------------------------------------- */
template <>
void FriendMaterial<MaterialElastic<1>>::testComputeTangentModuli() {
  Matrix<Real> tangent(1, 1);
  this->computeTangentModuliOnQuad(
      make_named_tuple("tangent_moduli"_n = tangent));
  EXPECT_NEAR(tangent(0, 0), E, 1e-14);
}

/* -------------------------------------------------------------------------- */
template <> void FriendMaterial<MaterialElastic<1>>::testCelerity() {
  auto wave_speed = this->getCelerity(Element());
  auto solution = std::sqrt(E / rho);
  EXPECT_NEAR(wave_speed, solution, 1e-14);
}

/* -------------------------------------------------------------------------- */
template <> void FriendMaterial<MaterialElastic<2>>::setParams() {
  Real E = 1.;
  Real nu = .3;
  Real rho = 2;
  set("E", E);
  set("nu", nu);
  set("rho", rho);
}

/* -------------------------------------------------------------------------- */
template <> void FriendMaterial<MaterialElastic<2>>::testComputeStress() {
  Real bulk_modulus_K = E / (3 * (1 - 2 * nu));
  Real shear_modulus_mu = E / (2 * (1 + nu));
  Real sigma_th = 0.;
  Matrix<Real> rotation_matrix = getRandomRotation();

  Matrix<Real, 2, 2> grad_u = this->getComposedStrain(1.).block<2, 2>(0, 0);

  auto grad_u_rot = this->applyRotation(grad_u, rotation_matrix);

  Matrix<Real, 2, 2> sigma_rot;
  this->computeStressOnQuad(make_named_tuple(
      "grad_u"_n = grad_u_rot, "sigma"_n = sigma_rot, "sigma_th"_n = sigma_th));

  auto sigma = this->reverseRotation(sigma_rot, rotation_matrix);

  auto identity = Matrix<Real, 2, 2>::Identity();

  Matrix<Real, 2, 2> strain = 0.5 * (grad_u + grad_u.transpose());
  Matrix<Real, 2, 2> deviatoric_strain =
      strain - 1. / 3. * strain.trace() * identity;

  Matrix<Real, 2, 2> sigma_expected =
      2 * shear_modulus_mu * deviatoric_strain +
      (sigma_th + 2. * bulk_modulus_K) * identity;

  auto diff = sigma - sigma_expected;
  Real stress_error =
      diff.lpNorm<Eigen::Infinity>() / sigma_expected.lpNorm<Eigen::Infinity>();
  EXPECT_NEAR(stress_error, 0., 1e-13);
}

/* -------------------------------------------------------------------------- */
template <> void FriendMaterial<MaterialElastic<2>>::testEnergyDensity() {
  Matrix<Real, 2, 2> sigma{{1, 2}, {2, 4}};
  Matrix<Real, 2, 2> eps{{1, 0}, {0, 1}};
  Real epot = 0;
  Real solution = 2.5;
  this->computePotentialEnergyOnQuad(
      make_named_tuple("grad_u"_n = eps, "sigma"_n = sigma), epot);
  EXPECT_NEAR(epot, solution, 1e-14);
}

/* -------------------------------------------------------------------------- */
template <>
void FriendMaterial<MaterialElastic<2>>::testComputeTangentModuli() {
  Matrix<Real, 3, 3> tangent;

  /* Plane Strain */
  // clang-format off
  Matrix<Real, 3,3> solution{
    {1 - nu, nu, 0},
    {nu, 1 - nu, 0},
    {0, 0, (1 - 2 * nu) / 2},
  };
  // clang-format on
  solution *= E / ((1 + nu) * (1 - 2 * nu));

  this->computeTangentModuliOnQuad(
      make_named_tuple("tangent_moduli"_n = tangent));
  Real tangent_error = (tangent - solution).lpNorm<Eigen::Infinity>();
  EXPECT_NEAR(tangent_error, 0, 1e-14);

  /* Plane Stress */
  this->plane_stress = true;
  this->updateInternalParameters();
  // clang-format off
  solution = Matrix<Real, 3,3>{
    {1, nu, 0},
    {nu, 1, 0},
    {0, 0, (1 - nu) / 2},
  };
  // clang-format on
  solution *= E / (1 - nu * nu);

  this->computeTangentModuliOnQuad(
      make_named_tuple("tangent_moduli"_n = tangent));
  tangent_error = (tangent - solution).lpNorm<Eigen::Infinity>();
  EXPECT_NEAR(tangent_error, 0, 1e-14);
}

/* -------------------------------------------------------------------------- */
template <> void FriendMaterial<MaterialElastic<2>>::testCelerity() {
  auto push_wave_speed = this->getPushWaveSpeed(Element());
  auto celerity = this->getCelerity(Element());

  Real K = E / (3 * (1 - 2 * nu));
  Real mu = E / (2 * (1 + nu));
  Real sol = std::sqrt((K + 4. / 3 * mu) / rho);

  EXPECT_NEAR(push_wave_speed, sol, 1e-14);
  EXPECT_NEAR(celerity, sol, 1e-14);

  auto shear_wave_speed = this->getShearWaveSpeed(Element());

  sol = std::sqrt(mu / rho);

  EXPECT_NEAR(shear_wave_speed, sol, 1e-14);
}

/* -------------------------------------------------------------------------- */
template <> void FriendMaterial<MaterialElastic<3>>::setParams() {
  Real E = 1.;
  Real nu = .3;
  Real rho = 2;
  set("E", E);
  set("nu", nu);
  set("rho", rho);
}

/* -------------------------------------------------------------------------- */
template <> void FriendMaterial<MaterialElastic<3>>::testComputeStress() {
  Real bulk_modulus_K = E / 3. / (1 - 2. * nu);
  Real shear_modulus_mu = 0.5 * E / (1 + nu);

  Matrix<Real, 3, 3> rotation_matrix = getRandomRotation();
  Matrix<Real, 3, 3> grad_u = this->getComposedStrain(1.);
  Matrix<Real, 3, 3> grad_u_rot = this->applyRotation(grad_u, rotation_matrix);

  Real sigma_th{0.};

  Matrix<Real, 3, 3> sigma_rot;
  this->computeStressOnQuad(make_named_tuple(
      "grad_u"_n = grad_u_rot, "sigma"_n = sigma_rot, "sigma_th"_n = sigma_th));

  Matrix<Real, 3, 3> sigma = this->reverseRotation(sigma_rot, rotation_matrix);

  Matrix<Real, 3, 3> identity = Matrix<Real, 3, 3>::Identity();

  Matrix<Real, 3, 3> strain = 0.5 * (grad_u + grad_u.transpose());
  Matrix<Real, 3, 3> deviatoric_strain =
      strain - 1. / 3. * strain.trace() * identity;

  Matrix<Real, 3, 3> sigma_expected =
      2 * shear_modulus_mu * deviatoric_strain +
      (sigma_th + 3. * bulk_modulus_K) * identity;

  auto diff = sigma - sigma_expected;
  Real stress_error = diff.lpNorm<Eigen::Infinity>();
  EXPECT_NEAR(stress_error, 0., 1e-14);
}

/* -------------------------------------------------------------------------- */
template <> void FriendMaterial<MaterialElastic<3>>::testEnergyDensity() {
  Matrix<Real, 3, 3> sigma{{1, 2, 3}, {2, 4, 5}, {3, 5, 6}};
  Matrix<Real, 3, 3> eps{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  Real epot = 0;
  Real solution = 5.5;
  this->computePotentialEnergyOnQuad(
      make_named_tuple("grad_u"_n = eps, "sigma"_n = sigma), epot);
  EXPECT_NEAR(epot, solution, 1e-14);
}

/* -------------------------------------------------------------------------- */
template <>
void FriendMaterial<MaterialElastic<3>>::testComputeTangentModuli() {
  Matrix<Real> tangent(6, 6);

  // clang-format off
  Matrix<Real,6,6> solution{
      {1 - nu, nu, nu, 0, 0, 0},
      {nu, 1 - nu, nu, 0, 0, 0},
      {nu, nu, 1 - nu, 0, 0, 0},
      {0, 0, 0, (1 - 2 * nu) / 2, 0, 0},
      {0, 0, 0, 0, (1 - 2 * nu) / 2, 0},
      {0, 0, 0, 0, 0, (1 - 2 * nu) / 2},
  };
  // clang-format on
  solution *= E / ((1 + nu) * (1 - 2 * nu));

  this->computeTangentModuliOnQuad(
      make_named_tuple("tangent_moduli"_n = tangent));
  Real tangent_error = (tangent - solution).norm();
  EXPECT_NEAR(tangent_error, 0, 1e-14);
}

/* -------------------------------------------------------------------------- */
template <> void FriendMaterial<MaterialElastic<3>>::testCelerity() {
  auto push_wave_speed = this->getPushWaveSpeed(Element());
  auto celerity = this->getCelerity(Element());

  Real K = E / (3 * (1 - 2 * nu));
  Real mu = E / (2 * (1 + nu));
  Real sol = std::sqrt((K + 4. / 3 * mu) / rho);

  EXPECT_NEAR(push_wave_speed, sol, 1e-14);
  EXPECT_NEAR(celerity, sol, 1e-14);

  auto shear_wave_speed = this->getShearWaveSpeed(Element());

  sol = std::sqrt(mu / rho);

  EXPECT_NEAR(shear_wave_speed, sol, 1e-14);
}

/* -------------------------------------------------------------------------- */
template <> void FriendMaterial<MaterialElasticOrthotropic<2>>::setParams() {
  // Note: for this test material and canonical basis coincide
  Vector<Real, 2> n1{1, 0};
  Vector<Real, 2> n2{0, 1};
  Real E1 = 1.;
  Real E2 = 2.;
  Real nu12 = 0.1;
  Real G12 = 2.;
  Real rho = 2.5;

  *this->dir_vecs[0] = n1;
  *this->dir_vecs[1] = n2;
  this->E1 = E1;
  this->E2 = E2;
  this->nu12 = nu12;
  this->G12 = G12;
  this->rho = rho;
}

/* -------------------------------------------------------------------------- */
template <>
void FriendMaterial<MaterialElasticOrthotropic<2>>::testComputeStress() {
  const Int Dim = 2;
  // material frame of reference is rotate by rotation_matrix starting from
  // canonical basis
  Matrix<Real> rotation_matrix = getRandomRotation();

  // canonical basis as expressed in the material frame of reference, as
  // required by MaterialElasticOrthotropic class (it is simply given by the
  // columns of the rotation_matrix; the lines give the material basis expressed
  // in the canonical frame of reference)
  *this->dir_vecs[0] = rotation_matrix(0);
  *this->dir_vecs[1] = rotation_matrix(1);

  // set internal Cijkl matrix expressed in the canonical frame of reference
  this->updateInternalParameters();

  // gradient in material frame of reference
  Matrix<Real, Dim, Dim> grad_u = this->getComposedStrain(2.).block<2, 2>(0, 0);

  // gradient in canonical basis (we need to rotate *back* to the canonical
  // basis)
  auto grad_u_rot = this->reverseRotation(grad_u, rotation_matrix);

  // stress in the canonical basis
  Matrix<Real, Dim, Dim> sigma_rot;
  this->computeStressOnQuad(
      make_named_tuple("grad_u"_n = grad_u_rot, "sigma"_n = sigma_rot));

  // stress in the material reference (we need to apply the rotation)
  auto sigma = this->applyRotation(sigma_rot, rotation_matrix);

  // construction of Cijkl engineering tensor in the *material* frame of
  // reference
  // ref: http://solidmechanics.org/Text/Chapter3_2/Chapter3_2.php#Sect3_2_13
  Real nu21 = nu12 * E2 / E1;
  Real gamma = 1 / (1 - nu12 * nu21);

  Matrix<Real, 2 * Dim, 2 * Dim> C_expected;
  C_expected.zero();

  C_expected(0, 0) = gamma * E1;
  C_expected(1, 1) = gamma * E2;
  C_expected(2, 2) = G12;

  C_expected(1, 0) = C_expected(0, 1) = gamma * E1 * nu21;

  // epsilon is computed directly in the *material* frame of reference
  Matrix<Real, Dim, Dim> epsilon = (grad_u + grad_u.transpose()) / 2.;

  // sigma_expected is computed directly in the *material* frame of reference
  Matrix<Real, Dim, Dim> sigma_expected;
  sigma_expected.zero();

  for (Int i = 0; i < Dim; ++i) {
    for (Int j = 0; j < Dim; ++j) {
      sigma_expected(i, i) += C_expected(i, j) * epsilon(j, j);
    }
  }

  sigma_expected(0, 1) = sigma_expected(1, 0) =
      C_expected(2, 2) * 2 * epsilon(0, 1);

  // sigmas are checked in the *material* frame of reference
  auto diff = sigma - sigma_expected;
  Real stress_error = diff.lpNorm<Eigen::Infinity>();
  EXPECT_NEAR(stress_error, 0., 1e-13);
}

/* -------------------------------------------------------------------------- */
template <>
void FriendMaterial<MaterialElasticOrthotropic<2>>::testEnergyDensity() {
  Matrix<Real, 2, 2> sigma{{1, 2}, {2, 4}};
  Matrix<Real, 2, 2> eps{{1, 0}, {0, 1}};
  Real epot = 0;
  Real solution = 2.5;
  this->computePotentialEnergyOnQuad(
      make_named_tuple("grad_u"_n = eps, "sigma"_n = sigma), epot);
  EXPECT_NEAR(epot, solution, 1e-14);
}

/* -------------------------------------------------------------------------- */
template <>
void FriendMaterial<MaterialElasticOrthotropic<2>>::testComputeTangentModuli() {
  // construction of Cijkl engineering tensor in the *material* frame of
  // reference
  // ref: http://solidmechanics.org/Text/Chapter3_2/Chapter3_2.php#Sect3_2_13
  Real nu21 = nu12 * E2 / E1;
  Real gamma = 1 / (1 - nu12 * nu21);

  Matrix<Real> C_expected(3, 3);
  C_expected.zero();

  C_expected(0, 0) = gamma * E1;
  C_expected(1, 1) = gamma * E2;
  C_expected(2, 2) = G12;

  C_expected(1, 0) = C_expected(0, 1) = gamma * E1 * nu21;

  Matrix<Real> tangent(3, 3);
  this->computeTangentModuliOnQuad(
      make_named_tuple("tangent_moduli"_n = tangent));

  Real tangent_error = (tangent - C_expected).norm();
  EXPECT_NEAR(tangent_error, 0, 1e-14);
}

/* -------------------------------------------------------------------------- */

template <> void FriendMaterial<MaterialElasticOrthotropic<2>>::testCelerity() {
  // construction of Cijkl engineering tensor in the *material* frame of
  // reference
  // ref: http://solidmechanics.org/Text/Chapter3_2/Chapter3_2.php#Sect3_2_13
  Real nu21 = nu12 * E2 / E1;
  Real gamma = 1 / (1 - nu12 * nu21);

  Matrix<Real, 3, 3> C_expected;
  C_expected.zero();

  C_expected(0, 0) = gamma * E1;
  C_expected(1, 1) = gamma * E2;
  C_expected(2, 2) = G12;

  C_expected(1, 0) = C_expected(0, 1) = gamma * E1 * nu21;

  Vector<Real, 3> eig_expected;
  C_expected.eig(eig_expected);

  auto celerity_expected = std::sqrt(eig_expected(0) / rho);

  auto celerity = this->getCelerity(Element());

  EXPECT_NEAR(celerity_expected, celerity, 1e-14);
}

/* -------------------------------------------------------------------------- */
template <> void FriendMaterial<MaterialElasticOrthotropic<3>>::setParams() {
  Vector<Real, 3> n1{1, 0, 0};
  Vector<Real, 3> n2{0, 1, 0};
  Vector<Real, 3> n3{0, 0, 1};
  Real E1 = 1.;
  Real E2 = 2.;
  Real E3 = 3.;
  Real nu12 = 0.1;
  Real nu13 = 0.2;
  Real nu23 = 0.3;
  Real G12 = 2.;
  Real G13 = 3.;
  Real G23 = 1.;
  Real rho = 2.3;

  *this->dir_vecs[0] = n1;
  *this->dir_vecs[1] = n2;
  *this->dir_vecs[2] = n3;
  this->E1 = E1;
  this->E2 = E2;
  this->E3 = E3;
  this->nu12 = nu12;
  this->nu13 = nu13;
  this->nu23 = nu23;
  this->G12 = G12;
  this->G13 = G13;
  this->G23 = G23;
  this->rho = rho;
}

/* -------------------------------------------------------------------------- */
template <>
void FriendMaterial<MaterialElasticOrthotropic<3>>::testComputeStress() {
  const Int Dim = 3;

  // material frame of reference is rotate by rotation_matrix starting from
  // canonical basis
  Matrix<Real, Dim, Dim> rotation_matrix = getRandomRotation();

  // canonical basis as expressed in the material frame of reference, as
  // required by MaterialElasticOrthotropic class (it is simply given by the
  // columns of the rotation_matrix; the lines give the material basis expressed
  // in the canonical frame of reference)
  *this->dir_vecs[0] = rotation_matrix(0);
  *this->dir_vecs[1] = rotation_matrix(1);
  *this->dir_vecs[2] = rotation_matrix(2);

  // set internal Cijkl matrix expressed in the canonical frame of reference
  this->updateInternalParameters();

  // gradient in material frame of reference
  Matrix<Real, Dim, Dim> grad_u = this->getComposedStrain(2.);

  // gradient in canonical basis (we need to rotate *back* to the canonical
  // basis)
  auto grad_u_rot = this->reverseRotation(grad_u, rotation_matrix);

  // stress in the canonical basis
  Matrix<Real> sigma_rot(3, 3);
  this->computeStressOnQuad(
      make_named_tuple("grad_u"_n = grad_u_rot, "sigma"_n = sigma_rot));

  // stress in the material reference (we need to apply the rotation)
  auto sigma = this->applyRotation(sigma_rot, rotation_matrix);

  // construction of Cijkl engineering tensor in the *material* frame of
  // reference
  // ref: http://solidmechanics.org/Text/Chapter3_2/Chapter3_2.php#Sect3_2_13
  Real nu21 = nu12 * E2 / E1;
  Real nu31 = nu13 * E3 / E1;
  Real nu32 = nu23 * E3 / E2;
  Real gamma = 1 / (1 - nu12 * nu21 - nu23 * nu32 - nu31 * nu13 -
                    2 * nu21 * nu32 * nu13);

  Matrix<Real, 2 * Dim, 2 * Dim> C_expected;
  C_expected.zero();

  C_expected(0, 0) = gamma * E1 * (1 - nu23 * nu32);
  C_expected(1, 1) = gamma * E2 * (1 - nu13 * nu31);
  C_expected(2, 2) = gamma * E3 * (1 - nu12 * nu21);

  C_expected(1, 0) = C_expected(0, 1) = gamma * E1 * (nu21 + nu31 * nu23);
  C_expected(2, 0) = C_expected(0, 2) = gamma * E1 * (nu31 + nu21 * nu32);
  C_expected(2, 1) = C_expected(1, 2) = gamma * E2 * (nu32 + nu12 * nu31);

  C_expected(3, 3) = G23;
  C_expected(4, 4) = G13;
  C_expected(5, 5) = G12;

  // epsilon is computed directly in the *material* frame of reference
  Matrix<Real, Dim, Dim> epsilon = 0.5 * (grad_u + grad_u.transpose());

  // sigma_expected is computed directly in the *material* frame of reference
  Matrix<Real, Dim, Dim> sigma_expected;
  sigma_expected.zero();

  for (Int i = 0; i < Dim; ++i) {
    for (Int j = 0; j < Dim; ++j) {
      sigma_expected(i, i) += C_expected(i, j) * epsilon(j, j);
    }
  }

  sigma_expected(0, 1) = C_expected(5, 5) * 2 * epsilon(0, 1);
  sigma_expected(0, 2) = C_expected(4, 4) * 2 * epsilon(0, 2);
  sigma_expected(1, 2) = C_expected(3, 3) * 2 * epsilon(1, 2);
  sigma_expected(1, 0) = sigma_expected(0, 1);
  sigma_expected(2, 0) = sigma_expected(0, 2);
  sigma_expected(2, 1) = sigma_expected(1, 2);

  // sigmas are checked in the *material* frame of reference
  auto diff = sigma - sigma_expected;
  Real stress_error = diff.lpNorm<Eigen::Infinity>();
  EXPECT_NEAR(stress_error, 0., 1e-13);
}

/* -------------------------------------------------------------------------- */
template <>
void FriendMaterial<MaterialElasticOrthotropic<3>>::testEnergyDensity() {
  Matrix<Real, 3, 3> sigma{{1, 2, 3}, {2, 4, 5}, {3, 5, 6}};
  Matrix<Real, 3, 3> eps{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  Real epot = 0;
  Real solution = 5.5;
  this->computePotentialEnergyOnQuad(
      make_named_tuple("grad_u"_n = eps, "sigma"_n = sigma), epot);
  EXPECT_NEAR(epot, solution, 1e-14);
}

/* -------------------------------------------------------------------------- */
template <>
void FriendMaterial<MaterialElasticOrthotropic<3>>::testComputeTangentModuli() {
  // Note: for this test material and canonical basis coincide
  UInt Dim = 3;

  // construction of Cijkl engineering tensor in the *material* frame of
  // reference
  // ref: http://solidmechanics.org/Text/Chapter3_2/Chapter3_2.php#Sect3_2_13
  auto nu21 = nu12 * E2 / E1;
  auto nu31 = nu13 * E3 / E1;
  auto nu32 = nu23 * E3 / E2;
  auto gamma = 1 / (1 - nu12 * nu21 - nu23 * nu32 - nu31 * nu13 -
                    2 * nu21 * nu32 * nu13);

  Matrix<Real> C_expected(2 * Dim, 2 * Dim);
  C_expected.zero();

  C_expected(0, 0) = gamma * E1 * (1 - nu23 * nu32);
  C_expected(1, 1) = gamma * E2 * (1 - nu13 * nu31);
  C_expected(2, 2) = gamma * E3 * (1 - nu12 * nu21);

  C_expected(1, 0) = C_expected(0, 1) = gamma * E1 * (nu21 + nu31 * nu23);
  C_expected(2, 0) = C_expected(0, 2) = gamma * E1 * (nu31 + nu21 * nu32);
  C_expected(2, 1) = C_expected(1, 2) = gamma * E2 * (nu32 + nu12 * nu31);

  C_expected(3, 3) = G23;
  C_expected(4, 4) = G13;
  C_expected(5, 5) = G12;

  Matrix<Real> tangent(6, 6);
  this->computeTangentModuliOnQuad(
      make_named_tuple("tangent_moduli"_n = tangent));

  Real tangent_error = (tangent - C_expected).norm();
  EXPECT_NEAR(tangent_error, 0, 1e-14);
}

/* -------------------------------------------------------------------------- */
template <> void FriendMaterial<MaterialElasticOrthotropic<3>>::testCelerity() {
  // Note: for this test material and canonical basis coincide
  UInt Dim = 3;
  // construction of Cijkl engineering tensor in the *material* frame of
  // reference
  // ref: http://solidmechanics.org/Text/Chapter3_2/Chapter3_2.php#Sect3_2_13
  auto nu21 = nu12 * E2 / E1;
  auto nu31 = nu13 * E3 / E1;
  auto nu32 = nu23 * E3 / E2;
  auto gamma = 1 / (1 - nu12 * nu21 - nu23 * nu32 - nu31 * nu13 -
                    2 * nu21 * nu32 * nu13);

  Matrix<Real> C_expected(2 * Dim, 2 * Dim);
  C_expected.zero();

  C_expected(0, 0) = gamma * E1 * (1 - nu23 * nu32);
  C_expected(1, 1) = gamma * E2 * (1 - nu13 * nu31);
  C_expected(2, 2) = gamma * E3 * (1 - nu12 * nu21);

  C_expected(1, 0) = C_expected(0, 1) = gamma * E1 * (nu21 + nu31 * nu23);
  C_expected(2, 0) = C_expected(0, 2) = gamma * E1 * (nu31 + nu21 * nu32);
  C_expected(2, 1) = C_expected(1, 2) = gamma * E2 * (nu32 + nu12 * nu31);

  C_expected(3, 3) = G23;
  C_expected(4, 4) = G13;
  C_expected(5, 5) = G12;

  Vector<Real> eig_expected(6);
  C_expected.eig(eig_expected);

  auto celerity_expected = std::sqrt(eig_expected(0) / rho);

  auto celerity = this->getCelerity(Element());

  EXPECT_NEAR(celerity_expected, celerity, 1e-14);
}

/* -------------------------------------------------------------------------- */
template <>
void FriendMaterial<MaterialElasticLinearAnisotropic<2>>::setParams() {
  Matrix<Real, 3, 3> C{{1.0, 0.3, 0.4}, {0.3, 2.0, 0.1}, {0.4, 0.1, 1.5}};

  Cprime.block<3, 3>(0, 0) = C;

  this->rho = 2.7;

  // material frame of reference is rotate by rotation_matrix starting from
  // canonical basis
  Matrix<Real> rotation_matrix = getRandomRotation();

  // canonical basis as expressed in the material frame of reference, as
  // required by MaterialElasticLinearAnisotropic class (it is simply given by
  // the columns of the rotation_matrix; the lines give the material basis
  // expressed in the canonical frame of reference)
  *this->dir_vecs[0] = rotation_matrix(0);
  *this->dir_vecs[1] = rotation_matrix(1);
}

/* -------------------------------------------------------------------------- */
template <>
void FriendMaterial<MaterialElasticLinearAnisotropic<2>>::testComputeStress() {
  Matrix<Real, 3, 3> C{{1.0, 0.3, 0.4}, {0.3, 2.0, 0.1}, {0.4, 0.1, 1.5}};

  Matrix<Real, 2, 2> rotation_matrix;

  rotation_matrix(0) = *this->dir_vecs[0];
  rotation_matrix(1) = *this->dir_vecs[1];

  // gradient in material frame of reference
  Matrix<Real, 2, 2> grad_u = this->getComposedStrain(1.).block<2, 2>(0, 0);

  // gradient in canonical basis (we need to rotate *back* to the canonical
  // basis)
  Matrix<Real, 2, 2> grad_u_rot =
      this->reverseRotation(grad_u, rotation_matrix);

  // stress in the canonical basis
  Matrix<Real, 2, 2> sigma_rot;
  this->computeStressOnQuad(
      make_named_tuple("grad_u"_n = grad_u_rot, "sigma"_n = sigma_rot));

  // stress in the material reference (we need to apply the rotation)
  auto sigma = this->applyRotation(sigma_rot, rotation_matrix);

  // epsilon is computed directly in the *material* frame of reference
  Matrix<Real, 2, 2> epsilon = (grad_u + grad_u.transpose()) / 2.;
  Vector<Real, 3> epsilon_voigt =
      VoigtHelper<2>::matrixToVoigtWithFactors(epsilon);

  // sigma_expected is computed directly in the *material* frame of reference
  Vector<Real, 3> sigma_voigt = C * epsilon_voigt;

  Matrix<Real, 2, 2> sigma_expected =
      VoigtHelper<2>::voigtToMatrix(sigma_voigt);

  // sigmas are checked in the *material* frame of reference
  auto diff = sigma - sigma_expected;
  Real stress_error = diff.lpNorm<Eigen::Infinity>();
  EXPECT_NEAR(stress_error, 0., 1e-13);
}

/* -------------------------------------------------------------------------- */
template <>
void FriendMaterial<MaterialElasticLinearAnisotropic<2>>::testEnergyDensity() {
  Matrix<Real> sigma{{1, 2}, {2, 4}};
  Matrix<Real> eps{{1, 0}, {0, 1}};
  Real epot = 0;
  Real solution = 2.5;
  this->computePotentialEnergyOnQuad(
      make_named_tuple("grad_u"_n = eps, "sigma"_n = sigma), epot);
  EXPECT_NEAR(epot, solution, 1e-14);
}

/* -------------------------------------------------------------------------- */
template <>
void FriendMaterial<
    MaterialElasticLinearAnisotropic<2>>::testComputeTangentModuli() {
  Matrix<Real> tangent(3, 3);
  this->computeTangentModuliOnQuad(
      make_named_tuple("tangent_moduli"_n = tangent));

  Real tangent_error = (tangent - C).norm();
  EXPECT_NEAR(tangent_error, 0, 1e-14);
}

/* -------------------------------------------------------------------------- */
template <>
void FriendMaterial<MaterialElasticLinearAnisotropic<2>>::testCelerity() {
  Vector<Real> eig_expected(3);
  C.eig(eig_expected);

  auto celerity_expected = std::sqrt(eig_expected(0) / this->rho);
  auto celerity = this->getCelerity(Element());

  EXPECT_NEAR(celerity_expected, celerity, 1e-14);
}

/* -------------------------------------------------------------------------- */
template <>
void FriendMaterial<MaterialElasticLinearAnisotropic<3>>::setParams() {
  // Note: for this test material and canonical basis coincide
  Matrix<Real, 6, 6> C{
      {1.0, 0.3, 0.4, 0.3, 0.2, 0.1}, {0.3, 2.0, 0.1, 0.2, 0.3, 0.2},
      {0.4, 0.1, 1.5, 0.1, 0.4, 0.3}, {0.3, 0.2, 0.1, 2.4, 0.1, 0.4},
      {0.2, 0.3, 0.4, 0.1, 0.9, 0.1}, {0.1, 0.2, 0.3, 0.4, 0.1, 1.2}};

  for (auto i = 0; i < C.rows(); ++i)
    for (auto j = 0; j < C.cols(); ++j)
      this->Cprime(i, j) = C(i, j);

  this->rho = 2.9;

  // material frame of reference is rotate by rotation_matrix starting from
  // canonical basis
  Matrix<Real> rotation_matrix = getRandomRotation();

  // canonical basis as expressed in the material frame of reference, as
  // required by MaterialElasticLinearAnisotropic class (it is simply given by
  // the columns of the rotation_matrix; the lines give the material basis
  // expressed in the canonical frame of reference)
  *this->dir_vecs[0] = rotation_matrix(0);
  *this->dir_vecs[1] = rotation_matrix(1);
  *this->dir_vecs[2] = rotation_matrix(2);
}

/* -------------------------------------------------------------------------- */
template <>
void FriendMaterial<MaterialElasticLinearAnisotropic<3>>::testComputeStress() {
  Matrix<Real, 6, 6> C{
      {1.0, 0.3, 0.4, 0.3, 0.2, 0.1}, {0.3, 2.0, 0.1, 0.2, 0.3, 0.2},
      {0.4, 0.1, 1.5, 0.1, 0.4, 0.3}, {0.3, 0.2, 0.1, 2.4, 0.1, 0.4},
      {0.2, 0.3, 0.4, 0.1, 0.9, 0.1}, {0.1, 0.2, 0.3, 0.4, 0.1, 1.2}};

  Matrix<Real, 3, 3> rotation_matrix;

  rotation_matrix(0) = *this->dir_vecs[0];
  rotation_matrix(1) = *this->dir_vecs[1];
  rotation_matrix(2) = *this->dir_vecs[2];

  // gradient in material frame of reference
  auto grad_u = this->getComposedStrain(2.);

  // gradient in canonical basis (we need to rotate *back* to the canonical
  // basis)
  auto grad_u_rot = this->reverseRotation(grad_u, rotation_matrix);

  // stress in the canonical basis
  Matrix<Real, 3, 3> sigma_rot;
  this->computeStressOnQuad(
      make_named_tuple("grad_u"_n = grad_u_rot, "sigma"_n = sigma_rot));

  // stress in the material reference (we need to apply the rotation)
  auto sigma = this->applyRotation(sigma_rot, rotation_matrix);

  // epsilon is computed directly in the *material* frame of reference
  Matrix<Real, 3, 3> epsilon = (grad_u + grad_u.transpose()) / 2.;
  Vector<Real, 6> epsilon_voigt =
      VoigtHelper<3>::matrixToVoigtWithFactors(epsilon);

  // sigma_expected is computed directly in the *material* frame of reference
  Vector<Real, 6> sigma_voigt = C * epsilon_voigt;
  Matrix<Real, 3, 3> sigma_expected =
      VoigtHelper<3>::voigtToMatrix(sigma_voigt);

  // sigmas are checked in the *material* frame of reference
  auto diff = sigma - sigma_expected;
  Real stress_error = diff.lpNorm<Eigen::Infinity>();
  EXPECT_NEAR(stress_error, 0., 1e-13);
}

/* -------------------------------------------------------------------------- */
template <>
void FriendMaterial<MaterialElasticLinearAnisotropic<3>>::testEnergyDensity() {
  Matrix<Real, 3, 3> sigma{{1, 2, 3}, {2, 4, 5}, {3, 5, 6}};
  Matrix<Real, 3, 3> eps{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  Real epot = 0;
  Real solution = 5.5;
  this->computePotentialEnergyOnQuad(
      make_named_tuple("grad_u"_n = eps, "sigma"_n = sigma), epot);
  EXPECT_NEAR(epot, solution, 1e-14);
}

/* -------------------------------------------------------------------------- */
template <>
void FriendMaterial<
    MaterialElasticLinearAnisotropic<3>>::testComputeTangentModuli() {
  Matrix<Real, 6, 6> tangent;
  this->computeTangentModuliOnQuad(
      make_named_tuple("tangent_moduli"_n = tangent));

  Real tangent_error = (tangent - C).norm();
  EXPECT_NEAR(tangent_error, 0, 1e-14);
}

/* -------------------------------------------------------------------------- */
template <>
void FriendMaterial<MaterialElasticLinearAnisotropic<3>>::testCelerity() {
  Vector<Real, 6> eig_expected;
  C.eig(eig_expected);

  auto celerity_expected = std::sqrt(eig_expected(0) / this->rho);

  auto celerity = this->getCelerity(Element());

  EXPECT_NEAR(celerity_expected, celerity, 1e-14);
}

/* -------------------------------------------------------------------------- */

namespace {

template <typename T>
class TestElasticMaterialFixture : public ::TestMaterialFixture<T> {};

TYPED_TEST_SUITE(TestElasticMaterialFixture, mat_types, );

TYPED_TEST(TestElasticMaterialFixture, ComputeStress) {
  this->material->testComputeStress();
}

TYPED_TEST(TestElasticMaterialFixture, EnergyDensity) {
  this->material->testEnergyDensity();
}

TYPED_TEST(TestElasticMaterialFixture, ComputeTangentModuli) {
  this->material->testComputeTangentModuli();
}

TYPED_TEST(TestElasticMaterialFixture, ComputeCelerity) {
  this->material->testCelerity();
}

} // namespace
