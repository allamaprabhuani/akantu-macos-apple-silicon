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
#include "dof_manager_default.hh"
#if defined(AKANTU_USE_PETSC)
#include "dof_manager_petsc.hh"
#endif
#include "patch_test_linear_fixture.hh"
#include "solid_mechanics_model.hh"
/* --------------------------------------------------------------------------
 */

#ifndef AKANTU_PATCH_TEST_LINEAR_SOLID_MECHANICS_FIXTURE_HH_
#define AKANTU_PATCH_TEST_LINEAR_SOLID_MECHANICS_FIXTURE_HH_

/* --------------------------------------------------------------------------
 */
template <typename tuple_>
class TestPatchTestSMMLinear
    : public TestPatchTestLinear<std::tuple_element_t<0, tuple_>,
                                 SolidMechanicsModel,
                                 std::tuple_element_t<2, tuple_>> {
  using parent =
      TestPatchTestLinear<std::tuple_element_t<0, tuple_>, SolidMechanicsModel,
                          std::tuple_element_t<2, tuple_>>;

public:
  static constexpr bool plane_strain = std::tuple_element_t<1, tuple_>::value;

  void applyBC() override {
    parent::applyBC();
    auto & displacement = this->model->getDisplacement();
    this->applyBConDOFs(displacement);
  }

  void checkForces() {
    auto & mat = this->model->getMaterial(0);
    auto & internal_forces = this->model->getInternalForce();
    auto & external_forces = this->model->getExternalForce();
    auto dim = this->dim;

    Matrix<Real> sigma =
        make_view(mat.getStress(this->type), dim, dim).begin()[0];

    external_forces.zero();
    if (dim > 1) {
      for (auto & eg : this->mesh->iterateElementGroups()) {
        this->model->applyBC(BC::Neumann::FromHigherDim(sigma), eg.getName());
      }
    } else {
      external_forces(0) = -sigma(0, 0);
      external_forces(1) = sigma(0, 0);
    }

    Real force_norm_inf = -std::numeric_limits<Real>::max();

    Vector<Real> total_force(dim);
    total_force.zero();

    for (auto && f : make_view(internal_forces, dim)) {
      total_force += f;
      force_norm_inf =
          std::max(force_norm_inf, f.template lpNorm<Eigen::Infinity>());
    }

    constexpr Real force_tol = 1e-8;
    EXPECT_NEAR(0,
                total_force.template lpNorm<Eigen::Infinity>() / force_norm_inf,
                force_tol);

    for (auto && [f_int, f_ext] : zip(make_view(internal_forces, dim),
                                      make_view(external_forces, dim))) {
      auto f = f_int + f_ext;
      EXPECT_NEAR(0, f.template lpNorm<Eigen::Infinity>() / force_norm_inf,
                  force_tol);
    }
  }

  void checkAll() {
    auto & displacement = this->model->getDisplacement();
    auto & mat = this->model->getMaterial(0);

    this->checkDOFs(displacement);
    this->checkGradient(mat.getGradU(this->type), displacement);
    this->checkResults(
        [&](const Matrix<Real> & pstrain) {
          Real nu = this->model->getMaterial(0).get("nu");
          Real E = this->model->getMaterial(0).get("E");

          Matrix<Real, parent::dim, parent::dim> strain =
              (pstrain + pstrain.transpose()) / 2.;
          auto trace = strain.trace();

          auto lambda = nu * E / ((1 + nu) * (1 - 2 * nu));
          auto mu = E / (2 * (1 + nu));

          if (not this->plane_strain) {
            lambda = nu * E / (1 - nu * nu);
          }

          Matrix<Real, parent::dim, parent::dim> stress;

          if (parent::dim == 1) {
            stress = E * strain;
          } else {
            stress = Matrix<Real, parent::dim, parent::dim>::Identity() *
                         lambda * trace +
                     2 * mu * strain;
          }

          return stress;
        },
        mat.getStress(this->type), displacement);
    this->checkForces();
  }
};

template <typename tuple_>
constexpr bool TestPatchTestSMMLinear<tuple_>::plane_strain;

template <typename T> struct invalid_plan_stress : std::true_type {};
template <typename type, typename bool_c, typename DM>
struct invalid_plan_stress<std::tuple<type, bool_c, DM>>
    : aka::bool_constant<ElementClass<type::value>::getSpatialDimension() !=
                             2 and
                         not bool_c::value> {};

using true_false =
    std::tuple<aka::bool_constant<true>, aka::bool_constant<false>>;

template <class NLS, class SST> struct TestSolverOptions {
  static constexpr NonLinearSolverType nls_type = NLS::value;
  static constexpr SparseSolverType ss_type = SST::value;
};

struct _non_linear_solver_auto
    : public std::integral_constant<NonLinearSolverType,
                                    NonLinearSolverType::_auto> {};
struct _non_linear_solver_petsc
    : public std::integral_constant<NonLinearSolverType,
                                    NonLinearSolverType::_petsc_snes> {};

struct _sparse_solver_eigen
    : public std::integral_constant<SparseSolverType,
                                    SparseSolverType::_eigen> {};
struct _sparse_solver_mumps
    : public std::integral_constant<SparseSolverType,
                                    SparseSolverType::_mumps> {};
struct _sparse_solver_petsc
    : public std::integral_constant<SparseSolverType,
                                    SparseSolverType::_petsc> {};

using solver_options = std::tuple<
    std::tuple<DOFManagerDefault,
               TestSolverOptions<_non_linear_solver_auto, _sparse_solver_eigen>>
#ifdef AKANTU_USE_MUMPS
    ,
    std::tuple<DOFManagerDefault,
               TestSolverOptions<_non_linear_solver_auto, _sparse_solver_mumps>>
#endif
#ifdef AKANTU_USE_PETSC
    ,
    std::tuple<DOFManagerPETSc, TestSolverOptions<_non_linear_solver_auto,
                                                  _sparse_solver_petsc>>,
    std::tuple<DOFManagerPETSc, TestSolverOptions<_non_linear_solver_petsc,
                                                  _sparse_solver_petsc>>
#endif
    >;
template <typename T> using valid_types = aka::negation<invalid_plan_stress<T>>;

using model_types = gtest_list_t<
    tuple_filter_t<valid_types, cross_product_t<TestElementTypes, true_false,
                                                solver_options>>>;

class NameGenerator {
public:
  template <typename T> static std::string GetName(int) {
    std::string name{};
    const auto element_type = std::tuple_element_t<0, T>::value;
    const bool plane_strain = std::tuple_element_t<1, T>::value;
    using options = std::tuple_element_t<2, T>;

    using dm_type = std::tuple_element_t<0, options>;

    const auto nls_type = std::tuple_element_t<1, options>::nls_type;
    const auto ss_type = std::tuple_element_t<1, options>::ss_type;

    name = std::to_string(element_type) + "/";

    name += plane_strain ? "plane_strain/" : "plane_stress/";

    if constexpr (std::is_same_v<DOFManagerDefault, dm_type>)
      name += "default/";
#ifdef AKANTU_USE_PETSC
    if constexpr (std::is_same_v<DOFManagerPETSc, dm_type>)
      name += "petsc/";
#endif

    if constexpr (nls_type == NonLinearSolverType::_auto)
      name += "auto";
    if constexpr (nls_type == NonLinearSolverType::_petsc_snes)
      name += "petsc_snes";
    if constexpr (nls_type == NonLinearSolverType::_newton_raphson)
      name += "newton_raphson";

    if constexpr (ss_type == SparseSolverType::_eigen)
      name += "[eigen]";
    if constexpr (ss_type == SparseSolverType::_mumps)
      name += "[mumps]";
    if constexpr (ss_type == SparseSolverType::_petsc)
      name += "[petsc]";

    return name;
  }
};

TYPED_TEST_SUITE(TestPatchTestSMMLinear, model_types, NameGenerator);

#endif /* AKANTU_PATCH_TEST_LINEAR_SOLID_MECHANICS_FIXTURE_HH_ */
