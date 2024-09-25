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
#include "element_group.hh"
#include "mesh_utils.hh"
#include "model.hh"
#include "test_gtest_utils.hh"
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>
#include <vector>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_PATCH_TEST_LINEAR_FIXTURE_HH_
#define AKANTU_PATCH_TEST_LINEAR_FIXTURE_HH_

// #define DEBUG_TEST

using namespace akantu;

template <typename element_type_, typename Model_, typename DOFManager_>
class TestPatchTestLinear : public ::testing::Test {
public:
  static constexpr ElementType type = element_type_::value;
  static constexpr Int dim = ElementClass<type>::getSpatialDimension();

  void SetUp() override {
    mesh = std::make_unique<Mesh>(dim);
    mesh->read(std::to_string(type) + ".msh");
    MeshUtils::buildFacets(*mesh);
    mesh->createBoundaryGroupFromGeometry();

    model =
        std::make_unique<Model_>(*mesh, _all_dimensions, std::to_string(type));

    if constexpr (std::is_same_v<DOFManager_, DOFManagerPETSc>) {
      model->initDOFManager("petsc");
    }
  }

  void TearDown() override {
    model.reset(nullptr);
    mesh.reset(nullptr);
  }

  virtual void initModel(const AnalysisMethod & method,
                         const std::string & material_file) {
    debug::setDebugLevel(dblError);
    getStaticParser().parse(material_file);
    this->model->initFull(_analysis_method = method);
    this->applyBC();

    if constexpr (std::is_same_v<DOFManager_, DOFManagerPETSc>) {
      auto & solver = model.getNonLinearSolver();
      auto & sparse_solver = solver.getSparseSolver();
      sparse_solver.set("max_iterations", 2);

      model->getSolver() initDOFManager("petsc");
    }

    if (method != _static) {
      this->model->setTimeStep(0.8 * this->model->getStableTimeStep());
    }
  }

  virtual void applyBC() {
    auto & boundary = this->model->getBlockedDOFs();

    for (auto & eg : mesh->iterateElementGroups()) {
      for (const auto & node : eg.getNodeGroup()) {
        for (Int s = 0; s < boundary.getNbComponent(); ++s) {
          boundary(node, s) = true;
        }
      }
    }
  }

  virtual void applyBConDOFs(Array<Real> & dofs) {
    const auto & coordinates = this->mesh->getNodes();
    for (auto & eg : this->mesh->iterateElementGroups()) {
      for (const auto & node : eg.getNodeGroup()) {
        this->setLinearDOF(dofs.begin(dofs.getNbComponent())[node],
                           coordinates.begin(this->dim)[node]);
      }
    }
  }

  template <typename V> Matrix<Real> prescribed_gradient(const V & dof) {
    Matrix<Real> gradient(dof.getNbComponent(), dim);

    for (Int i = 0; i < gradient.rows(); ++i) {
      for (Int j = 0; j < gradient.cols(); ++j) {
        gradient(i, j) = alpha(i, j + 1);
      }
    }
    return gradient;
  }

  template <typename Gradient, typename DOFs>
  void checkGradient(const Gradient & gradient, const DOFs & dofs) {
    auto pgrad = prescribed_gradient(dofs);

    for (auto & grad :
         make_view(gradient, gradient.getNbComponent() / dim, dim)) {
      auto diff = grad - pgrad;
      auto gradient_error = diff.template lpNorm<Eigen::Infinity>() /
                            grad.template lpNorm<Eigen::Infinity>();

      EXPECT_NEAR(0, gradient_error, gradient_tolerance);
    }
  }

  template <typename presult_func_t, typename Result, typename DOFs>
  void checkResults(presult_func_t && presult_func, const Result & results,
                    const DOFs & dofs) {
    Matrix<Real> presult = presult_func(prescribed_gradient(dofs));
    for (auto & result : make_view(results, presult.rows(), presult.cols())) {
      auto diff = result - presult;
      auto result_error = diff.template lpNorm<Eigen::Infinity>() /
                          presult.template lpNorm<Eigen::Infinity>();

      EXPECT_NEAR(0, result_error, result_tolerance);
    }
  }

  template <typename V1, typename V2>
  void setLinearDOF(V1 && dof, V2 && coord) {
    for (Int i = 0; i < dof.size(); ++i) {
      dof(i) = this->alpha(i, 0);
      for (Int j = 0; j < coord.size(); ++j) {
        dof(i) += this->alpha(i, j + 1) * coord(j);
      }
    }
  }

  template <typename V> void checkDOFs(V && dofs) {
    const auto & coordinates = mesh->getNodes();
    Vector<Real> ref_dof(dofs.getNbComponent());

    for (auto && [X, u] : zip(make_view(coordinates, dim),
                              make_view(dofs, dofs.getNbComponent()))) {
      setLinearDOF(ref_dof, X);
      auto diff = u - ref_dof;
      auto dofs_error = diff.template lpNorm<Eigen::Infinity>();

      EXPECT_NEAR(0, dofs_error, dofs_tolerance);
    }
  }

protected:
  std::unique_ptr<Mesh> mesh;
  std::unique_ptr<Model_> model;
  Matrix<Real> alpha{{0.01, 0.02, 0.03, 0.04},
                     {0.05, 0.06, 0.07, 0.08},
                     {0.09, 0.10, 0.11, 0.12}};

  Real gradient_tolerance{1e-13};
  Real result_tolerance{1e-13};
  Real dofs_tolerance{1e-15};
};

// template <typename type_, typename M>
// constexpr ElementType TestPatchTestLinear<type_, M>::type;

// template <typename tuple_, typename M>
// constexpr Int TestPatchTestLinear<tuple_, M>::dim;

#endif /* AKANTU_PATCH_TEST_LINEAR_FIXTURE_HH_ */
