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
#include "patch_test_linear_solid_mechanics_fixture.hh"
/* -------------------------------------------------------------------------- */
#include "non_linear_solver.hh"
/* -------------------------------------------------------------------------- */
constexpr Real ekin_tol = 1e-13;

using namespace akantu;

TYPED_TEST(TestPatchTestSMMLinear, Implicit) {
  std::string filename = "material_check_stress_plane_stress.dat";
  if (this->plane_strain)
    filename = "material_check_stress_plane_strain.dat";

  this->initModel(_implicit_dynamic, filename);

  const auto & coordinates = this->mesh->getNodes();
  auto & displacement = this->model->getDisplacement();
  // set the position of all nodes to the static solution
  for (auto && [X, u] : zip(make_view(coordinates, this->dim),
                            make_view(displacement, this->dim))) {
    this->setLinearDOF(u, X);
  }

  auto & solver = this->model->getNonLinearSolver();
  if (TestFixture::nls_type == NonLinearSolverType::_petsc_snes) {
    solver.set("snes_max_it", 2);
    solver.set("snes_rtol", 1e-10);
    solver.set("snes_atol", 1e-5);
  } else {
    solver.set("max_iterations", 2);
    solver.set("threshold", 1e-10);
    solver.set("convergence_type", SolveConvergenceCriteria::_residual);
  }

  for (Int s = 0; s < 100; ++s) {
    this->model->solveStep();
  }

  auto ekin = this->model->getEnergy("kinetic");
  EXPECT_NEAR(0, ekin, ekin_tol);

  this->checkAll();

#define debug 0
#if debug
  this->model->setBaseName(std::to_string(this->type) + "_implicit");
  this->model->addDumpField("stress");
  this->model->addDumpField("grad_u");
  this->model->addDumpFieldVector("internal_force");
  this->model->addDumpFieldVector("external_force");
  this->model->addDumpField("blocked_dofs");
  this->model->addDumpFieldVector("displacement");
  this->model->dump();
#endif
}

/* -------------------------------------------------------------------------- */
TYPED_TEST(TestPatchTestSMMLinear, Static) {
  std::string filename = "material_check_stress_plane_stress.dat";
  if (this->plane_strain)
    filename = "material_check_stress_plane_strain.dat";

  this->initModel(_static, filename);

  auto & solver = this->model->getNonLinearSolver();
  if (TestFixture::nls_type == NonLinearSolverType::_petsc_snes) {
    solver.set("snes_max_it", 2);
    solver.set("snes_rtol", 1e-10);
  } else {
    solver.set("max_iterations", 2);
    solver.set("threshold", 1e-10);
    solver.set("convergence_type", SolveConvergenceCriteria::_residual);
  }
  this->model->solveStep();

  this->checkAll();

#define debug 0
#if debug
  this->model->setBaseName(std::to_string(this->type) + "_static");
  this->model->addDumpField("stress");
  this->model->addDumpField("grad_u");
  this->model->addDumpFieldVector("internal_force");
  this->model->addDumpFieldVector("external_force");
  this->model->addDumpField("blocked_dofs");
  this->model->addDumpFieldVector("displacement");
  this->model->dump();
#endif
}

/* -------------------------------------------------------------------------- */
TYPED_TEST(TestPatchTestSMMLinear, ImplicitFiniteDeformation) {
  std::string filename =
      "material_check_stress_plane_stress_finite_deformation.dat";
  if (this->plane_strain)
    filename = "material_check_stress_plane_strain_finite_deformation.dat";
  else {
    SUCCEED();
    return;
  }

  this->initModel(_implicit_dynamic, filename);

  const auto & coordinates = this->mesh->getNodes();
  auto & displacement = this->model->getDisplacement();
  // set the position of all nodes to the static solution
  for (auto && tuple : zip(make_view(coordinates, this->dim),
                           make_view(displacement, this->dim))) {
    this->setLinearDOF(std::get<1>(tuple), std::get<0>(tuple));
  }
  for (Int s = 0; s < 100; ++s) {
    this->model->solveStep();
  }

  auto ekin = this->model->getEnergy("kinetic");
  EXPECT_NEAR(0, ekin, ekin_tol);

  this->checkAll();

#define debug 0
#if debug
  this->model->setBaseName(std::to_string(this->type) + "_implicit_finit_def");
  this->model->addDumpField("stress");
  this->model->addDumpField("grad_u");
  this->model->addDumpFieldVector("internal_force");
  this->model->addDumpFieldVector("external_force");
  this->model->addDumpField("blocked_dofs");
  this->model->addDumpFieldVector("displacement");
  this->model->dump();
#endif
}

/* -------------------------------------------------------------------------- */
TYPED_TEST(TestPatchTestSMMLinear, StaticFiniteDeformation) {
  std::string filename =
      "material_check_stress_plane_stress_finite_deformation.dat";
  if (this->plane_strain) {
    filename = "material_check_stress_plane_strain_finite_deformation.dat";
  } else {
    SUCCEED();
    return;
  }
  this->initModel(_static, filename);

  auto & solver = this->model->getNonLinearSolver();
  if (TestFixture::nls_type == NonLinearSolverType::_petsc_snes) {
    solver.set("snes_max_it", 2);
    solver.set("snes_rtol", 1e-10);
    solver.set("snes_atol", 1e-10);
  } else {
    solver.set("max_iterations", 2);
    solver.set("threshold", 1e-10);
    solver.set("convergence_type", SolveConvergenceCriteria::_residual);
  }

  this->model->solveStep();

  this->checkAll();
}
