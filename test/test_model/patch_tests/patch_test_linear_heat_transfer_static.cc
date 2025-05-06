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
#include "patch_test_linear_heat_transfer_fixture.hh"
/* -------------------------------------------------------------------------- */
#include "non_linear_solver.hh"
/* -------------------------------------------------------------------------- */

TYPED_TEST(TestPatchTestHTMLinear, Static) {
  this->initModel(_static, "heat_transfer_input.dat");

  auto & solver = this->model->getNonLinearSolver();
  if (TestFixture::nls_type == NonLinearSolverType::_petsc_snes) {
    solver.set("snes_max_it", 2);
    solver.set("snes_rtol", 2e-4);
  } else {
    solver.set("max_iterations", 2);
    solver.set("threshold", 2e-4);
    solver.set("convergence_type", SolveConvergenceCriteria::_residual);
  }

  this->model->solveStep();
  this->checkAll();
}
