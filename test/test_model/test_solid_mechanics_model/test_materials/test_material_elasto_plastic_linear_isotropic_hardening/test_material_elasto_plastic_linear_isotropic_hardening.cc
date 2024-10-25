/**
 * Copyright (©) 2014-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "non_linear_solver.hh"
#include "non_linear_solver_newton_raphson.hh"
#include "solid_mechanics_model.hh"
#if defined(AKANTU_USE_PETSC)
#include "sparse_solver_petsc.hh"
#endif

/* -------------------------------------------------------------------------- */
#include <iostream>
/* -------------------------------------------------------------------------- */

using namespace akantu;

/* -------------------------------------------------------------------------- */

int main(int argc, char * argv[]) {
  initialize("test_material_elasto_plastic_linear_isotropic_hardening.dat",
             argc, argv);

  const Int spatial_dimension = 2;
  const Real u_increment = 0.1;
  const Int steps = 20;

  Mesh mesh(spatial_dimension);
  mesh.read("test_material_elasto_plastic_linear_isotropic_hardening.msh");

  SolidMechanicsModel model(mesh);
  model.initFull(_analysis_method = _static);

  auto & solver = model.getNonLinearSolver("static");
  solver.set("max_iterations", 300);
  solver.set("threshold", 1e-5);

#if defined(AKANTU_USE_PETSC)
  if (aka::is_of_type<NonLinearSolverNewtonRaphson>(solver)) {
    auto & sparse_solver =
        aka::as_type<NonLinearSolverNewtonRaphson>(solver).getSparseSolver();

    if (aka::is_of_type<SparseSolverPETSc>(sparse_solver)) {
      sparse_solver.set("pc_type", "cholesky");
      sparse_solver.set("ksp_rtol", "1e-30");
    }
  }
#endif

  model.applyBC(BC::Dirichlet::FixedValue(0.0, _x), "left");
  model.applyBC(BC::Dirichlet::FixedValue(0.0, _y), "bottom");

  std::cout.precision(4);
  for (Int i = 0; i < steps; ++i) {

    model.applyBC(BC::Dirichlet::FixedValue(i * u_increment, _x), "right");

    try {
      model.solveStep();
    } catch (debug::NLSNotConvergedException & e) {
      std::cout << e.niter << " " << e.error << std::endl;
      throw;
    }
    Real strainxx = i * u_increment / 10.;

    const auto & edge_nodes =
        mesh.getElementGroup("right").getNodeGroup().getNodes();
    auto & residual = model.getInternalForce();
    Real reaction = 0;

    for (Int n = 0; n < edge_nodes.size(); n++) {
      reaction -= residual(edge_nodes(n), 0);
    }

    std::cout << strainxx << "," << reaction << std::endl;
  }

  finalize();
  return 0;
}
