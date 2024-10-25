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

#include "non_linear_solver.hh"
#include "non_linear_solver_newton_raphson.hh"
#include <solid_mechanics_model.hh>
#if defined(AKANTU_USE_PETSC)
#include "sparse_solver_petsc.hh"
#endif

using namespace akantu;

int main(int argc, char * argv[]) {
  initialize("test_multi_material_elastic.dat", argc, argv);

  const Int dim = 2;
  Mesh mesh(dim);

  mesh.read("test_multi_material_elastic.msh");

  SolidMechanicsModel model(mesh);

  auto && mat_sel = std::make_shared<MeshDataMaterialSelector<std::string>>(
      "physical_names", model);
  model.setMaterialSelector(mat_sel);

  model.initFull(_analysis_method = _static);

  model.applyBC(BC::Dirichlet::FlagOnly(_y), "ground");
  model.applyBC(BC::Dirichlet::FlagOnly(_x), "corner");

  Vector<Real> trac = Vector<Real>::Zero(dim);
  trac(_y) = 1.;
  model.applyBC(BC::Neumann::FromTraction(trac), "air");

  model.addDumpField("external_force");
  model.addDumpField("internal_force");
  model.addDumpField("blocked_dofs");
  model.addDumpField("displacement");
  model.addDumpField("stress");
  model.addDumpField("grad_u");

  model.dump();
  auto & solver = model.getNonLinearSolver("static");
  solver.set("max_iterations", 1);
  solver.set("threshold", 1e-8);
  solver.set("convergence_type", SolveConvergenceCriteria::_residual);

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

  model.solveStep();
  // model.dump();

  std::map<std::string, Matrix<Real>> ref_strain;
  ref_strain["strong"] = Matrix<Real, dim, dim>::Zero();
  ref_strain["strong"](_y, _y) = .5;

  ref_strain["weak"] = Matrix<Real, dim, dim>::Zero();
  ref_strain["weak"](_y, _y) = 1;

  Matrix<Real, dim, dim> ref_stress = Matrix<Real, dim, dim>::Zero();
  ref_stress(_y, _y) = 1.;

  std::vector<std::string> mats = {"strong", "weak"};

  auto check = [](auto && view, const Matrix<Real> & ref) -> bool {
    for (auto && mat : view) {
      Real dist = (mat - ref).norm();
      if (dist > 1e-10) {
        return false;
      }
    }
    return true;
  };

  for (const auto & type : mesh.elementTypes(dim)) {
    for (auto mat_id : mats) {
      auto & stress = model.getMaterial(mat_id).getStress(type);
      auto & grad_u = model.getMaterial(mat_id).getGradU(type);

      auto sview = make_view<dim, dim>(stress);
      auto gview = make_view<dim, dim>(grad_u);

      if (not check(sview, ref_stress)) {
        AKANTU_ERROR("The stresses are not correct");
      }
      if (not check(gview, ref_strain[mat_id])) {
        AKANTU_ERROR("The grad_u are not correct");
      }
    }
  }

  return 0;
}
