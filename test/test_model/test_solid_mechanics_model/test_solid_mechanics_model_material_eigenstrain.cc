/**
 * Copyright (©) 2011-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "mesh_utils.hh"
#include "non_linear_solver.hh"
#include "non_linear_solver_newton_raphson.hh"
#include "solid_mechanics_model.hh"
#if defined(AKANTU_USE_PETSC)
#include "sparse_solver_petsc.hh"
#endif

/* -------------------------------------------------------------------------- */

using namespace akantu;

Matrix<Real, 3, 4> alpha{{0.01, 0.02, 0.03, 0.04},
                         {0.05, 0.06, 0.07, 0.08},
                         {0.09, 0.10, 0.11, 0.12}};

/* -------------------------------------------------------------------------- */
template <ElementType type> static Matrix<Real> prescribed_strain() {
  Int spatial_dimension = ElementClass<type>::getSpatialDimension();
  Matrix<Real> strain(spatial_dimension, spatial_dimension);

  strain = alpha.block(0, 1, spatial_dimension, spatial_dimension);

  return strain;
}

template <ElementType type>
static Matrix<Real> prescribed_stress(Matrix<Real> prescribed_eigengradu) {
  Int spatial_dimension = ElementClass<type>::getSpatialDimension();
  Matrix<Real> stress(spatial_dimension, spatial_dimension);

  // plane strain in 2d
  Matrix<Real> strain(spatial_dimension, spatial_dimension);
  Matrix<Real> pstrain;
  pstrain = prescribed_strain<type>();
  Real nu = 0.3;
  Real E = 2.1e11;
  Real trace = 0;

  /// symetric part of the strain tensor
  strain = (pstrain + pstrain.transpose()) / 2.;

  // elastic strain is equal to elastic strain minus the eigenstrain
  strain -= prescribed_eigengradu;
  trace = strain.trace();

  Real lambda = nu * E / ((1 + nu) * (1 - 2 * nu));
  Real mu = E / (2 * (1 + nu));

  if (spatial_dimension == 1) {
    stress = E * strain;
  } else {
    auto I = Matrix<Real>::Identity(spatial_dimension, spatial_dimension);

    stress = I * trace * lambda + 2. * mu * strain;
  }

  return stress;
}

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
int main(int argc, char * argv[]) {
  initialize("material_elastic_plane_strain.dat", argc, argv);

  Int dim = 3;
  const ElementType element_type = _tetrahedron_4;

  Matrix<Real> prescribed_eigengradu(dim, dim);
  prescribed_eigengradu.set(0.1);

  /// load mesh
  Mesh mesh(dim);
  mesh.read("cube_3d_tet_4.msh");

  /// declaration of model
  SolidMechanicsModel model(mesh);
  /// model initialization
  model.initFull(_analysis_method = _static);

  // model.getNewSolver("static", TimeStepSolverType::_static,
  // NonLinearSolverType::_newton_raphson_modified);
  auto & solver = model.getNonLinearSolver("static");
  solver.set("threshold", 2e-4);
  solver.set("max_iterations", 2);
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

  const Array<Real> & coordinates = mesh.getNodes();
  Array<Real> & displacement = model.getDisplacement();
  Array<bool> & boundary = model.getBlockedDOFs();
  MeshUtils::buildFacets(mesh);

  mesh.createBoundaryGroupFromGeometry();

  // Loop over (Sub)Boundar(ies)
  for (auto & group : mesh.iterateElementGroups()) {
    for (const auto & n : group.getNodeGroup()) {
      std::cout << "Node " << n << std::endl;
      for (Int i = 0; i < dim; ++i) {
        displacement(n, i) = alpha(i, 0);
        for (Int j = 0; j < dim; ++j) {
          displacement(n, i) += alpha(i, j + 1) * coordinates(n, j);
        }
        boundary(n, i) = true;
      }
    }
  }

  /* ------------------------------------------------------------------------ */
  /* Apply eigenstrain in each element */
  /* ------------------------------------------------------------------------ */
  Array<Real> & eigengradu_vect =
      model.getMaterial(0).getInternal<Real>("eigen_grad_u")(element_type);
  auto eigengradu_it = eigengradu_vect.begin(dim, dim);
  auto eigengradu_end = eigengradu_vect.end(dim, dim);

  for (; eigengradu_it != eigengradu_end; ++eigengradu_it) {
    *eigengradu_it = prescribed_eigengradu;
  }

  /* ------------------------------------------------------------------------ */
  /* Static solve                                                             */
  /* ------------------------------------------------------------------------ */
  model.solveStep();

  std::cout << "Converged in " << Int(solver.get("nb_iterations")) << " ("
            << Real(solver.get("error")) << ")" << std::endl;

  /* ------------------------------------------------------------------------ */
  /* Checks                                                                   */
  /* ------------------------------------------------------------------------ */
  const Array<Real> & stress_vect =
      model.getMaterial(0).getStress(element_type);

  auto stress_it = stress_vect.begin(dim, dim);
  auto stress_end = stress_vect.end(dim, dim);

  Matrix<Real> presc_stress;
  presc_stress = prescribed_stress<element_type>(prescribed_eigengradu);

  Real stress_tolerance = 1e-13;

  for (; stress_it != stress_end; ++stress_it) {
    const auto & stress = *stress_it;
    Matrix<Real> diff(dim, dim);

    diff = stress;
    diff -= presc_stress;
    Real stress_error =
        diff.lpNorm<Eigen::Infinity>() / stress.lpNorm<Eigen::Infinity>();

    if (stress_error > stress_tolerance) {
      std::cerr << "stress error: " << stress_error << " > " << stress_tolerance
                << std::endl;
      std::cerr << "stress: " << stress << std::endl
                << "prescribed stress: " << presc_stress << std::endl;
      return EXIT_FAILURE;
    } // else {
    //   std::cerr << "stress error: " << stress_error
    //             << " < " << stress_tolerance << std::endl;
    // }
  }

  finalize();

  return EXIT_SUCCESS;
}
