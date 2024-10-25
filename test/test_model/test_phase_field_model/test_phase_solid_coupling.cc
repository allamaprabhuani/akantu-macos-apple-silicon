/**
 * Copyright (©) 2019-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_common.hh"
#include "material.hh"
#include "material_phasefield.hh"
#include "material_phasefield_anisotropic.hh"
#include "non_linear_solver.hh"
#include "non_linear_solver_newton_raphson.hh"
#include "phase_field_model.hh"
#include "solid_mechanics_model.hh"
#include "sparse_solver.hh"
#if defined(AKANTU_USE_PETSC)
#include "dof_manager_petsc.hh"
#endif
/* -------------------------------------------------------------------------- */
#include <fstream>
#include <iostream>
/* -------------------------------------------------------------------------- */

using namespace akantu;
const Int spatial_dimension = 2;

/* -------------------------------------------------------------------------- */
void applyDisplacement(SolidMechanicsModel &, Real &);
void computeStrainOnQuadPoints(SolidMechanicsModel &, PhaseFieldModel &,
                               GhostType);
void computeDamageOnQuadPoints(SolidMechanicsModel &, PhaseFieldModel &,
                               GhostType);
void gradUToEpsilon(const Matrix<Real> &, Matrix<Real> &);
/* -------------------------------------------------------------------------- */

int main(int argc, char * argv[]) {

  std::ofstream os("data.csv");
  os << "#strain stress damage analytical_sigma analytical_damage error_stress "
        "error_damage error_energy"
     << std::endl;

  initialize("material_coupling.dat", argc, argv);

  Mesh mesh(spatial_dimension);
  mesh.read("test_one_element.msh");

  SolidMechanicsModel model(mesh);
  auto & dof_manager = model.getDOFManager();
  model.initFull(_analysis_method = _static);

#if defined(AKANTU_USE_PETSC)
  if (aka::is_of_type<DOFManagerPETSc>(dof_manager)) {
    auto & solver = model.getNonLinearSolver();

    if (aka::is_of_type<NonLinearSolverNewtonRaphson>(solver)) {
      auto & sparse_solver =
          aka::as_type<NonLinearSolverNewtonRaphson>(solver).getSparseSolver();
      sparse_solver.set("pc_type", "cholesky");
      sparse_solver.set("ksp_rtol", "1e-30");
    }
  }
#endif

  PhaseFieldModel phase(mesh);
  auto && selector = std::make_shared<MeshDataPhaseFieldSelector<std::string>>(
      "physical_names", phase);
  phase.setPhaseFieldSelector(selector);
  phase.initFull(_analysis_method = _static);
  auto & dof_manager_phase_field = phase.getDOFManager();

#if defined(AKANTU_USE_PETSC)
  if (aka::is_of_type<DOFManagerPETSc>(dof_manager_phase_field)) {
    auto & solver = phase.getNonLinearSolver();

    if (aka::is_of_type<NonLinearSolverNewtonRaphson>(solver)) {
      auto & sparse_solver =
          aka::as_type<NonLinearSolverNewtonRaphson>(solver).getSparseSolver();
      sparse_solver.set("pc_type", "cholesky");
      sparse_solver.set("ksp_rtol", "1e-30");
    }
  }
#endif

  model.setBaseName("phase_solid");
  model.addDumpField("stress");
  model.addDumpField("grad_u");
  model.addDumpField("strain");
  model.addDumpFieldVector("displacement");
  model.addDumpField("damage");
  model.dump();

  Int nbSteps = 1000;
  Real increment = 1e-4;

  auto & stress = model.getMaterial(0).getArray<Real>("stress", _quadrangle_4);
  auto & damage = model.getMaterial(0).getArray<Real>("damage", _quadrangle_4);

  Real analytical_damage{0.};
  Real analytical_sigma{0.};
  Real analytical_energy{0.};

  auto & phasefield = phase.getPhaseField(0);

  const Real E = phasefield.getParam("E");
  const Real nu = phasefield.getParam("nu");
  Real c22 = E * (1 - nu) / ((1 + nu) * (1 - 2 * nu));

  const Real gc = phasefield.getParam("gc");
  const Real l0 = phasefield.getParam("l0");

  Real error_stress{0.};
  Real error_damage{0.};
  Real error_energy{0.};

  for (Int s = 0; s < nbSteps; ++s) {
    Real axial_strain = increment * s;
    applyDisplacement(model, axial_strain);

    model.solveStep("static");
    computeStrainOnQuadPoints(model, phase, _not_ghost);

    phase.solveStep();
    computeDamageOnQuadPoints(model, phase, _not_ghost);

    model.assembleInternalForces();

    analytical_damage = axial_strain * axial_strain * c22 /
                        (gc / l0 + axial_strain * axial_strain * c22);
    analytical_sigma =
        c22 * axial_strain * (1 - analytical_damage) * (1 - analytical_damage);

    analytical_energy = analytical_damage * analytical_damage * 0.5 * gc / l0;

    error_stress = std::abs(analytical_sigma - stress(0, 3)) / analytical_sigma;

    error_damage = std::abs(analytical_damage - damage(0)) / analytical_damage;

    error_energy =
        std::abs(analytical_energy - phase.getEnergy()) / analytical_energy;

    os << axial_strain << " " << stress(0, 3) << " " << damage(0) << " "
       << analytical_sigma << " " << analytical_damage << " " << error_stress
       << " " << error_damage << " " << error_energy << std::endl;

    if (error_damage > 1e-8 or error_stress > 1e-8) {
      std::cerr << std::left << std::setw(15) << "Step: " << s << std::endl;
      std::cerr << std::left << std::setw(15)
                << "Axial strain: " << axial_strain << std::endl;
      std::cerr << std::left << std::setw(15)
                << "Error damage: " << error_damage << std::endl;
      std::cerr << std::left << std::setw(15)
                << "Error stress: " << error_stress << std::endl;
      std::cerr << std::left << std::setw(15)
                << "Error energy: " << error_energy << std::endl;
      return EXIT_FAILURE;
    }

    model.dump();
  }

  os.close();
  finalize();

  return EXIT_SUCCESS;
}

/* -------------------------------------------------------------------------- */
void applyDisplacement(SolidMechanicsModel & model, Real & increment) {
  auto & displacement = model.getDisplacement();

  auto & positions = model.getMesh().getNodes();
  auto & blocked_dofs = model.getBlockedDOFs();

  for (Int n = 0; n < model.getMesh().getNbNodes(); ++n) {
    if (positions(n, 1) == -0.5) {
      displacement(n, 0) = 0;
      displacement(n, 1) = 0;
      blocked_dofs(n, 0) = true;
      blocked_dofs(n, 1) = true;
    } else {
      displacement(n, 0) = 0;
      displacement(n, 1) = increment;
      blocked_dofs(n, 0) = true;
      blocked_dofs(n, 1) = true;
    }
  }
}

/* -------------------------------------------------------------------------- */
void computeStrainOnQuadPoints(SolidMechanicsModel & solid,
                               PhaseFieldModel & phase, GhostType ghost_type) {

  auto & mesh = solid.getMesh();

  auto nb_materials = solid.getNbConstitutiveLaws();
  auto nb_phasefields = phase.getNbConstitutiveLaws();

  AKANTU_DEBUG_ASSERT(
      nb_phasefields == nb_materials,
      "The number of phasefields and materials should be equal");

  for (auto index : arange(nb_materials)) {
    auto & material = solid.getConstitutiveLaw(index);

    for (auto index2 : arange(nb_phasefields)) {
      auto & phasefield = phase.getConstitutiveLaw(index2);

      if (phasefield.getName() == material.getName()) {

        auto & strain_on_qpoints = phasefield.getStrain();
        auto & gradu_on_qpoints = material.getGradU();

        for (const auto & type :
             mesh.elementTypes(spatial_dimension, ghost_type)) {
          auto & strain_on_qpoints_vect = strain_on_qpoints(type, ghost_type);
          auto & gradu_on_qpoints_vect = gradu_on_qpoints(type, ghost_type);
          for (auto && values :
               zip(make_view(strain_on_qpoints_vect, spatial_dimension,
                             spatial_dimension),
                   make_view(gradu_on_qpoints_vect, spatial_dimension,
                             spatial_dimension))) {
            auto & strain = std::get<0>(values);
            auto & grad_u = std::get<1>(values);
            Material::gradUToEpsilon<spatial_dimension>(grad_u, strain);
          }
        }

        break;
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
void computeDamageOnQuadPoints(SolidMechanicsModel & solid,
                               PhaseFieldModel & phase, GhostType ghost_type) {

  auto & fem = phase.getFEEngine();
  auto & mesh = phase.getMesh();

  auto nb_materials = solid.getNbConstitutiveLaws();
  auto nb_phasefields = phase.getNbConstitutiveLaws();

  AKANTU_DEBUG_ASSERT(
      nb_phasefields == nb_materials,
      "The number of phasefields and materials should be equal");

  for (auto index : arange(nb_materials)) {
    auto & material = solid.getConstitutiveLaw(index);

    for (auto index2 : arange(nb_phasefields)) {
      auto & phasefield = phase.getConstitutiveLaw(index2);

      if (phasefield.getName() == material.getName()) {

        switch (spatial_dimension) {
        case 1: {
          auto & mat = dynamic_cast<MaterialDamage<1> &>(material);
          auto & solid_damage = mat.getDamage();

          for (const auto & type :
               mesh.elementTypes(spatial_dimension, ghost_type)) {
            auto & damage_on_qpoints_vect = solid_damage(type, ghost_type);

            fem.interpolateOnIntegrationPoints(
                phase.getDamage(), damage_on_qpoints_vect, 1, type, ghost_type);
          }

          break;
        }
        case 2: {
          auto & mat = dynamic_cast<MaterialDamage<2> &>(material);
          auto & solid_damage = mat.getDamage();

          for (const auto & type :
               mesh.elementTypes(spatial_dimension, ghost_type)) {
            auto & damage_on_qpoints_vect = solid_damage(type, ghost_type);

            fem.interpolateOnIntegrationPoints(
                phase.getDamage(), damage_on_qpoints_vect, 1, type, ghost_type);
          }
          break;
        }
        default:
          break;
        }
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
void gradUToEpsilon(const Matrix<Real> & grad_u, Matrix<Real> & epsilon) {
  for (Int i = 0; i < spatial_dimension; ++i) {
    for (Int j = 0; j < spatial_dimension; ++j)
      epsilon(i, j) = 0.5 * (grad_u(i, j) + grad_u(j, i));
  }
}
