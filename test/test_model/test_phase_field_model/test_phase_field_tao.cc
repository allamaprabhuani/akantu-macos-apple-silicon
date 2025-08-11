/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "material.hh"
#include "material_phasefield.hh"
#include "non_linear_solver.hh"
#include "non_linear_solver_tao.hh"
#include "phase_field_model.hh"
#include "solid_mechanics_model.hh"
#include "solver_vector_petsc.hh"
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

  initialize("material_linear_tao.dat", argc, argv);

  Mesh mesh(spatial_dimension);
  mesh.read("test_one_element.msh");

  SolidMechanicsModel model(mesh);
  model.initFull(_analysis_method = _static);

  auto & solver = model.getNonLinearSolver("static");
  solver.set("max_iterations", 100);
  solver.set("threshold", 1e-9);
  solver.set("convergence_type", SolveConvergenceCriteria::_residual);

  PhaseFieldModel phase(mesh);
  auto && selector = std::make_shared<MeshDataPhaseFieldSelector<std::string>>(
      "physical_names", phase);
  phase.setPhaseFieldSelector(selector);

  PetscOptionsSetValue(NULL, "-tao_type", "gpcg");
  PetscOptionsSetValue(NULL, "-tao_gatol", "1e-9");
  phase.initDOFManager("petsc");
  auto & DOFManager = aka::as_type<DOFManagerPETSc>(phase.getDOFManager());
  phase.initFull(_analysis_method = _static,
                 _solver_options = ModelSolverOptions{
                     .non_linear_solver_type = NonLinearSolverType::_petsc_tao,
                     .sparse_solver_type = SparseSolverType::_petsc});
  auto & tao_solver =
      aka::as_type<NonLinearSolverTAO>(phase.getNonLinearSolver());

  model.setBaseName("phase_solid");
  model.addDumpField("stress");
  model.addDumpField("grad_u");
  model.addDumpField("strain");
  model.addDumpFieldVector("displacement");
  model.addDumpField("damage");
  model.dump();

  Int nbSteps = 1000;
  Real increment = 1e-5;

  auto & stress = model.getMaterial(0).getArray<Real>("stress", _quadrangle_4);
  auto & damage = model.getMaterial(0).getArray<Real>("damage", _quadrangle_4);

  SolverVectorPETSc lower_bound(DOFManager, "lower_bound");
  lower_bound.resize();
  lower_bound.set(0.);
  SolverVectorPETSc upper_bound(DOFManager, "upper_bound");
  upper_bound.resize();
  upper_bound.set(1.);
  tao_solver.setBounds(lower_bound, upper_bound);

  auto & dof_damage = DOFManager.getDOFs("damage");

  Real dam{0.};
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
  Real tol = 1e-8;

  for (Int s = 0; s < nbSteps; ++s) {
    Real axial_strain{0.};
    if (s < 500) {
      axial_strain = increment * s;
    } else {
      axial_strain = (1000 - double(s)) * increment;
    }

    applyDisplacement(model, axial_strain);

    model.solveStep("static");
    computeStrainOnQuadPoints(model, phase, _not_ghost);

    phase.solveStep();
    computeDamageOnQuadPoints(model, phase, _not_ghost);

    model.assembleInternalForces();

    lower_bound.set(0.);
    upper_bound.set(1.);
    DOFManager.assembleToGlobalArray("damage", dof_damage, upper_bound, -1.);
    tao_solver.setBounds(lower_bound, upper_bound);

    dam = 1. - 3. * gc / (8. * l0 * axial_strain * axial_strain * c22);
    analytical_damage = std::max(dam, analytical_damage);
    analytical_sigma =
        c22 * axial_strain * (1 - analytical_damage) * (1 - analytical_damage);

    analytical_energy = analytical_damage * 3. * gc / (8. * l0);

    error_stress = std::abs(analytical_sigma - stress(0, 3)) / analytical_sigma;

    error_damage = std::abs(analytical_damage - damage(0));

    error_energy =
        std::abs(analytical_energy - phase.getEnergy()) / analytical_energy;

    os << axial_strain << " " << stress(0, 3) << " " << damage(0) << " "
       << analytical_sigma << " " << analytical_damage << " " << error_stress
       << " " << error_damage << " " << error_energy << std::endl;

    if (analytical_damage < 1e-8) {
      tol = 1e-5;
    } else {
      tol = 1e-8;
    }
    if ((error_damage > tol or error_stress > tol) and s > 0) {
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

  auto nb_materials = solid.getNbMaterials();
  auto nb_phasefields = phase.getNbPhaseFields();

  AKANTU_DEBUG_ASSERT(
      nb_phasefields == nb_materials,
      "The number of phasefields and materials should be equal");

  for (auto index : arange(nb_materials)) {
    auto & material = solid.getMaterial(index);

    for (auto index2 : arange(nb_phasefields)) {
      auto & phasefield = phase.getPhaseField(index2);

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

  auto nb_materials = solid.getNbMaterials();
  auto nb_phasefields = phase.getNbPhaseFields();

  AKANTU_DEBUG_ASSERT(
      nb_phasefields == nb_materials,
      "The number of phasefields and materials should be equal");

  for (auto index : arange(nb_materials)) {
    auto & material = solid.getMaterial(index);

    for (auto index2 : arange(nb_phasefields)) {
      auto & phasefield = phase.getPhaseField(index2);

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
