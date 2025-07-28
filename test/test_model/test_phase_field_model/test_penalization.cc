#include "aka_common.hh"
#include "coupler_solid_phasefield.hh"
#include "material.hh"
#include "material_phasefield.hh"
#include "non_linear_solver.hh"
#include "phase_field_model.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */
#include <cmath>
#include <fstream>
#include <iostream>
/* -------------------------------------------------------------------------- */

using namespace akantu;
const UInt spatial_dimension = 2;

/* -------------------------------------------------------------------------- */
void applyDisplacement(SolidMechanicsModel &, Real &);
/* -------------------------------------------------------------------------- */

int main(int argc, char * argv[]) {

  std::ofstream os("data.csv");
  os << "#strain stress damage analytical_sigma analytical_damage" << std::endl;

  initialize("material_penalization.dat", argc, argv);

  Mesh mesh(spatial_dimension);
  mesh.read("test_four_elements.msh");

  CouplerSolidPhaseField coupler(mesh);
  auto & model = coupler.getSolidMechanicsModel();
  auto & phase = coupler.getPhaseFieldModel();

  model.initFull(_analysis_method = _static);
  auto & solver = model.getNonLinearSolver("static");
  solver.set("max_iterations", 1000);
  solver.set("threshold", 1e-8);
  solver.set("convergence_type", SolveConvergenceCriteria::_residual);

  auto && selector = std::make_shared<MeshDataPhaseFieldSelector<std::string>>(
      "physical_names", phase);
  phase.setPhaseFieldSelector(selector);
  phase.initFull(_analysis_method = _static);
  auto & solver_phase = phase.getNonLinearSolver("static");
  solver_phase.set("max_iterations", 1000);
  solver_phase.set("threshold", 1e-6);
  solver_phase.set("convergence_type", SolveConvergenceCriteria::_solution);

  model.setBaseName("phase_solid");
  model.addDumpField("stress");
  model.addDumpField("grad_u");
  model.addDumpFieldVector("displacement");
  model.addDumpField("damage");
  model.dump();

  UInt nbSteps = 1000;
  Real increment = 1e-5;

  auto & stress = model.getMaterial(0).getArray<Real>("stress", _quadrangle_4);
  auto & damage = model.getMaterial(0).getArray<Real>("damage", _quadrangle_4);

  Real new_damage{0.};
  Real analytical_damage{0.};
  Real analytical_sigma{0.};
  Real axial_strain{0.};

  auto & phasefield = phase.getPhaseField(0);

  const Real E = phasefield.getParam("E");
  const Real nu = phasefield.getParam("nu");
  Real c22 = E * (1 - nu) / ((1 + nu) * (1 - 2 * nu));

  const Real gc = phasefield.getParam("gc");
  const Real l0 = phasefield.getParam("l0");

  Real error_stress{0.};

  Real error_damage{0.};

  for (UInt s = 0; s < nbSteps; ++s) {
    if (s < 500) {
      axial_strain += increment;
    } else {
      axial_strain -= increment;
    }
    applyDisplacement(model, axial_strain);

    coupler.solve("static", "static");
    phase.savePreviousState();

    new_damage = axial_strain * axial_strain * c22 /
                 (gc / l0 + axial_strain * axial_strain * c22);
    analytical_damage = std::max(analytical_damage, new_damage);
    analytical_sigma =
        c22 * axial_strain * (1 - analytical_damage) * (1 - analytical_damage);

    error_stress =
        std::abs(analytical_sigma - stress(0, 3)) / std::abs(analytical_sigma);

    error_damage = std::abs(analytical_damage - damage(0)) / analytical_damage;

    os << axial_strain << " " << stress(0, 3) << " " << damage(0) << " "
       << analytical_sigma << " " << analytical_damage << " " << error_stress
       << " " << error_damage << std::endl;

    if ((error_damage > 1e-3 or error_stress > 1e-3) and
        std::abs(axial_strain) > 1e-13) {
      std::cerr << std::left << std::setw(15) << "Step: " << s << std::endl;
      std::cerr << std::left << std::setw(15)
                << "Axial strain: " << axial_strain << std::endl;
      std::cerr << std::left << std::setw(15)
                << "An. damage: " << analytical_damage << std::endl;
      std::cerr << std::left << std::setw(15) << "Damage: " << damage(0)
                << std::endl;
      std::cerr << std::left << std::setw(15)
                << "Error damage: " << error_damage << std::endl;
      std::cerr << std::left << std::setw(15)
                << "Error stress: " << error_stress << std::endl;
      return EXIT_FAILURE;
    }

    model.dump();
  }

  os.close();
  finalize();

  return EXIT_SUCCESS;
}

/* --------------------------------------------------------------------------
 */
void applyDisplacement(SolidMechanicsModel & model, Real & increment) {
  auto & displacement = model.getDisplacement();

  auto & positions = model.getMesh().getNodes();
  auto & blocked_dofs = model.getBlockedDOFs();

  for (Idx n = 0; n < model.getMesh().getNbNodes(); ++n) {
    if (positions(n, 1) == -0.5) {
      displacement(n, 0) = 0;
      displacement(n, 1) = 0;
      blocked_dofs(n, 0) = true;
      blocked_dofs(n, 1) = true;
    } else if (positions(n, 1) == 0.5) {
      displacement(n, 0) = 0;
      displacement(n, 1) = increment;
      blocked_dofs(n, 0) = true;
      blocked_dofs(n, 1) = true;
    } else {
      displacement(n, 0) = 0;
      blocked_dofs(n, 0) = true;
    }
  }
}

/* --------------------------------------------------------------------------
 */
