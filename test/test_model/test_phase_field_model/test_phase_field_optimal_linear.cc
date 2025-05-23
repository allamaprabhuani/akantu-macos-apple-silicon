#include "aka_common.hh"
#include "coupler_solid_phasefield.hh"
#include "material.hh"
#include "material_phasefield.hh"
#include "non_linear_solver.hh"
#include "non_linear_solver_tao.hh"
#include "phase_field_model.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */
#include <cmath>
#include <fstream>
#include <iostream>
#include <ostream>
/* -------------------------------------------------------------------------- */

using namespace akantu;
const UInt spatial_dimension = 2;

/* -------------------------------------------------------------------------- */
Real computeOptimalDamage(Real const & x, Real const & l) {
  Real optimal_damage = std::pow((1. - std::abs(x) / (2. * l)), 2.) *
                        (x > -2. * l) * (x < 2. * l);
  return optimal_damage;
}
/* -------------------------------------------------------------------------- */

int main(int argc, char * argv[]) {

  std::ofstream os("data.csv");
  os << "position,damage,optimal_damage" << std::endl;

  initialize("material_optimal_linear.dat", argc, argv);

  Mesh mesh(spatial_dimension);
  mesh.read("test_optimal_damage.msh");
  // mesh.read("test_one_element.msh");

  PhaseFieldModel phase(mesh);

  const auto solver_type = NonLinearSolverType::_petsc_tao;
  if constexpr (solver_type == NonLinearSolverType::_petsc_snes) {
    phase.initDOFManager("petsc");
    phase.initFull(
        _analysis_method = _static,
        _solver_options = ModelSolverOptions{
            .non_linear_solver_type = NonLinearSolverType::_petsc_snes,
            .sparse_solver_type = SparseSolverType::_petsc});
    auto & solver_phase = phase.getNonLinearSolver("static");
    solver_phase.set("snes_max_it", 500);

    solver_phase.set("snes_rtol", 1e-10);
    solver_phase.set("snes_atol", 1e-6);
    solver_phase.set("snes_divergence_tolerance", 1e+10);
    solver_phase.set("pc_type", "lu");
    solver_phase.set("snes_linesearch_type", "none");
    solver_phase.set("snes_monitor", "stdout");
  } else if constexpr (solver_type == NonLinearSolverType::_petsc_tao) {
    phase.initDOFManager("petsc");
    phase.initFull(
        _analysis_method = _static,
        _solver_options = ModelSolverOptions{
            .non_linear_solver_type = NonLinearSolverType::_petsc_tao,
            .sparse_solver_type = SparseSolverType::_petsc});
    auto & tao_solver =
        aka::as_type<NonLinearSolverTAO>(phase.getNonLinearSolver());
    tao_solver.set("tao_type", "gpcg");
    tao_solver.set("tao_gatol", 1e-9);
    tao_solver.set("tao_monitor", "stdout");

    auto & dof_manager = aka::as_type<DOFManagerPETSc>(phase.getDOFManager());
    SolverVectorPETSc lower_bound(dof_manager, "lower_bound");
    SolverVectorPETSc upper_bound(dof_manager, "upper_bound");
    lower_bound.resize();
    upper_bound.resize();

    lower_bound.set(0.);
    upper_bound.set(1.);
    tao_solver.setBounds(lower_bound, upper_bound);
  } else {
    phase.initFull(_analysis_method = _static);
    auto & solver_phase = phase.getNonLinearSolver("static");
    solver_phase.set("max_iterations", 1000);
    solver_phase.set("threshold", 1e-6);
    solver_phase.set("convergence_type", SolveConvergenceCriteria::_residual);
    solver_phase.set("verbose_iteration", true);
  }

  auto & phasefield = phase.getPhaseField(0);
  auto & damage = phase.getDamage();

  Real analytical_damage{0.};

  // const Real gc = phasefield.getParam("gc");
  const Real l0 = phasefield.getParam("l0");
  // const Real L = 1000.;

  Real error_damage{0.};

  auto & positions = phase.getMesh().getNodes();
  //  auto & blocked_dofs = phase.getBlockedDOFs();

  phase.applyBC(BC::Dirichlet::FixedValue(1., _x), "blocked");

  phase.solveStep("static");
  phase.savePreviousState();

  for (Int n = 0; n < phase.getMesh().getNbNodes(); ++n) {
    // if (positions(n, 0) == 0.6) {
    analytical_damage = computeOptimalDamage(positions(n, 0), l0);
    // error_damage = std::abs(analytical_damage - damage(n)) /
    // analytical_damage;
    error_damage = std::abs(analytical_damage - damage(n));
    os << positions(n, 0) << "," << damage(n) << "," << analytical_damage
       << std::endl;
    if (error_damage > 1e-5) {
      std::cerr << "Position: " << positions(n, 0) << std::endl;
      std::cerr << "Optimal damage: " << analytical_damage << std::endl;
      std::cerr << "Damage: " << damage(n) << std::endl;
      return EXIT_FAILURE;
    }
    // }
  }

  os.close();
  finalize();

  return EXIT_SUCCESS;
}
