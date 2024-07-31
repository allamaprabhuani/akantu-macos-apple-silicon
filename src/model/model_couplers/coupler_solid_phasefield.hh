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
#include "boundary_condition.hh"
#include "data_accessor.hh"
#include "fe_engine.hh"
#include "material.hh"
#include "material_phasefield.hh"
#include "model.hh"
#include "phase_field_model.hh"
#include "phasefield.hh"
#include "solid_mechanics_model.hh"
#include "sparse_matrix.hh"
#include "time_step_solver.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_COUPLER_SOLID_PHASEFIELD_HH_
#define AKANTU_COUPLER_SOLID_PHASEFIELD_HH_

/* ------------------------------------------------------------------------ */
/* Coupling : Solid Mechanics / PhaseField                                  */
/* ------------------------------------------------------------------------ */
namespace akantu {
template <ElementKind kind, class IntegrationOrderFunctor>
class IntegratorGauss;
template <ElementKind kind> class ShapeLagrange;
class DOFManager;
} // namespace akantu

namespace akantu {

class CouplerSolidPhaseField
    : public Model,
      public DataAccessor<Element>,
      public DataAccessor<Idx>,
      public BoundaryCondition<CouplerSolidPhaseField> {

  /* ------------------------------------------------------------------------ */
  /*  Constructor/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

  using MyFEEngineType = FEEngineTemplate<IntegratorGauss, ShapeLagrange>;

public:
  CouplerSolidPhaseField(
      Mesh & mesh, Int dim = _all_dimensions,
      const ID & id = "coupler_solid_phasefield",
      ModelType model_type = ModelType::_coupler_solid_phasefield);

  ~CouplerSolidPhaseField() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
protected:
  /// initialize the complete model
  void initFullImpl(const ModelOptions & options) override;

  /// initialize the modelType
  void initModel() override;

  /// get some default values for derived classes
  std::tuple<ID, TimeStepSolverType>
  getDefaultSolverID(const AnalysisMethod & method) override;

  /* ------------------------------------------------------------------------ */
  /* Solver Interface                                                         */
  /* ------------------------------------------------------------------------ */
public:
  /// assembles the contact stiffness matrix
  virtual void assembleStiffnessMatrix();

  /// assembles the contant internal forces
  virtual void assembleInternalForces();

public:
  /// computes damage on quad points for solid mechanics model from
  /// damage array from phasefield model
  void computeDamageOnQuadPoints(GhostType ghost_type);

  /// computes strain on quadrature points for phasefield model from
  /// displacement gradient from solid mechanics model
  void computeStrainOnQuadPoints(GhostType ghost_type);

  /// solve the coupled model
  void solve(const ID & solid_solver_id = "", const ID & phase_solver_id = "");

private:
  /// test the convergence criteria
  bool checkConvergence(Array<Real> & /*u_new*/, Array<Real> & /*u_old*/,
                        Array<Real> & /*d_new*/, Array<Real> & /*d_old*/);

protected:
  /// callback for the solver, this adds f_{ext} - f_{int} to the residual
  void assembleResidual() override;

  /// callback for the solver, this adds f_{ext} or  f_{int} to the residual
  void assembleResidual(const ID & residual_part) override;
  bool canSplitResidual() const override { return true; }

  /// get the type of matrix needed
  MatrixType getMatrixType(const ID & matrix_id) const override;

  /// callback for the solver, this assembles different matrices
  void assembleMatrix(const ID & matrix_id) override;

  /// callback for the solver, this assembles the stiffness matrix
  void assembleLumpedMatrix(const ID & matrix_id) override;

  /// callback for the model to instantiate the matricess when needed
  void initSolver(TimeStepSolverType /*time_step_solver_type*/,
                  NonLinearSolverType /*non_linear_solver_type*/) override;

  /// callback for the solver, this is called at beginning of solve
  void predictor() override;

  /// callback for the solver, this is called at end of solve
  void corrector() override;

  /// callback for the solver, this is called at beginning of solve
  void beforeSolveStep() override;

  /// callback for the solver, this is called at end of solve
  void afterSolveStep(bool converged = true) override;

  /// solve the coupled model
  // void solveStep(const ID & solver_id = "") override;

  /// solve a step using a given pre instantiated time step solver and
  /// non linear solver with a user defined callback instead of the
  /// model itself /!\ This can mess up everything
  // void solveStep(SolverCallback & callback, const ID & solver_id = "")
  // override;

  /* ------------------------------------------------------------------------ */
  /* Mass matrix for solid mechanics model                                    */
  /* ------------------------------------------------------------------------ */
public:
  /// assemble the lumped mass matrix
  void assembleMassLumped();

  /// assemble the mass matrix for consistent mass resolutions
  void assembleMass();

protected:
  /// assemble the lumped mass matrix for local and ghost elements
  void assembleMassLumped(GhostType ghost_type);

  /// assemble the mass matrix for either _ghost or _not_ghost elements
  void assembleMass(GhostType ghost_type);

protected:
  /* --------------------------------------------------------------------------
   */
  TimeStepSolverType getDefaultSolverType() const override;
  /* --------------------------------------------------------------------------
   */
  ModelSolverOptions
  getDefaultSolverOptions(const TimeStepSolverType & type) const override;

public:
  bool isDefaultSolverExplicit() { return method == _explicit_lumped_mass; }

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  FEEngine & getFEEngineBoundary(const ID & name = "") override;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get the solid mechanics model
  AKANTU_GET_MACRO(SolidMechanicsModel, *solid, SolidMechanicsModel &);

  /// get the contact mechanics model
  AKANTU_GET_MACRO(PhaseFieldModel, *phase, PhaseFieldModel &);

  /* ------------------------------------------------------------------------ */
  /* Dumpable interface                                                       */
  /* ------------------------------------------------------------------------ */
public:
  std::shared_ptr<dumpers::Field>
  createNodalFieldReal(const std::string & field_name,
                       const std::string & group_name,
                       bool padding_flag) override;

  std::shared_ptr<dumpers::Field>
  createNodalFieldBool(const std::string & field_name,
                       const std::string & group_name,
                       bool padding_flag) override;

  std::shared_ptr<dumpers::Field>
  createElementalField(const std::string & field_name,
                       const std::string & group_name, bool padding_flag,
                       Int spatial_dimension, ElementKind kind) override;

  void dump(const std::string & dumper_name) override;
  void dump(const std::string & dumper_name, Int step) override;
  void dump(const std::string & dumper_name, Real time, Int step) override;

  void dump() override;
  void dump(Int step) override;
  void dump(Real time, Int step) override;

  /* ------------------------------------------------------------------------ */
  /* Members                                                                  */
  /* ------------------------------------------------------------------------ */
private:
  /// solid mechanics model
  std::unique_ptr<SolidMechanicsModel> solid;

  /// phasefield model
  std::unique_ptr<PhaseFieldModel> phase;

  Array<Real> * displacement{nullptr};

  ///
  Array<Real> * displacement_increment{nullptr};

  /// external forces array
  Array<Real> * external_force{nullptr};

  Array<Real> initial_mass;
  bool mass_initialized{false};
};

} // namespace akantu

#endif /* AKANTU_COUPLER_SOLID_PHASEFIELD_HH_ */
