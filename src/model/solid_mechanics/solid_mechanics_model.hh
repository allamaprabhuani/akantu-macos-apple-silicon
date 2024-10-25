/**
 * Copyright (©) 2010-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "constitutive_laws_handler.hh"
#include "fe_engine.hh"
#include "material.hh"
#include "model.hh"
#include "solid_mechanics_model_event_handler.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_SOLID_MECHANICS_MODEL_HH_
#define AKANTU_SOLID_MECHANICS_MODEL_HH_

namespace akantu {
class ConstitutiveLawSelector;
class DumperIOHelper;
template <ElementKind kind, class IntegrationOrderFunctor>
class IntegratorGauss;
template <ElementKind kind> class ShapeLagrange;
} // namespace akantu

/* -------------------------------------------------------------------------- */
namespace akantu {

/* -------------------------------------------------------------------------- */
class SolidMechanicsModel
    : public ConstitutiveLawsHandler<Material, Model>,
      public BoundaryCondition<SolidMechanicsModel>,
      public EventHandlerManager<SolidMechanicsModelEventHandler> {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  using MyFEEngineType = FEEngineTemplate<IntegratorGauss, ShapeLagrange>;

protected:
  using EventManager = EventHandlerManager<SolidMechanicsModelEventHandler>;
  using CLHParent = ConstitutiveLawsHandler<Material, Model>;

public:
  SolidMechanicsModel(Mesh & mesh, Int dim = _all_dimensions,
                      const ID & id = "solid_mechanics_model",
                      const std::shared_ptr<DOFManager> & dof_manager = nullptr,
                      ModelType model_type = ModelType::_solid_mechanics_model);

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
protected:
  /// initialize completely the model
  void initFullImpl(
      const ModelOptions & options = SolidMechanicsModelOptions()) override;

  void instantiateMaterials();

public:
  /// function to print the containt of the class
  void printself(std::ostream & stream, int indent = 0) const override;

  /// get some default values for derived classes
  [[nodiscard]] std::tuple<ID, TimeStepSolverType>
  getDefaultSolverID(const AnalysisMethod & method) override;

  /* ------------------------------------------------------------------------ */
  /* Solver interface                                                         */
  /* ------------------------------------------------------------------------ */
public:
  /// assembles the stiffness matrix,
  virtual void assembleStiffnessMatrix(bool need_to_reassemble = false);
  /// assembles the internal forces in the array internal_forces
  virtual void assembleInternalForces();

protected:
  /// callback for the solver, this adds f_{ext} - f_{int} to the residual
  void assembleResidual() override;

  /// callback for the solver, this adds f_{ext} or  f_{int} to the residual
  void assembleResidual(const ID & residual_part) override;
  [[nodiscard]] bool canSplitResidual() const override { return true; }

  /// get the type of matrix needed
  [[nodiscard]] MatrixType getMatrixType(const ID & matrix_id) const override;

  /// callback for the solver, this assembles different matrices
  void assembleMatrix(const ID & matrix_id) override;

  /// callback for the solver, this assembles the stiffness matrix
  void assembleLumpedMatrix(const ID & matrix_id) override;

  /// callback for the solver, this is called at beginning of solve
  void predictor() override;

  /// callback for the solver, this is called at end of solve
  void corrector() override;

  /// callback for the solver, this is called at beginning of solve
  void beforeSolveStep() override;
  /// callback for the solver, this is called at end of solve
  void afterSolveStep(bool converged = true) override;

  /// Callback for the model to instantiate the matricees when needed
  void initSolver(TimeStepSolverType time_step_solver_type) override;

public:
  /* ------------------------------------------------------------------------ */
  [[nodiscard]] TimeStepSolverType getDefaultSolverType() const override;

  [[nodiscard]] ModelSolverOptions
  getDefaultSolverOptions(const TimeStepSolverType & type) const override;

  bool isDefaultSolverExplicit() {
    return method == _explicit_lumped_mass ||
           method == _explicit_consistent_mass;
  }

protected:
  /// update the current position vector
  void updateCurrentPosition();

  /* ------------------------------------------------------------------------ */
  /* Materials (solid_mechanics_model_material.cc)                            */
  /* ------------------------------------------------------------------------ */
public:
  /// apply a constant eigen_grad_u on all quadrature points of a given material
  virtual void applyEigenGradU(const Matrix<Real> & prescribed_eigen_grad_u,
                               const ID & material_name,
                               GhostType ghost_type = _not_ghost);

  void computeNonLocalContribution(GhostType ghost_type) override;

  /* ------------------------------------------------------------------------ */
  /* Mass (solid_mechanics_model_mass.cc)                                     */
  /* ------------------------------------------------------------------------ */
public:
  /// force reassembly of the lumped mass matrix
  void forceReassembleMassLumped();

  /// assemble the lumped mass matrix
  void assembleMassLumped();

  /// assemble the mass matrix for consistent mass resolutions
  void assembleMass();

  /// assemble the lumped mass matrix for local and ghost elements
  void assembleMassLumped(GhostType ghost_type);

  /// assemble the mass matrix for either _ghost or _not_ghost elements
  void assembleMass(GhostType ghost_type);

  /// fill a vector of rho
  void computeRho(Array<Real> & rho, ElementType type, GhostType ghost_type);

  /* ------------------------------------------------------------------------ */
  /* Energies                                                                 */
  /* ------------------------------------------------------------------------ */
protected:
  /// compute the kinetic energy
  Real getKineticEnergy();

  [[deprecated("Use the interface with an Element")]] Real
  getKineticEnergy(ElementType type, Idx index) {
    return getKineticEnergy({type, index, _not_ghost});
  }

  Real getKineticEnergy(const Element & element);

  /// compute the external work (for impose displacement, the velocity should be
  /// given too)
  Real getExternalWork();

public:
  /// get the energies
  Real getEnergy(const std::string & energy_id);

  /// compute the energy for an element
  [[deprecated("Use the interface with an Element")]] Real
  getEnergy(const std::string & energy_id, ElementType type, Idx index) {
    return getEnergy(energy_id, Element{type, index, _not_ghost});
  };

  /// compute the energy for an element
  Real getEnergy(const std::string & energy_id, const Element & element);

  /// compute the energy for an element group
  Real getEnergy(const ID & energy_id, const ID & group_id);

  /* ------------------------------------------------------------------------ */
  /* Data Accessor inherited members                                          */
  /* ------------------------------------------------------------------------ */
public:
  [[nodiscard]] Int getNbData(const Array<Element> & elements,
                              const SynchronizationTag & tag) const override;

  void packData(CommunicationBuffer & buffer, const Array<Element> & elements,
                const SynchronizationTag & tag) const override;

  void unpackData(CommunicationBuffer & buffer, const Array<Element> & elements,
                  const SynchronizationTag & tag) override;

  [[nodiscard]] Int getNbData(const Array<Idx> & dofs,
                              const SynchronizationTag & tag) const override;

  void packData(CommunicationBuffer & buffer, const Array<Idx> & dofs,
                const SynchronizationTag & tag) const override;

  void unpackData(CommunicationBuffer & buffer, const Array<Idx> & dofs,
                  const SynchronizationTag & tag) override;

  /* ------------------------------------------------------------------------ */
  /* Mesh Event Handler inherited members                                     */
  /* ------------------------------------------------------------------------ */
protected:
  void onNodesAdded(const Array<Idx> & nodes_list,
                    const NewNodesEvent & event) override;
  void onNodesRemoved(const Array<Idx> & element_list,
                      const Array<Idx> & new_numbering,
                      const RemovedNodesEvent & event) override;

  /* ------------------------------------------------------------------------ */
  /* Dumpable interface (kept for convenience) and dumper relative functions  */
  /* ------------------------------------------------------------------------ */
public:
  virtual void onDump();

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
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// set the value of the time step
  void setTimeStep(Real time_step, const ID & solver_id = "") override;

  /// get the value of the conversion from forces/ mass to acceleration
  AKANTU_GET_MACRO(F_M2A, f_m2a, Real);

  /// set the value of the conversion from forces/ mass to acceleration
  AKANTU_SET_MACRO(F_M2A, f_m2a, Real);

  /// get the SolidMechanicsModel::displacement array
  AKANTU_GET_MACRO_DEREF_PTR_NOT_CONST(Displacement, displacement);
  /// get the SolidMechanicsModel::displacement array
  AKANTU_GET_MACRO_DEREF_PTR(Displacement, displacement);

  /// get the SolidMechanicsModel::previous_displacement array
  AKANTU_GET_MACRO_DEREF_PTR(PreviousDisplacement, previous_displacement);

  /// get the SolidMechanicsModel::current_position array
  const Array<Real> & getCurrentPosition();

  /// get  the SolidMechanicsModel::displacement_increment  array
  AKANTU_GET_MACRO_DEREF_PTR(Increment, displacement_increment);
  /// get  the SolidMechanicsModel::displacement_increment  array
  AKANTU_GET_MACRO_DEREF_PTR_NOT_CONST(Increment, displacement_increment);

  /// get the lumped SolidMechanicsModel::mass array
  AKANTU_GET_MACRO_DEREF_PTR(Mass, mass);

  /// get the SolidMechanicsModel::velocity array
  AKANTU_GET_MACRO_DEREF_PTR_NOT_CONST(Velocity, velocity);
  /// get the SolidMechanicsModel::velocity array
  AKANTU_GET_MACRO_DEREF_PTR(Velocity, velocity);

  /// get    the    SolidMechanicsModel::acceleration   array
  AKANTU_GET_MACRO_DEREF_PTR_NOT_CONST(Acceleration, acceleration);
  /// get    the    SolidMechanicsModel::acceleration   array
  AKANTU_GET_MACRO_DEREF_PTR(Acceleration, acceleration);

  /// get the SolidMechanicsModel::external_force array
  AKANTU_GET_MACRO_DEREF_PTR_NOT_CONST(ExternalForce, external_force);
  /// get the SolidMechanicsModel::external_force array
  AKANTU_GET_MACRO_DEREF_PTR(ExternalForce, external_force);

  /// get the SolidMechanicsModel::internal_force array (internal forces)
  AKANTU_GET_MACRO_DEREF_PTR_NOT_CONST(InternalForce, internal_force);
  /// get the SolidMechanicsModel::internal_force array (internal forces)
  AKANTU_GET_MACRO_DEREF_PTR(InternalForce, internal_force);

  /// get the SolidMechanicsModel::blocked_dofs array
  AKANTU_GET_MACRO_DEREF_PTR_NOT_CONST(BlockedDOFs, blocked_dofs);
  /// get the SolidMechanicsModel::blocked_dofs array
  AKANTU_GET_MACRO_DEREF_PTR(BlockedDOFs, blocked_dofs);

  AKANTU_GET_MACRO_AUTO(DisplacementRelease, displacement_release);
  AKANTU_GET_MACRO_AUTO(CurrentPositionRelease, current_position_release);

  /// get an iterable on the materials
  inline decltype(auto) getMaterials() { return this->getConstitutiveLaws(); }

  /// get an iterable on the materials
  [[nodiscard]] inline decltype(auto) getMaterials() const {
    return this->getConstitutiveLaws();
  }

  /// get a particular material (by numerical material index)
  inline Material & getMaterial(Idx mat_index) {
    return this->getConstitutiveLaw(mat_index);
  }

  /// get a particular material (by numerical material index)
  [[nodiscard]] inline const Material & getMaterial(Idx mat_index) const {
    return this->getConstitutiveLaw(mat_index);
  }

  /// get a particular material (by material name)
  inline Material & getMaterial(const std::string & name) {
    return this->getConstitutiveLaw(name);
  }

  /// get a particular material (by material name)
  [[nodiscard]] inline const Material &
  getMaterial(const std::string & name) const {
    return this->getConstitutiveLaw(name);
  }

  /// get a particular material (by material name)
  [[nodiscard]] inline const Material &
  getMaterial(const Element & element) const {
    return this->getConstitutiveLaw(element);
  }

  /// get a particular material id from is name
  [[nodiscard]] inline auto getMaterialIndex(const std::string & name) const {
    return this->getConstitutiveLawIndex(name);
  }

  /// give the number of materials
  [[nodiscard]] inline auto getNbMaterials() const {
    return this->getNbConstitutiveLaws();
  }

  void reassignMaterial() { this->reassignConstitutiveLaw(); }
  void registerNewMaterial(const ID & mat_name, const ID & mat_type,
                           const ID & opt_param);

  /// compute the stable time step
  Real getStableTimeStep();

  // this function is kept for backward compatinility
  [[nodiscard]] decltype(auto) getMaterialByElement() const {
    return this->getConstitutiveLawByElement();
  }

  // this function is kept for backward compatinility
  [[nodiscard]] decltype(auto) getMaterialLocalNumbering() const {
    return this->getConstitutiveLawLocalNumbering();
  }

  // this function is kept for backward compatinility
  [[nodiscard]] decltype(auto)
  getMaterialByElement(ElementType type,
                       GhostType ghost_type = _not_ghost) const {
    return this->getConstitutiveLawByElement(type, ghost_type);
  }

  // this function is kept for backward compatinility
  [[nodiscard]] decltype(auto)
  getMaterialLocalNumbering(ElementType type,
                            GhostType ghost_type = _not_ghost) const {
    return this->getConstitutiveLawLocalNumbering(type, ghost_type);
  }

  // this function is kept for backward compatinility
  decltype(auto) getMaterialSelector() {
    return this->getConstitutiveLawSelector();
  }

  // this function is kept for backward compatinility
  void setMaterialSelector(
      const std::shared_ptr<ConstitutiveLawSelector> & material_selector) {
    this->setConstitutiveLawSelector(material_selector);
  }

  /// get the FEEngine object to integrate or interpolate on the boundary
  FEEngine & getFEEngineBoundary(const ID & name = "") override;

protected:
  /// compute the stable time step
  Real getStableTimeStep(GhostType ghost_type);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// release version of the displacement array
  Int displacement_release{0};

  /// release version of the current_position array
  Int current_position_release{0};

  /// Check if materials need to recompute the mass array
  bool need_to_reassemble_lumped_mass{true};

  /// Check if materials need to recompute the mass matrix
  bool need_to_reassemble_mass{true};

protected:
  /// conversion coefficient form force/mass to acceleration
  Real f_m2a{1.0};

  /// displacements array
  std::unique_ptr<Array<Real>> displacement;

  /// displacements array at the previous time step (used in finite deformation)
  std::unique_ptr<Array<Real>> previous_displacement;

  /// increment of displacement
  std::unique_ptr<Array<Real>> displacement_increment;

  /// lumped mass array
  std::unique_ptr<Array<Real>> mass;

  /// velocities array
  std::unique_ptr<Array<Real>> velocity;

  /// accelerations array
  std::unique_ptr<Array<Real>> acceleration;

  /// external forces array
  std::unique_ptr<Array<Real>> external_force;

  /// internal forces array
  std::unique_ptr<Array<Real>> internal_force;

  /// array specifing if a degree of freedom is blocked or not
  std::unique_ptr<Array<bool>> blocked_dofs;

  /// array of current position used during update residual
  std::unique_ptr<Array<Real>> current_position;
};

/* -------------------------------------------------------------------------- */
namespace BC {
  namespace Neumann {
    using FromStress = FromHigherDim;
    using FromTraction = FromSameDim;
  } // namespace Neumann
} // namespace BC

} // namespace akantu

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
#include "material.hh"
#include "material_selector_tmpl.hh"
#include "parser.hh"
/* -------------------------------------------------------------------------- */

#endif /* AKANTU_SOLID_MECHANICS_MODEL_HH_ */
