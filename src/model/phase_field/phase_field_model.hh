/**
 * Copyright (©) 2018-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "model.hh"
#include "phasefield_selector.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_PHASE_FIELD_MODEL_HH_
#define AKANTU_PHASE_FIELD_MODEL_HH_

namespace akantu {
class PhaseField;
template <ElementKind kind, class IntegrationOrderFuntor> class IntegratorGauss;
template <ElementKind kind> class ShapeLagrange;
} // namespace akantu

/* -------------------------------------------------------------------------- */
namespace akantu {

/* -------------------------------------------------------------------------- */
class PhaseFieldModel : public ConstitutiveLawsHandler<PhaseField, Model>,
                        public BoundaryCondition<PhaseFieldModel> {
  using Parent = ConstitutiveLawsHandler<PhaseField, Model>;
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  using FEEngineType = FEEngineTemplate<IntegratorGauss, ShapeLagrange>;

  PhaseFieldModel(Mesh & mesh, Int dim = _all_dimensions,
                  const ID & id = "phase_field_model",
                  std::shared_ptr<DOFManager> dof_manager = nullptr,
                  ModelType model_type = ModelType::_phase_field_model);

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
protected:
  /// generic function to initialize everything ready for explicit dynamics
  void initFullImpl(const ModelOptions & options) override;

  /// allocate all vectors
  void initSolver(TimeStepSolverType /*unused*/) override;

  /// predictor
  void predictor() override;

  /// corrector
  void corrector() override;

  /// compute the heat flux
  void assembleResidual() override;

  /// get the type of matrix needed
  [[nodiscard]] MatrixType getMatrixType(const ID & /*unused*/) const override;

  /// callback to assemble a Matrix
  void assembleMatrix(const ID & /*unused*/) override;

  /// callback to assemble a lumped Matrix
  void assembleLumpedMatrix(const ID & /*unused*/) override;

protected:
  /* ------------------------------------------------------------------------ */
  [[nodiscard]] TimeStepSolverType getDefaultSolverType() const override;

  std::tuple<ID, TimeStepSolverType>
  getDefaultSolverID(const AnalysisMethod & method) override;

  [[nodiscard]] ModelSolverOptions
  getDefaultSolverOptions(const TimeStepSolverType & type) const override;

  /// set the element_id_by_phasefield and add the elements to the good
  /// phasefields
  void
  assignPhaseFieldToElements(const ElementTypeMapArray<Idx> * filter = nullptr);

  /* ------------------------------------------------------------------------ */
  /* Methods for static                                                       */
  /* ------------------------------------------------------------------------ */
public:
  /// assembles the phasefield stiffness matrix
  virtual void assembleStiffnessMatrix();

  /// compute the internal forces
  virtual void assembleInternalForces();

  // compute the internal forces
  void assembleInternalForces(GhostType ghost_type);

  // save fields history
  void savePreviousState();

  /* ------------------------------------------------------------------------ */
  /* Methods for dynamic                                                      */
  /* ------------------------------------------------------------------------ */
public:
  /// set the stable timestep
  void setTimeStep(Real time_step, const ID & solver_id = "") override;

protected:
  /// callback for the solver, this is called at beginning of solve
  void beforeSolveStep() override;
  /// callback for the solver, this is called at end of solve
  void afterSolveStep(bool converged = true) override;

  void computeNonLocalContribution(GhostType /*ghost_type*/) override{};
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

  [[nodiscard]] Int getNbData(const Array<Idx> & indexes,
                              const SynchronizationTag & tag) const override;

  void packData(CommunicationBuffer & buffer, const Array<Idx> & indexes,
                const SynchronizationTag & tag) const override;

  void unpackData(CommunicationBuffer & buffer, const Array<Idx> & indexes,
                  const SynchronizationTag & tag) override;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// return the damage array
  AKANTU_GET_MACRO_DEREF_PTR(Damage, damage);
  AKANTU_GET_MACRO_DEREF_PTR_NOT_CONST(Damage, damage);

  /// return the previous damage array
  AKANTU_GET_MACRO_DEREF_PTR(PreviousDamage, previous_damage);
  AKANTU_GET_MACRO_DEREF_PTR_NOT_CONST(PreviousDamage, previous_damage);

  /// get the PhaseFieldModel::internal_force vector (internal forces)
  AKANTU_GET_MACRO_DEREF_PTR(InternalForce, internal_force);
  AKANTU_GET_MACRO_DEREF_PTR_NOT_CONST(InternalForce, internal_force);

  /// get the PhaseFieldModel::external_force vector (external forces)
  AKANTU_GET_MACRO_DEREF_PTR(ExternalForce, external_force);
  AKANTU_GET_MACRO_DEREF_PTR_NOT_CONST(ExternalForce, external_force);

  /// get the PhaseFieldModel::blocked_dofs vector
  AKANTU_GET_MACRO_DEREF_PTR(BlockedDOFs, blocked_dofs);
  /// get the PhaseFieldModel::blocked_dofs vector
  AKANTU_GET_MACRO_DEREF_PTR_NOT_CONST(BlockedDOFs, blocked_dofs);

  /// Returns the total dissipated energy
  [[nodiscard]] virtual Real getEnergy(const ID & energy_id = "dissipated");

  /// Compute dissipated energy for an individual element
  [[nodiscard]] virtual Real getEnergy(const ID & energy_id,
                                       const Element & element);

  /// Compute dissipated energy for an element group
  [[nodiscard]] virtual Real getEnergy(const ID & energy_id,
                                       const ID & group_id);

  FEEngine & getFEEngineBoundary(const ID & name = "") override;

  void
  setPhaseFieldSelector(const std::shared_ptr<PhaseFieldSelector> & selector) {
    this->setConstitutiveLawSelector(selector);
  }

  PhaseField & getPhaseField(const ID & id) { return getConstitutiveLaw(id); }
  PhaseField & getPhaseField(Idx id) { return getConstitutiveLaw(id); }

  [[nodiscard]] const PhaseField & getPhaseField(const ID & id) const {
    return getConstitutiveLaw(id);
  }
  [[nodiscard]] const PhaseField & getPhaseField(Idx id) const {
    return getConstitutiveLaw(id);
  }

  /// give the number of phasefield materials
  [[nodiscard]] inline auto getNbPhaseFields() const {
    return this->getNbConstitutiveLaws();
  }

  /* ------------------------------------------------------------------------ */
  /* Dumpable Interface                                                       */
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

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// number of iterations
  Int n_iter;

  /// damage array
  std::unique_ptr<Array<Real>> damage;

  /// damage array at the previous time step
  std::unique_ptr<Array<Real>> previous_damage;

  /// boundary vector
  std::unique_ptr<Array<bool>> blocked_dofs;

  /// external force vector
  std::unique_ptr<Array<Real>> external_force;

  /// residuals array
  std::unique_ptr<Array<Real>> internal_force;
};

} // namespace akantu

#include "phasefield.hh"
/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
#include "phase_field_model_inline_impl.hh"
/* -------------------------------------------------------------------------- */

#endif
