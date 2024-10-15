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
#include "constitutive_laws_handler.hh"
#include "data_accessor.hh"
#include "fe_engine.hh"
#include "model.hh"
/* -------------------------------------------------------------------------- */
#include <array>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_DIFFUSION_MODEL_HH_
#define AKANTU_DIFFUSION_MODEL_HH_

namespace akantu {
template <ElementKind kind, class IntegrationOrderFunctor>
class IntegratorGauss;
template <ElementKind kind> class ShapeLagrange;
class DiffusionLaw;
} // namespace akantu

namespace akantu {

class DiffusionModel : public ConstitutiveLawsHandler<DiffusionLaw, Model> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
  using Parent = ConstitutiveLawsHandler<DiffusionLaw, Model>;

public:
  using FEEngineType = FEEngineTemplate<IntegratorGauss, ShapeLagrange>;

  DiffusionModel(Mesh & mesh, Int dim = _all_dimensions,
                 const ID & id = "diffusion",
                 const std::shared_ptr<DOFManager> & dof_manager = nullptr,
                 const ID & dof_name = "diffusion",
                 ModelType model_type = ModelType::_heat_transfer_model);

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
protected:
  /// allocate all vectors
  void initSolver(TimeStepSolverType time_step_solver_type) override;

  void predictor() override;
  void corrector() override;

  /// compute the heat flux
  void assembleResidual() override;

  /// callback to assemble a Matrix
  void assembleMatrix(const ID & matrix_id) override;

  /// callback to assemble a lumped Matrix
  void assembleLumpedMatrix(const ID & matrix_id) override;

  /// get the type of matrix needed
  [[nodiscard]] MatrixType getMatrixType(const ID & matrix_id) const override;

  [[nodiscard]] std::tuple<ID, TimeStepSolverType>
  getDefaultSolverID(const AnalysisMethod & method) override;

  [[nodiscard]] ModelSolverOptions
  getDefaultSolverOptions(const TimeStepSolverType & type) const override;
  /* ------------------------------------------------------------------------ */
  /* Methods for explicit                                                     */
  /* ------------------------------------------------------------------------ */
public:
  /// compute and get the stable time step
  Real getStableTimeStep();

  /// set the stable time step
  void setTimeStep(Real time_step, const ID & solver_id = "") override;

public:
  /// assemble the conductivity matrix
  void assembleDiffisivityMatrix();

  /// assemble the conductivity matrix
  void assembleCapacityMatrix();

  /// calculate the lumped capacity vector for heat transfer problem
  void assembleCapacityMatrixLumped();

private:
  /// compute the internal flow
  void assembleInternalFlow();

  /// calculate the lumped capacity vector for heat transfer problem (w
  /// ghost type)
  void assembleCapacityLumped(GhostType ghost_type);

  /// compute the thermal energy
  Real computeThermalEnergyByNode();

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

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get the boundary vector
  AKANTU_GET_MACRO_DEREF_PTR_NOT_CONST(BlockedDOFs, blocked_dofs);
  /// get the external heat rate vector
  AKANTU_GET_MACRO_DEREF_PTR(ExternalFlow, external_flow);
  /// get the external heat rate vector
  AKANTU_GET_MACRO_DEREF_PTR_NOT_CONST(ExternalFlow, external_flow);
  /// get the assembled heat flux
  AKANTU_GET_MACRO_DEREF_PTR(InternalFlow, internal_flow);
  /// get the assembled heat flux
  AKANTU_GET_MACRO_DEREF_PTR_NOT_CONST(InternalFlow, internal_flow);
  /// get the temperature
  AKANTU_GET_MACRO_DEREF_PTR(Diffusion, diffusion);
  /// get the temperature
  AKANTU_GET_MACRO_DEREF_PTR_NOT_CONST(Diffusion, diffusion);
  /// get the temperature derivative
  AKANTU_GET_MACRO_DEREF_PTR(DiffusionRate, diffusion_rate);

  AKANTU_GET_MACRO_AUTO(DOFName, dof_name);

  AKANTU_GET_MACRO_AUTO(DiffusionRelease, diffusion_release);

  /// get the energy denominated by thermal
  Real getEnergy(const ID & energy_id, const Element & element);
  /// get the energy denominated by thermal
  Real getEnergy(const ID & energy_id);

protected:
  /* ------------------------------------------------------------------------ */
  FEEngine & getFEEngineBoundary(const ID & name = "") override;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  ID dof_name;

  /// temperatures array
  std::unique_ptr<Array<Real>> diffusion;

  /// temperatures derivatives array
  std::unique_ptr<Array<Real>> diffusion_rate;

  /// increment array (@f$\delta \dot T@f$ or @f$\delta T@f$)
  std::unique_ptr<Array<Real>> increment;

  /// external flux vector
  std::unique_ptr<Array<Real>> external_flow;

  /// residuals array
  std::unique_ptr<Array<Real>> internal_flow;

  /// boundary vector
  std::unique_ptr<Array<bool>> blocked_dofs;

  bool need_to_reassemble_capacity{true};
  bool need_to_reassemble_capacity_lumped{true};
  Int diffusion_release{-1};
};

} // namespace akantu

#include "diffusion_law.hh"

#endif /* AKANTU_DIFFUSION_MODEL_HH_ */
