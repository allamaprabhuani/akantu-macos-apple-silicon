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
#include "contact_detector.hh"
#include "data_accessor.hh"
#include "fe_engine.hh"
#include "model.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_CONTACT_MECHANICS_MODEL_HH_
#define AKANTU_CONTACT_MECHANICS_MODEL_HH_

namespace akantu {
class Resolution;
template <ElementKind kind, class IntegrationOrderFunctor>
class IntegratorGauss;
template <ElementKind kind> class ShapeLagrange;
} // namespace akantu

/* -------------------------------------------------------------------------- */
namespace akantu {

/* -------------------------------------------------------------------------- */
class ContactMechanicsModel : public Model,
                              public DataAccessor<Element>,
                              public DataAccessor<Idx>,
                              public BoundaryCondition<ContactMechanicsModel> {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

  using MyFEEngineType = FEEngineTemplate<IntegratorGauss, ShapeLagrange>;

public:
  ContactMechanicsModel(
      Mesh & mesh, Int dim = _all_dimensions,
      const ID & id = "contact_mechanics_model",
      std::shared_ptr<DOFManager> dof_manager = nullptr,
      ModelType model_type = ModelType::_contact_mechanics_model);

  ~ContactMechanicsModel() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
protected:
  /// initialize completely the model
  void initFullImpl(const ModelOptions & options) override;

  /// allocate all vectors
  void initSolver(TimeStepSolverType /*unused*/) override;

  /// initialize all internal arrays for resolutions
  void initResolutions();

  /// initialize the modelType
  void initModel() override;

  /// call back for the solver, computes the force residual
  void assembleResidual() override;

  /// get the type of matrix needed
  MatrixType getMatrixType(const ID & matrix_id) const override;

  /// callback for the solver, this assembles different matrices
  void assembleMatrix(const ID & matrix_id) override;

  /// callback for the solver, this assembles the stiffness matrix
  void assembleLumpedMatrix(const ID & matrix_id) override;

  /// get some default values for derived classes
  std::tuple<ID, TimeStepSolverType>
  getDefaultSolverID(const AnalysisMethod & method) override;

  [[nodiscard]] ModelSolverOptions
  getDefaultSolverOptions(const TimeStepSolverType & type) const override;

  /// callback for the solver, this is called at beginning of solve
  void beforeSolveStep() override;

  /// callback for the solver, this is called at end of solve
  void afterSolveStep(bool converged = true) override;

  /// function to print the containt of the class
  void printself(std::ostream & stream, int indent = 0) const override;

  /* ------------------------------------------------------------------------ */
  /* Contact Detection                                                        */
  /* ------------------------------------------------------------------------ */
public:
  void search();

  void computeNodalAreas(GhostType ghost_type = _not_ghost);

  /* ------------------------------------------------------------------------ */
  /* Contact Resolution                                                       */
  /* ------------------------------------------------------------------------ */
public:
  /// register an empty contact resolution of a given type
  Resolution & registerNewResolution(const ID & res_name, const ID & res_type,
                                     const ID & opt_param);

protected:
  /// register a resolution in the dynamic database
  Resolution & registerNewResolution(const ParserSection & res_section);

  /// read the resolution files to instantiate all the resolutions
  void instantiateResolutions();

  /// save the parameters from previous state
  void savePreviousState();

  /* ------------------------------------------------------------------------ */
  /* Solver Interface                                                         */
  /* ------------------------------------------------------------------------ */
public:
  /// assembles the contact stiffness matrix
  virtual void assembleStiffnessMatrix();

  /// assembles the contant internal forces
  virtual void assembleInternalForces();

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  FEEngine & getFEEngineBoundary(const ID & name = "") override;

  /* ------------------------------------------------------------------------ */
  /* Dumpable interface                                                       */
  /* ------------------------------------------------------------------------ */
public:
  std::shared_ptr<dumpers::Field>
  createNodalFieldReal(const std::string & field_name,
                       const std::string & group_name,
                       bool padding_flag) override;

  std::shared_ptr<dumpers::Field>
  createNodalFieldInt(const std::string & field_name,
                      const std::string & group_name,
                      bool padding_flag) override;

  std::shared_ptr<dumpers::Field>
  createNodalFieldBool(const std::string & field_name,
                       const std::string & group_name,
                       bool padding_flag) override;

  /* ------------------------------------------------------------------------ */
  /* Data Accessor inherited members                                          */
  /* ------------------------------------------------------------------------ */
public:
  Int getNbData(const Array<Element> & elements,
                const SynchronizationTag & tag) const override;

  void packData(CommunicationBuffer & buffer, const Array<Element> & elements,
                const SynchronizationTag & tag) const override;

  void unpackData(CommunicationBuffer & buffer, const Array<Element> & elements,
                  const SynchronizationTag & tag) override;

  Int getNbData(const Array<Idx> & dofs,
                const SynchronizationTag & tag) const override;

  void packData(CommunicationBuffer & buffer, const Array<Idx> & dofs,
                const SynchronizationTag & tag) const override;

  void unpackData(CommunicationBuffer & buffer, const Array<Idx> & dofs,
                  const SynchronizationTag & tag) override;

protected:
  /// contact detection class
  friend class ContactDetector;

  /// contact resolution class
  friend class Resolution;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get the ContactMechanicsModel::displacement vector
  AKANTU_GET_MACRO(Displacement, *displacement, Array<Real> &);

  /// get  the ContactMechanicsModel::increment  vector \warn  only  consistent
  /// if ContactMechanicsModel::setIncrementFlagOn has been called before
  AKANTU_GET_MACRO(Increment, *displacement_increment, Array<Real> &);

  /// get the ContactMechanics::internal_force vector (internal forces)
  AKANTU_GET_MACRO(InternalForce, *internal_force, Array<Real> &);

  /// get the ContactMechanicsModel::external_force vector (external forces)
  AKANTU_GET_MACRO(ExternalForce, *external_force, Array<Real> &);

  /// get the ContactMechanics::normal_force vector (normal forces)
  AKANTU_GET_MACRO(NormalForce, *normal_force, Array<Real> &);

  /// get the ContactMechanics::tangential_force vector (friction forces)
  AKANTU_GET_MACRO(TangentialForce, *tangential_force, Array<Real> &);

  /// get the ContactMechanics::traction vector (friction traction)
  AKANTU_GET_MACRO(TangentialTractions, *tangential_tractions, Array<Real> &);

  /// get the ContactMechanics::previous_tangential_tractions vector
  AKANTU_GET_MACRO(PreviousTangentialTractions, *previous_tangential_tractions,
                   Array<Real> &);

  /// get the ContactMechanicsModel::force vector (external forces)
  Array<Real> & getForce() {
    AKANTU_DEBUG_WARNING("getForce was maintained for backward compatibility, "
                         "use getExternalForce instead");
    return *external_force;
  }

  /// get the ContactMechanics::blocked_dofs vector
  AKANTU_GET_MACRO(BlockedDOFs, *blocked_dofs, Array<bool> &);

  /// get the ContactMechanics::gaps (contact gaps)
  AKANTU_GET_MACRO(Gaps, *gaps, Array<Real> &);

  /// get the ContactMechanics::normals (normals on slave nodes)
  AKANTU_GET_MACRO(Normals, *normals, Array<Real> &);

  /// get the ContactMechanics::tangents (tangents on slave nodes)
  AKANTU_GET_MACRO(Tangents, *tangents, Array<Real> &);

  /// get the ContactMechanics::previous_tangents (tangents on slave nodes)
  AKANTU_GET_MACRO(PreviousTangents, *previous_tangents, Array<Real> &);

  /// get the ContactMechanics::areas (nodal areas)
  AKANTU_GET_MACRO(NodalArea, *nodal_area, Array<Real> &);

  /// get the ContactMechanics::previous_projections (previous_projections)
  AKANTU_GET_MACRO(PreviousProjections, *previous_projections, Array<Real> &);

  /// get the ContactMechanics::projections (projections)
  AKANTU_GET_MACRO(Projections, *projections, Array<Real> &);

  /// get the ContactMechanics::contact_state vector (no_contact/stick/slip
  /// state)
  AKANTU_GET_MACRO(ContactState, *contact_state, Array<ContactState> &);

  /// get the ContactMechanics::previous_master_elements
  AKANTU_GET_MACRO(PreviousMasterElements, *previous_master_elements,
                   Array<Element> &);

  /// get contact detector
  AKANTU_GET_MACRO_NOT_CONST(ContactDetector, *detector, ContactDetector &);

  /// get the contact elements
  inline const Array<ContactElement> & getContactElements() const {
    return contact_elements;
  }

  /// get the current positions of the nodes
  inline Array<Real> & getPositions() { return detector->getPositions(); }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// tells if the resolutions are instantiated
  bool are_resolutions_instantiated;

  /// displacements array
  std::unique_ptr<Array<Real>> displacement;

  /// increment of displacement
  std::unique_ptr<Array<Real>> displacement_increment;

  /// contact forces array
  std::unique_ptr<Array<Real>> internal_force;

  /// external forces array
  std::unique_ptr<Array<Real>> external_force;

  /// normal force array
  std::unique_ptr<Array<Real>> normal_force;

  /// friction force array
  std::unique_ptr<Array<Real>> tangential_force;

  /// friction traction array
  std::unique_ptr<Array<Real>> tangential_tractions;

  /// previous friction traction array
  std::unique_ptr<Array<Real>> previous_tangential_tractions;

  /// boundary vector
  std::unique_ptr<Array<bool>> blocked_dofs;

  /// array to store gap between slave and master
  std::unique_ptr<Array<Real>> gaps;

  /// array to store normals from master to slave
  std::unique_ptr<Array<Real>> normals;

  /// array to store tangents on the master element
  std::unique_ptr<Array<Real>> tangents;

  /// array to store previous tangents on the master element
  std::unique_ptr<Array<Real>> previous_tangents;

  /// array to store nodal areas
  std::unique_ptr<Array<Real>> nodal_area;

  /// array to store stick/slip state :
  std::unique_ptr<Array<ContactState>> contact_state;

  /// array to store previous projections in covariant basis
  std::unique_ptr<Array<Real>> previous_projections;

  // array to store projections in covariant basis
  std::unique_ptr<Array<Real>> projections;

  /// contact detection
  std::unique_ptr<ContactDetector> detector;

  /// list of contact resolutions
  std::vector<std::unique_ptr<Resolution>> resolutions;

  /// mapping between resolution name and resolution internal id
  std::map<std::string, UInt> resolutions_names_to_id;

  /// array to store contact elements
  Array<ContactElement> contact_elements;

  /// array to store previous master elements
  std::unique_ptr<Array<Element>> previous_master_elements;
};

} // namespace akantu

/* ------------------------------------------------------------------------ */
/* inline functions                                                         */
/* ------------------------------------------------------------------------ */
#include "parser.hh"
#include "resolution.hh"
/* ------------------------------------------------------------------------ */

#endif /* AKANTU_CONTACT_MECHANICS_MODEL_HH_ */
