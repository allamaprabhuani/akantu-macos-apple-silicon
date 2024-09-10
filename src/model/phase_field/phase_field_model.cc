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
#include "phase_field_model.hh"
#include "dumpable_inline_impl.hh"
#include "element_synchronizer.hh"
#include "fe_engine_template.hh"
#include "group_manager_inline_impl.hh"
#include "integrator_gauss.hh"
#include "shape_lagrange.hh"
/* -------------------------------------------------------------------------- */
#include "dumper_element_partition.hh"
#include "dumper_elemental_field.hh"
#include "dumper_internal_material_field.hh"
#include "dumper_iohelper_paraview.hh"
#include <algorithm>
#include <utility>
/* -------------------------------------------------------------------------- */
namespace akantu {

/* -------------------------------------------------------------------------- */
PhaseFieldModel::PhaseFieldModel(Mesh & mesh, Int dim, const ID & id,
                                 std::shared_ptr<DOFManager> dof_manager,
                                 ModelType model_type)
    : Parent(mesh, model_type, dim, id) {
  AKANTU_DEBUG_IN();

  this->initDOFManager(std::move(dof_manager));

  this->registerFEEngineObject<FEEngineType>("PhaseFieldFEEngine", mesh,
                                             Model::spatial_dimension);

  this->mesh.registerDumper<DumperParaview>("phase_field", id, true);
  this->mesh.addDumpMesh(mesh, Model::spatial_dimension, _not_ghost,
                         _ek_regular);

  if (this->mesh.isDistributed()) {
    auto & synchronizer = this->mesh.getElementSynchronizer();
    this->registerSynchronizer(synchronizer, SynchronizationTag::_pfm_damage);
    this->registerSynchronizer(synchronizer, SynchronizationTag::_for_dump);
  }

  this->parser_type = ParserType::_phasefield;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
MatrixType PhaseFieldModel::getMatrixType(const ID & matrix_id) const {
  if (matrix_id == "K") {
    return _symmetric;
  }
  return _mt_not_defined;
}

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::initFullImpl(const ModelOptions & options) {
  Parent::initFullImpl(options);

  this->initBC(*this, *damage, *external_force);
}

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::assembleMatrix(const ID & matrix_id) {

  if (matrix_id == "K") {
    this->assembleStiffnessMatrix();
  } else {
    AKANTU_ERROR("Unknown Matrix ID for PhaseFieldModel : " << matrix_id);
  }
}

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::predictor() {
  // AKANTU_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::corrector() {
  // AKANTU_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::initSolver(TimeStepSolverType time_step_solver_type,
                                 NonLinearSolverType /*unused*/) {
  DOFManager & dof_manager = this->getDOFManager();

  this->allocNodalField(this->damage, 1, "damage");
  this->allocNodalField(this->external_force, 1, "external_force");
  this->allocNodalField(this->internal_force, 1, "internal_force");
  this->allocNodalField(this->blocked_dofs, 1, "blocked_dofs");
  this->allocNodalField(this->previous_damage, 1, "previous_damage");

  if (!dof_manager.hasDOFs("damage")) {
    dof_manager.registerDOFs("damage", *this->damage, _dst_nodal);
    dof_manager.registerBlockedDOFs("damage", *this->blocked_dofs);
    dof_manager.registerDOFsPrevious("damage", *this->previous_damage);
  }

  if (time_step_solver_type == TimeStepSolverType::_dynamic) {
    AKANTU_TO_IMPLEMENT();
  }
}

/* -------------------------------------------------------------------------- */
FEEngine & PhaseFieldModel::getFEEngineBoundary(const ID & name) {
  return dynamic_cast<FEEngine &>(getFEEngineClassBoundary<FEEngineType>(name));
}

/* -------------------------------------------------------------------------- */
TimeStepSolverType PhaseFieldModel::getDefaultSolverType() const {
  return TimeStepSolverType::_static;
}

/* -------------------------------------------------------------------------- */
std::tuple<ID, TimeStepSolverType>
PhaseFieldModel::getDefaultSolverID(const AnalysisMethod & method) {
  switch (method) {
  case _explicit_lumped_mass: {
    return std::make_tuple("explicit_lumped",
                           TimeStepSolverType::_dynamic_lumped);
  }
  case _explicit_consistent_mass: {
    return std::make_tuple("explicit", TimeStepSolverType::_dynamic);
  }
  case _static: {
    return std::make_tuple("static", TimeStepSolverType::_static);
  }
  case _implicit_dynamic: {
    return std::make_tuple("implicit", TimeStepSolverType::_dynamic);
  }
  default:
    return std::make_tuple("unknown", TimeStepSolverType::_not_defined);
  }
}

/* -------------------------------------------------------------------------- */
ModelSolverOptions PhaseFieldModel::getDefaultSolverOptions(
    const TimeStepSolverType & type) const {
  ModelSolverOptions options;

  switch (type) {
  case TimeStepSolverType::_dynamic_lumped: {
    options.non_linear_solver_type = NonLinearSolverType::_lumped;
    options.integration_scheme_type["damage"] =
        IntegrationSchemeType::_central_difference;
    options.solution_type["damage"] = IntegrationScheme::_acceleration;
    break;
  }
  case TimeStepSolverType::_static: {
    options.non_linear_solver_type = NonLinearSolverType::_newton_raphson;
    options.integration_scheme_type["damage"] =
        IntegrationSchemeType::_pseudo_time;
    options.solution_type["damage"] = IntegrationScheme::_not_defined;
    break;
  }
  case TimeStepSolverType::_dynamic: {
    options.non_linear_solver_type = NonLinearSolverType::_newton_raphson;
    options.integration_scheme_type["damage"] =
        IntegrationSchemeType::_backward_euler;
    options.solution_type["damage"] = IntegrationScheme::_damage;
    break;
  }
  default:
    AKANTU_EXCEPTION(type << " is not a valid time step solver type");
  }

  return options;
}

/* -------------------------------------------------------------------------- */
Real PhaseFieldModel::getEnergy(const ID & energy_id) {
  AKANTU_DEBUG_IN();

  Real energy = 0.;
  for_each_constitutive_law([&energy, &energy_id](auto && phase_field) {
    energy += phase_field.getEnergy(energy_id);
  });

  /// reduction sum over all processors
  mesh.getCommunicator().allReduce(energy, SynchronizerOperation::_sum);

  AKANTU_DEBUG_OUT();
  return energy;
}

/* -------------------------------------------------------------------------- */
Real PhaseFieldModel::getEnergy(const ID & energy_id, const Element & element) {
  auto pf_element = element;
  auto phase_index = this->getConstitutiveLawByElement()(element);
  pf_element.element = this->getConstitutiveLawLocalNumbering()(element);
  Real energy =
      this->getConstitutiveLaw(phase_index).getEnergy(energy_id, pf_element);
  return energy;
}

/* -------------------------------------------------------------------------- */
Real PhaseFieldModel::getEnergy(const ID & energy_id, const ID & group_id) {
  auto && group = mesh.getElementGroup(group_id);
  auto energy = 0.;
  for (auto && type : group.elementTypes()) {
    for (auto el : group.getElementsIterable(type)) {
      energy += this->getEnergy(energy_id, el);
    }
  }

  /// reduction sum over all processors
  mesh.getCommunicator().allReduce(energy, SynchronizerOperation::_sum);

  return energy;
}

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::beforeSolveStep() {
  for_each_constitutive_law(
      [](auto && phasefield) { phasefield.beforeSolveStep(); });
}

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::afterSolveStep(bool converged) {
  if (not converged) {
    return;
  }

  // for (auto && values : zip(*damage, *previous_damage)) {
  //   auto & dam = std::get<0>(values);
  //   auto & prev_dam = std::get<1>(values);

  //   prev_dam = dam;
  // }

  for (auto && values : zip(*damage, *blocked_dofs)) {
    auto & dam = std::get<0>(values);
    auto & blocked = std::get<1>(values);

    dam = std::min(1., dam);
    if (!blocked && dam > 0.98) {
      blocked = Math::are_float_equal(dam, 1.);
    }
  }

  for_each_constitutive_law(
      [](auto && phasefield) { phasefield.afterSolveStep(); });
}

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::assembleStiffnessMatrix() {
  AKANTU_DEBUG_INFO("Assemble the new stiffness matrix");

  if (!this->getDOFManager().hasMatrix("K")) {
    this->getDOFManager().getNewMatrix("K", getMatrixType("K"));
  }

  this->getDOFManager().zeroMatrix("K");

  for_each_constitutive_law([](auto && phasefield) {
    phasefield.assembleStiffnessMatrix(_not_ghost);
  });
}

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::assembleResidual() {
  this->assembleInternalForces();

  this->getDOFManager().assembleToResidual("damage", *this->external_force, 1);
  this->getDOFManager().assembleToResidual("damage", *this->internal_force, 1);
}

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::assembleInternalForces() {
  AKANTU_DEBUG_INFO("Assemble the internal forces");

  this->internal_force->zero();

  for_each_constitutive_law([](auto && phasefield) {
    phasefield.computeAllDrivingForces(_not_ghost);
  });

  this->asynchronousSynchronize(SynchronizationTag::_pfm_damage);

  // assemble the forces due to local driving forces
  AKANTU_DEBUG_INFO("Assemble residual for local elements");
  for_each_constitutive_law([](auto && phasefield) {
    phasefield.assembleInternalForces(_not_ghost);
  });

  this->waitEndSynchronize(SynchronizationTag::_pfm_damage);
  // assemble the forces due to local driving forces
  AKANTU_DEBUG_INFO("Assemble residual for ghost elements");
  for_each_constitutive_law(
      [](auto && phasefield) { phasefield.assembleInternalForces(_ghost); });
}

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::savePreviousState() {
  for_each_constitutive_law(
      [](auto && phasefield) { phasefield.savePreviousState(); });
}

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::assembleLumpedMatrix(const ID & /*matrix_id*/) {}

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::setTimeStep(Real time_step, const ID & solver_id) {
  Model::setTimeStep(time_step, solver_id);

  this->mesh.getDumper("phase_field").setTimeStep(time_step);
}

/* -------------------------------------------------------------------------- */
Int PhaseFieldModel::getNbData(const Array<Element> & elements,
                               const SynchronizationTag & tag) const {
  Int size = 0;
  Int nb_nodes_per_element = 0;

  for (const Element & el : elements) {
    nb_nodes_per_element += Mesh::getNbNodesPerElement(el.type);
  }

  switch (tag) {
  case SynchronizationTag::_for_dump: {
    // damage
    size += nb_nodes_per_element * sizeof(Real);
    break;
  }
  // case SynchronizationTag::_pfm_damage: {
  //   size += nb_nodes_per_element * sizeof(Real);
  //   break;
  // }
  default: {
    // AKANTU_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }
  size += Parent::getNbData(elements, tag);

  return size;
}

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::packData(CommunicationBuffer & buffer,
                               const Array<Element> & elements,
                               const SynchronizationTag & tag) const {
  switch (tag) {
  case SynchronizationTag::_for_dump: {
    packNodalDataHelper(*damage, buffer, elements, mesh);
    break;
  }
  // case SynchronizationTag::_pfm_damage: {
  //   packNodalDataHelper(*damage, buffer, elements, mesh);
  //   break;
  // }
  default: {
    // AKANTU_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }

  Parent::packData(buffer, elements, tag);
}

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::unpackData(CommunicationBuffer & buffer,
                                 const Array<Element> & elements,
                                 const SynchronizationTag & tag) {
  switch (tag) {
  case SynchronizationTag::_for_dump: {
    unpackNodalDataHelper(*damage, buffer, elements, mesh);
    break;
  }
  // case SynchronizationTag::_pfm_damage: {
  //   unpackNodalDataHelper(*damage, buffer, elements, mesh);
  //   break;
  // }
  default: {
    // AKANTU_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }

  Parent::unpackData(buffer, elements, tag);
}

/* -------------------------------------------------------------------------- */
Int PhaseFieldModel::getNbData(const Array<Idx> & indexes,
                               const SynchronizationTag & tag) const {
  Int size = 0;
  Int nb_nodes = indexes.size();

  switch (tag) {
  case SynchronizationTag::_for_dump: {
    size += nb_nodes * sizeof(Real);
    break;
  }
  case SynchronizationTag::_pfm_damage: {
    size += nb_nodes * sizeof(Real);
    break;
  }
  default: {
    AKANTU_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }
  return size;
}

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::packData(CommunicationBuffer & buffer,
                               const Array<Idx> & indexes,
                               const SynchronizationTag & tag) const {
  switch (tag) {
  case SynchronizationTag::_for_dump: {
    packDOFDataHelper(*damage, buffer, indexes);
    break;
  }
  case SynchronizationTag::_pfm_damage: {
    packDOFDataHelper(*damage, buffer, indexes);
    break;
  }
  default: {
    AKANTU_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }
}

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::unpackData(CommunicationBuffer & buffer,
                                 const Array<Idx> & indexes,
                                 const SynchronizationTag & tag) {
  switch (tag) {
  case SynchronizationTag::_for_dump: {
    unpackDOFDataHelper(*damage, buffer, indexes);
    break;
  }
  case SynchronizationTag::_pfm_damage: {
    unpackDOFDataHelper(*damage, buffer, indexes);
    break;
  }
  default: {
    AKANTU_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }
}

/* -------------------------------------------------------------------------- */
std::shared_ptr<dumpers::Field>
PhaseFieldModel::createNodalFieldBool(const std::string & field_name,
                                      const std::string & group_name,
                                      bool /*unused*/) {
  std::map<std::string, Array<bool> *> uint_nodal_fields;
  uint_nodal_fields["blocked_dofs"] = blocked_dofs.get();

  return mesh.createNodalField(uint_nodal_fields[field_name], group_name);

  std::shared_ptr<dumpers::Field> field;
  return field;
}

/* -------------------------------------------------------------------------- */
std::shared_ptr<dumpers::Field>
PhaseFieldModel::createNodalFieldReal(const std::string & field_name,
                                      const std::string & group_name,
                                      bool /*unused*/) {
  std::map<std::string, Array<Real> *> real_nodal_fields;
  real_nodal_fields["damage"] = damage.get();
  real_nodal_fields["external_force"] = external_force.get();
  real_nodal_fields["internal_force"] = internal_force.get();

  return mesh.createNodalField(real_nodal_fields[field_name], group_name);

  std::shared_ptr<dumpers::Field> field;
  return field;
}

/* -------------------------------------------------------------------------- */
std::shared_ptr<dumpers::Field> PhaseFieldModel::createElementalField(
    const std::string & field_name, const std::string & group_name,
    bool /*unused*/, Int /*unused*/, ElementKind element_kind) {

  if (field_name == "partitions") {
    return mesh.createElementalField<Int, dumpers::ElementPartitionField>(
        mesh.getConnectivities(), group_name, this->spatial_dimension,
        element_kind);
  }

  std::shared_ptr<dumpers::Field> field;
  return field;
}

} // namespace akantu
