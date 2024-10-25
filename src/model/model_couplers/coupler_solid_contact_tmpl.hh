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
#include "coupler_solid_contact.hh"
#include "dumpable_inline_impl.hh"
/* -------------------------------------------------------------------------- */
#include "dumper_iohelper_paraview.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType>
void CouplerSolidContactTemplate<SolidMechanicsModelType>::initSolver(
    TimeStepSolverType time_step_solver_type) {
  auto & solid_model_solver = aka::as_type<ModelSolver>(*solid);
  solid_model_solver.initSolver(time_step_solver_type);

  auto & contact_model_solver = aka::as_type<ModelSolver>(*contact);
  contact_model_solver.initSolver(time_step_solver_type);
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType>
std::tuple<ID, TimeStepSolverType>
CouplerSolidContactTemplate<SolidMechanicsModelType>::getDefaultSolverID(
    const AnalysisMethod & method) {
  return solid->getDefaultSolverID(method);
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType>
TimeStepSolverType
CouplerSolidContactTemplate<SolidMechanicsModelType>::getDefaultSolverType()
    const {
  return solid->getDefaultSolverType();
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType>
ModelSolverOptions
CouplerSolidContactTemplate<SolidMechanicsModelType>::getDefaultSolverOptions(
    const TimeStepSolverType & type) const {
  return solid->getDefaultSolverOptions(type);
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType>
void CouplerSolidContactTemplate<SolidMechanicsModelType>::assembleResidual() {
  // computes the internal forces
  switch (method) {
  case _explicit_lumped_mass: {
    auto & current_positions = contact->getContactDetector().getPositions();
    current_positions.copy(solid->getCurrentPosition());
    contact->search();
    break;
  }
  default:
    break;
  }

  this->assembleInternalForces();

  auto & internal_force = solid->getInternalForce();
  auto & external_force = solid->getExternalForce();

  auto & contact_force = contact->getInternalForce();

  /* ------------------------------------------------------------------------ */
  this->getDOFManager().assembleToResidual("displacement", external_force, 1);
  this->getDOFManager().assembleToResidual("displacement", internal_force, 1);
  this->getDOFManager().assembleToResidual("displacement", contact_force, 1);
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType>
void CouplerSolidContactTemplate<SolidMechanicsModelType>::assembleResidual(
    const ID & residual_part) {
  AKANTU_DEBUG_IN();

  // contact->assembleInternalForces();

  auto & internal_force = solid->getInternalForce();
  auto & external_force = solid->getExternalForce();
  auto & contact_force = contact->getInternalForce();

  if ("external" == residual_part) {
    this->getDOFManager().assembleToResidual("displacement", external_force, 1);
    this->getDOFManager().assembleToResidual("displacement", contact_force, 1);
    AKANTU_DEBUG_OUT();
    return;
  }

  if ("internal" == residual_part) {
    this->getDOFManager().assembleToResidual("displacement", internal_force, 1);
    AKANTU_DEBUG_OUT();
    return;
  }

  AKANTU_CUSTOM_EXCEPTION(
      debug::SolverCallbackResidualPartUnknown(residual_part));

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType>
void CouplerSolidContactTemplate<SolidMechanicsModelType>::predictor() {
  auto & solid_model_solver = aka::as_type<ModelSolver>(*solid);
  solid_model_solver.predictor();
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType>
void CouplerSolidContactTemplate<SolidMechanicsModelType>::corrector() {
  auto & solid_model_solver = aka::as_type<ModelSolver>(*solid);
  solid_model_solver.corrector();

  switch (method) {
  case _static:
  case _implicit_dynamic: {
    auto & current_positions = contact->getContactDetector().getPositions();
    current_positions.copy(solid->getCurrentPosition());
    contact->search();
    break;
  }
  default:
    break;
  }
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType>
MatrixType CouplerSolidContactTemplate<SolidMechanicsModelType>::getMatrixType(
    const ID & matrix_id) const {
  if (matrix_id == "K") {
    return _symmetric;
  }
  if (matrix_id == "M") {
    return _symmetric;
  }
  return _mt_not_defined;
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType>
void CouplerSolidContactTemplate<SolidMechanicsModelType>::assembleMatrix(
    const ID & matrix_id) {
  if (matrix_id == "K") {
    this->assembleStiffnessMatrix();
  } else if (matrix_id == "M") {
    solid->assembleMass();
  }
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType>
void CouplerSolidContactTemplate<SolidMechanicsModelType>::assembleLumpedMatrix(
    const ID & matrix_id) {
  if (matrix_id == "M") {
    solid->assembleMassLumped();
  }
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType>
void CouplerSolidContactTemplate<SolidMechanicsModelType>::beforeSolveStep() {
  auto & solid_solver_callback = aka::as_type<SolverCallback>(*solid);
  solid_solver_callback.beforeSolveStep();

  auto & contact_solver_callback = aka::as_type<SolverCallback>(*contact);
  contact_solver_callback.beforeSolveStep();
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType>
void CouplerSolidContactTemplate<SolidMechanicsModelType>::afterSolveStep(
    bool converged) {
  auto & solid_solver_callback = aka::as_type<SolverCallback>(*solid);
  solid_solver_callback.afterSolveStep(converged);

  auto & contact_solver_callback = aka::as_type<SolverCallback>(*contact);
  contact_solver_callback.afterSolveStep(converged);
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType>
void CouplerSolidContactTemplate<
    SolidMechanicsModelType>::assembleInternalForces() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_INFO("Assemble the internal forces");

  solid->assembleInternalForces();
  contact->assembleInternalForces();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType>
void CouplerSolidContactTemplate<
    SolidMechanicsModelType>::assembleStiffnessMatrix() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_INFO("Assemble the new stiffness matrix");

  solid->assembleStiffnessMatrix(true);

  switch (method) {
  case _static:
  case _implicit_dynamic: {
    contact->assembleStiffnessMatrix();
    break;
  }
  default:
    break;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType>
void CouplerSolidContactTemplate<
    SolidMechanicsModelType>::assembleMassLumped() {
  solid->assembleMassLumped();
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType>
void CouplerSolidContactTemplate<SolidMechanicsModelType>::assembleMass() {
  solid->assembleMass();
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType>
void CouplerSolidContactTemplate<SolidMechanicsModelType>::assembleMassLumped(
    GhostType ghost_type) {
  solid->assembleMassLumped(ghost_type);
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType>
void CouplerSolidContactTemplate<SolidMechanicsModelType>::assembleMass(
    GhostType ghost_type) {
  solid->assembleMass(ghost_type);
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType>
std::shared_ptr<dumpers::Field>
CouplerSolidContactTemplate<SolidMechanicsModelType>::createElementalField(
    const std::string & field_name, const std::string & group_name,
    bool padding_flag, Int spatial_dimension, ElementKind kind) {

  std::shared_ptr<dumpers::Field> field;
  field = contact->createElementalField(field_name, group_name, padding_flag,
                                        spatial_dimension, kind);
  if (not field) {
    field = solid->createElementalField(field_name, group_name, padding_flag,
                                        spatial_dimension, kind);
  }
  return field;
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType>
std::shared_ptr<dumpers::Field>
CouplerSolidContactTemplate<SolidMechanicsModelType>::createNodalFieldReal(
    const std::string & field_name, const std::string & group_name,
    bool padding_flag) {

  std::shared_ptr<dumpers::Field> field;
  field = contact->createNodalFieldReal(field_name, group_name, padding_flag);
  if (not field) {
    field = solid->createNodalFieldReal(field_name, group_name, padding_flag);
  }

  return field;
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType>
std::shared_ptr<dumpers::Field>
CouplerSolidContactTemplate<SolidMechanicsModelType>::createNodalFieldInt(
    const std::string & field_name, const std::string & group_name,
    bool padding_flag) {
  std::shared_ptr<dumpers::Field> field;
  field = contact->createNodalFieldInt(field_name, group_name, padding_flag);
  if (not field) {
    field = solid->createNodalFieldInt(field_name, group_name, padding_flag);
  }

  return field;
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType>
std::shared_ptr<dumpers::Field>
CouplerSolidContactTemplate<SolidMechanicsModelType>::createNodalFieldBool(
    const std::string & field_name, const std::string & group_name,
    bool padding_flag) {
  std::shared_ptr<dumpers::Field> field;
  field = contact->createNodalFieldBool(field_name, group_name, padding_flag);
  if (not field) {
    field = solid->createNodalFieldBool(field_name, group_name, padding_flag);
  }
  return field;
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType>
void CouplerSolidContactTemplate<SolidMechanicsModelType>::dump(
    const std::string & dumper_name) {
  solid->onDump();
  Model::dump(dumper_name);
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType>
void CouplerSolidContactTemplate<SolidMechanicsModelType>::dump(
    const std::string & dumper_name, Int step) {
  solid->onDump();
  Model::dump(dumper_name, step);
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType>
void CouplerSolidContactTemplate<SolidMechanicsModelType>::dump(
    const std::string & dumper_name, Real time, Int step) {
  solid->onDump();
  Model::dump(dumper_name, time, step);
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType>
void CouplerSolidContactTemplate<SolidMechanicsModelType>::dump() {
  solid->onDump();
  Model::dump();
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType>
void CouplerSolidContactTemplate<SolidMechanicsModelType>::dump(Int step) {
  solid->onDump();
  Model::dump(step);
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType>
void CouplerSolidContactTemplate<SolidMechanicsModelType>::dump(Real time,
                                                                Int step) {
  solid->onDump();
  Model::dump(time, step);
}

} // namespace akantu
