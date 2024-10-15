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
#include "solid_mechanics_model.hh"
#include "integrator_gauss.hh"
#include "shape_lagrange.hh"

#include "communicator.hh"
#include "element_synchronizer.hh"
#include "sparse_matrix.hh"
#include "synchronizer_registry.hh"

#include "dumpable_inline_impl.hh" // NOLINT(unused-includes)
/* -------------------------------------------------------------------------- */
#include "dumper_iohelper_paraview.hh"
/* -------------------------------------------------------------------------- */
#include "material_non_local.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
/**
 * A solid mechanics model need a mesh  and a dimension to be created. the model
 * by it  self can not  do a lot,  the good init  functions should be  called in
 * order to configure the model depending on what we want to do.
 *
 * @param  mesh mesh  representing  the model  we  want to  simulate
 * @param dim spatial  dimension of the problem, if dim =  0 (default value) the
 * dimension of the problem is assumed to be the on of the mesh
 * @param id an id to identify the model
 * @param model_type this is an internal parameter for inheritance purposes
 */
SolidMechanicsModel::SolidMechanicsModel(
    Mesh & mesh, Int dim, const ID & id,
    const std::shared_ptr<DOFManager> & dof_manager, const ModelType model_type)
    : CLHParent(mesh, model_type, dim, id) {
  AKANTU_DEBUG_IN();

  this->initDOFManager(dof_manager);

  this->registerFEEngineObject<MyFEEngineType>("SolidMechanicsFEEngine", mesh,
                                               Model::spatial_dimension);

  this->mesh.registerDumper<DumperParaview>("solid_mechanics_model", id, true);
  this->mesh.addDumpMesh(mesh, Model::spatial_dimension, _not_ghost,
                         _ek_regular);

  if (this->mesh.isDistributed()) {
    auto & synchronizer = this->mesh.getElementSynchronizer();
    this->registerSynchronizer(synchronizer, SynchronizationTag::_smm_mass);
    this->registerSynchronizer(synchronizer, SynchronizationTag::_smm_stress);
    this->registerSynchronizer(synchronizer, SynchronizationTag::_smm_gradu);
    this->registerSynchronizer(synchronizer, SynchronizationTag::_for_dump);
  }

  this->parser_type = ParserType::_material;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::setTimeStep(Real time_step, const ID & solver_id) {
  Model::setTimeStep(time_step, solver_id);

  this->mesh.getDumper().setTimeStep(time_step);
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::instantiateMaterials() {}

/* -------------------------------------------------------------------------- */
/* Initialization                                                             */
/* -------------------------------------------------------------------------- */
/**
 * This function groups  many of the initialization in on  function. For most of
 * basics  case the  function should  be  enough. The  functions initialize  the
 * model, the internal  vectors, set them to 0, and  depending on the parameters
 * it also initialize the explicit or implicit solver.
 *
 * @param options
 * \parblock
 * contains the different options to initialize the model
 * \li \c analysis_method specify the type of solver to use
 * \endparblock
 */
void SolidMechanicsModel::initFullImpl(const ModelOptions & options) {
  CLHParent::initFullImpl(options);

  // // initialize the materials
  // if (not this->parser.getLastParsedFile().empty()) {
  //   this->instantiateMaterials();
  //   this->initConstitutiveLaws();
  // }

  this->initBC(*this, *displacement, *displacement_increment, *external_force);
}

/* -------------------------------------------------------------------------- */
TimeStepSolverType SolidMechanicsModel::getDefaultSolverType() const {
  return TimeStepSolverType::_dynamic_lumped;
}

/* -------------------------------------------------------------------------- */
ModelSolverOptions SolidMechanicsModel::getDefaultSolverOptions(
    const TimeStepSolverType & type) const {
  ModelSolverOptions options;

  switch (type) {
  case TimeStepSolverType::_dynamic_lumped: {
    options.non_linear_solver_type = NonLinearSolverType::_lumped;
    options.integration_scheme_type["displacement"] =
        IntegrationSchemeType::_central_difference;
    options.solution_type["displacement"] = IntegrationScheme::_acceleration;
    break;
  }
  case TimeStepSolverType::_static: {
    options.non_linear_solver_type = NonLinearSolverType::_newton_raphson;
    options.integration_scheme_type["displacement"] =
        IntegrationSchemeType::_pseudo_time;
    options.solution_type["displacement"] = IntegrationScheme::_not_defined;
    break;
  }
  case TimeStepSolverType::_dynamic: {
    if (this->method == _explicit_consistent_mass) {
      options.non_linear_solver_type = NonLinearSolverType::_newton_raphson;
      options.integration_scheme_type["displacement"] =
          IntegrationSchemeType::_central_difference;
      options.solution_type["displacement"] = IntegrationScheme::_acceleration;
    } else {
      options.non_linear_solver_type = NonLinearSolverType::_newton_raphson;
      options.integration_scheme_type["displacement"] =
          IntegrationSchemeType::_trapezoidal_rule_2;
      options.solution_type["displacement"] = IntegrationScheme::_displacement;
    }
    break;
  }
  default:
    AKANTU_EXCEPTION(type << " is not a valid time step solver type");
  }

  return options;
}

/* -------------------------------------------------------------------------- */
std::tuple<ID, TimeStepSolverType>
SolidMechanicsModel::getDefaultSolverID(const AnalysisMethod & method) {
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
void SolidMechanicsModel::initSolver(TimeStepSolverType time_step_solver_type) {
  auto & dof_manager = this->getDOFManager();

  /* ------------------------------------------------------------------------ */
  // for alloc type of solvers
  this->allocNodalField(this->displacement, spatial_dimension, "displacement");
  this->allocNodalField(this->previous_displacement, spatial_dimension,
                        "previous_displacement");
  this->allocNodalField(this->displacement_increment, spatial_dimension,
                        "displacement_increment");
  this->allocNodalField(this->internal_force, spatial_dimension,
                        "internal_force");
  this->allocNodalField(this->external_force, spatial_dimension,
                        "external_force");
  this->allocNodalField(this->blocked_dofs, spatial_dimension, "blocked_dofs");
  this->allocNodalField(this->current_position, spatial_dimension,
                        "current_position");

  // initialize the current positions
  this->current_position->copy(this->mesh.getNodes());

  /* ------------------------------------------------------------------------ */
  if (!dof_manager.hasDOFs("displacement")) {
    dof_manager.registerDOFs("displacement", *this->displacement, _dst_nodal);
    dof_manager.registerBlockedDOFs("displacement", *this->blocked_dofs);
    dof_manager.registerDOFsIncrement("displacement",
                                      *this->displacement_increment);
    dof_manager.registerDOFsPrevious("displacement",
                                     *this->previous_displacement);
  }

  /* ------------------------------------------------------------------------ */
  // for dynamic
  if (time_step_solver_type == TimeStepSolverType::_dynamic ||
      time_step_solver_type == TimeStepSolverType::_dynamic_lumped) {
    this->allocNodalField(this->velocity, spatial_dimension, "velocity");
    this->allocNodalField(this->acceleration, spatial_dimension,
                          "acceleration");

    if (!dof_manager.hasDOFsDerivatives("displacement", 1)) {
      dof_manager.registerDOFsDerivative("displacement", 1, *this->velocity);
      dof_manager.registerDOFsDerivative("displacement", 2,
                                         *this->acceleration);
    }
  }
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::assembleResidual() {
  AKANTU_DEBUG_IN();

  /* ------------------------------------------------------------------------ */
  // computes the internal forces
  this->assembleInternalForces();

  /* ------------------------------------------------------------------------ */
  this->getDOFManager().assembleToResidual("displacement",
                                           *this->external_force, 1);
  this->getDOFManager().assembleToResidual("displacement",
                                           *this->internal_force, 1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::assembleResidual(const ID & residual_part) {
  AKANTU_DEBUG_IN();

  if ("external" == residual_part) {
    this->getDOFManager().assembleToResidual("displacement",
                                             *this->external_force, 1);
    AKANTU_DEBUG_OUT();
    return;
  }

  if ("internal" == residual_part) {
    this->assembleInternalForces();
    this->getDOFManager().assembleToResidual("displacement",
                                             *this->internal_force, 1);
    AKANTU_DEBUG_OUT();
    return;
  }

  AKANTU_CUSTOM_EXCEPTION(
      debug::SolverCallbackResidualPartUnknown(residual_part));

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
MatrixType SolidMechanicsModel::getMatrixType(const ID & matrix_id) const {
  // \TODO check the materials to know what is the correct answer
  if (matrix_id == "C") {
    return _mt_not_defined;
  }

  if (matrix_id == "K") {
    auto matrix_type = _unsymmetric;

    for_each_constitutive_law([&](auto && material) {
      matrix_type = std::max(matrix_type, material.getMatrixType(matrix_id));
    });
  }

  return _symmetric;
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::assembleMatrix(const ID & matrix_id) {
  if (matrix_id == "K") {
    this->assembleStiffnessMatrix();
  } else if (matrix_id == "M") {
    this->assembleMass();
  }
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::assembleLumpedMatrix(const ID & matrix_id) {
  if (matrix_id == "M") {
    this->assembleMassLumped();
  }
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::beforeSolveStep() {
  for_each_constitutive_law(
      [&](auto && material) { material.beforeSolveStep(); });
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::afterSolveStep(bool converged) {
  for_each_constitutive_law(
      [&](auto && material) { material.afterSolveStep(converged); });
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::predictor() { ++displacement_release; }

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::corrector() { ++displacement_release; }

/* -------------------------------------------------------------------------- */
/**
 * This function computes the internal forces as \f$F_{int} = \int_{\Omega} N
 * \sigma d\Omega@\f$
 */
void SolidMechanicsModel::assembleInternalForces() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_INFO("Assemble the internal forces");

  this->internal_force->zero();

  // compute the stresses of local elements
  AKANTU_DEBUG_INFO("Compute local stresses");
  for_each_constitutive_law(
      [](auto && material) { material.computeAllStresses(_not_ghost); });

  /* ------------------------------------------------------------------------ */
  /* Computation of the non local part */
  if (this->isNonLocal()) {
    this->getNonLocalManager().computeAllNonLocalContribution();
  }

  // communicate the stresses
  AKANTU_DEBUG_INFO("Send data for residual assembly");
  this->asynchronousSynchronize(SynchronizationTag::_smm_stress);

  // assemble the forces due to local stresses
  AKANTU_DEBUG_INFO("Assemble residual for local elements");
  for_each_constitutive_law(
      [](auto && material) { material.assembleInternalForces(_not_ghost); });

  // finalize communications
  AKANTU_DEBUG_INFO("Wait distant stresses");
  this->waitEndSynchronize(SynchronizationTag::_smm_stress);

  // assemble the stresses due to ghost elements
  AKANTU_DEBUG_INFO("Assemble residual for ghost elements");
  for_each_constitutive_law(
      [](auto && material) { material.assembleInternalForces(_ghost); });

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::assembleStiffnessMatrix(bool need_to_reassemble) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_INFO("Assemble the new stiffness matrix.");

  if (not this->getDOFManager().hasMatrix("K")) {
    this->getDOFManager().getNewMatrix("K", this->getMatrixType("K"));
  }

  // Check if materials need to recompute the matrix
  for_each_constitutive_law([&](auto && material) {
    need_to_reassemble |= material.hasMatrixChanged("K");
  });

  if (need_to_reassemble) {
    this->getDOFManager().getMatrix("K").zero();

    // call compute stiffness matrix on each local elements
    for_each_constitutive_law(
        [](auto && material) { material.assembleStiffnessMatrix(_not_ghost); });
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::updateCurrentPosition() {
  if (this->current_position_release == this->displacement_release) {
    return;
  }

  this->current_position->copy(this->mesh.getNodes());

  for (auto && data : zip(make_view(*this->current_position, spatial_dimension),
                          make_view(*this->displacement, spatial_dimension))) {
    std::get<0>(data) += std::get<1>(data);
  }

  this->current_position_release = this->displacement_release;
}

/* -------------------------------------------------------------------------- */
const Array<Real> & SolidMechanicsModel::getCurrentPosition() {
  this->updateCurrentPosition();
  return *this->current_position;
}

/* -------------------------------------------------------------------------- */
/* Information                                                                */
/* -------------------------------------------------------------------------- */
Real SolidMechanicsModel::getStableTimeStep() {
  AKANTU_DEBUG_IN();

  Real min_dt = getStableTimeStep(_not_ghost);

  /// reduction min over all processors
  mesh.getCommunicator().allReduce(min_dt, SynchronizerOperation::_min);

  AKANTU_DEBUG_OUT();
  return min_dt;
}

/* -------------------------------------------------------------------------- */
Real SolidMechanicsModel::getStableTimeStep(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Real min_dt = std::numeric_limits<Real>::max();

  this->updateCurrentPosition();

  Element elem{_not_defined, 0, ghost_type};
  elem.ghost_type = ghost_type;

  for (auto type :
       mesh.elementTypes(Model::spatial_dimension, ghost_type, _ek_regular)) {
    elem.type = type;
    auto nb_nodes_per_element = mesh.getNbNodesPerElement(type);

    auto mat_indexes = this->getConstitutiveLawByElement(type, ghost_type);
    auto mat_loc_num = this->getConstitutiveLawLocalNumbering(type, ghost_type);

    Array<Real> X(0, nb_nodes_per_element * Model::spatial_dimension);
    FEEngine::extractNodalToElementField(mesh, *current_position, X, type,
                                         _not_ghost);

    for (auto && [X_el, mat_idx, el] :
         zip(make_view(X, spatial_dimension, nb_nodes_per_element),
             make_view(mat_indexes), make_view(mat_loc_num))) {
      elem.element = el;

      auto el_h = FEEngine::getElementInradius(X_el, type);
      auto el_c = this->getMaterial(mat_idx).getCelerity(elem);
      auto el_dt = el_h / el_c;

      min_dt = std::min(min_dt, el_dt);
    }
  }

  AKANTU_DEBUG_OUT();
  return min_dt;
}

/* -------------------------------------------------------------------------- */
Real SolidMechanicsModel::getKineticEnergy() {
  AKANTU_DEBUG_IN();

  Real ekin = 0.;
  auto nb_nodes = mesh.getNbNodes();

  if (this->getDOFManager().hasLumpedMatrix("M")) {
    this->assembleLumpedMatrix("M");

    auto m_it = this->mass->begin(Model::spatial_dimension);
    auto m_end = this->mass->end(Model::spatial_dimension);
    auto v_it = this->velocity->begin(Model::spatial_dimension);

    for (Int n = 0; m_it != m_end; ++n, ++m_it, ++v_it) {
      const auto & v = *v_it;
      const auto & m = *m_it;

      Real mv2 = 0.;
      auto is_local_node = mesh.isLocalOrMasterNode(n);
      // bool is_not_pbc_slave_node = !isPBCSlaveNode(n);
      auto count_node = is_local_node; // && is_not_pbc_slave_node;
      if (count_node) {
        for (Int i = 0; i < spatial_dimension; ++i) {
          if (m(i) > std::numeric_limits<Real>::epsilon()) {
            mv2 += v(i) * v(i) * m(i);
          }
        }
      }

      ekin += mv2;
    }
  } else if (this->getDOFManager().hasMatrix("M")) {
    this->assembleMatrix("M");

    Array<Real> Mv(nb_nodes, Model::spatial_dimension);
    this->getDOFManager().assembleMatMulVectToArray("displacement", "M",
                                                    *this->velocity, Mv);

    for (auto && data : zip(arange(nb_nodes), make_view(Mv, spatial_dimension),
                            make_view(*this->velocity, spatial_dimension))) {
      ekin += std::get<2>(data).dot(std::get<1>(data)) *
              static_cast<Real>(mesh.isLocalOrMasterNode(std::get<0>(data)));
    }
  } else {
    AKANTU_ERROR("No function called to assemble the mass matrix.");
  }

  mesh.getCommunicator().allReduce(ekin, SynchronizerOperation::_sum);

  AKANTU_DEBUG_OUT();
  return ekin / 2.;
}

/* -------------------------------------------------------------------------- */
Real SolidMechanicsModel::getKineticEnergy(const Element & element) {
  AKANTU_DEBUG_IN();

  auto nb_quadrature_points =
      getFEEngine().getNbIntegrationPoints(element.type);

  Array<Real> vel_on_quad(nb_quadrature_points, Model::spatial_dimension);
  Array<Idx> filter_element(1, 1, element.element);

  getFEEngine().interpolateOnIntegrationPoints(
      *velocity, vel_on_quad, Model::spatial_dimension, element.type,
      _not_ghost, filter_element);
  Vector<Real> rho_v2(nb_quadrature_points);
  Real rho = getConstitutiveLaw(element).getRho();

  for (auto && data : enumerate(make_view(vel_on_quad, spatial_dimension))) {
    auto && vel = std::get<1>(data);
    rho_v2(std::get<0>(data)) = rho * vel.dot(vel);
  }

  AKANTU_DEBUG_OUT();

  return .5 * getFEEngine().integrate(rho_v2, element);
}

/* -------------------------------------------------------------------------- */
Real SolidMechanicsModel::getExternalWork() {
  AKANTU_DEBUG_IN();

  Array<Real> * incrs_or_velos{nullptr};
  if (this->method == _static) {
    incrs_or_velos = this->displacement_increment.get();
  } else {
    incrs_or_velos = this->velocity.get();
  }

  Real work = 0.;

  auto nb_nodes = this->mesh.getNbNodes();

  for (auto && [ext_force, int_force, boun, incr_or_velo, n] :
       zip(make_view(*external_force, spatial_dimension),
           make_view(*internal_force, spatial_dimension),
           make_view(*blocked_dofs, spatial_dimension),
           make_view(*incrs_or_velos, spatial_dimension), arange(nb_nodes))) {
    auto is_local_node = this->mesh.isLocalOrMasterNode(n);
    // bool is_not_pbc_slave_node = !this->isPBCSlaveNode(n);
    auto count_node = is_local_node; // && is_not_pbc_slave_node;

    if (count_node) {
      for (Int i = 0; i < spatial_dimension; ++i) {
        if (boun(i)) {
          work -= int_force(i) * incr_or_velo(i);
        } else {
          work += ext_force(i) * incr_or_velo(i);
        }
      }
    }
  }

  mesh.getCommunicator().allReduce(work, SynchronizerOperation::_sum);

  if (this->method != _static) {
    work *= this->getTimeStep();
  }

  AKANTU_DEBUG_OUT();
  return work;
}

/* -------------------------------------------------------------------------- */
Real SolidMechanicsModel::getEnergy(const std::string & energy_id) {
  AKANTU_DEBUG_IN();

  if (energy_id == "kinetic") {
    return getKineticEnergy();
  }

  if (energy_id == "external work") {
    return getExternalWork();
  }

  Real energy = 0.;
  for_each_constitutive_law(
      [&](auto && material) { energy += material.getEnergy(energy_id); });

  /// reduction sum over all processors
  mesh.getCommunicator().allReduce(energy, SynchronizerOperation::_sum);

  AKANTU_DEBUG_OUT();
  return energy;
}

/* -------------------------------------------------------------------------- */
Real SolidMechanicsModel::getEnergy(const std::string & energy_id,
                                    const Element & element) {
  AKANTU_DEBUG_IN();

  if (energy_id == "kinetic") {
    return getKineticEnergy(element);
  }

  auto mat_element = element;
  mat_element.element = this->getConstitutiveLawLocalNumbering()(element);

  Real energy =
      this->getConstitutiveLaw(element).getEnergy(energy_id, mat_element);

  AKANTU_DEBUG_OUT();
  return energy;
}

/* -------------------------------------------------------------------------- */
Real SolidMechanicsModel::getEnergy(const ID & energy_id, const ID & group_id) {
  auto && group = mesh.getElementGroup(group_id);
  auto energy = 0.;
  for (auto && type : group.elementTypes()) {
    for (auto el : group.getElementsIterable(type)) {
      energy += getEnergy(energy_id, el);
    }
  }

  /// reduction sum over all processors
  mesh.getCommunicator().allReduce(energy, SynchronizerOperation::_sum);

  return energy;
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::onNodesAdded(const Array<Idx> & nodes_list,
                                       const NewNodesEvent & event) {
  AKANTU_DEBUG_IN();
  auto nb_nodes = mesh.getNbNodes();

  if (displacement) {
    displacement->resize(nb_nodes, 0.);
    ++displacement_release;
  }
  if (mass) {
    mass->resize(nb_nodes, 0.);
  }
  if (velocity) {
    velocity->resize(nb_nodes, 0.);
  }
  if (acceleration) {
    acceleration->resize(nb_nodes, 0.);
  }
  if (external_force) {
    external_force->resize(nb_nodes, 0.);
  }
  if (internal_force) {
    internal_force->resize(nb_nodes, 0.);
  }
  if (blocked_dofs) {
    blocked_dofs->resize(nb_nodes, false);
  }
  if (current_position) {
    current_position->resize(nb_nodes, 0.);
  }

  if (previous_displacement) {
    previous_displacement->resize(nb_nodes, 0.);
  }
  if (displacement_increment) {
    displacement_increment->resize(nb_nodes, 0.);
  }

  for_each_constitutive_law(
      [&](auto && material) { material.onNodesAdded(nodes_list, event); });

  need_to_reassemble_lumped_mass = true;
  need_to_reassemble_mass = true;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::onNodesRemoved(const Array<Idx> & /*element_list*/,
                                         const Array<Idx> & new_numbering,
                                         const RemovedNodesEvent & /*event*/) {
  if (displacement) {
    mesh.removeNodesFromArray(*displacement, new_numbering);
    ++displacement_release;
  }
  if (mass) {
    mesh.removeNodesFromArray(*mass, new_numbering);
  }
  if (velocity) {
    mesh.removeNodesFromArray(*velocity, new_numbering);
  }
  if (acceleration) {
    mesh.removeNodesFromArray(*acceleration, new_numbering);
  }
  if (internal_force) {
    mesh.removeNodesFromArray(*internal_force, new_numbering);
  }
  if (external_force) {
    mesh.removeNodesFromArray(*external_force, new_numbering);
  }
  if (blocked_dofs) {
    mesh.removeNodesFromArray(*blocked_dofs, new_numbering);
  }

  if (displacement_increment) {
    mesh.removeNodesFromArray(*displacement_increment, new_numbering);
  }

  if (previous_displacement) {
    mesh.removeNodesFromArray(*previous_displacement, new_numbering);
  }
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::printself(std::ostream & stream, int indent) const {
  std::string space(indent, AKANTU_INDENT);

  stream << space << "Solid Mechanics Model ["
         << "\n";
  stream << space << " + id                : " << id << "\n";
  stream << space << " + spatial dimension : " << Model::spatial_dimension
         << "\n";

  stream << space << " + fem ["
         << "\n";
  getFEEngine().printself(stream, indent + 2);
  stream << space << " ]"
         << "\n";

  stream << space << " + nodals information ["
         << "\n";
  displacement->printself(stream, indent + 2);
  if (velocity) {
    velocity->printself(stream, indent + 2);
  }
  if (acceleration) {
    acceleration->printself(stream, indent + 2);
  }
  if (mass) {
    mass->printself(stream, indent + 2);
  }
  external_force->printself(stream, indent + 2);
  internal_force->printself(stream, indent + 2);
  blocked_dofs->printself(stream, indent + 2);
  stream << space << " ]"
         << "\n";

  stream << space << " + materials ["
         << "\n";
  this->for_each_constitutive_law(
      [&](auto && material) { material.printself(stream, indent + 2); });
  stream << space << " ]"
         << "\n";

  stream << space << "]"
         << "\n";
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::computeNonLocalContribution(GhostType ghost_type) {
  for_each_constitutive_law([&](auto && material) {
    if (aka::is_of_type<MaterialNonLocalInterface>(material)) {
      auto & mat_non_local =
          dynamic_cast<MaterialNonLocalInterface &>(material);
      mat_non_local.computeNonLocalStresses(ghost_type);
    }
  });
}

/* -------------------------------------------------------------------------- */
FEEngine & SolidMechanicsModel::getFEEngineBoundary(const ID & name) {
  return getFEEngineClassBoundary<MyFEEngineType>(name);
}

/* -------------------------------------------------------------------------- */
Int SolidMechanicsModel::getNbData(const Array<Element> & elements,
                                   const SynchronizationTag & tag) const {

  Int size = 0;
  Int nb_nodes_per_element = 0;

  for (const Element & el : elements) {
    nb_nodes_per_element += Mesh::getNbNodesPerElement(el.type);
  }

  switch (tag) {
  case SynchronizationTag::_smm_mass: {
    size += spatial_dimension * nb_nodes_per_element * Int(sizeof(Real));
    break;
  }
  case SynchronizationTag::_smm_for_gradu: {
    size += spatial_dimension * nb_nodes_per_element * Int(sizeof(Real));
    break;
  }
  case SynchronizationTag::_for_dump: {
    size += 5 * spatial_dimension * nb_nodes_per_element * Int(sizeof(Real));
    break;
  }
  case SynchronizationTag::_smm_boundary: {
    size += 3 * spatial_dimension * nb_nodes_per_element * Int(sizeof(Real));
    break;
  }
  default: {
  }
  }

  size += CLHParent::getNbData(elements, tag);
  return size;
}

/* --------------------------------------------------------------------------
 */
void SolidMechanicsModel::packData(CommunicationBuffer & buffer,
                                   const Array<Element> & elements,
                                   const SynchronizationTag & tag) const {
  AKANTU_DEBUG_IN();

  switch (tag) {
  case SynchronizationTag::_smm_mass: {
    packNodalDataHelper(*mass, buffer, elements, mesh);
    break;
  }
  case SynchronizationTag::_smm_for_gradu: {
    packNodalDataHelper(*displacement, buffer, elements, mesh);
    break;
  }
  case SynchronizationTag::_for_dump: {
    packNodalDataHelper(*displacement, buffer, elements, mesh);
    packNodalDataHelper(*velocity, buffer, elements, mesh);
    packNodalDataHelper(*acceleration, buffer, elements, mesh);
    packNodalDataHelper(*internal_force, buffer, elements, mesh);
    packNodalDataHelper(*external_force, buffer, elements, mesh);
    break;
  }
  case SynchronizationTag::_smm_boundary: {
    packNodalDataHelper(*external_force, buffer, elements, mesh);
    packNodalDataHelper(*velocity, buffer, elements, mesh);
    packNodalDataHelper(*blocked_dofs, buffer, elements, mesh);
    break;
  }
  default: {
  }
  }

  CLHParent::packData(buffer, elements, tag);

  AKANTU_DEBUG_OUT();
}

/* --------------------------------------------------------------------------
 */
void SolidMechanicsModel::unpackData(CommunicationBuffer & buffer,
                                     const Array<Element> & elements,
                                     const SynchronizationTag & tag) {
  CLHParent::unpackData(buffer, elements, tag);
}

/* --------------------------------------------------------------------------
 */
Int SolidMechanicsModel::getNbData(const Array<Idx> & dofs,
                                   const SynchronizationTag & tag) const {
  Int size = 0;

  switch (tag) {
  case SynchronizationTag::_smm_uv: {
    size += Int(sizeof(Real)) * Model::spatial_dimension * 2;
    break;
  }
  case SynchronizationTag::_smm_res: /* FALLTHRU */
  case SynchronizationTag::_smm_mass: {
    size += Int(sizeof(Real)) * Model::spatial_dimension;
    break;
  }
  case SynchronizationTag::_for_dump: {
    size += Int(sizeof(Real)) * Model::spatial_dimension * 5;
    break;
  }
  default: {
    AKANTU_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }

  return size * dofs.size();
}

/* --------------------------------------------------------------------------
 */
void SolidMechanicsModel::packData(CommunicationBuffer & buffer,
                                   const Array<Idx> & dofs,
                                   const SynchronizationTag & tag) const {
  AKANTU_DEBUG_IN();

  switch (tag) {
  case SynchronizationTag::_smm_uv: {
    packDOFDataHelper(*displacement, buffer, dofs);
    packDOFDataHelper(*velocity, buffer, dofs);
    break;
  }
  case SynchronizationTag::_smm_res: {
    packDOFDataHelper(*internal_force, buffer, dofs);
    break;
  }
  case SynchronizationTag::_smm_mass: {
    packDOFDataHelper(*mass, buffer, dofs);
    break;
  }
  case SynchronizationTag::_for_dump: {
    packDOFDataHelper(*displacement, buffer, dofs);
    packDOFDataHelper(*velocity, buffer, dofs);
    packDOFDataHelper(*acceleration, buffer, dofs);
    packDOFDataHelper(*internal_force, buffer, dofs);
    packDOFDataHelper(*external_force, buffer, dofs);
    break;
  }
  default: {
    AKANTU_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }

  AKANTU_DEBUG_OUT();
}

/* --------------------------------------------------------------------------
 */
void SolidMechanicsModel::unpackData(CommunicationBuffer & buffer,
                                     const Array<Idx> & dofs,
                                     const SynchronizationTag & tag) {
  AKANTU_DEBUG_IN();

  switch (tag) {
  case SynchronizationTag::_smm_uv: {
    unpackDOFDataHelper(*displacement, buffer, dofs);
    unpackDOFDataHelper(*velocity, buffer, dofs);
    break;
  }
  case SynchronizationTag::_smm_res: {
    unpackDOFDataHelper(*internal_force, buffer, dofs);
    break;
  }
  case SynchronizationTag::_smm_mass: {
    unpackDOFDataHelper(*mass, buffer, dofs);
    break;
  }
  case SynchronizationTag::_for_dump: {
    unpackDOFDataHelper(*displacement, buffer, dofs);
    unpackDOFDataHelper(*velocity, buffer, dofs);
    unpackDOFDataHelper(*acceleration, buffer, dofs);
    unpackDOFDataHelper(*internal_force, buffer, dofs);
    unpackDOFDataHelper(*external_force, buffer, dofs);
    break;
  }
  default: {
    AKANTU_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }

  AKANTU_DEBUG_OUT();
}

/* --------------------------------------------------------------------------
 */
void SolidMechanicsModel::applyEigenGradU(
    const Matrix<Real> & prescribed_eigen_grad_u, const ID & material_name,
    const GhostType ghost_type) {
  AKANTU_DEBUG_ASSERT(prescribed_eigen_grad_u.size() ==
                          spatial_dimension * spatial_dimension,
                      "The prescribed grad_u is not of the good size");
  for_each_constitutive_law([&](auto && material) {
    if (material.getName() == material_name) {
      material.applyEigenGradU(prescribed_eigen_grad_u, ghost_type);
    }
  });
}

/* --------------------------------------------------------------------------
 */
void SolidMechanicsModel::registerNewMaterial(const ID & mat_name,
                                              const ID & mat_type,
                                              const ID & opt_param) {
  this->registerNewConstitutiveLaw(mat_name, mat_type, opt_param);
}

/* --------------------------------------------------------------------------
 */

} // namespace akantu
