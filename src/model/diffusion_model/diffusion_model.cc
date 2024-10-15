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
#include "heat_transfer_model.hh"
// #include "dumpable_inline_impl.hh"
#include "element_synchronizer.hh"
#include "fe_engine_template.hh"
// #include "generalized_trapezoidal.hh"
#include "diffusion_law.hh"
#include "group_manager_inline_impl.hh"
#include "integrator_gauss.hh"
#include "mesh.hh"
#include "parser.hh"
#include "shape_lagrange.hh"
/* -------------------------------------------------------------------------- */
#include "dumper_element_partition.hh"

// #include "dumper_elemental_field.hh"
#include "dumper_internal_material_field.hh"
#include "dumper_iohelper_paraview.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

namespace diffusion {
  namespace details {
    class ComputeRhoFunctor {
    public:
      ComputeRhoFunctor(const DiffusionModel & model) : model(model){};

      void operator()(Matrix<Real> & rho, const Element & element) {
        auto law_id = model.getConstitutiveLawByElement()(element);
        rho.array() = model.getConstitutiveLaw(law_id).getRho();
      }

    private:
      const DiffusionModel & model;
    };
  } // namespace details
} // namespace diffusion

/* -------------------------------------------------------------------------- */
DiffusionModel::DiffusionModel(Mesh & mesh, Int dim, const ID & id,
                               const std::shared_ptr<DOFManager> & dof_manager,
                               const ID & dof_name, ModelType model_type)
    : Parent(mesh, model_type, dim, id), dof_name(dof_name) {
  AKANTU_DEBUG_IN();

  this->initDOFManager(dof_manager);

  if (this->mesh.isDistributed()) {
    auto & synchronizer = this->mesh.getElementSynchronizer();
    this->registerSynchronizer(synchronizer, SynchronizationTag::_diffusion);
    this->registerSynchronizer(synchronizer,
                               SynchronizationTag::_diffusion_gradient);
  }

  registerFEEngineObject<FEEngineType>(id + ":fem", mesh, spatial_dimension);

  this->mesh.registerDumper<DumperParaview>(id, id, true);
  this->mesh.addDumpMesh(mesh, spatial_dimension, _not_ghost, _ek_regular);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
FEEngine & DiffusionModel::getFEEngineBoundary(const ID & name) {
  return aka::as_type<FEEngine>(getFEEngineClassBoundary<FEEngineType>(name));
}

/* -------------------------------------------------------------------------- */
MatrixType DiffusionModel::getMatrixType(const ID & matrix_id) const {
  if (matrix_id == "K" or matrix_id == "M") {
    return _symmetric;
  }

  return _mt_not_defined;
}

/* -------------------------------------------------------------------------- */
void DiffusionModel::assembleMatrix(const ID & matrix_id) {
  if (matrix_id == "K") {
    this->assembleDiffisivityMatrix();
  } else if (matrix_id == "M" and need_to_reassemble_capacity) {
    this->assembleCapacityMatrix();
  }
}

/* -------------------------------------------------------------------------- */
void DiffusionModel::assembleLumpedMatrix(const ID & matrix_id) {
  if (matrix_id == "M" and need_to_reassemble_capacity) {
    this->assembleCapacityMatrixLumped();
  }
}

/* -------------------------------------------------------------------------- */
void DiffusionModel::assembleResidual() {
  this->assembleInternalFlow();

  this->getDOFManager().assembleToResidual(dof_name, *this->external_flow, 1);
  this->getDOFManager().assembleToResidual(dof_name, *this->internal_flow, 1);
}

/* -------------------------------------------------------------------------- */
void DiffusionModel::predictor() { ++diffusion_release; }
/* -------------------------------------------------------------------------- */
void DiffusionModel::corrector() { ++diffusion_release; }

/* -------------------------------------------------------------------------- */
void DiffusionModel::initSolver(TimeStepSolverType time_step_solver_type) {
  DOFManager & dof_manager = this->getDOFManager();

  this->allocNodalField(this->diffusion, 1, dof_name);
  this->allocNodalField(this->external_flow, 1, "external_flow");
  this->allocNodalField(this->internal_flow, 1, "internal_flow");
  this->allocNodalField(this->blocked_dofs, 1, "blocked_dofs");

  if (not dof_manager.hasDOFs(dof_name)) {
    dof_manager.registerDOFs(dof_name, *this->diffusion, _dst_nodal);
    dof_manager.registerBlockedDOFs(dof_name, *this->blocked_dofs);
  }

  if (time_step_solver_type == TimeStepSolverType::_dynamic ||
      time_step_solver_type == TimeStepSolverType::_dynamic_lumped) {
    this->allocNodalField(this->diffusion_rate, 1, dof_name + "_rate");

    if (not dof_manager.hasDOFsDerivatives(dof_name, 1)) {
      dof_manager.registerDOFsDerivative(dof_name, 1, *this->diffusion_rate);
    }
  }
}

/* -------------------------------------------------------------------------- */
std::tuple<ID, TimeStepSolverType>
DiffusionModel::getDefaultSolverID(const AnalysisMethod & method) {
  switch (method) {
  case _explicit_lumped_mass: {
    return std::make_tuple("explicit_lumped",
                           TimeStepSolverType::_dynamic_lumped);
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
ModelSolverOptions
DiffusionModel::getDefaultSolverOptions(const TimeStepSolverType & type) const {
  ModelSolverOptions options;

  switch (type) {
  case TimeStepSolverType::_dynamic_lumped: {
    options.non_linear_solver_type = NonLinearSolverType::_lumped;
    options.integration_scheme_type[dof_name] =
        IntegrationSchemeType::_forward_euler;
    options.solution_type[dof_name] = IntegrationScheme::_temperature_rate;
    break;
  }
  case TimeStepSolverType::_static: {
    options.non_linear_solver_type = NonLinearSolverType::_newton_raphson;
    options.integration_scheme_type[dof_name] =
        IntegrationSchemeType::_pseudo_time;
    options.solution_type[dof_name] = IntegrationScheme::_not_defined;
    break;
  }
  case TimeStepSolverType::_dynamic: {
    if (this->method == _explicit_consistent_mass) {
      options.non_linear_solver_type = NonLinearSolverType::_newton_raphson;
      options.integration_scheme_type[dof_name] =
          IntegrationSchemeType::_forward_euler;
      options.solution_type[dof_name] = IntegrationScheme::_temperature_rate;
    } else {
      options.non_linear_solver_type = NonLinearSolverType::_newton_raphson;
      options.integration_scheme_type[dof_name] =
          IntegrationSchemeType::_backward_euler;
      options.solution_type[dof_name] = IntegrationScheme::_temperature;
    }
    break;
  }
  default:
    AKANTU_EXCEPTION(type << " is not a valid time step solver type");
  }

  return options;
}

/* -------------------------------------------------------------------------- */
void DiffusionModel::assembleDiffisivityMatrix() {
  AKANTU_DEBUG_ASSERT(this->getDOFManager().hasMatrix("K"),
                      "The K matrix has not been initialized yet.");

  this->getDOFManager().zeroMatrix("K");

  for_each_constitutive_law(
      [](auto && diffusion_law) { diffusion_law.assembleDiffusivityMatrix(); });
}

/* -------------------------------------------------------------------------- */
void DiffusionModel::assembleInternalFlow() {
  this->internal_flow->zero();

  this->synchronize(SynchronizationTag::_diffusion);

  for (auto ghost_type : ghost_types) {
    for_each_constitutive_law([&](auto && diffusion_law) {
      diffusion_law.computeDiffusivityGradU(ghost_type);
    });

    for_each_constitutive_law([&](auto && diffusion_law) {
      diffusion_law.assembleInternalFlow(ghost_type);
    });
  }
}

/* -------------------------------------------------------------------------- */
auto DiffusionModel::getStableTimeStep() -> Real {
  AKANTU_DEBUG_IN();

  Real el_size{};
  Real min_el_size = std::numeric_limits<Real>::max();

  for (auto && type : mesh.elementTypes(spatial_dimension, _not_ghost)) {
    auto nb_nodes_per_element = mesh.getNbNodesPerElement(type);
    Array<Real> coord(0, nb_nodes_per_element * spatial_dimension);
    FEEngine::extractNodalToElementField(mesh, mesh.getNodes(), coord, type,
                                         _not_ghost);

    for (auto && el_coord :
         make_view(coord, spatial_dimension, nb_nodes_per_element)) {
      el_size = FEEngine::getElementInradius(el_coord, type);
      min_el_size = std::min(min_el_size, el_size);
    }
  }

  Real min_dt = std::numeric_limits<Real>::max();
  for_each_constitutive_law([&](auto && diffusion_law) {
    min_dt = std::min(diffusion_law.getStableTimeStep(min_el_size), min_dt);
  });

  mesh.getCommunicator().allReduce(min_dt, SynchronizerOperation::_min);

  AKANTU_DEBUG_OUT();

  return min_dt;
}
/* -------------------------------------------------------------------------- */

void DiffusionModel::setTimeStep(Real time_step, const ID & solver_id) {
  Model::setTimeStep(time_step, solver_id);

  this->mesh.getDumper(id).setTimeStep(time_step);
}

/* -------------------------------------------------------------------------- */
void DiffusionModel::assembleCapacityMatrixLumped() {
  AKANTU_DEBUG_IN();

  if (!this->getDOFManager().hasLumpedMatrix("M")) {
    this->getDOFManager().getNewLumpedMatrix("M");
  }

  this->getDOFManager().zeroLumpedMatrix("M");

  for (auto && ghost_type : ghost_types) {
    auto & fem = getFEEngineClass<FEEngineType>();
    diffusion::details::ComputeRhoFunctor compute_rho(*this);

    for (auto && type :
         mesh.elementTypes(spatial_dimension, ghost_type, _ek_regular)) {
      fem.assembleFieldLumped(compute_rho, "M", dof_name, this->getDOFManager(),
                              type, ghost_type);
    }
  }

  need_to_reassemble_capacity_lumped = false;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DiffusionModel::assembleCapacityMatrix() {
  AKANTU_DEBUG_IN();
  auto ghost_type = _not_ghost;

  this->getDOFManager().zeroMatrix("M");

  auto & fem = getFEEngineClass<FEEngineType>();

  diffusion::details::ComputeRhoFunctor rho_functor(*this);

  for (auto && type :
       mesh.elementTypes(spatial_dimension, ghost_type, _ek_regular)) {
    fem.assembleFieldMatrix(rho_functor, "M", dof_name, this->getDOFManager(),
                            type, ghost_type);
  }

  need_to_reassemble_capacity = false;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
auto DiffusionModel::computeThermalEnergyByNode() -> Real {
  AKANTU_DEBUG_IN();
  auto time_step = this->getTimeStep();
  Real ethermal = 0.;

  for (auto && [n, flow] :
       enumerate(make_view(*internal_flow, internal_flow->getNbComponent()))) {
    Real E = 0.;
    bool count_node = mesh.isLocalOrMasterNode(n);

    if (count_node) {
      E += (flow.array() * time_step).sum();
    }

    ethermal += E;
  }

  mesh.getCommunicator().allReduce(ethermal, SynchronizerOperation::_sum);

  AKANTU_DEBUG_OUT();
  return ethermal;
}

/* -------------------------------------------------------------------------- */
auto DiffusionModel::getEnergy(const ID & energy_id) -> Real {
  AKANTU_DEBUG_IN();
  Real energy = 0;

  for_each_constitutive_law([&energy, &energy_id](auto && constitutive_law) {
    energy += constitutive_law.getEnergy(energy_id);
  });

  // reduction sum over all processors
  mesh.getCommunicator().allReduce(energy, SynchronizerOperation::_sum);

  AKANTU_DEBUG_OUT();
  return energy;
}

/* -------------------------------------------------------------------------- */
auto DiffusionModel::getEnergy(const ID & energy_id, const Element & element)
    -> Real {
  auto constitutive_law_element = element;
  constitutive_law_element.element =
      this->getConstitutiveLawLocalNumbering()(element);

  return this->getConstitutiveLaw(element).getEnergy(energy_id,
                                                     constitutive_law_element);
}

/* -------------------------------------------------------------------------- */
std::shared_ptr<dumpers::Field>
DiffusionModel::createNodalFieldBool(const std::string & field_name,
                                     const std::string & group_name,
                                     bool /*padding_flag*/) {

  std::map<std::string, Array<bool> *> uint_nodal_fields;
  uint_nodal_fields["blocked_dofs"] = blocked_dofs.get();

  auto field = mesh.createNodalField(uint_nodal_fields[field_name], group_name);
  return field;
}

/* -------------------------------------------------------------------------- */
std::shared_ptr<dumpers::Field> DiffusionModel::createNodalFieldReal(
    const std::string & field_name, const std::string & group_name,
    __attribute__((unused)) bool padding_flag) {

  if (field_name == "capacity_lumped") {
    AKANTU_EXCEPTION(
        "Capacity lumped is a nodal field now stored in the DOF manager."
        "Therefore it cannot be used by a dumper anymore");
  }

  std::map<std::string, Array<Real> *> real_nodal_fields;
  real_nodal_fields[dof_name] = diffusion.get();
  real_nodal_fields[dof_name + "_rate"] = diffusion_rate.get();
  real_nodal_fields["external_flow"] = external_flow.get();
  real_nodal_fields["internal_flow"] = internal_flow.get();
  real_nodal_fields["increment"] = increment.get();

  if (auto it = real_nodal_fields.find(field_name);
      it != real_nodal_fields.end()) {
    return mesh.createNodalField(it->second, group_name);
  }

  return nullptr;
}

/* -------------------------------------------------------------------------- */
std::shared_ptr<dumpers::Field> DiffusionModel::createElementalField(
    const std::string & field_name, const std::string & group_name,
    bool /*padding_flag*/, Int /*spatial_dimension*/,
    ElementKind element_kind) {

  std::shared_ptr<dumpers::Field> field;

  if (field_name == "partitions") {
    field = mesh.createElementalField<Int, dumpers::ElementPartitionField>(
        mesh.getConnectivities(), group_name, this->spatial_dimension,
        element_kind);
  }
  bool is_internal = this->isInternal(field_name, element_kind);

  if (is_internal) {
    auto nb_data_per_elem =
        this->getInternalDataPerElem(field_name, element_kind);
    auto & internal_flat = this->flattenInternal(field_name, element_kind);

    field = mesh.createElementalField<Real, dumpers::InternalMaterialField>(
        internal_flat, group_name, spatial_dimension, element_kind,
        nb_data_per_elem);

    // homogenize the field
    auto foo = dumpers::HomogenizerProxy::createHomogenizer(*field);

    field =
        dumpers::FieldComputeProxy::createFieldCompute(field, std::move(foo));
  };

  return field;
}

/* -------------------------------------------------------------------------- */
Int DiffusionModel::getNbData(const Array<Idx> & indexes,
                              const SynchronizationTag & tag) const {
  AKANTU_DEBUG_IN();

  Int size = 0;
  auto nb_nodes = indexes.size();

  switch (tag) {
  case SynchronizationTag::_diffusion: {
    size += nb_nodes * Int(sizeof(Real));
    break;
  }
  default: {
    AKANTU_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }

  AKANTU_DEBUG_OUT();
  return size;
}

/* -------------------------------------------------------------------------- */
void DiffusionModel::packData(CommunicationBuffer & buffer,
                              const Array<Idx> & indexes,
                              const SynchronizationTag & tag) const {
  AKANTU_DEBUG_IN();

  for (auto index : indexes) {
    switch (tag) {
    case SynchronizationTag::_diffusion: {
      buffer << (*diffusion)(index);
      break;
    }
    default: {
      AKANTU_ERROR("Unknown ghost synchronization tag : " << tag);
    }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DiffusionModel::unpackData(CommunicationBuffer & buffer,
                                const Array<Idx> & indexes,
                                const SynchronizationTag & tag) {
  AKANTU_DEBUG_IN();

  for (auto index : indexes) {
    switch (tag) {
    case SynchronizationTag::_diffusion: {
      buffer >> (*diffusion)(index);
      break;
    }
    default: {
      AKANTU_ERROR("Unknown ghost synchronization tag : " << tag);
    }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Int DiffusionModel::getNbData(const Array<Element> & elements,
                              const SynchronizationTag & tag) const {
  AKANTU_DEBUG_IN();

  Int size = 0;
  auto nb_nodes_per_element = 0;

  for (auto && el : elements) {
    nb_nodes_per_element += Mesh::getNbNodesPerElement(el.type);
  }

  if (tag == SynchronizationTag::_diffusion) {
    size += nb_nodes_per_element * Int(sizeof(Real)); // temperature
  }

  size += Parent::getNbData(elements, tag);

  AKANTU_DEBUG_OUT();
  return size;
}

/* -------------------------------------------------------------------------- */
void DiffusionModel::packData(CommunicationBuffer & buffer,
                              const Array<Element> & elements,
                              const SynchronizationTag & tag) const {
  if (tag == SynchronizationTag::_diffusion) {
    packNodalDataHelper(*diffusion, buffer, elements, mesh);
  }

  Parent::packData(buffer, elements, tag);
}

/* -------------------------------------------------------------------------- */
void DiffusionModel::unpackData(CommunicationBuffer & buffer,
                                const Array<Element> & elements,
                                const SynchronizationTag & tag) {
  if (tag == SynchronizationTag::_diffusion) {
    unpackNodalDataHelper(*diffusion, buffer, elements, mesh);
  }

  Parent::unpackData(buffer, elements, tag);
}

/* -------------------------------------------------------------------------- */
} // namespace akantu
