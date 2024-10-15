/**
 * Copyright (©) 2015-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "dof_manager.hh"
#include "communicator.hh"
#include "mesh.hh"
#include "mesh_utils.hh"
#include "node_group.hh"
#include "node_synchronizer.hh"
#include "non_linear_solver.hh"
#include "periodic_node_synchronizer.hh"
#include "time_step_solver.hh"
/* -------------------------------------------------------------------------- */
#include <memory>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
DOFManager::DOFManager(const ID & id)
    : id(id), dofs_flag(0, 1, std::string(id + ":dofs_type")),
      global_equation_number(0, 1, "global_equation_number"),
      communicator(Communicator::getStaticCommunicator()) {}

/* -------------------------------------------------------------------------- */
DOFManager::DOFManager(Mesh & mesh, const ID & id)
    : id(id), mesh(&mesh), dofs_flag(0, 1, std::string(id + ":dofs_type")),
      global_equation_number(0, 1, "global_equation_number"),
      communicator(mesh.getCommunicator()) {
  this->mesh->registerEventHandler(*this, _ehp_dof_manager);
}

/* -------------------------------------------------------------------------- */
DOFManager::~DOFManager() {
  if (this->mesh) {
    this->mesh->unregisterEventHandler(*this);
  }
}

/* -------------------------------------------------------------------------- */
std::vector<ID> DOFManager::getDOFIDs() const {
  std::vector<ID> keys;
  for (const auto & dof_data : this->dofs) {
    keys.push_back(dof_data.first);
  }

  return keys;
}

/* -------------------------------------------------------------------------- */
void DOFManager::assembleElementalArrayLocalArray(
    const Array<Real> & elementary_vect, Array<Real> & array_assembeled,
    ElementType type, GhostType ghost_type, Real scale_factor,
    const Array<Int> & filter_elements) {
  AKANTU_DEBUG_IN();

  Int nb_element;
  auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  auto nb_degree_of_freedom =
      elementary_vect.getNbComponent() / nb_nodes_per_element;

  const Idx * filter_it = nullptr;
  if (filter_elements != empty_filter) {
    nb_element = filter_elements.size();
    filter_it = filter_elements.data();
  } else {
    nb_element = this->mesh->getNbElement(type, ghost_type);
  }

  AKANTU_DEBUG_ASSERT(elementary_vect.size() == nb_element,
                      "The vector elementary_vect("
                          << elementary_vect.getID()
                          << ") has not the good size.");

  const auto & connectivity = this->mesh->getConnectivity(type, ghost_type);

  auto elem_it =
      make_view(elementary_vect, nb_degree_of_freedom, nb_nodes_per_element)
          .begin();
  auto assemble_it = make_view(array_assembeled, nb_degree_of_freedom).begin();

  for (Int el = 0; el < nb_element; ++el, ++elem_it) {
    auto element = el;
    if (filter_it != nullptr) {

      element = *filter_it;
    }

    // const Vector<UInt> & conn = *conn_it;
    const auto & elemental_val = *elem_it;
    for (Int n = 0; n < nb_nodes_per_element; ++n) {
      auto node = connectivity(element, n);

      auto && assemble = assemble_it[node];
      assemble += scale_factor * elemental_val(n);
    }

    if (filter_it != nullptr) {
      ++filter_it;
    }
    //    else
    //      ++conn_it;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DOFManager::assembleElementalArrayToResidual(
    const ID & dof_id, const Array<Real> & elementary_vect, ElementType type,
    GhostType ghost_type, Real scale_factor,
    const Array<Int> & filter_elements) {
  AKANTU_DEBUG_IN();

  auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  auto nb_degree_of_freedom =
      elementary_vect.getNbComponent() / nb_nodes_per_element;
  Array<Real> array_localy_assembeled(this->mesh->getNbNodes(),
                                      nb_degree_of_freedom);

  array_localy_assembeled.zero();

  this->assembleElementalArrayLocalArray(
      elementary_vect, array_localy_assembeled, type, ghost_type, scale_factor,
      filter_elements);

  this->assembleToResidual(dof_id, array_localy_assembeled, 1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DOFManager::assembleElementalArrayToLumpedMatrix(
    const ID & dof_id, const Array<Real> & elementary_vect,
    const ID & lumped_mtx, ElementType type, GhostType ghost_type,
    Real scale_factor, const Array<Int> & filter_elements) {
  AKANTU_DEBUG_IN();

  auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  auto nb_degree_of_freedom =
      elementary_vect.getNbComponent() / nb_nodes_per_element;
  Array<Real> array_localy_assembeled(this->mesh->getNbNodes(),
                                      nb_degree_of_freedom);

  array_localy_assembeled.zero();

  this->assembleElementalArrayLocalArray(
      elementary_vect, array_localy_assembeled, type, ghost_type, scale_factor,
      filter_elements);

  this->assembleToLumpedMatrix(dof_id, array_localy_assembeled, lumped_mtx, 1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DOFManager::assembleMatMulDOFsToResidual(const ID & A_id,
                                              Real scale_factor) {
  for (auto & pair : this->dofs) {
    const auto & dof_id = pair.first;
    auto & dof_data = *pair.second;

    this->assembleMatMulVectToResidual(dof_id, A_id, *dof_data.dof,
                                       scale_factor);
  }
}

/* -------------------------------------------------------------------------- */
void DOFManager::splitSolutionPerDOFs() {
  for (auto && [name, dof_data] : this->dofs) {
    dof_data->solution.resize(dof_data->dof->size() *
                              dof_data->dof->getNbComponent());
    this->getSolutionPerDOFs(name, dof_data->solution);
  }
}

/* -------------------------------------------------------------------------- */
void DOFManager::getSolutionPerDOFs(const ID & dof_id,
                                    Array<Real> & solution_array) {
  AKANTU_DEBUG_IN();
  this->getArrayPerDOFs(dof_id, this->getSolution(), solution_array);
  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
void DOFManager::getLumpedMatrixPerDOFs(const ID & dof_id,
                                        const ID & lumped_mtx,
                                        Array<Real> & lumped) {
  AKANTU_DEBUG_IN();
  this->getArrayPerDOFs(dof_id, this->getLumpedMatrix(lumped_mtx), lumped);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DOFManager::assembleToResidual(const ID & dof_id,
                                    Array<Real> & array_to_assemble,
                                    Real scale_factor) {
  AKANTU_DEBUG_IN();

  // this->makeConsistentForPeriodicity(dof_id, array_to_assemble);
  this->assembleToGlobalArray(dof_id, array_to_assemble, this->getResidual(),
                              scale_factor);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DOFManager::assembleToLumpedMatrix(const ID & dof_id,
                                        Array<Real> & array_to_assemble,
                                        const ID & lumped_mtx,
                                        Real scale_factor) {
  AKANTU_DEBUG_IN();

  // this->makeConsistentForPeriodicity(dof_id, array_to_assemble);
  auto & lumped = this->getLumpedMatrix(lumped_mtx);
  this->assembleToGlobalArray(dof_id, array_to_assemble, lumped, scale_factor);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
DOFManager::DOFData::DOFData(const ID & dof_id)
    : support_type(_dst_generic), group_support("__mesh__"),
      solution(0, 1, dof_id + ":solution"),
      local_equation_number(0, 1, dof_id + ":local_equation_number"),
      associated_nodes(0, 1, dof_id + "associated_nodes") {}

/* -------------------------------------------------------------------------- */
DOFManager::DOFData::~DOFData() = default;

/* -------------------------------------------------------------------------- */
template <typename Func>
auto DOFManager::countDOFsForNodes(const DOFData & dof_data, Int nb_nodes,
                                   Func && getNode) {
  auto nb_local_dofs = nb_nodes;
  decltype(nb_local_dofs) nb_pure_local = 0;
  for (auto n : arange(nb_nodes)) {
    UInt node = getNode(n);

    // http://www.open-std.org/jtc1/sc22/open/n2356/conv.html
    // bool are by convention casted to 0 and 1 when promoted to int
    nb_pure_local += this->mesh->isLocalOrMasterNode(node);
    nb_local_dofs -= this->mesh->isPeriodicSlave(node);
  }

  const auto & dofs_array = *dof_data.dof;
  nb_pure_local *= dofs_array.getNbComponent();
  nb_local_dofs *= dofs_array.getNbComponent();
  return std::make_pair(nb_local_dofs, nb_pure_local);
}

/* -------------------------------------------------------------------------- */
auto DOFManager::getNewDOFDataInternal(const ID & dof_id) -> DOFData & {
  auto it = this->dofs.find(dof_id);
  if (it != this->dofs.end()) {
    AKANTU_EXCEPTION("This dof array has already been registered");
  }

  std::unique_ptr<DOFData> dof_data_ptr = this->getNewDOFData(dof_id);
  DOFData & dof_data = *dof_data_ptr;

  this->dofs[dof_id] = std::move(dof_data_ptr);
  return dof_data;
}

/* -------------------------------------------------------------------------- */
void DOFManager::registerDOFs(const ID & dof_id, Array<Real> & dofs_array,
                              DOFSupportType support_type) {
  auto & dofs_storage = this->getNewDOFDataInternal(dof_id);
  dofs_storage.support_type = support_type;

  this->registerDOFsInternal(dof_id, dofs_array);

  resizeGlobalArrays();
}

/* -------------------------------------------------------------------------- */
void DOFManager::registerDOFs(const ID & dof_id, Array<Real> & dofs_array,
                              const ID & support_group) {
  auto & dofs_storage = this->getNewDOFDataInternal(dof_id);
  dofs_storage.support_type = _dst_nodal;
  dofs_storage.group_support = support_group;

  this->registerDOFsInternal(dof_id, dofs_array);

  resizeGlobalArrays();
}

/* -------------------------------------------------------------------------- */
std::tuple<Int, Int, Int>
DOFManager::registerDOFsInternal(const ID & dof_id, Array<Real> & dofs_array) {
  DOFData & dof_data = this->getDOFData(dof_id);
  dof_data.dof = &dofs_array;

  Int nb_local_dofs = 0;
  Int nb_pure_local = 0;

  const auto & support_type = dof_data.support_type;

  switch (support_type) {
  case _dst_nodal: {
    const auto & group = dof_data.group_support;

    std::function<Idx(Idx)> getNode;
    if (group == "__mesh__") {
      AKANTU_DEBUG_ASSERT(
          dofs_array.size() == this->mesh->getNbNodes(),
          "The array of dof is too short to be associated to nodes.");

      std::tie(nb_local_dofs, nb_pure_local) = countDOFsForNodes(
          dof_data, this->mesh->getNbNodes(), [](auto && n) { return n; });
    } else {
      const auto & node_group =
          this->mesh->getElementGroup(group).getNodeGroup().getNodes();

      AKANTU_DEBUG_ASSERT(
          dofs_array.size() == node_group.size(),
          "The array of dof is too short to be associated to nodes.");

      std::tie(nb_local_dofs, nb_pure_local) =
          countDOFsForNodes(dof_data, node_group.size(),
                            [&node_group](auto && n) { return node_group(n); });
    }

    break;
  }
  case _dst_generic: {
    nb_local_dofs = nb_pure_local =
        dofs_array.size() * dofs_array.getNbComponent();
    break;
  }
  default: {
    AKANTU_EXCEPTION("This type of dofs is not handled yet.");
  }
  }

  this->pure_local_system_size += nb_pure_local;
  this->local_system_size += nb_local_dofs;

  auto nb_total_pure_local = nb_pure_local;
  communicator.allReduce(nb_total_pure_local, SynchronizerOperation::_sum);

  this->system_size += nb_total_pure_local;

  // updating the dofs data after counting is finished
  switch (support_type) {
  case _dst_nodal: {
    const auto & group = dof_data.group_support;
    if (group != "__mesh__") {
      auto & support_nodes =
          this->mesh->getElementGroup(group).getNodeGroup().getNodes();
      this->updateDOFsData(
          dof_data, nb_local_dofs, nb_pure_local, support_nodes.size(),
          [&support_nodes](Idx node) -> Idx { return support_nodes[node]; });
    } else {
      this->updateDOFsData(dof_data, nb_local_dofs, nb_pure_local,
                           mesh->getNbNodes(),
                           [](Idx node) -> Idx { return node; });
    }
    break;
  }
  case _dst_generic: {
    this->updateDOFsData(dof_data, nb_local_dofs, nb_pure_local);
    break;
  }
  }

  return std::make_tuple(nb_local_dofs, nb_pure_local, nb_total_pure_local);
}

/* -------------------------------------------------------------------------- */
void DOFManager::registerDOFsPrevious(const ID & dof_id, Array<Real> & array) {
  DOFData & dof = this->getDOFData(dof_id);

  if (dof.previous != nullptr) {
    AKANTU_EXCEPTION("The previous dofs array for "
                     << dof_id << " has already been registered");
  }

  dof.previous = &array;
}

/* -------------------------------------------------------------------------- */
void DOFManager::registerDOFsIncrement(const ID & dof_id, Array<Real> & array) {
  DOFData & dof = this->getDOFData(dof_id);

  if (dof.increment != nullptr) {
    AKANTU_EXCEPTION("The dofs increment array for "
                     << dof_id << " has already been registered");
  }

  dof.increment = &array;
}

/* -------------------------------------------------------------------------- */
void DOFManager::registerDOFsDerivative(const ID & dof_id, Int order,
                                        Array<Real> & dofs_derivative) {
  auto & dof = this->getDOFData(dof_id);
  auto & derivatives = dof.dof_derivatives;

  if (Int(derivatives.size()) < order) {
    derivatives.resize(order, nullptr);
  } else {
    if (derivatives[order - 1] != nullptr) {
      AKANTU_EXCEPTION("The dof derivatives of order "
                       << order << " already been registered for this dof ("
                       << dof_id << ")");
    }
  }

  derivatives[order - 1] = &dofs_derivative;
}

/* -------------------------------------------------------------------------- */
void DOFManager::registerBlockedDOFs(const ID & dof_id,
                                     Array<bool> & blocked_dofs) {
  auto & dof = this->getDOFData(dof_id);

  if (dof.blocked_dofs != nullptr) {
    AKANTU_EXCEPTION("The blocked dofs array for "
                     << dof_id << " has already been registered");
  }

  dof.blocked_dofs = &blocked_dofs;
}

/* -------------------------------------------------------------------------- */
SparseMatrix &
DOFManager::registerSparseMatrix(const ID & matrix_id,
                                 std::unique_ptr<SparseMatrix> & matrix) {
  auto it = this->matrices.find(matrix_id);
  if (it != this->matrices.end()) {
    AKANTU_EXCEPTION("The matrix " << matrix_id << " already exists in "
                                   << this->id);
  }

  auto & ret = *matrix;
  this->matrices[matrix_id] = std::move(matrix);
  return ret;
}

/* -------------------------------------------------------------------------- */
/// Get an instance of a new SparseMatrix
SolverVector &
DOFManager::registerLumpedMatrix(const ID & matrix_id,
                                 std::unique_ptr<SolverVector> & matrix) {
  auto it = this->lumped_matrices.find(matrix_id);
  if (it != this->lumped_matrices.end()) {
    AKANTU_EXCEPTION("The lumped matrix " << matrix_id << " already exists in "
                                          << this->id);
  }

  auto & ret = *matrix;
  this->lumped_matrices[matrix_id] = std::move(matrix);
  ret.resize();
  return ret;
}

/* -------------------------------------------------------------------------- */
NonLinearSolver & DOFManager::registerNonLinearSolver(
    const ID & non_linear_solver_id,
    std::unique_ptr<NonLinearSolver> & non_linear_solver) {
  NonLinearSolversMap::const_iterator it =
      this->non_linear_solvers.find(non_linear_solver_id);
  if (it != this->non_linear_solvers.end()) {
    AKANTU_EXCEPTION("The non linear solver " << non_linear_solver_id
                                              << " already exists in "
                                              << this->id);
  }

  NonLinearSolver & ret = *non_linear_solver;
  this->non_linear_solvers[non_linear_solver_id] = std::move(non_linear_solver);

  return ret;
}

/* -------------------------------------------------------------------------- */
TimeStepSolver & DOFManager::registerTimeStepSolver(
    const ID & time_step_solver_id,
    std::unique_ptr<TimeStepSolver> & time_step_solver) {
  TimeStepSolversMap::const_iterator it =
      this->time_step_solvers.find(time_step_solver_id);
  if (it != this->time_step_solvers.end()) {
    AKANTU_EXCEPTION("The non linear solver " << time_step_solver_id
                                              << " already exists in "
                                              << this->id);
  }

  TimeStepSolver & ret = *time_step_solver;
  this->time_step_solvers[time_step_solver_id] = std::move(time_step_solver);
  return ret;
}

/* -------------------------------------------------------------------------- */
SparseMatrix & DOFManager::getMatrix(const ID & id) {
  ID matrix_id = this->id + ":mtx:" + id;
  SparseMatricesMap::const_iterator it = this->matrices.find(matrix_id);
  if (it == this->matrices.end()) {
    AKANTU_SILENT_EXCEPTION("The matrix " << matrix_id << " does not exists in "
                                          << this->id);
  }

  return *(it->second);
}

/* -------------------------------------------------------------------------- */
bool DOFManager::hasMatrix(const ID & id) const {
  ID mtx_id = this->id + ":mtx:" + id;
  auto it = this->matrices.find(mtx_id);
  return it != this->matrices.end();
}

/* -------------------------------------------------------------------------- */
SolverVector & DOFManager::getLumpedMatrix(const ID & id) {
  ID matrix_id = this->id + ":lumped_mtx:" + id;
  LumpedMatricesMap::const_iterator it = this->lumped_matrices.find(matrix_id);
  if (it == this->lumped_matrices.end()) {
    AKANTU_SILENT_EXCEPTION("The lumped matrix "
                            << matrix_id << " does not exists in " << this->id);
  }

  return *(it->second);
}

/* -------------------------------------------------------------------------- */
const SolverVector & DOFManager::getLumpedMatrix(const ID & id) const {
  ID matrix_id = this->id + ":lumped_mtx:" + id;
  auto it = this->lumped_matrices.find(matrix_id);
  if (it == this->lumped_matrices.end()) {
    AKANTU_SILENT_EXCEPTION("The lumped matrix "
                            << matrix_id << " does not exists in " << this->id);
  }

  return *(it->second);
}

/* -------------------------------------------------------------------------- */
bool DOFManager::hasLumpedMatrix(const ID & id) const {
  ID mtx_id = this->id + ":lumped_mtx:" + id;
  auto it = this->lumped_matrices.find(mtx_id);
  return it != this->lumped_matrices.end();
}

/* -------------------------------------------------------------------------- */
NonLinearSolver & DOFManager::getNonLinearSolver(const ID & id) {
  ID non_linear_solver_id = this->id + ":nls:" + id;
  NonLinearSolversMap::const_iterator it =
      this->non_linear_solvers.find(non_linear_solver_id);
  if (it == this->non_linear_solvers.end()) {
    AKANTU_EXCEPTION("The non linear solver " << non_linear_solver_id
                                              << " does not exists in "
                                              << this->id);
  }

  return *(it->second);
}

/* -------------------------------------------------------------------------- */
bool DOFManager::hasNonLinearSolver(const ID & id) const {
  ID solver_id = this->id + ":nls:" + id;
  auto it = this->non_linear_solvers.find(solver_id);
  return it != this->non_linear_solvers.end();
}

/* -------------------------------------------------------------------------- */
TimeStepSolver & DOFManager::getTimeStepSolver(const ID & id) {
  ID time_step_solver_id = this->id + ":tss:" + id;
  TimeStepSolversMap::const_iterator it =
      this->time_step_solvers.find(time_step_solver_id);
  if (it == this->time_step_solvers.end()) {
    AKANTU_EXCEPTION("The non linear solver " << time_step_solver_id
                                              << " does not exists in "
                                              << this->id);
  }

  return *(it->second);
}

/* -------------------------------------------------------------------------- */
bool DOFManager::hasTimeStepSolver(const ID & solver_id) const {
  ID time_step_solver_id = this->id + ":tss:" + solver_id;
  auto it = this->time_step_solvers.find(time_step_solver_id);
  return it != this->time_step_solvers.end();
}

/* -------------------------------------------------------------------------- */
void DOFManager::savePreviousDOFs(const ID & dofs_id) {
  this->getPreviousDOFs(dofs_id).copy(this->getDOFs(dofs_id));
}

/* -------------------------------------------------------------------------- */
void DOFManager::zeroResidual() { this->residual->zero(); }

/* -------------------------------------------------------------------------- */
void DOFManager::zeroMatrix(const ID & mtx) { this->getMatrix(mtx).zero(); }

/* -------------------------------------------------------------------------- */
void DOFManager::zeroLumpedMatrix(const ID & mtx) {
  this->getLumpedMatrix(mtx).zero();
}

/* -------------------------------------------------------------------------- */
/* Mesh Events                                                                */
/* -------------------------------------------------------------------------- */
std::pair<Int, Int> DOFManager::updateNodalDOFs(const ID & dof_id,
                                                const Array<Idx> & nodes_list) {
  auto & dof_data = this->getDOFData(dof_id);
  Int nb_new_local_dofs, nb_new_pure_local;

  std::tie(nb_new_local_dofs, nb_new_pure_local) =
      countDOFsForNodes(dof_data, nodes_list.size(),
                        [&nodes_list](auto && n) { return nodes_list(n); });

  this->pure_local_system_size += nb_new_pure_local;
  this->local_system_size += nb_new_local_dofs;

  auto nb_new_global = nb_new_pure_local;
  communicator.allReduce(nb_new_global, SynchronizerOperation::_sum);

  this->system_size += nb_new_global;

  dof_data.solution.resize(local_system_size);

  updateDOFsData(dof_data, nb_new_local_dofs, nb_new_pure_local,
                 nodes_list.size(),
                 [&nodes_list](auto pos) -> UInt { return nodes_list[pos]; });

  return std::make_pair(nb_new_local_dofs, nb_new_pure_local);
}

/* -------------------------------------------------------------------------- */
void DOFManager::resizeGlobalArrays() {
  // resize all relevant arrays
  this->residual->resize();
  this->solution->resize();
  this->data_cache->resize();

  for (auto & lumped_matrix : lumped_matrices) {
    lumped_matrix.second->resize();
  }

  for (auto & matrix : matrices) {
    matrix.second->clearProfile();
  }
}

/* -------------------------------------------------------------------------- */
void DOFManager::onNodesAdded(const Array<Idx> & nodes_list,
                              const NewNodesEvent &) {
  for (auto & pair : this->dofs) {
    const auto & dof_id = pair.first;
    auto & dof_data = this->getDOFData(dof_id);
    if (dof_data.support_type != _dst_nodal) {
      continue;
    }

    const auto & group = dof_data.group_support;

    if (group == "__mesh__") {
      this->updateNodalDOFs(dof_id, nodes_list);
    } else {
      const auto & node_group =
          this->mesh->getElementGroup(group).getNodeGroup();
      Array<Idx> new_nodes_list;
      for (const auto & node : nodes_list) {
        if (node_group.find(node) != Int(-1)) {
          new_nodes_list.push_back(node);
        }
      }

      this->updateNodalDOFs(dof_id, new_nodes_list);
    }
  }

  this->resizeGlobalArrays();
}

/* -------------------------------------------------------------------------- */
void DOFManager::onMeshIsDistributed(const MeshIsDistributedEvent & /*event*/) {
  AKANTU_DEBUG_ASSERT(this->mesh != nullptr, "The `Mesh` pointer is not set.");

  // check if the distributed state of the residual and the mesh are the same.
  if (this->mesh->isDistributed() != this->residual->isDistributed()) {
    // TODO: Allow to reallocate the internals, in that case one could actually
    // react on that event.
    auto is_or_is_not = [](bool q) {
      return ((q) ? std::string("is") : std::string("is not"));
    };
    AKANTU_EXCEPTION("There is an inconsistency about the distribution state "
                     "of the `DOFManager`."
                     " It seams that the `Mesh` "
                     << is_or_is_not(this->mesh->isDistributed())
                     << " distributed, but the `DOFManager`'s residual "
                     << is_or_is_not(this->residual->isDistributed())
                     << ", which is of type "
                     << debug::demangle(typeid(this->residual).name()) << ".");
  }
}

/* -------------------------------------------------------------------------- */
class GlobalDOFInfoDataAccessor : public DataAccessor<Idx> {
public:
  GlobalDOFInfoDataAccessor(DOFManager::DOFData & dof_data,
                            DOFManager & dof_manager)
      : dof_data(dof_data), dof_manager(dof_manager) {
    for (auto && pair :
         zip(dof_data.local_equation_number, dof_data.associated_nodes)) {
      Idx node;
      Idx dof;
      std::tie(dof, node) = pair;

      dofs_per_node[node].push_back(dof);
    }
  }

  Int getNbData(const Array<Idx> & nodes,
                const SynchronizationTag & tag) const override {
    if (tag == SynchronizationTag::_ask_nodes or
        tag == SynchronizationTag::_giu_global_conn) {
      return nodes.size() * dof_data.dof->getNbComponent() * sizeof(Idx);
    }

    return 0;
  }

  void packData(CommunicationBuffer & buffer, const Array<Idx> & nodes,
                const SynchronizationTag & tag) const override {
    if (tag == SynchronizationTag::_ask_nodes or
        tag == SynchronizationTag::_giu_global_conn) {
      for (const auto & node : nodes) {
        const auto & dofs = dofs_per_node.at(node);
        for (const auto & dof : dofs) {
          buffer << dof_manager.global_equation_number(dof);
        }
      }
    }
  }

  void unpackData(CommunicationBuffer & buffer, const Array<Idx> & nodes,
                  const SynchronizationTag & tag) override {
    if (tag == SynchronizationTag::_ask_nodes or
        tag == SynchronizationTag::_giu_global_conn) {
      for (const auto & node : nodes) {
        const auto & dofs = dofs_per_node[node];
        for (const auto & dof : dofs) {
          Idx global_dof;
          buffer >> global_dof;
          AKANTU_DEBUG_ASSERT(
              (dof_manager.global_equation_number(dof) == -1 or
               dof_manager.global_equation_number(dof) == global_dof),
              "This dof already had a global_dof_id which is different from "
              "the received one. "
                  << dof_manager.global_equation_number(dof)
                  << " != " << global_dof);
          dof_manager.global_equation_number(dof) = global_dof;
          dof_manager.global_to_local_mapping[global_dof] = dof;
        }
      }
    }
  }

protected:
  std::unordered_map<Idx, std::vector<Idx>> dofs_per_node;
  DOFManager::DOFData & dof_data;
  DOFManager & dof_manager;
};

/* -------------------------------------------------------------------------- */
auto DOFManager::computeFirstDOFIDs(Int nb_new_local_dofs,
                                    Int nb_new_pure_local) {
  // determine the first local/global dof id to use
  Int offset = 0;

  this->communicator.exclusiveScan(nb_new_pure_local, offset);

  auto first_global_dof_id = this->first_global_dof_id + offset;
  auto first_local_dof_id = this->local_system_size - nb_new_local_dofs;

  offset = nb_new_pure_local;
  this->communicator.allReduce(offset);
  this->first_global_dof_id += offset;

  return std::make_pair(first_local_dof_id, first_global_dof_id);
}

/* -------------------------------------------------------------------------- */
void DOFManager::updateDOFsData(DOFData & dof_data, Int nb_new_local_dofs,
                                Int nb_new_pure_local, Int nb_node,
                                const std::function<Idx(Idx)> & getNode) {
  auto nb_local_dofs_added = nb_node * dof_data.dof->getNbComponent();

  auto first_dof_pos = dof_data.local_equation_number.size();
  dof_data.local_equation_number.reserve(dof_data.local_equation_number.size() +
                                         nb_local_dofs_added);
  dof_data.associated_nodes.reserve(dof_data.associated_nodes.size() +
                                    nb_local_dofs_added);

  this->dofs_flag.resize(this->local_system_size, NodeFlag::_normal);
  this->global_equation_number.resize(this->local_system_size, -1);

  std::unordered_map<std::pair<Idx, Idx>, Idx> masters_dofs;

  // update per dof info
  Int local_eq_num, first_global_dof_id;
  std::tie(local_eq_num, first_global_dof_id) =
      computeFirstDOFIDs(nb_new_local_dofs, nb_new_pure_local);
  for (auto d : arange(nb_local_dofs_added)) {
    auto node = getNode(d / dof_data.dof->getNbComponent());
    auto dof_flag = this->mesh->getNodeFlag(node);

    dof_data.associated_nodes.push_back(node);
    auto is_local_dof = this->mesh->isLocalOrMasterNode(node);
    auto is_periodic_slave = this->mesh->isPeriodicSlave(node);
    auto is_periodic_master = this->mesh->isPeriodicMaster(node);

    if (is_periodic_slave) {
      dof_data.local_equation_number.push_back(-1);
      continue;
    }

    // update equation numbers
    this->dofs_flag(local_eq_num) = dof_flag;
    dof_data.local_equation_number.push_back(local_eq_num);

    if (is_local_dof) {
      this->global_equation_number(local_eq_num) = first_global_dof_id;
      this->global_to_local_mapping[first_global_dof_id] = local_eq_num;
      ++first_global_dof_id;
    } else {
      this->global_equation_number(local_eq_num) = -1;
    }

    if (is_periodic_master) {
      auto node = getNode(d / dof_data.dof->getNbComponent());
      auto dof = d % dof_data.dof->getNbComponent();
      masters_dofs.insert(
          std::make_pair(std::make_pair(node, dof), local_eq_num));
    }

    ++local_eq_num;
  }

  // correct periodic slave equation numbers
  if (this->mesh->isPeriodic()) {
    auto assoc_begin = dof_data.associated_nodes.begin();
    for (auto d : arange(nb_local_dofs_added)) {
      auto node = dof_data.associated_nodes(first_dof_pos + d);
      if (not this->mesh->isPeriodicSlave(node)) {
        continue;
      }

      auto master_node = this->mesh->getPeriodicMaster(node);
      auto dof = d % dof_data.dof->getNbComponent();
      dof_data.local_equation_number(first_dof_pos + d) =
          masters_dofs[std::make_pair(master_node, dof)];
    }
  }

  // synchronize the global numbering for slaves nodes
  if (this->mesh->isDistributed()) {
    GlobalDOFInfoDataAccessor data_accessor(dof_data, *this);

    if (this->mesh->isPeriodic()) {
      mesh->getPeriodicNodeSynchronizer().synchronizeOnce(
          data_accessor, SynchronizationTag::_giu_global_conn);
    }

    auto & node_synchronizer = this->mesh->getNodeSynchronizer();
    node_synchronizer.synchronizeOnce(data_accessor,
                                      SynchronizationTag::_ask_nodes);
  }
}

/* -------------------------------------------------------------------------- */
void DOFManager::updateDOFsData(DOFData & dof_data, Int nb_new_local_dofs,
                                Int nb_new_pure_local) {
  dof_data.local_equation_number.reserve(dof_data.local_equation_number.size() +
                                         nb_new_local_dofs);

  Int first_local_dof_id, first_global_dof_id;
  std::tie(first_local_dof_id, first_global_dof_id) =
      computeFirstDOFIDs(nb_new_local_dofs, nb_new_pure_local);

  this->dofs_flag.resize(this->local_system_size, NodeFlag::_normal);
  this->global_equation_number.resize(this->local_system_size, -1);

  // update per dof info
  for (auto _ [[gnu::unused]] : arange(nb_new_local_dofs)) {
    // update equation numbers
    this->dofs_flag(first_local_dof_id) = NodeFlag::_normal;

    dof_data.local_equation_number.push_back(first_local_dof_id);

    this->global_equation_number(first_local_dof_id) = first_global_dof_id;
    this->global_to_local_mapping[first_global_dof_id] = first_local_dof_id;

    ++first_global_dof_id;
    ++first_local_dof_id;
  }
}

/* -------------------------------------------------------------------------- */
void DOFManager::onNodesRemoved(const Array<Idx> &, const Array<Idx> &,
                                const RemovedNodesEvent &) {}

/* -------------------------------------------------------------------------- */
void DOFManager::onElementsAdded(const Array<Element> & /*unused*/,
                                 const NewElementsEvent & /*unused*/) {}

/* -------------------------------------------------------------------------- */
void DOFManager::onElementsRemoved(const Array<Element> &,
                                   const ElementTypeMapArray<Idx> &,
                                   const RemovedElementsEvent &) {}

/* -------------------------------------------------------------------------- */
void DOFManager::onElementsChanged(const Array<Element> &,
                                   const Array<Element> &,
                                   const ElementTypeMapArray<Idx> &,
                                   const ChangedElementsEvent &) {}

/* -------------------------------------------------------------------------- */
void DOFManager::updateGlobalBlockedDofs() {
  this->previous_global_blocked_dofs_indexes.copy(
      this->global_blocked_dofs_indexes);
  this->global_blocked_dofs_indexes.reserve(this->local_system_size, 0);
  this->previous_global_blocked_dofs_indexes_release =
      this->global_blocked_dofs_indexes_release;

  for (auto && [dof_id, dof_data] : dofs) {
    if (not this->hasBlockedDOFs(dof_id)) {
      continue;
    }

    for (auto && [dof, is_blocked] : zip(dof_data->getLocalEquationsNumbers(),
                                         make_view(*dof_data->blocked_dofs))) {
      if (is_blocked) {
        this->global_blocked_dofs_indexes.push_back(dof);
      }
    }
  }

  std::sort(this->global_blocked_dofs_indexes.begin(),
            this->global_blocked_dofs_indexes.end());
  auto last = std::unique(this->global_blocked_dofs_indexes.begin(),
                          this->global_blocked_dofs_indexes.end());
  this->global_blocked_dofs_indexes.resize(
      last - this->global_blocked_dofs_indexes.begin());

  auto are_equal = global_blocked_dofs_indexes.size() ==
                       previous_global_blocked_dofs_indexes.size() and
                   std::equal(global_blocked_dofs_indexes.begin(),
                              global_blocked_dofs_indexes.end(),
                              previous_global_blocked_dofs_indexes.begin());

  if (this->global_blocked_dofs_indexes.size() == 0) {
    global_blocked_dofs.resize(local_system_size);
    global_blocked_dofs.set(false);
  }

  if (not are_equal) {
    ++this->global_blocked_dofs_indexes_release;
  } else {
    return;
  }

  global_blocked_dofs.resize(local_system_size);
  global_blocked_dofs.set(false);
  for (const auto & dof : global_blocked_dofs_indexes) {
    global_blocked_dofs[dof] = true;
  }
}

/* -------------------------------------------------------------------------- */
void DOFManager::applyBoundary(const ID & matrix_id) {
  auto & J = this->getMatrix(matrix_id);

  if (this->jacobian_release == J.getRelease()) {
    if (this->hasBlockedDOFsChanged()) {
      J.applyBoundary();
    }

    previous_global_blocked_dofs_indexes.copy(global_blocked_dofs_indexes);
  } else {
    J.applyBoundary();
  }

  this->jacobian_release = J.getRelease();
  this->previous_global_blocked_dofs_indexes_release =
      this->global_blocked_dofs_indexes_release;
}

/* -------------------------------------------------------------------------- */
void DOFManager::assembleMatMulVectToGlobalArray(const ID & dof_id,
                                                 const ID & A_id,
                                                 const Array<Real> & x,
                                                 SolverVector & array,
                                                 Real scale_factor) {
  auto & A = this->getMatrix(A_id);

  data_cache->resize();
  data_cache->zero();
  this->assembleToGlobalArray(dof_id, x, *data_cache, 1.);

  A.matVecMul(*data_cache, array, scale_factor, 1.);
}

/* -------------------------------------------------------------------------- */
void DOFManager::assembleMatMulVectToResidual(const ID & dof_id,
                                              const ID & A_id,
                                              const Array<Real> & x,
                                              Real scale_factor) {
  assembleMatMulVectToGlobalArray(dof_id, A_id, x, *residual, scale_factor);
}
/* -------------------------------------------------------------------------- */

} // namespace akantu
