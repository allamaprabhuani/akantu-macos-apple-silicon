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
#ifndef AKANTU_MESH_HH_
#define AKANTU_MESH_HH_

/* -------------------------------------------------------------------------- */
#include "aka_array.hh"
#include "aka_bbox.hh"
#include "aka_event_handler_manager.hh"
#include "communicator.hh"
#include "dumpable.hh"
#include "element.hh"
#include "element_class.hh"
#include "element_type_map.hh"
#include "group_manager.hh"
#include "mesh_data.hh"
#include "mesh_events.hh"
/* -------------------------------------------------------------------------- */
#include <functional>
#include <set>
#include <unordered_map>
/* -------------------------------------------------------------------------- */

namespace akantu {
class ElementSynchronizer;
class NodeSynchronizer;
class PeriodicNodeSynchronizer;
class MeshGlobalDataUpdater;
} // namespace akantu

namespace akantu {

namespace {
  DECLARE_NAMED_ARGUMENT(communicator);
  DECLARE_NAMED_ARGUMENT(edge_weight_function);
  DECLARE_NAMED_ARGUMENT(vertex_weight_function);
} // namespace

/* -------------------------------------------------------------------------- */
/* Mesh                                                                       */
/* -------------------------------------------------------------------------- */

/**
 * @class  Mesh mesh.hh
 *
 * This class contaisn the coordinates of the nodes in the Mesh.nodes
 * akant::Array, and the connectivity. The connectivity are stored in by element
 * types.
 *
 * In order to loop on all element you have to loop on all types like this :
 * @code{.cpp}
 for(auto & type : mesh.elementTypes()) {
   Int nb_element  = mesh.getNbElement(type);
   const auto & conn = mesh.getConnectivity(type);
   for(Int e = 0; e < nb_element; ++e) {
     ...
   }
 }

 or

 for_each_element(mesh, [](Element & element) {
    std::cout << element << std::endl
  });
 @endcode
*/
class Mesh : public EventHandlerManager<MeshEventHandler>,
             public GroupManager,
             public MeshData,
             public Dumpable {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
private:
  /// default constructor used for chaining, the last parameter is just to
  /// differentiate constructors
  Mesh(Int spatial_dimension, const ID & id, Communicator & communicator);

public:
  /// constructor that create nodes coordinates array
  Mesh(Int spatial_dimension, const ID & id = "mesh");

  /// mesh not distributed and not using the default communicator
  Mesh(Int spatial_dimension, Communicator & communicator,
       const ID & id = "mesh");

  /**
   * constructor that use an existing nodes coordinates
   * array, by getting the vector of coordinates
   */
  Mesh(Int spatial_dimension, const std::shared_ptr<Array<Real>> & nodes,
       const ID & id = "mesh");

  ~Mesh() override;

  Mesh(const Mesh &) = delete;
  Mesh(Mesh &&) = delete;
  Mesh & operator=(const Mesh &) = delete;
  Mesh & operator=(Mesh &&) = delete;

  /// read the mesh from a file
  void read(const std::string & filename,
            const MeshIOType & mesh_io_type = _miot_auto);
  /// write the mesh to a file
  void write(const std::string & filename,
             const MeshIOType & mesh_io_type = _miot_auto);

protected:
  void makeReady();

private:
  /// initialize the connectivity to NULL and other stuff
  void init();

  /// function that computes the bounding box (fills xmin, xmax)
  void computeBoundingBox();

  /* ------------------------------------------------------------------------ */
  /* Distributed memory methods and accessors                                 */
  /* ------------------------------------------------------------------------ */
public:
protected:
  /// patitionate the mesh among the processors involved in their computation
  virtual void distributeImpl(
      Communicator & communicator,
      const std::function<Int(const Element &, const Element &)> &
          edge_weight_function,
      const std::function<Int(const Element &)> & vertex_weight_function);

public:
  /// with the arguments to pass to the partitionner
  template <typename... pack>
  std::enable_if_t<are_named_argument<pack...>::value>
  distribute(pack &&... _pack) {
    distributeImpl(
        OPTIONAL_NAMED_ARG(communicator, Communicator::getStaticCommunicator()),
        OPTIONAL_NAMED_ARG(edge_weight_function,
                           [](auto &&, auto &&) { return 1; }),
        OPTIONAL_NAMED_ARG(vertex_weight_function, [](auto &&) { return 1; }));
  }

  /// defines is the mesh is distributed or not
  inline bool isDistributed() const { return this->is_distributed; }

  /* ------------------------------------------------------------------------ */
  /* Periodicity methods and accessors                                        */
  /* ------------------------------------------------------------------------ */
public:
  /// set the periodicity in a given direction
  void makePeriodic(const SpatialDirection & direction);
  void makePeriodic(const SpatialDirection & direction, const ID & list_1,
                    const ID & list_2);

protected:
  void makePeriodic(const SpatialDirection & direction,
                    const Array<Idx> & list_1, const Array<Idx> & list_2);

  /// Removes the face that the mesh is periodic
  void wipePeriodicInfo();

  inline void addPeriodicSlave(Idx slave, Idx master);

  template <typename T>
  void synchronizePeriodicSlaveDataWithMaster(Array<T> & data);

  // update the periodic synchronizer (creates it if it does not exists)
  void updatePeriodicSynchronizer();

public:
  /// defines if the mesh is periodic or not
  inline bool isPeriodic() const { return this->is_periodic; }

  inline bool isPeriodic(const SpatialDirection & /*direction*/) const {
    return this->is_periodic;
  }

  class PeriodicSlaves;

  /// get the master node for a given slave nodes, except if node not a slave
  inline Idx getPeriodicMaster(Idx slave) const;

  /// get an iterable list of slaves for a given master node
  inline decltype(auto) getPeriodicSlaves(Idx master) const;

  /* ------------------------------------------------------------------------ */
  /* General Methods                                                          */
  /* ------------------------------------------------------------------------ */
public:
  /// function to print the containt of the class
  void printself(std::ostream & stream, int indent = 0) const override;

  /// extract coordinates of nodes from an element
  template <typename T, class Derived1, class Derived2,
            std::enable_if_t<aka::is_vector_v<Derived2>> * = nullptr>
  inline void extractNodalValuesFromElement(
      const Array<T> & nodal_values,
      Eigen::MatrixBase<Derived1> & elemental_values,
      const Eigen::MatrixBase<Derived2> & connectivity) const;

  /// extract coordinates of nodes from an element
  template <typename T>
  inline decltype(auto)
  extractNodalValuesFromElement(const Array<T> & nodal_values,
                                const Element & element) const;

  /// add a Array of connectivity for the given ElementType and GhostType .
  inline void addConnectivityType(ElementType type,
                                  GhostType ghost_type = _not_ghost);

  /* ------------------------------------------------------------------------ */
  template <class Event> void sendEvent(Event & event);

  /// prepare the  event to remove the elements listed
  void eraseElements(const Array<Element> & elements);

  /* ------------------------------------------------------------------------ */
  template <typename T>
  inline void removeNodesFromArray(Array<T> & vect,
                                   const Array<Int> & new_numbering);

  /// init facets' mesh
  Mesh & initMeshFacets(const ID & id = "mesh_facets");

  /// define parent mesh
  void defineMeshParent(const Mesh & mesh);

  /// get global connectivity array
  void getGlobalConnectivity(ElementTypeMapArray<Int> & global_connectivity);

public:
  void getAssociatedElements(const Array<Int> & node_list,
                             Array<Element> & elements);

  void getAssociatedElements(const Idx & node, Array<Element> & elements) const;

  inline decltype(auto) getAssociatedElements(const Idx & node) const;

public:
  /// fills the nodes_to_elements for given dimension elements
  void fillNodesToElements(Int dimension = _all_dimensions);

private:
  /// update the global ids, nodes type, ...
  std::tuple<Int, Int> updateGlobalData(NewNodesEvent & nodes_event,
                                        NewElementsEvent & elements_event);

  void registerGlobalDataUpdater(
      std::unique_ptr<MeshGlobalDataUpdater> && global_data_updater);
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get the id of the mesh
  AKANTU_GET_MACRO(ID, id, const ID &);

  /// get the spatial dimension of the mesh = number of component of the
  /// coordinates
  AKANTU_GET_MACRO(SpatialDimension, spatial_dimension, Int);

  /// get the nodes Array aka coordinates
  AKANTU_GET_MACRO(Nodes, *nodes, const Array<Real> &);
  AKANTU_GET_MACRO_NOT_CONST(Nodes, *nodes, Array<Real> &);

  /// get the number of nodes
  auto getNbNodes() const { return nodes->size(); }

  /// get the Array of global ids of the nodes (only used in parallel)
  AKANTU_GET_MACRO_AUTO(GlobalNodesIds, *nodes_global_ids);
  // AKANTU_GET_MACRO_NOT_CONST(GlobalNodesIds, *nodes_global_ids, Array<UInt>
  // &);

  /// get the global id of a node
  inline auto getNodeGlobalId(Idx local_id) const;

  /// get the global id of a node
  inline auto getNodeLocalId(Idx global_id) const;

  /// get the global number of nodes
  inline auto getNbGlobalNodes() const;

  /// get the nodes type Array
  AKANTU_GET_MACRO(NodesFlags, *nodes_flags, const Array<NodeFlag> &);

protected:
  AKANTU_GET_MACRO_NOT_CONST(NodesFlags, *nodes_flags, Array<NodeFlag> &);

public:
  inline NodeFlag getNodeFlag(Idx local_id) const;
  inline auto getNodePrank(Idx local_id) const;

  /// say if a node is a pure ghost node
  inline bool isPureGhostNode(Idx n) const;

  /// say if a node is pur local or master node
  inline bool isLocalOrMasterNode(Idx n) const;

  inline bool isLocalNode(Idx n) const;
  inline bool isMasterNode(Idx n) const;
  inline bool isSlaveNode(Idx n) const;

  inline bool isPeriodicSlave(Idx n) const;
  inline bool isPeriodicMaster(Idx n) const;

  const Vector<Real> & getLowerBounds() const { return bbox.getLowerBounds(); }
  const Vector<Real> & getUpperBounds() const { return bbox.getUpperBounds(); }
  AKANTU_GET_MACRO(BBox, bbox, const BBox &);

  const Vector<Real> & getLocalLowerBounds() const {
    return bbox_local.getLowerBounds();
  }
  const Vector<Real> & getLocalUpperBounds() const {
    return bbox_local.getUpperBounds();
  }
  AKANTU_GET_MACRO(LocalBBox, bbox_local, const BBox &);

  /// get the connectivity Array for a given type
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Connectivity, connectivities, Idx);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(Connectivity, connectivities, Idx);
  AKANTU_GET_MACRO(Connectivities, connectivities,
                   const ElementTypeMapArray<Idx> &);

  /// get the number of element of a type in the mesh
  inline auto getNbElement(ElementType type,
                           GhostType ghost_type = _not_ghost) const;

  /// get the number of element for a given ghost_type and a given dimension
  inline auto getNbElement(Int spatial_dimension = _all_dimensions,
                           GhostType ghost_type = _not_ghost,
                           ElementKind kind = _ek_not_defined) const;

  /// compute the barycenter of a given element
  template <class D, std::enable_if_t<aka::is_vector_v<D>> * = nullptr>
  inline void getBarycenter(const Element & element,
                            const Eigen::MatrixBase<D> & barycenter) const;

  inline Vector<Real> getBarycenter(const Element & element) const;

  void getBarycenters(Array<Real> & barycenter, ElementType type,
                      GhostType ghost_type) const;

  /// get the element connected to a subelement (element of lower dimension)
  decltype(auto) getElementToSubelement() const;

  /// get the element connected to a subelement
  const auto & getElementToSubelement(ElementType el_type,
                                      GhostType ghost_type = _not_ghost) const;

  /// get the elements connected to a subelement
  decltype(auto) getElementToSubelement(const Element & element) const;

  /// get the subelement (element of lower dimension) connected to a element
  decltype(auto) getSubelementToElement() const;

  /// get the subelement connected to an element
  const auto & getSubelementToElement(ElementType el_type,
                                      GhostType ghost_type = _not_ghost) const;

  /// get the subelement (element of lower dimension) connected to a element
  decltype(auto) getSubelementToElement(const Element & element) const;

  /// get connectivity of a given element
  inline decltype(auto) getConnectivity(const Element & element) const;
  inline decltype(auto)
  getConnectivityWithPeriodicity(const Element & element) const;

protected:
  /// get the element connected to a subelement (element of lower dimension)
  auto & getElementToSubelementNC();
  auto & getSubelementToElementNC();
  inline auto & getElementToSubelementNC(const Element & element);
  inline decltype(auto) getSubelementToElementNC(const Element & element);
  /// get the element connected to a subelement
  auto & getElementToSubelementNC(ElementType el_type,
                                  GhostType ghost_type = _not_ghost);
  /// get the subelement connected to an element
  auto & getSubelementToElementNC(ElementType el_type,
                                  GhostType ghost_type = _not_ghost);

  inline decltype(auto) getConnectivityNC(const Element & element);

public:
  /// get a name field associated to the mesh
  template <typename T>
  inline decltype(auto) getData(const ID & data_name, ElementType el_type,
                                GhostType ghost_type = _not_ghost) const;

  /// get a name field associated to the mesh
  template <typename T>
  inline decltype(auto) getData(const ID & data_name, ElementType el_type,
                                GhostType ghost_type = _not_ghost);

  /// get a name field associated to the mesh
  template <typename T>
  inline decltype(auto) getData(const ID & data_name) const;

  /// get a name field associated to the mesh
  template <typename T> inline decltype(auto) getData(const ID & data_name);

  template <typename T>
  inline decltype(auto) getData(const ID & data_name, Element element) const;

  template <typename T>
  inline decltype(auto) getData(const ID & data_name, Element element);

  template <typename T>
  auto getNbDataPerElem(ElementTypeMapArray<T> & array) -> ElementTypeMap<Int>;

  template <typename T>
  std::shared_ptr<dumpers::Field>
  createFieldFromAttachedData(const std::string & field_id,
                              const std::string & group_name,
                              ElementKind element_kind);

  /// templated getter returning the pointer to data in MeshData (modifiable)
  template <typename T>
  inline decltype(auto)
  getDataPointer(const std::string & data_name, ElementType el_type,
                 GhostType ghost_type = _not_ghost, Int nb_component = 1,
                 bool size_to_nb_element = true,
                 bool resize_with_parent = false);

  template <typename T>
  inline decltype(auto)
  getDataPointer(const ID & data_name, ElementType el_type,
                 GhostType ghost_type, Int nb_component,
                 bool size_to_nb_element, bool resize_with_parent,
                 const T & defaul_);

  /// Facets mesh accessor
  inline auto getMeshFacets() const -> const Mesh &;
  inline auto getMeshFacets() -> Mesh &;

  inline auto hasMeshFacets() const { return mesh_facets != nullptr; }

  /// Parent mesh accessor
  inline auto getMeshParent() const -> const Mesh &;

  inline auto isMeshFacets() const { return this->is_mesh_facets; }

  /// return the dumper from a group and and a dumper name
  auto getGroupDumper(const std::string & dumper_name,
                      const std::string & group_name) -> DumperIOHelper &;

  /* ------------------------------------------------------------------------ */
  /* Wrappers on ElementClass functions                                       */
  /* ------------------------------------------------------------------------ */
public:
  /// get the number of nodes per element for a given element type
  static inline constexpr auto getNbNodesPerElement(ElementType type) -> Int;

  /// get the number of nodes per element for a given element type considered as
  /// a first order element
  static inline constexpr auto getP1ElementType(ElementType type)
      -> ElementType;

  /// get the kind of the element type
  static inline constexpr auto getKind(ElementType type) -> ElementKind;

  /// get spatial dimension of a type of element
  static inline constexpr auto getSpatialDimension(ElementType type) -> Int;

  /// get the natural space dimension of a type of element
  static inline constexpr auto getNaturalSpaceDimension(ElementType type)
      -> Int;

  /// get number of facets of a given element type
  static inline constexpr auto getNbFacetsPerElement(ElementType type) -> Int;

  /// get number of facets of a given element type
  static inline constexpr auto getNbFacetsPerElement(ElementType type, Idx t)
      -> Int;

  /// get local connectivity of a facet for a given facet type
  static inline decltype(auto) getFacetLocalConnectivity(ElementType type,
                                                         Idx t = 0);

  /// get connectivity of facets for a given element
  inline auto getFacetConnectivity(const Element & element, Idx t = 0) const
      -> Matrix<Idx>;

  /// get the number of type of the surface element associated to a given
  /// element type
  static inline constexpr auto getNbFacetTypes(ElementType type, Idx t = 0)
      -> Int;

  /// get the type of the surface element associated to a given element
  static inline constexpr auto getFacetType(ElementType type, Idx t = 0)
      -> ElementType;

  /// get all the type of the surface element associated to a given element
  static inline decltype(auto) getAllFacetTypes(ElementType type);

  /// get the number of nodes in the given element list
  static inline auto getNbNodesPerElementList(const Array<Element> & elements)
      -> Int;

  /* ------------------------------------------------------------------------ */
  /* Element type Iterator                                                    */
  /* ------------------------------------------------------------------------ */
  using ElementTypesIteratorHelper =
      ElementTypeMapArray<Idx, ElementType>::ElementTypesIteratorHelper;

  template <typename... pack>
  auto elementTypes(pack &&... _pack) const -> ElementTypesIteratorHelper;

  AKANTU_GET_MACRO_DEREF_PTR(ElementSynchronizer, element_synchronizer);
  AKANTU_GET_MACRO_DEREF_PTR_NOT_CONST(ElementSynchronizer,
                                       element_synchronizer);
  AKANTU_GET_MACRO_DEREF_PTR(NodeSynchronizer, node_synchronizer);
  AKANTU_GET_MACRO_DEREF_PTR_NOT_CONST(NodeSynchronizer, node_synchronizer);
  AKANTU_GET_MACRO_DEREF_PTR(PeriodicNodeSynchronizer,
                             periodic_node_synchronizer);
  AKANTU_GET_MACRO_DEREF_PTR_NOT_CONST(PeriodicNodeSynchronizer,
                                       periodic_node_synchronizer);

  // AKANTU_GET_MACRO_NOT_CONST(Communicator, *communicator, StaticCommunicator
  // &);
  AKANTU_GET_MACRO_DEREF_PTR(Communicator, communicator);
  AKANTU_GET_MACRO_DEREF_PTR_NOT_CONST(Communicator, communicator);
  AKANTU_GET_MACRO_AUTO(PeriodicMasterSlaves, periodic_master_slave);

  /* ------------------------------------------------------------------------ */
  /* Private methods for friends                                              */
  /* ------------------------------------------------------------------------ */
private:
  friend class MeshAccessor;
  friend class MeshUtils;

  AKANTU_GET_MACRO_DEREF_PTR_NOT_CONST(NodesPointer, nodes);

  /// get a pointer to the nodes_global_ids Array<UInt> and create it if
  /// necessary
  inline auto getNodesGlobalIdsPointer() -> Array<Idx> &;

  /// get a pointer to the nodes_type Array<Int> and create it if necessary
  inline auto getNodesFlagsPointer() -> Array<NodeFlag> &;

  /// get a pointer to the connectivity Array for the given type and create it
  /// if necessary
  inline auto getConnectivityPointer(ElementType type,
                                     GhostType ghost_type = _not_ghost)
      -> Array<Idx> &;

  /// get the ghost element counter
  inline auto getGhostsCounters(ElementType type, GhostType ghost_type = _ghost)
      -> Array<Idx> & {
    AKANTU_DEBUG_ASSERT(ghost_type != _not_ghost,
                        "No ghost counter for _not_ghost elements");
    return ghosts_counters(type, ghost_type);
  }

  /// get a pointer to the element_to_subelement Array for the given type and
  /// create it if necessary
  inline decltype(auto)
  getElementToSubelementPointer(ElementType type,
                                GhostType ghost_type = _not_ghost);

  /// get a pointer to the subelement_to_element Array for the given type and
  /// create it if necessary
  inline decltype(auto)
  getSubelementToElementPointer(ElementType type,
                                GhostType ghost_type = _not_ghost);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  ID id;

  /// array of the nodes coordinates
  std::shared_ptr<Array<Real>> nodes;

  /// global node ids
  std::shared_ptr<Array<Idx>> nodes_global_ids;

  /// node flags (shared/periodic/...)
  std::shared_ptr<Array<NodeFlag>> nodes_flags;

  /// processor handling the node when not local or master
  std::unordered_map<Idx, Int> nodes_prank;

  /// global number of nodes;
  Int nb_global_nodes{0};

  /// all class of elements present in this mesh (for heterogenous meshes)
  ElementTypeMapArray<Idx> connectivities;

  /// count the references on ghost elements
  ElementTypeMapArray<Idx> ghosts_counters;

  /// the spatial dimension of this mesh
  Idx spatial_dimension{0};

  /// size covered by the mesh on each direction
  Vector<Real> size;

  /// global bounding box
  BBox bbox;

  /// local bounding box
  BBox bbox_local;

  /// Extra data loaded from the mesh file
  // MeshData mesh_data;

  /// facets' mesh
  std::unique_ptr<Mesh> mesh_facets;

  /// parent mesh (this is set for mesh_facets meshes)
  const Mesh * mesh_parent{nullptr};

  /// defines if current mesh is mesh_facets or not
  bool is_mesh_facets{false};

  /// defines if the mesh is centralized or distributed
  bool is_distributed{false};

  /// defines if the mesh is periodic
  bool is_periodic{false};

  /// Communicator on which mesh is distributed
  Communicator * communicator;

  /// Element synchronizer
  std::unique_ptr<ElementSynchronizer> element_synchronizer;

  /// Node synchronizer
  std::unique_ptr<NodeSynchronizer> node_synchronizer;

  /// Node synchronizer for periodic nodes
  std::unique_ptr<PeriodicNodeSynchronizer> periodic_node_synchronizer;

  using NodesToElements = std::vector<std::unique_ptr<std::set<Element>>>;

  /// class to update global data using external knowledge
  std::unique_ptr<MeshGlobalDataUpdater> global_data_updater;

  /// This info is stored to simplify the dynamic changes
  NodesToElements nodes_to_elements;

  /// periodicity local info
  std::unordered_map<Idx, Idx> periodic_slave_master;
  std::unordered_multimap<Idx, Idx> periodic_master_slave;
};

/// standard output stream operator
inline std::ostream & operator<<(std::ostream & stream, const Mesh & _this) {
  _this.printself(stream);
  return stream;
}

/* -------------------------------------------------------------------------- */
inline constexpr auto Mesh::getNbNodesPerElement(ElementType type) -> Int {
  return tuple_dispatch_with_default<AllElementTypes>(
      [](auto && enum_type) {
        constexpr ElementType type = aka::decay_v<decltype(enum_type)>;
        return ElementClass<type>::getNbNodesPerElement();
      },
      [](auto && /*type*/) { return 0; }, type);
}

/* -------------------------------------------------------------------------- */
inline auto Mesh::getNbElement(ElementType type, GhostType ghost_type) const {
  try {
    const auto & conn = connectivities(type, ghost_type);
    return conn.size();
  } catch (...) {
    return 0;
  }
}

/* -------------------------------------------------------------------------- */
inline auto Mesh::getNbElement(const Int spatial_dimension,
                               GhostType ghost_type, ElementKind kind) const {
  AKANTU_DEBUG_ASSERT(spatial_dimension <= 3 || spatial_dimension == Int(-1),
                      "spatial_dimension is " << spatial_dimension
                                              << " and is greater than 3 !");
  Int nb_element = 0;

  for (auto type : elementTypes(spatial_dimension, ghost_type, kind)) {
    nb_element += getNbElement(type, ghost_type);
  }

  return nb_element;
}
/* -------------------------------------------------------------------------- */

} // namespace akantu

/* -------------------------------------------------------------------------- */
/* Inline functions                                                           */
/* -------------------------------------------------------------------------- */
#include "element_type_map_tmpl.hh"
#include "mesh_inline_impl.hh"

/* -------------------------------------------------------------------------- */
#include "element_group.hh"
#include "node_group.hh"
/* -------------------------------------------------------------------------- */

#endif /* AKANTU_MESH_HH_ */
