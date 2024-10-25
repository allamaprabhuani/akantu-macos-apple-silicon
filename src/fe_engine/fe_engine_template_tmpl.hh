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
#include "aka_common.hh"
#include "dof_manager.hh"
#include "fe_engine_template.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::FEEngineTemplate(
    Mesh & mesh, Int spatial_dimension, const ID & id, bool do_not_precompute)
    : FEEngine(mesh, spatial_dimension, id),
      integrator(mesh, spatial_dimension, id),
      shape_functions(mesh, spatial_dimension, id) {
  if (not do_not_precompute) {
    initShapeFunctions(_not_ghost);
    initShapeFunctions(_ghost);
  }
}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::~FEEngineTemplate() =
    default;

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::
    gradientOnIntegrationPoints(const Array<Real> & u, Array<Real> & nablauq,
                                Int nb_degree_of_freedom, ElementType type,
                                GhostType ghost_type,
                                const Array<Int> & filter_elements) const {
  AKANTU_DEBUG_IN();

  auto nb_element = mesh.getNbElement(type, ghost_type);
  if (filter_elements != empty_filter) {
    nb_element = filter_elements.size();
  }
  auto nb_points =
      shape_functions.getIntegrationPoints(type, ghost_type).cols();

#ifndef AKANTU_NDEBUG

  auto element_dimension = mesh.getSpatialDimension(type);

  AKANTU_DEBUG_ASSERT(u.size() == mesh.getNbNodes(),
                      "The vector u(" << u.getID()
                                      << ") has not the good size.");
  AKANTU_DEBUG_ASSERT(u.getNbComponent() == nb_degree_of_freedom,
                      "The vector u("
                          << u.getID()
                          << ") has not the good number of component.");

  AKANTU_DEBUG_ASSERT(
      nablauq.getNbComponent() == nb_degree_of_freedom * element_dimension,
      "The vector nablauq(" << nablauq.getID()
                            << ") has not the good number of component.");
#endif

  nablauq.resize(nb_element * nb_points);

  auto call = [&](auto && integral_type) {
    constexpr ElementType type = std::decay_t<decltype(integral_type)>::value;
    if (element_dimension == ElementClass<type>::getSpatialDimension())
      shape_functions.template gradientOnIntegrationPoints<type>(
          u, nablauq, nb_degree_of_freedom, ghost_type, filter_elements);
  };

  tuple_dispatch<ElementTypes_t<kind>>(call, type);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::initShapeFunctions(
    GhostType ghost_type) {
  initShapeFunctions(mesh.getNodes(), ghost_type);
}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::initShapeFunctions(
    const Array<Real> & nodes, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  for (const auto & type :
       mesh.elementTypes(element_dimension, ghost_type, kind)) {
    integrator.initIntegrator(nodes, type, ghost_type);
    const auto & control_points = getIntegrationPoints(type, ghost_type);
    shape_functions.initShapeFunctions(nodes, control_points, type, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::integrate(
    const ElementTypeMapArray<Real> & f, ElementTypeMapArray<Real> & intf,
    const ElementTypeMapArray<Idx> * filter_elements) const {
  AKANTU_DEBUG_IN();

  const Array<Idx> * filter = nullptr;

  for (auto ghost_type : ghost_types) {
    auto && types = f.elementTypes(_all_dimensions, ghost_type, kind);
    for (const auto & type : types) {

      if (filter_elements != nullptr) {
        filter = &((*filter_elements)(type, ghost_type));
      } else {
        filter = &empty_filter;
      }

      const Array<Real> & fq = f(type, ghost_type);
      Array<Real> & intf_el = intf(type, ghost_type);

      integrate(fq, intf_el, fq.getNbComponent(), type, ghost_type, *filter);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::integrate(
    const Array<Real> & f, Array<Real> & intf, Int nb_degree_of_freedom,
    ElementType type, GhostType ghost_type,
    const Array<Idx> & filter_elements) const {

  auto nb_element = mesh.getNbElement(type, ghost_type);
  if (filter_elements != empty_filter) {
    nb_element = filter_elements.size();
  }
#ifndef AKANTU_NDEBUG

  auto nb_quadrature_points = getNbIntegrationPoints(type);

  AKANTU_DEBUG_ASSERT(f.size() == nb_element * nb_quadrature_points,
                      "The vector f(" << f.getID() << " size " << f.size()
                                      << ") has not the good size ("
                                      << nb_element << ").");
  AKANTU_DEBUG_ASSERT(f.getNbComponent() == nb_degree_of_freedom,
                      "The vector f("
                          << f.getID()
                          << ") has not the good number of component.");
  AKANTU_DEBUG_ASSERT(intf.getNbComponent() == nb_degree_of_freedom,
                      "The vector intf("
                          << intf.getID()
                          << ") has not the good number of component.");
#endif

  intf.resize(nb_element);

  auto && call = [&](auto && enum_type) {
    constexpr ElementType type = std ::decay_t<decltype(enum_type)>::value;
    integrator.template integrate<ElementType(type)>(
        f, intf, nb_degree_of_freedom, ghost_type, filter_elements);
  };
  tuple_dispatch<ElementTypes_t<kind>>(call, type);
}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
Real FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::integrate(
    const Array<Real> & f, ElementType type, GhostType ghost_type,
    const Array<Int> & filter_elements) const {
  AKANTU_DEBUG_IN();

#ifndef AKANTU_NDEBUG
  auto nb_element = mesh.getNbElement(type, ghost_type);
  if (filter_elements != empty_filter) {
    nb_element = filter_elements.size();
  }

  auto nb_quadrature_points = getNbIntegrationPoints(type, ghost_type);

  AKANTU_DEBUG_ASSERT(
      f.size() == nb_element * nb_quadrature_points,
      "The vector f(" << f.getID() << ") has not the good size. (" << f.size()
                      << "!=" << nb_quadrature_points * nb_element << ")");
  AKANTU_DEBUG_ASSERT(f.getNbComponent() == 1,
                      "The vector f("
                          << f.getID()
                          << ") has not the good number of component.");
#endif

  auto && call = [&](auto && enum_type) {
    constexpr ElementType type = std ::decay_t<decltype(enum_type)>::value;
    return integrator.template integrate<type>(f, ghost_type, filter_elements);
  };
  return tuple_dispatch<ElementTypes_t<kind>>(call, type);
}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
Real FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::integrate(
    const Ref<const VectorXr> f, ElementType type, Int index,
    GhostType ghost_type) const {
  auto && call = [&](auto && enum_type) {
    constexpr ElementType type = std ::decay_t<decltype(enum_type)>::value;
    return integrator.template integrate<type>(f, index, ghost_type);
  };
  return tuple_dispatch<ElementTypes_t<kind>>(call, type);
}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::
    interpolateOnIntegrationPoints(const Array<Real> & u, Array<Real> & uq,
                                   Int nb_degree_of_freedom, ElementType type,
                                   GhostType ghost_type,
                                   const Array<Int> & filter_elements) const {
  AKANTU_DEBUG_IN();

  auto nb_points =
      shape_functions.getIntegrationPoints(type, ghost_type).cols();
  auto nb_element = mesh.getNbElement(type, ghost_type);
  if (filter_elements != empty_filter) {
    nb_element = filter_elements.size();
  }

  AKANTU_DEBUG_ASSERT(u.size() == mesh.getNbNodes(),
                      "The vector u(" << u.getID()
                                      << ") has not the good size.");
  AKANTU_DEBUG_ASSERT(u.getNbComponent() == nb_degree_of_freedom,
                      "The vector u("
                          << u.getID()
                          << ") has not the good number of component.");
  AKANTU_DEBUG_ASSERT(uq.getNbComponent() == nb_degree_of_freedom,
                      "The vector uq("
                          << uq.getID()
                          << ") has not the good number of component.");

  uq.resize(nb_element * nb_points);

  auto && call = [&](auto && enum_type) {
    constexpr ElementType type = std ::decay_t<decltype(enum_type)>::value;
    shape_functions.template interpolateOnIntegrationPoints<type>(
        u, uq, nb_degree_of_freedom, ghost_type, filter_elements);
  };
  tuple_dispatch<ElementTypes_t<kind>>(call, type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::
    interpolateOnIntegrationPoints(
        const Array<Real> & u, ElementTypeMapArray<Real> & uq,
        const ElementTypeMapArray<Idx> * filter_elements) const {
  AKANTU_DEBUG_IN();

  const Array<Idx> * filter = nullptr;

  for (auto ghost_type : ghost_types) {
    auto && types = uq.elementTypes(_all_dimensions, ghost_type, kind);
    for (const auto & type : types) {
      auto nb_quad_per_element = getNbIntegrationPoints(type, ghost_type);

      Int nb_element = 0;

      if (filter_elements != nullptr) {
        filter = &((*filter_elements)(type, ghost_type));
        nb_element = filter->size();
      } else {
        filter = &empty_filter;
        nb_element = mesh.getNbElement(type, ghost_type);
      }

      auto nb_tot_quad = nb_quad_per_element * nb_element;

      Array<Real> & quad = uq(type, ghost_type);
      quad.resize(nb_tot_quad);

      interpolateOnIntegrationPoints(u, quad, quad.getNbComponent(), type,
                                     ghost_type, *filter);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
inline void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::computeBtD(
    const Array<Real> & Ds, Array<Real> & BtDs, ElementType type,
    GhostType ghost_type, const Array<Idx> & filter_elements) const {

  auto && call = [&](auto && enum_type) {
    constexpr ElementType type = std ::decay_t<decltype(enum_type)>::value;
    shape_functions.template computeBtD<type>(Ds, BtDs, ghost_type,
                                              filter_elements);
  };
  tuple_dispatch<ElementTypes_t<kind>>(call, type);
}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
inline void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::computeBtDB(
    const Array<Real> & Ds, Array<Real> & BtDBs, Int order_d, ElementType type,
    GhostType ghost_type, const Array<Idx> & filter_elements) const {
  auto && call = [&](auto && enum_type) {
    constexpr ElementType type = std ::decay_t<decltype(enum_type)>::value;
    shape_functions.template computeBtDB<type>(Ds, BtDBs, order_d, ghost_type,
                                               filter_elements);
  };
  tuple_dispatch<ElementTypes_t<kind>>(call, type);
}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
inline void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::computeNtbN(
    const Array<Real> & bs, Array<Real> & NtbNs, ElementType type,
    GhostType ghost_type, const Array<Idx> & filter_elements) const {
  auto && call = [&](auto && enum_type) {
    constexpr ElementType type = std ::decay_t<decltype(enum_type)>::value;
    shape_functions.template computeNtbN<type>(bs, NtbNs, ghost_type,
                                               filter_elements);
  };
  tuple_dispatch<ElementTypes_t<kind>>(call, type);
}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
inline void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::computeNtb(
    const Array<Real> & bs, Array<Real> & Ntbs, ElementType type,
    GhostType ghost_type, const Array<Idx> & filter_elements) const {
  auto && call = [&](auto && enum_type) {
    constexpr ElementType type = std ::decay_t<decltype(enum_type)>::value;
    shape_functions.template computeNtb<type>(bs, Ntbs, ghost_type,
                                              filter_elements);
  };
  tuple_dispatch<ElementTypes_t<kind>>(call, type);
}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
inline void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::
    computeIntegrationPointsCoordinates(
        ElementTypeMapArray<Real> & quadrature_points_coordinates,
        const ElementTypeMapArray<Idx> * filter_elements) const {

  const Array<Real> & nodes_coordinates = mesh.getNodes();

  interpolateOnIntegrationPoints(
      nodes_coordinates, quadrature_points_coordinates, filter_elements);
}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
inline void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::
    computeIntegrationPointsCoordinates(
        Array<Real> & quadrature_points_coordinates, ElementType type,
        GhostType ghost_type, const Array<Idx> & filter_elements) const {
  const Array<Real> & nodes_coordinates = mesh.getNodes();

  auto spatial_dimension = mesh.getSpatialDimension();

  interpolateOnIntegrationPoints(
      nodes_coordinates, quadrature_points_coordinates, spatial_dimension, type,
      ghost_type, filter_elements);
}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
inline void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::
    initElementalFieldInterpolationFromIntegrationPoints(
        const ElementTypeMapArray<Real> & interpolation_points_coordinates,
        ElementTypeMapArray<Real> & interpolation_points_coordinates_matrices,
        ElementTypeMapArray<Real> & quad_points_coordinates_inv_matrices,
        const ElementTypeMapArray<Idx> * element_filter) const {
  AKANTU_DEBUG_IN();

  auto spatial_dimension = this->mesh.getSpatialDimension();

  ElementTypeMapArray<Real> quadrature_points_coordinates(
      "quadrature_points_coordinates_for_interpolation", getID());

  quadrature_points_coordinates.initialize(*this,
                                           _nb_component = spatial_dimension);

  computeIntegrationPointsCoordinates(quadrature_points_coordinates,
                                      element_filter);
  shape_functions.initElementalFieldInterpolationFromIntegrationPoints(
      interpolation_points_coordinates,
      interpolation_points_coordinates_matrices,
      quad_points_coordinates_inv_matrices, quadrature_points_coordinates,
      element_filter);
}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
inline void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::
    interpolateElementalFieldFromIntegrationPoints(
        const ElementTypeMapArray<Real> & field,
        const ElementTypeMapArray<Real> & interpolation_points_coordinates,
        ElementTypeMapArray<Real> & result, const GhostType ghost_type,
        const ElementTypeMapArray<Idx> * element_filter) const {

  ElementTypeMapArray<Real> interpolation_points_coordinates_matrices(
      "interpolation_points_coordinates_matrices", id);
  ElementTypeMapArray<Real> quad_points_coordinates_inv_matrices(
      "quad_points_coordinates_inv_matrices", id);

  initElementalFieldInterpolationFromIntegrationPoints(
      interpolation_points_coordinates,
      interpolation_points_coordinates_matrices,
      quad_points_coordinates_inv_matrices, element_filter);

  interpolateElementalFieldFromIntegrationPoints(
      field, interpolation_points_coordinates_matrices,
      quad_points_coordinates_inv_matrices, result, ghost_type, element_filter);
}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
inline void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::
    interpolateElementalFieldFromIntegrationPoints(
        const ElementTypeMapArray<Real> & field,
        const ElementTypeMapArray<Real> &
            interpolation_points_coordinates_matrices,
        const ElementTypeMapArray<Real> & quad_points_coordinates_inv_matrices,
        ElementTypeMapArray<Real> & result, const GhostType ghost_type,
        const ElementTypeMapArray<Idx> * element_filter) const {
  shape_functions.interpolateElementalFieldFromIntegrationPoints(
      field, interpolation_points_coordinates_matrices,
      quad_points_coordinates_inv_matrices, result, ghost_type, element_filter);
}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
template <ElementKind kind_, typename D1, typename D2, typename D3,
          std::enable_if_t<aka::are_vectors<D1, D3>::value and
                           kind_ == _ek_regular> *>
inline void
FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::interpolateImpl(
    const Eigen::MatrixBase<D1> & real_coords,
    const Eigen::MatrixBase<D2> & nodal_values,
    Eigen::MatrixBase<D3> & interpolated, const Element & element) const {
  /// add sfinea to call only on _ek_regular
  auto && call = [&](auto && enum_type) {
    constexpr ElementType type = std::decay_t<decltype(enum_type)>::value;
    shape_functions.template interpolate<type>(real_coords, element.element,
                                               nodal_values, interpolated,
                                               element.ghost_type);
  };
  tuple_dispatch<ElementTypes_t<kind>>(call, element.type);
}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
inline void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::interpolate(
    const Ref<const VectorXr> real_coords,
    const Ref<const MatrixXr> nodal_values, Ref<VectorXr> interpolated,
    const Element & element) const {
  interpolateImpl(real_coords, nodal_values, interpolated, element);
}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::
    computeNormalsOnIntegrationPoints(GhostType ghost_type) {
  computeNormalsOnIntegrationPoints(mesh.getNodes(), ghost_type);
}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::
    computeNormalsOnIntegrationPoints(const Array<Real> & field,
                                      GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  //  Real * coord = mesh.getNodes().data();
  auto spatial_dimension = mesh.getSpatialDimension();

  // allocate the normal arrays
  normals_on_integration_points.initialize(
      *this, _nb_component = spatial_dimension,
      _spatial_dimension = element_dimension, _ghost_type = ghost_type,
      _element_kind = kind);

  // loop over the type to build the normals
  for (const auto & type :
       mesh.elementTypes(element_dimension, ghost_type, kind)) {
    auto & normals_on_quad = normals_on_integration_points(type, ghost_type);
    computeNormalsOnIntegrationPoints(field, normals_on_quad, type, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
template <ElementType type, ElementKind kind_,
          std::enable_if_t<kind_ == _ek_regular and type != _point_1> *>
void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::
    computeNormalsOnIntegrationPoints(const Array<Real> & field,
                                      Array<Real> & normal,
                                      GhostType ghost_type) const {
  AKANTU_DEBUG_IN();

  auto spatial_dimension = mesh.getSpatialDimension();
  constexpr auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  auto nb_points = getNbIntegrationPoints(type, ghost_type);

  auto nb_element = mesh.getConnectivity(type, ghost_type).size();
  normal.resize(nb_element * nb_points);

  Array<Real> f_el(0, spatial_dimension * nb_nodes_per_element);
  FEEngine::extractNodalToElementField(mesh, field, f_el, type, ghost_type);

  const auto & quads =
      integrator.template getIntegrationPoints<type>(ghost_type);

  for (auto && data : zip(make_view(normal, spatial_dimension, nb_points),
                          make_view<Eigen::Dynamic, nb_nodes_per_element>(
                              f_el, spatial_dimension, nb_nodes_per_element))) {
    ElementClass<type>::computeNormalsOnNaturalCoordinates(
        quads, std::get<1>(data), std::get<0>(data));
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
template <ElementType type, ElementKind kind_,
          std::enable_if_t<kind_ == _ek_regular and type == _point_1> *>
inline void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::
    computeNormalsOnIntegrationPoints(const Array<Real> & /*field*/,
                                      Array<Real> & normal,
                                      GhostType ghost_type) const {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(mesh.getSpatialDimension() == 1,
                      "Mesh dimension must be 1 to compute normals on points!");

  auto spatial_dimension = mesh.getSpatialDimension();
  // Int nb_nodes_per_element  = Mesh::getNbNodesPerElement(type);
  auto nb_points = getNbIntegrationPoints(type, ghost_type);
  const auto & connectivity = mesh.getConnectivity(type, ghost_type);
  auto nb_element = connectivity.size();

  normal.resize(nb_element * nb_points);
  auto normals_on_quad =
      make_view(normal, spatial_dimension, nb_points).begin();
  const auto & segments = mesh.getElementToSubelement(type, ghost_type);
  const auto & coords = mesh.getNodes();

  const Mesh * mesh_segment;
  if (mesh.isMeshFacets()) {
    mesh_segment = &(mesh.getMeshParent());
  } else {
    mesh_segment = &mesh;
  }

  for (Idx elem = 0; elem < nb_element; ++elem) {
    auto nb_segment = segments(elem).size();
    AKANTU_DEBUG_ASSERT(
        nb_segment > 0,
        "Impossible to compute a normal on a point connected to 0 segments");

    Real normal_value = 1;
    if (nb_segment == 1) {
      auto point = connectivity(elem);
      const auto segment = segments(elem)[0];
      const auto & segment_connectivity =
          mesh_segment->getConnectivity(segment.type, segment.ghost_type);
      Vector<Idx> segment_points = segment_connectivity.begin(
          Mesh::getNbNodesPerElement(segment.type))[segment.element];
      Real difference;
      if (segment_points(0) == point) {
        difference = coords(elem) - coords(segment_points(1));
      } else {
        difference = coords(elem) - coords(segment_points(0));
      }

      normal_value = difference / std::abs(difference);
    }

    for (Idx n(0); n < nb_points; ++n) {
      (*normals_on_quad)(0, n) = normal_value;
    }
    ++normals_on_quad;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::
    computeNormalsOnIntegrationPoints(const Array<Real> & field,
                                      Array<Real> & normal, ElementType type,
                                      GhostType ghost_type) const {
  auto && call = [&](auto && enum_type) {
    constexpr ElementType type = std ::decay_t<decltype(enum_type)>::value;
    this->computeNormalsOnIntegrationPoints<type>(field, normal, ghost_type);
  };
  tuple_dispatch<ElementTypes_t<kind>>(call, type);
}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
inline void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::inverseMap(
    const Ref<const VectorXr> real_coords, Int element, ElementType type,
    Ref<VectorXr> natural_coords, GhostType ghost_type) const {
  /// need sfinea to avoid structural
  auto && call = [&](auto && enum_type) {
    constexpr ElementType type = std ::decay_t<decltype(enum_type)>::value;
    shape_functions.template inverseMap<type>(real_coords, element,
                                              natural_coords, ghost_type);
  };
  tuple_dispatch<ElementTypes_t<kind>>(call, type);
}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
inline bool FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::contains(
    const Vector<Real> & real_coords, Int element, ElementType type,
    GhostType ghost_type) const {
  auto && call = [&](auto && enum_type) {
    constexpr ElementType type = std ::decay_t<decltype(enum_type)>::value;
    return shape_functions.template contains<type>(real_coords, element,
                                                   ghost_type);
  };
  return tuple_dispatch<ElementTypes_t<kind>>(call, type);
}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
template <ElementKind kind_, typename D1, typename D2,
          std::enable_if_t<aka::are_vectors<D1, D2>::value and
                           kind_ != _ek_cohesive> *>
inline void
FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::computeShapesImpl(
    const Eigen::MatrixBase<D1> & real_coords, Idx element, ElementType type,
    Eigen::MatrixBase<D2> & shapes, GhostType ghost_type) const {

  auto && call = [&](auto && enum_type) {
    constexpr ElementType type = std ::decay_t<decltype(enum_type)>::value;
    this->shape_functions.template computeShapes<type>(real_coords, element,
                                                       shapes, ghost_type);
  };
  tuple_dispatch<ElementTypes_t<kind>>(call, type);
}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
template <ElementKind kind_, typename D1, typename D2,
          std::enable_if_t<aka::is_vector_v<D1> and kind_ != _ek_cohesive> *>
inline void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::
    computeShapeDerivativesImpl(const Eigen::MatrixBase<D1> & real_coords,
                                Int element, ElementType type,
                                Eigen::MatrixBase<D2> & shape_derivatives,
                                GhostType ghost_type) const {
  MatrixProxy<const Real> coords_mat(real_coords.derived().data(),
                                     shape_derivatives.rows(), 1);
  Tensor3Proxy<Real> shapesd_tensor(shape_derivatives.derived().data(),
                                    shape_derivatives.rows(),
                                    shape_derivatives.cols(), 1);
  auto && call = [&](auto && enum_type) {
    constexpr ElementType type = std ::decay_t<decltype(enum_type)>::value;
    shape_functions.template computeShapeDerivatives<type>(
        coords_mat, element, shapesd_tensor, ghost_type);
  };
  tuple_dispatch<ElementTypes_t<kind>>(call, type);
}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
inline Int
FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::getNbIntegrationPoints(
    ElementType type, GhostType ghost_type) const {
  return tuple_dispatch<ElementTypes_t<kind>>(
      [&](auto && enum_type) {
        constexpr ElementType type = std ::decay_t<decltype(enum_type)>::value;
        return integrator.template getNbIntegrationPoints<type>(ghost_type);
      },
      type);
}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
inline const Array<Real> &
FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::getShapes(
    ElementType type, GhostType ghost_type, Int /*id*/) const {
  return shape_functions.getShapes(type, ghost_type);
}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
inline const Array<Real> &
FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::getShapesDerivatives(
    ElementType type, GhostType ghost_type, Int /*id*/) const {
  return shape_functions.getShapesDerivatives(type, ghost_type);
}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
inline const Matrix<Real> &
FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::getIntegrationPoints(
    ElementType type, GhostType ghost_type) const {

  return tuple_dispatch<ElementTypes_t<kind>>(
      [&](auto && enum_type) -> const Matrix<Real> & {
        constexpr ElementType type = std ::decay_t<decltype(enum_type)>::value;
        return (integrator.template getIntegrationPoints<type>(ghost_type));
      },
      type);
}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::printself(
    std::ostream & stream, int indent) const {
  std::string space(indent, AKANTU_INDENT);

  stream << space << "FEEngineTemplate [" << std::endl;
  stream << space << " + parent [" << std::endl;
  FEEngine::printself(stream, indent + 3);
  stream << space << "   ]" << std::endl;
  stream << space << " + shape functions [" << std::endl;
  shape_functions.printself(stream, indent + 3);
  stream << space << "   ]" << std::endl;
  stream << space << " + integrator [" << std::endl;
  integrator.printself(stream, indent + 3);
  stream << space << "   ]" << std::endl;
  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::onElementsAdded(
    const Array<Element> & new_elements, const NewElementsEvent & /*unused*/) {
  integrator.onElementsAdded(new_elements);

  const auto & points = integrator.getIntegrationPoints();
  // for each distinct type add missing integration points to shape_functions
  for (auto ghost_type : ghost_types) {
    for (const auto & type : points.elementTypes(_ghost_type = ghost_type)) {
      tuple_dispatch<ElementTypes_t<kind>>(
          [&](auto && enum_type) {
            constexpr ElementType type =
                std ::decay_t<decltype(enum_type)>::value;
            shape_functions.template setIntegrationPointsByType<type>(
                points(type, ghost_type), ghost_type);
          },
          type);
    }
  }

  shape_functions.onElementsAdded(new_elements);
}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::onElementsRemoved(
    const Array<Element> & removed_elements,
    const ElementTypeMapArray<Idx> & new_numbering,
    const RemovedElementsEvent & /*event*/) {
  integrator.onElementsRemoved(removed_elements, new_numbering);
  shape_functions.onElementsRemoved(removed_elements, new_numbering);
}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::onElementsChanged(
    const Array<Element> & /*unused*/, const Array<Element> & /*unused*/,
    const ElementTypeMapArray<Idx> & /*unused*/,
    const ChangedElementsEvent & /*unused*/) {}

} // namespace akantu
