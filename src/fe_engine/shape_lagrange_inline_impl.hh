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
#include "aka_iterators.hh"
#include "aka_voigthelper.hh"
#include "fe_engine.hh"
// #include "shape_lagrange.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_SHAPE_LAGRANGE_INLINE_IMPL_HH_
#define AKANTU_SHAPE_LAGRANGE_INLINE_IMPL_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <typename D>
inline void ShapeLagrange<kind>::initShapeFunctions(
    const Array<Real> & /*nodes*/,
    const Eigen::MatrixBase<D> & /*integration_points*/, ElementType /*type*/,
    GhostType /*ghost_type*/) {}

/* -------------------------------------------------------------------------- */
template <>
template <typename D>
inline void ShapeLagrange<_ek_regular>::initShapeFunctions(
    const Array<Real> & nodes, const Eigen::MatrixBase<D> & integration_points,
    ElementType type, GhostType ghost_type) {
  tuple_dispatch<ElementTypes_t<_ek_regular>>(
      [&](auto && enum_type) {
        constexpr ElementType type = std::decay_t<decltype(enum_type)>::value;
        this->setIntegrationPointsByType<type>(integration_points, ghost_type);
        this->precomputeShapesOnIntegrationPoints<type>(nodes, ghost_type);
        if (ElementClass<type>::getNaturalSpaceDimension() ==
            mesh.getSpatialDimension()) {
          this->precomputeShapeDerivativesOnIntegrationPoints<type>(nodes,
                                                                    ghost_type);
        }
      },
      type);
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type, typename D1, typename D2>
inline void ShapeLagrange<kind>::computeShapeDerivativesOnCPointsByElement(
    const Eigen::MatrixBase<D1> & node_coords,
    const Eigen::MatrixBase<D2> & natural_coords,
    Tensor3Base<Real> & shapesd) const {
  AKANTU_DEBUG_IN();

  // compute dnds
  Tensor3<Real> dnds(node_coords.rows(), node_coords.cols(),
                     natural_coords.cols());
  ElementClass<type>::computeDNDS(natural_coords, dnds);
  // compute jacobian
  Tensor3<Real> J(node_coords.rows(), natural_coords.rows(),
                  natural_coords.cols());
  ElementClass<type>::computeJMat(dnds, node_coords, J);

  // compute dndx
  ElementClass<type>::computeShapeDerivatives(J, dnds, shapesd);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type, typename D1, typename D2>
void ShapeLagrange<kind>::inverseMap(
    const Eigen::MatrixBase<D1> & real_coords, Int elem,
    const Eigen::MatrixBase<D2> & natural_coords_, GhostType ghost_type) const {
  AKANTU_DEBUG_IN();

  // as advised by the Eigen developers even though this is a UB
  auto & natural_coords = const_cast<Eigen::MatrixBase<D2> &>(natural_coords_);

  auto nodes_coord = mesh.extractNodalValuesFromElement(
      mesh.getNodes(), Element{type, elem, ghost_type});

  ElementClass<type>::inverseMap(real_coords, nodes_coord, natural_coords);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type, typename D, std::enable_if_t<aka::is_vector_v<D>> *>
bool ShapeLagrange<kind>::contains(const Eigen::MatrixBase<D> & real_coords,
                                   Idx elem, GhostType ghost_type) const {
  auto spatial_dimension = mesh.getSpatialDimension();
  Vector<Real> natural_coords(spatial_dimension);

  inverseMap<type>(real_coords, elem, natural_coords, ghost_type);
  return ElementClass<type>::contains(natural_coords);
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type, typename D1, typename D2, typename D3,
          std::enable_if_t<aka::are_vectors<D1, D3>::value> *>
void ShapeLagrange<kind>::interpolate(
    const Eigen::MatrixBase<D1> & real_coords, Idx elem,
    const Eigen::MatrixBase<D2> & nodal_values,
    Eigen::MatrixBase<D3> & interpolated, GhostType ghost_type) const {
  constexpr auto nb_shapes = ElementClass<type>::getShapeSize();
  Vector<Real, nb_shapes> shapes;
  computeShapes<type>(real_coords, elem, shapes, ghost_type);
  ElementClass<type>::interpolate(nodal_values, shapes, interpolated);
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type, typename D1, typename D2>
void ShapeLagrange<kind>::computeShapes(
    const Eigen::MatrixBase<D1> & real_coords, Idx elem,
    Eigen::MatrixBase<D2> & shapes, GhostType ghost_type) const {
  AKANTU_DEBUG_IN();

  auto spatial_dimension = mesh.getSpatialDimension();
  Vector<Real> natural_coords(spatial_dimension);

  inverseMap<type>(real_coords, elem, natural_coords, ghost_type);
  ElementClass<type>::computeShapes(natural_coords, shapes);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type, typename D>
void ShapeLagrange<kind>::computeShapeDerivatives(
    const Eigen::MatrixBase<D> & real_coords, Idx elem,
    Tensor3Base<Real> & shapesd, GhostType ghost_type) const {
  AKANTU_DEBUG_IN();

  auto spatial_dimension = mesh.getSpatialDimension();
  auto nb_points = real_coords.cols();

#if !defined(AKANTU_NDEBUG)
  const auto nb_nodes_per_element =
      ElementClass<type>::getNbNodesPerInterpolationElement();
  AKANTU_DEBUG_ASSERT(mesh.getSpatialDimension() == shapesd.size(0) &&
                          nb_nodes_per_element == shapesd.size(1),
                      "Shape size doesn't match");
  AKANTU_DEBUG_ASSERT(nb_points == shapesd.size(2),
                      "Number of points doesn't match shapes size");
#endif

  Matrix<Real> natural_coords(spatial_dimension, nb_points);

  // Creates the matrix of natural coordinates
  for (Int i = 0; i < nb_points; i++) {
    inverseMap<type>(real_coords(i), elem, natural_coords(i), ghost_type);
  }

  auto nodes_coord = mesh.extractNodalValuesFromElement(
      mesh.getNodes(), Element{type, elem, ghost_type});

  computeShapeDerivativesOnCPointsByElement<type>(nodes_coord, natural_coords,
                                                  shapesd);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
ShapeLagrange<kind>::ShapeLagrange(const Mesh & mesh, Int spatial_dimension,
                                   const ID & id)
    : ShapeLagrangeBase(mesh, spatial_dimension, kind, id) {}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type, typename D>
void ShapeLagrange<kind>::computeShapeDerivativesOnIntegrationPoints(
    const Array<Real> & nodes, const Eigen::MatrixBase<D> & integration_points,
    Array<Real> & shape_derivatives, GhostType ghost_type,
    const Array<Idx> & filter_elements) const {
  AKANTU_DEBUG_IN();

  auto spatial_dimension = mesh.getSpatialDimension();

  auto nb_nodes_per_element =
      ElementClass<type>::getNbNodesPerInterpolationElement();

  auto nb_points = integration_points.cols();
  auto nb_element = mesh.getConnectivity(type, ghost_type).size();

  auto size_of_shapesd = ElementClass<type>::getShapeDerivativesSize();
  AKANTU_DEBUG_ASSERT(shape_derivatives.getNbComponent() == size_of_shapesd,
                      "The shapes_derivatives array does not have the correct "
                          << "number of component");
  shape_derivatives.resize(nb_element * nb_points);

  Array<Real> x_el(0, spatial_dimension * nb_nodes_per_element);
  FEEngine::extractNodalToElementField(mesh, nodes, x_el, type, ghost_type,
                                       filter_elements);

  auto * shapesd_val = shape_derivatives.data();
  auto x_it = x_el.begin(spatial_dimension, nb_nodes_per_element);

  if (filter_elements != empty_filter) {
    nb_element = filter_elements.size();
  }

  for (Int elem = 0; elem < nb_element; ++elem, ++x_it) {
    if (filter_elements != empty_filter) {
      shapesd_val = shape_derivatives.data() +
                    filter_elements(elem) * size_of_shapesd * nb_points;
    }

    auto & X = *x_it;
    Tensor3Proxy<Real> B(shapesd_val, spatial_dimension, nb_nodes_per_element,
                         nb_points);
    computeShapeDerivativesOnCPointsByElement<type>(X, integration_points, B);

    if (filter_elements == empty_filter) {
      shapesd_val += size_of_shapesd * nb_points;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
void ShapeLagrange<kind>::computeShapeDerivativesOnIntegrationPoints(
    const Array<Real> & nodes, const Ref<const MatrixXr> integration_points,
    Array<Real> & shape_derivatives, ElementType type, GhostType ghost_type,
    const Array<Idx> & filter_elements) const {
  tuple_dispatch<ElementTypes_t<_ek_regular>>(
      [&](auto && enum_type) {
        constexpr ElementType type = aka::decay_v<decltype(enum_type)>;
        this->computeShapeDerivativesOnIntegrationPoints<type>(
            nodes, integration_points, shape_derivatives, ghost_type,
            filter_elements);
      },
      type);
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeLagrange<kind>::precomputeShapesOnIntegrationPoints(
    const Array<Real> & nodes, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  InterpolationType itp_type = ElementClassProperty<type>::interpolation_type;
  Matrix<Real> & natural_coords = integration_points(type, ghost_type);
  auto size_of_shapes = ElementClass<type>::getShapeSize();

  Array<Real> & shapes_tmp =
      shapes.alloc(0, size_of_shapes, itp_type, ghost_type);

  this->computeShapesOnIntegrationPoints<type>(nodes, natural_coords,
                                               shapes_tmp, ghost_type);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeLagrange<kind>::precomputeShapeDerivativesOnIntegrationPoints(
    const Array<Real> & nodes, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  InterpolationType itp_type = ElementClassProperty<type>::interpolation_type;
  Matrix<Real> & natural_coords = integration_points(type, ghost_type);
  auto size_of_shapesd = ElementClass<type>::getShapeDerivativesSize();

  Array<Real> & shapes_derivatives_tmp =
      shapes_derivatives.alloc(0, size_of_shapesd, itp_type, ghost_type);

  this->computeShapeDerivativesOnIntegrationPoints<type>(
      nodes, natural_coords, shapes_derivatives_tmp, ghost_type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeLagrange<kind>::interpolateOnIntegrationPoints(
    const Array<Real> & in_u, Array<Real> & out_uq, Int nb_degree_of_freedom,
    const Array<Real> & shapes, GhostType ghost_type,
    const Array<Idx> & filter_elements) const {
  AKANTU_DEBUG_IN();

  auto nb_nodes_per_element =
      ElementClass<type>::getNbNodesPerInterpolationElement();

  Array<Real> u_el(0, nb_degree_of_freedom * nb_nodes_per_element);
  FEEngine::extractNodalToElementField(mesh, in_u, u_el, type, ghost_type,
                                       filter_elements);

  this->interpolateElementalFieldOnIntegrationPoints<type>(
      u_el, out_uq, ghost_type, shapes, filter_elements);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeLagrange<kind>::interpolateOnIntegrationPoints(
    const Array<Real> & in_u, Array<Real> & out_uq, Int nb_degree_of_freedom,
    GhostType ghost_type, const Array<Idx> & filter_elements) const {
  AKANTU_DEBUG_IN();

  InterpolationType itp_type = ElementClassProperty<type>::interpolation_type;
  AKANTU_DEBUG_ASSERT(shapes.exists(itp_type, ghost_type),
                      "No shapes for the type "
                          << shapes.printType(itp_type, ghost_type));

  this->interpolateOnIntegrationPoints<type>(in_u, out_uq, nb_degree_of_freedom,
                                             shapes(itp_type, ghost_type),
                                             ghost_type, filter_elements);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeLagrange<kind>::gradientOnIntegrationPoints(
    const Array<Real> & in_u, Array<Real> & out_nablauq,
    Int nb_degree_of_freedom, GhostType ghost_type,
    const Array<Idx> & filter_elements) const {
  AKANTU_DEBUG_IN();

  InterpolationType itp_type = ElementClassProperty<type>::interpolation_type;
  AKANTU_DEBUG_ASSERT(
      shapes_derivatives.exists(itp_type, ghost_type),
      "No shapes derivatives for the type "
          << shapes_derivatives.printType(itp_type, ghost_type));

  auto nb_nodes_per_element =
      ElementClass<type>::getNbNodesPerInterpolationElement();

  Array<Real> u_el(0, nb_degree_of_freedom * nb_nodes_per_element);
  FEEngine::extractNodalToElementField(mesh, in_u, u_el, type, ghost_type,
                                       filter_elements);

  this->gradientElementalFieldOnIntegrationPoints<type>(
      u_el, out_nablauq, ghost_type, shapes_derivatives(itp_type, ghost_type),
      filter_elements);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeLagrange<kind>::computeBtD(const Array<Real> & Ds, Array<Real> & BtDs,
                                     GhostType ghost_type,
                                     const Array<Idx> & filter_elements) const {
  BtDs.resize(Ds.size());

  auto itp_type = ElementClassProperty<type>::interpolation_type;
  const auto & shapes_derivatives =
      this->shapes_derivatives(itp_type, ghost_type);

  auto spatial_dimension = mesh.getSpatialDimension();
  auto nb_nodes_per_element = mesh.getNbNodesPerElement(type);

  Array<Real> shapes_derivatives_filtered(0,
                                          shapes_derivatives.getNbComponent());
  auto view = make_const_view(shapes_derivatives, spatial_dimension,
                              nb_nodes_per_element);

  if (filter_elements != empty_filter) {
    FEEngine::filterElementalData(this->mesh, shapes_derivatives,
                                  shapes_derivatives_filtered, type, ghost_type,
                                  filter_elements);
    view = make_const_view(shapes_derivatives_filtered, spatial_dimension,
                           nb_nodes_per_element);
  }

  for (auto && [B, D, Bt_D] :
       zip(view,
           make_view(Ds, Ds.getNbComponent() / spatial_dimension,
                     spatial_dimension),
           make_view(BtDs, BtDs.getNbComponent() / nb_nodes_per_element,
                     nb_nodes_per_element))) {
    // transposed due to the storage layout of B
    Bt_D.noalias() = D * B;
  }
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <
    ElementType type,
    std::enable_if_t<ElementClass<type>::getNaturalSpaceDimension() != 0> *>
void ShapeLagrange<kind>::computeBtDB(
    const Array<Real> & Ds, Array<Real> & BtDBs, Int order_d,
    GhostType ghost_type, const Array<Idx> & filter_elements) const {
  BtDBs.resize(Ds.size());

  auto itp_type = ElementClassProperty<type>::interpolation_type;
  const auto & shapes_derivatives =
      this->shapes_derivatives(itp_type, ghost_type);

  constexpr auto dim = ElementClass<type>::getSpatialDimension();
  auto nb_nodes_per_element = mesh.getNbNodesPerElement(type);

  Array<Real> shapes_derivatives_filtered(0,
                                          shapes_derivatives.getNbComponent());
  auto && view = make_const_view(shapes_derivatives, dim, nb_nodes_per_element);

  if (filter_elements != empty_filter) {
    FEEngine::filterElementalData(this->mesh, shapes_derivatives,
                                  shapes_derivatives_filtered, type, ghost_type,
                                  filter_elements);
    view =
        make_const_view(shapes_derivatives_filtered, dim, nb_nodes_per_element);
  }

  if (order_d == 4) {
    auto tangent_size = VoigtHelper<dim>::size;
    Matrix<Real> B(tangent_size, dim * nb_nodes_per_element);

    for (auto && values : zip(view, make_view(Ds, tangent_size, tangent_size),
                              make_view(BtDBs, dim * nb_nodes_per_element,
                                        dim * nb_nodes_per_element))) {
      const auto & Bfull = std::get<0>(values);
      const auto & D = std::get<1>(values);
      auto & Bt_D_B = std::get<2>(values);

      VoigtHelper<dim>::transferBMatrixToSymVoigtBMatrix(Bfull, B,
                                                         nb_nodes_per_element);
      Bt_D_B.noalias() = B.transpose() * D * B;
    }
  } else if (order_d == 2) {
    for (auto && values :
         zip(view, make_view(Ds, dim, dim),
             make_view(BtDBs, nb_nodes_per_element, nb_nodes_per_element))) {
      const auto & B = std::get<0>(values);
      const auto & D = std::get<1>(values);
      auto & Bt_D_B = std::get<2>(values);
      Bt_D_B.noalias() = B.transpose() * D * B;
    }
  }
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeLagrange<kind>::computeNtbN(
    const Array<Real> & bs, Array<Real> & NtbNs, GhostType ghost_type,
    const Array<Idx> & filter_elements) const {
  NtbNs.resize(bs.size());

  auto itp_type = ElementClassProperty<type>::interpolation_type;
  auto size_of_shapes = ElementClass<type>::getShapeSize();
  auto nb_degree_of_freedom = bs.getNbComponent();

  auto nb_nodes_per_element = mesh.getNbNodesPerElement(type);
  Array<Real> shapes_filtered(0, size_of_shapes);

  auto && view =
      make_const_view(shapes(itp_type, ghost_type), 1, size_of_shapes);

  if (filter_elements != empty_filter) {
    FEEngine::filterElementalData(this->mesh, shapes(itp_type, ghost_type),
                                  shapes_filtered, type, ghost_type,
                                  filter_elements);
    view = make_const_view(shapes_filtered, 1, size_of_shapes);
  }

  for (auto && values :
       zip(view, make_view(bs, nb_degree_of_freedom, 1),
           make_view(NtbNs, nb_nodes_per_element, nb_nodes_per_element))) {
    const auto & N = std::get<0>(values);
    const auto & b = std::get<1>(values);
    auto & Nt_b_N = std::get<2>(values);

    Nt_b_N.noalias() = N.transpose() * b * N;
  }
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeLagrange<kind>::computeNtb(const Array<Real> & bs, Array<Real> & Ntbs,
                                     GhostType ghost_type,
                                     const Array<Idx> & filter_elements) const {
  AKANTU_DEBUG_IN();

  Ntbs.resize(bs.size());

  auto size_of_shapes = ElementClass<type>::getShapeSize();
  auto itp_type = ElementClassProperty<type>::interpolation_type;
  auto nb_degree_of_freedom = bs.getNbComponent();

  Array<Real> shapes_filtered(0, size_of_shapes);
  auto && view =
      make_const_view(shapes(itp_type, ghost_type), 1, size_of_shapes);

  if (filter_elements != empty_filter) {
    FEEngine::filterElementalData(this->mesh, shapes(itp_type, ghost_type),
                                  shapes_filtered, type, ghost_type,
                                  filter_elements);
    view = make_const_view(shapes_filtered, 1, size_of_shapes);
  }

  for (auto && values :
       zip(make_view(bs, nb_degree_of_freedom, 1), view,
           make_view(Ntbs, nb_degree_of_freedom, size_of_shapes))) {
    const auto & b = std::get<0>(values);
    const auto & N = std::get<1>(values);
    auto & Ntb = std::get<2>(values);

    Ntb.noalias() = b * N;
  }

  AKANTU_DEBUG_OUT();
}

} // namespace akantu

#endif /* AKANTU_SHAPE_LAGRANGE_INLINE_IMPL_HH_ */
