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
#include "shape_lagrange_base.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_SHAPE_LAGRANGE_HH_
#define AKANTU_SHAPE_LAGRANGE_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
template <class Shape> class ShapeCohesive;
class ShapeIGFEM;

template <ElementKind kind> class ShapeLagrange : public ShapeLagrangeBase {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  ShapeLagrange(const Mesh & mesh, Int spatial_dimension,
                const ID & id = "shape_lagrange");
  ~ShapeLagrange() override = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// initialization function for structural elements not yet implemented
  template <typename D>
  inline void
  initShapeFunctions(const Array<Real> & nodes,
                     const Eigen::MatrixBase<D> & integration_points,
                     ElementType type, GhostType ghost_type);

  /// computes the shape functions derivatives for given interpolation points
  template <ElementType type, typename D>
  void computeShapeDerivativesOnIntegrationPoints(
      const Array<Real> & nodes,
      const Eigen::MatrixBase<D> & integration_points,
      Array<Real> & shape_derivatives, GhostType ghost_type,
      const Array<Idx> & filter_elements = empty_filter) const;

  void computeShapeDerivativesOnIntegrationPoints(
      const Array<Real> & nodes, const Ref<const MatrixXr> integration_points,
      Array<Real> & shape_derivatives, ElementType type, GhostType ghost_type,
      const Array<Idx> & filter_elements) const override;

  /// pre compute all shapes on the element integration points from natural
  /// coordinates
  template <ElementType type>
  void precomputeShapesOnIntegrationPoints(const Array<Real> & nodes,
                                           GhostType ghost_type);

  /// pre compute all shape derivatives on the element integration points from
  /// natural coordinates
  template <ElementType type>
  void precomputeShapeDerivativesOnIntegrationPoints(const Array<Real> & nodes,
                                                     GhostType ghost_type);

  /// interpolate nodal values on the integration points
  template <ElementType type>
  void interpolateOnIntegrationPoints(
      const Array<Real> & u, Array<Real> & uq, Int nb_degree_of_freedom,
      GhostType ghost_type = _not_ghost,
      const Array<Idx> & filter_elements = empty_filter) const;

  template <ElementType type>
  void interpolateOnIntegrationPoints(
      const Array<Real> & in_u, Array<Real> & out_uq, Int nb_degree_of_freedom,
      const Array<Real> & shapes, GhostType ghost_type = _not_ghost,
      const Array<Idx> & filter_elements = empty_filter) const;

  /// interpolate on physical point
  template <ElementType type, typename D1, typename D2, typename D3,
            std::enable_if_t<aka::are_vectors<D1, D3>::value> * = nullptr>
  void interpolate(const Eigen::MatrixBase<D1> & real_coords, Idx elem,
                   const Eigen::MatrixBase<D2> & nodal_values,
                   Eigen::MatrixBase<D3> & interpolated,
                   GhostType ghost_type) const;

  /// compute the gradient of u on the integration points
  template <ElementType type>
  void gradientOnIntegrationPoints(
      const Array<Real> & u, Array<Real> & nablauq, Int nb_degree_of_freedom,
      GhostType ghost_type = _not_ghost,
      const Array<Idx> & filter_elements = empty_filter) const;

  template <ElementType type>
  void computeBtD(const Array<Real> & Ds, Array<Real> & BtDs,
                  GhostType ghost_type,
                  const Array<Idx> & filter_elements) const;

  template <ElementType type,
            std::enable_if_t<ElementClass<type>::getNaturalSpaceDimension() !=
                             0> * = nullptr>
  void computeBtDB(const Array<Real> & Ds, Array<Real> & BtDBs, Int order_d,
                   GhostType ghost_type,
                   const Array<Idx> & filter_elements) const;

  template <ElementType type,
            std::enable_if_t<ElementClass<type>::getNaturalSpaceDimension() ==
                             0> * = nullptr>
  void computeBtDB(const Array<Real> & /*Ds*/, Array<Real> & /*BtDBs*/,
                   Int /*order_d*/, GhostType /*ghost_type*/,
                   const Array<Idx> & /*filter_elements*/) const {
    AKANTU_TO_IMPLEMENT();
  }

  /// multiply a field by shape functions  @f$ fts_{ij} = f_i * \varphi_j @f$
  template <ElementType type>
  void computeNtb(const Array<Real> & bs, Array<Real> & Ntbs,
                  GhostType ghost_type,
                  const Array<Idx> & filter_elements = empty_filter) const;

  template <ElementType type>
  void computeNtbN(const Array<Real> & bs, Array<Real> & NtbNs,
                   GhostType ghost_type,
                   const Array<Idx> & filter_elements) const;

  /// find natural coords from real coords provided an element
  template <ElementType type, typename D1, typename D2>
  void inverseMap(const Eigen::MatrixBase<D1> & real_coords, Idx element,
                  const Eigen::MatrixBase<D2> & natural_coords,
                  GhostType ghost_type = _not_ghost) const;

  /// return true if the coordinates provided are inside the element, false
  /// otherwise
  template <ElementType type, typename D>
  bool contains(const Eigen::MatrixBase<D> & real_coords, Idx elem,
                GhostType ghost_type) const;

  /// compute the shape on a provided point
  template <ElementType type, typename D1, typename D2>
  void computeShapes(const Eigen::MatrixBase<D1> & real_coords, Idx elem,
                     Eigen::MatrixBase<D2> & shapes,
                     GhostType ghost_type) const;

  /// compute the shape derivatives on a provided point
  template <ElementType type, typename D>
  void computeShapeDerivatives(const Eigen::MatrixBase<D> & real_coords,
                               Idx elem, Tensor3Base<Real> & shapes,
                               GhostType ghost_type) const;

protected:
  /// compute the shape derivatives on integration points for a given element
  template <ElementType type, typename D1, typename D2>
  inline void computeShapeDerivativesOnCPointsByElement(
      const Eigen::MatrixBase<D1> & node_coords,
      const Eigen::MatrixBase<D2> & natural_coords,
      Tensor3Base<Real> & shapesd) const;
};

} // namespace akantu

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
#include "shape_lagrange_inline_impl.hh"

#endif /* AKANTU_SHAPE_LAGRANGE_HH_ */
