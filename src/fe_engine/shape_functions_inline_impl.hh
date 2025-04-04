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
#include "fe_engine.hh"
#include "shape_functions.hh"
/* -------------------------------------------------------------------------- */

//#ifndef AKANTU_SHAPE_FUNCTIONS_INLINE_IMPL_HH_
//#define AKANTU_SHAPE_FUNCTIONS_INLINE_IMPL_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
inline const Array<Real> &
ShapeFunctions::getShapes(ElementType el_type, GhostType ghost_type) const {
  return shapes(FEEngine::getInterpolationType(el_type), ghost_type);
}

/* -------------------------------------------------------------------------- */
inline const Array<Real> &
ShapeFunctions::getShapesDerivatives(ElementType el_type,
                                     GhostType ghost_type) const {
  return shapes_derivatives(FEEngine::getInterpolationType(el_type),
                            ghost_type);
}

/* -------------------------------------------------------------------------- */
inline Int ShapeFunctions::getShapeSize(ElementType type) {
  return tuple_dispatch<AllElementTypes>(
      [](auto type) { return ElementClass<type.value>::getShapeSize(); }, type);
}

/* -------------------------------------------------------------------------- */
inline Int ShapeFunctions::getShapeDerivativesSize(ElementType type) {
  return tuple_dispatch<AllElementTypes>(
      [](auto type) {
        return ElementClass<type.value>::getShapeDerivativesSize();
      },
      type);
}

/* -------------------------------------------------------------------------- */
template <ElementType type, class D1>
void ShapeFunctions::setIntegrationPointsByType(
    const Eigen::MatrixBase<D1> & points, GhostType ghost_type) {
  this->integration_points(type, ghost_type) = points;
}

/* -------------------------------------------------------------------------- */
template <typename D1, typename D2>
inline void ShapeFunctions::buildInterpolationMatrix(
    const Eigen::MatrixBase<D1> & coordinates,
    Eigen::MatrixBase<D2> & coordMatrix, Int integration_order) const {
  switch (integration_order) {
  case 1: {
    for (Int i = 0; i < coordinates.cols(); ++i) {
      coordMatrix(i, 0) = 1;
    }
    break;
  }
  case 2: {
    auto nb_quadrature_points = coordMatrix.cols();

    for (Int i = 0; i < coordinates.cols(); ++i) {
      coordMatrix(i, 0) = 1;
      for (Int j = 1; j < nb_quadrature_points; ++j) {
        coordMatrix(i, j) = coordinates(j - 1, i);
      }
    }
    break;
  }
  default: {
    AKANTU_TO_IMPLEMENT();
    break;
  }
  }
}

/* -------------------------------------------------------------------------- */
template <ElementType type> struct BuildElementalFieldInterpolationMatrix {
  template <typename ShapeFunction, typename D1, typename D2>
  static inline void
  call(ShapeFunction && func, const Eigen::MatrixBase<D1> & coordinates,
       Eigen::MatrixBase<D2> & coord_matrix, Int integration_order) {
    func.buildInterpolationMatrix(coordinates, coord_matrix, integration_order);
  }
};

/**
 * @todo Write a more efficient interpolation for quadrangles by
 * dropping unnecessary quadrature points
 *
 */
/* -------------------------------------------------------------------------- */
template <> struct BuildElementalFieldInterpolationMatrix<_quadrangle_4> {
  template <typename ShapeFunction, typename D1, typename D2>
  static inline void
  call(ShapeFunction && /*func*/, const Eigen::MatrixBase<D1> & coordinates,
       Eigen::MatrixBase<D2> & coord_matrix, Int integration_order) {

    if (integration_order !=
        ElementClassProperty<_quadrangle_4>::polynomial_degree) {
      AKANTU_TO_IMPLEMENT();
    } else {
      for (Int i = 0; i < coordinates.cols(); ++i) {
        auto x = coordinates(0, i);
        auto y = coordinates(1, i);

        coord_matrix(i, 0) = 1;
        coord_matrix(i, 1) = x;
        coord_matrix(i, 2) = y;
        coord_matrix(i, 3) = x * y;
      }
    }
  }
};

/* -------------------------------------------------------------------------- */
template <> struct BuildElementalFieldInterpolationMatrix<_quadrangle_8> {
  template <typename ShapeFunction, typename D1, typename D2>
  static inline void
  call(ShapeFunction && /*func*/, const Eigen::MatrixBase<D1> & coordinates,
       Eigen::MatrixBase<D2> & coordMatrix, Int integration_order) {

    if (integration_order !=
        ElementClassProperty<_quadrangle_8>::polynomial_degree) {
      AKANTU_TO_IMPLEMENT();
    } else {
      for (Int i = 0; i < coordinates.cols(); ++i) {
        // Int j = 0;
        auto x = coordinates(0, i);
        auto y = coordinates(1, i);

        coordMatrix(i, 0) = 1;
        coordMatrix(i, 1) = x;
        coordMatrix(i, 2) = y;
        coordMatrix(i, 3) = x * y;
      }
    }
  }
};

/* -------------------------------------------------------------------------- */
template <ElementType type, typename D1, typename D2>
inline void ShapeFunctions::buildElementalFieldInterpolationMatrix(
    const Eigen::MatrixBase<D1> & coordinates,
    Eigen::MatrixBase<D2> & coordMatrix, Int integration_order) const {
  BuildElementalFieldInterpolationMatrix<type>::call(
      *this, coordinates, coordMatrix, integration_order);
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
inline void ShapeFunctions::interpolateElementalFieldFromIntegrationPoints(
    const Array<Real> & field,
    const Array<Real> & interpolation_points_coordinates_matrices,
    const Array<Real> & quad_points_coordinates_inv_matrices,
    ElementTypeMapArray<Real> & result, GhostType ghost_type,
    const Array<Int> & element_filter) const {
  AKANTU_DEBUG_IN();

  constexpr auto nb_quad_per_element =
      GaussIntegrationElement<type>::getNbQuadraturePoints();
  auto nb_interpolation_points_per_elem =
      interpolation_points_coordinates_matrices.getNbComponent() /
      nb_quad_per_element;

  if (not result.exists(type, ghost_type)) {
    auto nb_element = this->mesh.getNbElement(type, ghost_type);
    result.alloc(nb_element * nb_interpolation_points_per_elem,
                 field.getNbComponent(), type, ghost_type);
  }

  AKANTU_DEBUG_ASSERT(element_filter != empty_filter,
                      "This function does not work without an element_filter");
  //  auto nb_element = element_filter.size();

  Matrix<Real> coefficients(nb_quad_per_element, field.getNbComponent());

  auto result_begin =
      make_view(result(type, ghost_type), field.getNbComponent(),
                nb_interpolation_points_per_elem)
          .begin();

  /// loop over the elements of the current filter and element type
  for (auto && [el, field, itp_coord, inv_quad_coord] :
       zip(element_filter,
           make_view(field, field.getNbComponent(), nb_quad_per_element),
           make_view(interpolation_points_coordinates_matrices,
                     nb_interpolation_points_per_elem, nb_quad_per_element),
           make_view(quad_points_coordinates_inv_matrices, nb_quad_per_element,
                     nb_quad_per_element))) {

    coefficients = inv_quad_coord * field.transpose();
    result_begin[el] = coefficients.transpose() * itp_coord.transpose();
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
inline void ShapeFunctions::interpolateElementalFieldOnIntegrationPoints(
    const Array<Real> & u_el, Array<Real> & uq, GhostType ghost_type,
    const Array<Real> & shapes, const Array<Int> & filter_elements) const {
  auto nb_element = mesh.getNbElement(type, ghost_type);

  if (nb_element == 0) {
    return;
  }

  auto nb_nodes_per_element = ElementClass<type>::getShapeSize();
  auto nb_points = shapes.size() / nb_element;
  auto nb_degree_of_freedom = u_el.getNbComponent() / nb_nodes_per_element;

  auto N_view = make_view(shapes, nb_nodes_per_element, nb_points);
  std::unique_ptr<Array<Real>> filtered_N;
  if (filter_elements != empty_filter) {
    nb_element = filter_elements.size();
    filtered_N = std::make_unique<Array<Real>>(0, shapes.getNbComponent());
    FEEngine::filterElementalData(mesh, shapes, *filtered_N, type, ghost_type,
                                  filter_elements);
    N_view = make_const_view(*filtered_N, nb_nodes_per_element, nb_points);
  }

  uq.resize(nb_element * nb_points);

  for (auto && data :
       zip(N_view, make_view(uq, nb_degree_of_freedom, nb_points),
           make_view(u_el, nb_degree_of_freedom, nb_nodes_per_element))) {
    const auto & u = std::get<2>(data);
    const auto & N = std::get<0>(data);
    auto & uq = std::get<1>(data);

    uq.noalias() = u * N;
  }
}

/* -------------------------------------------------------------------------- */
template <
    ElementType type,
    std::enable_if_t<ElementClass<type>::getNaturalSpaceDimension() != 0> *>
void ShapeFunctions::gradientElementalFieldOnIntegrationPoints(
    const Array<Real> & u_el, Array<Real> & out_nablauq, GhostType ghost_type,
    const Array<Real> & shapes_derivatives,
    const Array<Int> & filter_elements) const {
  AKANTU_DEBUG_IN();

  constexpr auto nb_nodes_per_element =
      ElementClass<type>::getNbNodesPerInterpolationElement();
  constexpr auto element_dimension =
      ElementClass<type>::getNaturalSpaceDimension();

  auto nb_points = integration_points(type, ghost_type).cols();
  auto nb_element = mesh.getNbElement(type, ghost_type);
  auto nb_degree_of_freedom = u_el.getNbComponent() / nb_nodes_per_element;

  auto B_it = make_const_view<element_dimension, nb_nodes_per_element>(
                  shapes_derivatives)
                  .begin();

  Array<Real> * filtered_B = nullptr;
  if (filter_elements != empty_filter) {
    nb_element = filter_elements.size();
    filtered_B = new Array<Real>(0, shapes_derivatives.getNbComponent());
    FEEngine::filterElementalData(mesh, shapes_derivatives, *filtered_B, type,
                                  ghost_type, filter_elements);
    B_it = make_const_view<element_dimension, nb_nodes_per_element>(*filtered_B)
               .begin();
  }

  out_nablauq.resize(nb_element * nb_points);
  auto u_it = make_view<Eigen::Dynamic, nb_nodes_per_element>(
                  u_el, nb_degree_of_freedom, nb_nodes_per_element)
                  .begin();
  auto nabla_u_it = make_view<Eigen::Dynamic, element_dimension>(
                        out_nablauq, nb_degree_of_freedom, element_dimension)
                        .begin();

  for (Int el = 0; el < nb_element; ++el, ++u_it) {
    const auto & u = *u_it;
    for (Int q = 0; q < nb_points; ++q, ++B_it, ++nabla_u_it) {
      const auto & B = *B_it;
      auto & nabla_u = *nabla_u_it;
      nabla_u.noalias() = u * B.transpose();
    }
  }

  delete filtered_B;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
} // namespace akantu

// #endif /* AKANTU_SHAPE_FUNCTIONS_INLINE_IMPL_HH_ */
