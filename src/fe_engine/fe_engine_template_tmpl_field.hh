/**
 * Copyright (©) 2017-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
// #include "fe_engine_template.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_FE_ENGINE_TEMPLATE_TMPL_FIELD_HH_
#define AKANTU_FE_ENGINE_TEMPLATE_TMPL_FIELD_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
/* Matrix lumping functions                                                   */
/* -------------------------------------------------------------------------- */
namespace fe_engine::details {
template <class Functor>
void fillField(const Functor & field_funct, Array<Real> & field, Int nb_element,
               Int nb_integration_points, ElementType type,
               GhostType ghost_type) {
  auto nb_degree_of_freedom = field.getNbComponent();
  field.resize(nb_integration_points * nb_element);

  Matrix<Real> mat(nb_degree_of_freedom, nb_integration_points);

  Element el{type, 0, ghost_type};
  for (auto && [e, f] : enumerate(
           make_view(field, nb_degree_of_freedom, nb_integration_points))) {
    el.element = e;
    mat.zero();
    field_funct(mat, el);
    f = mat;
  }
}
} // namespace fe_engine::details

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IOF>
void FEEngineTemplate<I, S, kind, IOF>::assembleFieldLumped(
    const std::function<void(Matrix<Real> &, const Element &)> & field_funct,
    const ID & matrix_id, const ID & dof_id, DOFManager & dof_manager,
    ElementType type, GhostType ghost_type) const {
  tuple_dispatch<ElementTypes_t<kind>>(
      [&](auto && enum_type) {
        constexpr ElementType type = aka::decay_v<decltype(enum_type)>;
        this->assembleFieldLumped<type>(field_funct, matrix_id, dof_id,
                                        dof_manager, ghost_type);
      },
      type);
}

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
template <ElementType type>
void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::assembleFieldLumped(
    const std::function<void(Matrix<Real> &, const Element &)> & field_funct,
    const ID & matrix_id, const ID & dof_id, DOFManager & dof_manager,
    GhostType ghost_type) const {
  auto nb_degree_of_freedom = dof_manager.getDOFs(dof_id).getNbComponent();
  auto nb_element = mesh.getNbElement(type, ghost_type);
  auto nb_integration_points = this->getNbIntegrationPoints(type);

  Array<Real> field(0, nb_degree_of_freedom);
  fe_engine::details::fillField(field_funct, field, nb_element,
                                nb_integration_points, type, ghost_type);

  switch (type) {
  case _triangle_6:
  case _quadrangle_8:
  case _tetrahedron_10:
  case _hexahedron_20:
  case _pentahedron_15:
    this->template assembleLumpedDiagonalScaling<type>(field, matrix_id, dof_id,
                                                       dof_manager, ghost_type);
    break;
  default:
    this->template assembleLumpedRowSum<type>(field, matrix_id, dof_id,
                                              dof_manager, ghost_type);
  }
}

/* -------------------------------------------------------------------------- */
/**
 * @f$ \tilde{M}_{i} = \sum_j M_{ij} = \sum_j \int \rho \varphi_i \varphi_j dV =
 * \int \rho \varphi_i dV @f$
 */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
template <ElementType type>
void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::
    assembleLumpedRowSum(const Array<Real> & field, const ID & matrix_id,
                         const ID & dof_id, DOFManager & dof_manager,
                         GhostType ghost_type) const {
  auto shapes_size = ElementClass<type>::getShapeSize();
  auto nb_degree_of_freedom = field.getNbComponent();

  auto field_times_shapes =
      std::make_shared<Array<Real>>(0, shapes_size * nb_degree_of_freedom);

  shape_functions.template computeNtb<type>(field, *field_times_shapes,
                                            ghost_type);

  auto nb_element = mesh.getNbElement(type, ghost_type);
  auto int_field_times_shapes = std::make_shared<Array<Real>>(
      nb_element, shapes_size * nb_degree_of_freedom, "inte_rho_x_shapes");

  integrator.template integrate<type>(
      *field_times_shapes, *int_field_times_shapes,
      nb_degree_of_freedom * shapes_size, ghost_type, empty_filter);

  dof_manager.assembleElementalArrayToLumpedMatrix(
      dof_id, *int_field_times_shapes, matrix_id, type, ghost_type);
}

/* -------------------------------------------------------------------------- */
/**
 * @f$ \tilde{M}_{i} = c * M_{ii} = \int_{V_e} \rho dV @f$
 */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
template <ElementType type>
void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::
    assembleLumpedDiagonalScaling(const Array<Real> & field,
                                  const ID & matrix_id, const ID & dof_id,
                                  DOFManager & dof_manager,
                                  GhostType ghost_type) const {
  const auto & type_p1 = ElementClass<type>::getP1ElementType();
  auto nb_nodes_per_element_p1 = Mesh::getNbNodesPerElement(type_p1);
  auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  auto nb_degree_of_freedom = field.getNbComponent();
  auto nb_element = mesh.getNbElement(type, ghost_type);

  Vector<Real> nodal_factor(nb_nodes_per_element);

  auto assign_weight_to_nodes = [&nodal_factor, nb_nodes_per_element,
                                 nb_nodes_per_element_p1](Real corner,
                                                          Real mid) {
    for (Int n = 0; n < nb_nodes_per_element_p1; n++) {
      nodal_factor(n) = corner;
    }
    for (Int n = nb_nodes_per_element_p1; n < nb_nodes_per_element; n++) {
      nodal_factor(n) = mid;
    }
  };

  if (type == _triangle_6) {
    assign_weight_to_nodes(1. / 12., 1. / 4.);
  }
  if (type == _tetrahedron_10) {
    assign_weight_to_nodes(1. / 32., 7. / 48.);
  }
  if (type == _quadrangle_8) {
    assign_weight_to_nodes(
        3. / 76.,
        16. / 76.); /** coeff. derived by scaling
                     * the diagonal terms of the corresponding
                     * consistent mass computed with 3x3 gauss points;
                     * coeff. are (1./36., 8./36.) for 2x2 gauss points */
  }
  if (type == _hexahedron_20) {
    assign_weight_to_nodes(
        7. / 248., 16. / 248.); /** coeff. derived by scaling
                                 * the diagonal terms of the corresponding
                                 * consistent mass computed with 3x3x3 gauss
                                 * points; coeff. are (1./40.,
                                 * 1./15.) for 2x2x2 gauss points */
  }
  if (type == _pentahedron_15) {
    // coefficients derived by scaling the diagonal terms of the corresponding
    // consistent mass computed with 8 gauss points;
    for (Int n = 0; n < nb_nodes_per_element_p1; n++) {
      nodal_factor(n) = 51. / 2358.;
    }

    Real mid_triangle = 192. / 2358.;
    Real mid_quadrangle = 300. / 2358.;

    nodal_factor(6) = mid_triangle;
    nodal_factor(7) = mid_triangle;
    nodal_factor(8) = mid_triangle;
    nodal_factor(9) = mid_quadrangle;
    nodal_factor(10) = mid_quadrangle;
    nodal_factor(11) = mid_quadrangle;
    nodal_factor(12) = mid_triangle;
    nodal_factor(13) = mid_triangle;
    nodal_factor(14) = mid_triangle;
  }

  if (nb_element == 0) {
    return;
  }

  /// compute @f$ \int \rho dV = \rho V @f$ for each element
  auto int_field = std::make_unique<Array<Real>>(
      field.size(), nb_degree_of_freedom, "inte_rho_x");
  integrator.template integrate<type>(field, *int_field, nb_degree_of_freedom,
                                      ghost_type, empty_filter);

  /// distribute the mass of the element to the nodes
  auto lumped_per_node = std::make_unique<Array<Real>>(
      nb_element, nb_degree_of_freedom * nb_nodes_per_element, "mass_per_node");
  auto int_field_it = int_field->begin(nb_degree_of_freedom);
  auto lumped_per_node_it =
      lumped_per_node->begin(nb_degree_of_freedom, nb_nodes_per_element);

  for (Int e = 0; e < nb_element; ++e) {
    for (Int n = 0; n < nb_nodes_per_element; ++n) {
      auto && l = (*lumped_per_node_it)(n);
      l = *int_field_it * nodal_factor(n);
    }
    ++int_field_it;
    ++lumped_per_node_it;
  }

  dof_manager.assembleElementalArrayToLumpedMatrix(dof_id, *lumped_per_node,
                                                   matrix_id, type, ghost_type);
}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IOF>
template <ElementKind kind_, std::enable_if_t<kind_ != _ek_cohesive> *>
void FEEngineTemplate<I, S, kind, IOF>::assembleFieldMatrixImpl(
    const std::function<void(Matrix<Real> &, const Element &)> & field_funct,
    const ID & matrix_id, const ID & dof_id, DOFManager & dof_manager,
    ElementType type, GhostType ghost_type) const {
  tuple_dispatch<ElementTypes_t<kind>>(
      [&](auto && enum_type) {
        constexpr ElementType type = aka::decay_v<decltype(enum_type)>;
        this->assembleFieldMatrix<type>(field_funct, matrix_id, dof_id,
                                        dof_manager, ghost_type);
      },
      type);
}

template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IOF>
template <ElementKind kind_, std::enable_if_t<kind_ == _ek_cohesive> *>
void FEEngineTemplate<I, S, kind, IOF>::assembleFieldMatrixImpl(
    const std::function<void(Matrix<Real> &,
                             const Element &)> & /*field_funct*/,
    const ID & /*matrix_id*/, const ID & /*dof_id*/,
    DOFManager & /*dof_manager*/, ElementType /*type*/,
    GhostType /*ghost_type*/) const {
  AKANTU_TO_IMPLEMENT();
}

namespace fe_engine {
namespace details {
template <ElementKind kind> struct ShapesForMassHelper {
  template <ElementType type, class ShapeFunctions>
  static auto getShapes(ShapeFunctions & shape_functions,
                        const Matrix<Real> & integration_points,
                        const Array<Real> & nodes, Int & nb_degree_of_freedom,
                        Int nb_element, GhostType ghost_type) {

    auto shapes_size = ElementClass<type>::getShapeSize();
    Array<Real> shapes(0, shapes_size);

    shape_functions.template computeShapesOnIntegrationPoints<type>(
        nodes, integration_points, shapes, ghost_type);

    auto nb_integration_points = integration_points.cols();
    auto vect_size = nb_integration_points * nb_element;
    auto lmat_size = nb_degree_of_freedom * shapes_size;

    // Extending the shape functions
    /// \todo move this in the shape functions as Voigt format shapes to
    /// have the code in common with the structural elements
    auto shapes_voigt = std::make_unique<Array<Real>>(
        vect_size, lmat_size * nb_degree_of_freedom, 0.);
    auto mshapes_it = shapes_voigt->begin(nb_degree_of_freedom, lmat_size);
    auto shapes_it = shapes.begin(shapes_size);

    for (Int q = 0; q < vect_size; ++q, ++mshapes_it, ++shapes_it) {
      for (Int d = 0; d < nb_degree_of_freedom; ++d) {
        for (Int s = 0; s < shapes_size; ++s) {
          (*mshapes_it)(d, s * nb_degree_of_freedom + d) = (*shapes_it)(s);
        }
      }
    }

    return shapes_voigt;
  }
};

#if defined(AKANTU_STRUCTURAL_MECHANICS)
template <> struct ShapesForMassHelper<_ek_structural> {
  template <ElementType type, class ShapeFunctions>
  static auto getShapes(ShapeFunctions & shape_functions,
                        const Matrix<Real> & integration_points,
                        const Array<Real> & nodes, Int & nb_degree_of_freedom,
                        Int /*nb_element*/, GhostType ghost_type) {
    static_assert(ElementClass<type>::getKind() == _ek_structural,
                  "getShapes for structural elements called with non "
                  "strutral element type");

    auto nb_unknown = ElementClass<type>::getNbStressComponents();
    auto nb_degree_of_freedom_ = ElementClass<type>::getNbDegreeOfFreedom();
    auto nb_nodes_per_element = ElementClass<type>::getNbNodesPerElement();
    auto shapes = std::make_unique<Array<Real>>(
        0, nb_unknown * nb_nodes_per_element * nb_degree_of_freedom_);
    nb_degree_of_freedom = nb_unknown;
    shape_functions.template computeShapesMassOnIntegrationPoints<type>(
        nodes, integration_points, *shapes, ghost_type);

    return shapes;
  }
};
#endif
} // namespace details
} // namespace fe_engine
  //
/* -------------------------------------------------------------------------- */
/**
 * @f$ \tilde{M}_{i} = \sum_j M_{ij} = \sum_j \int \rho \varphi_i \varphi_j dV =
 * \int \rho \varphi_i dV @f$
 */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
template <ElementType type>
void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::assembleFieldMatrix(
    const std::function<void(Matrix<Real> &, const Element &)> & field_funct,
    const ID & matrix_id, const ID & dof_id, DOFManager & dof_manager,
    GhostType ghost_type) const {
  auto nb_degree_of_freedom = dof_manager.getDOFs(dof_id).getNbComponent();
  auto nb_element = mesh.getNbElement(type, ghost_type);

  // \int N * N  so degree 2 * degree of N
  const auto polynomial_degree =
      2 * ElementClassProperty<type>::polynomial_degree;

  // getting the integration points
  auto integration_points =
      integrator.template getIntegrationPoints<type, polynomial_degree>();

  auto nb_integration_points = integration_points.cols();
  auto vect_size = nb_integration_points * nb_element;

  // getting the shapes on the integration points
  auto shapes_voigt =
      fe_engine::details::ShapesForMassHelper<kind>::template getShapes<type>(
          shape_functions, integration_points, mesh.getNodes(),
          nb_degree_of_freedom, nb_element, ghost_type);

  auto lmat_size = shapes_voigt->getNbComponent() / nb_degree_of_freedom;

  // getting the value to assemble on the integration points
  Array<Real> field(vect_size, nb_degree_of_freedom);
  fe_engine::details::fillField(field_funct, field, nb_element,
                                integration_points.cols(), type, ghost_type);

  Array<Real> local_mat(vect_size, lmat_size * lmat_size);

  // computing \rho * N
  for (auto && [mat, N, rho] :
       zip(make_view(local_mat, lmat_size, lmat_size),
           make_view(*shapes_voigt, nb_degree_of_freedom, lmat_size),
           make_view(field, nb_degree_of_freedom))) {
    Matrix<Real> Nt = N.transpose();
    for (Int d = 0; d < Nt.cols(); ++d) {
      Nt(d) *= rho(d);
    }

    mat = Nt * N;
  }

  // integrate the elemental values
  Array<Real> int_field_times_shapes(nb_element, lmat_size * lmat_size,
                                     "inte_rho_x_shapes");
  this->integrator.template integrate<type, polynomial_degree>(
      local_mat, int_field_times_shapes, lmat_size * lmat_size, ghost_type);

  // assemble the elemental values to the matrix
  dof_manager.assembleElementalMatricesToMatrix(
      matrix_id, dof_id, int_field_times_shapes, type, ghost_type);
}

} // namespace akantu

#endif /* AKANTU_FE_ENGINE_TEMPLATE_TMPL_FIELD_HH_ */
