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
/* -------------------------------------------------------------------------- */
#include <type_traits>
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
#ifndef AKANTU_FE_ENGINE_TEMPLATE_HH_
#define AKANTU_FE_ENGINE_TEMPLATE_HH_

namespace akantu {
class Integrator;
class ShapeFunctions;
} // namespace akantu

namespace akantu {
class DOFManager;
namespace fe_engine {
  namespace details {
    template <ElementKind> struct AssembleLumpedTemplateHelper;
    template <ElementKind> struct AssembleFieldMatrixHelper;
  } // namespace details
} // namespace fe_engine

template <ElementKind, typename> struct AssembleFieldMatrixStructHelper;

struct DefaultIntegrationOrderFunctor {
  template <ElementType type> static inline constexpr int getOrder() {
    return ElementClassProperty<type>::polynomial_degree;
  }
};

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind = _ek_regular,
          class IntegrationOrderFunctor = DefaultIntegrationOrderFunctor>
class FEEngineTemplate : public FEEngine {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  using Integ = I<kind, IntegrationOrderFunctor>;
  using Shape = S<kind>;

  FEEngineTemplate(Mesh & mesh, Int spatial_dimension = _all_dimensions,
                   const ID & id = "fem", bool do_not_precompute = false);

  ~FEEngineTemplate() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// pre-compute all the shape functions, their derivatives and the jacobians
  void initShapeFunctions(GhostType ghost_type = _not_ghost) override;
  void initShapeFunctions(const Array<Real> & nodes,
                          GhostType ghost_type = _not_ghost);

  /* ------------------------------------------------------------------------ */
  /* Integration method bridges                                               */
  /* ------------------------------------------------------------------------ */
  /// integrate f for all elements
  void
  integrate(const ElementTypeMapArray<Real> & f,
            ElementTypeMapArray<Real> & intf,
            const ElementTypeMapArray<Idx> * filter_elements = nullptr) const;

  /// integrate f for all elements of type "type"
  void
  integrate(const Array<Real> & f, Array<Real> & intf, Int nb_degree_of_freedom,
            ElementType type, GhostType ghost_type = _not_ghost,
            const Array<Idx> & filter_elements = empty_filter) const override;

  /// integrate a scalar value on all elements of type "type"
  [[nodiscard]] Real
  integrate(const Array<Real> & f, ElementType type,
            GhostType ghost_type = _not_ghost,
            const Array<Idx> & filter_elements = empty_filter) const override;

  /// integrate one element scalar value on all elements of type "type"
  [[nodiscard]] Real
  integrate(const Ref<const VectorXr> f, ElementType type, Int index,
            GhostType ghost_type = _not_ghost) const override;

private:
  template <ElementKind kind_ = kind, typename D1, typename D2, typename D3,
            std::enable_if_t<aka::are_vectors<D1, D3>::value and
                             kind_ == _ek_regular> * = nullptr>
  inline void interpolateImpl(const Eigen::MatrixBase<D1> & real_coords,
                              const Eigen::MatrixBase<D2> & nodal_values,
                              Eigen::MatrixBase<D3> & interpolated,
                              const Element & element) const;

  template <ElementKind kind_ = kind, typename D1, typename D2, typename D3,
            std::enable_if_t<aka::are_vectors<D1, D3>::value and
                             kind_ != _ek_regular> * = nullptr>
  inline void interpolateImpl(const Eigen::MatrixBase<D1> & /*real_coords*/,
                              const Eigen::MatrixBase<D2> & /*nodal_values*/,
                              Eigen::MatrixBase<D3> & /*interpolated*/,
                              const Element & /*element*/) const {
    AKANTU_TO_IMPLEMENT();
  }

public:
  /// interpolate on a phyiscal point inside an element
  void interpolate(const Ref<const VectorXr> real_coords,
                   const Ref<const MatrixXr> nodal_values,
                   Ref<VectorXr> interpolated,
                   const Element & element) const override;

  /// get the number of integration points
  [[nodiscard]] Int
  getNbIntegrationPoints(ElementType type,
                         GhostType ghost_type = _not_ghost) const override;

  /// get shapes precomputed
  [[nodiscard]] const Array<Real> & getShapes(ElementType type,
                                              GhostType ghost_type = _not_ghost,
                                              Int id = 0) const override;

  /// get the derivatives of shapes
  [[nodiscard]] const Array<Real> &
  getShapesDerivatives(ElementType type, GhostType ghost_type = _not_ghost,
                       Int id = 0) const override;

  /// get integration points
  [[nodiscard]] inline const Matrix<Real> &
  getIntegrationPoints(ElementType type,
                       GhostType ghost_type = _not_ghost) const override;

  /* ------------------------------------------------------------------------ */
  /* Shape method bridges                                                     */
  /* ------------------------------------------------------------------------ */

  /// compute the gradient of a nodal field on the integration points
  void gradientOnIntegrationPoints(
      const Array<Real> & u, Array<Real> & nablauq, Int nb_degree_of_freedom,
      ElementType type, GhostType ghost_type = _not_ghost,
      const Array<Idx> & filter_elements = empty_filter) const override;

  /// interpolate a nodal field on the integration points
  void interpolateOnIntegrationPoints(
      const Array<Real> & u, Array<Real> & uq, Int nb_degree_of_freedom,
      ElementType type, GhostType ghost_type = _not_ghost,
      const Array<Idx> & filter_elements = empty_filter) const override;

  /// interpolate a nodal field on the integration points given a
  /// by_element_type
  void interpolateOnIntegrationPoints(
      const Array<Real> & u, ElementTypeMapArray<Real> & uq,
      const ElementTypeMapArray<Idx> * filter_elements =
          nullptr) const override;

  /// pre multiplies a tensor by the shapes derivaties
  void
  computeBtD(const Array<Real> & Ds, Array<Real> & BtDs, ElementType type,
             GhostType ghost_type,
             const Array<Idx> & filter_elements = empty_filter) const override;

  /// left and right  multiplies a tensor by the shapes derivaties
  void
  computeBtDB(const Array<Real> & Ds, Array<Real> & BtDBs, Int order_d,
              ElementType type, GhostType ghost_type,
              const Array<Idx> & filter_elements = empty_filter) const override;

  /// left multiples a vector by the shape functions
  void computeNtb(const Array<Real> & bs, Array<Real> & Ntbs, ElementType type,
                  GhostType ghost_type,
                  const Array<Idx> & filter_elements) const override;

  /// left and right  multiplies a tensor by the shapes
  void
  computeNtbN(const Array<Real> & bs, Array<Real> & NtbNs, ElementType type,
              GhostType ghost_type,
              const Array<Idx> & filter_elements = empty_filter) const override;

  /// compute the position of integration points given by an element_type_map
  /// from nodes position
  inline void computeIntegrationPointsCoordinates(
      ElementTypeMapArray<Real> & quadrature_points_coordinates,
      const ElementTypeMapArray<Idx> * filter_elements =
          nullptr) const override;

  /// compute the position of integration points from nodes position
  inline void computeIntegrationPointsCoordinates(
      Array<Real> & quadrature_points_coordinates, ElementType type,
      GhostType ghost_type = _not_ghost,
      const Array<Idx> & filter_elements = empty_filter) const override;

  /// interpolate field at given position (interpolation_points) from given
  /// values of this field at integration points (field)
  inline void interpolateElementalFieldFromIntegrationPoints(
      const ElementTypeMapArray<Real> & field,
      const ElementTypeMapArray<Real> & interpolation_points_coordinates,
      ElementTypeMapArray<Real> & result, GhostType ghost_type,
      const ElementTypeMapArray<Idx> * element_filter) const override;

  /// Interpolate field at given position from given values of this field at
  /// integration points (field)
  /// using matrices precomputed with
  /// initElementalFieldInterplationFromIntegrationPoints
  inline void interpolateElementalFieldFromIntegrationPoints(
      const ElementTypeMapArray<Real> & field,
      const ElementTypeMapArray<Real> &
          interpolation_points_coordinates_matrices,
      const ElementTypeMapArray<Real> & quad_points_coordinates_inv_matrices,
      ElementTypeMapArray<Real> & result, GhostType ghost_type,
      const ElementTypeMapArray<Idx> * element_filter) const override;

  /// Build pre-computed matrices for interpolation of field form integration
  /// points at other given positions (interpolation_points)
  inline void initElementalFieldInterpolationFromIntegrationPoints(
      const ElementTypeMapArray<Real> & interpolation_points_coordinates,
      ElementTypeMapArray<Real> & interpolation_points_coordinates_matrices,
      ElementTypeMapArray<Real> & quad_points_coordinates_inv_matrices,
      const ElementTypeMapArray<Idx> * element_filter = nullptr) const override;

  /// find natural coords from real coords provided an element
  void inverseMap(const Ref<const VectorXr> real_coords, Int element,
                  ElementType type, Ref<VectorXr> natural_coords,
                  GhostType ghost_type = _not_ghost) const;

  /// return true if the coordinates provided are inside the element, false
  /// otherwise
  inline bool contains(const Vector<Real> & real_coords, Int element,
                       ElementType type,
                       GhostType ghost_type = _not_ghost) const;

private:
  template <ElementKind kind_ = kind, typename D1, typename D2,
            std::enable_if_t<aka::are_vectors<D1, D2>::value and
                             kind_ != _ek_cohesive> * = nullptr>
  inline void computeShapesImpl(const Eigen::MatrixBase<D1> & real_coords,
                                Int element, ElementType type,
                                Eigen::MatrixBase<D2> & shapes,
                                GhostType ghost_type = _not_ghost) const;

  template <ElementKind kind_ = kind, typename D1, typename D2,
            std::enable_if_t<aka::are_vectors<D1, D2>::value and
                             kind_ == _ek_cohesive> * = nullptr>
  inline void computeShapesImpl(const Eigen::MatrixBase<D1> & /*real_coords*/,
                                Int /*element*/, ElementType /*type*/,
                                Eigen::MatrixBase<D2> & /*shapes*/,
                                GhostType /*ghost_type*/ = _not_ghost) const {
    AKANTU_TO_IMPLEMENT();
  }

  template <ElementKind kind_ = kind, typename D1, typename D2,
            std::enable_if_t<aka::is_vector_v<D1> and kind_ != _ek_cohesive> * =
                nullptr>
  inline void
  computeShapeDerivativesImpl(const Eigen::MatrixBase<D1> & real_coords,
                              Int element, ElementType type,
                              Eigen::MatrixBase<D2> & shape_derivatives,
                              GhostType ghost_type = _not_ghost) const;

  template <ElementKind kind_ = kind, typename D1, typename D2,
            std::enable_if_t<aka::is_vector_v<D1> and kind_ == _ek_cohesive> * =
                nullptr>
  inline void
  computeShapeDerivativesImpl(const Eigen::MatrixBase<D1> & /*real_coords*/,
                              Int /*element*/, ElementType /*type*/,
                              Eigen::MatrixBase<D2> & /*shape_derivatives*/,
                              GhostType /*ghost_type*/ = _not_ghost) const {
    AKANTU_TO_IMPLEMENT();
  }

public:
  /// compute the shape on a provided point
  inline void computeShapes(const Ref<const VectorXr> real_coords, Int element,
                            ElementType type, Ref<VectorXr> shapes,
                            GhostType ghost_type = _not_ghost) const override {
    this->template computeShapesImpl(real_coords, element, type, shapes,
                                     ghost_type);
  }

  /// compute the shape derivatives on a provided point
  inline void
  computeShapeDerivatives(const Ref<const VectorXr> real_coords, Int element,
                          ElementType type, Ref<MatrixXr> shape_derivatives,
                          GhostType ghost_type = _not_ghost) const override {
    this->template computeShapeDerivativesImpl<kind>(
        real_coords, element, type, shape_derivatives, ghost_type);
  }

  /* ------------------------------------------------------------------------ */
  /* Other methods                                                            */
  /* ------------------------------------------------------------------------ */
  /// pre-compute normals on integration points
  void
  computeNormalsOnIntegrationPoints(GhostType ghost_type = _not_ghost) override;
  void
  computeNormalsOnIntegrationPoints(const Array<Real> & field,
                                    GhostType ghost_type = _not_ghost) override;
  void computeNormalsOnIntegrationPoints(
      const Array<Real> & field, Array<Real> & normal, ElementType type,
      GhostType ghost_type = _not_ghost) const override;

  template <ElementType type, ElementKind kind_ = kind,
            std::enable_if_t<kind_ != _ek_regular> * = nullptr>
  void computeNormalsOnIntegrationPoints(const Array<Real> & /*field*/,
                                         Array<Real> & /*normal*/,
                                         GhostType /*ghost_type*/) const {
    AKANTU_TO_IMPLEMENT();
  }

  template <
      ElementType type, ElementKind kind_ = kind,
      std::enable_if_t<kind_ == _ek_regular and type != _point_1> * = nullptr>
  void computeNormalsOnIntegrationPoints(const Array<Real> & field,
                                         Array<Real> & normal,
                                         GhostType ghost_type) const;

  template <
      ElementType type, ElementKind kind_ = kind,
      std::enable_if_t<kind_ == _ek_regular and type == _point_1> * = nullptr>
  void computeNormalsOnIntegrationPoints(const Array<Real> & field,
                                         Array<Real> & normal,
                                         GhostType ghost_type) const;

public:
  /// function to print the contain of the class
  void printself(std::ostream & stream, int indent = 0) const override;

  void assembleFieldLumped(
      const std::function<void(Matrix<Real> &, const Element &)> & field_funct,
      const ID & matrix_id, const ID & dof_id, DOFManager & dof_manager,
      ElementType type, GhostType ghost_type) const override;

private:
  template <ElementKind kind_ = kind,
            std::enable_if_t<kind_ != _ek_cohesive> * = nullptr>
  void assembleFieldMatrixImpl(
      const std::function<void(Matrix<Real> &, const Element &)> & field_funct,
      const ID & matrix_id, const ID & dof_id, DOFManager & dof_manager,
      ElementType type, GhostType ghost_type) const;

  template <ElementKind kind_ = kind,
            std::enable_if_t<kind_ == _ek_cohesive> * = nullptr>
  void assembleFieldMatrixImpl(
      const std::function<void(Matrix<Real> &, const Element &)> & field_funct,
      const ID & matrix_id, const ID & dof_id, DOFManager & dof_manager,
      ElementType type, GhostType ghost_type) const;

public:
  /// assemble a field as a matrix (ex. rho to mass matrix)
  void assembleFieldMatrix(
      const std::function<void(Matrix<Real> &, const Element &)> & field_funct,
      const ID & matrix_id, const ID & dof_id, DOFManager & dof_manager,
      ElementType type, GhostType ghost_type) const override {
    this->assembleFieldMatrixImpl(field_funct, matrix_id, dof_id, dof_manager,
                                  type, ghost_type);
  }

private:
  friend struct fe_engine::details::AssembleLumpedTemplateHelper<kind>;
  friend struct fe_engine::details::AssembleFieldMatrixHelper<kind>;
  friend struct AssembleFieldMatrixStructHelper<kind, void>;

  /// templated function to compute the scaling to assemble a lumped matrix
  template <ElementType type>
  void assembleFieldLumped(
      const std::function<void(Matrix<Real> &, const Element &)> & field_funct,
      const ID & matrix_id, const ID & dof_id, DOFManager & dof_manager,
      GhostType ghost_type) const;

  /// @f$ \tilde{M}_{i} = \sum_j M_{ij} = \sum_j \int \rho \varphi_i \varphi_j
  /// dV = \int \rho \varphi_i dV @f$
  template <ElementType type>
  void assembleLumpedRowSum(const Array<Real> & field, const ID & matrix_id,
                            const ID & dof_id, DOFManager & dof_manager,
                            GhostType ghost_type) const;

  /// @f$ \tilde{M}_{i} = c * M_{ii} = \int_{V_e} \rho dV @f$
  template <ElementType type>
  void assembleLumpedDiagonalScaling(const Array<Real> & field,
                                     const ID & matrix_id, const ID & dof_id,
                                     DOFManager & dof_manager,
                                     GhostType ghost_type) const;

  /// assemble a field as a matrix (ex. rho to mass matrix)
  template <ElementType type>
  void assembleFieldMatrix(
      const std::function<void(Matrix<Real> &, const Element &)> & field_funct,
      const ID & matrix_id, const ID & dof_id, DOFManager & dof_manager,
      GhostType ghost_type) const;

#ifdef AKANTU_STRUCTURAL_MECHANICS

  /// assemble a field as a matrix for structural elements (ex. rho to mass
  /// matrix)
  template <ElementType type>
  void assembleFieldMatrix(const Array<Real> & field_1,
                           Int nb_degree_of_freedom, SparseMatrix & M,
                           Array<Real> * n,
                           ElementTypeMapArray<Real> & rotation_mat,
                           GhostType ghost_type) const;
#endif

  /* ------------------------------------------------------------------------ */
  /* Mesh Event Handler interface                                             */
  /* ------------------------------------------------------------------------ */
public:
  void onElementsAdded(const Array<Element> & /*new_elements*/,
                       const NewElementsEvent & /*unused*/) override;
  void onElementsRemoved(const Array<Element> & /*unused*/,
                         const ElementTypeMapArray<Idx> & /*unused*/,
                         const RemovedElementsEvent & /*unused*/) override;
  void onElementsChanged(const Array<Element> & /*unused*/,
                         const Array<Element> & /*unused*/,
                         const ElementTypeMapArray<Idx> & /*unused*/,
                         const ChangedElementsEvent & /*unused*/) override;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get the shape class (probably useless: see getShapeFunction)
  const ShapeFunctions & getShapeFunctionsInterface() const override {
    return shape_functions;
  };
  /// get the shape class
  const Shape & getShapeFunctions() const { return shape_functions; };

  /// get the integrator class (probably useless: see getIntegrator)
  const Integrator & getIntegratorInterface() const override {
    return integrator;
  };
  /// get the integrator class
  const Integ & getIntegrator() const { return integrator; };

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  Integ integrator;
  Shape shape_functions;
};

} // namespace akantu

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
#include "fe_engine_template_tmpl.hh"
#include "fe_engine_template_tmpl_field.hh"
/* -------------------------------------------------------------------------- */
/* Shape Linked specialization                                                */
/* -------------------------------------------------------------------------- */
#if defined(AKANTU_STRUCTURAL_MECHANICS)
#include "fe_engine_template_tmpl_struct.hh"
#endif
/* -------------------------------------------------------------------------- */
/* Shape IGFEM specialization                                                 */
/* -------------------------------------------------------------------------- */
#if defined(AKANTU_IGFEM)
#include "fe_engine_template_tmpl_igfem.hh"
#endif

#endif /* AKANTU_FE_ENGINE_TEMPLATE_HH_ */
