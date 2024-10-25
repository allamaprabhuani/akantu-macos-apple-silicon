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
#include "element_type_map.hh"
#include "mesh_events.hh"
/* -------------------------------------------------------------------------- */
#include <functional>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_FE_ENGINE_HH_
#define AKANTU_FE_ENGINE_HH_

namespace akantu {
class Mesh;
class Integrator;
class ShapeFunctions;
class DOFManager;
class Element;
} // namespace akantu

/* -------------------------------------------------------------------------- */
namespace akantu {
/* -------------------------------------------------------------------------- */

/**
 * The  generic  FEEngine class  derived  in  a  FEEngineTemplate class
 * containing  the
 * shape functions and the integration method
 */
class FEEngine : public MeshEventHandler {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  FEEngine(Mesh & mesh, Int element_dimension = _all_dimensions,
           const ID & id = "fem");

  ~FEEngine() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// pre-compute all the shape functions, their derivatives and the jacobians
  virtual void initShapeFunctions(GhostType ghost_type = _not_ghost) = 0;

  /// extract the nodal values and store them per element
  template <typename T>
  static void
  extractNodalToElementField(const Mesh & mesh, const Array<T> & nodal_f,
                             Array<T> & elemental_f, ElementType type,
                             GhostType ghost_type = _not_ghost,
                             const Array<Int> & filter_elements = empty_filter);

  /// filter a field
  template <typename T>
  static void
  filterElementalData(const Mesh & mesh, const Array<T> & quad_f,
                      Array<T> & filtered_f, ElementType type,
                      GhostType ghost_type = _not_ghost,
                      const Array<Idx> & filter_elements = empty_filter);

  /* ------------------------------------------------------------------------ */
  /* Integration method bridges                                               */
  /* ------------------------------------------------------------------------ */
  /// integrate f for all elements
  virtual void integrate(
      const ElementTypeMapArray<Real> & f, ElementTypeMapArray<Real> & intf,
      const ElementTypeMapArray<Idx> * filter_elements = nullptr) const = 0;

  /// integrate f for all elements of type "type"
  virtual void
  integrate(const Array<Real> & f, Array<Real> & intf, Int nb_degree_of_freedom,
            ElementType type, GhostType ghost_type = _not_ghost,
            const Array<Idx> & filter_elements = empty_filter) const = 0;

  /// integrate a scalar value f on all elements of type "type"
  [[nodiscard]] virtual Real
  integrate(const Array<Real> & f, ElementType type,
            GhostType ghost_type = _not_ghost,
            const Array<Idx> & filter_elements = empty_filter) const = 0;

  /// integrate one element scalar value on all elements of type "type"
  [[nodiscard]] Real integrate(const Ref<const VectorXr> f,
                               const Element & element) const {
    return integrate(f, element.type, element.element, element.ghost_type);
  }

private:
  [[nodiscard]] virtual Real
  integrate(const Ref<const VectorXr> f, ElementType type, Idx index,
            GhostType ghost_type = _not_ghost) const = 0;

  /* ------------------------------------------------------------------------ */
  /* compatibility with old FEEngine fashion */
  /* ------------------------------------------------------------------------ */
public:
  /// get the number of integration points
  [[nodiscard]] virtual Int
  getNbIntegrationPoints(ElementType type,
                         GhostType ghost_type = _not_ghost) const = 0;

  /// get the precomputed shapes
  [[nodiscard]] const virtual Array<Real> &
  getShapes(ElementType type, GhostType ghost_type = _not_ghost,
            Idx id = 0) const = 0;

  /// get the derivatives of shapes
  [[nodiscard]] virtual const Array<Real> &
  getShapesDerivatives(ElementType type, GhostType ghost_type = _not_ghost,
                       Idx id = 0) const = 0;

  /// get integration points
  [[nodiscard]] virtual const MatrixXr &
  getIntegrationPoints(ElementType type,
                       GhostType ghost_type = _not_ghost) const = 0;

  /* ------------------------------------------------------------------------ */
  /* Shape method bridges                                                     */
  /* ------------------------------------------------------------------------ */
  /// Compute the gradient nablauq on the integration points of an element type
  /// from nodal values u
  virtual void gradientOnIntegrationPoints(
      const Array<Real> & u, Array<Real> & nablauq, Int nb_degree_of_freedom,
      ElementType type, GhostType ghost_type = _not_ghost,
      const Array<Idx> & filter_elements = empty_filter) const = 0;

  /// Interpolate a nodal field u at the integration points of an element type
  /// -> uq
  virtual void interpolateOnIntegrationPoints(
      const Array<Real> & u, Array<Real> & uq, Int nb_degree_of_freedom,
      ElementType type, GhostType ghost_type = _not_ghost,
      const Array<Idx> & filter_elements = empty_filter) const = 0;

  /// Interpolate a nodal field u at the integration points of many element
  /// types -> uq
  virtual void interpolateOnIntegrationPoints(
      const Array<Real> & u, ElementTypeMapArray<Real> & uq,
      const ElementTypeMapArray<Idx> * filter_elements = nullptr) const = 0;

  /// pre multiplies a tensor by the shapes derivaties
  virtual void
  computeBtD(const Array<Real> & Ds, Array<Real> & BtDs, ElementType type,
             GhostType ghost_type = _not_ghost,
             const Array<Idx> & filter_elements = empty_filter) const = 0;

  /// left and right  multiplies a tensor by the shapes derivaties
  virtual void
  computeBtDB(const Array<Real> & Ds, Array<Real> & BtDBs, Int order_d,
              ElementType type, GhostType ghost_type = _not_ghost,
              const Array<Idx> & filter_elements = empty_filter) const = 0;

  /// left multiples a vector by the shape functions
  virtual void
  computeNtb(const Array<Real> & bs, Array<Real> & Ntbs, ElementType type,
             GhostType ghost_type = _not_ghost,
             const Array<Idx> & filter_elements = empty_filter) const = 0;

  /// left and right  multiplies a tensor by the shapes
  virtual void
  computeNtbN(const Array<Real> & bs, Array<Real> & NtbNs, ElementType type,
              GhostType ghost_type = _not_ghost,
              const Array<Idx> & filter_elements = empty_filter) const = 0;

  /// Compute the interpolation point position in the global coordinates for
  /// many element types
  virtual void computeIntegrationPointsCoordinates(
      ElementTypeMapArray<Real> & integration_points_coordinates,
      const ElementTypeMapArray<Idx> * filter_elements = nullptr) const = 0;

  /// Compute the interpolation point position in the global coordinates for an
  /// element type
  virtual void computeIntegrationPointsCoordinates(
      Array<Real> & integration_points_coordinates, ElementType type,
      GhostType ghost_type = _not_ghost,
      const Array<Idx> & filter_elements = empty_filter) const = 0;

  /// Build pre-computed matrices for interpolation of field form integration
  /// points at other given positions (interpolation_points)
  virtual void initElementalFieldInterpolationFromIntegrationPoints(
      const ElementTypeMapArray<Real> & interpolation_points_coordinates,
      ElementTypeMapArray<Real> & interpolation_points_coordinates_matrices,
      ElementTypeMapArray<Real> & integration_points_coordinates_inv_matrices,
      const ElementTypeMapArray<Idx> * element_filter) const = 0;

  /// interpolate field at given position (interpolation_points) from given
  /// values of this field at integration points (field)
  virtual void interpolateElementalFieldFromIntegrationPoints(
      const ElementTypeMapArray<Real> & field,
      const ElementTypeMapArray<Real> & interpolation_points_coordinates,
      ElementTypeMapArray<Real> & result, const GhostType ghost_type,
      const ElementTypeMapArray<Idx> * element_filter) const = 0;

  /// Interpolate field at given position from given values of this field at
  /// integration points (field)
  /// using matrices precomputed with
  /// initElementalFieldInterplationFromIntegrationPoints
  virtual void interpolateElementalFieldFromIntegrationPoints(
      const ElementTypeMapArray<Real> & field,
      const ElementTypeMapArray<Real> &
          interpolation_points_coordinates_matrices,
      const ElementTypeMapArray<Real> &
          integration_points_coordinates_inv_matrices,
      ElementTypeMapArray<Real> & result, const GhostType ghost_type,
      const ElementTypeMapArray<Idx> * element_filter) const = 0;

  /// interpolate on a phyiscal point inside an element
  virtual void interpolate(const Ref<const VectorXr> real_coords,
                           const Ref<const MatrixXr> nodal_values,
                           Ref<VectorXr> interpolated,
                           const Element & element) const = 0;

  /// compute the shape on a provided point
  virtual void computeShapes(const Ref<const VectorXr> real_coords, Int elem,
                             ElementType type, Ref<VectorXr> shapes,
                             GhostType ghost_type = _not_ghost) const = 0;

  /// compute the shape derivatives on a provided point
  virtual void
  computeShapeDerivatives(const Ref<const VectorXr> real_coords, Int element,
                          ElementType type, Ref<MatrixXr> shape_derivatives,
                          GhostType ghost_type = _not_ghost) const = 0;

  /// assembles the lumped version of @f[ \int N^t rho N @f]
  virtual void assembleFieldLumped(
      const std::function<void(Matrix<Real> &, const Element &)> & field_funct,
      const ID & matrix_id, const ID & dof_id, DOFManager & dof_manager,
      ElementType type, GhostType ghost_type = _not_ghost) const = 0;

  /// assembles the matrix @f[ \int N^t rho N @f]
  virtual void assembleFieldMatrix(
      const std::function<void(Matrix<Real> &, const Element &)> & field_funct,
      const ID & matrix_id, const ID & dof_id, DOFManager & dof_manager,
      ElementType type, GhostType ghost_type = _not_ghost) const = 0;

  /* ------------------------------------------------------------------------ */
  /* Other methods                                                            */
  /* ------------------------------------------------------------------------ */

  /// pre-compute normals on integration points
  virtual void
  computeNormalsOnIntegrationPoints(GhostType ghost_type = _not_ghost) = 0;

  /// pre-compute normals on integration points
  virtual void
  computeNormalsOnIntegrationPoints(const Array<Real> & /*field*/,
                                    GhostType /*ghost_type*/ = _not_ghost) {
    AKANTU_TO_IMPLEMENT();
  }

  /// pre-compute normals on integration points
  virtual void computeNormalsOnIntegrationPoints(
      const Array<Real> & /*field*/, Array<Real> & /*normal*/,
      ElementType /*type*/, GhostType /*ghost_type*/ = _not_ghost) const {
    AKANTU_TO_IMPLEMENT();
  }

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

private:
  /// initialise the class
  void init();

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  using ElementTypesIteratorHelper =
      ElementTypeMapArray<Real, ElementType>::ElementTypesIteratorHelper;

  [[nodiscard]] ElementTypesIteratorHelper
  elementTypes(Int dim = _all_dimensions, GhostType ghost_type = _not_ghost,
               ElementKind kind = _ek_regular) const;

  /// get the dimension of the element handeled by this fe_engine object
  AKANTU_GET_MACRO_AUTO(ElementDimension, element_dimension);

  /// get the mesh contained in the fem object
  AKANTU_GET_MACRO_AUTO(Mesh, mesh);
  /// get the mesh contained in the fem object
  AKANTU_GET_MACRO_NOT_CONST(Mesh, mesh, Mesh &);

  /// get the in-radius of an element
  template <class Derived>
  [[nodiscard]] static inline Real
  getElementInradius(const Eigen::MatrixBase<Derived> & coord,
                     ElementType type);

  [[nodiscard]] inline Real getElementInradius(const Element & element) const;

  /// get the normals on integration points
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(NormalsOnIntegrationPoints,
                                         normals_on_integration_points, Real);

  /// get cohesive element type for a given facet type
  [[nodiscard]] static inline constexpr ElementType
  getCohesiveElementType(ElementType type_facet);

  /// get igfem element type for a given regular type
  [[nodiscard]] static inline Vector<ElementType>
  getIGFEMElementTypes(ElementType type);

  /// get the interpolation element associated to an element type
  [[nodiscard]] static inline constexpr auto
  getInterpolationType(ElementType el_type);

  /// get the shape function class (probably useless: see getShapeFunction in
  /// fe_engine_template.hh)
  [[nodiscard]] virtual const ShapeFunctions &
  getShapeFunctionsInterface() const = 0;
  /// get the integrator class (probably useless: see getIntegrator in
  /// fe_engine_template.hh)
  [[nodiscard]] virtual const Integrator & getIntegratorInterface() const = 0;

  AKANTU_GET_MACRO(ID, id, ID);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  ID id;

  /// spatial dimension of the problem
  Int element_dimension;

  /// the mesh on which all computation are made
  Mesh & mesh;

  /// normals at integration points
  ElementTypeMapArray<Real> normals_on_integration_points;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

/// standard output stream operator
inline std::ostream & operator<<(std::ostream & stream,
                                 const FEEngine & _this) {
  _this.printself(stream);
  return stream;
}

} // namespace akantu

#include "fe_engine_inline_impl.hh"
#include "fe_engine_template.hh"

#endif /* AKANTU_FE_ENGINE_HH_ */
