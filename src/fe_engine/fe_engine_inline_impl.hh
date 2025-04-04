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
#include "element_class.hh"
// #include "fe_engine.hh"
#include "mesh.hh"
/* -------------------------------------------------------------------------- */
#include "element_type_conversion.hh"
/* -------------------------------------------------------------------------- */

// #ifndef __AKANTU_FE_ENGINE_INLINE_IMPL_CC__
// #define __AKANTU_FE_ENGINE_INLINE_IMPL_CC__

namespace akantu {

/* -------------------------------------------------------------------------- */
template <class Derived>
inline Real
FEEngine::getElementInradius(const Eigen::MatrixBase<Derived> & coord,
                             ElementType type) {
  return tuple_dispatch<AllElementTypes>(
      [&coord](auto type) -> Real {
        return ElementClass<type.value>::getInradius(coord);
      },
      type);
}

/* -------------------------------------------------------------------------- */
inline Real FEEngine::getElementInradius(const Element & element) const {
  auto spatial_dimension = mesh.getSpatialDimension();
  auto positions = make_view(mesh.getNodes(), spatial_dimension).begin();
  auto connectivity = mesh.getConnectivities().get(element);

  Matrix<Real> coords(spatial_dimension, connectivity.size());

  for (auto && data : zip(connectivity, coords)) {
    std::get<1>(data) = positions[std::get<0>(data)];
  }

  return getElementInradius(coords, element.type);
}

/* -------------------------------------------------------------------------- */
inline constexpr auto FEEngine::getInterpolationType(ElementType type) {
  return convertType<ElementType, InterpolationType>(type);
}

/* -------------------------------------------------------------------------- */
/// @todo rewrite this function in order to get the cohesive element
/// type directly from the facet
#if defined(AKANTU_COHESIVE_ELEMENT)
inline constexpr ElementType
FEEngine::getCohesiveElementType(ElementType type) {
  return tuple_dispatch_with_default<ElementTypes_t<_ek_regular>>(
      [](auto type) -> ElementType {
        return CohesiveFacetProperty<type.value>::cohesive_type;
      },
      [](auto /*type*/) { return _not_defined; }, type);
}
#endif

/* -------------------------------------------------------------------------- */
#if defined(AKANTU_IGFEM)
} // akantu
#include "igfem_helper.hh"
namespace akantu {

inline Vector<ElementType> FEEngine::getIGFEMElementTypes(ElementType type) {
  tuple_dispatch<ElementTypes_t<_ek_regular>>(
      [](auto type) { return IGFEMHelper::getIGFEMElementTypes<type.value>(); },
      type);
}
#endif

/* -------------------------------------------------------------------------- */
template <typename T>
void FEEngine::extractNodalToElementField(const Mesh & mesh,
                                          const Array<T> & nodal_f,
                                          Array<T> & elemental_f,
                                          ElementType type,
                                          GhostType ghost_type,
                                          const Array<Int> & filter_elements) {
  AKANTU_DEBUG_IN();

  auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  auto nb_degree_of_freedom = nodal_f.getNbComponent();
  auto nb_element = mesh.getNbElement(type, ghost_type);
  auto * conn_val = mesh.getConnectivity(type, ghost_type).data();

  if (filter_elements != empty_filter) {
    nb_element = filter_elements.size();
  }

  elemental_f.resize(nb_element);

  const T * nodal_f_val = nodal_f.data();
  T * f_val = elemental_f.data();

  const Idx * el_conn;
  for (Int el = 0; el < nb_element; ++el) {
    if (filter_elements != empty_filter) {
      el_conn = conn_val + filter_elements(el) * nb_nodes_per_element;
    } else {
      el_conn = conn_val + el * nb_nodes_per_element;
    }

    for (Int n = 0; n < nb_nodes_per_element; ++n) {
      auto node = *(el_conn + n);
      std::copy(nodal_f_val + node * nb_degree_of_freedom,
                nodal_f_val + (node + 1) * nb_degree_of_freedom, f_val);
      f_val += nb_degree_of_freedom;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <typename T>
void FEEngine::filterElementalData(const Mesh & mesh, const Array<T> & elem_f,
                                   Array<T> & filtered_f, ElementType type,
                                   GhostType ghost_type,
                                   const Array<Int> & filter_elements) {
  AKANTU_DEBUG_IN();

  auto nb_element = mesh.getNbElement(type, ghost_type);
  if (nb_element == 0) {
    filtered_f.resize(0);
    return;
  }

  auto nb_degree_of_freedom = elem_f.getNbComponent();
  auto nb_data_per_element = elem_f.size() / nb_element;

  if (filter_elements != empty_filter) {
    nb_element = filter_elements.size();
  }

  filtered_f.resize(nb_element * nb_data_per_element);

  const T * elem_f_val = elem_f.data();
  T * f_val = filtered_f.data();

  UInt el_offset;
  for (Idx el = 0; el < nb_element; ++el) {
    if (filter_elements != empty_filter) {
      el_offset = filter_elements(el);
    } else {
      el_offset = el;
    }

    std::copy(elem_f_val +
                  el_offset * nb_data_per_element * nb_degree_of_freedom,
              elem_f_val +
                  (el_offset + 1) * nb_data_per_element * nb_degree_of_freedom,
              f_val);
    f_val += nb_degree_of_freedom * nb_data_per_element;
  }

  AKANTU_DEBUG_OUT();
}

} // namespace akantu

// #endif /* __AKANTU_FE_ENGINE_INLINE_IMPL_CC__ */
