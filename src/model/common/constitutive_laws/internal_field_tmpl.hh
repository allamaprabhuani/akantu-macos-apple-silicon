/**
 * Copyright (©) 2013-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "constitutive_law.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_INTERNAL_FIELD_TMPL_HH_
#define AKANTU_INTERNAL_FIELD_TMPL_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
template <typename T>
InternalField<T>::InternalField(
    const ID & id, ConstitutiveLawInternalHandler & constitutive_law, Int dim,
    const ID & fem_id, const ElementTypeMapArray<Idx> & element_filter)
    : InternalFieldBase(id),
      ElementTypeMapArray<T>(id, constitutive_law.getID()),
      constitutive_law(constitutive_law),
      fem(constitutive_law.getFEEngine(fem_id)), element_filter(element_filter),
      spatial_dimension(dim) {}

/* -------------------------------------------------------------------------- */
template <typename T>
InternalField<T>::InternalField(
    const ID & id, ConstitutiveLawInternalHandler & constitutive_law)
    : InternalField(id, constitutive_law,
                    constitutive_law.getSpatialDimension(), "",
                    constitutive_law.getElementFilter()) {}

/* -------------------------------------------------------------------------- */
template <typename T>
InternalField<T>::InternalField(
    const ID & id, ConstitutiveLawInternalHandler & constitutive_law,
    const ID & fem_id, const ElementTypeMapArray<Idx> & element_filter)
    : InternalField(id, constitutive_law,
                    constitutive_law.getSpatialDimension(), fem_id,
                    element_filter) {}

/* -------------------------------------------------------------------------- */
template <typename T>
InternalField<T>::InternalField(const ID & id, const InternalField<T> & other)
    : InternalFieldBase(id),
      ElementTypeMapArray<T>(id, other.constitutive_law.getID()),
      constitutive_law(other.constitutive_law), fem(other.fem),
      element_filter(other.element_filter), default_value(other.default_value),
      spatial_dimension(other.spatial_dimension),
      element_kind(other.element_kind), nb_component(other.nb_component) {
  AKANTU_DEBUG_ASSERT(other.is_init,
                      "Cannot create a copy of a non initialized field");
  this->internalInitialize(this->nb_component);
}

/* -------------------------------------------------------------------------- */
template <typename T>
void InternalField<T>::setElementKind(ElementKind element_kind) {
  this->element_kind = element_kind;
}

/* -------------------------------------------------------------------------- */
template <typename T> void InternalField<T>::initialize(Int nb_component) {
  internalInitialize(nb_component);
}

/* -------------------------------------------------------------------------- */
template <typename T> void InternalField<T>::initializeHistory() {
  if (not previous_values) {
    previous_values = std::shared_ptr<InternalField<T>>(
        new InternalField<T>("previous_" + this->getID(), *this));
  }
}

/* -------------------------------------------------------------------------- */
template <typename T> void InternalField<T>::resize() {
  if (not this->is_init) {
    return;
  }

  ElementTypeMap<Int> old_sizes;
  for (auto ghost_type : ghost_types) {
    for (const auto & type : this->filterTypes(ghost_type)) {
      if (this->exists(type, ghost_type)) {
        old_sizes(type, ghost_type) = this->operator()(type, ghost_type).size();
      } else {
        old_sizes(type, ghost_type) = 0;
      }
    }
  }

  ElementTypeMapArray<T>::initialize(
      fem, _element_filter = &element_filter, _element_kind = element_kind,
      _nb_component = nb_component, _with_nb_element = true,
      _do_not_default = true);

  for (auto ghost_type : ghost_types) {
    for (const auto & type : this->elementTypes(ghost_type)) {
      auto & vect = this->operator()(type, ghost_type);
      auto old_size = old_sizes(type, ghost_type);
      auto new_size = vect.size();
      this->setArrayValues(vect.data() + old_size * vect.getNbComponent(),
                           vect.data() + new_size * vect.getNbComponent());

      this->releases(type, ghost_type) += 1;
    }
  }

  if (this->previous_values) {
    this->previous_values->resize();
  }
}

/* -------------------------------------------------------------------------- */
template <typename T> void InternalField<T>::setDefaultValue(const T & value) {
  this->default_value = value;
  this->reset();
}

/* -------------------------------------------------------------------------- */
template <typename T> void InternalField<T>::reset() {
  for (auto ghost_type : ghost_types) {
    for (const auto & type : this->elementTypes(ghost_type)) {
      auto & vect = (*this)(type, ghost_type);
      this->setArrayValues(vect.data(),
                           vect.data() + vect.size() * vect.getNbComponent());
      this->releases(type, ghost_type) += 1;
    }
  }
}

/* -------------------------------------------------------------------------- */
template <typename T>
void InternalField<T>::internalInitialize(Int nb_component) {
  if (not this->is_init) {
    this->nb_component = nb_component;
    this->is_init = true;
  } else {
    resize();
    return;
  }

  for (auto ghost_type : ghost_types) {
    for (const auto & type : this->filterTypes(ghost_type)) {
      this->releases(type, ghost_type) = -1;
    }
  }

  ElementTypeMapArray<T>::initialize(
      fem, _element_filter = &element_filter, _element_kind = element_kind,
      _nb_component = nb_component, _with_nb_element = true,
      _do_not_default = true);

  this->reset();

  if (this->previous_values) {
    this->previous_values->internalInitialize(nb_component);
  }
}

/* -------------------------------------------------------------------------- */
template <typename T>
void InternalField<T>::setArrayValues(T * begin, T * end) {
  for (; begin < end; ++begin) {
    *begin = this->default_value;
  }
}

/* -------------------------------------------------------------------------- */
template <typename T> void InternalField<T>::saveCurrentValues() {
  AKANTU_DEBUG_ASSERT(this->previous_values != nullptr,
                      "The history of the internal "
                          << this->getID() << " has not been activated");
  if (not this->is_init) {
    return;
  }

  for (auto ghost_type : ghost_types) {
    for (const auto & type : this->elementTypes(ghost_type)) {
      (*this->previous_values)(type, ghost_type)
          .copy((*this)(type, ghost_type));
    }
  }
}

/* -------------------------------------------------------------------------- */
template <typename T> void InternalField<T>::restorePreviousValues() {
  AKANTU_DEBUG_ASSERT(this->previous_values != nullptr,
                      "The history of the internal "
                          << this->getID() << " has not been activated");

  if (not this->is_init) {
    return;
  }

  for (auto ghost_type : ghost_types) {
    for (const auto & type : this->elementTypes(ghost_type)) {
      (*this)(type, ghost_type)
          .copy((*this->previous_values)(type, ghost_type));
    }
  }
}

/* -------------------------------------------------------------------------- */
template <typename T>
void InternalField<T>::removeIntegrationPoints(
    const ElementTypeMapArray<Idx> & new_numbering) {
  for (auto ghost_type : ghost_types) {
    for (auto type : new_numbering.elementTypes(_all_dimensions, ghost_type,
                                                _ek_not_defined)) {
      if (not this->exists(type, ghost_type)) {
        continue;
      }

      Array<T> & vect = (*this)(type, ghost_type);
      if (vect.empty()) {
        continue;
      }

      const auto & renumbering = new_numbering(type, ghost_type);

      auto nb_quad_per_elem = fem.getNbIntegrationPoints(type, ghost_type);
      auto nb_component = vect.getNbComponent();

      Array<T> tmp(renumbering.size() * nb_quad_per_elem, nb_component);

      AKANTU_DEBUG_ASSERT(
          tmp.size() == vect.size(),
          "Something strange append some mater was created from nowhere!!");

      AKANTU_DEBUG_ASSERT(
          tmp.size() == vect.size(),
          "Something strange append some mater was created or disappeared in "
              << vect.getID() << "(" << vect.size() << "!=" << tmp.size()
              << ") "
                 "!!");

      Int new_size = 0;
      for (Int i = 0; i < renumbering.size(); ++i) {
        auto new_i = renumbering(i);
        if (new_i != Int(-1)) {
          memcpy(tmp.data() + new_i * nb_component * nb_quad_per_elem,
                 vect.data() + i * nb_component * nb_quad_per_elem,
                 nb_component * nb_quad_per_elem * sizeof(T));
          ++new_size;
        }
      }
      tmp.resize(new_size * nb_quad_per_elem);
      vect.copy(tmp);
    }
  }
}

/* -------------------------------------------------------------------------- */
template <typename T>
void InternalField<T>::printself(std::ostream & stream,
                                 int indent [[gnu::unused]]) const {
  stream << "InternalField [ " << this->getID();
#if !defined(AKANTU_NDEBUG)
  if (AKANTU_DEBUG_TEST(dblDump)) {
    stream << "\n";
    ElementTypeMapArray<T>::printself(stream, indent + 3);
  } else {
#endif
    stream << " {" << this->getData(_not_ghost).size() << " types - "
           << this->getData(_ghost).size() << " ghost types"
           << "}";
#if !defined(AKANTU_NDEBUG)
  }
#endif
  stream << " ]";
}

/* -------------------------------------------------------------------------- */
template <>
inline void
ParameterTyped<InternalField<Real>>::setAuto(const ParserParameter & in_param) {
  Parameter::setAuto(in_param);
  Real r = in_param;
  param.setDefaultValue(r);
}
/* -------------------------------------------------------------------------- */

template <>
inline void ParameterTyped<InternalField<Real>>::set(std::any value) {
  param.setDefaultValue(std::any_cast<Real>(value));
}

/* -------------------------------------------------------------------------- */
template <typename T> inline InternalField<T>::operator T() const {
  return default_value;
}

} // namespace akantu

#endif /* AKANTU_INTERNAL_FIELD_TMPL_HH_ */
