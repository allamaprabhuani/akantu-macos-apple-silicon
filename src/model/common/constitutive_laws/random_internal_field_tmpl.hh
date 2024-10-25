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
#include "aka_common.hh"
#include "aka_random_generator.hh"
#include "internal_field_tmpl.hh"

/* -------------------------------------------------------------------------- */

#ifndef AKANTU_RANDOM_INTERNAL_FIELD_TMPL_HH_
#define AKANTU_RANDOM_INTERNAL_FIELD_TMPL_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
template <typename T, template <typename> class BaseField,
          template <typename> class Generator>
void RandomInternalField<T, BaseField, Generator>::initialize(
    Int nb_component) {
  this->internalInitialize(nb_component);
}

/* ------------------------------------------------------------------------ */
template <typename T, template <typename> class BaseField,
          template <typename> class Generator>
void RandomInternalField<T, BaseField, Generator>::setDefaultValue(
    const T & value) {
  random_parameter.setBaseValue(value);
  this->reset();
}

/* ------------------------------------------------------------------------ */
template <typename T, template <typename> class BaseField,
          template <typename> class Generator>
void RandomInternalField<T, BaseField, Generator>::setRandomDistribution(
    const RandomParameter<T> & param) {
  random_parameter = param;
  this->reset();
}

/* ------------------------------------------------------------------------ */
template <typename T, template <typename> class BaseField,
          template <typename> class Generator>
void RandomInternalField<T, BaseField, Generator>::printself(
    std::ostream & stream, int indent [[gnu::unused]]) const {
  stream << "RandomInternalField [ ";
  random_parameter.printself(stream);
  stream << " ]";
#if !defined(AKANTU_NDEBUG)
  if (AKANTU_DEBUG_TEST(dblDump)) {
    stream << std::endl;
    BaseField<T>::printself(stream, indent);
  }
#endif
}

/* -------------------------------------------------------------------------- */
template <typename T, template <typename> class BaseField,
          template <typename> class Generator>
void RandomInternalField<T, BaseField, Generator>::setArrayValues(T * begin,
                                                                  T * end) {
  random_parameter.template setValues<Generator>(begin, end);
}

/* -------------------------------------------------------------------------- */
template <typename T, template <typename> class BaseField,
          template <typename> class Generator>
inline RandomInternalField<T, BaseField, Generator>::operator Real() const {
  return random_parameter.getBaseValue();
}

/* -------------------------------------------------------------------------- */
template <>
inline void ParameterTyped<RandomInternalField<Real>>::setAuto(
    const ParserParameter & in_param) {
  Parameter::setAuto(in_param);
  RandomParameter<Real> r = in_param;
  param.setRandomDistribution(r);
}

/* -------------------------------------------------------------------------- */
template <>
inline void ParameterTyped<RandomInternalField<Real>>::set(std::any value) {
  auto r = std::any_cast<RandomParameter<Real>>(value);
  param.setRandomDistribution(r);
}
/* -------------------------------------------------------------------------- */

} // namespace akantu

#endif /* AKANTU_RANDOM_INTERNAL_FIELD_TMPL_HH_ */
