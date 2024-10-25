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

#ifndef AKANTU_COHESIVE_INTERNAL_FIELD_TMPL_HH_
#define AKANTU_COHESIVE_INTERNAL_FIELD_TMPL_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
template <typename T>
CohesiveInternalField<T>::CohesiveInternalField(
    const ID & id, ConstitutiveLawInternalHandler & constitutive_law, Int dim,
    const ID & fem_id, const ElementTypeMapArray<Idx> & element_filter)
    : InternalField<T>(id, constitutive_law, dim, fem_id, element_filter) {
  this->element_kind = _ek_cohesive;
}

template <typename T>
void CohesiveInternalField<T>::initialize(Int nb_component) {
  this->internalInitialize(nb_component);
}

/* -------------------------------------------------------------------------- */
template <>
inline CohesiveInternalField<Real> &
ConstitutiveLawInternalHandler::registerInternal<Real, CohesiveInternalField>(
    const ID & id, Int nb_component) {
  return this->registerInternal<Real, CohesiveInternalField>(
      id, nb_component, "CohesiveFEEngine");
}
/* -------------------------------------------------------------------------- */
template <>
inline CohesiveRandomInternalField<Real> &
ConstitutiveLawInternalHandler::registerInternal<
    Real, CohesiveRandomInternalField>(const ID & id, Int nb_component) {
  return this->registerInternal<Real, CohesiveRandomInternalField>(
      id, nb_component, "CohesiveFEEngine");
}

/* -------------------------------------------------------------------------- */
template <typename T>
FacetInternalField<T>::FacetInternalField(
    const ID & id, ConstitutiveLawInternalHandler & constitutive_law, Int dim,
    const ID & fem_id, const ElementTypeMapArray<Idx> & element_filter)
    : InternalField<T>(id, constitutive_law, dim, fem_id, element_filter) {
  this->spatial_dimension -= 1;
  this->element_kind = _ek_regular;
}

template <typename T> void FacetInternalField<T>::initialize(Int nb_component) {
  this->internalInitialize(nb_component);
}

/* -------------------------------------------------------------------------- */
template <>
inline FacetInternalField<Real> &
ConstitutiveLawInternalHandler::registerInternal<Real, FacetInternalField>(
    const ID & id, Int nb_component) {
  return this->registerInternal<Real, FacetInternalField>(
      id, nb_component, "FacetsFEEngine",
      aka::as_type<MaterialCohesive>(*this).getFacetFilter());
}
/* -------------------------------------------------------------------------- */
template <>
inline FacetRandomInternalField<Real> &
ConstitutiveLawInternalHandler::registerInternal<
    Real, FacetRandomInternalField>(const ID & id, Int nb_component) {
  return this->registerInternal<Real, FacetRandomInternalField>(
      id, nb_component, "FacetsFEEngine",
      aka::as_type<MaterialCohesive>(*this).getFacetFilter());
}

/* -------------------------------------------------------------------------- */
template <>
inline void ParameterTyped<FacetRandomInternalField<Real>>::setAuto(
    const ParserParameter & in_param) {
  Parameter::setAuto(in_param);
  RandomParameter<Real> r = in_param;
  param.setRandomDistribution(r);
}

/* -------------------------------------------------------------------------- */

template <>
inline void
ParameterTyped<FacetRandomInternalField<Real>>::set(std::any value) {
  auto r = std::any_cast<RandomParameter<Real>>(value);
  param.setRandomDistribution(r);
}

/* -------------------------------------------------------------------------- */
template <>
inline void ParameterTyped<CohesiveRandomInternalField<Real>>::setAuto(
    const ParserParameter & in_param) {
  Parameter::setAuto(in_param);
  RandomParameter<Real> r = in_param;
  param.setRandomDistribution(r);
}

/* -------------------------------------------------------------------------- */
template <>
inline void
ParameterTyped<CohesiveRandomInternalField<Real>>::set(std::any value) {
  auto r = std::any_cast<RandomParameter<Real>>(value);
  param.setRandomDistribution(r);
}

} // namespace akantu

#endif /* AKANTU_COHESIVE_INTERNAL_FIELD_TMPL_HH_ */
