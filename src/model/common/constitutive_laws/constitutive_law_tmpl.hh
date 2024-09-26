/**
 * @file   constitutive_law_tmpl.hh *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Tue Jul 27 2010
 * @date last modification: Wed Feb 21 2018
 *
 * @brief  Implementation of the templated part of the constitutive law class
 *
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "constitutive_law.hh" // NOLINT
#include "constitutive_laws_handler.hh"
#include "fe_engine.hh"
#include "internal_field.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_CONSTITUTIVE_LAW_TMPL_HH
#define AKANTU_CONSTITUTIVE_LAW_TMPL_HH

namespace akantu {

/* -------------------------------------------------------------------------- */
template <typename T, template <typename Type> class InternalFieldType>
inline InternalFieldType<T> &
ConstitutiveLawInternalHandler::registerInternal(const ID & id,
                                                 Int nb_component) {
  return this->registerInternal<T, InternalFieldType>(
      id, nb_component, this->default_fe_engine_id);
}

/* -------------------------------------------------------------------------- */
template <typename T, template <typename Type> class InternalFieldType>
inline InternalFieldType<T> & ConstitutiveLawInternalHandler::registerInternal(
    const ID & id, Int nb_component, const ID & fe_engine_id) {
  return this->registerInternal<T, InternalFieldType>(
      id, nb_component, fe_engine_id, this->getElementFilter());
}
/* -------------------------------------------------------------------------- */
template <typename T, template <typename Type> class InternalFieldType>
inline InternalFieldType<T> & ConstitutiveLawInternalHandler::registerInternal(
    const ID & id, Int nb_component, const ID & fe_engine_id,
    const ElementTypeMapArray<Idx> & element_filter) {
  auto && internal =
      std::shared_ptr<InternalFieldType<T>>(new InternalFieldType<T>(
          id, *this, this->spatial_dimension, fe_engine_id, element_filter));
  internal->initialize(nb_component);
  internal_vectors[internal->getRegisterID()] = internal;
  return *internal;
}

/* -------------------------------------------------------------------------- */
inline void ConstitutiveLawInternalHandler::unregisterInternal(const ID & id) {
  internal_vectors.erase(id);
}

/* -------------------------------------------------------------------------- */
inline void ConstitutiveLawInternalHandler::savePreviousState() {
  for (auto && [_, internal] : internal_vectors) {
    if (internal->hasHistory()) {
      internal->saveCurrentValues();
    }
  }
}

/* -------------------------------------------------------------------------- */
inline void ConstitutiveLawInternalHandler::restorePreviousState() {
  for (auto && [_, internal] : internal_vectors) {
    if (internal->hasHistory()) {
      internal->restorePreviousValues();
    }
  }
}

/* -------------------------------------------------------------------------- */
inline void ConstitutiveLawInternalHandler::resizeInternals() {
  for (auto && [_, internal] : internal_vectors) {
    internal->resize();
  }
}

/* -------------------------------------------------------------------------- */
template <typename T>
const InternalField<T> &
ConstitutiveLawInternalHandler::getInternal(const ID & id) const {
  auto it = internal_vectors.find(id);
  if (it != internal_vectors.end() and
      aka::is_of_type<InternalField<T>>(*it->second)) {
    return aka::as_type<InternalField<T>>(*it->second);
  }

  AKANTU_SILENT_EXCEPTION("The constitutive law "
                          << name << "(" << getID()
                          << ") does not contain an internal " << id << " ("
                          << (getID() + ":" + id) << ")");
}

/* -------------------------------------------------------------------------- */
template <typename T>
InternalField<T> & ConstitutiveLawInternalHandler::getInternal(const ID & id) {
  if (auto it = internal_vectors.find(id);
      it != internal_vectors.end() and
      aka::is_of_type<InternalField<T>>(*it->second)) {
    return aka::as_type<InternalField<T>>(*it->second);
  }

  AKANTU_SILENT_EXCEPTION("The constitutive law "
                          << name << "(" << getID()
                          << ") does not contain an internal " << id);
}

/* -------------------------------------------------------------------------- */
template <typename T, template <typename Type> class InternalFieldType>
std::shared_ptr<InternalFieldType<T>>
ConstitutiveLawInternalHandler::getSharedPtrInternal(const ID & id) {
  if (auto it = this->internal_vectors.find(id);
      it != internal_vectors.end() and
      aka::is_of_type<InternalFieldType<T>>(*it->second)) {
    return std::dynamic_pointer_cast<InternalFieldType<T>>(it->second);
  }

  AKANTU_SILENT_EXCEPTION("The constitutive law "
                          << name << "(" << getID()
                          << ") does not contain an internal " << id);
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline bool ConstitutiveLawInternalHandler::isInternal(
    const ID & id, const ElementKind & element_kind) const {
  auto it = internal_vectors.find(id);

  return ((it != internal_vectors.end()) and
          aka::is_of_type<InternalField<T>>(*it->second) and
          aka::as_type<InternalField<T>>(*it->second).getElementKind() ==
              element_kind);
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
template <typename T>
const Array<T> &
ConstitutiveLawInternalHandler::getArray(const ID & vect_id, ElementType type,
                                         GhostType ghost_type) const {
  if (isInternal<T>(vect_id, Mesh::getKind(type))) {
    auto && internal = this->template getInternal<T>(vect_id);
    if (internal.exists(type, ghost_type)) {
      return internal(type, ghost_type);
    }

    AKANTU_SILENT_EXCEPTION(
        "The internal " << vect_id << " in the constitutive law " << name
                        << " (" << getID() << ") does not contain the type ["
                        << type << ":" << ghost_type << "]");
  }

  AKANTU_SILENT_EXCEPTION("The constitutive law "
                          << name << " (" << getID()
                          << ") does not contain an internal field "
                          << vect_id);
}

/* -------------------------------------------------------------------------- */
template <typename T>
Array<T> & ConstitutiveLawInternalHandler::getArray(const ID & vect_id,
                                                    ElementType type,
                                                    GhostType ghost_type) {
  if (isInternal<T>(vect_id, Mesh::getKind(type))) {
    auto && internal = this->template getInternal<T>(vect_id);
    if (internal.exists(type, ghost_type)) {
      return internal(type, ghost_type);
    }
    AKANTU_SILENT_EXCEPTION(
        "The internal " << vect_id << " in the constitutive law " << name
                        << " (" << getID() << ") does not contain the type ["
                        << type << ":" << ghost_type << "]");
  }

  AKANTU_SILENT_EXCEPTION("The constitutive law "
                          << name << " (" << getID()
                          << ") does not contain an internal field "
                          << vect_id);
}

/* -------------------------------------------------------------------------- */
inline void ConstitutiveLawInternalHandler::removeIntegrationPoints(
    ElementTypeMapArray<Idx> & new_numbering) {
  for (auto && [_, internal] : internal_vectors) {
    internal->removeIntegrationPoints(new_numbering);
  }
}
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawsHandler_>
ConstitutiveLaw<ConstitutiveLawsHandler_>::ConstitutiveLaw(
    ConstitutiveLawsHandler_ & handler, const ID & id, Int spatial_dimension,
    ElementKind element_kind, const ID & fe_engine_id)
    : ConstitutiveLawInternalHandler(id, spatial_dimension, fe_engine_id),
      Parsable(handler.getConstitutiveLawParserType(), id), handler(handler) {

  /// for each connectivity types allocate the element filer array of
  /// the constitutive law
  this->getElementFilter().initialize(
      handler.getMesh(), _spatial_dimension = handler.getSpatialDimension(),
      _element_kind = element_kind);

  this->initialize();
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawsHandler_>
void ConstitutiveLaw<ConstitutiveLawsHandler_>::initialize() {
  registerParam("name", name, std::string(), _pat_parsable | _pat_readable);
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawsHandler_>
void ConstitutiveLaw<ConstitutiveLawsHandler_>::initConstitutiveLaw() {
  this->resizeInternals();
  this->updateInternalParameters();
  is_init = true;
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawsHandler_>
void ConstitutiveLaw<ConstitutiveLawsHandler_>::addElements(
    const Array<Element> & elements_to_add) {
  AKANTU_DEBUG_IN();

  Int law_id = handler.getConstitutiveLawIndex(name);
  for (const auto & element : elements_to_add) {
    auto index = this->addElement(element);
    handler.constitutive_law_index(element) = law_id;
    handler.constitutive_law_local_numbering(element) = index;
  }

  this->resizeInternals();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawsHandler_>
void ConstitutiveLaw<ConstitutiveLawsHandler_>::removeElements(
    const Array<Element> & elements_to_remove) {
  AKANTU_DEBUG_IN();

  auto el_begin = elements_to_remove.begin();
  auto el_end = elements_to_remove.end();

  if (elements_to_remove.empty()) {
    return;
  }

  auto & mesh = handler.getMesh();

  ElementTypeMapArray<Idx> constitutive_law_local_new_numbering(
      "remove constitutive law filter elem", id);

  constitutive_law_local_new_numbering.initialize(
      mesh, _element_filter = &getElementFilter(),
      _element_kind = _ek_not_defined, _with_nb_element = true);

  ElementTypeMapArray<Idx> element_filter_tmp("element_filter_tmp", id);

  element_filter_tmp.initialize(mesh, _element_filter = &getElementFilter(),
                                _element_kind = _ek_not_defined);

  ElementTypeMap<Idx> new_ids, element_ids;

  for_each_element(
      mesh,
      [&](auto && el) {
        if (not new_ids(el.type, el.ghost_type)) {
          element_ids(el.type, el.ghost_type) = 0;
        }

        auto & element_id = element_ids(el.type, el.ghost_type);
        auto l_el = Element{el.type, element_id, el.ghost_type};
        if (std::find(el_begin, el_end, el) != el_end) {
          constitutive_law_local_new_numbering(l_el) = Idx(-1);
          return;
        }

        element_filter_tmp(el.type, el.ghost_type).push_back(el.element);
        if (not new_ids(el.type, el.ghost_type)) {
          new_ids(el.type, el.ghost_type) = 0;
        }

        auto & new_id = new_ids(el.type, el.ghost_type);

        constitutive_law_local_new_numbering(l_el) = new_id;
        handler.constitutive_law_local_numbering(el) = new_id;

        ++new_id;
        ++element_id;
      },
      _element_filter = &getElementFilter(), _element_kind = _ek_not_defined);

  for (auto ghost_type : ghost_types) {
    for (const auto & type : getElementFilter().elementTypes(
             _ghost_type = ghost_type, _element_kind = _ek_not_defined)) {
      getElementFilter(type, ghost_type)
          .copy(element_filter_tmp(type, ghost_type));
    }
  }

  this->removeIntegrationPoints(constitutive_law_local_new_numbering);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawsHandler_>
void ConstitutiveLaw<ConstitutiveLawsHandler_>::onElementsAdded(
    const Array<Element> & /*unused*/, const NewElementsEvent & /*unused*/) {
  this->resizeInternals();
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawsHandler_>
void ConstitutiveLaw<ConstitutiveLawsHandler_>::onElementsRemoved(
    const Array<Element> & element_list,
    const ElementTypeMapArray<Idx> & new_numbering,
    const RemovedElementsEvent & /*event*/) {

  auto my_num = handler.getInternalIndexFromID(getID());

  ElementTypeMapArray<Idx> constitutive_law_local_new_numbering(
      "remove constitutive law filter elem", getID());

  auto el_begin = element_list.begin();
  auto el_end = element_list.end();

  for (auto && gt : ghost_types) {
    for (auto && type :
         new_numbering.elementTypes(_all_dimensions, gt, _ek_not_defined)) {

      if (not getElementFilter().exists(type, gt) or
          getElementFilter(type, gt).empty()) {
        continue;
      }

      auto & elem_filter = getElementFilter(type, gt);
      auto & law_indexes = this->handler.constitutive_law_index(type, gt);
      auto & law_loc_num =
          this->handler.constitutive_law_local_numbering(type, gt);
      auto nb_element = this->handler.getMesh().getNbElement(type, gt);

      // all constitutive laws will resized to the same size...
      law_indexes.resize(nb_element);
      law_loc_num.resize(nb_element);

      if (not constitutive_law_local_new_numbering.exists(type, gt)) {
        constitutive_law_local_new_numbering.alloc(elem_filter.size(), 1, type,
                                                   gt);
      }

      auto & law_renumbering = constitutive_law_local_new_numbering(type, gt);
      const auto & renumbering = new_numbering(type, gt);
      Array<Idx> elem_filter_tmp;
      Int ni = 0;
      Element el{type, 0, gt};

      for (auto && [i, el_id] : enumerate(elem_filter)) {
        el.element = el_id;

        if (std::find(el_begin, el_end, el) == el_end) {
          auto new_el = renumbering(el_id);
          AKANTU_DEBUG_ASSERT(new_el != -1,
                              "A not removed element as been badly renumbered");
          elem_filter_tmp.push_back(new_el);
          law_renumbering(i) = ni;
          law_indexes(new_el) = my_num;
          law_loc_num(new_el) = ni;
          ++ni;
        } else {
          law_renumbering(i) = -1;
        }
      }

      elem_filter.resize(elem_filter_tmp.size());
      elem_filter.copy(elem_filter_tmp);
    }
  }

  this->removeIntegrationPoints(constitutive_law_local_new_numbering);
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawsHandler_>
template <typename T>
inline void ConstitutiveLaw<ConstitutiveLawsHandler_>::packInternalFieldHelper(
    const InternalField<T> & data_to_pack, CommunicationBuffer & buffer,
    const Array<Element> & elements) const {
  DataAccessor::packElementalDataHelper<T>(data_to_pack, buffer, elements,
                                           data_to_pack.getFEEngine());
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawsHandler_>
template <typename T>
inline void
ConstitutiveLaw<ConstitutiveLawsHandler_>::unpackInternalFieldHelper(
    InternalField<T> & data_to_unpack, CommunicationBuffer & buffer,
    const Array<Element> & elements) {
  DataAccessor::unpackElementalDataHelper<T>(data_to_unpack, buffer, elements,
                                             data_to_unpack.getFEEngine());
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawsHandler_>
inline Element ConstitutiveLaw<ConstitutiveLawsHandler_>::convertToLocalElement(
    const Element & global_element) const {
#ifndef AKANTU_NDEBUG
  auto model_law_index = handler.getConstitutiveLawByElement()(global_element);
  auto law_index = handler.getConstitutiveLawIndex(this->name);
  AKANTU_DEBUG_ASSERT(model_law_index == law_index,
                      "Conversion of a global  element in a local element for "
                      "the wrong constitutive law "
                          << this->name << std::endl);
#endif
  auto local_element = global_element;
  local_element.element =
      handler.getConstitutiveLawLocalNumbering()(global_element);

  return local_element;
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawsHandler_>
inline Element
ConstitutiveLaw<ConstitutiveLawsHandler_>::convertToGlobalElement(
    const Element & local_element) const {
  auto global_element = local_element;
  global_element.element = this->getElementFilter()(local_element);
  return global_element;
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawsHandler_>
inline Idx
ConstitutiveLaw<ConstitutiveLawsHandler_>::addElement(const Element & element) {
  auto & el_filter = this->getElementFilter(element.type, element.ghost_type);
  el_filter.push_back(element.element);
  return el_filter.size() - 1;
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawsHandler_>
inline const Parameter &
ConstitutiveLaw<ConstitutiveLawsHandler_>::getParam(const ID & param) const {
  try {
    return get(param);
  } catch (...) {
    AKANTU_EXCEPTION("No parameter " << param << " in the constitutive law"
                                     << getID());
  }
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawsHandler_>
template <typename T>
inline void
ConstitutiveLaw<ConstitutiveLawsHandler_>::setParam(const ID & param, T value) {
  try {
    set(param, value);
  } catch (...) {
    AKANTU_EXCEPTION("No parameter " << param << " in the constitutive law "
                                     << getID());
  }
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawsHandler_>
template <typename T>
inline ElementTypeMap<Int>
ConstitutiveLaw<ConstitutiveLawsHandler_>::getInternalDataPerElem(
    const ID & field_id, ElementKind element_kind) const {

  if (not this->template isInternal<T>(field_id, element_kind)) {
    AKANTU_EXCEPTION("Cannot find internal field "
                     << id << " in the constitutive law " << this->name);
  }

  const auto & internal_field = this->template getInternal<T>(field_id);
  const FEEngine & fe_engine = internal_field.getFEEngine();
  auto nb_data_per_quad = internal_field.getNbComponent();

  ElementTypeMap<Int> res;
  for (auto ghost_type : ghost_types) {
    for (auto && type : internal_field.elementTypes(ghost_type)) {
      auto nb_quadrature_points =
          fe_engine.getNbIntegrationPoints(type, ghost_type);
      res(type, ghost_type) = nb_data_per_quad * nb_quadrature_points;
    }
  }

  return res;
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawsHandler_>
template <typename T>
void ConstitutiveLaw<ConstitutiveLawsHandler_>::flattenInternal(
    const std::string & field_id, ElementTypeMapArray<T> & internal_flat,
    const GhostType ghost_type, ElementKind element_kind) const {

  if (!this->template isInternal<T>(field_id, element_kind)) {
    AKANTU_EXCEPTION("Cannot find internal field "
                     << id << " in the constitutive law " << this->name);
  }

  const auto & internal_field = this->template getInternal<T>(field_id);

  const auto & fe_engine = internal_field.getFEEngine();
  const auto & mesh = fe_engine.getMesh();

  for (auto && type : internal_field.filterTypes(ghost_type)) {
    const auto & src_vect = internal_field(type, ghost_type);
    const auto & filter = internal_field.getFilter(type, ghost_type);

    // total number of elements in the corresponding mesh
    auto nb_element_dst = mesh.getNbElement(type, ghost_type);
    // number of element in the internal field
    auto nb_element_src = filter.size();
    // number of quadrature points per elem
    auto nb_quad_per_elem = fe_engine.getNbIntegrationPoints(type);
    // number of data per quadrature point
    auto nb_data_per_quad = internal_field.getNbComponent();

    if (not internal_flat.exists(type, ghost_type)) {
      internal_flat.alloc(nb_element_dst * nb_quad_per_elem, nb_data_per_quad,
                          type, ghost_type);
    }

    if (nb_element_src == 0) {
      continue;
    }

    // number of data per element
    auto nb_data = nb_quad_per_elem * nb_data_per_quad;

    auto & dst_vect = internal_flat(type, ghost_type);
    dst_vect.resize(nb_element_dst * nb_quad_per_elem);

    auto it_dst = make_view(dst_vect, nb_data).begin();

    for (auto && [element, src_value] :
         zip(filter, make_view(src_vect, nb_data))) {
      it_dst[element] = src_value;
    }
  }
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawsHandler_>
template <typename T>
void ConstitutiveLaw<ConstitutiveLawsHandler_>::inflateInternal(
    const std::string & field_id, const ElementTypeMapArray<T> & field,
    GhostType ghost_type, ElementKind element_kind) {
  if (not this->template isInternal<T>(field_id, element_kind)) {
    AKANTU_EXCEPTION("Cannot find internal field " << id << " in material "
                                                   << this->name);
  }

  auto & internal_field = this->template getInternal<T>(field_id);
  const auto & fe_engine = internal_field.getFEEngine();

  for (auto && type : field.elementTypes(_ghost_type = ghost_type)) {
    if (not internal_field.exists(type, ghost_type)) {
      continue;
    }
    const auto & filter = internal_field.getFilter(type, ghost_type);

    const auto & src_array = field(type, ghost_type);
    auto & dest_array = internal_field(type, ghost_type);

    auto nb_quad_per_elem = fe_engine.getNbIntegrationPoints(type);
    auto nb_component = src_array.getNbComponent();

    AKANTU_DEBUG_ASSERT(
        field.size() == fe_engine.getMesh().getNbElement(type, ghost_type) *
                            nb_quad_per_elem,
        "The ElementTypeMapArray to inflate is not of the proper size");
    AKANTU_DEBUG_ASSERT(dest_array.getNbComponent() == nb_component,
                        "The ElementTypeMapArray has not the proper "
                        "number of components");

    auto src =
        make_view(field(type, ghost_type), nb_component, nb_quad_per_elem)
            .begin();
    for (auto && [el, dest] :
         zip(filter, make_view(dest_array, nb_component, nb_quad_per_elem))) {
      dest = src[el];
    }
  }
}

} // namespace akantu

#endif // AKANTU_CONSTITUTIVE_LAW_TMPL_HH
