/**
 * Copyright (©) 2016-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include <utility>

#include "parameter_registry.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

Parameter::Parameter() = default;

/* -------------------------------------------------------------------------- */
Parameter::Parameter(std::string name, std::string description,
                     ParameterAccessType param_type)
    : name(std::move(name)), description(std::move(description)),
      param_type(param_type) {}

/* -------------------------------------------------------------------------- */
bool Parameter::isWritable() const { return (param_type & _pat_writable) != 0; }

/* -------------------------------------------------------------------------- */
bool Parameter::isReadable() const { return (param_type & _pat_readable) != 0; }

/* -------------------------------------------------------------------------- */
bool Parameter::isInternal() const { return (param_type & _pat_internal) != 0; }

/* -------------------------------------------------------------------------- */
bool Parameter::isParsable() const { return (param_type & _pat_parsable) != 0; }

/* -------------------------------------------------------------------------- */
void Parameter::setAccessType(ParameterAccessType ptype) {
  this->param_type = ptype;
}

/* -------------------------------------------------------------------------- */
void Parameter::printself(std::ostream & stream) const {
  stream << " ";
  if (isInternal()) {
    stream << "iii";
  } else {
    if (isReadable()) {
      stream << "r";
    } else {
      stream << "-";
    }

    if (isWritable()) {
      stream << "w";
    } else {
      stream << "-";
    }

    if (isParsable()) {
      stream << "p";
    } else {
      stream << "-";
    }
  }
  stream << " ";

  std::stringstream sstr;
  sstr << name;
  UInt width = std::max(int(10 - sstr.str().length()), 0);
  sstr.width(width);

  if (not description.empty()) {
    sstr << " [" << description << "]";
  }

  stream << sstr.str();
  width = std::max(int(50 - sstr.str().length()), 0);
  stream.width(width);

  stream << " : ";
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
ParameterRegistry::ParameterRegistry() = default;

/* -------------------------------------------------------------------------- */
ParameterRegistry::~ParameterRegistry() = default;

/* -------------------------------------------------------------------------- */
void ParameterRegistry::printself(std::ostream & stream, int indent) const {
  std::string space(indent, AKANTU_INDENT);

  for (auto && [_, param] : params) {
    stream << space;
    param->printself(stream);
  }

  for (auto && [_, registry] : sub_registries) {
    stream << space << "Registry [" << std::endl;
    registry.get().printself(stream, indent + 1);
    stream << space << "]";
  }
}

/* -------------------------------------------------------------------------- */
void ParameterRegistry::registerSubRegistry(const ID & id,
                                            ParameterRegistry & registry) {
  sub_registries.insert_or_assign(id, std::ref(registry));
}

/* -------------------------------------------------------------------------- */
void ParameterRegistry::setParameterAccessType(const std::string & name,
                                               ParameterAccessType ptype) {
  auto it = params.find(name);
  if (it == params.end()) {
    AKANTU_CUSTOM_EXCEPTION(debug::ParameterUnexistingException(name, *this));
  }
  Parameter & param = *(it->second);
  param.setAccessType(ptype);
}

/* -------------------------------------------------------------------------- */
void ParameterRegistry::set(const std::string & name, std::any value) {
  this->setMixed(name, value);
  updateInternalParameters();
}
/* -------------------------------------------------------------------------- */

void ParameterRegistry::setMixed(const std::string & name, std::any value) {
  auto it = params.find(name);
  if (it == params.end()) {
    if (consisder_sub) {
      for (auto && [_, registry] : sub_registries) {
        registry.get().setMixed(name, value);
      }
    } else {
      AKANTU_CUSTOM_EXCEPTION(debug::ParameterUnexistingException(name, *this));
    }
  } else {
    Parameter & param = *(it->second);
    param.set(value);
  }
}
/* -------------------------------------------------------------------------- */

} // namespace akantu
