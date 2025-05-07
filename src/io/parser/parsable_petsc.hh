/**
 * Copyright (©) 2012-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "parsable.hh"
/* -------------------------------------------------------------------------- */
#include <petscsys.h>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_PARSABLE_PETSC_HH
#define AKANTU_PARSABLE_PETSC_HH

namespace akantu {

template <class ParsableParent, class PETScType>
class ParsablePETSc : public ParsableParent {
public:
  template <class... Args>
  ParsablePETSc(
      PETScType & petsc_object,
      std::function<PetscErrorCode(PETScType)> set_object_from_options,
      Args &&... args)
      : ParsableParent(std::forward<Args>(args)...), petsc_object(petsc_object),
        set_object_from_options(std::move(set_object_from_options)) {
    PetscOptionsCreate(&options);
  }
  ~ParsablePETSc() override = default;

  /* ------------------------------------------------------------------------ */
  void updateInternalParameters() override {
    if (not petsc_object) {
      return;
    }

    // set the options from parser
    PetscOptionsPush(options);
    set_object_from_options(petsc_object);
    PetscOptionsPop();

    // set the options from command line and environment
    set_object_from_options(petsc_object);
  }

  /* ------------------------------------------------------------------------ */
  void parseSection(const ParserSection & section) override {
    auto && [params_begin, params_end] = section.getParameters();
    for (auto && param : range(params_begin, params_end)) {
      exceptIfWrongParameter(param);
    }

    Parsable::parseSection(section);

    auto parameters = section.getParameters();
    for (auto && param : range(parameters.first, parameters.second)) {
      PetscOptionsSetValue(options, ("-" + param.getName()).c_str(),
                           param.getValue().c_str());
    }
    updateInternalParameters();
  }

  /* ------------------------------------------------------------------------ */
  void set(const std::string & name, std::any value) override {
    exceptIfWrongParameter(name);

    if (this->hasParameter(name)) {
      Parsable::set(name, value);
    } else {
      std::string option{};
      if (value.type() == typeid(Real)) {
        std::stringstream sstr;
        sstr << std::any_cast<Real>(value);
        option = sstr.str();
      } else if (value.type() == typeid(Int)) {
        option = std::to_string(std::any_cast<Int>(value));
      } else if (value.type() == typeid(std::string)) {
        option = std::any_cast<std::string>(value);
      } else if (value.type() == typeid(const char *)) {
        option = std::any_cast<const char *>(value);
      }

      if (not option.empty()) {
        PetscOptionsSetValue(options, ("-" + name).c_str(), option.c_str());
      }
    }
    updateInternalParameters();
  }

private:
  void exceptIfWrongParameter(const std::string param) const {
    if (auto it = akantu_to_petsc_options.find(param);
        it != akantu_to_petsc_options.end()) {

      const char * type{nullptr};
      const char * name{nullptr};
      auto object = reinterpret_cast<PetscObject>(petsc_object);
      PetscObjectGetType(object, &type);
      PetscObjectGetName(object, &name);
      AKANTU_EXCEPTION(
          "Cannot set the parameter \""
          << it->first << "\" on object \"" << name << "\" of type \"" << type
          << "\" use the PETSc option \"" << it->second << "\" instead.");
    }
  }

  PetscOptions options{};
  PETScType & petsc_object;
  std::function<PetscErrorCode(PETScType)> set_object_from_options;

protected:
  std::map<ID, ID> akantu_to_petsc_options;
};

} // namespace akantu

#endif // AKANTU_PARSABLE_PETSC_HH
