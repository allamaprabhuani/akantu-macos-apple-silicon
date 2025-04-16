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
    for (auto && [akantu, petsc] : akantu_to_petsc_options) {
      auto & value = this->get(akantu);
      PetscOptionsSetValue(options, ("-" + petsc).c_str(),
                           value.to_string().c_str());
    }

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
    Parsable::parseSection(section);

    auto parameters = section.getParameters();
    for (auto && param : range(parameters.first, parameters.second)) {
      PetscOptionsSetValue(options, ("-" + param.getName()).c_str(),
                           param.getValue().c_str());
    }
  }

  /* ------------------------------------------------------------------------ */
  void set(const std::string & name, std::any value) override {
    if (this->hasParameter(name)) {
      Parsable::set(name, value);
    } else {
      try {
        std::string option = std::any_cast<const char *>(value);
        PetscOptionsSetValue(options, ("-" + name).c_str(), option.c_str());
      } catch (std::bad_any_cast & c) {
        AKANTU_DEBUG_WARNING("\"" << name
                                  << "\" is not a known option of the object "
                                  << this->pid);
      }
    }
    updateInternalParameters();
  }

private:
  PetscOptions options{};
  PETScType & petsc_object;
  std::function<PetscErrorCode(PETScType)> set_object_from_options;

protected:
  std::map<ID, ID> akantu_to_petsc_options;
};

} // namespace akantu

#endif // AKANTU_PARSABLE_PETSC_HH
