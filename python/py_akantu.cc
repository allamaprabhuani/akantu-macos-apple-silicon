/**
 * Copyright (©) 2018-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_config.hh"
// for NLSNotConvergedException
#include "non_linear_solver.hh"
/* -------------------------------------------------------------------------- */
#include "py_aka_common.hh"
#include "py_aka_error.hh"
#include "py_boundary_conditions.hh"
#include "py_constitutive_law.hh"
#include "py_constitutive_law_selector.hh"
#include "py_dof_manager.hh"
#include "py_dumpable.hh"
#include "py_fe_engine.hh"
#include "py_group_manager.hh"
#include "py_integration_scheme.hh"
#include "py_mesh.hh"
#include "py_model.hh"
#include "py_parser.hh"
#include "py_solver.hh"

#if defined(AKANTU_SOLID_MECHANICS)
#include "py_material.hh"
#include "py_solid_mechanics_model.hh"
#endif

#if defined(AKANTU_DIFFUSION)
#include "py_heat_transfer_model.hh"
#endif

#if defined(AKANTU_COHESIVE_ELEMENT)
#include "py_fragment_manager.hh"
#include "py_solid_mechanics_model_cohesive.hh"
#endif

#if defined(AKANTU_CONTACT_MECHANICS)
#include "py_contact_mechanics_model.hh"
#include "py_model_couplers.hh"
#endif

#if defined(AKANTU_PHASE_FIELD)
#include "py_phase_field_model.hh"
#endif

#if defined(AKANTU_STRUCTURAL_MECHANICS)
#include "py_structural_mechanics_model.hh"
#endif

/* -------------------------------------------------------------------------- */
#include <aka_error.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/pybind11.h>
/* -------------------------------------------------------------------------- */
#include <iostream>
/* -------------------------------------------------------------------------- */

namespace py = pybind11;

namespace akantu {
void register_all(pybind11::module & mod) {
  register_initialize(mod);
  register_enums(mod);
  register_error(mod);
  register_functions(mod);
  register_parser(mod);
  register_solvers(mod);

  register_dumpable(mod);
  register_group_manager(mod);
  register_mesh(mod);

  register_fe_engine(mod);

  register_integration_schemes(mod);
  register_dof_manager(mod);

  register_boundary_conditions(mod);
  register_model(mod);
  register_constitutive_law_selector(mod);
  register_constitutive_law_internal_handler(mod);

#if defined(AKANTU_DIFFUSION)
  register_heat_transfer_model(mod);
#endif

#if defined(AKANTU_SOLID_MECHANICS)
  register_solid_mechanics_model(mod);
  register_material(mod);
#endif

#if defined(AKANTU_COHESIVE_ELEMENT)
  register_solid_mechanics_model_cohesive(mod);
  register_fragment_manager(mod);
#endif

#if defined(AKANTU_STRUCTURAL_MECHANICS)
  register_structural_mechanics_model(mod);
#endif

#if defined(AKANTU_CONTACT_MECHANICS)
  register_contact_mechanics_model(mod);
  register_model_couplers(mod);
#endif

#if defined(AKANTU_PHASE_FIELD)
  register_phase_field_model(mod);
  register_phase_field_coupler(mod);
#endif
}
} // namespace akantu

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
PYBIND11_MODULE(py11_akantu, mod) {
  mod.doc() = "Akantu python interface";

  static py::exception<akantu::debug::Exception> akantu_exception(mod,
                                                                  "Exception");

  static py::exception<akantu::debug::NLSNotConvergedException>
      akantu_exception_nls_not_converged(mod, "NLSNotConvergedException");

  py::register_exception_translator([](std::exception_ptr ptr) {
    try {
      if (ptr) {
        std::rethrow_exception(ptr);
      }
    } catch (akantu::debug::NLSNotConvergedException & e) {
#if PYBIND11_VERSION_HEX >= 0x020C0000
      py::set_error(akantu_exception_nls_not_converged, e.info().c_str());
#else
      akantu_exception_nls_not_converged(e.info().c_str());
#endif
    } catch (akantu::debug::Exception & e) {
      if (akantu::debug::debugger.printBacktrace()) {
        akantu::debug::printBacktrace();
      }
#if PYBIND11_VERSION_HEX >= 0x020C0000
      py::set_error(akantu_exception, e.info().c_str());
#else
      akantu_exception(e.info().c_str());
#endif
    }
  });

  akantu::register_all(mod);

  mod.def("has_mpi",
          []() {
#if defined(AKANTU_USE_MPI)
            return true;
#else
            return false;
#endif
          })
      .def("getVersion", &akantu::getVersion);

} // Module akantu
