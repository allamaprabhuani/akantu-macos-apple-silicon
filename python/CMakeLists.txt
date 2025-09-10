#===============================================================================
# Copyright (©) 2014-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# This file is part of Akantu
# 
# Akantu is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
# 
# Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
# 
# You should have received a copy of the GNU Lesser General Public License along
# with Akantu. If not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================


if(NOT SKBUILD)
  package_get_all_include_directories(
    AKANTU_LIBRARY_INCLUDE_DIRS
    )

  package_get_all_external_informations(
    PRIVATE_INCLUDE AKANTU_PRIVATE_EXTERNAL_INCLUDE_DIR
    INTERFACE_INCLUDE AKANTU_INTERFACE_EXTERNAL_INCLUDE_DIR
    LIBRARIES AKANTU_EXTERNAL_LIBRARIES
    )
endif()

set(PYAKANTU_SRCS
  py_aka_common.cc
  py_aka_error.cc
  py_akantu.cc
  py_boundary_conditions.cc
  py_constitutive_law.cc
  py_constitutive_law_selector.cc
  py_dof_manager.cc
  py_dumpable.cc
  py_fe_engine.cc
  py_group_manager.cc
  py_integration_scheme.cc
  py_mesh.cc
  py_model.cc
  py_parser.cc
  py_solver.cc
  )

package_is_activated(solid_mechanics _is_activated)
if (_is_activated)
  list(APPEND PYAKANTU_SRCS
    py_solid_mechanics_model.cc
    py_material.cc
    )
endif()

package_is_activated(cohesive_element _is_activated)
if (_is_activated)
  list(APPEND PYAKANTU_SRCS
    py_solid_mechanics_model_cohesive.cc
    py_fragment_manager.cc
    )
endif()

package_is_activated(diffusion _is_activated)
if (_is_activated)
  list(APPEND PYAKANTU_SRCS
    py_heat_transfer_model.cc
    )
endif()


package_is_activated(contact_mechanics _is_activated)
if(_is_activated)
  list(APPEND PYAKANTU_SRCS
    py_contact_mechanics_model.cc
    py_model_couplers.cc
    )
endif()

package_is_activated(phase_field _is_activated)
if (_is_activated)
  list(APPEND PYAKANTU_SRCS
    py_phase_field_model.cc
    py_phasefield.cc
    )
endif()

package_is_activated(structural_mechanics _is_activated)
if (_is_activated)
  list(APPEND PYAKANTU_SRCS
    py_structural_mechanics_model.cc
    )
endif()

pybind11_add_module(py11_akantu ${PYAKANTU_SRCS})

# to avoid compilation warnings from pybind11
target_include_directories(py11_akantu
  SYSTEM BEFORE
  PRIVATE ${PYBIND11_INCLUDE_DIR}
  PRIVATE ${pybind11_INCLUDE_DIR}
  PRIVATE ${Python_INCLUDE_DIRS})

target_link_libraries(py11_akantu PUBLIC akantu)
set_target_properties(py11_akantu PROPERTIES
  DEBUG_POSTFIX ""
  LIBRARY_OUTPUT_DIRECTORY akantu)

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  target_compile_options(py11_akantu PUBLIC -fsized-deallocation)
endif()


file(COPY akantu DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

if(NOT Python_MAJOR)
  set(Python_VERSION_MAJOR ${PYTHON_VERSION_MAJOR})
  set(Python_VERSION_MINOR ${PYTHON_VERSION_MINOR})
endif()

if(NOT SKBUILD)
  set(_python_install_dir
    ${CMAKE_INSTALL_LIBDIR}/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages/akantu)
else()
  set(_python_install_dir python/akantu)
endif()

install(TARGETS py11_akantu
  LIBRARY DESTINATION ${_python_install_dir})

if(NOT SKBUILD)
  install(DIRECTORY akantu
    DESTINATION ${_python_install_dir}
    FILES_MATCHING PATTERN "*.py")
endif()
