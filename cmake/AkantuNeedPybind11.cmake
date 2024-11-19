#===============================================================================
# Copyright (©) 2021-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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


if(DEFINED AKANTU_NEED_PYBIND11_LOADED)
  return()
endif()
set(AKANTU_NEED_PYBIND11_LOADED TRUE)

set(AKANTU_USE_SYSTEM_PYBIND11 AUTO CACHE STRING
  "Should akantu compile the third-party: pybind11")
mark_as_advanced(AKANTU_USE_SYSTEM_PYBIND11)
set_property(CACHE AKANTU_USE_SYSTEM_PYBIND11 PROPERTY STRINGS ON OFF AUTO)

set(PYBIND11_PYTHON_VERSION ${AKANTU_PREFERRED_PYTHON_VERSION} CACHE INTERNAL "")

set(AKANTU_PYBIND11_VERSION 2.13.6)

if (AKANTU_USE_SYSTEM_PYBIND11 MATCHES "(ON|AUTO)")
  if(CMAKE_VERSION VERSION_GREATER 3.27)
    cmake_policy(SET CMP0148 OLD)
  endif()

  find_package(pybind11 ${AKANTU_PYBIND11_VERSION} QUIET)
else()
  set(pybind11_FOUND FALSE)
endif()

mark_as_advanced(
  pybind11_DIR
  PYBIND11_PYTHONLIBS_OVERWRITE
  )

if (NOT pybind11_FOUND)
  set(PYBIND11_VERSION "v${AKANTU_PYBIND11_VERSION}")
  set(PYBIND11_GIT "https://github.com/pybind/pybind11.git")
  set(CMAKE_CXX_STANDARD ${AKANTU_CXX_STANDARD}) # Otherwhy pybind11 default to cxx14
  include(${PROJECT_SOURCE_DIR}/third-party/cmake/pybind11.cmake)
else()
  message(STATUS "Found pybind11: ${pybind11_INCLUDE_DIRS} (${pybind11_VERSION})")
endif()
