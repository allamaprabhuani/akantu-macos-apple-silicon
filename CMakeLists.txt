#===============================================================================
# Copyright (©) 2010-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
# @section DESCRIPTION
#-------------------------------------------------------------------------------
#                      _               _
#                     | |             | |
#                 __ _| | ____ _ _ __ | |_ _   _
#                / _` | |/ / _` | '_ \| __| | | |
#               | (_| |   < (_| | | | | |_| |_| |
#                \__,_|_|\_\__,_|_| |_|\__|\__,_|
#
#===============================================================================


#===============================================================================
# CMake Project
#===============================================================================
cmake_minimum_required(VERSION 3.16.3)

# add this options before PROJECT keyword
set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

cmake_policy(SET CMP0074 NEW)
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.30)
  # FindBOOST deprecated to old due to manylinux build
  cmake_policy(SET CMP0167 OLD)
endif()

set(AKANTU_COPYRIGHT "2010-2021, EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)")
set(AKANTU_MAINTAINER "Nicolas Richart")
set(AKANTU_HOMEPAGE_URL "https://akantu.ch")

project(Akantu
  HOMEPAGE_URL "https://akantu.ch")

enable_language(CXX)

set(ENV{RUNNING_IN_CMAKE} True)  # Used by sphinx but could be used by any script

#===============================================================================
# Misc. config for cmake
#===============================================================================
set(AKANTU_CMAKE_DIR "${PROJECT_SOURCE_DIR}/cmake")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/Modules")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE BOOL
  "Enable/Disable output of compile commands during generation" FORCE)

mark_as_advanced(BUILD_SHARED_LIBS)

set(AKANTU_CXX_STANDARD 17 CACHE INTERNAL "")

if(NOT AKANTU_TARGETS_EXPORT)
  set(AKANTU_TARGETS_EXPORT AkantuTargets)
endif()

include(CMakeVersionGenerator)
include(CMakePackagesSystem)
include(CMakeFlagsHandling)

include(AkantuPackagesSystem)
include(AkantuMacros)
include(AkantuCleaning)

#cmake_activate_debug_message()

include(GNUInstallDirs)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

#===============================================================================
# Version Number
#===============================================================================
# AKANTU version number.
define_project_version()

#===============================================================================
# Options
#===============================================================================
option(AKANTU_EXAMPLES "Activate examples" OFF)
option(AKANTU_TESTS "Activate tests" OFF)
option(AKANTU_SHARED "Build Akantu as a shared library" ON)
option(AKANTU_POSITION_INDEPENDENT "Build with -fPIC when static" ON)

option(AKANTU_RUN_IN_DOCKER "Set the approriate flage tu run in docker" OFF)
set(AKANTU_PREFERRED_PYTHON_VERSION 3
  CACHE STRING "Preferred version for python related things")
mark_as_advanced(
  AKANTU_PREFERRED_PYTHON_VERSION
  AKANTU_RUN_IN_DOCKER
  AKANTU_POSITION_INDEPENDENT
  AKANTU_SHARED
  )

if (AKANTU_SHARED)
  set(BUILD_SHARED_LIBS ON CACHE BOOL "Build shared libraries.")
else()
  set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build shared libraries.")
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

include(AkantuExtraCompilationProfiles)
#===============================================================================
# Dependencies
#===============================================================================
declare_akantu_types()

package_list_packages(${PROJECT_SOURCE_DIR}/packages
  EXTRA_PACKAGES_FOLDER ${PROJECT_SOURCE_DIR}/extra_packages
  NO_AUTO_COMPILE_FLAGS)

#===============================================================================
# Akantu library
#===============================================================================
if (NOT AKANTU_BYPASS_AKANTU_TARGET)
  add_subdirectory(src)
else()
  find_package(Akantu REQUIRED)
  if (Akantu_FOUND)
    get_target_property(_lib akantu INTERFACE_LINK_LIBRARIES)
    message(STATUS "Found Akantu: ${_lib} (found version ${AKANTU_VERSION})")
  endif()
endif()

#===============================================================================
# Documentation
#===============================================================================
if(AKANTU_DOCUMENTATION OR AKANTU_DOCUMENTATION_MANUAL)
  file(GLOB_RECURSE AKANTU_EXAMPLES_DOC_SRCS CONFIGURE_DEPENDS "examples/*.rst")
  add_subdirectory(doc)
else()
  set(AKANTU_DOC_EXCLUDE_FILES "${PROJECT_SOURCE_DIR}/doc/manual" CACHE INTERNAL "")
endif()

#===============================================================================
# Python interface
#===============================================================================
package_is_activated(python_interface _python_act)
if(_python_act)
  include(AkantuNeedPybind11)
  if(IS_ABSOLUTE "${CMAKE_INSTALL_PREFIX}")
    set(AKANTU_PYTHON_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
  else()
    set(AKANTU_PYTHON_INSTALL_PREFIX "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_PREFIX}")
  endif()
  add_subdirectory(python)
endif()

#===============================================================================
# Examples and tests
#===============================================================================
include(AkantuTestsMacros)
include(AkantuExampleMacros)

if(AKANTU_TESTS)
  include(AkantuNeedPybind11)
  option(AKANTU_BUILD_ALL_TESTS "Build all tests" ON)
  find_package(GMSH REQUIRED)
endif()

# tests
add_test_tree(test)

if(AKANTU_EXAMPLES)
  if(AKANTU_TESTS)
    option(AKANTU_TEST_EXAMPLES "Run the examples" ON)
  endif()
  find_package(GMSH REQUIRED)
  add_subdirectory(examples)
endif()

#===============================================================================
# Install and Packaging
#===============================================================================
if (NOT AKANTU_BYPASS_AKANTU_TARGET)
  include(AkantuInstall)

  option(AKANTU_DISABLE_CPACK
    "This option commands the generation of extra info for the \"make package\" target" ON)
  mark_as_advanced(AKANTU_DISABLE_CPACK)
  if(NOT AKANTU_DISABLE_CPACK)
    include(AkantuCPack)
  endif()
endif()
