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
#===============================================================================


#===============================================================================
# Package Management
#===============================================================================
package_get_all_source_files(
  AKANTU_LIBRARY_SRCS
  AKANTU_LIBRARY_PUBLIC_HDRS
  AKANTU_LIBRARY_PRIVATE_HDRS
  )

package_get_all_include_directories(
  AKANTU_LIBRARY_INCLUDE_DIRS
  )

package_get_all_external_informations(
  PRIVATE_INCLUDE AKANTU_PRIVATE_EXTERNAL_INCLUDE_DIR
  INTERFACE_INCLUDE AKANTU_INTERFACE_EXTERNAL_INCLUDE_DIR
  LIBRARIES AKANTU_EXTERNAL_LIBRARIES
  )

package_get_all_compilation_flags(CXX _cxx_flags)
set(AKANTU_EXTRA_CXX_FLAGS
  "${_cxx_flags}" CACHE STRING "Extra flags defined by loaded packages" FORCE)
mark_as_advanced(AKANTU_EXTRA_CXX_FLAGS)

foreach(src_ ${AKANTU_SPIRIT_SOURCES})
  set_property(SOURCE ${src_} PROPERTY COMPILE_FLAGS "-g0 -Werror")
endforeach()

list(APPEND AKANTU_LIBRARY_INCLUDE_DIRS "${CMAKE_CURRENT_BINARY_DIR}")
set(AKANTU_INCLUDE_DIRS
  ${CMAKE_CURRENT_BINARY_DIR} ${AKANTU_LIBRARY_INCLUDE_DIRS}
  CACHE INTERNAL "Internal include directories to link with Akantu as a subproject")

#===========================================================================
# configurations
#===========================================================================
configure_file(common/aka_config.hh.in
  "${CMAKE_CURRENT_BINARY_DIR}/aka_config.hh" @ONLY)

list(APPEND AKANTU_LIBRARY_PUBLIC_HDRS
  "${CMAKE_CURRENT_BINARY_DIR}/aka_config.hh")

configure_file(common/aka_config.cc.in
  "${CMAKE_CURRENT_BINARY_DIR}/aka_config.cc" @ONLY)

list(APPEND AKANTU_LIBRARY_SRCS
  "${CMAKE_CURRENT_BINARY_DIR}/aka_config.cc")

#===============================================================================
# Debug infos
#===============================================================================
set(AKANTU_GDB_DIR ${PROJECT_SOURCE_DIR}/cmake)
if(UNIX AND NOT APPLE)
  string(TOUPPER "${CMAKE_BUILD_TYPE}" _u_build_type)
  if(_u_build_type MATCHES "DEBUG" OR _u_build_type MATCHES "RELWITHDEBINFO")
    configure_file(${PROJECT_SOURCE_DIR}/cmake/libakantu-gdb.py.in
      "${PROJECT_BINARY_DIR}/libakantu-gdb.py"
      @ONLY)
    configure_file(${PROJECT_SOURCE_DIR}/cmake/akantu-debug.cc.in
      "${PROJECT_BINARY_DIR}/akantu-debug.cc" @ONLY)

    list(APPEND AKANTU_LIBRARY_SRCS ${PROJECT_BINARY_DIR}/akantu-debug.cc)
  endif()
else()
  find_program(GDB_EXECUTABLE gdb)

  if(GDB_EXECUTABLE)
    execute_process(COMMAND
      ${GDB_EXECUTABLE} --batch -x "${PROJECT_SOURCE_DIR}/cmake/gdb_python_path"
      OUTPUT_VARIABLE AKANTU_PYTHON_GDB_DIR
      ERROR_QUIET
      RESULT_VARIABLE _res)

    if(_res EQUAL 0 AND UNIX)
      set(GDB_USER_CONFIG $ENV{HOME}/.gdb/auto-load)
      file(MAKE_DIRECTORY ${GDB_USER_CONFIG})

      configure_file(${PROJECT_SOURCE_DIR}/cmake/libakantu-gdb.py.in
        "${GDB_USER_CONFIG}/${CMAKE_SHARED_LIBRARY_PREFIX}akantu${CMAKE_SHARED_LIBRARY_SUFFIX}.${AKANTU_VERSION}-gdb.py"
        @ONLY)
    endif()
  endif()
endif()

#===============================================================================
# GIT info
#===============================================================================

#include(AkantuBuildOptions)
#git_version_info(akantu_git_info ALL
#  OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/akantu_git_info.hh
#  )
#list(APPEND AKANTU_LIBRARY_SRCS ${CMAKE_CURRENT_BINARY_DIR}/akantu_git_info.hh)
#===============================================================================
# Library generation
#===============================================================================
add_library(akantu ${AKANTU_LIBRARY_SRCS})

target_include_directories(akantu
  PRIVATE   $<BUILD_INTERFACE:${AKANTU_INCLUDE_DIRS}>
  INTERFACE $<INSTALL_INTERFACE:include/akantu>
  )

# small trick for build includes in public
set_property(TARGET akantu APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES
  $<BUILD_INTERFACE:${AKANTU_INCLUDE_DIRS}>)

target_include_directories(akantu SYSTEM
  PUBLIC ${AKANTU_INTERFACE_EXTERNAL_INCLUDE_DIR}
  )

target_include_directories(akantu SYSTEM
  PRIVATE ${AKANTU_PRIVATE_EXTERNAL_INCLUDE_DIR}
  )

target_link_libraries(akantu PUBLIC ${AKANTU_EXTERNAL_LIBRARIES})

separate_arguments(_separated_cxx_flags UNIX_COMMAND ${_cxx_flags})
target_compile_options(akantu PUBLIC ${_separated_cxx_flags})

set_property(TARGET akantu PROPERTY JOB_POOL_COMPILE akantu_compile)
set_property(TARGET akantu PROPERTY JOB_POOL_LINK akantu_link)

if(AKANTU_LIBRARY_PROPERTIES)
  set_target_properties(akantu
    PROPERTIES ${AKANTU_LIBRARY_PROPERTIES} # this contains the version
    )
endif()

if(NOT AKANTU_SHARED)
  set_property(TARGET akantu PROPERTY
    POSITION_INDEPENDENT_CODE ${AKANTU_POSITION_INDEPENDENT})
endif()
  
if(AKANTU_LIBRARY_PUBLIC_HDRS)
  set_property(TARGET akantu PROPERTY PUBLIC_HEADER ${AKANTU_LIBRARY_PUBLIC_HDRS})
endif()

if(AKANTU_LIBRARY_PRIVATE_HDRS)
  set_property(TARGET akantu PROPERTY PRIVATE_HEADER ${AKANTU_LIBRARY_PRIVATE_HDRS})
endif()

if(AKANTU_COMPILER_LAUNCHER)
  set_property(TARGET akantu PROPERTY CXX_COMPILER_LAUNCHER "${AKANTU_COMPILER_LAUNCHER}")
endif()

set_property(TARGET akantu PROPERTY CXX_EXTENSIONS OFF)

target_compile_features(akantu PUBLIC cxx_std_${AKANTU_CXX_STANDARD})
set_target_properties(akantu
  PROPERTIES CXX_STANDARD ${AKANTU_CXX_STANDARD})

if (AKANTU_COMPILER_LAUNCHER)
  set_target_properties(akantu
    PROPERTIES CXX_COMPILER_LAUNCHER ${AKANTU_COMPILER_LAUNCHER})
endif()

package_get_all_extra_dependencies(_extra_target_dependencies)
if(_extra_target_dependencies)
  # This only adding todo: find a solution for when a dependency was add the is removed...
  add_dependencies(akantu ${_extra_target_dependencies})
endif()

if(AKANTU_SPLIT_DWARF)
  target_compile_options(akantu
    PRIVATE $<$<CONFIG:DEBUG>:-gsplit-dwarf -gdwarf-5> $<$<CONFIG:SANITIZEDEBUG>:-gsplit-dwarf -gdwarf-5>
  )
  target_link_options(akantu
    PRIVATE -fuse-ld=gold $<$<CONFIG:DEBUG>:-Wl,--gdb-index> $<$<CONFIG:SANITIZEDEBUG>:-Wl,--gdb-index>
  )
endif()

package_get_all_export_list(AKANTU_EXPORT_LIST)
list(APPEND AKANTU_EXPORT_LIST akantu)


# TODO separate public from private headers
install(TARGETS akantu
  EXPORT ${AKANTU_TARGETS_EXPORT}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
          COMPONENT Akantu_runtime
  #         NAMELINK_ONLY
  # LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  #         COMPONENT Akantu_development
  #         NAMELINK_SKIP Akantu_development
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
          COMPONENT Akantu_runtime
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
          COMPONENT Akantu_runtime
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/akantu/
          COMPONENT Akantu_development
  )

if("${AKANTU_TARGETS_EXPORT}" STREQUAL "AkantuTargets")
  install(EXPORT AkantuTargets
    DESTINATION lib/cmake/${PROJECT_NAME}
    COMPONENT dev)

  #Export for build tree
  export(EXPORT AkantuTargets
    FILE "${CMAKE_BINARY_DIR}/AkantuTargets.cmake")
  export(PACKAGE Akantu)
endif()

#===============================================================================
# Adding module names for debug
package_get_all_packages(_pkg_list)
foreach(_pkg ${_pkg_list})
  _package_get_real_name(${_pkg} _pkg_name)
  _package_get_source_files(${_pkg} _srcs _public_hdrs _private_hdrs)
  string(TOLOWER "${_pkg_name}" _l_package_name)
  set_property(SOURCE ${_srcs} ${_public_hdrs} ${_private_hdrs} 
    APPEND PROPERTY COMPILE_DEFINITIONS AKANTU_MODULE=${_l_package_name})
endforeach()

register_target_to_tidy(akantu)
register_tidy_all(${AKANTU_LIBRARY_SRCS})
register_code_to_format(
  ${AKANTU_LIBRARY_SRCS}
  ${AKANTU_LIBRARY_PUBLIC_HDRS}
  ${AKANTU_LIBRARY_PRIVATE_HDRS}
  )

cmake_policy(SET CMP0069 NEW)

include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_supported OUTPUT ipo_error LANGUAGES CXX)
if(ipo_supported)
  string(TOUPPER "${CMAKE_BUILD_TYPE}" _u_build_type)
  if(NOT _u_build_type MATCHES "DEBUG")
    message(STATUS "IPO/LTO enabled")
    set_property(TARGET akantu PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
  else()
    message(STATUS "IPO/LTO supported but deactivated in debug")
    set_property(TARGET akantu PROPERTY INTERPROCEDURAL_OPTIMIZATION FALSE)
  endif()
else()
    message(STATUS "IPO/LTO not supported: <${ipo_error}>")
endif()
