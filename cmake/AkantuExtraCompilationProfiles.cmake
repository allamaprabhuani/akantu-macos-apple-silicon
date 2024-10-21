#===============================================================================
# Copyright (©) 2016-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
option (FORCE_COLORED_OUTPUT "Always produce ANSI-colored output (GNU/Clang only)." FALSE)
mark_as_advanced(FORCE_COLORED_OUTPUT)

set(AKANTU_CXX_EXTRA_FLAGS)
if(FORCE_COLORED_OUTPUT)
  if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    list(APPEND AKANTU_CXX_EXTRA_FLAGS "-fcolor-diagnostics")
  else()
    list(APPEND AKANTU_CXX_EXTRA_FLAGS "-fdiagnostics-color=always")
  endif()
endif()

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL "11")
  list(APPEND AKANTU_CXX_EXTRA_FLAGS "-flarge-source-files")
endif()

if (AKANTU_CXX_EXTRA_FLAGS)
  set(CMAKE_CXX_FLAGS_RELEASE
    "${CMAKE_CXX_FLAGS_RELEASE_INIT} -DAKANTU_NDEBUG ${AKANTU_CXX_EXTRA_FLAGS}"
    CACHE STRING "Flags used by the compiler during Release builds" FORCE)
  set(CMAKE_CXX_FLAGS_DEBUG
    "${CMAKE_CXX_FLAGS_DEBUG_INIT} -g3 -ggdb3 ${AKANTU_CXX_EXTRA_FLAGS}"
    CACHE STRING "Flags used by the compiler during Debug builds" FORCE)
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO
    "${CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT} -g3 -ggdb3 ${AKANTU_CXX_EXTRA_FLAGS}"
    CACHE STRING "Flags used by the compiler during RelWithDebInfo builds" FORCE)
endif()

function(declare_compilation_profile name)
  include(CMakeParseArguments)

  cmake_parse_arguments(_args
    "" "COMPILER;LINKER;DOC" "" ${ARGN})

  string(TOUPPER "${name}" _u_name)

  if(NOT _args_DOC)
    string(TOLOWER "${name}" _args_DOC)
  endif()

  if(NOT _args_COMPILER)
    message(FATAL_ERROR "declare_compilation_profile: you should at least give COMPILER flags")
  endif()

  if(NOT _args_LINKER)
    set(_args_LINKER ${_args_COMPILER})
  endif()

  foreach(_flag CXX C Fortran SHARED_LINKER EXE_LINKER)
    set(_stage "compiler")
    set(_flags ${_args_COMPILER})
    if(_stage MATCHES ".*LINKER")
      set(_stage "linker")
      set(_flags ${_args_LINKER})
    endif()
    set(CMAKE_${_flag}_FLAGS_${_u_name} ${_flags}
      CACHE STRING "Flags used by the ${_stage} during coverage builds" FORCE)
    mark_as_advanced(CMAKE_${_flag}_FLAGS_${_u_name})
  endforeach()
endfunction()


# Profiling
declare_compilation_profile(PROFILING
  COMPILER "-g -ggdb3 -pg -DNDEBUG -DAKANTU_NDEBUG -O3")

# Valgrind
declare_compilation_profile(VALGRIND
  COMPILER "-g -ggdb3 -DNDEBUG -DAKANTU_NDEBUG -O3")

# Coverage
declare_compilation_profile(COVERAGE
  COMPILER "-g -ggdb3 -DNDEBUG -DAKANTU_NDEBUG -O2 --coverage")

# Sanitize the code
if ((CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "5.2") OR
    CMAKE_CXX_COMPILER_ID STREQUAL "Clang")

  if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(_blacklist " -fsanitize-blacklist=${PROJECT_SOURCE_DIR}/cmake/sanitize-blacklist.txt")
  endif()

  declare_compilation_profile(SANITIZE
    COMPILER "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -fsanitize=address -fsanitize=leak -fsanitize=undefined -fno-omit-frame-pointer${_blacklist}")

  declare_compilation_profile(SANITIZEDEBUG
    COMPILER "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address -fsanitize=leak -fsanitize=undefined -fno-omit-frame-pointer${_blacklist}")
endif()

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  declare_compilation_profile(SANITIZEMEMORY
    COMPILER "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -fPIE -fsanitize=memory -fsanitize-memory-track-origins -fsanitize-recover=all -fno-omit-frame-pointer -fsanitize-blacklist=${PROJECT_SOURCE_DIR}/cmake/sanitize-blacklist.txt"
    DOC "\"sanitize memory\"")
endif()

string(TOLOWER "${CMAKE_BUILD_TYPE}" _cmake_build_type_lower)
if (_cmake_build_type_lower MATCHES "valgrind")
  find_program(VALGRIND_EXECUTABLE valgrind)
endif()

find_program(CCACHE_EXECUTABLE ccache)
if(CCACHE_EXECUTABLE)
  option(AKANTU_USE_CCACHE "Use ccache if available to build akantu" ON)
  mark_as_advanced(AKANTU_USE_CCACHE)
endif()
option(AKANTU_SPLIT_DWARF "Split the debug symbols in separate DWO files" OFF)
mark_as_advanced(AKANTU_SPLIT_DWARF)
if (CCACHE_EXECUTABLE AND AKANTU_USE_CCACHE)
  set(AKANTU_COMPILER_LAUNCHER "${CCACHE_EXECUTABLE}")
endif()

include(ProcessorCount)
ProcessorCount(N)

set_property(GLOBAL PROPERTY JOB_POOLS
  akantu_compile=${N}
  akantu_link=1
  akantu_link_test=${N}
  akantu_compile_heavy=2)
