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


set(_MUMPS_COMPONENTS "sequential" "parallel" "double" "float" "complex_double" "complex_float")

if(NOT Mumps_FIND_COMPONENTS)
  set(Mumps_FIND_COMPONENTS "parallel" "double" "float" "complex_double" "complex_float")
endif()

#===============================================================================
enable_language(Fortran)

option(MUMPS_DETECT_DEBUG "Helps to debug mumps detection problems" OFF)
mark_as_advanced(MUMPS_DETECT_DEBUG)

set(MUMPS_PRECISIONS)
set(MUMPS_PLAT)

foreach(_comp ${Mumps_FIND_COMPONENTS})
  if("${_comp}" STREQUAL "sequential")
    set(MUMPS_PLAT _seq) #default plat on debian based distribution
  endif()

  if("${_comp}" STREQUAL "float")
    list(APPEND MUMPS_PRECISIONS s)
  endif()
  if("${_comp}" STREQUAL "double")
    list(APPEND MUMPS_PRECISIONS d)
  endif()
  if("${_comp}" STREQUAL "complex_float")
    list(APPEND MUMPS_PRECISIONS c)
  endif()
  if("${_comp}" STREQUAL "complex_double")
    list(APPEND MUMPS_PRECISIONS z)
  endif()
endforeach()

if(NOT MUMPS_PRECISIONS)
  set(MUMPS_PRECISIONS s d c z)
endif()

list(GET MUMPS_PRECISIONS 0 _first_precision)

string(TOUPPER "${_first_precision}" _u_first_precision)

find_path(MUMPS_INCLUDE_DIR ${_first_precision}mumps_c.h
  PATHS "${MUMPS_DIR}" ENV MUMPS_DIR
  PATH_SUFFIXES include
  )
mark_as_advanced(MUMPS_INCLUDE_DIR)

set(_mumps_required_vars)
foreach(_precision ${MUMPS_PRECISIONS})
  string(TOUPPER "${_precision}" _u_precision)

  if(DEFINED MUMPS_LIBRARY_${_u_precision}MUMPS AND
      (NOT "${Mumps_FIND_COMPONENTS}" STREQUAL "${Mumps_FIND_COMPONENTS_SAVE}"))
    set(MUMPS_LIBRARY_${_u_precision}MUMPS NOTFOUND CACHE PATH "" FORCE)
  endif()

  
  find_library(MUMPS_LIBRARY_${_u_precision}MUMPS
    NAMES ${_precision}mumps${MUMPS_PLAT} ${_precision}mumps
    PATHS "${MUMPS_DIR}" ENV MUMPS_DIR
    PATH_SUFFIXES lib
    )
  mark_as_advanced(MUMPS_LIBRARY_${_u_precision}MUMPS)
  list(APPEND _mumps_required_vars MUMPS_LIBRARY_${_u_precision}MUMPS)

  list(APPEND MUMPS_LIBRARIES_ALL ${MUMPS_LIBRARY_${_u_precision}MUMPS})
endforeach()

if(MUMPS_LIBRARY_${_u_first_precision}MUMPS MATCHES ".*${_first_precision}mumps.*${CMAKE_STATIC_LIBRARY_SUFFIX}")
  # Assuming mumps was compiled as a static library
  set(MUMPS_LIBRARY_TYPE STATIC CACHE INTERNAL "" FORCE)
else()
  set(MUMPS_LIBRARY_TYPE SHARED CACHE INTERNAL "" FORCE)
endif()


macro(find_mpiseq)
  find_library(MUMPS_LIBRARY_MPISEQ mpiseq${MUMPS_PLAT} mpiseq
    PATHS "${MUMPS_DIR}" ENV MUMPS_DIR
    PATH_SUFFIXES lib
    )
  if (NOT MUMPS_LIBRARY_MPISEQ)
    message("Could not find libmpiseq for sequential version of MUMPS, was "
      "MUMPS compiled in sequential ?")
  endif()
  set(${_libs} ${MUMPS_LIBRARY_MPISEQ} PARENT_SCOPE)
  mark_as_advanced(MUMPS_LIBRARY_MPISEQ)
endmacro()

macro(debug_message)
  if(MUMPS_DETECT_DEBUG)
    message(${ARGN})
  endif()
endmacro()

function(mumps_add_dependency _pdep _libs _incs)
  string(TOUPPER ${_pdep} _u_pdep)
  if(_pdep STREQUAL "mumps_common")
    find_library(MUMPS_LIBRARY_COMMON
      NAMES mumps_common${MUMPS_PLAT} mumps_common
      PATHS "${MUMPS_DIR}" ENV MUMPS_DIR
      PATH_SUFFIXES lib
      )
    set(${_libs} ${MUMPS_LIBRARY_COMMON} PARENT_SCOPE)
    mark_as_advanced(MUMPS_LIBRARY_COMMON)
  elseif(_pdep STREQUAL "pord")
    find_library(MUMPS_LIBRARY_PORD
      NAMES pord${MUMPS_PLAT} pord
      PATHS "${MUMPS_DIR}" ENV MUMPS_DIR
      PATH_SUFFIXES lib
      )
    set(${_libs} ${MUMPS_LIBRARY_PORD} PARENT_SCOPE)
    mark_as_advanced(MUMPS_LIBRARY_PORD)
  elseif(_pdep MATCHES "Scotch")
    find_package(Scotch REQUIRED ${ARGN} QUIET)
    if(ARGN)
      list(GET ARGN 1 _comp)
      string(TOUPPER ${_comp} _u_comp)
      set(${_libs} ${SCOTCH_LIBRARY_${_u_comp}} PARENT_SCOPE)
    else()
      set(${_libs} ${${_u_pdep}_LIBRARIES} PARENT_SCOPE)
    endif()
  elseif(_pdep MATCHES "MPI")
    if(MUMPS_PLAT STREQUAL "_seq")
      find_mpiseq()
    else()
      if(NOT CMAKE_Fortran_COMPILER_LOADED)
        enable_language(Fortran)
      endif()
      find_package(MPI REQUIRED QUIET
        COMPONENTS C Fortran)
      set(${_libs} ${MPI_C_LIBRARIES} ${MPI_Fortran_LIBRARIES} PARENT_SCOPE)
      set(${_incs}
        ${MPI_C_INCLUDE_PATH} # deprecated
        ${MPI_C_INCLUDE_DIRS}
        ${MPI_Fortran_INCLUDE_PATH} # deprecated
        ${MPI_Fortran_INCLUDE_DIRS}
        PARENT_SCOPE)
    endif()
  elseif(_pdep MATCHES "Threads")
    find_package(Threads REQUIRED QUIET)
    set(${_libs} Threads::Threads PARENT_SCOPE)
    set(_target_to_import Threads PARENT_SCOPE)
  elseif(_pdep MATCHES "OpenMP")
    find_package(OpenMP REQUIRED QUIET)
    set(${_libs} OpenMP::OpenMP_C PARENT_SCOPE)
    set(_target_to_import OpenMP PARENT_SCOPE)
  elseif(_pdep MATCHES "Math")
    set(${_libs} m PARENT_SCOPE)
  elseif(_pdep MATCHES "ScaLAPACK")
    if(MUMPS_PLAT STREQUAL "_seq")
      # ScaLAPACK symbols are in mpiseq form 5.20+
      find_mpiseq()
    else()
      find_package(ScaLAPACK REQUIRED QUIET)
      set(${_libs} ${SCALAPACK_LIBRARIES} PARENT_SCOPE)
    endif()
  elseif(_pdep MATCHES "gfortran")
    if(NOT CMAKE_Fortran_COMPILER_LOADED)
      enable_language(Fortran)
    endif()
    set(${_libs} gfortran PARENT_SCOPE)
  else()
    find_package(${_pdep} REQUIRED QUIET)
    set(${_libs} ${${_u_pdep}_LIBRARIES} ${${_u_pdep}_LIBRARY} PARENT_SCOPE)
  endif()
endfunction()

function(mumps_find_dependencies)
  set(_libraries_all m ${MUMPS_LIBRARIES_ALL})
  set(_include_dirs ${MUMPS_INCLUDE_DIR})
  set(_targets_to_import)

  set(_mumps_test_dir "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}")
  file(READ ${CMAKE_CURRENT_LIST_DIR}/CheckFindMumps.c _output)
  file(WRITE "${_mumps_test_dir}/mumps_test_code.c"
    "#include <${_first_precision}mumps_c.h>
${_u_first_precision}MUMPS_STRUC_C id;

#define mumps_c ${_first_precision}mumps_c
#define Real ${_u_first_precision}MUMPS_REAL
")

  if(MUMPS_PLAT STREQUAL _seq)
    file(APPEND "${_mumps_test_dir}/mumps_test_code.c"
      "#define MUMPS_SEQ
")
  else()
    file(APPEND "${_mumps_test_dir}/mumps_test_code.c"
      "// #undef MUMPS_SEQ
")
  endif()

  file(APPEND "${_mumps_test_dir}/mumps_test_code.c" "${_output}")

  #===============================================================================
  # ADD here the symbols needed to compile
  set(_mumps_dep_compile_MPI mpi.h)
  # ADD here the symbols needed to link
  set(_mumps_dep_link_MPI mpi_send mpi_type_free mpi_allreduce MPI_)
  set(_mumps_dep_link_BLAS ${_first_precision}gemm)
  set(_mumps_dep_link_ScaLAPACK numroc)
  set(_mumps_dep_link_LAPACK ilaenv)
  set(_mumps_dep_link_Scotch SCOTCH_graphInit scotchfstratexit)
  set(_mumps_dep_link_Scotch_ptscotch scotchfdgraphexit)
  set(_mumps_dep_link_Scotch_esmumps esmumps)
  set(_mumps_dep_link_mumps_common mumps_abort mumps_get_perlu)
  set(_mumps_dep_link_pord SPACE_ordering)
  set(_mumps_dep_link_METIS metis_nodend)
  set(_mumps_dep_link_Threads pthread_create)
  set(_mumps_dep_link_OpenMP omp_ GOMP_loop_end_nowait)
  set(_mumps_dep_link_gfortran gfortran)
  set(_mumps_dep_link_Math lround)
  set(_mumps_dep_link_ParMETIS ParMETIS_V3_NodeND)

  # ADD here the symbols needed to run
  set(_mumps_dep_run_mumps_common mumps_fac_descband)
  set(_mumps_dep_run_MPI mpi_bcast)
  set(_mumps_dep_run_ScaLAPACK idamax numroc)
  set(_mumps_dep_run_Scotch_ptscotch scotchfdgraphbuild)

  set(_mumps_dep_comp_Scotch_ptscotch COMPONENTS ptscotch)
  set(_mumps_dep_comp_Scotch_esmumps COMPONENTS esmumps)

  set(_mumps_potential_dependencies
    mumps_common pord
    MPI Threads OpenMP
    BLAS LAPACK ScaLAPACK
    Scotch Scotch_ptscotch Scotch_esmumps
    METIS ParMETIS
    gfortran Math)
  #===============================================================================

  set(_retry_try_run TRUE)
  set(_retry_count 0)

  set(ENV{LANG} C) # Trying to get the error messages in English


  # trying only as long as we add dependencies to avoid infinite loop in case of
  # an unknown dependency
  while (_retry_try_run AND _retry_count LESS 100)
    try_run(_mumps_run _mumps_compiles
      "${_mumps_test_dir}"
      "${_mumps_test_dir}/mumps_test_code.c"
      CMAKE_FLAGS "-DINCLUDE_DIRECTORIES:STRING=${_include_dirs}"
      LINK_LIBRARIES ${_libraries_all} ${_libraries_all}
      RUN_OUTPUT_VARIABLE _run
      COMPILE_OUTPUT_VARIABLE _out)

    set(_retry_compile FALSE)

    debug_message("COMPILATION outputs: \n${_out} \n RUN OUTPUT \n${_run}")

    if(_mumps_compiles AND NOT (_mumps_run STREQUAL "FAILED_TO_RUN"))
      break()
    endif()

    # On macOS with homebrew, MUMPS might be compiled with MPI even if we want sequential
    # The test might fail to run but compilation works, so we accept that
    if(APPLE AND _mumps_compiles AND _retry_count GREATER 5)
      message(STATUS "MUMPS test compiled successfully on macOS, accepting despite runtime issues")
      break()
    endif()

    if(_retry_count EQUAL 0 AND
        (NOT _mumps_compiles OR _mumps_run STREQUAL "FAILED_TO_RUN"))
      message(STATUS "Searching for MUMPS link dependencies")
    endif()
    
    foreach(_pdep ${_mumps_potential_dependencies})
      set(_libs)
      set(_incs)
      set(_add_pdep FALSE)
      
      debug_message("Trying to add: ${_pdep} as a dependency")
      
      if (NOT _mumps_compiles)
        if(DEFINED _mumps_dep_link_${_pdep})
          foreach (_link_dep ${_mumps_dep_link_${_pdep}})
            debug_message(" - test ${_link_dep}")
            if(_out MATCHES "undefined reference to.*${_link_dep}" OR
                _out MATCHES "${_link_dep}.*referenced from")
              set(_add_pdep TRUE)
              debug_message(" - ${_pdep} is a link dependency")
            endif()
          endforeach()
        endif()
        if (DEFINED _mumps_dep_compile_${_pdep})
          foreach (_compile_dep ${_mumps_dep_compile_${_pdep}})
            if(_out MATCHES "${_compile_dep}.*(No such file|file not found)")
              set(_add_pdep TRUE)
              debug_message(" - ${_pdep} is a compile dependency")
            endif()
          endforeach()
        endif()
      elseif(_mumps_run STREQUAL "FAILED_TO_RUN" AND
          DEFINED _mumps_dep_run_${_pdep})
        foreach(_run_dep ${_mumps_dep_run_${_pdep}})
          if(_run MATCHES "${_run_dep}")
            set(_add_pdep TRUE)
            debug_message(" - ${_pdep} is a run dependency")
          endif()
        endforeach()
      endif()

      if(_add_pdep)
        set(_target_to_import)
        mumps_add_dependency(${_pdep} _libs _incs ${_mumps_dep_comp_${_pdep}})
        debug_message(" - Found: ${_pdep} (${_libs})")

        if(_target_to_import)
          list(APPEND _targets_to_import ${_target_to_import})
        endif()

        if(NOT _libs)
          message(FATAL_ERROR "MUMPS depends on ${_pdep} but no libraries where found")
        else()
          message(STATUS "  Found MUMPS dependency ${_pdep} (${_libs})")
        endif()
        list(APPEND _libraries_all ${_libs})
        if(_incs)
          list(APPEND _include_dirs ${_incs})
        endif()
        set(_retry_try_run TRUE)
      endif()
    endforeach()

    math(EXPR _retry_count "${_retry_count} + 1")
  endwhile()

  if(_retry_count GREATER 10)
    message(FATAL_ERROR "Do not know what to do to link with mumps on your system, I give up!"
      "Last compilation outputs: \n${_out} \n And last run output \n${_run}")
  endif()

  if(APPLE)
    # in doubt add some stuff because mumps was perhaps badly compiled
    mumps_add_dependency(pord _libs _incs)
    list(APPEND _libraries_all ${_libs})
    mumps_add_dependency(mumps_common _libs _incs)
    list(APPEND _libraries_all ${_libs})
    mumps_add_dependency(BLAS _libs _incs)
    list(APPEND _libraries_all ${_libs})
    mumps_add_dependency(LAPACK _libs _incs)
    list(APPEND _libraries_all ${_libs})
    mumps_add_dependency(ScaLAPACK _libs _incs)
    list(APPEND _libraries_all ${_libs})
    mumps_add_dependency(gfortran _libs _incs)
    list(APPEND _libraries_all ${_libs})
  endif()

  set(MUMPS_LIBRARIES_ALL ${_libraries_all} PARENT_SCOPE)

  if(_targets_to_import)
    set(MUMPS_TARGETS_TO_IMPORT ${_targets_to_import} CACHE INTERNAL "" FORCE)
  endif()
endfunction()

mumps_find_dependencies()

set(MUMPS_LIBRARIES ${MUMPS_LIBRARIES_ALL} CACHE INTERNAL "" FORCE)
#===============================================================================
include(FindPackageHandleStandardArgs)
if(MUMPS_INCLUDE_DIR)
  file(STRINGS ${MUMPS_INCLUDE_DIR}/dmumps_c.h _versions
    REGEX "^#define MUMPS_VERSION .*")
  foreach(_ver ${_versions})
    string(REGEX MATCH "MUMPS_VERSION *\"([0-9.]+)\"" _tmp "${_ver}")
    set(_mumps_VERSION ${CMAKE_MATCH_1})
  endforeach()
  set(MUMPS_VERSION "${_mumps_VERSION}" CACHE INTERNAL "")
endif()

find_package_handle_standard_args(Mumps
  REQUIRED_VARS ${_mumps_required_vars}
  MUMPS_INCLUDE_DIR
  VERSION_VAR MUMPS_VERSION
)

if(Mumps_FOUND)
  set(Mumps_FIND_COMPONENTS_SAVE "${Mumps_FIND_COMPONENTS}" CACHE INTERNAL "")
endif()
