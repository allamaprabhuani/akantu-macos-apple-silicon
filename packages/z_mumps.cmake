#===============================================================================
# Copyright (©) 2011-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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


package_declare(Mumps EXTERNAL
  DESCRIPTION "Add Mumps support in akantu"
  )

package_declare_sources(Mumps
  solver/sparse_solver_mumps.cc
  solver/sparse_solver_mumps.hh
  )

set(_mumps_float_type ${AKANTU_FLOAT_TYPE})

if(AKANTU_FLOAT_TYPE STREQUAL "float" OR
    AKANTU_FLOAT_TYPE STREQUAL "double")
  set(_mumps_components ${AKANTU_FLOAT_TYPE})
else()
  if(DEFINED AKANTU_FLOAT_TYPE)
    message(FATAL_ERROR "MUMPS does not support floating point type \"${AKANTU_FLOAT_TYPE}\"")
  endif()
endif()

package_get_option_name(parallel _par_option)
if(${_par_option})
  list(APPEND _mumps_components "parallel")

  package_set_package_system_dependency(Mumps deb libmumps)
  package_set_package_system_dependency(Mumps deb-src libmumps-dev)
  package_add_dependencies(Mumps PRIVATE parallel)
else()
  list(APPEND _mumps_components "sequential")
  package_set_package_system_dependency(Mumps deb libmumps-seq)
  package_set_package_system_dependency(Mumps deb-src libmumps-seq-dev)
  package_remove_dependencies(Mumps PRIVATE parallel)
endif()

package_set_find_package_extra_options(Mumps ARGS COMPONENTS "${_mumps_components}")
package_declare_extra_files_to_package(Mumps
  PROJECT
    cmake/Modules/FindMumps.cmake
  )
