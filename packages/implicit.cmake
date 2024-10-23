#===============================================================================
# Copyright (©) 2012-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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


package_declare(implicit META
  DESCRIPTION "Add support for implicit time scheme")


set(AKANTU_IMPLICIT_SOLVER "Mumps+PETSc"
  CACHE STRING "Solver activated in Akantu")

set_property(CACHE AKANTU_IMPLICIT_SOLVER PROPERTY STRINGS
  Eigen
  Mumps
  PETSc
  Mumps+PETSc
)

package_declare_sources(implicit
  model/common/non_linear_solver/non_linear_solver_newton_raphson.cc
  model/common/non_linear_solver/non_linear_solver_newton_raphson.hh
)

package_is_activated(parallel _is_parallel)

if(_is_parallel AND AKANTU_IMPLICIT_SOLVER MATCHES "Eigen")
  message(WARNING "The Eigen solver does not work in parallel")
endif()

if(AKANTU_IMPLICIT_SOLVER MATCHES "Mumps")
  package_add_dependencies(implicit PRIVATE Mumps)
else()
  package_remove_dependencies(implicit Mumps)
  set(AKANTU_USE_MUMPS OFF CACHE BOOL "" FORCE)
endif()

if(AKANTU_IMPLICIT_SOLVER MATCHES "PETSc")
  package_add_dependencies(implicit
    PRIVATE PETSc)
else()
  package_remove_dependency(implicit PETSc)
  set(AKANTU_USE_PETSC OFF CACHE BOOL "" FORCE)
endif()
