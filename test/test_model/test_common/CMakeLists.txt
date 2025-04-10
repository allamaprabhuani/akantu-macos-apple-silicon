#===============================================================================
# Copyright (©) 2019-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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


add_akantu_test(test_model_solver "Test for the solvers")
add_akantu_test(test_non_local_toolbox "Test of the functionalities in the non-local toolbox")

add_mesh(dof_manager_mesh mesh.geo 3 1)
register_gtest_sources(
  SOURCES test_dof_manager.cc
  PACKAGE core)

register_gtest_sources(
  SOURCES test_dof_mesh_distribute.cc
  PACKAGE parallel)

register_gtest_test(test_dof_manager DEPENDS dof_manager_mesh
  PARALLEL)
