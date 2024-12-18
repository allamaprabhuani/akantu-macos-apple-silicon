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


register_test(test_phasefield_selector
  SOURCES test_phasefield_selector.cc
  FILES_TO_COPY phasefield_selector.dat phasefield_selector.msh
  PACKAGE phase_field
  )

register_test(test_phase_solid_coupling
  SOURCES test_phase_solid_coupling.cc
  FILES_TO_COPY material_coupling.dat test_one_element.msh
  PACKAGE phase_field
  )

register_test(test_phase_field_anisotropic
  SOURCES test_phase_field_anisotropic.cc
  FILES_TO_COPY material_hybrid.dat test_one_element.msh material_penalization.dat
  PACKAGE phase_field
  )

register_test(test_phase_field_deviatoric_split
  SOURCES test_phase_field_deviatoric_split.cc
  FILES_TO_COPY material_deviatoric.dat test_one_element.msh
  PACKAGE phase_field
  )

register_test(test_phase_field_optimal
  SOURCES test_phase_field_optimal.cc
  FILES_TO_COPY material_optimal.dat test_optimal_damage.msh
  PACKAGE phase_field
  )

register_test(test_phase_field_optimal_linear
  SOURCES test_phase_field_optimal_linear.cc
  FILES_TO_COPY material_optimal_linear.dat test_optimal_damage.msh
  PACKAGE phase_field
  )

register_test(test_phase_solid_explicit
  SOURCES test_phase_solid_explicit.cc
  FILES_TO_COPY material_coupling.dat test_one_element.msh
  PACKAGE phase_field
  )

register_test(test_multi_material
  SOURCES test_multi_material.cc
  FILES_TO_COPY material_multiple.dat test_two_element.msh
  PACKAGE phase_field
  UNSTABLE
  )

register_test(test_penalization
  SOURCES test_penalization.cc
  FILES_TO_COPY material_penalization.dat test_four_elements.msh
  PACKAGE phase_field
  )

register_test(test_phase_field_linear
  SOURCES test_phase_field_linear.cc
  FILES_TO_COPY material_linear.dat test_one_element.msh
  PACKAGE phase_field
  )

register_test(test_phase_field_tao
  SOURCES test_phase_field_tao.cc
  FILES_TO_COPY material_linear_tao.dat test_one_element.msh
  PACKAGE phase_field
  )
