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


package_declare(phase_field DEFAULT ON
  DEPENDS model_couplers
  DESCRIPTION "Use Phase Field package of Akantu")

package_declare_sources(phase_field
  model/phase_field/phasefield.cc
  model/phase_field/phasefield.hh
  model/phase_field/phasefield_inline_impl.hh
  model/phase_field/phasefield_selector.hh
  model/phase_field/phasefield_selector_tmpl.hh

  model/phase_field/phasefields/phasefield_quadratic.hh
  model/phase_field/phasefields/phasefield_quadratic.cc
  model/phase_field/phasefields/phasefield_quadratic_inline_impl.hh
  
  model/phase_field/phasefields/phasefield_linear.hh
  model/phase_field/phasefields/phasefield_linear.cc
  model/phase_field/phasefields/phasefield_linear_inline_impl.hh
  
  model/phase_field/phase_field_model.cc
  model/phase_field/phase_field_model.hh
  model/phase_field/phase_field_model_inline_impl.hh

  model/solid_mechanics/materials/material_damage/material_phasefield.cc
  model/solid_mechanics/materials/material_damage/material_phasefield.hh
  model/solid_mechanics/materials/material_damage/material_phasefield_inline_impl.hh
  model/model_couplers/coupler_solid_phasefield.hh
  model/model_couplers/coupler_solid_phasefield.cc

  model/phase_field/phase_field_element_filter.hh

  model/phase_field/energy_split.hh
  model/phase_field/energy_split.cc

  model/phase_field/energy_splits/volumetric_deviatoric_split.hh
  model/phase_field/energy_splits/volumetric_deviatoric_split_inline_impl.hh

  model/phase_field/energy_splits/no_energy_split.hh
  model/phase_field/energy_splits/no_energy_split_inline_impl.hh
  )
