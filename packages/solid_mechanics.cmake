#===============================================================================
# Copyright (©) 2017-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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


package_declare(solid_mechanics DEFAULT ON
  DESCRIPTION "Solid mechanics model"
  DEPENDS core
  )

package_declare_sources(solid_mechanics
  model/solid_mechanics/material.cc
  model/solid_mechanics/material.hh
  model/solid_mechanics/material_inline_impl.hh
  model/solid_mechanics/material_selector.hh
  model/solid_mechanics/material_selector_tmpl.hh
  model/solid_mechanics/solid_mechanics_model.cc
  model/solid_mechanics/solid_mechanics_model.hh
  model/solid_mechanics/solid_mechanics_model_io.cc
  model/solid_mechanics/solid_mechanics_model_mass.cc
  model/solid_mechanics/solid_mechanics_model_event_handler.hh
  model/solid_mechanics/materials/plane_stress_toolbox.hh
  model/solid_mechanics/materials/plane_stress_toolbox_tmpl.hh

  model/solid_mechanics/materials/material_elastic.cc
  model/solid_mechanics/materials/material_elastic.hh
  model/solid_mechanics/materials/material_elastic_inline_impl.hh
  model/solid_mechanics/materials/material_thermal.cc
  model/solid_mechanics/materials/material_thermal.hh
  model/solid_mechanics/materials/material_elastic_linear_anisotropic.cc
  model/solid_mechanics/materials/material_elastic_linear_anisotropic.hh
  model/solid_mechanics/materials/material_elastic_linear_anisotropic_inline_impl.hh
  model/solid_mechanics/materials/material_elastic_orthotropic.cc
  model/solid_mechanics/materials/material_elastic_orthotropic.hh
  model/solid_mechanics/materials/material_damage/material_anisotropic_damage.hh
  model/solid_mechanics/materials/material_damage/material_anisotropic_damage.cc
  model/solid_mechanics/materials/material_damage/material_anisotropic_damage_tmpl.hh
  model/solid_mechanics/materials/material_damage/material_damage.hh
  model/solid_mechanics/materials/material_damage/material_damage_tmpl.hh
  model/solid_mechanics/materials/material_damage/material_marigo.cc
  model/solid_mechanics/materials/material_damage/material_marigo.hh
  model/solid_mechanics/materials/material_damage/material_marigo_inline_impl.hh
  model/solid_mechanics/materials/material_damage/material_mazars.cc
  model/solid_mechanics/materials/material_damage/material_mazars.hh
  model/solid_mechanics/materials/material_damage/material_mazars_inline_impl.hh
  model/solid_mechanics/materials/material_finite_deformation/material_neohookean.cc
  model/solid_mechanics/materials/material_finite_deformation/material_neohookean.hh
  model/solid_mechanics/materials/material_finite_deformation/material_neohookean_inline_impl.hh
  model/solid_mechanics/materials/material_plastic/material_plastic.cc
  model/solid_mechanics/materials/material_plastic/material_plastic.hh
  model/solid_mechanics/materials/material_plastic/material_plastic_inline_impl.hh
  model/solid_mechanics/materials/material_plastic/material_drucker_prager.cc
  model/solid_mechanics/materials/material_plastic/material_drucker_prager.hh
  model/solid_mechanics/materials/material_plastic/material_drucker_prager_inline_impl.hh
  model/solid_mechanics/materials/material_plastic/material_linear_isotropic_hardening.cc
  model/solid_mechanics/materials/material_plastic/material_linear_isotropic_hardening.hh
  model/solid_mechanics/materials/material_plastic/material_linear_isotropic_hardening_inline_impl.hh
  model/solid_mechanics/materials/material_damage/material_von_mises_mazars.cc
  model/solid_mechanics/materials/material_damage/material_von_mises_mazars.hh
  model/solid_mechanics/materials/material_damage/material_von_mises_mazars_inline_impl.hh
  model/solid_mechanics/materials/material_viscoelastic/material_standard_linear_solid_deviatoric.cc
  model/solid_mechanics/materials/material_viscoelastic/material_standard_linear_solid_deviatoric.hh
  model/solid_mechanics/materials/material_viscoelastic/material_viscoelastic_maxwell.cc
  model/solid_mechanics/materials/material_viscoelastic/material_viscoelastic_maxwell.hh

  model/solid_mechanics/materials/material_non_local.hh
  )
