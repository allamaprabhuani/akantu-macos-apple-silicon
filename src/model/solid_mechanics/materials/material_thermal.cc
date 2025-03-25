/**
 * Copyright (©) 2010-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * This file is part of Akantu
 *
 * Akantu is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 */

/* -------------------------------------------------------------------------- */
#include "material_thermal.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int dim>
MaterialThermal<dim>::MaterialThermal(SolidMechanicsModel & model,
                                      const ID & id, const ID & fe_engine_id)
    : Parent(model, id, fe_engine_id),
      delta_T(this->registerInternal("delta_T", 1)),
      sigma_th(this->registerInternal("sigma_th", 1)),
      epsilon_th(this->registerInternal("epsilon_th", 1)) {
  sigma_th.initializeHistory();
  epsilon_th.initializeHistory();

  this->registerParam("E", E, Real(0.), _pat_parsable | _pat_modifiable,
                      "Young's modulus");
  this->registerParam("nu", nu, Real(0.5), _pat_parsable | _pat_modifiable,
                      "Poisson's ratio");
  this->registerParam("alpha", alpha, Real(0.), _pat_parsable | _pat_modifiable,
                      "Thermal expansion coefficient");
  this->registerParam("delta_T", delta_T, _pat_parsable | _pat_modifiable,
                      "Uniform temperature field");
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialThermal<dim>::computeStress(ElementType el_type,
                                         GhostType ghost_type) {
  auto && arguments = getArguments(el_type, ghost_type);
  for (auto && args : arguments) {
    computeStressOnQuad(args);
  }
}

/* -------------------------------------------------------------------------- */
template <>
void MaterialThermal<2>::computeStress(ElementType el_type,
                                       GhostType ghost_type) {
  auto && arguments = getArguments(el_type, ghost_type);
  if (this->plane_stress) {
    for (auto && args : arguments) {
      computeStressOnQuadPlaneStress(args);
    }
  } else {
    for (auto && args : arguments) {
      computeStressOnQuad(args);
    }
  }
}

/* -------------------------------------------------------------------------- */
template class MaterialThermal<1>;
template class MaterialThermal<2>;
template class MaterialThermal<3>;

} // namespace akantu
