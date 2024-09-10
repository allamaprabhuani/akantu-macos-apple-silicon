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
#include "material_phasefield.hh"
#include "aka_common.hh"
#include "solid_mechanics_model.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int dim>
MaterialPhaseField<dim>::MaterialPhaseField(SolidMechanicsModel & model,
                                            const ID & id)
    : Parent(model, id),
      effective_damage(this->registerInternal("effective_damage", 1)) {
  this->registerParam("eta", eta, Real(0.), _pat_parsable,
                      "Residual rigidity factor");
  this->registerParam("is_hybrid", is_hybrid, false,
                      _pat_parsable | _pat_readable, "Use hybrid formulation");
  this->registerParam("degrade_mass", degrade_mass, false,
                      _pat_parsable | _pat_readable,
                      "Degrade mass with damage");
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialPhaseField<dim>::computeStress(ElementType el_type,
                                            GhostType ghost_type) {

  if (this->is_hybrid) {
    computeEffectiveDamage(el_type, ghost_type);

    for (auto && args : getArguments(el_type, ghost_type)) {
      auto && dam = args["effective_damage"_n];
      computeStressOnQuad(tuple::replace(args, "damage"_n = dam));
    }
  } else {
    for (auto && args : getArguments(el_type, ghost_type)) {
      computeStressOnQuad(args);
    }
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialPhaseField<dim>::computeTangentModuli(ElementType el_type,
                                                   Array<Real> & tangent_matrix,
                                                   GhostType ghost_type) {
  computeEffectiveDamage(el_type, ghost_type);

  if (this->is_hybrid) {
    computeEffectiveDamage(el_type, ghost_type);

    for (auto && args :
         getArgumentsTangent(tangent_matrix, el_type, ghost_type)) {
      auto && dam = args["effective_damage"_n];
      computeTangentModuliOnQuad(tuple::replace(args, "damage"_n = dam));
    }
  } else {
    for (auto && args :
         getArgumentsTangent(tangent_matrix, el_type, ghost_type)) {
      computeTangentModuliOnQuad(args);
    }
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialPhaseField<dim>::computeEffectiveDamage(ElementType el_type,
                                                     GhostType ghost_type) {

  for (auto && args : getArguments(el_type, ghost_type)) {
    computeEffectiveDamageOnQuad(args);
  }
}

/* -------------------------------------------------------------------------- */
template class MaterialPhaseField<1>;
template class MaterialPhaseField<2>;
template class MaterialPhaseField<3>;

const bool material_is_allocated_phasefield [[maybe_unused]] =
    instantiateMaterial<MaterialPhaseField>("phasefield");

} // namespace akantu
