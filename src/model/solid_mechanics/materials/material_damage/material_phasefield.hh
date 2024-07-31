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
#include "aka_common.hh"
#include "material.hh"
#include "material_damage.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MATERIAL_PHASEFIELD_HH_
#define AKANTU_MATERIAL_PHASEFIELD_HH_

namespace akantu {

template <Int dim> class MaterialPhaseField : public MaterialDamage<dim> {
  using Parent = MaterialDamage<dim>;
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MaterialPhaseField(SolidMechanicsModel & model, const ID & id = "");

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// constitutive law for all element of a type
  void computeStress(ElementType el_type,
                     GhostType ghost_type = _not_ghost) override;

  /// compute the tangent stiffness matrix for an element type
  void computeTangentModuli(ElementType el_type, Array<Real> & tangent_matrix,
                            GhostType ghost_type = _not_ghost) override;

  /// get mass density degraded by damage
  Vector<Real> getRho(const Element & element) const override;

  bool hasMassMatrixChanged() override { return true; };

  /* ------------------------------------------------------------------------ */
  decltype(auto) getArguments(ElementType el_type,
                              GhostType ghost_type = _not_ghost) {
    return zip_append(Parent::getArguments(el_type, ghost_type),
                      "effective_damage"_n = make_view(
                          this->effective_damage(el_type, ghost_type)));
  }

  decltype(auto) getArgumentsTangent(Array<Real> & tangent_matrix,
                                     ElementType el_type,
                                     GhostType ghost_type) {
    return zip_append(
        Parent::getArgumentsTangent(tangent_matrix, el_type, ghost_type),
        "effective_damage"_n =
            make_view(this->effective_damage(el_type, ghost_type)));
  }

protected:
  /// constitutive law for a given quadrature point
  template <class Args> inline void computeStressOnQuad(Args && args);

  /// compute the tangent stiffness matrix for a given quadrature point
  template <class Args> inline void computeTangentModuliOnQuad(Args && args);

  /// Compute the effective damage
  void computeEffectiveDamage(ElementType el_type,
                              GhostType ghost_type = _not_ghost);

  template <class Args> inline void computeEffectiveDamageOnQuad(Args && args);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  Real eta{0.};

  /// Phasefield isotropic
  bool is_hybrid;

  // effective damage to conserve stiffness in compression
  InternalField<Real> & effective_damage;
};

} // namespace akantu

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "material_phasefield_inline_impl.hh"

#endif /* AKANTU_MATERIAL_PHASEFIELD_HH_ */
