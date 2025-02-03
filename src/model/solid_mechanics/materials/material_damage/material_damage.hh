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
#include "material_elastic.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MATERIAL_DAMAGE_HH_
#define AKANTU_MATERIAL_DAMAGE_HH_

namespace akantu {
template <Int dim, template <Int> class Parent = MaterialElastic>
class MaterialDamage : public Parent<dim> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MaterialDamage(SolidMechanicsModel & model, const ID & id = "");

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// compute the tangent stiffness matrix for an element type
  void computeTangentModuli(ElementType el_type, Array<Real> & tangent_matrix,
                            GhostType ghost_type = _not_ghost) override;

  auto hasStiffnessMatrixChanged() -> bool override { return true; }

protected:
  /// update the dissipated energy, must be called after the stress have been
  /// computed
  void updateEnergies(ElementType el_type) override;

  /// compute the tangent stiffness matrix for a given quadrature point
  template <class Args>
  inline void computeTangentModuliOnQuad(Args && arguments);

  auto getDissipatedEnergy() const -> Real;

  /* ------------------------------------------------------------------------ */
  /* DataAccessor inherited members                                           */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  decltype(auto) getArguments(ElementType el_type,
                              GhostType ghost_type = _not_ghost) {
    return zip_append(Parent<dim>::getArguments(el_type, ghost_type),
                      "damage"_n = this->damage(el_type, ghost_type),
                      "previous_damage"_n =
                          this->damage.previous(el_type, ghost_type));
  }

  decltype(auto) getArgumentsTangent(Array<Real> & tangent_matrix,
                                     ElementType el_type,
                                     GhostType ghost_type) {
    return zip_append(
        Parent<dim>::getArgumentsTangent(tangent_matrix, el_type, ghost_type),
        "damage"_n = this->damage(el_type, ghost_type),
        "previous_damage"_n = this->damage.previous(el_type, ghost_type));
  }

  Real getEnergy(const std::string & type) override;

  // AKANTU_GET_MACRO_AUTO_NOT_CONST(Damage, damage);
  // AKANTU_GET_MACRO_AUTO(Damage, damage);
  // InternalField<Real> getDamage() { return damage; }
  virtual inline InternalField<Real> & getDamage() { return damage; }
  inline InternalField<Real> & getDamage() const { return damage; }

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Damage, damage, Real);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// damage internal variable
  InternalField<Real> & damage;

  /// dissipated energy
  InternalField<Real> & dissipated_energy;

  /// contain the current value of @f$ \int_0^{\epsilon}\sigma(\omega)d\omega
  /// @f$ the dissipated energy
  InternalField<Real> & int_sigma;
};

} // namespace akantu

#include "material_damage_tmpl.hh"

#endif /* AKANTU_MATERIAL_DAMAGE_HH_ */
