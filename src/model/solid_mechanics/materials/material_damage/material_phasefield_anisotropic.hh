/**
 * @file   material_phasefield_anisotropic.hh
 *
 * @author Shad Durussel <shad.durussel@epfl.ch>
 *
 * @date creation: Mon Mar 27 2023
 * @date last modification: Mon Mar 27 2023
 *
 * @brief  Phasefield anisotropic damage law
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
 *
 */

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "energy_split.hh"
#include "material.hh"
#include "material_damage.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_PHASEFIELD_ANISOTROPIC_HH__
#define __AKANTU_MATERIAL_PHASEFIELD_ANISOTROPIC_HH__

namespace akantu {

template <Int dim>
class MaterialPhaseFieldAnisotropic : public MaterialDamage<dim> {
  using Parent = MaterialDamage<dim>;
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MaterialPhaseFieldAnisotropic(SolidMechanicsModel & model,
                                const ID & id = "");
  ~MaterialPhaseFieldAnisotropic() override = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void initMaterial() override;
  /// constitutive law for all element of a type
  void computeStress(ElementType el_type,
                     GhostType ghost_type = _not_ghost) override;

  /// compute the tangent stiffness matrix for an element type
  void computeTangentModuli(ElementType el_type, Array<Real> & tangent_matrix,
                            GhostType ghost_type = _not_ghost) override;

  /* ------------------------------------------------------------------------ */
  // decltype(auto) getArguments(ElementType el_type,
  //                             GhostType ghost_type = _not_ghost) {
  //   return zip_append(Parent::getArguments(el_type, ghost_type),
  //                     "effective_damage"_n = make_view(
  //                         this->effective_damage(el_type, ghost_type)));
  // }

  // decltype(auto) getArgumentsTangent(Array<Real> & tangent_matrix,
  //                                    ElementType el_type,
  //                                    GhostType ghost_type) {
  //   return zip_append(
  //       Parent::getArgumentsTangent(tangent_matrix, el_type, ghost_type),
  //       "effective_damage"_n =
  //           make_view(this->effective_damage(el_type, ghost_type)));
  // }

  /// get mass density degraded by damage
  Vector<Real> getRho(const Element & element) const override;

  bool hasMassMatrixChanged() override { return degrade_mass; };

  /* ------------------------------------------------------------------------ */
  /* DataAccessor inherited members                                           */
  /* ------------------------------------------------------------------------ */
public:
  [[nodiscard]] inline Int
  getNbData(const Array<Element> & elements,
            const SynchronizationTag & tag) const override;

  inline void packData(CommunicationBuffer & buffer,
                       const Array<Element> & elements,
                       const SynchronizationTag & tag) const override;

  inline void unpackData(CommunicationBuffer & buffer,
                         const Array<Element> & elements,
                         const SynchronizationTag & tag) override;

protected:
  /// constitutive law for a given quadrature point
  template <class Args> inline void computeStressOnQuad(Args && args);

  /// compute the tangent stiffness matrix for a given quadrature point
  template <class Args> inline void computeTangentModuliOnQuad(Args && args);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// Residual stiffness parameter
  Real eta;

  /// Phasefield isotropic
  bool is_isotropic;

  // Dimension considered in volumetric-deviatoric split
  Int dev_dim;

  bool degrade_mass;

  std::shared_ptr<EnergySplit> energy_split{nullptr};
};

} // namespace akantu

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "material_phasefield_anisotropic_inline_impl.hh"

#endif /* __AKANTU_MATERIAL_PHASEFIELD_HH__ */
