/**
 * @file   material_phasefield.hh
 *
 * @author Shad Durussel <shad.durussel@epfl.ch>
 *
 * @date creation: Mon Mar 27 2023
 * @date last modification: Mon Mar 27 2023
 *
 * @brief  Phasefield damage law
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
#include "no_energy_split.hh"
#include "volumetric_deviatoric_split.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_PHASEFIELD_HH__
#define __AKANTU_MATERIAL_PHASEFIELD_HH__

namespace akantu {

/* ------------------------------------------------------------------------ */

template <Int dim, template <Int> class EnergySplit_>
requires CanComputeSigmaAndTangent<dim, EnergySplit_>
class MaterialPhaseField : public MaterialDamage<dim> {
  using Parent = MaterialDamage<dim>;
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MaterialPhaseField(SolidMechanicsModel & model, const ID & id = "");
  ~MaterialPhaseField() override = default;

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
  /// get mass density degraded by damage
  void getRho(Ref<Vector<Real>> rhos, const Element & element) const override;

  bool hasMassMatrixChanged() override { return degrade_mass; };

  /* ------------------------------------------------------------------------ */
  /* DataAccessor inherited members                                           */
  /* ------------------------------------------------------------------------ */
public:
  [[nodiscard]] Int getNbData(const Array<Element> & elements,
                              const SynchronizationTag & tag) const override;

  void packData(CommunicationBuffer & buffer, const Array<Element> & elements,
                const SynchronizationTag & tag) const override;

  void unpackData(CommunicationBuffer & buffer, const Array<Element> & elements,
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

  /// energy split
  std::shared_ptr<EnergySplit_<dim>> energy_split;
};

} // namespace akantu
/* -------------------------------------------------------------------------- */
namespace akantu {
namespace {
template <template <Int, template <Int> class> class Mat>
bool instantiateMaterialPhaseField(const ID & id) {
  return MaterialFactory::getInstance().registerAllocator(
      id,
      [](Int dim, const ID & energy_split, SolidMechanicsModel & model,
         const ID & id) -> std::unique_ptr<Material> {
        return tuple_dispatch<AllSpatialDimensions>(
            [&](auto && _) -> std::unique_ptr<Material> {
              constexpr auto && dim_ = aka::decay_v<decltype(_)>;

              if (energy_split.empty() or energy_split == "no_split") {
                return std::make_unique<Mat<dim_, NoEnergySplit>>(model, id);
              }
              if (energy_split == "volumetric_deviatoric") {
                return std::make_unique<Mat<dim_, VolumetricDeviatoricSplit>>(
                    model, id);
              }
              AKANTU_ERROR("Unknown energy split type: " << energy_split);
              return nullptr;
            },
            dim);
      });
}
} // namespace
} // namespace akantu

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "material_phasefield_inline_impl.hh"

#endif /* __AKANTU_MATERIAL_PHASEFIELD_HH__ */
