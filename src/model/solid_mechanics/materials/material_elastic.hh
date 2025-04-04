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

#ifndef AKANTU_MATERIAL_ELASTIC_HH_
#define AKANTU_MATERIAL_ELASTIC_HH_

namespace akantu {

/**
 * Material elastic isotropic
 *
 * parameters in the material files :
 *   - E   : Young's modulus (default: 0)
 *   - nu  : Poisson's ratio (default: 1/2)
 *   - Plane_Stress : if 0: plane strain, else: plane stress (default: 0)
 */
template <Int dim> class MaterialElastic : public MaterialThermal<dim> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
private:
  using Parent = MaterialThermal<dim>;

public:
  MaterialElastic(SolidMechanicsModel & model, const ID & id = "",
                  const ID & fe_engine_id = "");

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

  /// compute the elastic potential energy
  void computePotentialEnergy(ElementType el_type) override;

  [[deprecated("Use the interface with an Element")]] void
  computePotentialEnergyByElement(ElementType type, Int index,
                                  Vector<Real> & epot_on_quad_points) override {
    computePotentialEnergyByElement({type, index, _not_ghost},
                                    epot_on_quad_points);
  }

  void
  computePotentialEnergyByElement(const Element & element,
                                  Vector<Real> & epot_on_quad_points) override;

  /// compute the p-wave speed in the material
  [[nodiscard]] auto getPushWaveSpeed(const Element & element) const
      -> Real override;

  /// compute the s-wave speed in the material
  [[nodiscard]] auto getShearWaveSpeed(const Element & element) const
      -> Real override;

protected:
  /// constitutive law for a given quadrature point
  template <typename Args> inline void computeStressOnQuad(Args && args) const;

  /// compute the tangent stiffness matrix for an element
  template <typename Args>
  inline void computeTangentModuliOnQuad(Args && args) const;

  /// recompute the lame coefficient if E or nu changes
  void updateInternalParameters() override;

  template <class Args>
  static inline void computePotentialEnergyOnQuad(Args && args, Real & epot);

  auto hasStiffnessMatrixChanged() -> bool override {
    return (not was_stiffness_assembled);
  }

  auto getTangentType() -> MatrixType override { return _symmetric; }

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get first Lame constant
  AKANTU_GET_MACRO(Lambda, lambda, Real);

  /// get second Lame constant
  AKANTU_GET_MACRO(Mu, mu, Real);

  /// get bulk modulus
  AKANTU_GET_MACRO(Kappa, kpa, Real);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// First Lamé coefficient
  Real lambda{0.};

  /// Second Lamé coefficient (shear modulus)
  Real mu{0.};

  /// Bulk modulus
  Real kpa{0.};

  /// defines if the stiffness was computed
  bool was_stiffness_assembled{false};
};

} // namespace akantu

#include "material_elastic_inline_impl.hh"

#endif /* AKANTU_MATERIAL_ELASTIC_HH_ */
