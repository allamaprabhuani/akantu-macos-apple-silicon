/**
 * Copyright (©) 2020-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_factory.hh"
#include "constitutive_law.hh"
#include "energy_split.hh"
#include "no_energy_split.hh"
#include "volumetric_deviatoric_split.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_PHASEFIELD_HH_
#define AKANTU_PHASEFIELD_HH_

/* -------------------------------------------------------------------------- */
namespace akantu {
class PhaseFieldModel;
class PhaseField;
} // namespace akantu

namespace akantu {

using PhaseFieldFactory =
    Factory<PhaseField, ID, Int, const ID &, PhaseFieldModel &, const ID &>;

class PhaseField : public ConstitutiveLaw<PhaseFieldModel> {
  using Parent = ConstitutiveLaw<PhaseFieldModel>;
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  PhaseField(PhaseFieldModel & model, const ID & id = "",
             const ID & fe_engine_id = "");

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// initialize the phasefield computed parameter
  virtual void initPhaseField() { Parent::initConstitutiveLaw(); }
  void initConstitutiveLaw() override { this->initPhaseField(); }

  ///
  virtual void beforeSolveStep();

  virtual void afterSolveStep();

  /// assemble the residual for this phasefield
  virtual void assembleInternalForces(GhostType ghost_type);

  /// assemble the stiffness matrix for this phasefield
  virtual void assembleStiffnessMatrix(GhostType ghost_type);

  /// compute the driving force for this phasefield
  virtual void computeAllDrivingForces(GhostType ghost_type = _not_ghost);

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
  /// compute the dissipated energy by element
  void computeDissipatedEnergyByElements();

  /// function called to updatet the internal parameters when the
  /// modifiable parameters are modified
  void updateInternalParameters() override;

  // constitutive law for driving force
  virtual void computeDrivingForce(ElementType /* el_type */,
                                   GhostType /* ghost_type */ = _not_ghost) {
    AKANTU_TO_IMPLEMENT();
  }

  /// compute the dissiapted energy
  virtual void computeDissipatedEnergy(ElementType el_type);

  /// compute the dissipated energy for an element
  virtual void
  computeDissipatedEnergyByElement(const Element & /*element*/,
                                   Vector<Real> & /*edis_on_quad_points*/) {
    AKANTU_TO_IMPLEMENT();
  }

  /// compute the dissipated energy for an element
  virtual void
  computeDissipatedEnergyByElement(ElementType /*type*/, Idx /*index*/,
                                   Vector<Real> & /*edis_on_quad_points*/) {
    AKANTU_TO_IMPLEMENT();
  }

  template <Int dim>
  decltype(auto) getArguments(ElementType el_type, GhostType ghost_type) {
    using namespace tuple;
    return zip(
        "phi"_n = phi(el_type, ghost_type),
        "previous_phi"_n = phi.previous(el_type, ghost_type),
        "strain"_n = make_view<dim, dim>(strain(el_type, ghost_type)),
        "driving_force"_n = driving_force(el_type, ghost_type),
        "damage_energy_density"_n = damage_energy_density(el_type, ghost_type),
        "damage"_n = damage_on_qpoints(el_type, ghost_type),
        "previous_damage"_n = damage_on_qpoints.previous(el_type, ghost_type),
        "g_c"_n = g_c(el_type, ghost_type),
        "gradd"_n = make_view<dim>(gradd(el_type, ghost_type)),
        "driving_energy"_n =
            make_view<dim>(driving_energy(el_type, ghost_type)),
        "damage_energy"_n =
            make_view<dim, dim>(damage_energy(el_type, ghost_type)),
        "dissipated_energy"_n = dissipated_energy(el_type, ghost_type));
  }

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// static method to reteive the material factory
  static PhaseFieldFactory & getFactory();

  /// return the damage energyfor the subset of elements contained
  /// by the phasefield
  [[nodiscard]] Real getEnergy(const ID & energy_id) override;

  /// Compute dissipated energy for an individual element
  [[nodiscard]] Real getEnergy(const ID & energy_id,
                               const Element & element) override;

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Strain, strain, Real);
  AKANTU_GET_MACRO_AUTO(Strain, strain);
  AKANTU_GET_MACRO_AUTO_NOT_CONST(Strain, strain);

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Damage, damage_on_qpoints, Real);
  AKANTU_GET_MACRO_AUTO_NOT_CONST(Damage, damage_on_qpoints);
  AKANTU_GET_MACRO_AUTO(Damage, damage_on_qpoints);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// length scale parameter
  Real l0{0.};

  /// Young's modulus
  Real E{0.};

  /// Poisson ratio
  Real nu{0.};

  /// Finite deformation
  bool plane_stress{false};

  /// Finite deformation
  bool finite_deformation{false};

  /// Use history
  bool use_history{true};

  /// Use penalization
  bool use_penalization{false};

  /// Use TAO
  bool use_tao{false};

  /// Way of enforcing irreversibility
  std::string irreversibility_type{"history"};

  /// Lame's first parameter
  Real lambda{0.};

  /// Lame's second paramter
  Real mu{0.};

  /// critical energy release rate
  // Real g_c;
  DefaultRandomInternalField<Real> & g_c;

  /// damage arrays ordered by element types
  InternalField<Real> & damage_on_qpoints;

  /// grad_d arrays ordered by element types
  InternalField<Real> & gradd;

  /// phi arrays ordered by element types
  InternalField<Real> & phi;

  /// strain arrays ordered by element types
  InternalField<Real> & strain;

  /// driving force ordered by element types
  InternalField<Real> & driving_force;

  /// driving energy ordered by element types
  InternalField<Real> & driving_energy;

  /// damage energy ordered by element types
  InternalField<Real> & damage_energy;

  /// damage energy density ordered by element types
  InternalField<Real> & damage_energy_density;

  /// dissipated energy by element
  InternalField<Real> & dissipated_energy;
};

} // namespace akantu

#include "phasefield_inline_impl.hh"

#include "internal_field_tmpl.hh"
#include "random_internal_field_tmpl.hh"

namespace akantu {
namespace {
template <template <Int, template <Int> class> class PF>
bool instantiatePhaseField(const ID & id) {
  return PhaseFieldFactory::getInstance().registerAllocator(
      id, [](Int dim, const ID & energy_split, PhaseFieldModel & model,
             const ID & id) {
        return tuple_dispatch<AllSpatialDimensions>(
            [&](auto && _) -> std::unique_ptr<PhaseField> {
              constexpr auto && dim_ = aka::decay_v<decltype(_)>;
              if (energy_split.empty() or energy_split == "no_split") {
                return std::make_unique<PF<dim_, NoEnergySplit>>(model, id);
              } else if (energy_split == "volumetric_deviatoric") {
                return std::make_unique<PF<dim_, VolumetricDeviatoricSplit>>(
                    model, id);
              } else {
                AKANTU_ERROR("Unknown energy split type: " << energy_split);
                return nullptr;
              }
            },
            dim);
      });
}
} // namespace
} // namespace akantu

#endif
