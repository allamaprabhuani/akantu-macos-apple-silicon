/**
 * @file   material_phasefield.cc
 *
 * @author Shad Durussel <shad.durussel@epfl.ch>
 *
 * @date creation: Mon Mar 27 2023
 * @date last modification: Mon Mar 27 2023
 *
 * @brief  Specialization of the material class for the phasefield material
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
#include "material_phasefield.hh"
#include "solid_mechanics_model.hh"
#include "volumetric_deviatoric_split.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int dim, template <Int> class EnergySplit_>
  requires CanComputeSigmaAndTangent<dim, EnergySplit_>
MaterialPhaseField<dim, EnergySplit_>::MaterialPhaseField(
    SolidMechanicsModel & model, const ID & id)
    : Parent(model, id) {
  this->registerParam("eta", eta, Real(0.), _pat_parsable, "eta");
  this->registerParam("is_isotropic", is_isotropic, false,
                      _pat_parsable | _pat_readable,
                      "Use isotropic formulation");
  this->registerParam("degrade_mass", degrade_mass, false,
                      _pat_parsable | _pat_readable,
                      "Degrade mass with damage");
}

/* -------------------------------------------------------------------------- */
template <Int dim, template <Int> class EnergySplit_>
  requires CanComputeSigmaAndTangent<dim, EnergySplit_>
void MaterialPhaseField<dim, EnergySplit_>::computeStress(
    ElementType el_type, GhostType ghost_type) {

  MaterialThermal<dim>::computeStress(el_type, ghost_type);
  auto && arguments = Parent::getArguments(el_type, ghost_type);

  if (not this->finite_deformation) {
    for (auto && args : arguments) {
      this->computeStressOnQuad(args);
    }
  } else {
    for (auto && args : arguments) {
      auto && E = this->template gradUToE<dim>(args["grad_u"_n]);
      this->computeStressOnQuad(tuple::replace(args, "grad_u"_n = E));
    }
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim, template <Int> class EnergySplit_>
  requires CanComputeSigmaAndTangent<dim, EnergySplit_>
void MaterialPhaseField<dim, EnergySplit_>::initMaterial() {
  MaterialDamage<dim>::initMaterial();

  this->energy_split = std::make_shared<EnergySplit_<dim>>();
  if constexpr (dim == 2) {
    this->energy_split->updateMaterialProperties(this->E, this->nu,
                                                 this->plane_stress);
  } else {
    this->energy_split->updateMaterialProperties(this->E, this->nu, false);
  }
}

/* -------------------------------------------------------------------------- */
// template <template <Int> class EnergySplit_>
// void MaterialPhaseField<2, EnergySplit_>::initMaterial() {
//   MaterialDamage<2>::initMaterial();
//
//   this->energy_split = std::make_shared<EnergySplit_<2>>();
//   this->energy_split->updateMaterialProperties(this->E, this->nu,
//                                                this->plane_stress);
// }

/* -------------------------------------------------------------------------- */
template <Int dim, template <Int> class EnergySplit_>
  requires CanComputeSigmaAndTangent<dim, EnergySplit_>
void MaterialPhaseField<dim, EnergySplit_>::computeTangentModuli(
    ElementType el_type, Array<Real> & tangent_matrix, GhostType ghost_type) {

  auto && arguments =
      Parent::getArgumentsTangent(tangent_matrix, el_type, ghost_type);

  if (not this->finite_deformation) {
    for (auto && args : arguments) {
      this->computeTangentModuliOnQuad(args);
    }
  } else {
    for (auto && args : arguments) {
      auto && E = this->template gradUToE<dim>(args["grad_u"_n]);
      this->computeTangentModuliOnQuad(tuple::replace(args, "grad_u"_n = E));
    }
  }

  // for (auto && args :
  //      Parent::getArgumentsTangent(tangent_matrix, el_type, ghost_type)) {
  //   computeTangentModuliOnQuad(args);
  // }
}

/* -------------------------------------------------------------------------- */
template class MaterialPhaseField<1, VolumetricDeviatoricSplit>;
template class MaterialPhaseField<1, NoEnergySplit>;
template class MaterialPhaseField<2, NoEnergySplit>;
template class MaterialPhaseField<3, NoEnergySplit>;
template class MaterialPhaseField<2, VolumetricDeviatoricSplit>;
template class MaterialPhaseField<3, VolumetricDeviatoricSplit>;

static bool material_is_allocated_phasefield [[maybe_unused]] =
    instantiateMaterialPhaseField("phasefield");

} // namespace akantu
