/**
 * @file   material_phasefield_anisotropic.cc
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
#include "material_phasefield_anisotropic.hh"
#include "aka_common.hh"
#include "solid_mechanics_model.hh"
#include "volumetric_deviatoric_split.hh"
#include <memory>

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int dim>
MaterialPhaseFieldAnisotropic<dim>::MaterialPhaseFieldAnisotropic(
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

template <Int dim>
void MaterialPhaseFieldAnisotropic<dim>::computeStress(ElementType el_type,
                                                       GhostType ghost_type) {

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
template <Int dim> void MaterialPhaseFieldAnisotropic<dim>::initMaterial() {
  MaterialDamage<dim>::initMaterial();
  this->energy_split =
      std::make_shared<VolumetricDeviatoricSplit<dim>>(this->E, this->nu);
}

/* -------------------------------------------------------------------------- */
template <> void MaterialPhaseFieldAnisotropic<2>::initMaterial() {
  MaterialDamage<2>::initMaterial();

  this->dev_dim = 2;
  if (!this->plane_stress) {
    this->dev_dim = 3;
  }
  this->energy_split = std::make_shared<VolumetricDeviatoricSplit<2>>(
      this->E, this->nu, this->plane_stress);
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialPhaseFieldAnisotropic<dim>::computeTangentModuli(
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
template class MaterialPhaseFieldAnisotropic<1>;
template class MaterialPhaseFieldAnisotropic<2>;
template class MaterialPhaseFieldAnisotropic<3>;

static bool material_is_allocated_phasefield =
    instantiateMaterial<MaterialPhaseFieldAnisotropic>(
        "phasefield_anisotropic");

} // namespace akantu
