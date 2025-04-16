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

#include "material_phasefield.hh"
#include "solid_mechanics_model.hh"
#include <algorithm>

#ifndef AKANTU_MATERIAL_PHASEFIELD_INLINE_IMPL_HH_
#define AKANTU_MATERIAL_PHASEFIELD_INLINE_IMPL_HH_
/* -------------------------------------------------------------------------- */
namespace akantu {
template <Int dim>
template <class Args>
inline void MaterialPhaseField<dim>::computeStressOnQuad(Args && args) {

  auto && dam = args["damage"_n];
  args["sigma"_n] *= (1 - dam) * (1 - dam) + eta;
}

/* -------------------------------------------------------------------------- */
template <Int dim>
template <class Args>
void MaterialPhaseField<dim>::computeTangentModuliOnQuad(Args && args) {

  auto dam = args["damage"_n];
  args["tangent_moduli"_n] *= (1 - dam) * (1 - dam) + eta;
}

/* -------------------------------------------------------------------------- */
template <Int dim>
template <class Args>
inline void
MaterialPhaseField<dim>::computeEffectiveDamageOnQuad(Args && args) {
  // using Mat = Matrix<Real, dim, dim>;

  auto strain = Material::gradUToEpsilon<dim>(args["grad_u"_n]);

  Real trace = strain.trace();
  Real trace_plus = std::max(Real(0.), trace);
  Real trace_minus = std::min(Real(0.), trace);

  Matrix<Real> strain_dev(dim, dim);
  strain_dev = strain - trace / Real(dim) * Matrix<Real>::Identity(dim, dim);

  Real kappa = this->lambda + 2. / Real(dim) * this->mu;

  Real strain_energy_plus = 0.5 * kappa * trace_plus * trace_plus +
                            this->mu * strain_dev.doubleDot(strain_dev);
  Real strain_energy_minus = 0.5 * kappa * trace_minus * trace_minus;

  args["effective_damage"_n] =
      args["damage"_n] * (strain_energy_minus < strain_energy_plus);
}

/* -------------------------------------------------------------------------- */
template <Int dim>
inline void MaterialPhaseField<dim>::getRho(Ref<Vector<Real>> rhos,
                                            const Element & element) const {
  Parent::getRho(rhos, element);

  if (not degrade_mass) {
    return;
  }

  auto damage_it = this->damage(element.type, element.ghost_type).begin();

  auto & fem = this->getFEEngine();
  auto nb_quadrature_points =
      fem.getNbIntegrationPoints(element.type, element.ghost_type);

  auto ratio = rhos.size() / nb_quadrature_points;
  damage_it += element.element * nb_quadrature_points;
  auto damage_end = damage_it + nb_quadrature_points;

  auto rho_base = Parent::getRho();
  for (auto && [rho, d] :
       zip(MatrixProxy<Real>(rhos.data(), rhos.size() / ratio,
                             nb_quadrature_points)
               .colwise(),
           range(damage_it, damage_end))) {
    rho *= (1 - d) * (1 - d) + eta;
  }

  for (auto & r : rhos) {
    r = std::min(rho_base, r);
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim>
inline Int
MaterialPhaseField<dim>::getNbData(const Array<Element> & elements,
                                   const SynchronizationTag & tag) const {

  if (tag == SynchronizationTag::_smm_density && degrade_mass) {
    return Int(sizeof(Real)) *
           this->getHandler().getNbIntegrationPoints(elements);
  }

  return Parent::getNbData(elements, tag);
}

/* -------------------------------------------------------------------------- */
template <Int dim>
inline void
MaterialPhaseField<dim>::packData(CommunicationBuffer & buffer,
                                  const Array<Element> & elements,
                                  const SynchronizationTag & tag) const {
  Parent::packData(buffer, elements, tag);

  if (tag == SynchronizationTag::_smm_density && degrade_mass) {
    this->packInternalFieldHelper(this->damage, buffer, elements);
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim>
inline void
MaterialPhaseField<dim>::unpackData(CommunicationBuffer & buffer,
                                    const Array<Element> & elements,
                                    const SynchronizationTag & tag) {
  Parent::unpackData(buffer, elements, tag);

  if (tag == SynchronizationTag::_smm_density && degrade_mass) {
    this->unpackInternalFieldHelper(this->damage, buffer, elements);
  }
}

} // namespace akantu
#endif
