/**
 * @file   material_phasefield_anisotropic_inline_impl.cc
 *
 * @author Shad Durussel <shad.durussel@epfl.ch>
 *
 * @date creation: Mon Mar 27 2023
 * @date last modification: Mon Mar 27 2023
 *
 * @brief  Implementation of the inline functions of the material phasefield
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

#include "material_phasefield_anisotropic.hh"
#include <algorithm>

#ifndef AKANTU_MATERIAL_PHASEFIELD_ANISOTROPIC_INLINE_IMPL_HH_
#define AKANTU_MATERIAL_PHASEFIELD_ANISOTROPIC_INLINE_IMPL_HH_
/* -------------------------------------------------------------------------- */
namespace akantu {
template <Int dim>
template <class Args>
inline void
MaterialPhaseFieldAnisotropic<dim>::computeStressOnQuad(Args && args) {

  auto && dam = args["damage"_n];

  auto sigma = args["sigma"_n];
  auto strain = Material::gradUToEpsilon<dim>(args["grad_u"_n]);

  auto && sigma_th = args["sigma_th"_n];

  Real g_d = (1 - dam) * (1 - dam) + eta;

  Matrix<Real> sigma_plus = Matrix<Real>::Zero(dim, dim);
  Matrix<Real> sigma_minus = Matrix<Real>::Zero(dim, dim);
  this->energy_split->computeSigmaOnQuad(strain, sigma_th, sigma_plus,
                                         sigma_minus);

  sigma = g_d * sigma_plus.topLeftCorner(dim, dim) +
          sigma_minus.topLeftCorner(dim, dim);
}

/* -------------------------------------------------------------------------- */
template <Int dim>
template <class Args>
void MaterialPhaseFieldAnisotropic<dim>::computeTangentModuliOnQuad(
    Args && args) {

  auto && dam = args["damage"_n];

  auto strain = Material::gradUToEpsilon<dim>(args["grad_u"_n]);

  Real g_d = (1 - dam) * (1 - dam) + eta;

  auto && tangent = args["tangent_moduli"_n];
  tangent.zero();

  constexpr auto n = Material::getTangentStiffnessVoigtSize(dim);
  Matrix<Real> tmp = Matrix<Real>::Zero(n, n);

  this->energy_split->computeTangentCoefsOnQuad(strain, g_d, tmp);
  tangent = tmp;
}

/* -------------------------------------------------------------------------- */
template <Int dim>
inline Vector<Real>
MaterialPhaseFieldAnisotropic<dim>::getRho(const Element & element) const {
  Vector<Real> rhos = Parent::getRho(element);

  if (not degrade_mass) {
    return rhos;
  }

  auto damage_it = this->damage(element.type, element.ghost_type).begin();

  auto gradu_view = make_view<dim, dim>(this->gradu(element.type));

  auto & fem = this->getFEEngine();
  UInt nb_quadrature_points =
      fem.getNbIntegrationPoints(element.type, element.ghost_type);

  auto gradu_it = gradu_view.begin() + element.element * nb_quadrature_points;
  auto gradu_end = gradu_it + nb_quadrature_points;

  damage_it += element.element * nb_quadrature_points;

  Real rho_base = Parent::getRho();
  for (auto & rho : rhos) {
    Real trace = gradu_it->trace();
    if (trace > 0) {
      rho *= (1 - *damage_it) * (1 - *damage_it) + eta;
      rho *= (1 - *damage_it) * (1 - *damage_it) + eta;
      rho = std::min(rho_base, rho);
    }
    ++damage_it;
    ++gradu_it;
  }
  return rhos;
}
/* -------------------------------------------------------------------------- */
template <Int dim>
inline Int
MaterialPhaseFieldAnisotropic<dim>::getNbData(const Array<Element> & elements,
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
MaterialPhaseFieldAnisotropic<dim>::packData(CommunicationBuffer & buffer,
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
MaterialPhaseFieldAnisotropic<dim>::unpackData(CommunicationBuffer & buffer,
                                    const Array<Element> & elements,
                                    const SynchronizationTag & tag) {
  Parent::unpackData(buffer, elements, tag);

  if (tag == SynchronizationTag::_smm_density && degrade_mass) {
    this->unpackInternalFieldHelper(this->damage, buffer, elements);
  }
}


} // namespace akantu
#endif
