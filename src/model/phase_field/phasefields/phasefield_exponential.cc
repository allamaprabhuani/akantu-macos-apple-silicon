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
#include "phasefield_exponential.hh"
#include "aka_common.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int dim>
PhaseFieldExponential<dim>::PhaseFieldExponential(PhaseFieldModel & model,
                                                  const ID & id)
    : PhaseField(model, id) {}

/* -------------------------------------------------------------------------- */
template <Int dim> void PhaseFieldExponential<dim>::initPhaseField() {
  PhaseField::initPhaseField();

  this->dev_dim = dim;
  if (dim == 2 && !this->plane_stress) {
    this->dev_dim = 3;
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim> void PhaseFieldExponential<dim>::updateInternalParameters() {
  PhaseField::updateInternalParameters();

  for (const auto & type : getElementFilter().elementTypes(dim, _not_ghost)) {
    for (auto && [dam, gc] :
         zip(make_view<dim, dim>(this->damage_energy(type, _not_ghost)),
             this->g_c(type, _not_ghost))) {
      dam = Matrix<Real, dim, dim>::Identity() * gc * this->l0;
    }
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void PhaseFieldExponential<dim>::computeDrivingForce(ElementType el_type,
                                                     GhostType ghost_type) {

  auto && arguments = zip(
      this->phi(el_type, ghost_type), this->phi.previous(el_type, ghost_type),
      make_view(this->strain(el_type, ghost_type), dim, dim));
  if (this->isotropic) {
    for (auto && [phi_quad, phi_hist_quad, strain] : arguments) {
      computePhiIsotropicOnQuad(strain, phi_quad, phi_hist_quad);
    }
  } else {
    for (auto && [phi_quad, phi_hist_quad, strain] : arguments) {
      computePhiOnQuad(strain, phi_quad, phi_hist_quad);
    }
  }

  for (auto && [phi_quad, driving_force_quad, dam_energy_density_quad,
                dam_on_quad, driving_energy_quad, damage_energy_quad,
                gradd_quad, g_c_quad] :
       zip(this->phi(el_type, ghost_type),
           this->driving_force(el_type, ghost_type),
           this->damage_energy_density(el_type, ghost_type),
           this->damage_on_qpoints(el_type, _not_ghost),
           make_view(this->driving_energy(el_type, ghost_type), dim),
           make_view(this->damage_energy(el_type, ghost_type), dim, dim),
           make_view(this->gradd(el_type, ghost_type), dim),
           this->g_c(el_type, ghost_type))) {
    computeDamageEnergyDensityOnQuad(phi_quad, dam_energy_density_quad,
                                     g_c_quad);

    driving_force_quad = dam_on_quad * dam_energy_density_quad - 2 * phi_quad;
    driving_energy_quad = damage_energy_quad * gradd_quad;
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void PhaseFieldExponential<dim>::computeDissipatedEnergy(ElementType el_type) {
  AKANTU_DEBUG_IN();

  for (auto && [dis_energy, damage, grad_d, g_c] :
       zip(this->dissipated_energy(el_type, _not_ghost),
           this->damage_on_qpoints(el_type, _not_ghost),
           make_view<dim>(this->gradd(el_type, _not_ghost)),
           this->g_c(el_type, _not_ghost))) {
    this->computeDissipatedEnergyOnQuad(damage, grad_d, dis_energy, g_c);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void PhaseFieldExponential<dim>::computeDissipatedEnergyByElement(
    ElementType type, Idx index, Vector<Real> & edis_on_quad_points) {
  auto && fem = this->getFEEngine();
  auto nb_quadrature_points = fem.getNbIntegrationPoints(type);

  auto quad = index * nb_quadrature_points;

  auto gradd_it = make_view<dim>(this->gradd(type)).begin() + quad;
  auto gradd_end = gradd_it + nb_quadrature_points;

  auto damage_it = this->damage_on_qpoints(type).begin() + quad;
  auto g_c_it = this->g_c(type).begin() + quad;

  auto edis_quad = edis_on_quad_points.begin();

  for (; gradd_it != gradd_end; ++gradd_it, ++damage_it, ++edis_quad) {
    this->computeDissipatedEnergyOnQuad(*damage_it, *gradd_it, *edis_quad,
                                        *g_c_it);
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void PhaseFieldExponential<dim>::computeDissipatedEnergyByElement(
    const Element & element, Vector<Real> & edis_on_quad_points) {
  computeDissipatedEnergyByElement(element.type, element.element,
                                   edis_on_quad_points);
}

/* -------------------------------------------------------------------------- */
template class PhaseFieldExponential<1>;
template class PhaseFieldExponential<2>;
template class PhaseFieldExponential<3>;

const bool phase_field_exponential_is_allocated [[maybe_unused]] =
    instantiatePhaseField<PhaseFieldExponential>("exponential");

} // namespace akantu
