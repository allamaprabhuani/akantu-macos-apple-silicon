/**
 * @file   phasefield_exponential.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Fri Jun 19 2020
 * @date last modification: Wed Jun 23 2021
 *
 * @brief  Specialization of the phasefield law class for exponential type
 * law
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2018-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "phasefield_quadratic.hh"
#include "aka_common.hh"
#include <cmath>

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int dim, template <Int> class EnergySplit_>
requires CanComputePhi<dim, EnergySplit_>
PhaseFieldQuadratic<dim, EnergySplit_>::PhaseFieldQuadratic(
    PhaseFieldModel & model, const ID & id)
    : PhaseField(model, id) {
  registerParam("irreversibility_tol", tol_ir, Real(1e-2),
                _pat_parsable | _pat_readable, "Irreversibility tolerance");
}

/* -------------------------------------------------------------------------- */
template <Int dim, template <Int> class EnergySplit_>
requires CanComputePhi<dim, EnergySplit_>
void PhaseFieldQuadratic<dim, EnergySplit_>::initPhaseField() {
  PhaseField::initPhaseField();

  this->gamma = Real(this->g_c) / this->l0 * (1. / (tol_ir * tol_ir) - 1.);

  this->energy_split = std::make_shared<EnergySplit_<dim>>();
  this->energy_split->updateMaterialProperties(this->E, this->nu,
                                               this->plane_stress);
}

/* -------------------------------------------------------------------------- */
template <Int dim, template <Int> class EnergySplit_>
requires CanComputePhi<dim, EnergySplit_>
void PhaseFieldQuadratic<dim, EnergySplit_>::updateInternalParameters() {
  PhaseField::updateInternalParameters();

  for (const auto & type : getElementFilter().elementTypes(dim, _not_ghost)) {
    auto && arguments = PhaseField::getArguments<dim>(type, _not_ghost);
    for (auto && args : arguments) {
      auto & dam_energy_quad = args["damage_energy"_n];
      auto & gc_quad = args["g_c"_n];
      auto & phi_hist_quad = args["previous_phi"_n];
      Matrix<Real> d(dim, dim);
      // eye g_c * l0
      d.eye(gc_quad * this->l0);
      dam_energy_quad = d;
      phi_hist_quad = 0;
    }
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim, template <Int> class EnergySplit_>
requires CanComputePhi<dim, EnergySplit_>
void PhaseFieldQuadratic<dim, EnergySplit_>::computeDrivingForce(
    ElementType el_type, GhostType ghost_type) {

  auto && arguments = PhaseField::getArguments<dim>(el_type, ghost_type);

  for (auto && args : arguments) {
    auto & phi_quad = args["phi"_n];
    auto & strain = args["strain"_n];
    this->energy_split->computePhiOnQuad(strain, phi_quad);
  }

  if (this->use_history) {
    for (auto && args : arguments) {
      auto & phi_quad = args["phi"_n];
      auto & phi_hist_quad = args["previous_phi"_n];
      if (phi_quad < phi_hist_quad) {
        phi_quad = phi_hist_quad;
      }
    }
  }

  for (auto && args : arguments) {
    auto & phi_quad = args["phi"_n];
    auto & driving_force_quad = args["driving_force"_n];
    auto & dam_energy_density_quad = args["damage_energy_density"_n];
    auto & driving_energy_quad = args["driving_energy"_n];
    auto & g_c_quad = args["g_c"_n];

    computeDamageEnergyDensityOnQuad(phi_quad, dam_energy_density_quad,
                                     g_c_quad);
    driving_force_quad = -2 * phi_quad;
    driving_energy_quad *= 0;
  }

  computeResidual(el_type, ghost_type);

  if (this->use_penalization) {
    applyPenalization(el_type, ghost_type);
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim, template <Int> class EnergySplit_>
requires CanComputePhi<dim, EnergySplit_>
void PhaseFieldQuadratic<dim, EnergySplit_>::computeResidual(
    ElementType el_type, GhostType ghost_type) {
  auto && arguments = PhaseField::getArguments<dim>(el_type, ghost_type);
  for (auto && args : arguments) {
    auto & dam_energy_density_quad = args["damage_energy_density"_n];
    auto & driving_force_quad = args["driving_force"_n];
    auto & driving_energy_quad = args["driving_energy"_n];
    auto & dam_on_quad = args["damage"_n];
    auto & gradd_quad = args["gradd"_n];
    auto & damage_energy_quad = args["damage_energy"_n];

    driving_force_quad += dam_on_quad * dam_energy_density_quad;
    driving_energy_quad += damage_energy_quad * gradd_quad;
  }
}

template <Int dim, template <Int> class EnergySplit_>
requires CanComputePhi<dim, EnergySplit_>
void PhaseFieldQuadratic<dim, EnergySplit_>::applyPenalization(
    ElementType el_type, GhostType ghost_type) {
  auto && arguments = PhaseField::getArguments<dim>(el_type, ghost_type);
  for (auto && args : arguments) {
    auto & dam_on_quad = args["damage"_n];
    auto & dam_prev_quad = args["previous_damage"_n];
    auto & dam_energy_density_quad = args["damage_energy_density"_n];
    auto & driving_force_quad = args["driving_force"_n];

    Real penalization =
        this->gamma * std::min(Real(0.), dam_on_quad - dam_prev_quad);

    driving_force_quad += penalization;
    dam_energy_density_quad += this->gamma * (dam_on_quad < dam_prev_quad);
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim, template <Int> class EnergySplit_>
requires CanComputePhi<dim, EnergySplit_>
void PhaseFieldQuadratic<dim, EnergySplit_>::computeDissipatedEnergy(
    ElementType el_type) {
  AKANTU_DEBUG_IN();

  auto && arguments = PhaseField::getArguments<dim>(el_type, _not_ghost);

  for (auto && args : arguments) {
    auto & dam_on_quad = args["damage"_n];
    auto & gradd_quad = args["gradd"_n];
    auto & g_c_quad = args["g_c"_n];
    auto & edis_quad = args["dissipated_energy"_n];

    computeDissipatedEnergyOnQuad(dam_on_quad, gradd_quad, edis_quad, g_c_quad);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim, template <Int> class EnergySplit_>
requires CanComputePhi<dim, EnergySplit_>
void PhaseFieldQuadratic<dim, EnergySplit_>::computeDissipatedEnergyByElement(
    ElementType type, Idx index, Vector<Real> & edis_on_quad_points) {
  auto gradd_it = this->gradd(type).begin(dim);
  auto gradd_end = this->gradd(type).begin(dim);
  auto damage_it = this->damage_on_qpoints(type).begin();
  auto g_c_it = this->g_c(type).begin();

  auto & fem = this->getFEEngine();
  UInt nb_quadrature_points = fem.getNbIntegrationPoints(type);

  gradd_it += index * nb_quadrature_points;
  gradd_end += (index + 1) * nb_quadrature_points;
  damage_it += index * nb_quadrature_points;
  g_c_it += index * nb_quadrature_points;

  Real * edis_quad = edis_on_quad_points.data();

  for (; gradd_it != gradd_end; ++gradd_it, ++damage_it, ++edis_quad) {
    this->computeDissipatedEnergyOnQuad(*damage_it, *gradd_it, *edis_quad,
                                        *g_c_it);
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim, template <Int> class EnergySplit_>
requires CanComputePhi<dim, EnergySplit_>
void PhaseFieldQuadratic<dim, EnergySplit_>::computeDissipatedEnergyByElement(
    const Element & element, Vector<Real> & edis_on_quad_points) {
  computeDissipatedEnergyByElement(element.type, element.element,
                                   edis_on_quad_points);
}

/* -------------------------------------------------------------------------- */
const bool phase_field_quadratic_is_allocated [[maybe_unused]] =
    instantiatePhaseField<PhaseFieldQuadratic>("quadratic");

} // namespace akantu
