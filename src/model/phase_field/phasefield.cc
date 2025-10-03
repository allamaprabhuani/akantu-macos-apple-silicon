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
#include "phasefield.hh"
#include "aka_common.hh"
#include "no_energy_split.hh"
#include "phase_field_model.hh"
#include "random_internal_field.hh"
#include "volumetric_deviatoric_split.hh"
#include <string>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
PhaseField::PhaseField(PhaseFieldModel & model, const ID & id,
                       const ID & fe_engine_id)
    : Parent(model, id, model.getSpatialDimension(), _ek_regular, fe_engine_id),
      g_c(this->registerInternal<Real, DefaultRandomInternalField>(
          "g_c", 1, fe_engine_id)),
      damage_on_qpoints(this->registerInternal("damage", 1, fe_engine_id)),
      gradd(this->registerInternal("grad_d", spatial_dimension, fe_engine_id)),
      phi(this->registerInternal("phi", 1, fe_engine_id)),
      phi_minus(this->registerInternal("phi_minus", 1, fe_engine_id)),
      phi_plus(this->registerInternal("phi_plus", 1, fe_engine_id)),
      strain(this->registerInternal(
          "strain", spatial_dimension * spatial_dimension, fe_engine_id)),
      driving_force(this->registerInternal("driving_force", 1, fe_engine_id)),
      driving_energy(this->registerInternal("driving_energy", spatial_dimension,
                                            fe_engine_id)),
      damage_energy(this->registerInternal(
          "damage_energy", spatial_dimension * spatial_dimension,
          fe_engine_id)),
      damage_energy_density(
          this->registerInternal("damage_energy_density", 1, fe_engine_id)),
      dissipated_energy(
          this->registerInternal("dissipated_energy", 1, fe_engine_id)) {

  this->phi.setDefaultValue(0.);
  this->phi_minus.setDefaultValue(0.);
  this->phi_plus.setDefaultValue(0.);
  this->damage_on_qpoints.setDefaultValue(0.);
  this->phi.initializeHistory();
  this->damage_on_qpoints.initializeHistory();

  this->registerParam("l0", l0, Real(0.), _pat_parsable | _pat_readable,
                      "length scale parameter");
  this->registerParam("gc", g_c, _pat_parsable | _pat_readable,
                      "critical local fracture energy density");
  this->registerParam("E", E, _pat_parsable | _pat_readable, "Young's modulus");
  this->registerParam("nu", nu, _pat_parsable | _pat_readable, "Poisson ratio");
  this->registerParam("irreversibility_type", irreversibility_type,
                      std::string("history"), _pat_parsable | _pat_readable,
                      "Type of irreversibility");
  this->registerParam("Plane_Stress", plane_stress, false,
                      _pat_parsable | _pat_readable, "Is plane stress");
}

/* -------------------------------------------------------------------------- */
void PhaseField::updateInternalParameters() {
  lambda = nu * E / ((1 + nu) * (1 - 2 * nu));
  if (plane_stress) {
    lambda = nu * E / ((1 + nu) * (1 - nu));
  }
  mu = E / (2 * (1 + nu));

  if (irreversibility_type == "history") {
    use_history = true;
  } else if (irreversibility_type == "penalization") {
    use_penalization = true;
    use_history = false;
  } else if (irreversibility_type == "tao") {
    use_tao = true;
    use_history = false;
  } else {
    AKANTU_EXCEPTION("Unknown irreversibility type: " << irreversibility_type);
  }

  Parent::updateInternalParameters();
}

/* -------------------------------------------------------------------------- */
void PhaseField::computeAllDrivingForces(GhostType ghost_type) {
  auto & damage = handler.getDamage();
  auto & fem = this->getFEEngine();

  for (const auto & type : this->getElementFilter().elementTypes(
           this->spatial_dimension, ghost_type)) {
    auto & elem_filter = this->getElementFilter(type, ghost_type);
    if (elem_filter.empty()) {
      continue;
    }

    // compute the damage on quadrature points
    auto & damage_interpolated = damage_on_qpoints(type, ghost_type);
    fem.interpolateOnIntegrationPoints(damage, damage_interpolated, 1, type,
                                       ghost_type);

    auto & gradd_vect = gradd(type, _not_ghost);
    /// compute @f$\nabla u@f$
    fem.gradientOnIntegrationPoints(damage, gradd_vect, 1, type, ghost_type,
                                    elem_filter);

    computeDrivingForce(type, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void PhaseField::assembleInternalForces(GhostType ghost_type) {
  Array<Real> & internal_force = handler.getInternalForce();
  auto & fem = this->getFEEngine();

  for (auto type : getElementFilter().elementTypes(_ghost_type = ghost_type)) {
    auto & elem_filter = getElementFilter(type, ghost_type);
    if (elem_filter.empty()) {
      continue;
    }

    auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);

    auto & driving_force_vect = driving_force(type, ghost_type);

    Array<Real> nt_driving_force(0, nb_nodes_per_element);
    fem.computeNtb(driving_force_vect, nt_driving_force, type, ghost_type,
                   elem_filter);

    Array<Real> int_nt_driving_force(0, nb_nodes_per_element);
    fem.integrate(nt_driving_force, int_nt_driving_force, nb_nodes_per_element,
                  type, ghost_type, elem_filter);

    handler.getDOFManager().assembleElementalArrayLocalArray(
        int_nt_driving_force, internal_force, type, ghost_type, -1,
        elem_filter);

    // damage_energy_on_qpoints = gc*l0 = scalar
    auto & driving_energy_vect = driving_energy(type, ghost_type);

    Array<Real> bt_driving_energy(0, nb_nodes_per_element);
    fem.computeBtD(driving_energy_vect, bt_driving_energy, type, ghost_type,
                   elem_filter);

    Array<Real> int_bt_driving_energy(0, nb_nodes_per_element);
    fem.integrate(bt_driving_energy, int_bt_driving_energy,
                  nb_nodes_per_element, type, ghost_type, elem_filter);

    handler.getDOFManager().assembleElementalArrayLocalArray(
        int_bt_driving_energy, internal_force, type, ghost_type, -1,
        elem_filter);
  }
}

/* -------------------------------------------------------------------------- */
void PhaseField::assembleStiffnessMatrix(GhostType ghost_type) {
  AKANTU_DEBUG_INFO("Assemble the new stiffness matrix");
  auto & fem = this->getFEEngine();

  for (auto type :
       getElementFilter().elementTypes(spatial_dimension, ghost_type)) {
    auto & elem_filter = getElementFilter(type, ghost_type);
    if (elem_filter.empty()) {
      return;
    }

    auto nb_element = elem_filter.size();
    auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
    auto nb_quadrature_points = fem.getNbIntegrationPoints(type, ghost_type);

    auto nt_b_n = std::make_unique<Array<Real>>(
        nb_element * nb_quadrature_points,
        nb_nodes_per_element * nb_nodes_per_element, "N^t*b*N");

    auto bt_d_b = std::make_unique<Array<Real>>(
        nb_element * nb_quadrature_points,
        nb_nodes_per_element * nb_nodes_per_element, "B^t*D*B");

    // damage_energy_density_on_qpoints = gc/l0 + phi = scalar
    auto & damage_energy_density_vect = damage_energy_density(type, ghost_type);

    // damage_energy_on_qpoints = gc*l0 = scalar
    auto & damage_energy_vect = damage_energy(type, ghost_type);

    fem.computeBtDB(damage_energy_vect, *bt_d_b, 2, type, ghost_type,
                    elem_filter);

    fem.computeNtbN(damage_energy_density_vect, *nt_b_n, type, ghost_type,
                    elem_filter);

    /// compute @f$ K_{\grad d} = \int_e \mathbf{N}^t * \mathbf{w} *
    /// \mathbf{N}@f$
    auto K_n = std::make_unique<Array<Real>>(
        nb_element, nb_nodes_per_element * nb_nodes_per_element, "K_n");

    fem.integrate(*nt_b_n, *K_n, nb_nodes_per_element * nb_nodes_per_element,
                  type, ghost_type, elem_filter);

    handler.getDOFManager().assembleElementalMatricesToMatrix(
        "K", "damage", *K_n, type, _not_ghost, _symmetric, elem_filter);

    /// compute @f$ K_{\grad d} = \int_e \mathbf{B}^t * \mathbf{W} *
    /// \mathbf{B}@f$
    auto K_b = std::make_unique<Array<Real>>(
        nb_element, nb_nodes_per_element * nb_nodes_per_element, "K_b");

    fem.integrate(*bt_d_b, *K_b, nb_nodes_per_element * nb_nodes_per_element,
                  type, ghost_type, elem_filter);

    handler.getDOFManager().assembleElementalMatricesToMatrix(
        "K", "damage", *K_b, type, _not_ghost, _symmetric, elem_filter);
  }
}

/* -------------------------------------------------------------------------- */
void PhaseField::computeDissipatedEnergyByElements() {
  const Array<Real> & damage = handler.getDamage();
  auto & fem = this->getFEEngine();

  for (auto type :
       getElementFilter().elementTypes(spatial_dimension, _not_ghost)) {
    Array<Idx> & elem_filter = getElementFilter(type, _not_ghost);
    if (elem_filter.empty()) {
      continue;
    }

    Array<Real> & damage_interpolated = damage_on_qpoints(type, _not_ghost);

    // compute the damage on quadrature points
    fem.interpolateOnIntegrationPoints(damage, damage_interpolated, 1, type,
                                       _not_ghost);

    Array<Real> & gradd_vect = gradd(type, _not_ghost);

    /// compute @f$\nabla u@f$
    fem.gradientOnIntegrationPoints(damage, gradd_vect, 1, type, _not_ghost,
                                    elem_filter);

    computeDissipatedEnergy(type);
  }
}

/* -------------------------------------------------------------------------- */
void PhaseField::computeDissipatedEnergy(ElementType /*unused*/) {
  AKANTU_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
PhaseFieldFactory & PhaseField::getFactory() {
  return PhaseFieldFactory::getInstance();
}

/* -------------------------------------------------------------------------- */
Real PhaseField::getEnergy(const ID & energy_id) {
  if (energy_id != "dissipated") {
    return 0.;
  }

  Real edis = 0.;
  auto & fem = this->getFEEngine();

  computeDissipatedEnergyByElements();

  /// integrate the dissipated energy for each type of elements
  for (auto type :
       getElementFilter().elementTypes(spatial_dimension, _not_ghost)) {
    edis += fem.integrate(dissipated_energy(type, _not_ghost), type, _not_ghost,
                          getElementFilter(type, _not_ghost));
  }

  return edis;
}

/* -------------------------------------------------------------------------- */
Real PhaseField::getEnergy(const ID & energy_id, const Element & element) {
  if (energy_id != "dissipated") {
    return 0.;
  }

  auto & fem = this->getFEEngine();
  Vector<Real> edis_on_quad_points(fem.getNbIntegrationPoints(element.type));
  computeDissipatedEnergyByElement(element.type, element.element,
                                   edis_on_quad_points);
  return fem.integrate(edis_on_quad_points, element);
}

/* -------------------------------------------------------------------------- */
void PhaseField::beforeSolveStep() {
  // this->savePreviousState();
  this->computeAllDrivingForces(_not_ghost);
  this->computeAllDrivingForces(_ghost);
}

/* -------------------------------------------------------------------------- */
void PhaseField::afterSolveStep() {}

} // namespace akantu
