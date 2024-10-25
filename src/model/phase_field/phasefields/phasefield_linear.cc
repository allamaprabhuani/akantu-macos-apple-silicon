/* -------------------------------------------------------------------------- */
#include "phasefield_linear.hh"
#include "aka_common.hh"
#include <tuple>

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int dim>
PhaseFieldLinear<dim>::PhaseFieldLinear(PhaseFieldModel & model, const ID & id)
    : PhaseField(model, id) {
  registerParam("irreversibility_tol", tol_ir, Real(1e-2),
                _pat_parsable | _pat_readable, "Irreversibility tolerance");
  registerParam("recovery_tol", tol_rec, Real(1e-2),
                _pat_parsable | _pat_readable, "Recovery tolerance");
  registerParam("non_linear", non_linear, true, _pat_parsable | _pat_readable,
                "Non linear assembly");
}

/* -------------------------------------------------------------------------- */
template <Int dim> void PhaseFieldLinear<dim>::initPhaseField() {
  PhaseField::initPhaseField();

  // initiate irreversibility param
  this->gamma = Real(this->g_c) / this->l0 * 27. / (64. * tol_ir * tol_ir);

  // initiate recovery param
  Real diam = 0.;
  BBox bbox = this->getHandler().getMesh().getBBox();
  for (auto i : arange(dim)) {
    diam += bbox.size(SpatialDirection(i)) * bbox.size(SpatialDirection(i));
  }
  diam = std::sqrt(diam);
  this->rho_rec = Real(this->g_c) / this->l0 * 9. * (diam / this->l0 - 2.) /
                  (64. * tol_rec);

  this->dev_dim = dim;
  if (dim == 2 && !this->plane_stress) {
    this->dev_dim = 3;
  }

  if (this->isotropic) {
    this->energy_split = std::make_shared<NoEnergySplit<dim>>(
        this->E, this->nu, this->plane_stress);
  } else {
    this->energy_split = std::make_shared<VolumetricDeviatoricSplit<dim>>(
        this->E, this->nu, this->plane_stress);
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim> void PhaseFieldLinear<dim>::updateInternalParameters() {
  PhaseField::updateInternalParameters();

  for (const auto & type : getElementFilter().elementTypes(dim, _not_ghost)) {
    for (auto && tuple :
         zip(make_view(this->damage_energy(type, _not_ghost), dim, dim),
             this->g_c(type, _not_ghost))) {
      Matrix<Real> d(dim, dim);
      // eye g_c * l0
      d.eye(3. / 4. * std::get<1>(tuple) * this->l0);
      std::get<0>(tuple) = d;
    }
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void PhaseFieldLinear<dim>::computeDrivingForce(ElementType el_type,
                                                GhostType ghost_type) {

  for (auto && tuple :
       zip(this->phi(el_type, ghost_type),
           make_view(this->strain(el_type, ghost_type), dim, dim))) {
    auto & phi_quad = std::get<0>(tuple);
    auto & strain = std::get<1>(tuple);
    // computePhiOnQuad(strain, phi_quad);
    this->energy_split->computePhiOnQuad(strain, phi_quad);
  }

  for (auto && tuple :
       zip(this->phi(el_type, ghost_type),
           this->driving_force(el_type, ghost_type),
           this->damage_energy_density(el_type, ghost_type),
           this->damage_on_qpoints(el_type, _not_ghost),
           make_view(this->driving_energy(el_type, ghost_type), dim),
           make_view(this->damage_energy(el_type, ghost_type), dim, dim),
           make_view(this->gradd(el_type, ghost_type), dim),
           this->g_c(el_type, ghost_type),
           this->damage_on_qpoints.previous(el_type, ghost_type))) {
    auto & phi_quad = std::get<0>(tuple);
    auto & driving_force_quad = std::get<1>(tuple);
    auto & dam_energy_density_quad = std::get<2>(tuple);
    auto & dam_on_quad = std::get<3>(tuple);
    auto & driving_energy_quad = std::get<4>(tuple);
    auto & damage_energy_quad = std::get<5>(tuple);
    auto & gradd_quad = std::get<6>(tuple);
    auto & g_c_quad = std::get<7>(tuple);
    auto & dam_prev_quad = std::get<8>(tuple);

    if (non_linear) {
      computeDamageEnergyDensityOnQuad(phi_quad, dam_energy_density_quad);
      Real penalization_ir =
          this->gamma * std::min(Real(0.), dam_on_quad - dam_prev_quad);
      Real penalization_rec = this->rho_rec * std::min(Real(0.), dam_on_quad);

      driving_force_quad = dam_on_quad * dam_energy_density_quad -
                           2 * phi_quad + 3 * g_c_quad / (8 * this->l0) +
                           penalization_ir + penalization_rec;
      driving_energy_quad = damage_energy_quad * gradd_quad;

      dam_energy_density_quad += this->gamma * (dam_on_quad < dam_prev_quad);
      dam_energy_density_quad += this->rho_rec * (dam_on_quad < 0);
    } else {
      computeDamageEnergyDensityOnQuad(phi_quad, dam_energy_density_quad);
      driving_force_quad = -2 * phi_quad + 3 * g_c_quad / (8 * this->l0);
      driving_energy_quad = 0 * gradd_quad;
    }
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void PhaseFieldLinear<dim>::computeDissipatedEnergy(ElementType el_type) {
  AKANTU_DEBUG_IN();

  for (auto && tuple :
       zip(this->dissipated_energy(el_type, _not_ghost),
           this->damage_on_qpoints(el_type, _not_ghost),
           this->damage_on_qpoints.previous(el_type, _not_ghost),
           make_view(this->gradd(el_type, _not_ghost), dim),
           this->g_c(el_type, _not_ghost))) {

    this->computeDissipatedEnergyOnQuad(std::get<1>(tuple), std::get<2>(tuple),
                                        std::get<3>(tuple), std::get<0>(tuple),
                                        std::get<4>(tuple));
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void PhaseFieldLinear<dim>::computeDissipatedEnergyByElement(
    ElementType type, Idx index, Vector<Real> & edis_on_quad_points) {
  auto gradd_it = this->gradd(type).begin(dim);
  auto gradd_end = this->gradd(type).begin(dim);
  auto damage_it = this->damage_on_qpoints(type).begin();
  auto damage_prev_it = this->damage_on_qpoints.previous(type).begin();
  auto g_c_it = this->g_c(type).begin();

  auto & fem = this->getFEEngine();
  UInt nb_quadrature_points = fem.getNbIntegrationPoints(type);

  gradd_it += index * nb_quadrature_points;
  gradd_end += (index + 1) * nb_quadrature_points;
  damage_it += index * nb_quadrature_points;
  damage_prev_it += index * nb_quadrature_points;
  g_c_it += index * nb_quadrature_points;

  Real * edis_quad = edis_on_quad_points.data();

  for (; gradd_it != gradd_end; ++gradd_it, ++damage_it, ++edis_quad) {
    this->computeDissipatedEnergyOnQuad(*damage_it, *damage_prev_it, *gradd_it,
                                        *edis_quad, *g_c_it);
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void PhaseFieldLinear<dim>::computeDissipatedEnergyByElement(
    const Element & element, Vector<Real> & edis_on_quad_points) {
  computeDissipatedEnergyByElement(element.type, element.element,
                                   edis_on_quad_points);
}

/* -------------------------------------------------------------------------- */
template <Int dim> void PhaseFieldLinear<dim>::afterSolveStep() {
  // clamp negative damage to 0
  // for (auto & dam : this->getHandler().getDamage()) {
  //   dam = std::max(Real(0.), dam);
  // }
}

/* -------------------------------------------------------------------------- */
template class PhaseFieldLinear<1>;
template class PhaseFieldLinear<2>;
template class PhaseFieldLinear<3>;

const bool phase_field_linear_is_allocated [[maybe_unused]] =
    instantiatePhaseField<PhaseFieldLinear>("linear");

} // namespace akantu
