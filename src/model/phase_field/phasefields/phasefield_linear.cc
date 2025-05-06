/* -------------------------------------------------------------------------- */
#include "phasefield_linear.hh"
#include "aka_common.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int dim>
PhaseFieldLinear<dim>::PhaseFieldLinear(PhaseFieldModel & model, const ID & id)
    : PhaseField(model, id) {
  registerParam("irreversibility_tol", tol_ir, Real(1e-2),
                _pat_parsable | _pat_readable, "Irreversibility tolerance");
  registerParam("recovery_tol", tol_rec, Real(1e-2),
                _pat_parsable | _pat_readable, "Recovery tolerance");
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
    auto && arguments = PhaseField::getArguments<dim>(type, _not_ghost);
    for (auto && args : arguments) {
      auto & dam_energy_quad = args["damage_energy"_n];
      auto & gc_quad = args["g_c"_n];
      auto & phi_hist_quad = args["previous_phi"_n];
      Matrix<Real> d(dim, dim);
      // eye g_c * l0
      d.eye(3. / 4. * gc_quad * this->l0);
      dam_energy_quad = d;
      phi_hist_quad = 3 * gc_quad / (16 * this->l0);
    }
  }
}

/* -------------------------------------------------------------------------- */
// template<Int dim, class EnergySplit>
template <Int dim>
void PhaseFieldLinear<dim>::computeDrivingForce(ElementType el_type,
                                                GhostType ghost_type) {
  auto && arguments = PhaseField::getArguments<dim>(el_type, ghost_type);

  for (auto && args : arguments) {
    auto & phi_quad = args["phi"_n];
    auto & strain = args["strain"_n];
    // EnergySplit::computePhiOnQuad(strain, phi_quad);
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

    computeDamageEnergyDensityOnQuad(phi_quad, dam_energy_density_quad);
    driving_force_quad = -2 * phi_quad + 3 * g_c_quad / (8 * this->l0);
    driving_energy_quad *= 0;
  }

  if (not this->use_tao) {
    computeResidual(el_type, ghost_type);
  }

  if (this->use_penalization) {
    applyPenalization(el_type, ghost_type);
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void PhaseFieldLinear<dim>::computeResidual(ElementType el_type,
                                            GhostType ghost_type) {
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

template <Int dim>
void PhaseFieldLinear<dim>::applyPenalization(ElementType el_type,
                                              GhostType ghost_type) {
  auto && arguments = PhaseField::getArguments<dim>(el_type, ghost_type);
  for (auto && args : arguments) {
    auto & dam_on_quad = args["damage"_n];
    auto & dam_prev_quad = args["previous_damage"_n];
    auto & dam_energy_density_quad = args["damage_energy_density"_n];
    auto & driving_force_quad = args["driving_force"_n];

    Real penalization_ir =
        this->gamma * std::min(Real(0.), dam_on_quad - dam_prev_quad);
    Real penalization_rec = this->rho_rec * std::min(Real(0.), dam_on_quad);

    driving_force_quad += penalization_ir + penalization_rec;
    dam_energy_density_quad += this->gamma * (dam_on_quad < dam_prev_quad);
    dam_energy_density_quad += this->rho_rec * (dam_on_quad < 0);
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void PhaseFieldLinear<dim>::computeDissipatedEnergy(ElementType el_type) {
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

/* --------------------------------------------------------------------------
 */
template <Int dim>
void PhaseFieldLinear<dim>::computeDissipatedEnergyByElement(
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

/* --------------------------------------------------------------------------
 */
template <Int dim>
void PhaseFieldLinear<dim>::computeDissipatedEnergyByElement(
    const Element & element, Vector<Real> & edis_on_quad_points) {
  computeDissipatedEnergyByElement(element.type, element.element,
                                   edis_on_quad_points);
}

/* --------------------------------------------------------------------------
 */
template <Int dim> void PhaseFieldLinear<dim>::afterSolveStep() {
  // clamp negative damage to 0
  // for (auto & dam : this->getHandler().getDamage()) {
  //   dam = std::max(Real(0.), dam);
  // }
}

/* --------------------------------------------------------------------------
 */
template class PhaseFieldLinear<1>;
template class PhaseFieldLinear<2>;
template class PhaseFieldLinear<3>;

const bool phase_field_linear_is_allocated [[maybe_unused]] =
    instantiatePhaseField<PhaseFieldLinear>("linear");

} // namespace akantu
