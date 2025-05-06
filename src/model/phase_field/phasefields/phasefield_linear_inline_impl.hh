#include "phasefield_linear.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int dim>
inline void PhaseFieldLinear<dim>::computeDissipatedEnergyOnQuad(
    const Real & dam, const Vector<Real> & grad_d, Real & edis,
    Real & g_c_quad) {

  edis = 0.;
  for (auto i : arange(dim)) {
    edis += 3. / 8 * g_c_quad * this->l0 * grad_d[i] * grad_d[i];
  }

  edis += 3 * g_c_quad * dam / (8 * this->l0);
}

/* -------------------------------------------------------------------------- */
template <Int dim>
inline void PhaseFieldLinear<dim>::computeDrivingForceOnQuad(
    const Real & phi_quad, Real & driving_force_quad, const Real & g_c_quad) {
  driving_force_quad = phi_quad - 3 * g_c_quad / (16 * this->l0);
}

/* -------------------------------------------------------------------------- */
template <Int dim>
inline void PhaseFieldLinear<dim>::computeDamageEnergyDensityOnQuad(
    const Real & phi_quad, Real & dam_energy_quad) {
  dam_energy_quad = 2 * phi_quad;
}

} // namespace akantu
