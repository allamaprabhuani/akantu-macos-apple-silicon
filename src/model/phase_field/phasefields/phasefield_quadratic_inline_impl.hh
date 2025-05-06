#include "phasefield_quadratic.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int dim>
inline void PhaseFieldQuadratic<dim>::computeDissipatedEnergyOnQuad(
    const Real & dam, const Vector<Real> & grad_d, Real & edis,
    Real & g_c_quad) {

  edis = 0.;
  for (auto i : arange(dim)) {
    edis += 0.5 * g_c_quad * this->l0 * grad_d[i] * grad_d[i];
  }

  edis += g_c_quad * dam * dam / (2. * this->l0);
}

/* -------------------------------------------------------------------------- */
template <Int dim>
inline void PhaseFieldQuadratic<dim>::computeDamageEnergyDensityOnQuad(
    const Real & phi_quad, Real & dam_energy_quad, const Real & g_c_quad) {
  dam_energy_quad = 2.0 * phi_quad + g_c_quad / this->l0;
}

} // namespace akantu
