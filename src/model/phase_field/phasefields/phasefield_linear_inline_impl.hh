#include "phasefield_linear.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int dim>
inline void PhaseFieldLinear<dim>::computeDissipatedEnergyOnQuad(
    const Real & dam, const Real & /*dam_prev*/, const Vector<Real> & grad_d,
    Real & edis, Real & g_c_quad) {

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

/* -------------------------------------------------------------------------- */
template <Int dim>
inline void
PhaseFieldLinear<dim>::computePhiOnQuad(const Matrix<Real> & strain_quad,
                                        Real & phi_quad) {
  Real trace = strain_quad.trace();
  Real trace_plus = std::max(Real(0.), trace);

  Matrix<Real> strain_dev = Matrix<Real>::Zero(dev_dim, dev_dim);
  Matrix<Real> strain_tmp = Matrix<Real>::Zero(dev_dim, dev_dim);
  strain_tmp.topLeftCorner(dim, dim) = strain_quad;

  strain_dev = strain_tmp -
               trace / Real(dev_dim) * Matrix<Real>::Identity(dev_dim, dev_dim);

  Real kpa = this->lambda + 2. * this->mu / Real(dev_dim);

  phi_quad = 0.5 * kpa * trace_plus * trace_plus +
             this->mu * strain_dev.doubleDot(strain_dev);
}

/* -------------------------------------------------------------------------- */
template <Int dim>
inline void PhaseFieldLinear<dim>::computePhiIsotropicOnQuad(
    const Matrix<Real> & strain_quad, Real & phi_quad) {
  Real trace = strain_quad.trace();

  phi_quad = 0.5 * this->lambda * trace * trace +
             this->mu * strain_quad.doubleDot(strain_quad);
}

} // namespace akantu
