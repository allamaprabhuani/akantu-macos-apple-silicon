
/* -------------------------------------------------------------------------- */
#include "volumetric_deviatoric_split.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int dim>
inline void VolumetricDeviatoricSplit<dim>::computePhiOnQuad(
    const Matrix<Real> & strain_quad, Real & phi_quad) {
  Real trace = strain_quad.trace();
  Real trace_plus = std::max(Real(0.), trace);

  Matrix<Real> strain_dev =
      strain_quad - trace / 3. * Matrix<Real>::Identity(dim, dim);

  Real kappa = this->lambda + 2. / 3. * this->mu;

  phi_quad = 0.5 * kappa * trace_plus * trace_plus +
             this->mu * strain_dev.doubleDot(strain_dev);
}

/* -------------------------------------------------------------------------- */
template <Int dim>
inline void VolumetricDeviatoricSplit<dim>::computePhiMinusOnQuad(
    const Matrix<Real> & strain_quad, Real & phi_minus_quad) {
  Real trace = strain_quad.trace();
  Real trace_minus = std::min(Real(0.), trace);

  Real kappa = this->lambda + 2. / 3. * this->mu;

  phi_minus_quad = 0.5 * kappa * trace_minus * trace_minus;
}

/* -------------------------------------------------------------------------- */
template <Int dim>
inline void VolumetricDeviatoricSplit<dim>::computeSigmaOnQuad(
    const Matrix<Real> & strain_quad, const Real & epsilon_th,
    Matrix<Real> & sigma_plus, Matrix<Real> & sigma_minus) {
  auto strain = strain_quad - epsilon_th * Matrix<Real>::Identity(dim, dim);
  Real trace = strain.trace();
  Real trace_plus = std::max(Real(0.), trace);
  Real trace_minus = std::min(Real(0.), trace);

  Matrix<Real> strain_dev =
      strain - trace / 3 * Matrix<Real>::Identity(dim, dim);

  Real kappa = this->lambda + 2. / 3. * this->mu;

  sigma_plus = (kappa * trace_plus) * Matrix<Real>::Identity(dim, dim) +
               2. * this->mu * strain_dev;
  sigma_minus = (kappa * trace_minus) * Matrix<Real>::Identity(dim, dim);
}

/* -------------------------------------------------------------------------- */
template <Int dim>
inline void VolumetricDeviatoricSplit<dim>::computeTangentCoefsOnQuad(
    const Matrix<Real> & strain_quad, const Real & g_d,
    Matrix<Real> & tangent) {

  constexpr auto n = (dim * (dim - 1) / 2 + dim);

  Real kappa = this->lambda + 2. / 3. * this->mu;

  Real g_d_hyd = strain_quad.trace() > 0 ? g_d : 1;

  if constexpr (dim == 1) {
    tangent(0, 0) = this->E * g_d_hyd;
    return;
  }

  auto Miiii = g_d_hyd * kappa + g_d * 2. * this->mu * (1. - 1. / 3.);
  [[maybe_unused]] auto Miijj = g_d_hyd * kappa - g_d * 2. * this->mu / 3.;
  [[maybe_unused]] auto Mijij = g_d * this->mu;

  tangent(0, 0) = Miiii;

  if constexpr (dim >= 2) {
    tangent(1, 1) = Miiii;
    tangent(0, 1) = Miijj;
    tangent(1, 0) = Miijj;

    tangent(n - 1, n - 1) = Mijij;
  }

  if constexpr (dim == 3) {
    tangent(2, 2) = Miiii;
    tangent(0, 2) = Miijj;
    tangent(2, 0) = Miijj;
    tangent(1, 2) = Miijj;
    tangent(2, 1) = Miijj;

    tangent(3, 3) = Mijij;
    tangent(4, 4) = Mijij;
  }
}

/* -------------------------------------------------------------------------- */
template class VolumetricDeviatoricSplit<1>;
template class VolumetricDeviatoricSplit<2>;
template class VolumetricDeviatoricSplit<3>;

} // namespace akantu
