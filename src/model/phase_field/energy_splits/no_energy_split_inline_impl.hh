
/* -------------------------------------------------------------------------- */
#include "no_energy_split.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int dim>
inline void
NoEnergySplit<dim>::computePhiOnQuad(const Matrix<Real> & strain_quad,
                                     Real & phi_quad) {
  Real trace = strain_quad.trace();
  phi_quad = 0.5 * this->lambda * trace * trace +
             this->mu * strain_quad.doubleDot(strain_quad);
}

/* -------------------------------------------------------------------------- */
template <Int dim>
inline void
NoEnergySplit<dim>::computePhiMinusOnQuad(const Matrix<Real> & /*strain_quad*/,
                                          Real & phi_minus_quad) {
  phi_minus_quad = 0.0;
}

/* -------------------------------------------------------------------------- */
template <Int dim>
inline void NoEnergySplit<dim>::computeSigmaOnQuad(
    const Matrix<Real> & strain_quad, const Real & sigma_th,
    Matrix<Real> & sigma_plus, Matrix<Real> & /*sigma_minus*/) {
  auto strain = strain_quad - sigma_th * Matrix<Real>::Identity(dim, dim);
  Real trace = strain.trace();
  sigma_plus = 2. * this->mu * strain +
               this->lambda * trace * Matrix<Real>::Identity(dim, dim);
}

/* -------------------------------------------------------------------------- */
template <Int dim>
inline void NoEnergySplit<dim>::computeTangentCoefsOnQuad(
    const Matrix<Real> & /*strain_quad*/, const Real & g_d,
    Matrix<Real> & tangent) {

  constexpr auto n = (dim * (dim - 1) / 2 + dim);

  if constexpr (dim == 1) {
    tangent(0, 0) = this->E * g_d;
    return;
  }

  auto Miiii = g_d * (this->lambda + 2. * this->mu);
  [[maybe_unused]] auto Miijj = g_d * this->lambda;
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
template class NoEnergySplit<1>;
template class NoEnergySplit<2>;
template class NoEnergySplit<3>;

} // namespace akantu
