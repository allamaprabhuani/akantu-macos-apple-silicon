
/* -------------------------------------------------------------------------- */
#include "energy_split.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_NO_ENERGY_SPLIT_HH_
#define AKANTU_NO_ENERGY_SPLIT_HH_

/* -------------------------------------------------------------------------- */
namespace akantu {

template <Int dim> class NoEnergySplit : public EnergySplit {
public:
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
  NoEnergySplit() : EnergySplit() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
  // compute strain energy density on quad
  inline void computePhiOnQuad(const Matrix<Real> & strain_quad,
                               Real & phi_quad);

  // compute stress on quad
  inline void computeSigmaOnQuad(const Matrix<Real> & strain_quad,
                                 const Real & epsilon_th,
                                 Matrix<Real> & sigma_plus,
                                 Matrix<Real> & /*sigma_minus*/);

  // compute tangent moduli coefficients on quad
  inline void computeTangentCoefsOnQuad(const Matrix<Real> & /*strain_quad*/,
                                        const Real & g_d,
                                        Matrix<Real> & tangent);
};

} // namespace akantu

#include "no_energy_split_inline_impl.hh"

#endif
