
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
  NoEnergySplit(Real E, Real nu, bool plane_stress = false);

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  // compute strain energy density on quad
  void computePhiOnQuad(const Matrix<Real> & strain_quad,
                        Real & phi_quad) override;

  // compute stress on quad
  void computeSigmaOnQuad(const Matrix<Real> & strain_quad,
                          const Real & sigma_th, Matrix<Real> & sigma_plus,
                          Matrix<Real> & /*sigma_minus*/) override;

  // compute tangent moduli coefficients on quad
  void computeTangentCoefsOnQuad(const Matrix<Real> & /*strain_quad*/,
                                 const Real & g_d,
                                 Matrix<Real> & tangent) override;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
};

} // namespace akantu

#endif
