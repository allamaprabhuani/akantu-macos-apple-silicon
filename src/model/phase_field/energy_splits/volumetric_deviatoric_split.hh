
/* -------------------------------------------------------------------------- */
#include "energy_split.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_VOLUMETRIC_DEVIATORIC_SPLIT_HH_
#define AKANTU_VOLUMETRIC_DEVIATORIC_SPLIT_HH_

/* -------------------------------------------------------------------------- */
namespace akantu {

template <Int dim> class VolumetricDeviatoricSplit : public EnergySplit {
public:
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------s------------------------------------------
   */
  VolumetricDeviatoricSplit(Real E, Real nu, bool plane_stress = false);

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
                          Matrix<Real> & sigma_minus) override;

  // compute tangent moduli coefficients on quad
  void computeTangentCoefsOnQuad(const Matrix<Real> & strain_quad,
                                 const Real & g_d,
                                 Matrix<Real> & tangent) override;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  // dimension to consider in deviatoric split
  Int dev_dim;
};

} // namespace akantu

#endif
