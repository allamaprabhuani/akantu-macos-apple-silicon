
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
  VolumetricDeviatoricSplit() : EnergySplit() {
    // dev_dim = plane_stress ? dim : 3;
    dev_dim = dim;
  };

  /* ------------------------------------------------------------------------*/
  /* Methods */
  /* ------------------------------------------------------------------------*/
public:
  // compute strain energy density on quad
  inline void computePhiOnQuad(const Matrix<Real> & strain_quad,
                               Real & phi_quad);

  // compute stress on quad
  inline void computeSigmaOnQuad(const Matrix<Real> & strain_quad,
                                 const Real & epsilon_th,
                                 Matrix<Real> & sigma_plus,
                                 Matrix<Real> & sigma_minus);

  // compute tangent moduli coefficients on quad
  inline void computeTangentCoefsOnQuad(const Matrix<Real> & strain_quad,
                                        const Real & g_d,
                                        Matrix<Real> & tangent);

  /* ------------------------------------------------------------------------
   */
  /* Class Members */
  /* ------------------------------------------------------------------------
   */
private:
  // dimension to consider in deviatoric split
  Int dev_dim;
};

} // namespace akantu

#include "volumetric_deviatoric_split_inline_impl.hh"

#endif
