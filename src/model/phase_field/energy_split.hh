
/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "aka_types.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_ENERGY_SPLIT_HH_
#define AKANTU_ENERGY_SPLIT_HH_

/* -------------------------------------------------------------------------- */
namespace akantu {

class EnergySplit {
public:
  EnergySplit(const EnergySplit &) = default;
  EnergySplit(EnergySplit &&) = delete;
  auto operator=(const EnergySplit &) -> EnergySplit & = default;
  auto operator=(EnergySplit &&) -> EnergySplit & = delete;

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
  EnergySplit(){};
  virtual ~EnergySplit() = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
  // update material properties
  virtual void updateMaterialProperties(Real E, Real nu, bool plane_stress);

  // // compute strain energy density on quad
  // virtual void computePhiOnQuad(const Matrix<Real> & /*strain_quad*/,
  //                               Real & /*phi_quad*/) = 0;

  // // compute stress on quad
  // virtual void computeSigmaOnQuad(const Matrix<Real> & /*strain_quad*/,
  //                                 const Real & /*sigma_th*/,
  //                                 Matrix<Real> & /*sigma_plus*/,
  //                                 Matrix<Real> & /*sigma_minus*/) = 0;

  // // compute tangent moduli coefficients on quad
  // virtual void computeTangentCoefsOnQuad(const Matrix<Real> & /*strain_quad*/,
  //                                        const Real & /*g_d*/,
  //                                        Matrix<Real> & /*tangent*/) = 0;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// Young's modulus
  Real E{0.};

  /// Poisson ratio
  Real nu{0.};

  /// Finite deformation
  bool plane_stress{false};

  /// Lame's first parameter
  Real lambda{0.};

  /// Lame's second paramter
  Real mu{0.};
};
} // namespace akantu

#endif
