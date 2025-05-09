
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
  EnergySplit() {};
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
  // virtual void computeTangentCoefsOnQuad(const Matrix<Real> &
  // /*strain_quad*/,
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

/* -------------------------------------------------------------------------- */
template <Int dim, template <Int> class EnergySplit_>
concept CanComputePhi =
    requires(EnergySplit_<dim> energy_split, const Matrix<Real> & strain_quad,
             Real & phi_quad) {
      {
        energy_split.computePhiOnQuad(strain_quad, phi_quad)
      } -> std::same_as<void>;
    };

/* ------------------------------------------------------------------------ */
template <Int dim, template <Int> class EnergySplit_>
concept CanComputeSigma =
    requires(EnergySplit_<dim> energy_split, const Matrix<Real> & strain_quad,
             const Real & sigma_th, Matrix<Real> & sigma_plus,
             Matrix<Real> & sigma_minus) {
      {
        energy_split.computeSigmaOnQuad(strain_quad, sigma_th, sigma_plus,
                                        sigma_minus)
      } -> std::same_as<void>;
    };

/* -------------------------------------------------------------------------- */
template <Int dim, template <Int> class EnergySplit_>
concept CanComputeTangentCoefs =
    requires(EnergySplit_<dim> energy_split, const Matrix<Real> & strain_quad,
             const Real & g_d, Matrix<Real> & tangent) {
      {
        energy_split.computeTangentCoefsOnQuad(strain_quad, g_d, tangent)
      } -> std::same_as<void>;
    };

/* -------------------------------------------------------------------------- */
template <Int dim, template <Int> class EnergySplit_>
concept CanComputeSigmaAndTangent = CanComputeSigma<dim, EnergySplit_> &&
                                    CanComputeTangentCoefs<dim, EnergySplit_>;

} // namespace akantu
#endif
