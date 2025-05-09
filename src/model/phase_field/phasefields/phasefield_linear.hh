/* -------------------------------------------------------------------------- */
#include "energy_split.hh"
#include "phasefield.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_PHASEFIELD_LINEAR_HH__
#define __AKANTU_PHASEFIELD_LINEAR_HH__

namespace akantu {

template <Int dim, template <Int> class EnergySplit_>
  requires CanComputePhi<dim, EnergySplit_>
class PhaseFieldLinear : public PhaseField {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  PhaseFieldLinear(PhaseFieldModel & model, const ID & id = "");

  ~PhaseFieldLinear() override = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// compute the dissiapted energy
  void computeDissipatedEnergy(ElementType el_type) override;

  void
  computeDissipatedEnergyByElement(const Element & element,
                                   Vector<Real> & edis_on_quad_points) override;

  void afterSolveStep() override;

protected:
  void
  computeDissipatedEnergyByElement(ElementType type, Idx index,
                                   Vector<Real> & edis_on_quad_points) override;

  void computeDrivingForce(ElementType /*el_type*/,
                           GhostType /*ghost_type*/) override;

  void computeResidual(ElementType /*el_type*/, GhostType /*ghost_type*/);

  void applyPenalization(ElementType /*el_type*/, GhostType /*ghost_type*/);

  inline void computeDrivingForceOnQuad(const Real & /*phi_quad*/,
                                        Real & /*driving_force_quad*/,
                                        const Real & /*g_c_quad*/);

  inline void computeDamageEnergyDensityOnQuad(const Real & /*phi_quad*/,
                                               Real & /*dam_energy_quad*/);

  inline void
  computeDissipatedEnergyOnQuad(const Real & /*dam_quad*/,
                                const Vector<Real> & /*grad_d_quad */,
                                Real & /*energy*/, Real & /*g_c_quad*/);

public:
  void updateInternalParameters() override;

  void initPhaseField() override;

private:
  // irreversibility tolerance
  Real tol_ir;

  // recovery tolerance
  Real tol_rec;

  // irreversibility penalization parameter
  Real gamma;

  // recovery penalization parameter
  Real rho_rec;

  /// energy split
  std::shared_ptr<EnergySplit_<dim>> energy_split;
};

} // namespace akantu

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
#include "phasefield_linear_inline_impl.hh"

#endif
