
/* -------------------------------------------------------------------------- */
#include "energy_split.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
void EnergySplit::updateMaterialProperties(Real E, Real nu, bool plane_stress) {
  this->E = E;
  this->nu = nu;
  this->plane_stress = plane_stress;
  lambda = nu * E / ((1 + nu) * (1 - 2 * nu));
  if (plane_stress) {
    lambda = nu * E / ((1 + nu) * (1 - nu));
  }
  mu = E / (2 * (1 + nu));
}

} // namespace akantu
