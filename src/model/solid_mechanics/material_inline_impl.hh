/**
 * Copyright (©) 2010-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * This file is part of Akantu
 *
 * Akantu is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 */

/* -------------------------------------------------------------------------- */
#include "integration_point.hh"
#include "material.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

// #ifndef __AKANTU_MATERIAL_INLINE_IMPL_CC__
// #define __AKANTU_MATERIAL_INLINE_IMPL_CC__

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int dim, typename D1, typename D2>
constexpr inline void Material::gradUToF(const Eigen::MatrixBase<D1> & grad_u,
                                         Eigen::MatrixBase<D2> & F) {
  assert(F.size() >= grad_u.size() && grad_u.size() == dim * dim &&
         "The dimension of the tensor F should be greater or "
         "equal to the dimension of the tensor grad_u.");

  F.setIdentity();
  F.template block<dim, dim>(0, 0) += grad_u;
}

/* -------------------------------------------------------------------------- */
template <Int dim, typename D1>
constexpr inline decltype(auto)
Material::gradUToF(const Eigen::MatrixBase<D1> & grad_u) {
  Matrix<Real, dim, dim> F;
  gradUToF<dim>(grad_u, F);
  return F;
}

/* -------------------------------------------------------------------------- */
template <Int dim, typename D1, typename D2, typename D3>
constexpr inline void Material::StoCauchy(const Eigen::MatrixBase<D1> & F,
                                          const Eigen::MatrixBase<D2> & S,
                                          Eigen::MatrixBase<D3> & sigma,
                                          const Real & C33) {
  Real J = F.determinant() * std::sqrt(C33);

  Matrix<Real, dim, dim> F_S;
  F_S = F * S;
  Real constant = (J < std::numeric_limits<Real>::epsilon()) ? 0. : 1. / J;
  sigma = constant * F_S * F.transpose();
}

/* -------------------------------------------------------------------------- */
template <Int dim, typename D1, typename D2>
constexpr inline decltype(auto)
Material::StoCauchy(const Eigen::MatrixBase<D1> & F,
                    const Eigen::MatrixBase<D2> & S, const Real & C33) {
  Matrix<Real, dim, dim> sigma;
  Material::StoCauchy<dim>(F, S, sigma, C33);
  return sigma;
}
/* -------------------------------------------------------------------------- */
template <typename D1, typename D2>
constexpr inline void Material::rightCauchy(const Eigen::MatrixBase<D1> & F,
                                            Eigen::MatrixBase<D2> & C) {
  C = F.transpose() * F;
}

/* -------------------------------------------------------------------------- */
template <Int dim, typename D>
constexpr inline decltype(auto)
Material::rightCauchy(const Eigen::MatrixBase<D> & F) {
  Matrix<Real, dim, dim> C;
  rightCauchy(F, C);
  return C;
}

/* -------------------------------------------------------------------------- */
template <typename D1, typename D2>
constexpr inline void Material::leftCauchy(const Eigen::MatrixBase<D1> & F,
                                           Eigen::MatrixBase<D2> & B) {
  B = F * F.transpose();
}

/* -------------------------------------------------------------------------- */
template <Int dim, typename D>
constexpr inline decltype(auto)
Material::leftCauchy(const Eigen::MatrixBase<D> & F) {
  Matrix<Real, dim, dim> B;
  rightCauchy(F, B);
  return B;
}

/* -------------------------------------------------------------------------- */
template <Int dim, typename D1, typename D2>
constexpr inline void
Material::gradUToEpsilon(const Eigen::MatrixBase<D1> & grad_u,
                         Eigen::MatrixBase<D2> & epsilon) {
  epsilon = .5 * (grad_u.transpose() + grad_u);
}

/* -------------------------------------------------------------------------- */
template <Int dim, typename D1>
inline decltype(auto) constexpr Material::gradUToEpsilon(
    const Eigen::MatrixBase<D1> & grad_u) {
  Matrix<Real, dim, dim> epsilon;
  Material::gradUToEpsilon<dim>(grad_u, epsilon);
  return epsilon;
}

/* -------------------------------------------------------------------------- */
template <Int dim, typename D1, typename D2>
constexpr inline void Material::gradUToE(const Eigen::MatrixBase<D1> & grad_u,
                                         Eigen::MatrixBase<D2> & E) {
  E = (grad_u.transpose() * grad_u + grad_u.transpose() + grad_u) / 2.;
}

/* -------------------------------------------------------------------------- */
template <Int dim, typename D1>
constexpr inline decltype(auto)
Material::gradUToE(const Eigen::MatrixBase<D1> & grad_u) {
  Matrix<Real, dim, dim> E;
  gradUToE<dim>(grad_u, E);
  return E;
}

/* -------------------------------------------------------------------------- */
template <typename D1>
inline Real Material::stressToVonMises(const Eigen::MatrixBase<D1> & stress) {
  // compute deviatoric stress
  auto dim = stress.cols();
  auto && deviatoric_stress =
      stress - Matrix<Real>::Identity(dim, dim) * stress.trace() / 3.;

  // return Von Mises stress
  return std::sqrt(3. * deviatoric_stress.doubleDot(deviatoric_stress) / 2.);
}

/* -------------------------------------------------------------------------- */
template <Int dim, typename D1, typename D2>
constexpr inline void
Material::setCauchyStressMatrix(const Eigen::MatrixBase<D1> & S_t,
                                Eigen::MatrixBase<D2> & sigma) {
  sigma.zero();

  /// see Finite ekement formulations for large deformation dynamic analysis,
  /// Bathe et al. IJNME vol 9, 1975, page 364 ^t \f$\tau\f$
  for (Int i = 0; i < dim; ++i) {
    for (Int m = 0; m < dim; ++m) {
      for (Int n = 0; n < dim; ++n) {
        sigma(i * dim + m, i * dim + n) = S_t(m, n);
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
inline Int Material::getNbData(const Array<Element> & elements,
                               const SynchronizationTag & tag) const {
  if (tag == SynchronizationTag::_smm_stress) {
    return (this->isFiniteDeformation() ? 3 : 1) * spatial_dimension *
           spatial_dimension * Int(sizeof(Real)) *
           this->getModel().getNbIntegrationPoints(elements);
  }
  if (tag == SynchronizationTag::_smm_gradu) {
    return spatial_dimension * spatial_dimension * Int(sizeof(Real)) *
           this->getModel().getNbIntegrationPoints(elements);
  }
  return 0;
}

/* -------------------------------------------------------------------------- */
inline void Material::packData(CommunicationBuffer & buffer,
                               const Array<Element> & elements,
                               const SynchronizationTag & tag) const {
  if (tag == SynchronizationTag::_smm_stress) {
    if (this->isFiniteDeformation()) {
      packInternalFieldHelper(*piola_kirchhoff_2, buffer, elements);
      packInternalFieldHelper(gradu, buffer, elements);
    }
    packInternalFieldHelper(stress, buffer, elements);
  }

  if (tag == SynchronizationTag::_smm_gradu) {
    packInternalFieldHelper(gradu, buffer, elements);
  }
}

/* -------------------------------------------------------------------------- */
inline void Material::unpackData(CommunicationBuffer & buffer,
                                 const Array<Element> & elements,
                                 const SynchronizationTag & tag) {
  if (tag == SynchronizationTag::_smm_stress) {
    if (this->isFiniteDeformation()) {
      unpackInternalFieldHelper(*piola_kirchhoff_2, buffer, elements);
      unpackInternalFieldHelper(gradu, buffer, elements);
    }
    unpackInternalFieldHelper(stress, buffer, elements);
  }

  if (tag == SynchronizationTag::_smm_gradu) {
    unpackInternalFieldHelper(gradu, buffer, elements);
  }
}

/* -------------------------------------------------------------------------- */
/// return rho per quad point of element
Vector<Real> Material::getRho(Element & element) {
  auto & fem = getFEEngine();
  auto nb_quadrature_points = fem.getNbIntegrationPoints(element.type, element.ghost_type);

  Vector<Real> rhos(nb_quadrature_points);
  rhos.set(rho);
  return rhos;
}

} // namespace akantu

// #endif /* __AKANTU_MATERIAL_INLINE_IMPL_CC__ */
