/**
 * Copyright (©) 2017-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_common.hh"
#include "aka_named_argument.hh"
#include "model_solver.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MODEL_OPTIONS_HH_
#define AKANTU_MODEL_OPTIONS_HH_

namespace akantu {

namespace {
  DECLARE_NAMED_ARGUMENT(analysis_method);
  DECLARE_NAMED_ARGUMENT(solver_options);
} // namespace

class ModelOptions {
public:
  explicit ModelOptions(
      AnalysisMethod analysis_method = _explicit_lumped_mass,
      ModelSolverOptions solver_options = detail::_default_solver_options)
      : analysis_method(analysis_method),
        solver_options(std::move(solver_options)) {}

  template <typename... pack>
  ModelOptions(use_named_args_t /*unused*/, pack &&... _pack)
      : ModelOptions(OPTIONAL_NAMED_ARG(analysis_method, _explicit_lumped_mass),
                     OPTIONAL_NAMED_ARG(solver_options,
                                        detail::_default_solver_options)) {}

  virtual ~ModelOptions() = default;

  AnalysisMethod analysis_method;
  ModelSolverOptions solver_options;
};

#ifdef AKANTU_SOLID_MECHANICS
/* -------------------------------------------------------------------------- */
using SolidMechanicsModelOptions = ModelOptions;
#endif

/* -------------------------------------------------------------------------- */
#ifdef AKANTU_COHESIVE_ELEMENT
namespace {
  DECLARE_NAMED_ARGUMENT(is_extrinsic);
}
/* -------------------------------------------------------------------------- */
struct SolidMechanicsModelCohesiveOptions : public SolidMechanicsModelOptions {
  SolidMechanicsModelCohesiveOptions(
      AnalysisMethod analysis_method = _explicit_lumped_mass,
      ModelSolverOptions solver_options = detail::_default_solver_options,
      bool extrinsic = false)
      : SolidMechanicsModelOptions(analysis_method, solver_options),
        is_extrinsic(extrinsic) {}

  template <typename... pack>
  SolidMechanicsModelCohesiveOptions(use_named_args_t /*unused*/,
                                     pack &&... _pack)
      : SolidMechanicsModelCohesiveOptions(
            OPTIONAL_NAMED_ARG(analysis_method, _explicit_lumped_mass),
            OPTIONAL_NAMED_ARG(solver_options, detail::_default_solver_options),
            OPTIONAL_NAMED_ARG(is_extrinsic, false)) {}

  bool is_extrinsic{false};
};

#endif

#ifdef AKANTU_DIFFUSION
/* -------------------------------------------------------------------------- */
using HeatTransferModelOptions = ModelOptions;
using DiffusionModelOptions = ModelOptions;
#endif

#ifdef AKANTU_PHASE_FIELD
/* -------------------------------------------------------------------------- */
using PhaseFieldModelOptions = ModelOptions;
#endif

/* -------------------------------------------------------------------------- */
#ifdef AKANTU_EMBEDDED

namespace {
  DECLARE_NAMED_ARGUMENT(init_intersections, non_linear_solver);
}
/* -------------------------------------------------------------------------- */
struct EmbeddedInterfaceModelOptions : SolidMechanicsModelOptions {
  /**
   * @brief Constructor for EmbeddedInterfaceModelOptions
   * @param analysis_method see SolidMechanicsModelOptions
   * @param init_intersections compute intersections
   */
  EmbeddedInterfaceModelOptions(
      AnalysisMethod analysis_method = _explicit_lumped_mass,
      ModelSolverOptions solver_options = detail::_default_solver_options,
      bool init_intersections = true)
      : SolidMechanicsModelOptions(analysis_method, solver_options),
        has_intersections(init_intersections) {}

  template <typename... pack>
  EmbeddedInterfaceModelOptions(use_named_args_t /*unused*/, pack &&... _pack)
      : EmbeddedInterfaceModelOptions(
            OPTIONAL_NAMED_ARG(analysis_method, _explicit_lumped_mass),
            OPTIONAL_NAMED_ARG(solver_options, detail::_default_solver_options),
            OPTIONAL_NAMED_ARG(init_intersections, true)) {}

  /// Should consider reinforcements
  bool has_intersections;
};
#endif

#ifdef AKANTU_CONTACT_MECHANICS
/* -------------------------------------------------------------------------- */
using ContactMechanicsModelOptions = ModelOptions;
#endif

#ifdef AKANTU_MODEL_COUPLERS
/* -------------------------------------------------------------------------- */
using CouplerSolidContactOptions = ModelOptions;

#ifdef AKANTU_COHESIVE_ELEMENT
/* -------------------------------------------------------------------------- */
struct CouplerSolidCohesiveContactOptions : public ModelOptions {

  CouplerSolidCohesiveContactOptions(
      AnalysisMethod analysis_method = _explicit_lumped_mass,
      ModelSolverOptions solver_options = detail::_default_solver_options,
      bool extrinsic = false)
      : ModelOptions(analysis_method, solver_options), is_extrinsic(extrinsic) {
  }

  template <typename... pack>
  CouplerSolidCohesiveContactOptions(use_named_args_t /*unused*/,
                                     pack &&... _pack)
      : CouplerSolidCohesiveContactOptions(
            OPTIONAL_NAMED_ARG(analysis_method, _explicit_lumped_mass),
            OPTIONAL_NAMED_ARG(solver_options, detail::_default_solver_options),
            OPTIONAL_NAMED_ARG(is_extrinsic, false)) {}

  bool is_extrinsic{false};
};
#endif

#endif

} // namespace akantu

#endif /* AKANTU_MODEL_OPTIONS_HH_ */
