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
#include "aka_common.hh"
#include "integration_scheme.hh"
#include "non_linear_solver.hh"
#include "parsable.hh"
#include "solver_callback.hh"
#include "synchronizer_registry.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MODEL_SOLVER_HH_
#define AKANTU_MODEL_SOLVER_HH_

namespace akantu {
class Mesh;
class DOFManager;
class TimeStepSolver;
struct ModelSolverOptions;
} // namespace akantu

namespace akantu {

class ModelSolver : public Parsable,
                    public SolverCallback,
                    public SynchronizerRegistry {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  ModelSolver(Mesh & mesh, const ModelType & type, const ID & id);

  /// initialize the dof manager based on solver type passed in the input file
  std::shared_ptr<DOFManager>
  initDOFManager(const std::shared_ptr<DOFManager> & dof_manager = nullptr);
  /// initialize the dof manager based on the used chosen solver type
  std::shared_ptr<DOFManager> initDOFManager(const ID & dof_manager_type);

protected:
  /// initialize the dof manager based on the used chosen solver type
  std::shared_ptr<DOFManager> initDOFManager(const ParserSection & section,
                                             const ID & solver_type);

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// Callback for the model to instantiate the matrices when needed
  virtual void initSolver(TimeStepSolverType /*time_step_solver_type*/) {}

  /// get the section in the input file (if it exsits) corresponding to this
  /// model
  std::tuple<ParserSection, bool> getParserSection();

  /// solve a step using a given pre instantiated time step solver and
  /// non linear solver
  virtual void solveStep(const ID & solver_id = "");

  /// solve a step using a given pre instantiated time step solver and
  /// non linear solver with a user defined callback instead of the
  /// model itself /!\ This can mess up everything
  virtual void solveStep(SolverCallback & callback, const ID & solver_id = "");

  /// Initialize a time solver that can be used afterwards with its id
  void getNewSolver(const ID & solver_id, ModelSolverOptions solver_options);

  /// set an integration scheme for a given dof and a given solver
  void
  setIntegrationScheme(const ID & solver_id, const ID & dof_id,
                       const IntegrationSchemeType & integration_scheme_type,
                       IntegrationScheme::SolutionType solution_type =
                           IntegrationScheme::_not_defined);

  /// set an externally instantiated integration scheme
  void
  setIntegrationScheme(const ID & solver_id, const ID & dof_id,
                       std::unique_ptr<IntegrationScheme> & integration_scheme,
                       IntegrationScheme::SolutionType solution_type =
                           IntegrationScheme::_not_defined);

  /* ------------------------------------------------------------------------ */
  /* SolverCallback interface                                                 */
  /* ------------------------------------------------------------------------ */
public:
  /// Predictor interface for the callback
  void predictor() override;

  /// Corrector interface for the callback
  void corrector() override;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// Default time step solver to instantiate for this model
  [[nodiscard]] virtual TimeStepSolverType getDefaultSolverType() const;

  /// Default configurations for a given time step solver
  [[nodiscard]] virtual ModelSolverOptions
  getDefaultSolverOptions(const TimeStepSolverType & type) const;

  /// get access to the internal dof manager
  DOFManager & getDOFManager() { return *this->dof_manager; }

  /// get the time step of a given solver
  [[nodiscard]] Real getTimeStep(const ID & solver_id = "") const;
  /// set the time step of a given solver
  virtual void setTimeStep(Real time_step, const ID & solver_id = "");

  /// answer to the question "does the solver exists ?"
  [[nodiscard]] bool hasSolver(const ID & solver_id) const;

  /// changes the current default solver
  void setDefaultSolver(const ID & solver_id);

  /// is a default solver defined
  [[nodiscard]] bool hasDefaultSolver() const;

  /// is an integration scheme set for a given solver and a given dof
  [[nodiscard]] bool hasIntegrationScheme(const ID & solver_id,
                                          const ID & dof_id) const;

  [[nodiscard]] TimeStepSolver & getTimeStepSolver(const ID & solver_id = "");
  [[nodiscard]] NonLinearSolver & getNonLinearSolver(const ID & solver_id = "");

  [[nodiscard]] const TimeStepSolver &
  getTimeStepSolver(const ID & solver_id = "") const;
  [[nodiscard]] const NonLinearSolver &
  getNonLinearSolver(const ID & solver_id = "") const;

  /// get id of model
  AKANTU_GET_MACRO(ID, id, const ID &)

private:
  [[nodiscard]] TimeStepSolver & getSolver(const ID & solver_id);
  [[nodiscard]] const TimeStepSolver & getSolver(const ID & solver_id) const;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  ModelType model_type;

  /// Underlying dof_manager (the brain...)
  std::shared_ptr<DOFManager> dof_manager;

  ID id;

private:
  /// Underlying mesh
  Mesh & mesh;

  /// Default time step solver to use
  ID default_solver_id;
};

struct ModelSolverOptions {
  TimeStepSolverType timestep_solver_type{TimeStepSolverType::_not_defined};
  NonLinearSolverType non_linear_solver_type{NonLinearSolverType::_auto};
  SparseSolverType sparse_solver_type{SparseSolverType::_auto};
  std::map<ID, IntegrationSchemeType> integration_scheme_type{};
  std::map<ID, IntegrationScheme::SolutionType> solution_type{};

  void update(const ModelSolverOptions & other);
};

namespace detail {
const inline ModelSolverOptions _default_solver_options;
}

} // namespace akantu

#endif /* AKANTU_MODEL_SOLVER_HH_ */
