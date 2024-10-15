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
#include "aka_named_argument.hh"
#include "fe_engine.hh"
#include "mesh.hh"
#include "model_options.hh"
#include "model_solver.hh"

/* -------------------------------------------------------------------------- */
#include <typeindex>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MODEL_HH_
#define AKANTU_MODEL_HH_

namespace akantu {
class Parser;
class DumperIOHelper;
class DOFManager;
} // namespace akantu

/* -------------------------------------------------------------------------- */
namespace akantu {

class Model : public ModelSolver, public MeshEventHandler {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  /// Normal constructor where the DOFManager is created internally
  Model(Mesh & mesh, const ModelType & type, Int dim = _all_dimensions,
        const ID & id = "model");

  ~Model() override;

  /* ------------------------------------------------------------------------ */
  Model(const Model & /*other*/) = delete;
  Model(Model && /*other*/) = delete;
  Model & operator=(const Model & /*other*/) = delete;
  Model & operator=(Model && /*other*/) = delete;
  /* ------------------------------------------------------------------------ */

  using FEEngineMap = std::map<std::string, std::unique_ptr<FEEngine>>;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
protected:
  virtual void initFullImpl(const ModelOptions & options);

public:
  template <typename... pack>
  std::enable_if_t<are_named_argument<pack...>::value>
  initFull(pack &&... _pack) {
    switch (this->model_type) {
#ifdef AKANTU_SOLID_MECHANICS
    case ModelType::_solid_mechanics_model:
      this->initFullImpl(SolidMechanicsModelOptions{
          use_named_args, std::forward<decltype(_pack)>(_pack)...});
      break;
#endif
#ifdef AKANTU_COHESIVE_ELEMENT
    case ModelType::_solid_mechanics_model_cohesive:
      this->initFullImpl(SolidMechanicsModelCohesiveOptions{
          use_named_args, std::forward<decltype(_pack)>(_pack)...});
      break;
#endif
#ifdef AKANTU_DIFFUSION
    case ModelType::_heat_transfer_model:
      this->initFullImpl(HeatTransferModelOptions{
          use_named_args, std::forward<decltype(_pack)>(_pack)...});
      break;
#endif
#ifdef AKANTU_PHASE_FIELD
    case ModelType::_phase_field_model:
      this->initFullImpl(PhaseFieldModelOptions{
          use_named_args, std::forward<decltype(_pack)>(_pack)...});
      break;
#endif
#ifdef AKANTU_EMBEDDED
    case ModelType::_embedded_model:
      this->initFullImpl(EmbeddedInterfaceModelOptions{
          use_named_args, std::forward<decltype(_pack)>(_pack)...});
      break;
#endif
#ifdef AKANTU_CONTACT_MECHANICS
    case ModelType::_contact_mechanics_model:
      this->initFullImpl(ContactMechanicsModelOptions{
          use_named_args, std::forward<decltype(_pack)>(_pack)...});
      break;
#endif
#ifdef AKANTU_MODEL_COUPLERS
    case ModelType::_coupler_solid_contact:
      this->initFullImpl(CouplerSolidContactOptions{
          use_named_args, std::forward<decltype(_pack)>(_pack)...});
      break;
    case ModelType::_coupler_solid_cohesive_contact:
      this->initFullImpl(CouplerSolidCohesiveContactOptions{
          use_named_args, std::forward<decltype(_pack)>(_pack)...});
      break;
#endif
    default:
      this->initFullImpl(ModelOptions{use_named_args,
                                      std::forward<decltype(_pack)>(_pack)...});
    }
  }

  template <typename... pack>
  std::enable_if_t<not are_named_argument<pack...>::value>
  initFull(pack &&... _pack) {
    this->initFullImpl(std::forward<decltype(_pack)>(_pack)...);
  }

  /// initialize a new solver if needed
  void initNewSolver(const AnalysisMethod & method,
                     ModelSolverOptions options = ModelSolverOptions());

protected:
  /// get some default values for derived classes
  virtual std::tuple<ID, TimeStepSolverType>
  getDefaultSolverID(const AnalysisMethod & /*method*/) {
    return {"none", TimeStepSolverType::_not_defined};
  }

  virtual void initModel() {}

  virtual void initFEEngineBoundary();

  /// function to print the containt of the class
  void printself(std::ostream & /*stream*/, int /*indent*/ = 0) const override {
  }

public:
  /* ------------------------------------------------------------------------ */
  /* Access to the dumpable interface of the boundaries                       */
  /* ------------------------------------------------------------------------ */
  /// Dump the data for a given group
  void dumpGroup(const std::string & group_name);
  void dumpGroup(const std::string & group_name,
                 const std::string & dumper_name);
  /// Dump the data for all boundaries
  void dumpGroup();
  /// Set the directory for a given group
  void setGroupDirectory(const std::string & directory,
                         const std::string & group_name);
  /// Set the directory for all boundaries
  void setGroupDirectory(const std::string & directory);
  /// Set the base name for a given group
  void setGroupBaseName(const std::string & basename,
                        const std::string & group_name);
  /// Get the internal dumper of a given group
  DumperIOHelper & getGroupDumper(const std::string & group_name);

  /* ------------------------------------------------------------------------ */
  /* Function for non local capabilities                                      */
  /* ------------------------------------------------------------------------ */
  virtual void
  updateDataForNonLocalCriterion(ElementTypeMapReal & /*criterion*/) {
    AKANTU_TO_IMPLEMENT();
  }

protected:
  template <typename T>
  void allocNodalField(std::unique_ptr<Array<T>> & array, Int nb_component,
                       const ID & name) const;

  /* ------------------------------------------------------------------------ */
  /* Accessors */
  /* ------------------------------------------------------------------------ */
public:
  /// get the number of surfaces
  AKANTU_GET_MACRO(Mesh, mesh, Mesh &)

  /// synchronize the boundary in case of parallel run
  virtual void synchronizeBoundaries(){};

  /// return the fem object associated with a provided name
  [[nodiscard]] inline FEEngine & getFEEngine(const ID & name = "") const;

  /// return the fem boundary object associated with a provided name
  virtual FEEngine & getFEEngineBoundary(const ID & name = "");

  inline bool hasFEEngineBoundary(const ID & name = "");

  /// register a fem object associated with name
  template <typename FEEngineClass>
  inline void registerFEEngineObject(const ID & name, Mesh & mesh,
                                     Int spatial_dimension = _all_dimensions,
                                     bool do_not_precompute = false);

  /// unregister a fem object associated with name
  inline void unRegisterFEEngineObject(const ID & name);

  /// return the synchronizer registry
  SynchronizerRegistry & getSynchronizerRegistry();

  /// return the fem object associated with a provided name
  template <typename FEEngineClass>
  inline FEEngineClass & getFEEngineClass(const ID & name = "") const;

  /// return the fem boundary object associated with a provided name
  template <typename FEEngineClass>
  inline FEEngineClass & getFEEngineClassBoundary(const ID & name = "");

  /// Get the type of analysis method used
  AKANTU_GET_MACRO(AnalysisMethod, method, AnalysisMethod);

  /// return the dimension of the system space
  AKANTU_GET_MACRO(SpatialDimension, Model::spatial_dimension, Int);

  /* ------------------------------------------------------------------------ */
  /* Pack and unpack hexlper functions */
  /* ------------------------------------------------------------------------ */
public:
  [[nodiscard]] inline Int
  getNbIntegrationPoints(const Array<Element> & elements,
                         const ID & fem_id = ID()) const;

  /* ------------------------------------------------------------------------ */
  /* Dumpable interface (kept for convenience) and dumper relative functions  */
  /* ------------------------------------------------------------------------ */
  void setTextModeToDumper();

  virtual void addDumpGroupFieldToDumper(const std::string & field_id,
                                         std::shared_ptr<dumpers::Field> field,
                                         DumperIOHelper & dumper);

  virtual void addDumpField(const std::string & field_id);

  virtual void addDumpFieldVector(const std::string & field_id);

  virtual void addDumpFieldToDumper(const std::string & dumper_name,
                                    const std::string & field_id);

  virtual void addDumpFieldVectorToDumper(const std::string & dumper_name,
                                          const std::string & field_id);

  virtual void addDumpFieldTensorToDumper(const std::string & dumper_name,
                                          const std::string & field_id);

  virtual void addDumpFieldTensor(const std::string & field_id);

  virtual void setBaseName(const std::string & field_id);

  virtual void setBaseNameToDumper(const std::string & dumper_name,
                                   const std::string & basename);

  virtual void addDumpGroupField(const std::string & field_id,
                                 const std::string & group_name);

  virtual void addDumpGroupFieldToDumper(const std::string & dumper_name,
                                         const std::string & field_id,
                                         const std::string & group_name,
                                         ElementKind element_kind,
                                         bool padding_flag);

  virtual void addDumpGroupFieldToDumper(const std::string & dumper_name,
                                         const std::string & field_id,
                                         const std::string & group_name,
                                         Int spatial_dimension,
                                         ElementKind element_kind,
                                         bool padding_flag);

  virtual void removeDumpGroupField(const std::string & field_id,
                                    const std::string & group_name);
  virtual void removeDumpGroupFieldFromDumper(const std::string & dumper_name,
                                              const std::string & field_id,
                                              const std::string & group_name);

  virtual void addDumpGroupFieldVector(const std::string & field_id,
                                       const std::string & group_name);

  virtual void addDumpGroupFieldVectorToDumper(const std::string & dumper_name,
                                               const std::string & field_id,
                                               const std::string & group_name);

  virtual std::shared_ptr<dumpers::Field>
  createNodalFieldReal(const std::string & /*field_name*/,
                       const std::string & /*group_name*/,
                       bool /*padding_flag*/) {
    return nullptr;
  }

  virtual std::shared_ptr<dumpers::Field>
  createNodalFieldInt(const std::string & /*field_name*/,
                      const std::string & /*group_name*/,
                      bool /*padding_flag*/) {
    return nullptr;
  }

  virtual std::shared_ptr<dumpers::Field>
  createNodalFieldBool(const std::string & /*field_name*/,
                       const std::string & /*group_name*/,
                       bool /*padding_flag*/) {
    return nullptr;
  }

  virtual std::shared_ptr<dumpers::Field> createElementalField(
      const std::string & /*field_name*/, const std::string & /*group_name*/,
      bool /*padding_flag*/, Int /*spatial_dimension*/, ElementKind /*kind*/) {
    return nullptr;
  }

  void setDirectory(const std::string & directory);
  void setDirectoryToDumper(const std::string & dumper_name,
                            const std::string & directory);

  /* ------------------------------------------------------------------------ */
  virtual void dump(const std::string & dumper_name);
  virtual void dump(const std::string & dumper_name, Int step);
  virtual void dump(const std::string & dumper_name, Real time, Int step);
  /* ------------------------------------------------------------------------ */
  virtual void dump();
  virtual void dump(Int step);
  virtual void dump(Real time, Int step);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  friend std::ostream & operator<<(std::ostream & /*stream*/,
                                   const Model & /*this*/);

  /// analysis method check the list in akantu::AnalysisMethod
  AnalysisMethod method;

  /// Mesh
  Mesh & mesh;

  /// Spatial dimension of the problem
  Int spatial_dimension;

  /// the main fem object present in all  models
  FEEngineMap fems;

  /// default fem object
  std::string default_fem;

  /// parser to the pointer to use
  Parser & parser;

  /// default ElementKind for dumper
  ElementKind dumper_default_element_kind{_ek_regular};
};

/// standard output stream operator
inline std::ostream & operator<<(std::ostream & stream, const Model & _this) {
  _this.printself(stream);
  return stream;
}

} // namespace akantu

#include "model_inline_impl.hh"

#endif /* AKANTU_MODEL_HH_ */
