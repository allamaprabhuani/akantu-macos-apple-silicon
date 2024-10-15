/**
 * Copyright (©) 2011-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_named_argument.hh"
#include "model.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_STRUCTURAL_MECHANICS_MODEL_HH_
#define AKANTU_STRUCTURAL_MECHANICS_MODEL_HH_

/* -------------------------------------------------------------------------- */
namespace akantu {
class Material;
class MaterialSelector;
class DumperIOHelper;
class NonLocalManager;
template <ElementKind kind, class IntegrationOrderFunctor>
class IntegratorGauss;
template <ElementKind kind> class ShapeStructural;
} // namespace akantu

namespace akantu {

struct StructuralMaterial {
  Real E{0};
  Real A{1};
  Real I{0};
  Real Iz{0};
  Real Iy{0};
  Real GJ{0};
  Real rho{0};
  Real t{0};
  Real nu{0};
};

class StructuralMechanicsModel : public Model {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  using MyFEEngineType =
      FEEngineTemplate<IntegratorGauss, ShapeStructural, _ek_structural>;

  StructuralMechanicsModel(Mesh & mesh, Int dim = _all_dimensions,
                           const ID & id = "structural_mechanics_model");

  ~StructuralMechanicsModel() override;

  /// Init full model
  void initFullImpl(const ModelOptions & options) override;

  /// Init boundary FEEngine
  void initFEEngineBoundary() override;

  /* ------------------------------------------------------------------------ */
  /* Virtual methods from SolverCallback                                      */
  /* ------------------------------------------------------------------------ */
  /// get the type of matrix needed
  MatrixType getMatrixType(const ID & matrix_id) const override;

  /// callback to assemble a Matrix
  void assembleMatrix(const ID & matrix_id) override;

  /// callback to assemble a lumped Matrix
  void assembleLumpedMatrix(const ID & matrix_id) override;

  /// callback to assemble the residual (rhs)
  void assembleResidual() override;

  void assembleResidual(const ID & residual_part) override;

  bool canSplitResidual() const override { return true; }

  void afterSolveStep(bool converged) override;

  /// compute kinetic energy
  Real getKineticEnergy();

  /// compute potential energy
  Real getPotentialEnergy();

  /// compute the specified energy
  Real getEnergy(const ID & energy);

  /**
   * \brief This function computes the an approximation of the lumped mass.
   *
   * The mass is computed by looping over all beams and computing their mass.
   * The mass of a single beam is computed by the (initial) length of the beam,
   * its cross sectional area and its density.
   * The beam mass is then equaly distributed among the two nodes.
   *
   * For computing the rotational inertia, the function assumes that the mass of
   * a node is uniformaly distributed inside a disc (2D) or a sphere (3D). The
   * size of that disc, depends on the volume of the beam.
   *
   * Note that the computation of the mass is not unambigius.
   * The reason for this is, that the units of `StructralMaterial::rho` are not
   * clear. By default the function assumes that its unit are 'Mass per Volume'.
   * However, this makes the computed mass different than the consistent mass,
   * which seams to assume that its units are 'mass per unit length'.
   * The main difference between thge two are not the values, but that the
   * first version depends on `StructuralMaterial::A` while the later does not.
   * By defining the macro `AKANTU_STRUCTURAL_MECHANICS_CONSISTENT_LUMPED_MASS`
   * the function will compute the mass in a way that is consistent with the
   * consistent mass matrix.
   *
   * \note	The lumped mass is not stored inside the DOFManager.
   *
   * \param  ghost_type 	Should ghost types be computed.
   */
  void assembleLumpedMassMatrix();

  /* ------------------------------------------------------------------------ */
  /* Virtual methods from Model                                               */
  /* ------------------------------------------------------------------------ */
protected:
  /// get some default values for derived classes
  std::tuple<ID, TimeStepSolverType>
  getDefaultSolverID(const AnalysisMethod & method) override;

  ModelSolverOptions
  getDefaultSolverOptions(const TimeStepSolverType & type) const override;

  static UInt getNbDegreeOfFreedom(ElementType type);

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
  void initSolver(TimeStepSolverType time_step_solver_type) override;

  /// initialize the model
  void initModel() override;

  /// compute the stresses per elements
  void computeStresses();

  /// compute the nodal forces
  void assembleInternalForce();

  /// compute the nodal forces for an element type
  void assembleInternalForce(ElementType type, GhostType gt);

  /// assemble the stiffness matrix
  void assembleStiffnessMatrix();

  /// assemble the mass matrix for consistent mass resolutions
  void assembleMassMatrix();

protected:
  /// assemble the mass matrix for either _ghost or _not_ghost elements
  void assembleMassMatrix(GhostType ghost_type);

  /// computes rho
  void computeRho(Array<Real> & rho, ElementType type, GhostType ghost_type);

  /// finish the computation of residual to solve in increment
  void updateResidualInternal();

  /* ------------------------------------------------------------------------ */
private:
  template <ElementType type> void assembleStiffnessMatrix();
  template <ElementType type> void computeStressOnQuad();
  template <ElementType type>
  void computeTangentModuli(Array<Real> & tangent_moduli);

  /* ------------------------------------------------------------------------ */
  /* Dumpable interface                                                       */
  /* ------------------------------------------------------------------------ */
public:
  std::shared_ptr<dumpers::Field>
  createNodalFieldReal(const std::string & field_name,
                       const std::string & group_name,
                       bool padding_flag) override;

  std::shared_ptr<dumpers::Field>
  createNodalFieldBool(const std::string & field_name,
                       const std::string & group_name,
                       bool padding_flag) override;

  std::shared_ptr<dumpers::Field>
  createElementalField(const std::string & field_name,
                       const std::string & group_name, bool padding_flag,
                       Int spatial_dimension, ElementKind kind) override;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// set the value of the time step
  void setTimeStep(Real time_step, const ID & solver_id = "") override;

  /// get the StructuralMechanicsModel::displacement vector
  AKANTU_GET_MACRO(Displacement, *displacement_rotation, Array<Real> &);

  /// get the StructuralMechanicsModel::velocity vector
  AKANTU_GET_MACRO(Velocity, *velocity, Array<Real> &);

  /// get    the    StructuralMechanicsModel::acceleration    vector,   updated
  /// by
  /// StructuralMechanicsModel::updateAcceleration
  AKANTU_GET_MACRO(Acceleration, *acceleration, Array<Real> &);

  /// get the StructuralMechanicsModel::external_force vector
  AKANTU_GET_MACRO(ExternalForce, *external_force, Array<Real> &);

  /// get the StructuralMechanicsModel::internal_force vector (boundary forces)
  AKANTU_GET_MACRO(InternalForce, *internal_force, Array<Real> &);

  /// get the StructuralMechanicsModel::boundary vector
  AKANTU_GET_MACRO(BlockedDOFs, *blocked_dofs, Array<bool> &);

  /**
   * Returns a const reference to the array that stores the lumped mass.
   *
   * The returned array has dimension `N x d` where `N` is the number of nodes
   * and `d`, is the number of degrees of freedom per node.
   */
  inline const Array<Real> & getLumpedMass() const {
    if (this->mass == nullptr) {
      AKANTU_EXCEPTION("The pointer to the mass was not allocated.");
    };
    return *(this->mass);
  };

  // These function is an alias, for compability with the solid mechanics
  inline const Array<Real> & getMass() const { return this->getLumpedMass(); }

  // Creates the array for storing the mass
  bool allocateLumpedMassArray();

  /**
   * Tests if *this has a lumped mass pointer.
   */
  inline bool hasLumpedMass() const { return (this->mass != nullptr); };

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(RotationMatrix, rotation_matrix, Real);

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Stress, stress, Real);

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(ElementMaterial, element_material, UInt);

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(Set_ID, set_ID, UInt);

  /**
   * \brief This function adds the `StructuralMaterial` material to the list of
   * materials managed by *this.
   *
   * It is important that this function might invalidate all references to
   * structural materials, that were previously optained by `getMaterial()`.
   *
   * \param  material The new material.
   * \param  name the name associated to this material
   *
   * \return The ID of the material that was added.
   *
   * \note The return type is is new.
   */
  UInt addMaterial(StructuralMaterial & material, const ID & name = "");

  const StructuralMaterial &
  getMaterialByElement(const Element & element) const;

  /**
   * \brief Returns the ith material of *this.
   * \param material_index The index of requested material
   */
  const StructuralMaterial & getMaterial(UInt material_index) const;

  const StructuralMaterial & getMaterial(const ID & name) const;

  /**
   * \brief Returns the number of the different materials inside *this.
   */
  UInt getNbMaterials() const { return materials.size(); }

  /* ------------------------------------------------------------------------ */
  /* Boundaries (structural_mechanics_model_boundary.cc)                      */
  /* ------------------------------------------------------------------------ */
public:
  /// Compute Linear load function set in global axis
  void computeForcesByGlobalTractionArray(const Array<Real> & traction_global,
                                          ElementType type);

  /// Compute Linear load function set in local axis
  void computeForcesByLocalTractionArray(const Array<Real> & tractions,
                                         ElementType type);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// time step
  Real time_step;

  /// conversion coefficient form force/mass to acceleration
  Real f_m2a;

  /// displacements array
  std::unique_ptr<Array<Real>> displacement_rotation;

  /// velocities array
  std::unique_ptr<Array<Real>> velocity;

  /// accelerations array
  std::unique_ptr<Array<Real>> acceleration;

  /// forces array
  std::unique_ptr<Array<Real>> internal_force;

  /// forces array
  std::unique_ptr<Array<Real>> external_force;

  /**
   * \brief	This is the "lumped" mass array.
   *
   * It is a bit special, since it is not a one dimensional array, bit it is
   * actually a matrix. The number of rows equals the number of nodes. The
   * number of colums equals the number of degrees of freedoms per node. This
   * layout makes the thing a bit more simple.
   *
   * Note that it is only allocated in case, the "Lumped" mode is enabled.
   */
  std::unique_ptr<Array<Real>> mass;

  /// boundaries array
  std::unique_ptr<Array<bool>> blocked_dofs;

  /// stress array
  ElementTypeMapArray<Real> stress;

  ElementTypeMapArray<UInt> element_material;

  // Define sets of beams
  ElementTypeMapArray<UInt> set_ID;

  /// number of degre of freedom
  Int nb_degree_of_freedom;

  // Rotation matrix
  ElementTypeMapArray<Real> rotation_matrix;

  // /// analysis method check the list in akantu::AnalysisMethod
  // AnalysisMethod method;

  /// flag defining if the increment must be computed or not
  bool increment_flag;

  bool need_to_reassemble_mass{true};
  bool need_to_reassemble_stiffness{true};
  bool need_to_reassemble_lumped_mass{true};

  /* ------------------------------------------------------------------------ */
  std::vector<StructuralMaterial> materials;
  std::map<std::string, UInt> materials_names_to_id;
};

} // namespace akantu

#include "structural_mechanics_model_inline_impl.hh"

#endif /* AKANTU_STRUCTURAL_MECHANICS_MODEL_HH_ */
