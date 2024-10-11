/**
 * Copyright (©) 2012-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "mesh_iterators.hh"
#include "solid_mechanics_model_cohesive.hh"
/* -------------------------------------------------------------------------- */
#include <iostream>
/* -------------------------------------------------------------------------- */

using namespace akantu;

static void updateDisplacement(SolidMechanicsModelCohesive & model,
                               const ElementGroup & group, Real increment);

/* -------------------------------------------------------------------------- */
int main(int argc, char * argv[]) {
  initialize("material.dat", argc, argv);

  const Int spatial_dimension = 2;
  const Int max_steps = 350;

  Mesh mesh(spatial_dimension);
  mesh.read("triangle.msh");

  SolidMechanicsModelCohesive model(mesh);
  // To restric the insertion to the range [-0.26, -0.24] in the x direction
  model.getElementInserter().setLimit(_x, -0.26, -0.24);

  /// model initialization
  // _is_extrinsic = false for intrinsic
  model.initFull(_analysis_method = _explicit_lumped_mass,
                 _is_extrinsic = false);

  Real time_step = model.getStableTimeStep() * 0.8;
  model.setTimeStep(time_step);
  std::cout << "Time step: " << time_step << "\n";

  Array<bool> & boundary = model.getBlockedDOFs();

  // Int nb_nodes = mesh.getNbNodes();

  /// boundary conditions
  boundary.set(true);

  model.setBaseName("intrinsic");
  model.addDumpFieldVector("displacement");
  model.addDumpField("velocity");
  model.addDumpField("acceleration");
  model.addDumpField("grad_u");
  model.addDumpField("external_force");
  model.addDumpField("internal_force");
  model.dump();

  /// update displacement
  auto && elements = mesh.createElementGroup("diplacement");
  Vector<Real> barycenter(spatial_dimension);

  for_each_element(
      mesh,
      [&](auto && el) {
        mesh.getBarycenter(el, barycenter);
        if (barycenter(_x) > -0.25) {
          elements.add(el, true);
        }
      },
      _element_kind = _ek_regular);

  Real increment = 0.01;

  updateDisplacement(model, elements, increment);

  /// Main loop
  for (Int s = 1; s <= max_steps; ++s) {
    model.solveStep();

    updateDisplacement(model, elements, increment);
    if (s % 1 == 0) {
      model.dump();
      std::cout << "passing step " << s << "/" << max_steps << "\n";
    }
  }

  Real Ed = model.getEnergy("dissipated");
  Real Edt = 2 * sqrt(2);

  std::cout << Ed << " " << Edt << "\n";

  if (Ed < Edt * 0.999 || Ed > Edt * 1.001 || std::isnan(Ed)) {
    std::cout << "The dissipated energy is incorrect"
              << "\n";
    return -1;
  }

  return 0;
}

/* -------------------------------------------------------------------------- */
static void updateDisplacement(SolidMechanicsModelCohesive & model,
                               const ElementGroup & group, Real increment) {
  Array<Real> & displacement = model.getDisplacement();

  for (auto && node : group.getNodeGroup().getNodes()) {
    displacement(node, 0) += increment;
  }
}
