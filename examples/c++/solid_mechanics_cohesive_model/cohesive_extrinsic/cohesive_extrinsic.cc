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
#include "solid_mechanics_model_cohesive.hh"
/* -------------------------------------------------------------------------- */
#include <iostream>
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char * argv[]) {
  initialize("material.dat", argc, argv);

  const Int dim = 2;
  const Int max_steps = 1000;

  Mesh mesh(dim);
  mesh.read("triangle.msh");

  SolidMechanicsModelCohesive model(mesh);

  /// model initialization
  model.initFull(_analysis_method = _explicit_lumped_mass,
                 _is_extrinsic = true);

  Real time_step = model.getStableTimeStep() * 0.05;
  model.setTimeStep(time_step);
  std::cout << "Time step: " << time_step << "\n";

  // updated the insertion limits
  CohesiveElementInserter & inserter = model.getElementInserter();
  inserter.setLimit(_y, 0.30, 0.20);
  model.updateAutomaticInsertion();

  Array<Real> & position = mesh.getNodes();
  Array<Real> & velocity = model.getVelocity();
  Array<bool> & boundary = model.getBlockedDOFs();
  Array<Real> & displacement = model.getDisplacement();

  // Int nb_nodes = mesh.getNbNodes();

  /// boundary conditions
  for (auto && [pos, boun] :
       zip(make_view(position, dim), make_view(boundary, dim))) {
    if (pos(_y) > 0.99 or pos(_y) < -0.99) {
      boun(_y) = true;
    }
    if (pos(_x) > 0.99 or pos(_x) < -0.99) {
      boun(_x) = true;
    }
  }

  model.setBaseName("extrinsic");
  model.addDumpFieldVector("displacement");
  model.addDumpField("velocity");
  model.addDumpField("acceleration");
  model.addDumpField("internal_force");
  model.addDumpField("stress");
  model.addDumpField("grad_u");
  model.dump();

  /// initial conditions
  Real loading_rate = 0.5;
  Real disp_update = loading_rate * time_step;
  for (auto && [pos, vel] :
       zip(make_view(position, dim), make_view(velocity, dim))) {
    vel(_y) = loading_rate * pos(_y);
  }

  /// Main loop
  for (Int s = 1; s <= max_steps; ++s) {

    /// update displacement on extreme nodes
    for (auto && [pos, disp] :
         zip(make_view(position, dim), make_view(displacement, dim))) {
      if (pos(_y) > 0.99 or pos(_y) < -0.99) {
        disp(_y) += disp_update * pos(_y);
      }
    }

    // check wether cohesive elements should be inserted
    model.checkCohesiveStress();
    model.solveStep();

    if (s % 10 == 0) {
      model.dump();

      std::cout << "passing step " << s << "/" << max_steps << "\n";
    }
  }

  Real Ed = model.getEnergy("dissipated");

  Real Edt = 200 * std::sqrt(2);
  std::cout << Ed << " " << Edt << "\n";
  if (Ed < Edt * 0.999 || Ed > Edt * 1.001 || std::isnan(Ed)) {
    std::cout << "The dissipated energy is incorrect"
              << "\n";
    return -1;
  }

  return 0;
}
