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
#include "dof_synchronizer.hh"
#include "element_synchronizer.hh"
#include "mesh.hh"
#include "mesh_accessor.hh"
#include "mesh_partition_scotch.hh"
#include "sparse_matrix_aij.hh"
#include "sparse_solver_eigen.hh"
#include "terms_to_assemble.hh"
/* -------------------------------------------------------------------------- */
#include <iostream>
/* -------------------------------------------------------------------------- */

using namespace akantu;

/* -------------------------------------------------------------------------- */
void genMesh(Mesh & mesh, Int nb_nodes);
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
int main(int argc, char * argv[]) {
  initialize(argc, argv);
  const Int spatial_dimension = 1;
  const Int nb_global_dof = 11;
  const auto & comm = Communicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();

  Mesh mesh(spatial_dimension);

  if (prank == 0) {
    genMesh(mesh, nb_global_dof);
    RandomGenerator<Idx>::seed(1496137735);
  } else {
    RandomGenerator<Idx>::seed(2992275470);
  }

  mesh.distribute();
  UInt node = 0;
  for (auto pos : mesh.getNodes()) {
    std::cout << prank << " " << node << " pos: " << pos << " ["
              << mesh.getNodeGlobalId(node) << "] " << mesh.getNodeFlag(node)
              << std::endl;
    ++node;
  }
  Int nb_nodes = mesh.getNbNodes();

  DOFManagerDefault dof_manager(mesh, "test_dof_manager");

  Array<Real> x(nb_nodes);
  x.set(-1);
  dof_manager.registerDOFs("x", x, _dst_nodal);

  const auto & local_equation_number =
      dof_manager.getLocalEquationsNumbers("x");

  auto & A = dof_manager.getNewMatrix("A", _symmetric);

  Array<Real> b(nb_nodes);
  TermsToAssemble terms("x", "x");

  for (Int i = 0; i < nb_nodes; ++i) {
    if (dof_manager.isLocalOrMasterDOF(i)) {
      auto li = local_equation_number(i);
      auto gi = dof_manager.localToGlobalEquationNumber(li);
      terms(i, i) = 1. / (1. + gi);
    }
  }

  dof_manager.assemblePreassembledMatrix("A", terms);

  std::stringstream str;
  str << "Matrix_" << prank << ".mtx";
  A.saveMatrix(str.str());

  for (Int n = 0; n < nb_nodes; ++n) {
    b(n) = 1.;
  }

  SparseSolverEigen solver(dof_manager, "A");
  dof_manager.assembleToResidual("x", b);

  solver.solve();
  Array<Real> & y = dof_manager.getSolution("x");

  auto && check = [&](auto && xs) {
    debug::setDebugLevel(dblTest);
    std::cout << xs << std::endl;
    debug::setDebugLevel(dblWarning);

    Int d = 1.;
    for (auto x : xs) {
      if (std::abs(x - d) / d > 1e-15)
        AKANTU_EXCEPTION("Error in the solution: " << x << " != " << d << " ["
                                                   << (std::abs(x - d) / d)
                                                   << "].");
      ++d;
    }
  };

  if (psize > 1) {
    auto & sync =
        dynamic_cast<DOFManagerDefault &>(dof_manager).getSynchronizer();

    if (prank == 0) {
      Array<Real> x_gathered(dof_manager.getSystemSize());
      sync.gather(y, x_gathered);
      check(x_gathered);
    } else {
      sync.gather(x);
    }
  } else {
    check(y);
  }

  finalize();

  return 0;
}

/* -------------------------------------------------------------------------- */
void genMesh(Mesh & mesh, Int nb_nodes) {
  MeshAccessor mesh_accessor(mesh);
  Array<Real> & nodes = mesh_accessor.getNodes();
  Array<Idx> & conn = mesh_accessor.getConnectivity(_segment_2);

  nodes.resize(nb_nodes);

  for (Int n = 0; n < nb_nodes; ++n) {
    nodes(n, _x) = n * (1. / (nb_nodes - 1));
  }

  conn.resize(nb_nodes - 1);
  for (Int n = 0; n < nb_nodes - 1; ++n) {
    conn(n, 0) = n;
    conn(n, 1) = n + 1;
  }

  mesh_accessor.makeReady();
}
