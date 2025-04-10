/**
 * Copyright (©) 2019-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#include "test_gtest_utils.hh"
/* -------------------------------------------------------------------------- */
#include <aka_common.hh>
#include <algorithm>
#include <communicator.hh>
#include <dof_manager.hh>
#include <element.hh>
#include <mesh.hh>
#include <mesh_partition_scotch.hh>
#include <mesh_utils.hh>
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>
#include <memory>
#include <numeric>
#include <string>
#include <type_traits>
#include <unordered_map>
/* -------------------------------------------------------------------------- */

namespace akantu {
enum DOFManagerType { _dmt_default, _dmt_petsc };
} // namespace akantu

AKANTU_ENUM_HASH(DOFManagerType)

using namespace akantu;

// defined as struct to get there names in gtest outputs
struct dof_manager_default_t
    : public std::integral_constant<DOFManagerType, _dmt_default> {};
struct dof_manager_petsc_t
    : public std::integral_constant<DOFManagerType, _dmt_petsc> {};

using dof_manager_types = ::testing::Types<
#ifdef AKANTU_USE_PETSC
    dof_manager_petsc_t,
#endif
    dof_manager_default_t>;

namespace std {

auto to_string(const DOFManagerType & type) -> std::string {
  std::unordered_map<DOFManagerType, std::string> map{
#ifdef AKANTU_USE_PETSC
      {_dmt_petsc, "petsc"},
#endif
      {_dmt_default, "default"},
  };
  return map.at(type);
}

} // namespace std

/* -------------------------------------------------------------------------- */
using namespace akantu;
/* -------------------------------------------------------------------------- */

template <class T> class DOFManagerFixture : public ::testing::Test {
public:
  constexpr static DOFManagerType type = T::value;
  constexpr static Int dim = 3;
  void SetUp() override {
    mesh = std::make_unique<Mesh>(this->dim);

    auto & communicator = Communicator::getStaticCommunicator();
    prank = communicator.whoAmI();
    psize = communicator.getNbProc();

    if (prank == 0) {
      mesh->read("mesh.msh");
    }
    mesh->distribute();

    nb_nodes = this->mesh->getNbNodes();
    nb_total_nodes = this->mesh->getNbGlobalNodes();

    auto && range_nodes = arange(nb_nodes);
    nb_pure_local =
        std::accumulate(range_nodes.begin(), range_nodes.end(), 0,
                        [&](auto && init, auto && val) {
                          return init + mesh->isLocalOrMasterNode(val);
                        });

    dof_manager = this->allocDOFManager();
  }
  void TearDown() override {
    dof_manager.reset();
    mesh.reset();
    dof1.reset();
    dof2.reset();
  }

  auto allocDOFManager() -> decltype(auto) {
    return DOFManagerFactory::getInstance().allocate(std::to_string(type),
                                                     *mesh, "dof_manager");
  }

  void registerDOFs(DOFSupportType dst1, DOFSupportType dst2) {
    auto n1 = dst1 == _dst_nodal ? nb_nodes : nb_pure_local;
    this->dof1 = std::make_unique<Array<Real>>(n1, 3);

    this->dof_manager->registerDOFs("dofs1", *this->dof1, dst1);

    EXPECT_EQ(dof_manager->getResidual().size(), nb_total_nodes * 3);

    auto n2 = dst2 == _dst_nodal ? nb_nodes : nb_pure_local;
    this->dof2 = std::make_unique<Array<Real>>(n2, 5);

    this->dof_manager->registerDOFs("dofs2", *this->dof2, dst2);

    EXPECT_EQ(dof_manager->getResidual().size(), nb_total_nodes * 8);
  }

protected:
  Int nb_nodes{0}, nb_total_nodes{0}, nb_pure_local{0};
  int prank{0}, psize{1};
  std::unique_ptr<Mesh> mesh;
  std::unique_ptr<Array<Real>> dof1;
  std::unique_ptr<Array<Real>> dof2;

  std::unique_ptr<DOFManager> dof_manager;
};

TYPED_TEST_SUITE(DOFManagerFixture, dof_manager_types, );

/* -------------------------------------------------------------------------- */
TYPED_TEST(DOFManagerFixture, Construction) {
  // Construction in SetUp
}

/* -------------------------------------------------------------------------- */
TYPED_TEST(DOFManagerFixture, DoubleConstruction) {
  this->dof_manager = this->allocDOFManager();
}

/* -------------------------------------------------------------------------- */
TYPED_TEST(DOFManagerFixture, RegisterGenericDOF1) {
  Array<Real> dofs(this->nb_pure_local, 3);

  this->dof_manager->registerDOFs("dofs1", dofs, _dst_generic);
  EXPECT_GE(this->dof_manager->getResidual().size(), this->nb_total_nodes * 3);
}

/* -------------------------------------------------------------------------- */
TYPED_TEST(DOFManagerFixture, RegisterNodalDOF1) {
  Array<Real> dofs(this->nb_nodes, 3);
  this->dof_manager->registerDOFs("dofs1", dofs, _dst_nodal);
  EXPECT_GE(this->dof_manager->getResidual().size(), this->nb_total_nodes * 3);
}

/* -------------------------------------------------------------------------- */
TYPED_TEST(DOFManagerFixture, RegisterGenericDOF2) {
  this->registerDOFs(_dst_generic, _dst_generic);
}

/* -------------------------------------------------------------------------- */
TYPED_TEST(DOFManagerFixture, RegisterNodalDOF2) {
  this->registerDOFs(_dst_nodal, _dst_nodal);
}

/* -------------------------------------------------------------------------- */
TYPED_TEST(DOFManagerFixture, RegisterMixedDOF) {
  this->registerDOFs(_dst_nodal, _dst_generic);
}

/* -------------------------------------------------------------------------- */
TYPED_TEST(DOFManagerFixture, splitArrayPerDOFs) {
  Array<Real> dofs(this->nb_nodes, 1);

  auto && range = arange(this->nb_nodes);
  std::transform(range.begin(), range.end(), dofs.begin(), [&](auto n) {
    Idx gid =
        this->mesh->isLocalOrMasterNode(n) ? this->mesh->getNodeGlobalId(n) : 0;
    return gid;
  });

  this->dof_manager->registerDOFs("dofs1", dofs, _dst_nodal);

  auto & nodes = this->dof_manager->getNewLumpedMatrix("nodes");

  this->dof_manager->assembleToGlobalArray("dofs1", dofs, nodes, 1.);

  dofs.set(0.);

  this->dof_manager->getArrayPerDOFs("dofs1", nodes, dofs);

  for (auto && [i, node] : enumerate(dofs)) {
    if (not this->mesh->isPureGhostNode(i)) {
      EXPECT_EQ(node, this->mesh->getNodeGlobalId(i));
    }
  }
}

/* -------------------------------------------------------------------------- */
TYPED_TEST(DOFManagerFixture, AssembleVector) {
  this->registerDOFs(_dst_nodal, _dst_generic);

  this->dof_manager->getResidual().zero();

  for (auto && [n, l] :
       enumerate(make_view(*this->dof1, this->dof1->getNbComponent()))) {
    l.set(1. * this->mesh->isLocalOrMasterNode(n));
  }

  this->dof2->set(2.);

  this->dof_manager->assembleToResidual("dofs1", *this->dof1);
  this->dof_manager->assembleToResidual("dofs2", *this->dof2);

  this->dof1->set(0.);
  this->dof2->set(0.);

  const auto ref1 = Vector<Real>{1., 1., 1.};

  this->dof_manager->getArrayPerDOFs("dofs1", this->dof_manager->getResidual(),
                                     *this->dof1);
  for (auto && [node, dof] :
       enumerate(make_view(*this->dof1, this->dof1->getNbComponent()))) {
    if (not this->mesh->isPureGhostNode(node)) {
      auto e = (dof - ref1).norm() / ref1.norm();
      EXPECT_NEAR(e, 0., 1e-14)
          << "[" << this->prank << "/" << this->psize << "] DOF: " << dof
          << " should be: " << ref1 << " for node: " << node;
    }
  }

  const auto ref2 = Vector<Real>{2., 2., 2., 2., 2.};
  this->dof_manager->getArrayPerDOFs("dofs2", this->dof_manager->getResidual(),
                                     *this->dof2);
  for (auto && [node, dof] :
       enumerate(make_view(*this->dof2, this->dof2->getNbComponent()))) {
    if (not this->mesh->isPureGhostNode(node)) {
      auto e = (dof - ref2).norm() / ref2.norm();
      EXPECT_NEAR(e, 0., 1e-14)
          << "[" << this->prank << "/" << this->psize << "] DOF: " << dof
          << " should be: " << ref2 << " for node: " << node;
    }
  }
}

/* -------------------------------------------------------------------------- */
TYPED_TEST(DOFManagerFixture, AssembleMatrixNodal) {
  this->registerDOFs(_dst_nodal, _dst_nodal);

  auto && K = this->dof_manager->getNewMatrix("K", _symmetric);
  K.zero();

  auto && elemental_matrix = std::make_unique<Array<Real>>(
      this->mesh->getNbElement(this->dim), 8 * 3 * 8 * 3);

  for (auto && m : make_view(*elemental_matrix, 8 * 3, 8 * 3)) {
    m.set(1.);
  }

  this->dof_manager->assembleElementalMatricesToMatrix(
      "K", "dofs1", *elemental_matrix, _hexahedron_8);

  elemental_matrix = std::make_unique<Array<Real>>(
      this->mesh->getNbElement(this->dim), 8 * 5 * 8 * 5);

  for (auto && m : make_view(*elemental_matrix, 8 * 5, 8 * 5)) {
    m.set(2.);
  }

  this->dof_manager->assembleElementalMatricesToMatrix(
      "K", "dofs2", *elemental_matrix, _hexahedron_8);

  K.saveMatrix("K_" + std::to_string(this->type) + ".mtx");
}
