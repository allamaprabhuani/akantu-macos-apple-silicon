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
 * @section DESCRIPTION
 *
 * This code is computing the gradient of a linear field and check that it gives
 * a constant result.  It also compute the gradient the  coordinates of the mesh
 * and check that it gives the identity
 *
 */

/* -------------------------------------------------------------------------- */
#include "test_fe_engine_fixture.hh"
/* -------------------------------------------------------------------------- */
#include <cstdlib>
#include <iostream>
/* -------------------------------------------------------------------------- */

using namespace akantu;

struct Support {
  const Array<Idx> & elem_filter;
  FEEngine & fem;
  int spatial_dimension;
  ElementType type;
  GhostType ghost_type;
};

struct TensorField {
  TensorField(const std::string & name, Support & support)
      : name(name), support(support) {}

  TensorField(const TensorField & f) = default;
  virtual void eval(Array<Real> & output) = 0;
  std::string name;
  Support & support;
  int getNbComponent() { throw; };

  TensorField & transpose() { throw; };
  TensorField & operator*(TensorField & f);
};

struct NodalTensorField : public TensorField {
  NodalTensorField(const std::string & name, Support & support)
      : TensorField(name, support) {}

  void eval(Array<Real> & output) override {
    support.fem.interpolateOnIntegrationPoints(
        nodal_field, output, nodal_field.getNbComponent(), support.type,
        support.ghost_type);
  }
  Array<Real> nodal_field;
};

struct IntegrationPointTensorField : public TensorField {
  IntegrationPointTensorField(const std::string & name, Support & support)
      : TensorField(name, support) {}

  void eval(Array<Real> &) override { throw; };
};

struct DotField : public TensorField {

  DotField(TensorField & f1, TensorField & f2)
      : TensorField(f1.name + "." + f2.name, f1.support), field1(f1),
        field2(f2) {}
  void eval(Array<Real> & output) override;

  TensorField & field1;
  TensorField & field2;
};

void DotField::eval(Array<Real> & output) {
  Array<Real> o1, o2;
  field1.eval(o1);
  field2.eval(o2);
  for (auto && [o, o1, o2] : zip(output, o1, o2)) {
    o = o1 * o2;
  }
};

TensorField & TensorField::operator*(TensorField & f) {
  return *new DotField{*this, f};
}

struct GradientOperator : public TensorField {

  GradientOperator(Support & support) : TensorField("gradient", support) {}
  void eval(Array<Real> & output) {
    const auto & shapes_derivatives =
        support.fem.getShapesDerivatives(support.type, support.ghost_type);

    output = shapes_derivatives;
  }
};

struct FieldIntegrator {

  static Array<Real> integrate(TensorField & field, Support & support) {
    auto nb_element = support.elem_filter.size();
    auto nb_nodes_per_element = Mesh::getNbNodesPerElement(support.type);
    Array<Real> res(nb_element,
                    nb_nodes_per_element * support.spatial_dimension,
                    "int_sigma_x_dphi_/_dX");

    auto nb_component = field.getNbComponent();

    Array<Real> field_eval;
    field.eval(field_eval);

    support.fem.integrate(field_eval, res, nb_component, support.type,
                          support.ghost_type, support.elem_filter);
    return res;
  }
};
/* -------------------------------------------------------------------------- */

TYPED_TEST(TestFEMFixture, evaluate_field_on_support) {
  this->fem->initShapeFunctions();
  const auto dim = TestFixture::dim;
  const auto type = TestFixture::type;

  Support support{empty_filter, *this->fem, dim, type, _not_ghost};
  NodalTensorField u("disp", support);

  Array<Real> u_d;
  // u.eval(u_d); // u evaluated on quadrature points of the support
  // std::cout << u_d << std::endl;
}

TYPED_TEST(TestFEMFixture, consistent_force) {
  this->fem->initShapeFunctions();
  const auto dim = TestFixture::dim;
  const auto type = TestFixture::type;

  Support support{empty_filter, *this->fem, dim, type, _not_ghost};
  GradientOperator grad{support};
  IntegrationPointTensorField stress("stress", support);

  // auto f = FieldIntegrator::integrate(stress * grad, support);
  // std::cout << f << std::endl;
}

TYPED_TEST(TestFEMFixture, stiffness_matrix) {
  this->fem->initShapeFunctions();
  const auto dim = TestFixture::dim;
  const auto type = TestFixture::type;

  Support support{empty_filter, *this->fem, dim, type, _not_ghost};
  GradientOperator grad{support};
  IntegrationPointTensorField elastic_tensor("C", support);

  // auto f = FieldIntegrator::integrate(grad.transpose() * elastic_tensor *
  // grad, support);
  // std::cout << f << std::endl;
}
