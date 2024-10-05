/**
 * Copyright (©) 2018-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include <aka_common.hh>
#include <aka_random_generator.hh>
/* -------------------------------------------------------------------------- */
#include <boost/preprocessor.hpp>
/* -------------------------------------------------------------------------- */
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
namespace py = pybind11;

namespace akantu {

/* -------------------------------------------------------------------------- */
#define PY_AKANTU_PP_VALUE(s, data, elem)                                      \
  .value(BOOST_PP_STRINGIZE(elem), BOOST_PP_CAT(data, elem))

#define PY_AKANTU_REGISTER_ENUM_(type_name, list, prefix, mod)                 \
  py::enum_<type_name>(mod, BOOST_PP_STRINGIZE(type_name))                     \
                                BOOST_PP_SEQ_FOR_EACH(PY_AKANTU_PP_VALUE,      \
                                                      prefix, list)            \
                                    .export_values()

#define PY_AKANTU_REGISTER_CLASS_ENUM(type_name, list, mod)                    \
  PY_AKANTU_REGISTER_ENUM_(type_name, list, type_name::_, mod)

#define PY_AKANTU_REGISTER_ENUM(type_name, list, mod)                          \
  PY_AKANTU_REGISTER_ENUM_(type_name, list, , mod)

/* -------------------------------------------------------------------------- */

std::vector<std::string> python_argv_str;
std::vector<char *> python_argv_char_ptr;
char ** python_argv = nullptr;
int python_argc;

void register_initialize(py::module & mod) {
  mod.def("__initialize", []() {
    auto && sys = py::module::import("sys");
    auto && argv = sys.attr("argv");
    python_argv_str = py::cast<std::vector<std::string>>(argv);
    for (auto && v : python_argv_str) {
      python_argv_char_ptr.push_back((char *)v.c_str());
    }
    python_argc = python_argv_str.size();
    python_argv = &python_argv_char_ptr[0];
    initialize(python_argc, python_argv);
  });
}

void register_enums(py::module & mod) {
  py::enum_<SpatialDirection>(mod, "SpatialDirection")
      .value("_x", _x)
      .value("_y", _y)
      .value("_z", _z)
      .export_values();

  py::enum_<AnalysisMethod>(mod, "AnalysisMethod")
      .value("_static", _static)
      .value("_implicit_dynamic", _implicit_dynamic)
      .value("_explicit_lumped_mass", _explicit_lumped_mass)
      .value("_explicit_lumped_capacity", _explicit_lumped_capacity)
      .value("_explicit_consistent_mass", _explicit_consistent_mass)
      .value("_explicit_contact", _explicit_contact)
      .value("_implicit_contact", _implicit_contact)
      .export_values();

  PY_AKANTU_REGISTER_CLASS_ENUM(ModelType, AKANTU_MODEL_TYPES, mod);
  PY_AKANTU_REGISTER_CLASS_ENUM(NonLinearSolverType,
                                AKANTU_NON_LINEAR_SOLVER_TYPES, mod);
  PY_AKANTU_REGISTER_CLASS_ENUM(TimeStepSolverType,
                                AKANTU_TIME_STEP_SOLVER_TYPE, mod);
  PY_AKANTU_REGISTER_CLASS_ENUM(IntegrationSchemeType,
                                AKANTU_INTEGRATION_SCHEME_TYPE, mod);
  PY_AKANTU_REGISTER_CLASS_ENUM(SolveConvergenceCriteria,
                                AKANTU_SOLVE_CONVERGENCE_CRITERIA, mod);

  py::enum_<CohesiveMethod>(mod, "CohesiveMethod")
      .value("_intrinsic", _intrinsic)
      .value("_extrinsic", _extrinsic)
      .export_values();

  py::enum_<GhostType>(mod, "GhostType")
      .value("_not_ghost", _not_ghost)
      .value("_ghost", _ghost)
      .value("_casper", _casper)
      .export_values();

  py::enum_<MeshIOType>(mod, "MeshIOType")
      .value("_miot_auto", _miot_auto)
      .value("_miot_gmsh", _miot_gmsh)
      .value("_miot_gmsh_struct", _miot_gmsh_struct)
      .value("_miot_diana", _miot_diana)
      .value("_miot_abaqus", _miot_abaqus)
      .export_values();

  py::enum_<MatrixType>(mod, "MatrixType")
      .value("_unsymmetric", _unsymmetric)
      .value("_symmetric", _symmetric)
      .export_values();

  py::module mod_rdt = mod.def_submodule("RandomDistributionType");
  py::enum_<RandomDistributionType>(mod_rdt, "RandomDistributionTypes")
      .value("_uniform", RandomDistributionType::_uniform)
      .value("_exponential", RandomDistributionType::_exponential)
      .value("_gamma", RandomDistributionType::_gamma)
      .value("_weibull", RandomDistributionType::_weibull)
      .value("_extreme_value", RandomDistributionType::_extreme_value)
      .value("_normal", RandomDistributionType::_normal)
      .value("_lognormal", RandomDistributionType::_lognormal)
      .value("_chi_squared", RandomDistributionType::_chi_squared)
      .value("_cauchy", RandomDistributionType::_cauchy)
      .value("_fisher_f", RandomDistributionType::_fisher_f)
      .value("_student_t", RandomDistributionType::_student_t)
      .export_values();

  PY_AKANTU_REGISTER_ENUM(ElementType, AKANTU_ALL_ELEMENT_TYPE(_not_defined),
                          mod);
  PY_AKANTU_REGISTER_ENUM(ElementKind, AKANTU_ELEMENT_KIND(_ek_not_defined),
                          mod);
}

/* -------------------------------------------------------------------------- */
#define AKANTU_PP_STR_TO_TYPE2(s, data, elem) ({BOOST_PP_STRINGIZE(elem), elem})

void register_functions(py::module & mod) {

  mod.def("getElementTypes", []() {
    std::map<std::string, akantu::ElementType> element_types{
        BOOST_PP_SEQ_FOR_EACH_I(
            AKANTU_PP_ENUM, BOOST_PP_SEQ_SIZE(AKANTU_ek_regular_ELEMENT_TYPE),
            BOOST_PP_SEQ_TRANSFORM(AKANTU_PP_STR_TO_TYPE2, akantu,
                                   AKANTU_ek_regular_ELEMENT_TYPE))};

    return element_types;
  });
}

#undef AKANTU_PP_STR_TO_TYPE2

} // namespace akantu
