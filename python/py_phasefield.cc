/* -------------------------------------------------------------------------- */
#include "py_aka_array.hh"
#include "py_constitutive_law.hh"
/* -------------------------------------------------------------------------- */
#include <constitutive_law.hh>
#include <phase_field_model.hh>
#include <phasefield_linear.hh>
#include <phasefield_quadratic.hh>
#include <phasefield_selector.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
namespace py = pybind11;
/* -------------------------------------------------------------------------- */

namespace akantu {

namespace {
template <typename _PhaseField> class PyPhaseField : public _PhaseField {
public:
  /* Inherit the constructors */
  using _PhaseField::_PhaseField;

  void initPhaseField() override {
    // NOLINTNEXTLINE
    PYBIND11_OVERRIDE(void, _PhaseField, initPhaseField, );
  };

  void computeDrivingForce(ElementType el_type,
                           GhostType ghost_type = _not_ghost) override {
    // NOLINTNEXTLINE
    PYBIND11_OVERRIDE(void, _PhaseField, computeDrivingForce, el_type,
                      ghost_type);
  }

  void computeDissipatedEnergy(ElementType el_type) override {
    // NOLINTNEXTLINE
    PYBIND11_OVERRIDE(void, _PhaseField, computeDissipatedEnergy, el_type);
  }
};

/* ------------------------------------------------------------------------ */
template <typename _PhaseField>
void register_phasefield_classes(py::module & mod, const std::string & name) {
  py::class_<_PhaseField, PhaseField, PyPhaseField<_PhaseField>>(
      mod, name.c_str(), py::multiple_inheritance())
      .def(py::init<PhaseFieldModel &, const ID &>());
}

} // namespace

/* -------------------------------------------------------------------------- */
void register_phasefield(py::module & mod) {
  register_constitutive_law<PhaseFieldModel>(mod);

  py::class_<PhaseField, PyPhaseField<PhaseField>,
             ConstitutiveLaw<PhaseFieldModel>>(mod, "PhaseField",
                                               py::multiple_inheritance())
      .def(py::init<PhaseFieldModel &, const ID &>())
      .def(
          "getDamage",
          [](PhaseField & self, ElementType el_type,
             GhostType ghost_type = _not_ghost) -> decltype(auto) {
            return self.getDamage(el_type, ghost_type);
          },
          py::arg("el_type"), py::arg("ghost_type") = _not_ghost,
          py::return_value_policy::reference)

      .def(
          "getDamage",
          [](PhaseField & self) -> decltype(auto) { return self.getDamage(); },
          py::return_value_policy::reference)
      .def(
          "getStrain",
          [](PhaseField & self, ElementType el_type,
             GhostType ghost_type = _not_ghost) -> decltype(auto) {
            return self.getStrain(el_type, ghost_type);
          },
          py::arg("el_type"), py::arg("ghost_type") = _not_ghost,
          py::return_value_policy::reference)

      .def(
          "getStrain",
          [](PhaseField & self) -> decltype(auto) { return self.getStrain(); },
          py::return_value_policy::reference)
      .def(
          "getEnergy",
          [](PhaseField & self, ID energy_id) -> Real {
            return self.getEnergy(energy_id);
          },
          py::arg("energy_id"))
      .def(
          "getEnergy",
          [](PhaseField & self, ID energy_id, Element element) -> Real {
            return self.getEnergy(energy_id, element);
          },
          py::arg("energy_id"), py::arg("element"));

  py::class_<PhaseFieldFactory>(mod, "PhaseFieldFactory")
      .def_static(
          "getInstance",
          []() -> PhaseFieldFactory & { return PhaseField::getFactory(); },
          py::return_value_policy::reference)
      .def("registerAllocator", [](PhaseFieldFactory & self,
                                   const std::string id, py::function func) {
        self.registerAllocator(
            id,
            [func, id](Int dim, const ID & energy_split,
                       PhaseFieldModel & model,
                       const ID & option) -> std::unique_ptr<PhaseField> {
              py::object obj = func(dim, energy_split, std::ref(model), id);
              auto & ptr = py::cast<PhaseField &>(obj);
              obj.release();
              return std::unique_ptr<PhaseField>(&ptr);
            });
      });
}

} // namespace akantu
