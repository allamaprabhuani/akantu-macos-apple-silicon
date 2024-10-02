#include "py_aka_array.hh"
#include "py_akantu_pybind11_compatibility.hh"
/* -------------------------------------------------------------------------- */
#include <constitutive_law.hh>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_PY_CONSTITUTIVE_LAW_HH_
#define AKANTU_PY_CONSTITUTIVE_LAW_HH_

namespace akantu {

void register_constitutive_law_internal_handler(pybind11::module & mod);

template <class ConstitutiveLawHandler>
class PyConstitutiveLaw : public ConstitutiveLaw<ConstitutiveLawHandler> {
public:
  using Parent = ConstitutiveLaw<ConstitutiveLawHandler>;
  /* Inherit the constructors */
  using Parent::Parent;

  ~PyConstitutiveLaw() override = default;

  void initConstitutiveLaw() override {
    // NOLINTNEXTLINE
    PYBIND11_OVERRIDE(void, Parent, initConstitutiveLaw, );
  }

  Real getEnergy(const ID & energy_id) override {
    // NOLINTNEXTLINE
    PYBIND11_OVERRIDE(Real, Parent, getEnergy, energy_id);
  }
  Real getEnergy(const ID & energy_id, const Element & element) override {
    // NOLINTNEXTLINE
    PYBIND11_OVERRIDE(Real, Parent, getEnergy, energy_id, element);
  }

  // methods need to be defined to be able to do the python interface
  Int getNbData(const Array<Element> & /*elements*/,
                const SynchronizationTag & /*tag*/) const override {
    return 0;
  }
  void packData(CommunicationBuffer & /*buffer*/,
                const Array<Element> & /*element*/,
                const SynchronizationTag & /*tag*/) const override {}
  void unpackData(CommunicationBuffer & /*buffer*/,
                  const Array<Element> & /*element*/,
                  const SynchronizationTag & /*tag*/) override{};
};

/* ------------------------------------------------------------------------ */
template <class ConstitutiveLawsHandler_>
void register_constitutive_law(py::module & mod) {
  using CL = ConstitutiveLaw<ConstitutiveLawsHandler_>;
  const std::string & name =
      "ConstitutiveLaw" + debug::demangle<ConstitutiveLawsHandler_>();

  py::class_<CL, ConstitutiveLawInternalHandler, Parsable,
             PyConstitutiveLaw<ConstitutiveLawsHandler_>>(
      mod, name.c_str(), py::multiple_inheritance())
      .def(py::template init<ConstitutiveLawsHandler_ &, const ID &, Int,
                             ElementKind, const ID &>())
      /*
       * These functions override the `Parsable` interface.
       * This ensure that the `updateInternalParameters()` function is called.
       */
      .def(
          "setReal",
          [](CL & self, const ID & name, const Real value) -> void {
            self.set(name, value);
            return;
          },
          py::arg("name"), py::arg("value"))
      .def(
          "setBool",
          [](CL & self, const ID & name, const bool value) -> void {
            self.set(name, value);
            return;
          },
          py::arg("name"), py::arg("value"))
      .def(
          "setString",
          [](CL & self, const ID & name, const std::string & value) -> void {
            self.set(name, value);
            return;
          },
          py::arg("name"), py::arg("value"))
      .def(
          "setInt",
          [](CL & self, const ID & name, const int value) -> void {
            self.set(name, value);
            return;
          },
          py::arg("name"), py::arg("value"));
}

} // namespace akantu

#endif // AKANTU_PY_CONSTITUTIVE_LAW_HH_
