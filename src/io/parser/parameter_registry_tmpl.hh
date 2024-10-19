/**
 * Copyright (©) 2016-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_error.hh"
#include "aka_iterators.hh"
#include "aka_types.hh"
#include "parser.hh"
/* -------------------------------------------------------------------------- */
#include <algorithm>
#include <set>
#include <string>
#include <vector>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_PARAMETER_REGISTRY_TMPL_HH_
#define AKANTU_PARAMETER_REGISTRY_TMPL_HH_

namespace akantu {

namespace debug {
  class ParameterException : public Exception {
  public:
    ParameterException(const std::string & name, const std::string & message)
        : Exception(message), name(name) {}
    const std::string & name;
  };

  class ParameterUnexistingException : public ParameterException {
  public:
    ParameterUnexistingException(const std::string & name,
                                 const ParameterRegistry & registry)
        : ParameterException(name, "Parameter " + name +
                                       " does not exists in this scope") {
      auto && params = registry.listParameters();
      this->_info =
          std::accumulate(params.begin(), params.end(),
                          this->_info + "\n Possible parameters are: ",
                          [](auto && str, auto && param) {
                            static auto first = true;
                            auto ret = str + (first ? " " : ", ") + param;
                            first = false;
                            return ret;
                          });
    }
  };

  class ParameterAccessRightException : public ParameterException {
  public:
    ParameterAccessRightException(const std::string & name,
                                  const std::string & perm)
        : ParameterException(name, "Parameter " + name + " is not " + perm) {}
  };

  class ParameterWrongTypeException : public ParameterException {
  public:
    ParameterWrongTypeException(const std::string & name,
                                const std::type_info & wrong_type,
                                const std::type_info & type)
        : ParameterException(name, "Parameter " + name +
                                       " type error, cannot convert " +
                                       debug::demangle(type.name()) + " to " +
                                       debug::demangle(wrong_type.name())) {}
  };
} // namespace debug
/* -------------------------------------------------------------------------- */
template <typename T>
const ParameterTyped<T> & Parameter::getParameterTyped() const {
  try {
    const auto & tmp = aka::as_type<ParameterTyped<T>>(*this);
    return tmp;
  } catch (std::bad_cast &) {
    AKANTU_CUSTOM_EXCEPTION(
        debug::ParameterWrongTypeException(name, typeid(T), this->type()));
  }
}

/* -------------------------------------------------------------------------- */
template <typename T> ParameterTyped<T> & Parameter::getParameterTyped() {
  try {
    auto & tmp = aka::as_type<ParameterTyped<T>>(*this);
    return tmp;
  } catch (std::bad_cast &) {
    AKANTU_CUSTOM_EXCEPTION(
        debug::ParameterWrongTypeException(name, typeid(T), this->type()));
  }
}

/* ------------------------------------------------------------------------ */
inline void Parameter::setAuto(const ParserParameter & /*value*/) {
  if (not isParsable()) {
    AKANTU_CUSTOM_EXCEPTION(
        debug::ParameterAccessRightException(name, "parsable"));
  }
}

/* -------------------------------------------------------------------------- */
template <typename T> const T & Parameter::get() const {
  if (not isReadable()) {
    AKANTU_CUSTOM_EXCEPTION(
        debug::ParameterAccessRightException(name, "readable"));
  }
  const ParameterTyped<T> & typed_param = getParameterTyped<T>();
  return typed_param.getTyped();
}

/* -------------------------------------------------------------------------- */
template <typename T> T & Parameter::get() {
  ParameterTyped<T> & typed_param = getParameterTyped<T>();
  if (not isReadable() or not this->isWritable()) {
    AKANTU_CUSTOM_EXCEPTION(
        debug::ParameterAccessRightException(name, "accessible"));
  }
  return typed_param.getTyped();
}

/* -------------------------------------------------------------------------- */
template <typename T> inline Parameter::operator T() const {
  return this->get<T>();
}

/* -------------------------------------------------------------------------- */
template <typename T>
ParameterTyped<T>::ParameterTyped(const std::string & name,
                                  const std::string & description,
                                  ParameterAccessType param_type, T & param)
    : Parameter(name, description, param_type), param(param) {}

/* -------------------------------------------------------------------------- */
template <typename T> inline void ParameterTyped<T>::set(std::any value) {
  try {
    param = std::any_cast<T>(value);
  } catch (std::bad_any_cast & e) {
    AKANTU_EXCEPTION("Cannot cast std::any to parameter\n" << e.what());
  }
}

/* -------------------------------------------------------------------------- */
template <typename T> inline std::string ParameterTyped<T>::to_string() {
  if constexpr (std::is_same_v<std::string, T>)
    return param;
  else {
    return std::to_string(param);
  }
}
/* -------------------------------------------------------------------------- */

template <typename T>
inline void ParameterTyped<T>::setAuto(const ParserParameter & value) {
  Parameter::setAuto(value);
  param = static_cast<T>(value);
}

/* -------------------------------------------------------------------------- */
template <>
inline void
ParameterTyped<std::string>::setAuto(const ParserParameter & value) {
  Parameter::setAuto(value);
  param = value.getValue();
}

/* -------------------------------------------------------------------------- */
template <typename T> const T & ParameterTyped<T>::getTyped() const {
  return param;
}
/* -------------------------------------------------------------------------- */
template <typename T> T & ParameterTyped<T>::getTyped() { return param; }

/* -------------------------------------------------------------------------- */
template <typename T>
inline void ParameterTyped<T>::printself(std::ostream & stream) const {
  Parameter::printself(stream);
  stream << param << "\n";
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
template <typename T, Int m, Int n>
class ParameterTyped<Eigen::Matrix<T, m, n>> : public Parameter {
public:
  ParameterTyped(const std::string & name, const std::string & description,
                 ParameterAccessType param_type, Eigen::Matrix<T, m, n> & param)
      : Parameter(name, description, param_type), param(param) {}

  std::string to_string() override {
    std::stringstream sstr;
    sstr << param;
    return sstr.str();
  }

  void set(std::any value) {
    param = std::any_cast<Eigen::Matrix<T, m, n>>(value);
  }
  void setAuto(const ParserParameter & value) override {
    Parameter::setAuto(value);
    if constexpr (n == 1) {
      Vector<Real> tmp = value;
      if (param.rows() != 0) {
        param = tmp.block(0, 0, param.rows(), 1);
      } else {
        param = tmp;
      }
    } else {
      Matrix<Real> tmp = value;
      if (param.rows() != 0 and param.cols() != 0) {
        param = tmp.block(0, 0, param.rows(), param.cols());
      } else {
        param = tmp;
      }
    }
  }

  Eigen::Matrix<T, m, n> & getTyped() { return param; }
  const Eigen::Matrix<T, m, n> & getTyped() const { return param; }

  void printself(std::ostream & stream) const override {
    Parameter::printself(stream);

    stream << "[ ";
    for (int i : arange(param.rows())) {
      if (param.cols() != 1) {
        stream << "[ ";
      }
      for (int j : arange(param.cols())) {
        stream << param(i, j) << " ";
      }
      if (param.cols() != 1) {
        stream << "] ";
      }
    }
    stream << "]\n";
  }

  inline const std::type_info & type() const override {
    return typeid(Eigen::Matrix<T, m, n>);
  }

private:
  /// Value of parameter
  Eigen::Matrix<T, m, n> & param;
};

/* -------------------------------------------------------------------------- */
template <typename T> class ParameterTyped<std::vector<T>> : public Parameter {
public:
  ParameterTyped(const std::string & name, const std::string & description,
                 ParameterAccessType param_type, std::vector<T> & param)
      : Parameter(name, description, param_type), param(param) {}

  /* ------------------------------------------------------------------------ */

  std::string to_string() override {
    std::stringstream sstr;
    sstr << param;
    return sstr.str();
  }

  void set(std::any value) { param = value; }
  void setAuto(const ParserParameter & value) override {
    Parameter::setAuto(value);
    param.zero();
    const std::vector<T> & tmp = value;
    for (auto && z : tmp) {
      param.emplace_back(z);
    }
  }

  std::vector<T> & getTyped() { return param; }
  const std::vector<T> & getTyped() const { return param; }

  void printself(std::ostream & stream) const override {
    Parameter::printself(stream);
    stream << "[ ";
    for (auto && v : param) {
      stream << v << " ";
    }
    stream << "]\n";
  }

  inline const std::type_info & type() const override {
    return typeid(std::vector<T>);
  }

private:
  /// Value of parameter
  std::vector<T> & param;
};

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
template <typename T> class ParameterTyped<std::set<T>> : public Parameter {
public:
  ParameterTyped(const std::string & name, const std::string & description,
                 ParameterAccessType param_type, std::set<T> & param)
      : Parameter(name, description, param_type), param(param) {}

  /* ------------------------------------------------------------------------ */
  std::string to_string() override {
    std::stringstream sstr;
    sstr << "set[";
    for (auto && el : param) {
      sstr << el << ",";
    }
    sstr << "]";
    return sstr.str();
  }

  void set(std::any value) { param = std::any_cast<std::set<T>>(value); }
  void setAuto(const ParserParameter & value) {
    Parameter::setAuto(value);
    param.clear();
    const std::set<T> & tmp = value;
    for (auto && z : tmp) {
      param.emplace(z);
    }
  }

  std::set<T> & getTyped() { return param; }
  const std::set<T> & getTyped() const { return param; }

  void printself(std::ostream & stream) const {
    Parameter::printself(stream);
    stream << "[ ";
    for (auto && v : param) {
      stream << v << " ";
    }
    stream << "]\n";
  }

  inline const std::type_info & type() const { return typeid(std::set<T>); }

private:
  /// Value of parameter
  std::set<T> & param;
};

/* --------------------------------------------------------------------------
 */
template <>
inline void ParameterTyped<bool>::printself(std::ostream & stream) const {
  Parameter::printself(stream);
  stream << std::boolalpha << param << "\n";
}

/* --------------------------------------------------------------------------
 */
template <typename T>
void ParameterRegistry::registerParam(const std::string & name, T & variable,
                                      ParameterAccessType type,
                                      const std::string & description) {
  auto it = params.find(name);
  if (it != params.end()) {
    AKANTU_CUSTOM_EXCEPTION(debug::ParameterException(
        name, "Parameter named " + name + " already registered."));
  }
  auto && param =
      std::make_shared<ParameterTyped<T>>(name, description, type, variable);
  params[name] = std::move(param);
}

/* --------------------------------------------------------------------------
 */
template <typename T>
void ParameterRegistry::registerParam(const std::string & name, T & variable,
                                      const T & default_value,
                                      ParameterAccessType type,
                                      const std::string & description) {
  variable = default_value;
  registerParam(name, variable, type, description);
}

/* --------------------------------------------------------------------------
 */
template <typename T> T & ParameterRegistry::get_(const std::string & name) {
  auto it = params.find(name);
  if (it == params.end()) {
    if (consisder_sub) {
      for (auto && [_, registery] : sub_registries) {
        try {
          return registery.get().template get_<T>(name);
        } catch (...) {
        }
      }
    }

    // nothing was found not even in sub registries
    AKANTU_CUSTOM_EXCEPTION(debug::ParameterUnexistingException(name, *this));
  }

  Parameter & param = *(it->second);
  return param.get<T>();
}

/* --------------------------------------------------------------------------
 */
const Parameter & ParameterRegistry::get(const std::string & name) const {
  auto it = params.find(name);
  if (it == params.end()) {
    if (consisder_sub) {
      for (auto && [_, registery] : sub_registries) {
        try {
          return registery.get().get(name);
        } catch (...) {
        }
      }
    }

    // nothing was found not even in sub registries
    AKANTU_CUSTOM_EXCEPTION(debug::ParameterUnexistingException(name, *this));
  }

  Parameter & param = *(it->second);
  return param;
}

/* --------------------------------------------------------------------------
 */
Parameter & ParameterRegistry::get(const std::string & name) {
  auto it = params.find(name);
  if (it == params.end()) {
    if (consisder_sub) {
      for (auto && [_, registery] : sub_registries) {
        try {
          return registery.get().get(name);
        } catch (...) {
        }
      }
    }

    // nothing was found not even in sub registries
    AKANTU_CUSTOM_EXCEPTION(debug::ParameterUnexistingException(name, *this));
  }

  Parameter & param = *(it->second);
  return param;
}

/* --------------------------------------------------------------------------
 */
namespace details {
  template <class T, class R, class Enable = void> struct CastHelper {
    static R convert(const T & /*unused*/) { throw std::bad_cast(); }
  };

  template <class T, class R>
  struct CastHelper<T, R, std::enable_if_t<std::is_convertible<T, R>::value>> {
    static R convert(const T & val) { return val; }
  };
} // namespace details

template <typename T> inline ParameterTyped<T>::operator Real() const {
  if (not isReadable()) {
    AKANTU_CUSTOM_EXCEPTION(
        debug::ParameterAccessRightException(name, "accessible"));
  }
  return details::CastHelper<T, Real>::convert(param);
}

} // namespace akantu

#endif /* AKANTU_PARAMETER_REGISTRY_TMPL_HH_ */
