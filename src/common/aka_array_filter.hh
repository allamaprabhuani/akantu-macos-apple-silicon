/**
 * Copyright (©) 2014-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
// IWYU pragma: private, include "aka_array.hh"

/* -------------------------------------------------------------------------- */
#include "aka_array.hh"
#include "aka_common.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_ARRAY_FILTER_HH
#define AKANTU_ARRAY_FILTER_HH

namespace akantu {

/* ------------------------------------------------------------------------ */
/* ArrayFilter                                                              */
/* ------------------------------------------------------------------------ */
template <typename T> class ArrayFilter {
public:
  /// const iterator
  template <class original_iterator, typename filter_iterator>
  class const_iterator {
  public:
    auto getCurrentIndex() -> Idx {
      return ((*this->filter_it * this->nb_item_per_elem) +
              this->sub_element_counter);
    }

    const_iterator() = default;
    const_iterator(original_iterator origin_it, filter_iterator filter_it,
                   Int nb_item_per_elem)
        : origin_it(std::move(origin_it)), filter_it(std::move(filter_it)),
          nb_item_per_elem(nb_item_per_elem), sub_element_counter(0) {};

    auto operator!=(const_iterator & other) const -> bool {
      return !((*this) == other);
    }
    auto operator==(const_iterator & other) const -> bool {
      return (this->origin_it == other.origin_it &&
              this->filter_it == other.filter_it &&
              this->sub_element_counter == other.sub_element_counter);
    }

    auto operator!=(const const_iterator & other) const -> bool {
      return !((*this) == other);
    }
    auto operator==(const const_iterator & other) const -> bool {
      return (this->origin_it == other.origin_it &&
              this->filter_it == other.filter_it &&
              this->sub_element_counter == other.sub_element_counter);
    }

    auto operator++() -> const_iterator & {
      ++sub_element_counter;
      if (sub_element_counter == nb_item_per_elem) {
        sub_element_counter = 0;
        ++filter_it;
      }
      return *this;
    };

    auto operator*() -> decltype(auto) {
      return origin_it[(nb_item_per_elem * (*filter_it)) + sub_element_counter];
    };

  private:
    original_iterator origin_it;
    filter_iterator filter_it;
    /// the number of item per element
    Int nb_item_per_elem;
    /// counter for every sub element group
    Int sub_element_counter;
  };

  using vector_iterator = void; // iterator<Vector<T>>;
  using array_type = Array<T>;
  using const_vector_iterator =
      const_iterator<typename array_type::const_vector_iterator,
                     Array<Idx>::const_scalar_iterator>;
  using value_type = typename array_type::value_type;

private:
  /* ---------------------------------------------------------------------- */
  /* Constructors/Destructors                                               */
  /* ---------------------------------------------------------------------- */
public:
  ArrayFilter(const Array<T> & array, const Array<Idx> & filter,
              Int nb_item_per_elem)
      : array(array), filter(filter), nb_item_per_elem(nb_item_per_elem) {};

  auto begin_reinterpret(Int n, Int new_size) const -> decltype(auto) {
    Int new_nb_item_per_elem = this->nb_item_per_elem;
    if (new_size != 0 && n != 0) {
      new_nb_item_per_elem = this->array.getNbComponent() *
                             this->filter.size() * this->nb_item_per_elem /
                             (n * new_size);
    }

    return const_vector_iterator(make_view(this->array, n).begin(),
                                 this->filter.begin(), new_nb_item_per_elem);
  };

  auto end_reinterpret(Int n, Int new_size) const -> decltype(auto) {
    Int new_nb_item_per_elem = this->nb_item_per_elem;
    if (new_size != 0 && n != 0) {
      new_nb_item_per_elem = this->array.getNbComponent() *
                             this->filter.size() * this->nb_item_per_elem /
                             (n * new_size);
    }

    return const_vector_iterator(make_view(this->array, n).begin(),
                                 this->filter.end(), new_nb_item_per_elem);
  };

  /// return the size of the filtered array which is the filter size
  [[nodiscard]] auto size() const -> Int {
    return this->filter.size() * this->nb_item_per_elem;
  };
  /// the number of components of the filtered array
  [[nodiscard]] auto getNbComponent() const -> Int {
    return this->array.getNbComponent();
  };

  /// tells if the container is empty
  [[nodiscard]] auto empty() const -> bool { return (size() == 0); }
  /* ---------------------------------------------------------------------- */
  /* Class Members                                                          */
  /* ---------------------------------------------------------------------- */
private:
  /// reference to array of data
  const Array<T> & array;
  /// reference to the filter used to select elements
  const Array<Idx> & filter;
  /// the number of item per element
  Int nb_item_per_elem;
};

} // namespace akantu

#endif // AKANTU_ARRAY_FILTER_HH
