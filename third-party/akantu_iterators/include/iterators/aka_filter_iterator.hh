/**
 * @file   aka_filter_iterator.hh
 *
 * @author Nicolas Richart
 *
 * @date creation  jeu déc 12 2019
 *
 * @brief A Documented file.
 *
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */
/* -------------------------------------------------------------------------- */
#include <iterator>
#include <utility>
/* -------------------------------------------------------------------------- */

#ifndef AKA_FILTER_ITERATOR_HH
#define AKA_FILTER_ITERATOR_HH

#ifndef AKANTU_ITERATORS_NAMESPACE
#define AKANTU_ITERATORS_NAMESPACE akantu
#endif

namespace AKANTU_ITERATORS_NAMESPACE {

/* -------------------------------------------------------------------------- */
namespace iterators AKA_ITERATOR_EXPORT_NAMESPACE {
  template <class filter_iterator_t, class container_iterator_t>
  class FilterIterator
      : public details::CopyAssignmentEnabler<aka::conjunction<
            std::is_copy_assignable<filter_iterator_t>,
            std::is_copy_constructible<filter_iterator_t>,
            std::is_copy_assignable<container_iterator_t>,
            std::is_copy_constructible<container_iterator_t>>::value>,
        public details::MoveAssignmentEnabler<aka::conjunction<
            std::is_move_assignable<filter_iterator_t>,
            std::is_move_constructible<filter_iterator_t>,
            std::is_move_assignable<container_iterator_t>,
            std::is_move_constructible<container_iterator_t>>::value> {
  public:
    using value_type = typename std::decay_t<container_iterator_t>::value_type;
    using difference_type =
        typename std::decay_t<filter_iterator_t>::difference_type;
    using pointer = std::decay_t<value_type> *;
    using reference = typename std::decay_t<container_iterator_t>::reference;
    using iterator_category =
        typename std::decay_t<filter_iterator_t>::iterator_category;

    FilterIterator(filter_iterator_t filter_it,
                   const container_iterator_t & container_begin)
        : filter_it(std::move(filter_it)), container_begin(container_begin),
          container_it(container_begin) {}

    auto operator++() -> FilterIterator & {
      ++filter_it;
      return *this;
    }

    auto operator++(int) -> FilterIterator {
      auto cpy = *this;
      this->operator++();
      return cpy;
    }

    auto operator*() -> reference {
      container_it = this->container_begin + *this->filter_it;
      return *container_it;
    }

    auto operator!=(const FilterIterator & other) const -> bool {
      return (filter_it != other.filter_it) or
             (container_begin != other.container_begin);
    }

    template <class iterator_category_ = iterator_category,
              std::enable_if_t<aka::is_iterator_category_at_least<
                  iterator_category_,
                  std::bidirectional_iterator_tag>::value> * = nullptr>
    auto operator--() -> FilterIterator & {
      --filter_it;
      return *this;
    }

    template <class iterator_category_ = iterator_category,
              std::enable_if_t<aka::is_iterator_category_at_least<
                  iterator_category_,
                  std::bidirectional_iterator_tag>::value> * = nullptr>
    auto operator--(int) -> FilterIterator {
      auto cpy = *this;
      this->operator--();
      return cpy;
    }

    template <class iterator_category_ = iterator_category,
              std::enable_if_t<aka::is_iterator_category_at_least<
                  iterator_category_,
                  std::random_access_iterator_tag>::value> * = nullptr>
    auto operator-(const FilterIterator & other) -> difference_type {
      return filter_it - other.index;
    }

    // random iterator it[idx]
    template <class iterator_category_ = iterator_category,
              std::enable_if_t<aka::is_iterator_category_at_least<
                  iterator_category_,
                  std::random_access_iterator_tag>::value> * = nullptr>
    auto operator[](difference_type idx) -> reference {
      container_it = this->container_begin + *(this->filter_it + idx);
      return *container_it;
    }

    // random iterator it + n
    template <class iterator_category_ = iterator_category,
              std::enable_if_t<aka::is_iterator_category_at_least<
                  iterator_category_,
                  std::random_access_iterator_tag>::value> * = nullptr>
    auto operator+(difference_type n) -> FilterIterator {
      auto it = FilterIterator(filter_it + n, container_begin);
      return it;
    }

    // random iterator it - n
    template <class iterator_category_ = iterator_category,
              std::enable_if_t<aka::is_iterator_category_at_least<
                  iterator_category_,
                  std::random_access_iterator_tag>::value> * = nullptr>
    auto operator-(difference_type n) -> FilterIterator {
      auto it = FilterIterator(filter_it - n, container_begin);
      return it;
    }

    template <
        class iterator_category_ = iterator_category,
        std::enable_if_t<aka::is_iterator_category_at_least<
            iterator_category_, std::forward_iterator_tag>::value> * = nullptr>
    auto operator==(const FilterIterator & other) const -> bool {
      return (filter_it == other.filter_it) and
             (container_begin == other.container_begin);
    }

  private:
    filter_iterator_t filter_it;
    container_iterator_t container_begin;
    container_iterator_t container_it;
  };

  template <class filter_iterator_t, class container_iterator_t>
  auto make_filter_iterator(filter_iterator_t && filter_it,
                            container_iterator_t && container_begin)
      -> decltype(auto) {
    return FilterIterator<filter_iterator_t, container_iterator_t>(
        std::forward<filter_iterator_t>(filter_it),
        std::forward<container_iterator_t>(container_begin));
  }

  template <class container_iterator_t, class Predicate>
  class FilterIfIterator
      : public details::CopyAssignmentEnabler<aka::conjunction<
            std::is_copy_assignable<Predicate>,
            std::is_copy_constructible<Predicate>,
            std::is_copy_assignable<container_iterator_t>,
            std::is_copy_constructible<container_iterator_t>>::value>,
        public details::MoveAssignmentEnabler<aka::conjunction<
            std::is_move_assignable<Predicate>,
            std::is_move_constructible<Predicate>,
            std::is_move_assignable<container_iterator_t>,
            std::is_move_constructible<container_iterator_t>>::value> {
  public:
    using value_type = typename std::decay_t<container_iterator_t>::value_type;
    using difference_type =
        typename std::decay_t<container_iterator_t>::difference_type;
    using pointer = std::decay_t<value_type> *;
    using reference = typename std::decay_t<container_iterator_t>::reference;
    using iterator_category = std::common_type_t<
        std::forward_iterator_tag,
        typename std::decay_t<container_iterator_t>::iterator_category>;

    FilterIfIterator(container_iterator_t && container_begin,
                     container_iterator_t && container_end,
                     Predicate && predicate)
        : container_it(container_begin), container_end(container_end),
          predicate(predicate) {
      if (container_it != container_end and not predicate(*container_it)) {
        this->operator++();
      }
    }

    auto operator++() -> FilterIfIterator & {
      do {
        ++container_it;
      } while (container_it != container_end and not predicate(*container_it));
      return *this;
    }

    auto operator++(int) -> FilterIfIterator {
      auto cpy = *this;
      this->operator++();
      return cpy;
    }

    auto operator*() -> decltype(auto) { return (*container_it); }

    auto operator!=(const FilterIfIterator & other) const -> bool {
      return (container_it != other.container_it);
    }

    template <
        class iterator_category_ = iterator_category,
        std::enable_if_t<aka::is_iterator_category_at_least<
            iterator_category_, std::forward_iterator_tag>::value> * = nullptr>
    auto operator==(const FilterIfIterator & other) const -> bool {
      return (container_it == other.container_it);
    }

  private:
    container_iterator_t container_it, container_end;
    Predicate predicate;
  };

  template <class container_iterator_t, class Predicate>
  auto make_filter_if_iterator(container_iterator_t && container_begin,
                               container_iterator_t && container_end,
                               Predicate && predicate) -> decltype(auto) {
    return FilterIfIterator<container_iterator_t, Predicate>(
        std::forward<container_iterator_t>(container_begin),
        std::forward<container_iterator_t>(container_end),
        std::forward<Predicate>(predicate));
  }
} // namespace iterators AKA_ITERATOR_EXPORT_NAMESPACE

namespace containers AKA_ITERATOR_EXPORT_NAMESPACE {
  template <class filter_t, class Container> class FilterAdaptor {
  public:
    using size_type = typename std::decay_t<Container>::size_type;

    FilterAdaptor(filter_t && filter, Container && container)
        : filter(std::forward<filter_t>(filter)),
          container(std::forward<Container>(container)) {
      static_assert(
          std::is_same<typename decltype(container.begin())::iterator_category,
                       std::random_access_iterator_tag>::value,
          "Containers must all have random iterators");
    }

    auto begin() const -> decltype(auto) {
      return iterators::make_filter_iterator(filter.begin(), container.begin());
    }

    auto begin() -> decltype(auto) {
      return iterators::make_filter_iterator(filter.begin(), container.begin());
    }

    auto end() const -> decltype(auto) {
      return iterators::make_filter_iterator(filter.end(), container.begin());
    }

    auto end() -> decltype(auto) {
      return iterators::make_filter_iterator(filter.end(), container.begin());
    }

  private:
    filter_t filter;
    Container container;
  };

  template <class Container, class Predicate> class FilterIfAdaptor {
  public:
    using size_type = typename std::decay_t<Container>::size_type;

    FilterIfAdaptor(Container && container, Predicate && predicate)
        : container(container), predicate(predicate) {}

    auto begin() const -> decltype(auto) {
      return iterators::make_filter_if_iterator(container.begin(),
                                                container.end(), predicate);
    }

    auto end() const -> decltype(auto) {
      return iterators::make_filter_if_iterator(container.end(),
                                                container.end(), predicate);
    }

    auto begin() -> decltype(auto) {
      return iterators::make_filter_if_iterator(container.begin(),
                                                container.end(), predicate);
    }

    auto end() -> decltype(auto) {
      return iterators::make_filter_if_iterator(container.end(),
                                                container.end(), predicate);
    }

  private:
    Container container;
    Predicate predicate;
  };

} // namespace containers AKA_ITERATOR_EXPORT_NAMESPACE

template <class filter_t, class Container>
auto filter(filter_t && filter, Container && container) -> decltype(auto) {
  return containers::FilterAdaptor<filter_t, Container>(
      std::forward<filter_t>(filter), std::forward<Container>(container));
}

template <class Container, class Predicate>
auto filter_if(Container && container, Predicate && predicate)
    -> decltype(auto) {
  return containers::FilterIfAdaptor<Container, Predicate>(
      std::forward<Container>(container), std::forward<Predicate>(predicate));
}

} // namespace AKANTU_ITERATORS_NAMESPACE

namespace std {
template <class filter_iterator_t, class container_iterator_t>
struct iterator_traits<::AKANTU_ITERATORS_NAMESPACE::iterators::FilterIterator<
    filter_iterator_t, container_iterator_t>> {
private:
  using iterator_type =
      typename ::AKANTU_ITERATORS_NAMESPACE::iterators::FilterIterator<
          filter_iterator_t, container_iterator_t>;

public:
  using iterator_category = typename iterator_type::iterator_category;
  using value_type = typename iterator_type::value_type;
  using difference_type = typename iterator_type::difference_type;
  using pointer = typename iterator_type::pointer;
  using reference = typename iterator_type::reference;
};
} // namespace std

#endif // AKA_FILTER_ITERATOR_HH
