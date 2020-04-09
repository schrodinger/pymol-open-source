#pragma once
#include <cassert>
#include <iostream>
#include <iterator>

namespace pymol
{

/**
 * Non-owning view into PyMOL's singly linked list.
 */

template <typename T> class ListAdapter
{
  /**
   * Iterator for pymol::ListAdapter
   * Not really  aiming to satisfy any of the standard
   * named requirements for iterators, but to only implement
   * ones that will likely be used.
   * Uses IsConst pattern to avoid having to define non-const and const iterator
   * versions
   */

  template <bool IsConst> class ListAdapterIterator
  {
  public:
    using iterator_category = std::forward_iterator_tag;
    using difference_type = std::ptrdiff_t;
    using value_type = typename std::conditional<IsConst, const T, T>::type;
    using pointer = value_type*;
    using reference = value_type&;

    explicit ListAdapterIterator(pointer ptr)
        : m_ptr(ptr)
    {
    }

    /**
     * Creates a const_iterator from iterator but not the other way around
     * @param other non const iterator
     */
    template <bool InnerIsConst = IsConst,
        typename = typename std::enable_if<InnerIsConst>>
    ListAdapterIterator(const ListAdapterIterator<false>& other)
        : m_ptr(other.m_ptr)
    {
    }

    reference operator*() { return *m_ptr; }

    pointer operator->() { return m_ptr; }

    ListAdapterIterator& operator++()
    {
      assert(m_ptr);
      m_ptr = m_ptr->next;
      return *this;
    }

    ListAdapterIterator operator++(int)
    {
      ListAdapterIterator tmp(*this);
      operator++();
      return tmp;
    }

    bool operator==(const ListAdapterIterator& iter) const
    {
      return m_ptr == iter.m_ptr;
    }

    bool operator!=(const ListAdapterIterator& iter) const
    {
      return !(m_ptr == iter.m_ptr);
    }

    friend class ListAdapterIterator<true>;

  private:
    pointer m_ptr;
  };

public:
  using type = T;
  using iterator = ListAdapterIterator<false>;
  using const_iterator = ListAdapterIterator<true>;
  explicit ListAdapter(T* list)
      : m_list(list)
  {
  }
  iterator begin() noexcept { return iterator(m_list); }
  iterator end() noexcept { return iterator(nullptr); }
  const_iterator begin() const noexcept { return const_iterator(m_list); }
  const_iterator end() const noexcept { return const_iterator(nullptr); }
  const_iterator cbegin() const noexcept { return const_iterator(m_list); }
  const_iterator cend() const noexcept { return const_iterator(nullptr); }

private:
  type* m_list;
};

template <typename T> ListAdapter<T> make_list_adapter(T *list) {
  return ListAdapter<T>(list);
}

} // namespace pymol
