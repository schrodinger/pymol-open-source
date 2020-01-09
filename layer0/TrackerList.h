#pragma once

#include "Tracker.h"

namespace pymol
{

/**
 * Manages a Tracker list provided from outside function
 */

template <typename T> class TrackerAdapter
{
  /**
   * Iterator for pymol::ListAdapter
   * Not really  aiming to satisfy any of the standard
   * named requirements for iterators, but to only implement
   * ones that will likely be used.
   * Uses IsConst pattern to avoid having to define non-const and const iterator
   * versions
   */

  template <bool IsConst> class TrackerAdapterIterator
  {
  public:
    using value_type = typename std::conditional<IsConst, const T, T>::type;
    using pointer = value_type*;
    using reference = value_type&;

    TrackerAdapterIterator()
        : m_ptr(nullptr)
    {
    }

    explicit TrackerAdapterIterator(CTracker* tracker, int list_id)
        : m_tracker(tracker)
        , m_iter_id(TrackerNewIter(tracker, 0, list_id))
    {
      ++*this;
    }

    ~TrackerAdapterIterator()
    {
      // safe to call for m_iter_id == -1
      TrackerDelIter(m_tracker, m_iter_id);
    }

    // move only because it owns m_iter_id
    TrackerAdapterIterator(TrackerAdapterIterator&& other)
        : m_tracker(other.m_tracker)
        , m_ptr(other.m_ptr)
    {
      std::swap(m_iter_id, other.m_iter_id);
    }

    reference operator*() { return *m_ptr; }

    pointer operator->() { return m_ptr; }

    TrackerAdapterIterator& operator++()
    {
      m_ptr = nullptr;
      TrackerIterNextCandInList(
          m_tracker, m_iter_id, (TrackerRef**) (void*) &m_ptr);
      return *this;
    }

    bool operator==(const TrackerAdapterIterator& iter) const
    {
      return m_ptr == iter.m_ptr;
    }

    bool operator!=(const TrackerAdapterIterator& iter) const
    {
      return !(m_ptr == iter.m_ptr);
    }

  private:
    int m_iter_id = -1;
    CTracker* m_tracker;
    pointer m_ptr;
  };

public:
  using type = T;
  using iterator = TrackerAdapterIterator<false>;
  using const_iterator = TrackerAdapterIterator<true>;

  // Takes ownership of tracker id: list_id
  explicit TrackerAdapter(CTracker* tracker, int list_id)
      : m_tracker(tracker)
      , m_list_id(list_id)
  {
  }
  ~TrackerAdapter()
  {
    TrackerDelList(m_tracker, m_list_id);
  }
  iterator begin() { return iterator(m_tracker, m_list_id); }
  const_iterator begin() const { return const_iterator(m_tracker, m_list_id); }
  iterator end() { return iterator(); }
  const_iterator end() const { return const_iterator(); }

private:
  CTracker* m_tracker;
  int m_list_id;
};

} // namespace pymol
