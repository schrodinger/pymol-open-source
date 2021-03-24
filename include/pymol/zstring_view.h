/**
 * pymol::zstring_view
 * pymol::null_safe_zstring_view
 *
 * (c) Schrodinger, Inc.
 */

#pragma once

#include <cassert>
#include <cstring>
#include <ostream>

namespace pymol
{

/**
 * Describes a null-terminated read-only string.
 *
 * The length is not stored. Calling `size()` or `end()` calculates the length
 * on the fly in `O(size())`. All other methods behave exactly like
 * `std::string_view` with the same complexity.
 *
 * Constructing from a nullptr is undefined. See `null_safe_zstring_view` for a
 * null-safe version.
 */
class zstring_view
{
  typedef char CharT;

  const CharT* m_data = nullptr;

public:
  typedef const CharT* const_pointer;
  typedef const CharT* const_iterator;
  typedef const CharT& const_reference;
  typedef size_t size_type;
  typedef CharT value_type;

  enum : size_type { npos = size_type(-1) };

  constexpr zstring_view() = default;
  constexpr zstring_view(const CharT* s) : m_data(s) {}

  /**
   * Construct from any type which has a `c_str()` method
   */
  template <typename T, typename = decltype(&T::c_str)>
  constexpr zstring_view(const T& t) : m_data(t.c_str())
  {
  }

  constexpr const_iterator begin() const noexcept { return m_data; }
  /**
   * Like `std::string_view::end` but has complexity `O(size())`
   */
  const_iterator end() const noexcept { return m_data + size(); }

  constexpr const_reference operator[](size_type pos) const
  {
    return m_data[pos];
  }

  constexpr const_pointer data() const noexcept { return m_data; }
  constexpr const_pointer c_str() const noexcept { return m_data; }

  /**
   * True if not viewing a NULL pointer
   */
  constexpr explicit operator bool() const noexcept { return m_data; }

  /**
   * Position of the null-terminator. Complexity `O(size())`
   */
  size_type size() const noexcept { return std::strlen(m_data); }

  constexpr bool empty() const noexcept { return !m_data[0]; }

  void remove_prefix(size_type n) { m_data += n; }

  // SKIP swap

  size_type copy(CharT* dest, size_type count, size_type pos = 0) const
  {
    assert(pos <= size()); // std::out_of_range with std::string_view
    size_type i = 0;
    for (; i != count && m_data[i + pos]; ++i)
      dest[i] = m_data[i + pos];
    return i;
  }

  constexpr zstring_view substr(size_type pos) const noexcept
  {
    return m_data + pos;
  }

  int compare(zstring_view v) const noexcept
  {
    return std::strcmp(m_data, v.m_data);
  }

  bool starts_with(zstring_view x) const noexcept
  {
    for (auto it_x = x.begin(), it = begin(); *it_x; ++it_x, ++it) {
      if (*it_x != *it)
        return false;
    }
    return true;
  }
  constexpr bool starts_with(CharT x) const noexcept { return m_data[0] == x; }

  bool ends_with(zstring_view x) const noexcept
  {
    auto n = size(), n_x = x.size();
    return n >= n_x && std::strcmp(m_data + n - n_x, x.m_data) == 0;
  }
  bool ends_with(CharT x) const noexcept
  {
    auto n = size();
    return n > 0 && m_data[n - 1] == x;
  }

  size_type find(zstring_view v, size_type pos = 0) const noexcept
  {
    const_pointer p = std::strstr(m_data + pos, v.m_data);
    return p ? p - m_data : npos;
  }
  size_type find(CharT ch, size_type pos = 0) const noexcept
  {
    const_pointer p = ch ? std::strchr(m_data + pos, ch) : nullptr;
    return p ? p - m_data : npos;
  }

  size_type find_first_of(zstring_view v, size_type pos = 0) const noexcept
  {
    for (auto p = m_data + pos; *p; ++p) {
      if (std::strchr(v.m_data, *p) != nullptr)
        return p - m_data;
    }
    return npos;
  }
  size_type find_first_of(CharT ch, size_type pos = 0) const noexcept
  {
    return find(ch, pos);
  }

  size_type find_first_not_of(zstring_view v, size_type pos = 0) const noexcept
  {
    for (auto p = m_data + pos; *p; ++p) {
      if (std::strchr(v.m_data, *p) == nullptr)
        return p - m_data;
    }
    return npos;
  }
  size_type find_first_not_of(CharT ch, size_type pos = 0) const noexcept
  {
    for (auto p = m_data + pos; *p; ++p) {
      if (ch != *p)
        return p - m_data;
    }
    return npos;
  }

  // SKIP find_last_of
  // SKIP find_last_not_of
};

inline bool operator==(zstring_view lhs, zstring_view rhs) noexcept
{
  return lhs.compare(rhs) == 0;
}
inline bool operator!=(zstring_view lhs, zstring_view rhs) noexcept
{
  return lhs.compare(rhs) != 0;
}
inline bool operator<(zstring_view lhs, zstring_view rhs) noexcept
{
  return lhs.compare(rhs) < 0;
}
inline bool operator>(zstring_view lhs, zstring_view rhs) noexcept
{
  return lhs.compare(rhs) > 0;
}
inline bool operator<=(zstring_view lhs, zstring_view rhs) noexcept
{
  return lhs.compare(rhs) <= 0;
}
inline bool operator>=(zstring_view lhs, zstring_view rhs) noexcept
{
  return lhs.compare(rhs) >= 0;
}

/**
 * Describes a null-terminated read-only string which is never nullptr. If
 * constructred from a nullptr, it will point to an empty string instead.
 */
class null_safe_zstring_view : public zstring_view
{
public:
  /**
   * Construct from the empty string.
   */
  constexpr null_safe_zstring_view() : zstring_view("") {}

  /**
   * Copy constructor
   */
  constexpr null_safe_zstring_view(const null_safe_zstring_view& v)
      : zstring_view(v.data())
  {
  }

  /**
   * Construct from a null-terminated character string or from the empty string
   * if `s` is nullptr.
   */
  constexpr null_safe_zstring_view(const char* s) : zstring_view(s ? s : "") {}

  /**
   * Construct from any type which has a `c_str()` method under the assumption
   * that `c_str()` never returns nullptr. Assumption is not true for
   * `zstring_view`, see special overload.
   */
  template <typename T, typename = decltype(&T::c_str)>
  constexpr null_safe_zstring_view(const T& t) : zstring_view(t.c_str())
  {
  }

  /**
   * Construct from a `zstring_view` or from the empty string if `v` points to
   * nullptr.
   */
  constexpr null_safe_zstring_view(zstring_view v)
      : zstring_view(v.data() ? v.data() : "")
  {
  }
};

} // namespace pymol

namespace std
{
inline ostream& operator<<(ostream& os, pymol::zstring_view v)
{
  os << v.c_str();
  return os;
}

template <> struct hash<pymol::zstring_view> {
  size_t operator()(pymol::zstring_view v) const
  {
    // OVLexicon hash function (derived from djb2)
    size_t i = 0;
    size_t hashval = v[0] << 7;
    for (; v[i]; ++i) {
      hashval += (hashval << 5) + v[i];
    }
    return hashval ^ i;
  }
};
} // namespace std
