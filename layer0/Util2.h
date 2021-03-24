/*
 * Utility functions
 *
 * (c) 2015 Schrodinger, Inc.
 */

#pragma once

#include <algorithm>
#include <array>
#include <string>
#include <vector>

#include <string.h>
#include <sstream>

#include "pymol/type_traits.h"
#include "pymol/zstring_view.h"

std::vector<std::string> strsplit(const std::string &s, char delim=0);

bool cstrlessnat(const char * a, const char * b);
bool strlessnat(const std::string& a, const std::string& b);

/**
 * C string comparison class
 */
struct cstrless_t {
  bool operator()(const char * a, const char * b) const {
    return strcmp(a, b) < 0;
  }
};

bool p_strstartswith(const char * s, const char * prefix);
bool p_strcasestartswith(const char * s, const char * prefix);

namespace pymol
{
namespace join_to_string_detail
{
inline void join_to_string_impl(std::ostringstream&){};

template <typename T, typename... OtherTs>
void join_to_string_impl(std::ostringstream& stream, T&& t, OtherTs&&... ts)
{
  stream << std::forward<T>(t);
  join_to_string_impl(stream, std::forward<OtherTs>(ts)...);
}
} // namespace join_to_string_detail

/**
 * Joins data into std::string.
 * @param ts printable types to be joined into a string
 * @return joined string
 * Note: Types must be supported be std::ostringstream::operator<<
 */

template <typename... PrintableTs>
std::string join_to_string(PrintableTs&&... ts)
{
  std::ostringstream stream;
  join_to_string_detail::join_to_string_impl(
      stream, std::forward<PrintableTs>(ts)...);
  return stream.str();
}

namespace string_format_detail
{
template <typename T> const T& fwdArgs(const T& t)
{
#if defined(__clang__) || !defined(__GNUC__) || __GNUC__ >= 5
  // Not available in GCC 4.8
  static_assert(std::is_trivially_copyable<T>::value, "");
#endif
  return t;
}

inline const char* fwdArgs(const std::string& t)
{
  return t.c_str();
}

inline const char* fwdArgs(const pymol::zstring_view& t)
{
  return t.c_str();
}

template <typename... FmtArgs>
std::string string_format_impl(const char* const fmt, FmtArgs&&... fmtargs)
{
  auto size = snprintf(nullptr, 0, fmt, fmtargs...);
  std::string tmp(size, ' ');
  snprintf(&tmp[0], size + 1, fmt, fmtargs...);
  return tmp;
}

} // namespace string_format_detail

/**
 * C++ version of sprintf
 * @param fmt formatting string
 * @param fmtargs formatting string arguments
 * @return result as std::string
 *
 * Note: For std::string arguments, the underlying null-terminated string will
 * be used.
 */

template <std::size_t N, typename... FmtArgs>
std::string string_format(const char (&fmt)[N], FmtArgs&&... fmtargs)
{
  static_assert(N > 1, "Format string must not be empty");
  return string_format_detail::string_format_impl(
      fmt, string_format_detail::fwdArgs(std::forward<FmtArgs>(fmtargs))...);
}

/**
 * C++20's std::string::start_with
 * @param str string whose contents will be checked
 * @param pre candidate prefix string
 * @return true if str begins with pre
 *
 */
inline bool starts_with(pymol::zstring_view str, pymol::zstring_view pre)
{
  return str.starts_with(pre);
}

double pretty_f2d(float v);

/**
 * Compares two strings with consideration of case sensitivity
 * @param str1 first string
 * @param str2 second string
 * @param case_insensitive determines whether the comparison should not consider
 * case
 */

bool string_equal_case(pymol::zstring_view str1, pymol::zstring_view str2,
    bool case_insensitive = true);

template <typename T> struct cache_value {
  using value_type = T;

  value_type value = value_type{};

  cache_value() = default;

  cache_value(T t)
      : value(std::move(t))
  {
  }

  cache_value(const cache_value&) {}
  cache_value& operator=(const cache_value&) { return *this; }
  cache_value(cache_value&& other)
  {
    std::swap(value, other.value);
  }
  cache_value& operator=(cache_value&& other)
  {
    std::swap(value, other.value);
    return *this;
  }

  operator T&() { return value; }
};

template <typename T, std::size_t N>
using cache_array = std::array<cache_value<T>, N>;

} // namespace pymol
