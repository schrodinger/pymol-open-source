/*
 * Utility functions
 *
 * (c) 2015 Schrodinger, Inc.
 */

#pragma once

#include <string>
#include <vector>

#include <string.h>
#include <sstream>

#include "pymol/type_traits.h"

std::vector<std::string> strsplit(const std::string &s, char delim=0);

bool cstrlessnat(const char * a, const char * b);
bool strlessnat(const std::string& a, const std::string& b);

/*
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
template <typename T,
    enable_if_t<std::is_same<remove_cvref_t<T>, std::string>::value>* = nullptr>
const char* fwdArgs(T&& t)
{
  return t.c_str();
}

template <typename T,
    enable_if_t<!std::is_same<remove_cvref_t<T>, std::string>::value>* = nullptr>
T&& fwdArgs(T&& t)
{
  return std::forward<T>(t);
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
  static_assert(N > 0, "Format string must not be empty");
  return string_format_detail::string_format_impl(
      fmt, string_format_detail::fwdArgs(std::forward<FmtArgs>(fmtargs))...);
}

double pretty_f2d(float v);

} // namespace pymol
