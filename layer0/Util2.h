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
} // namespace pymol

