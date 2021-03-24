/*
 * Utility functions
 *
 * (c) 2015 Schrodinger, Inc.
 */

#include <cstdio>
#include <cctype>
#include <cmath>
#include <sstream>
#include <string>
#include <vector>

#include "Util2.h"
#include "pymol/algorithm.h"

/**
 * strsplit: split string `s` by character `delim`
 *
 * If `delim` is null, then do whitespace split.
 */
std::vector<std::string> strsplit(const std::string &s, char delim) {
  std::vector<std::string> elems;
  std::istringstream iss(s);
  std::string item;

  if (delim) {
    // token split
    while (std::getline(iss, item, delim))
      elems.push_back(item);
  } else {
    // whitespace split
    while (iss >> item)
      elems.push_back(item);
  }

  return elems;
}

/**
 * Natural string compare: F1 < F2 < F10
 *
 * Return true if a < b
 */
bool cstrlessnat(const char * a, const char * b) {
  if (!b[0])
    return false;
  if (!a[0])
    return true;

  bool a_digit = isdigit(a[0]);
  bool b_digit = isdigit(b[0]);

  if (a_digit && !b_digit)
    return true;
  if (!a_digit && b_digit)
    return false;

  if (!a_digit && !b_digit) {
    if (a[0] != b[0])
      return (a[0] < b[0]);

    return cstrlessnat(a + 1, b + 1);
  }

  int ia, ib, na, nb;
  sscanf(a, "%d%n", &ia, &na);
  sscanf(b, "%d%n", &ib, &nb);

  if (ia != ib)
    return ia < ib;

  return cstrlessnat(a + na, b + nb);
}

/**
 * Natural string compare: F1 < F2 < F10
 */
bool strlessnat(const std::string& a, const std::string& b) {
  return cstrlessnat(a.c_str(), b.c_str());
}

/**
 * Return true if s starts with the specified prefix, false otherwise.
 */
bool p_strstartswith(const char * s, const char * prefix) {
  while (*prefix)
    if (*s++ != *prefix++)
      return false;
  return true;
}

/**
 * case-insensitive version of p_strstartswith
 */
bool p_strcasestartswith(const char * s, const char * prefix) {
  for (; *prefix; ++s, ++prefix)
    if (*s != *prefix && tolower(*s) != tolower(*prefix))
      return false;
  return true;
}

namespace pymol
{

/**
 * Convert float to double, rounded to decimal precision.
 */
double pretty_f2d(float f)
{
  if (f == 0.0f) {
    // don't call log10(0.0)
    return 0.0;
  }

  int digits = std::ceil(std::log10(std::fabs(f)));
  auto factor = std::pow(10.0L, 7 - digits);
  return std::round(f * factor) / factor;
}

bool string_equal_case(
    pymol::zstring_view str1, pymol::zstring_view str2, bool case_insensitive)
{
  return pymol::ranges::equal(str1, str2,
      [case_insensitive](char c1, char c2) {
        return case_insensitive ? std::tolower(c1) == std::tolower(c2)
                                : c1 == c2;
      });
}

} // namespace pymol
