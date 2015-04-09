/*
 * Utility functions
 *
 * (c) 2015 Schrodinger, Inc.
 */

#include <cstdio>
#include <cctype>
#include <sstream>
#include <string>
#include <vector>

#include "Util2.h"

/*
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

/*
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

/*
 * Natural string compare: F1 < F2 < F10
 */
bool strlessnat(const std::string& a, const std::string& b) {
  return cstrlessnat(a.c_str(), b.c_str());
}


