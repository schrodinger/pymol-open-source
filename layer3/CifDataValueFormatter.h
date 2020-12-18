/*
 * Cif export helper stuff
 *
 * (c) 2016 Schrodinger, Inc.
 */

#pragma once

#include <string>
#include <vector>

#include "os_std.h"

/**
 * Callable class to format a CIF (STAR) data value. If the string is
 * a "simple data vlaue", then return it as-is. Otherwise return a quoted
 * copy of the string, pointing to an internal memory buffer.
 */
class CifDataValueFormatter {
  int m_i;

  std::string & nextbuf();
  const char * quoted(const char * s);

public:
  // circular buffer for quoted strings
  std::vector<std::string> m_buf;

  CifDataValueFormatter(size_t size=1) : m_i(0), m_buf(size) {}

  const char * operator() (const char * s, const char * d=".");
  const char * operator() (char c, const char * d=".");
};
