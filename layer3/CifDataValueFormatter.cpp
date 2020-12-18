/*
 * Cif export helper stuff
 *
 * (c) 2016 Schrodinger, Inc.
 */

#include "os_std.h"
#include "CifDataValueFormatter.h"
#include "strcasecmp.h"

/**
 * Simplified whitespace test. Returns true if string `s` contains any character
 * between 0x01 and 0x20 (this range contains all whitespace characters and
 * no printable characters)
 */
static
bool has_whitespace(const char * s) {
  for (; *s; ++s)
    if (*s <= ' ')
      return true;
  return false;
}

/**
 * Return true if `s` contains the given quote character followed by white space
 */
static
bool has_quotespace(const char * s, char quote) {
  for (; (s = strchr(s, quote)); ++s)
    if (s[1] && s[1] <= ' ')
      return true;
  return false;
}

/**
 * Return true if `s` is a "simple data value" according to the CIF
 * syntax specification
 */
static
bool cif_is_simpledatavalue(const char * s) {
  return (
      // first character is special
      !strchr("_#$'\"[];", s[0]) &&
      // whitespace
      !has_whitespace(s) &&
      // special values '.' (inapplicable) and '?' (unknown)
      !((s[0] == '.' || s[0] == '?') && !s[1]) &&
      // prefix is special
      strncasecmp("data_", s, 5) &&
      strncasecmp("save_", s, 5) &&
      // reserved words
      strcasecmp("loop_", s) &&
      strcasecmp("stop_", s) &&
      strcasecmp("global_", s));
}

std::string & CifDataValueFormatter::nextbuf() {
  // advance circular pointer
  m_i = (m_i + 1) % m_buf.size();
  return m_buf[m_i];
}

const char * CifDataValueFormatter::quoted(const char * s) {
  const char * quote = nullptr;

  if (!strchr(s, '\n')) {
    if (!has_quotespace(s, '\'')) {
      quote = "'";
    } else if (!has_quotespace(s, '"')) {
      quote = "\"";
    }
  }

  if (!quote) {
    quote = "\n;";
    if (strstr(s, quote)) {
      printf(" CIF-Warning: data value contains unquotable <newline><semicolon>\n");
      return "<UNQUOTABLE>";
    }
  }

  return nextbuf().assign(quote).append(s).append(quote).c_str();
}

const char * CifDataValueFormatter::operator() (const char * s, const char * d) {
  if (!s[0])
    return d;

  if (cif_is_simpledatavalue(s))
    return s;

  return quoted(s);
}

const char * CifDataValueFormatter::operator() (char c, const char * d) {
  return (*this)(nextbuf().assign(1, c).c_str(), d);
}
