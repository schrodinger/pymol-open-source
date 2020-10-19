/*
 * CIF tokenizer
 *
 * All keys are canonicalized to lowercase
 *
 * (c) 2014 Schrodinger, Inc.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cassert>
#include <iostream>
#include <string>
#include <vector>

#include "CifFile.h"
#include "File.h"
#include "MemoryDebug.h"
#include "strcasecmp.h"

namespace pymol {
namespace _cif_detail {

template <> const char* raw_to_typed(const char* s) { return s; }
template <> std::string raw_to_typed(const char* s) { return s; }
template <> char        raw_to_typed(const char* s) { return s[0]; }
template <> int         raw_to_typed(const char* s) { return atoi(s); }

/**
 * Convert to floating point number, ignores uncertainty notation
 * 1.23(45)e2 -> 1.23e2
 */
template <> double raw_to_typed(const char* s)
{
  const char *close, *open = strchr(s, '(');
  if (open && (close = strchr(open, ')'))) {
    return atof(std::string(s, open - s).append(close + 1).c_str());
  }
  return atof(s);
}

template <> float raw_to_typed(const char* s)
{
  return static_cast<float>(raw_to_typed<double>(s));
}

} // namespace _cif_detail

// basic IO and string handling

// Return true if "c" is whitespace or null
static bool iswhitespace0(char c) {
  return strchr(" \t\r\n", c) ? true : false;
}

// Return true if "c" is whitespace
static bool iswhitespace(char c) {
  return (c && iswhitespace0(c));
}

// Return true if "c" is line feed or carriage return
static bool islinefeed(char c) {
  return (c == '\r' || c == '\n');
}

// Return true if "c" is line feed or carriage return or null
static bool islinefeed0(char c) {
  return (!c || islinefeed(c));
}

// Return true if "c" is double or single quote
static bool isquote(char c) {
  return (c == '"' || c == '\'');
}

// FreeBSD name conflict
#ifdef isspecial
#undef isspecial
#endif

// Return true if token is a STAR keyword
static bool isspecial(const char *token) {
  return (token[0] == '_'
      || strncasecmp("data_", token, 5) == 0
      || strncasecmp("save_", token, 5) == 0
      || strcasecmp("loop_", token) == 0
      || strcasecmp("stop_", token) == 0
      || strcasecmp("global_", token) == 0);
}

// convert all chars to lowercase
static void tolowerinplace(char *p) {
  for (; *p; p++) {
    if (*p <= 'Z' && *p >= 'A')
      *p -= 'Z' - 'z';
  }
}

// CIF stuff

static const cif_array EMPTY_ARRAY(nullptr);

/*
 * Class to store CIF loops. Only for parsing, do not use in any higher level
 * reading functions.
 */
class cif_loop {
public:
  int ncols;
  int nrows;
  const char **values;

  // methods
  const char * get_value_raw(int row, int col) const;
};

// get table value, return NULL if indices out of bounds
const char * cif_loop::get_value_raw(int row, int col) const {
  if (row >= nrows)
    return nullptr;
  return values[row * ncols + col];
}

// get the number of elements in this array
unsigned cif_array::size() const {
  return (col == NOT_IN_LOOP) ? 1 : pointer.loop->nrows;
}

/// Get array value, return NULL if `pos >= size()` or value in ['.', '?']
const char* cif_array::get_value_raw(unsigned pos) const
{
  if (col == NOT_IN_LOOP)
    return (pos > 0) ? nullptr : pointer.value;
  return pointer.loop->get_value_raw(pos, col);
}

// true if all values in ['.', '?']
bool cif_array::is_missing_all() const {
  for (unsigned i = 0, n = size(); i != n; ++i) {
    if (!is_missing(i))
      return false;
  }

  return true;
}

/**
 * Get a pointer to array or NULL if not found
 *
 * Can lookup different aliases, the first one found is returned.
 * Also supports an alias shortcut for the trivial case where mmCIF uses
 * a colon and CIF uses an underscore: `get_arr("_foo?bar")` is identical to
 * `get_arr("_foo.bar", "_foo_bar")`
 *
 * @param key data name, must be lower case
 */
const cif_array * cif_data::get_arr(const char * key) const {
  const char* p = strchr(key, '?');
  decltype(m_dict)::const_iterator it;

#ifndef NDEBUG
  for (const char* q = key; *q; ++q) {
    assert("key must be lower case" && !('Z' >= *q && *q >= 'A'));
  }
#endif

  // support alias shortcut: '?' matches '.' and '_'
  if (p != nullptr) {
    std::string tmp(key);
    // replace '?' by '.' or '_'
    tmp[p - key] = '.';
    if ((it = m_dict.find(tmp.c_str())) != m_dict.end())
      return &it->second;
    tmp[p - key] = '_';
    if ((it = m_dict.find(tmp.c_str())) != m_dict.end())
      return &it->second;
  } else {
    if ((it = m_dict.find(key)) != m_dict.end())
      return &it->second;
  }

  return nullptr;
}

const cif_array* cif_data::empty_array() {
  return &EMPTY_ARRAY;
}

const cif_data* cif_data::get_saveframe(const char* code) const {
  auto it = m_saveframes.find(code);
  if (it != m_saveframes.end())
    return &it->second;
  return nullptr;
}

bool cif_file::parse_file(const char* filename) {
  char* contents = FileGetContents(filename, nullptr);

  if (!contents) {
    error(std::string("failed to read file ").append(filename).c_str());
    return false;
  }

  return parse(std::move(contents));
}

bool cif_file::parse_string(const char* contents) {
  return parse(std::move(mstrdup(contents)));
}

void cif_file::error(const char* msg) {
  std::cout << "ERROR " << msg << std::endl;
}

// constructor
cif_file::cif_file(const char* filename, const char* contents_) {
  if (contents_) {
    parse_string(contents_);
  } else if (filename) {
    parse_file(filename);
  }
}

// constructor
cif_file::cif_file() = default;
cif_file::cif_file(cif_file&&) = default;

// move assignment
cif_file& cif_file::operator=(cif_file&&) = default;

// destructor
cif_file::~cif_file() = default;

bool cif_file::parse(char*&& p) {
  m_datablocks.clear();
  m_tokens.clear();
  m_contents.reset(p);

  if (!p) {
    error("parse(nullptr)");
    return false;
  }

  auto& tokens = m_tokens;
  char quote;
  char prev = '\0';

  std::vector<bool> keypossible;

  // tokenize
  while (true) {
    while (iswhitespace(*p))
      prev = *(p++);

    if (!*p)
      break;

    if (*p == '#') {
      while (!(islinefeed0(*++p)));
      prev = *p;
    } else if (isquote(*p)) { // will NULL the closing quote
      quote = *p;
      keypossible.push_back(false);
      tokens.push_back(p + 1);
      while (*++p && !(*p == quote && iswhitespace0(p[1])));
      if (*p)
        *(p++) = 0;
      prev = *p;
    } else if (*p == ';' && islinefeed(prev)) {
      // multi-line tokens start with ";" and end with "\n;"
      // multi-line tokens cannot be keys, only values.
      keypossible.push_back(false);
      tokens.push_back(p + 1);
      // advance until `\n;`
      while (*++p && !(islinefeed(*p) && p[1] == ';'));
      // step to next line and null the line feed
      if (*p) {
        *p = 0;
        // \r\n on Windows)
        if (p - 1 > tokens.back() && *(p - 1) == '\r') {
          *(p - 1) = 0;
        }
        p += 2;
      }
      prev = ';';
    } else { // will null the whitespace
      char * q = p++;
      while (!iswhitespace0(*p)) ++p;
      prev = *p;
      if (p - q == 1 && (*q == '?' || *q == '.')) {
        // store values '.' (inapplicable) and '?' (unknown) as null-pointers
        q = nullptr;
        keypossible.push_back(false);
      } else {
        if (*p)
          *(p++) = 0;
        keypossible.push_back(true);
      }
      tokens.push_back(q);
    }
  }

  cif_data* current_frame = nullptr;
  std::vector<cif_data*> frame_stack;
  std::unique_ptr<cif_data> global_block;
  decltype(m_datablocks) datablocksnew;

  // parse into dictionary
  for (unsigned int i = 0, n = tokens.size(); i < n; i++) {
    if (!keypossible[i]) {
      error("expected key (1)");
      return false;
    } else if (tokens[i][0] == '_') {
      if (!current_frame) {
        error("missing data_ (unexpected data name)");
        return false;
      }

      if (i + 1 == n) {
        error("truncated");
        return false;
      }

      tolowerinplace(tokens[i]);
      current_frame->m_dict[tokens[i]].set_value(tokens[i + 1]);

      i++;
    } else if (strcasecmp("loop_", tokens[i]) == 0) {
      if (!current_frame) {
        error("missing data_ (unexpected loop)");
        return false;
      }

      int ncols = 0;
      int nrows = 0;
      cif_loop *loop = nullptr;

      // loop data
      loop = new cif_loop;
      current_frame->m_loops.emplace_back(loop);

      // columns
      while (++i < n && keypossible[i] && tokens[i][0] == '_') {
        tolowerinplace(tokens[i]);

        current_frame->m_dict[tokens[i]].set_loop(loop, ncols);

        ncols++;
      }

      if (loop) {
        // loop data
        loop->values = (const char **) &tokens[i];
        loop->ncols = ncols;
      }

      // rows
      while (i < n && !(keypossible[i] && isspecial(tokens[i]))) {
        i += ncols;

        if (i > n) {
          error("truncated loop");
          return false;
        }

        nrows++;
      }

      // loop data
      if (loop) {
        loop->nrows = nrows;
      }

      i--;

    } else if (strncasecmp("data_", tokens[i], 5) == 0) {
      datablocksnew.emplace_back();
      current_frame = &datablocksnew.back();
      current_frame->m_code = tokens[i] + 5;
      frame_stack = {current_frame};

    } else if (strncasecmp("global_", tokens[i], 5) == 0) {
      // STAR feature, not supported in CIF
      current_frame = new cif_data;
      global_block.reset(current_frame);
      frame_stack = {current_frame};

    } else if (strncasecmp("save_", tokens[i], 5) == 0) {
      if (tokens[i][5]) {
        // begin
        if (!current_frame) {
          error("top-level save_");
          return false;
        }

        const char * key(tokens[i] + 5);
        current_frame = &current_frame->m_saveframes[key];
        frame_stack.push_back(current_frame);
      } else {
        // end
        if (frame_stack.size() < 2) {
          error("unexpected save_");
          return false;
        }

        frame_stack.pop_back();
        current_frame = frame_stack.back();
      }
    } else {
      error("expected key (2)");
      return false;
    }
  }

  m_datablocks = std::move(datablocksnew);

  return true;
}

} // namespace pymol

// vi:sw=2:ts=2
