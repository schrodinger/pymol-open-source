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

#include <vector>
#include <iostream>

#include "CifFile.h"
#include "File.h"
#include "MemoryDebug.h"
#include "strcasecmp.h"

// basic IO and string handling

/*
 * atof which ignores uncertainty notation
 * 1.23(45)e2 -> 1.23e2
 */
double scifloat(const char *str) {
  const char *close, *open = strchr(str, '(');
  if (open && (close = strchr(open, ')'))) {
    double value;
    char *copy = strdup(str);
    strcpy(copy + (open - str), close + 1);
    value = atof(copy);
    free(copy);
    return value;
  }
  return atof(str);
}

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

static const char * EMPTY_STRING = "";
static cif_array EMPTY_ARRAY(NULL);

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
    return NULL;
  return values[row * ncols + col];
}

// get the number of elements in this array
int cif_array::get_nrows() const {
  return (col < 0) ? 1 : pointer.loop->nrows;
}

// get array value, return NULL if row-index out of bounds
// or value in ['.', '?']
const char * cif_array::get_value(int row) const {
  if (col < 0)
    return (row > 0) ? NULL : pointer.value;
  return pointer.loop->get_value_raw(row, col);
}

// get array value, return an empty string if missing
const char * cif_array::as_s(int row) const {
  const char * s = get_value(row);
  return s ? s : EMPTY_STRING;
}

// get array value as integer, return d (default 0) if missing
int cif_array::as_i(int row, int d) const {
  const char * s = get_value(row);
  return s ? atoi(s) : d;
}

// get array value as double, return d (default 0.0) if missing
double cif_array::as_d(int row, double d) const {
  const char * s = get_value(row);
  return s ? scifloat(s) : d;
}

// true if all values in ['.', '?']
bool cif_array::is_missing_all() const {
  int n = get_nrows();

  for (int i = 0; i < n; ++i) {
    if (!is_missing(i))
      return false;
  }

  return true;
}

// templated getters
template <> const char* cif_array::as<const char* >(int row) const { return get_value(row); }
template <> std::string cif_array::as<std::string >(int row) const { return as_s(row); }
template <> int         cif_array::as<int         >(int row) const { return as_i(row); }
template <> double      cif_array::as<double      >(int row) const { return as_d(row); }
template <> float       cif_array::as<float       >(int row) const { return as_d(row); }

/*
 * Get a pointer to array or NULL if not found
 *
 * Can lookup up to 3 different aliases, the first one found is returned.
 * Also supports an alias shortcut for the trivial case where mmCIF uses
 * a colon and CIF uses an underscore: (key="_foo?bar") is identical to
 * (key="_foo.bar", alias1="_foo_bar")
 */
const cif_array * cif_data::get_arr(const char * key, const char * alias1, const char * alias2) const {
  const char * p;
  const char * aliases[] = {alias1, alias2, NULL};
  m_str_cifarray_t::const_iterator it;

  for (int j = 0; key; key = aliases[j++]) {
    // support alias shortcut: '?' matches '.' and '_'
    if ((p = strchr(key, '?'))) {
      std::string tmp(key);
      for (const char * d = "._"; *d; ++d) {
        // replace '?' by '.' or '_'
        tmp[p - key] = *d;
        if ((it = dict.find(tmp.c_str())) != dict.end())
          return &it->second;
      }
    } else {
      if ((it = dict.find(key)) != dict.end())
        return &it->second;
    }
  }

  return NULL;
}

// Get a pointer to array or to a default value if not found
const cif_array * cif_data::get_opt(const char * key, const char * alias1, const char * alias2) const {
  const cif_array * arr = get_arr(key, alias1, alias2);
  if (arr == NULL)
    return &EMPTY_ARRAY;
  return arr;
}

// constructor
cif_file::cif_file(const char* filename, const char* contents_) {
  if (contents_) {
    contents = mstrdup(contents_);
  } else {
    contents = FileGetContents(filename, NULL);
    if (!contents)
      std::cerr << "ERROR: Failed to load file '" << filename << "'" << std::endl;
  }

  if (contents)
    parse();
}

// destructor
cif_file::~cif_file() {
  for (auto& datablock : datablocks)
    delete datablock.second;

  if (contents)
    mfree(contents);
}

// destructor
cif_data::~cif_data() {
  for (auto& saveframe : saveframes)
    delete saveframe.second;

  for (auto& loop : loops)
    delete loop;
}

// parse CIF contents
bool cif_file::parse() {
  char *p = contents;
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
    } else if (*p == ';' && islinefeed(prev)) { // will NULL the line feed before the closing semicolon
      keypossible.push_back(false);
      tokens.push_back(p + 1);
      while (*++p && !(islinefeed(*p) && p[1] == ';'));
      if (*p) {
        *p = 0;
        p += 2;
      }
      prev = ';';
    } else { // will null the whitespace
      char * q = p++;
      while (!iswhitespace0(*p)) ++p;
      prev = *p;
      if (p - q == 1 && (*q == '?' || *q == '.')) {
        // store values '.' (inapplicable) and '?' (unknown) as null-pointers
        q = NULL;
        keypossible.push_back(false);
      } else {
        if (*p)
          *(p++) = 0;
        keypossible.push_back(true);
      }
      tokens.push_back(q);
    }
  }

  cif_data *current_data = NULL, *current_frame = NULL, *global_block = NULL;

  // parse into dictionary
  for (unsigned int i = 0, n = tokens.size(); i < n; i++) {
    if (!keypossible[i]) {
      std::cout << "ERROR" << std::endl;
      break;
    } else if (tokens[i][0] == '_') {
      if (i + 1 == n) {
        std::cout << "ERROR truncated" << std::endl;
        break;
      }

      if (current_frame) {
        tolowerinplace(tokens[i]);
        current_frame->dict[tokens[i]].set_value(tokens[i + 1]);
      }

      i++;
    } else if (strcasecmp("loop_", tokens[i]) == 0) {
      int ncols = 0;
      int nrows = 0;
      cif_loop *loop = NULL;

      // loop data
      if (current_frame) {
        loop = new cif_loop;

        // add to loops list
        current_frame->loops.push_back(loop);
      }

      // columns
      while (++i < n && keypossible[i] && tokens[i][0] == '_') {
        tolowerinplace(tokens[i]);

        if (current_frame) {
          current_frame->dict[tokens[i]].set_loop(loop, ncols);
        }

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
          std::cout << "ERROR truncated loop" << std::endl;
          break;
        }

        nrows++;
      }

      // loop data
      if (loop) {
        loop->nrows = nrows;
      }

      i--;

    } else if (strncasecmp("data_", tokens[i], 5) == 0) {
      const char * key(tokens[i] + 5);
      datablocks[key] = current_data = current_frame = new cif_data;

    } else if (strncasecmp("global_", tokens[i], 5) == 0) {
      // STAR feature, not supported in CIF
      global_block = current_data = current_frame = new cif_data;

    } else if (strncasecmp("save_", tokens[i], 5) == 0) {
      if (tokens[i][5] && current_data) {
        // begin
        const char * key(tokens[i] + 5);
        current_data->saveframes[key] = current_frame = new cif_data;
      } else {
        // end
        current_frame = current_data;
      }
    } else {
      std::cout << "ERROR" << std::endl;
      break;
    }
  }

  if (global_block)
    delete global_block;

  return true;
}

// vi:sw=2:ts=2
