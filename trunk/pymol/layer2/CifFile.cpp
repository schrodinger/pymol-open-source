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
#include <stdexcept>

#include "CifFile.h"
#include "File.h"

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

// get table value, return NULL if indices out of bounds
const char * cif_loop::get_value_raw(int row, int col) const {
  if (row >= nrows)
    return NULL;
  return values[row * ncols + col];
}

// get array value, return NULL if row-index out of bounds
const char * cif_array::get_value_raw(int row) const {
  if (col < 0)
    return (row > 0) ? NULL : pointer.value;
  return pointer.loop->get_value_raw(row, col);
};

// get array value, return NULL if value in ['.', '?']
const char * cif_array::get_value(int row) const {
  const char * s = get_value_raw(row);
  return (s && (s[0] == '?' || s[0] == '.') && !s[1]) ? NULL : s;
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

/*
 * Get a pointer to array or NULL if not found
 *
 * Can lookup up to 3 different aliases, the first one found is returned.
 * Also supports an alias shortcut for the trivial case where mmCIF uses
 * a colon and CIF uses an underscore: (key="_foo?bar") is identical to
 * (key="_foo.bar", alias1="_foo_bar")
 */
cif_array * cif_data::get_arr(const char * key, const char * alias1, const char * alias2) {
  std::string tmp1, tmp2;
  const char * p;

  // support alias shortcut: '?' matches '.' and '_'
  if (!alias1 && (p = strchr(key, '?'))) {
    int i = p - key;

    // '.' version
    tmp1 = key;
    tmp1[i] = '.';
    key = tmp1.c_str();

    // '_' version
    tmp2 = key;
    tmp2[i] = '_';
    alias1 = tmp2.c_str();
  }

  // dict lookup, return first hit
  m_str_cifarray_t::iterator arr, end = dict.end();
  if ((arr = dict.find(key)) != end ||
      alias1 && (arr = dict.find(alias1)) != end ||
      alias2 && (arr = dict.find(alias2)) != end)
    return &arr->second;
  return NULL;
}

// Get a pointer to array or to a default value if not found
cif_array * cif_data::get_opt(const char * key, const char * alias1, const char * alias2) {
  cif_array * arr = get_arr(key, alias1, alias2);
  if (arr == NULL)
    return &EMPTY_ARRAY;
  return arr;
}

// constructor - mode = 0: filename, mode = 1: content (free), mode = 2: content (don't free)
void cif_file::init(char * source, int mode) {
  owncontents = true;

  switch (mode) {
    case 0:
      contents = FileGetContents(source, NULL);

      if (!contents)
        std::cerr << "ERROR: Failed to load file '" << source << "'" << std::endl;

      break;
    case 2:
      owncontents = false;
    default:
      contents = source;
      break;
  }

  if (contents)
    parse();
}

// destructor
cif_file::~cif_file() {
  for (m_str_cifdatap_t::iterator it = datablocks.begin(),
      it_end = datablocks.end(); it != it_end; ++it)
    delete it->second;

  if (contents && owncontents)
    free(contents);
}

// destructor
cif_data::~cif_data() {
  for (m_str_cifdatap_t::iterator it = saveframes.begin(),
      it_end = saveframes.end(); it != it_end; ++it)
    delete it->second;

  for (v_cifloopp_t::iterator it = loops.begin(),
      it_end = loops.end(); it != it_end; ++it)
    delete *it;
}

// parse CIF contents
bool cif_file::parse() {
  char *p = contents;
  char quote;
  char prev = '\0';

  std::vector<char> codes;

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
      codes.push_back('Q');
      tokens.push_back(p + 1);
      while (*++p && !(*p == quote && iswhitespace0(p[1])));
      if (*p)
        *(p++) = 0;
      prev = *p;
    } else if (*p == ';' && islinefeed(prev)) { // will NULL the line feed before the closing semicolon
      codes.push_back('Q');
      tokens.push_back(p + 1);
      while (*++p && !(islinefeed(*p) && p[1] == ';'));
      if (*p) {
        *p = 0;
        p += 2;
      }
      prev = ';';
    } else { // will null the whitespace
      codes.push_back('R');
      tokens.push_back(p);
      while (!iswhitespace0(*p)) ++p;
      prev = *p;
      if (*p)
        *(p++) = 0;
    }
  }

  cif_data *current_data = NULL, *current_frame = NULL;

  // parse into dictionary
  for (unsigned int i = 0, n = tokens.size(); i < n; i++) {
    if (codes[i] == 'Q') {
      std::cout << "ERROR" << std::endl;
      break;
    } else if (tokens[i][0] == '_') {
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
      while (codes[++i] != 'Q' && tokens[i][0] == '_') {
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
      while (i < n && (codes[i] == 'Q' || !isspecial(tokens[i]))) {
        i += ncols;
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

    } else if (strncasecmp("save_", tokens[i], 5) == 0) {
      if (tokens[i][5]) {
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

  return true;
}

// vi:sw=2:ts=2
