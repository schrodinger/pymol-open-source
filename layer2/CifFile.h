/*
 * CIF tokenizer
 *
 * (c) 2014 Schrodinger, Inc.
 */

#ifndef _H_CIFFILE
#define _H_CIFFILE

#include <vector>
#include <map>
#include <stdexcept>

#include <string.h>

#ifdef WIN32
  #define strcasecmp(s1, s2) _stricmp(s1, s2)
  #define strncasecmp(s1, s2, n) _strnicmp(s1, s2, n)
#endif

/*
 * C string comparison class
 */
struct strless2_t {
  bool operator()(const char * a, const char * b) const {
    return strcmp(a, b) < 0;
  }
};

// cif data types
class cif_data;
class cif_loop;
class cif_array;
typedef std::vector<cif_loop*> v_cifloopp_t;
typedef std::map<const char*, cif_data*, strless2_t> m_str_cifdatap_t;
typedef std::map<const char*, cif_array, strless2_t> m_str_cifarray_t;

// atof with uncertanty notation handling
double scifloat(const char *);

/*
 * Class for reading CIF files.
 * Parses the entire file and exposes its data blocks.
 */
class cif_file {
public:
  m_str_cifdatap_t datablocks;

  // constructors & destructor
  cif_file(const char* filename, const char* contents=NULL);
  ~cif_file();

private:
  char * contents;

  std::vector<char*> tokens;

  // methods
  bool parse();
};

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

/*
 * High-level access to CIF arrays
 */
class cif_array {
private:
  // column index, -1 if not in loop
  short col;

  // pointer to either loop or single value
  union {
    const cif_loop * loop;
    const char * value;
  } pointer;

  // methods
  const char * get_value_raw(int row = 0) const;
  const char * get_value(int row = 0) const;

public:
  // point this array to a loop (only for parsing)
  void set_loop(const cif_loop * loop, short col_) {
    col = col_;
    pointer.loop = loop;
  };

  // point this array to a single value (only for parsing)
  void set_value(const char * value) {
    col = -1;
    pointer.value = value;
  };

  // constructor
  cif_array() {
  };

  // constructor (only needed for EMPTY_ARRAY)
  cif_array(const char * value) {
    set_value(value);
  };

  // get the number of elements in this array
  int get_nrows() const {
    return (col < 0) ? 1 : pointer.loop->nrows;
  };

  // get element as string, integer or double. If index is out of bounds,
  // then return a default value.
  const char * as_s(int row = 0) const;
  int          as_i(int row = 0, int d = 0) const;
  double       as_d(int row = 0, double d = 0.0) const;

  // true if value in ['.', '?']
  bool is_missing(int row = 0) const {
    return !get_value(row);
  }

  // templated getter
  template <typename T> T as(int row = 0) const;

  // get a copy of the entire array
  template <typename T> std::vector<T> to_vector() const {
    int n = get_nrows();
    std::vector<T> v;
    v.reserve(n);
    for (int i = 0; i < n; ++i)
      v.push_back(as<T>(i));
    return v;
  }
};

/*
 * High-level access to CIF data blocks
 */
class cif_data {
  friend class cif_file;

private:
  m_str_cifarray_t dict;
  m_str_cifdatap_t saveframes;

  // only needed for freeing
  v_cifloopp_t loops;

public:
  // Get a pointer to array or NULL if not found
  const cif_array * get_arr(const char * key, const char * alias1=NULL, const char * alias2=NULL) const;

  // Get a pointer to array or to a default value if not found
  const cif_array * get_opt(const char * key, const char * alias1=NULL, const char * alias2=NULL) const;

  // destructor
  ~cif_data();
};

#endif
// vi:sw=2:ts=2
