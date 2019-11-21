/*
 * Copyright (c) Schrodinger, LLC.
 *
 * Datastructure for `bond_dict`, parsed from CIF files with
 * `chem_comp_bond` records.
 *
 * This file includes a generic reference implementation (1), and an
 * optimized fast and memory efficient implementation (2) with key
 * string length limitations.
 */

#ifndef CIFBONDDICT_H
#define CIFBONDDICT_H

#include <map>
#include <algorithm>
#include <string>

#include "os_std.h"
#include "PyMOLGlobals.h"

#if 1
// (2)
// Optimized implementation, tuned on speed and memory usage.
// Only supports resn length <= 8 and atom name length <= 4.
// Working with this limitation is currently reasonable in PyMOL,
// since cResnLen=5 and cAtomNameLen=4 (see AtomInfo.h).

#include <cstdint>
#include <unordered_map>
#include <string.h>

// mapped type of bond_dict_t
class res_bond_dict_t {
  using key_type = std::int_fast64_t;
  using halfkey_t = std::int32_t;
  using mapped_type = signed char;

  static halfkey_t make_halfkey(const char* name) {
    union { halfkey_t i; char s[sizeof(halfkey_t)]; };
    strncpy(s, name, sizeof(halfkey_t));
    return i;
  }

  static key_type make_key(const char * name1, const char * name2) {
    auto i1 = make_halfkey(name1);
    auto i2 = make_halfkey(name2);

    if (i1 > i2)
      std::swap(i1, i2);

    // construct 8-byte key from the two 4-byte names
    return (((key_type)i1) << 32) | i2;
  }

  std::unordered_map<key_type, mapped_type> m_map;
  std::unordered_map<halfkey_t, std::string> m_alt_names;

  const char* get_standard_name(const char* name) const {
    auto it = m_alt_names.find(make_halfkey(name));
    if (it != m_alt_names.end()) {
      return it->second.c_str();
    }
    return name;
  }

public:
  void add_alt_name(const char* name, const char* alt) {
    m_alt_names[make_halfkey(alt)] = name;
  }

  void set(const char * name1, const char * name2, mapped_type order) {
    m_map[make_key(name1, name2)] = order;
  }

  mapped_type get(const char * name1, const char * name2) const {
    name1 = get_standard_name(name1);
    name2 = get_standard_name(name2);
    auto it = m_map.find(make_key(name1, name2));
    if (it == m_map.end())
      return -1;
    return it->second;
  }
};

// type for mapping: resn -> ((name1, name2) -> bond order)
class bond_dict_t {
  using key_type = std::int_fast64_t;
  using mapped_type = res_bond_dict_t;

  static key_type make_key(const char * s_) {
    union { key_type i; char s[sizeof(key_type)]; };
    strncpy(s, s_, sizeof(s));
    return i;
  }

  std::map<key_type, mapped_type> m_map;
  std::set<key_type> unknown_resn;

public:
  bool empty() const { return m_map.empty(); }

  mapped_type& operator[](const char* resn) { return m_map[make_key(resn)]; }

  void set_unknown(const char * resn) {
    unknown_resn.insert(make_key(resn));
  }

  const mapped_type * get(PyMOLGlobals *, const char * resn, bool try_download=true);
};

#else
// (1)
// Generic reference implementation. Supports arbitrary key lengths,
// but is slower and less memory efficient. Equivalent to:
//
// typedef std::map<std::string,
//         std::map<std::string,
//         std::map<std::string, int> > > bond_dict_t;
//
// with added get/set methods.

#include <string>

class res_bond_dict_t {
  std::map<std::string, std::map<std::string, int> > m_data;

public:
  void set(const char * name1, const char * name2, int order) {
    if (strcmp(name1, name2) < 0)
      std::swap(name1, name2);
    m_data[name1][name2] = order;
  }

  int get(const char * name1, const char * name2) const {
    if (strcmp(name1, name2) < 0)
      std::swap(name1, name2);
    try {
      return m_data.at(name1).at(name2);
    } catch (const std::out_of_range& e) {
      return -1;
    }
  }
};

class bond_dict_t {
  std::map<std::string, res_bond_dict_t> m_data;

public:
  void set(const char * resn, const char * name1, const char * name2, int order) {
    m_data[resn].set(name1, name2, order);
  }

  const res_bond_dict_t * get(const char * resn) const {
    auto it = m_data.find(resn);
    if (it == m_data.end())
      return NULL;
    return &it->second;
  }

  bool empty() const {
    return m_data.empty();
  }
};

#endif
#endif
