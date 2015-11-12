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

#include "os_std.h"
#include "PyMOLGlobals.h"

#if 1
// (2)
// Optimized implementation, tuned on speed and memory usage.
// Only supports resn length <= 8 and atom name length <= 4.
// Working with this limitation is currently reasonable in PyMOL,
// since cResnLen=5 and cAtomNameLen=4 (see AtomInfo.h).

#ifdef _PYMOL_NO_CXX11
#define unordered_map map
#else
#include <unordered_map>
#endif

#include <stdint.h>
#include <string.h>

// mapped type of bond_dict_t
class res_bond_dict_t : std::unordered_map<int64_t, signed char> {
  static key_type make_key(const char * name1, const char * name2) {
    union { char s[4]; int32_t i; } u1, u2;

    strncpy(u1.s, name1, 4);
    strncpy(u2.s, name2, 4);
    if (u1.i > u2.i)
      std::swap(u1.i, u2.i);

    // construct 8-byte key from the two 4-byte names
    return (((key_type)u1.i) << 32) | u2.i;
  }

public:
  void set(const char * name1, const char * name2, mapped_type order) {
    (*this)[make_key(name1, name2)] = order;
  }

  mapped_type get(const char * name1, const char * name2) const {
    auto it = find(make_key(name1, name2));
    if (it == end())
      return -1;
    return it->second;
  }
};

// type for mapping: resn -> ((name1, name2) -> bond order)
class bond_dict_t : public std::map<int64_t, res_bond_dict_t> {
  static key_type make_key(const char * s_) {
    union { key_type i; char s[sizeof(key_type)]; };
    strncpy(s, s_, sizeof(s));
    return i;
  }

  std::set<key_type> unknown_resn;

public:
  void set_unknown(const char * resn) {
    unknown_resn.insert(make_key(resn));
  }

  void set(const char * resn, const char * name1, const char * name2, int order) {
    (*this)[make_key(resn)].set(name1, name2, order);
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
