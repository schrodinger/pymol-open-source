/*
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* copyright 1998-2000 by Warren Lyford Delano of DeLano Scientific.
D* -------------------------------------------------------------------
E* It is unlawful to modify or remove this copyright notice.
F* -------------------------------------------------------------------
G* Please see the accompanying LICENSE file for further information.
H* -------------------------------------------------------------------
I* Additional authors of this source file include:
-*
-*
-*
Z* -------------------------------------------------------------------
*/
#pragma once

#include "PyMOLGlobals.h"

#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <variant>

#include "PConv.h"
#include "OVLexicon.h"
#include "Result.h"

enum class PropertyType {
  Auto = -1,
  Blank = 0,
  Boolean = 1,
  Int = 2,
  Float = 3,
  // Unused = 4,
  Color = 5,
  String = 6,
};

void PropertyUniqueDetachChain(PyMOLGlobals* G, int index);

struct PropertyUniqueEntry {
  int prop_id;
  PropertyType type;
  std::variant<int, double> value;
  int next; /* for per- property lists & memory management */
};

struct CPropertyUnique {
  OVLexicon* propnames{};     /* Property Name Lookup */
  OVLexicon* string_values{}; /* String lookup */
  std::unordered_map<int, int> id2offset;
  std::vector<PropertyUniqueEntry> entry;
  int next_free{};

  int NextUniqueID = 1;
  std::unordered_set<int> ActiveIDs{};
  constexpr static int InitNumEntries = 10;
};
int PropertyGetNewUniqueID(PyMOLGlobals* G);

void PropertySetromString(
    PyMOLGlobals* G, int unique_id, const char* prop_name, const char* value);
void PropertySet(
    PyMOLGlobals* G, int unique_id, const char* prop_name, const char* value);
void PropertyUniqueSet_b(
    PyMOLGlobals* G, int unique_id, const char* prop_name, int value);
void PropertySet(
    PyMOLGlobals* G, int unique_id, const char* prop_name, int value);
void PropertySet(
    PyMOLGlobals* G, int unique_id, const char* prop_name, double value);
void PropertyUniqueSet_color(
    PyMOLGlobals* G, int unique_id, const char* prop_name, int value);
int PropertyUniqueSetTypedValue(PyMOLGlobals* G, int unique_id,
    const char* prop_name, PropertyType prop_type, const void* value);

void PropertyInit(PyMOLGlobals* G);
void PropertyFree(PyMOLGlobals* G);

int PropertyCopyProperties(
    PyMOLGlobals* G, int src_unique_id, int dst_unique_id);

CPythonVal* PropertyGetPropertyImpl(
    PyMOLGlobals* G, int prop_id, const char* propname, short returnList);
int PropertySetPropertyImpl(PyMOLGlobals* G, int prop_id, const char* propname,
    CPythonVal* value, PropertyType proptype);

int PropertyFromPyList(PyMOLGlobals* G, PyObject* list);

int PropertyAddAllToDictItem(
    PyMOLGlobals* G, PyObject* settingdict, int src_unique_id);
PyObject* PropertyAsPyList(PyMOLGlobals* G, int unique_id, int include_type);
PyObject* PropertyGetNamesAsPyList(PyMOLGlobals* G, int unique_id);

float PropertyGetAsFloat(PyMOLGlobals* G, int prop_id, const char* prop_name);

struct PropertyFmtOptions {
  std::string_view decFmt = "%d";
  std::string_view boolFmt = "%d";
  std::string_view floatFmt = "%.5f";
};
const char* PropertyGetAsString(PyMOLGlobals* G, int prop_id, const char* prop_name,
    char* buffer, PropertyFmtOptions fmt = {});
pymol::Result<std::string> PropertyGetAsString(PyMOLGlobals* G, int prop_id,
    std::string_view prop_name, PropertyFmtOptions fmt = {});
PropertyUniqueEntry* PropertyFindPropertyUniqueEntry(
    PyMOLGlobals* G, int prop_id, int propname_id);

// -------------------------------------------------------------
// templated functions to handle any type with a "prop_id" field
// -------------------------------------------------------------

template <typename T> inline int PropertyCheckUniqueID(PyMOLGlobals* G, T* x)
{
  if (!x->prop_id)
    x->prop_id = PropertyGetNewUniqueID(G);
  return x->prop_id;
}

template <typename T, typename V>
inline int PropertySet(PyMOLGlobals* G, T* x, const char* propname, V value)
{
  PropertyCheckUniqueID(G, x);
  PropertySet(G, x->prop_id, propname, value);
  return true;
}

// -------------------------------------------------------------
// Python
// -------------------------------------------------------------

// derive property type from Python type
int PropertySet(PyMOLGlobals* G, int id, const char* propname, PyObject* value);

// handles auto and color
template <typename T>
inline int PropertySet(PyMOLGlobals* G, T* x, const char* propname,
    PyObject* value, PropertyType proptype)
{
  PropertyCheckUniqueID(G, x);
  return PropertySetPropertyImpl(G, x->prop_id, propname, value, proptype);
}

// handles three cases... (1) list of names (2) value (3) [type, value]
template <typename T>
PyObject* PropertyGetPyObject(
    PyMOLGlobals* G, T* x, const char* propname, bool returnList = true)
{
  PyObject* r = nullptr;
  if (x->prop_id)
    r = PropertyGetPropertyImpl(G, x->prop_id, propname, returnList);
  return r;
}
