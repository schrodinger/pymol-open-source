

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

#include "Property.h"

#include "os_std.h"

#include "Base.h"
#include "PConv.h"
#include "Err.h"
#include "Feedback.h"
#include "MemoryDebug.h"
#include "OVContext.h"
#include "Ortho.h"
#include "PConv.h"
#include "Util.h"
#include "strcasecmp.h"

static PropertyUniqueEntry* PropertyFindPropertyUniqueEntry(
    PyMOLGlobals* G, int prop_id, const char* propname);

/**
 * Delete all properties for unique_id (e.g. AtomInfoType::prop_id)
 */
void PropertyUniqueDetachChain(PyMOLGlobals* G, int unique_id)
{
  CPropertyUnique* I = G->PropertyUnique;

  auto offsetIt = I->id2offset.find(unique_id);
  if (offsetIt != I->id2offset.end()) {
    int offset = offsetIt->second;
    int next;
    I->id2offset.erase(offsetIt);
    {
      while (offset) {
        auto entry = &I->entry[offset];
        OVLexicon_DecRef(I->propnames, entry->prop_id);
        if (entry->type == PropertyType::String) {
          OVLexicon_DecRef(I->string_values, std::get<int>(entry->value));
        }
        next = entry->next;
        entry->next = I->next_free;
        I->next_free = offset;
        offset = next;
      }
    }
  } else {
    /* uncaught error */
  }
  I->ActiveIDs.erase(unique_id);
}

/**
 * Get an unused unique id and mark it used (in ActiveIDs)
 */
int PropertyGetNewUniqueID(PyMOLGlobals* G)
{
  CPropertyUnique* I = G->PropertyUnique;
  int result = 0;
  while (true) {
    result = I->NextUniqueID++;
    if (result) { /* skip zero */
      if (I->ActiveIDs.count(result) == 0) {
        I->ActiveIDs.emplace(result);
        break;
      }
    }
  }
  return result;
}

void PropertyInit(PyMOLGlobals* G)
{
  G->PropertyUnique = new CPropertyUnique();
  auto I = G->PropertyUnique;

  I->entry.resize(CPropertyUnique::InitNumEntries, PropertyUniqueEntry{});
  /* note: intentially skip index 0  */
  for (int a = 2; a < I->entry.size(); a++) {
    I->entry[a].next = a - 1; /* 1-based linked list with 0 as sentinel */
  }
  I->next_free = I->entry.size() - 1;
  OVContext* C = G->Context;
  I->propnames = OVLexicon_New(C->heap);
  I->string_values = OVLexicon_New(C->heap);
}

void PropertyFree(PyMOLGlobals* G)
{
  CPropertyUnique* I = G->PropertyUnique;
  OVLexicon_DEL_AUTO_NULL(I->propnames);
  OVLexicon_DEL_AUTO_NULL(I->string_values);
  DeleteP(G->PropertyUnique);
}

static void PropertyUniqueExpand(PyMOLGlobals* G)
{
  CPropertyUnique* I = G->PropertyUnique;

  if (!I->next_free) {
    int new_n_alloc = (I->entry.size() * 3) / 2;
    auto sizeBefore = I->entry.size();
    I->entry.resize(new_n_alloc, PropertyUniqueEntry{});
    for (int a = sizeBefore; a < I->entry.size(); a++) {
      I->entry[a].next = I->next_free;
      I->next_free = a;
    }
  }
}

/**
 * Guess type, in order: int, float, bool (true/yes/false/no), string
 */
void PropertySetromString(
    PyMOLGlobals* G, int unique_id, const char* prop_name, const char* value)
{
  char nextchar;
  union {
    double fval;
    int ival;
  };

  if (1 == sscanf(value, "%i %c", &ival, &nextchar)) {
    PropertySet(G, unique_id, prop_name, ival);
    return;
  }

  if (1 == sscanf(value, "%lf %c", &fval, &nextchar)) {
    PropertySet(G, unique_id, prop_name, fval);
    return;
  }

  if (strcasecmp(value, "true") == 0 || strcasecmp(value, "yes") == 0) {
    PropertyUniqueSet_b(G, unique_id, prop_name, 1);
    return;
  }

  if (strcasecmp(value, "false") == 0 || strcasecmp(value, "no") == 0) {
    PropertyUniqueSet_b(G, unique_id, prop_name, 0);
    return;
  }

  PropertySet(G, unique_id, prop_name, value);
}

void PropertySet(
    PyMOLGlobals* G, int unique_id, const char* prop_name, const char* value)
{
  int val;
  CPropertyUnique* I = G->PropertyUnique;
  OVreturn_word result;
  if (OVreturn_IS_OK(
          (result = OVLexicon_GetFromCString(I->string_values, value)))) {
    val = result.word;
    PropertyUniqueSetTypedValue(
        G, unique_id, prop_name, PropertyType::String, &val);
  }
}

void PropertyUniqueSet_b(
    PyMOLGlobals* G, int unique_id, const char* prop_name, int value)
{
  PropertyUniqueSetTypedValue(
      G, unique_id, prop_name, PropertyType::Boolean, &value);
}

void PropertySet(
    PyMOLGlobals* G, int unique_id, const char* prop_name, int value)
{
  PropertyUniqueSetTypedValue(
      G, unique_id, prop_name, PropertyType::Int, &value);
}

void PropertySet(
    PyMOLGlobals* G, int unique_id, const char* prop_name, double value)
{
  PropertyUniqueSetTypedValue(
      G, unique_id, prop_name, PropertyType::Float, &value);
}

void PropertyUniqueSet_color(
    PyMOLGlobals* G, int unique_id, const char* prop_name, int value)
{
  PropertyUniqueSetTypedValue(
      G, unique_id, prop_name, PropertyType::Color, &value);
}

int PropertyUniqueSetTypedValue(PyMOLGlobals* G, int unique_id,
    const char* prop_name, PropertyType prop_type, const void* value)
{
  CPropertyUnique* I = G->PropertyUnique;

  int isset = false;
  int prop_id = OVLexicon_GetFromCString(I->propnames, prop_name).word;
  auto offsetIt = I->id2offset.find(unique_id);
  if (offsetIt != I->id2offset.end()) { /* setting list exists for atom */
    int offset = offsetIt->second;
    int prev = 0;
    int found = false;
    while (offset) {
      auto entry = &I->entry[offset];
      if (entry->prop_id == prop_id) {
        found = true; /* this setting is already defined */
        if (value) {  /* if redefining value */
          switch (prop_type) {
          case PropertyType::Float:
            if (std::get<double>(entry->value) != *(double*) value ||
                entry->type != prop_type) {
              if (entry->type == PropertyType::String &&
                  std::holds_alternative<int>(entry->value)) {
                OVLexicon_DecRef(I->string_values, std::get<int>(entry->value));
              }
              entry->value = *(double*) value;
              entry->type = prop_type;
              isset = true;
            }
            break;
          default:
            if (std::get<int>(entry->value) != *(int*) value ||
                entry->type != prop_type) {
              if (entry->type == PropertyType::String &&
                  std::holds_alternative<int>(entry->value)) {
                OVLexicon_DecRef(I->string_values, std::get<int>(entry->value));
              }
              entry->value = *(int*) value;
              entry->type = prop_type;
              isset = true;
            }
          }
        } else {       /* or nullptr value means delete this setting */
          if (!prev) { /* if first entry in list */
            I->id2offset.erase(offsetIt);
            if (entry->next) { /* set new list start */
              I->id2offset[unique_id] = entry->next;
            }
          } else { /* otherwise excise from middle or end */
            I->entry[prev].next = entry->next;
          }
          entry->next = I->next_free;
          I->next_free = offset;
          isset = true;
        }
        break;
      }
      prev = offset;
      offset = entry->next;
    }
    if (found) {
      OVLexicon_DecRef(I->propnames, prop_id);
    }
    if ((!found) &&
        value) { /* setting not found in existing list, so append new value */
      if (!I->next_free)
        PropertyUniqueExpand(G);
      if (I->next_free) {
        offset = I->next_free;
        {
          auto entry = &I->entry[offset];
          short should_set = false;
          I->next_free = entry->next;
          entry->next = 0;

          if (prev) { /* append onto existing list */
            I->entry[prev].next = offset;
            should_set = true;
          } else {
            I->id2offset[unique_id] = offset;
            should_set = true;
          }
          if (should_set) {
            entry->type = prop_type;
            switch (prop_type) {
            case PropertyType::Float:
              entry->value = *(double*) value;
              break;
            default:
              entry->value = *(int*) value;
            }
            entry->prop_id = prop_id;
            isset = true;
          }
        }
      }
    }
  } else if (value) { /* new setting list for atom */
    if (!I->next_free)
      PropertyUniqueExpand(G);
    if (I->next_free) {
      int offset = I->next_free;
      auto entry = &I->entry[offset];

      I->id2offset[unique_id] = offset;
      I->next_free = entry->next;
      entry->type = prop_type;
      switch (prop_type) {
      case PropertyType::Float:
        entry->value = *(double*) value;
        break;
      default:
        entry->value = *(int*) value;
      }
      entry->prop_id = prop_id;
      entry->next = 0;
      isset = true;
    }
  }
  return isset;
}

/**
 * For copying atoms, coordsets, etc...
 */
int PropertyCopyProperties(
    PyMOLGlobals* G, int src_unique_id, int dst_unique_id)
{
  int ok = true;
  CPropertyUnique* I = G->PropertyUnique;

  auto dstOffsetIt = I->id2offset.find(dst_unique_id);
  if (dstOffsetIt != I->id2offset.end()) { /* setting list exists for atom */
    /* note, this code path not yet tested...doesn't occur */

    auto srcOffsetIt = I->id2offset.find(src_unique_id);
    if (srcOffsetIt != I->id2offset.end()) {
      int src_offset = srcOffsetIt->second;
      int dst_offset = dstOffsetIt->second;
      PropertyUniqueEntry* src_entry;
      while (src_offset) {
        src_entry = &I->entry[src_offset];

        {
          int prop_id = src_entry->prop_id;
          auto property_type = src_entry->type;
          PropertyUniqueEntry* property_value = src_entry;
          int prev = 0;
          int found = false;
          while (dst_offset) {
            auto dst_entry = &I->entry[dst_offset];
            if (dst_entry->prop_id == prop_id) {
              found = true; /* this setting is already defined */
              dst_entry->value = property_value->value;
              dst_entry->type = property_type;
              break;
            }
            prev = dst_offset;
            dst_offset = dst_entry->next;
          }
          if (!found) { /* setting not found in existing list, so append new
                           value */
            if (!I->next_free)
              PropertyUniqueExpand(G);
            if (I->next_free) {
              dst_offset = I->next_free;
              {
                auto dst_entry = &I->entry[dst_offset];
                I->next_free = dst_entry->next;
                dst_entry->next = 0;
                if (prev) { /* append onto existing list */
                  I->entry[prev].next = dst_offset;
                  dst_entry->type = property_type;
                  dst_entry->value = property_value->value;
                  dst_entry->prop_id = prop_id;
                  // if copying, need to increment refs in lexicons
                  OVLexicon_IncRef(I->propnames, dst_entry->prop_id);
                  if (dst_entry->type == PropertyType::String) {
                    OVLexicon_IncRef(
                        I->string_values, std::get<int>(dst_entry->value));
                  }
                } else {
                  I->id2offset[dst_unique_id] = dst_offset;
                  /* create new list */
                  dst_entry->type = property_type;
                  dst_entry->value = property_value->value;
                  dst_entry->prop_id = prop_id;
                  // if copying, need to increment refs in lexicons
                  OVLexicon_IncRef(I->propnames, dst_entry->prop_id);
                  if (dst_entry->type == PropertyType::String) {
                    OVLexicon_IncRef(
                        I->string_values, std::get<int>(dst_entry->value));
                  }
                }
              }
            }
          }
        }
        src_offset =
            I->entry[src_offset]
                .next; /* src_entry invalid, since I->entry may have changed */
      }
    }
  } else { /* new setting list for atom */
    auto srcOffsetIt = I->id2offset.find(src_unique_id);
    if (srcOffsetIt != I->id2offset.end()) {
      int src_offset = srcOffsetIt->second;
      int prev = 0;
      PropertyUniqueEntry* src_entry;
      while (ok && src_offset) {
        if (!I->next_free)
          PropertyUniqueExpand(G);
        {
          src_entry = &I->entry[src_offset];
          {
            int prop_id = src_entry->prop_id;
            auto property_type = src_entry->type;
            PropertyUniqueEntry* property_value = src_entry;
            if (I->next_free) {
              int dst_offset = I->next_free;
              auto dst_entry = &I->entry[dst_offset];
              I->next_free = dst_entry->next;

              if (!prev) {
                I->id2offset[dst_unique_id] = dst_offset;
              } else {
                I->entry[prev].next = dst_offset;
              }

              if (ok) {
                dst_entry->type = property_type;
                dst_entry->value = property_value->value;
                dst_entry->prop_id = prop_id;
                dst_entry->next = 0;
                OVLexicon_IncRef(I->propnames, dst_entry->prop_id);
                if (dst_entry->type == PropertyType::String) {
                  OVLexicon_IncRef(
                      I->string_values, std::get<int>(dst_entry->value));
                }
              }
              prev = dst_offset;
            }
          }
        }
        src_offset =
            I->entry[src_offset]
                .next; /* src_entry invalid, since I->entry may have changed */
      }
    }
  }

  return ok;
}

// -------------------------------------------------------------------------------------
// //

/**
 * Return a list of all of the property names for this object
 */
PyObject* PropertyGetNamesAsPyList(PyMOLGlobals* G, int unique_id)
{
  CPropertyUnique* I = G->PropertyUnique;
  PyObject* list = PyList_New(0);

  auto offsetIt = I->id2offset.find(unique_id);
  if (offsetIt != I->id2offset.end()) {
    int src_offset = offsetIt->second;
    auto src_entry = &I->entry[src_offset];
    while (src_offset) {
      char* s = OVLexicon_FetchCString(I->propnames, src_entry->prop_id);
      PyObject* item = PyString_FromString(s);
      PyList_Append(list, item);
      Py_DECREF(item);
      src_offset = src_entry->next;
      src_entry = &I->entry[src_offset];
    }
  }

  return list;
}

/**
 * Handles three cases:
 * 1. prop_name == nullptr:       return list of all keys
 * 2. returnList = false:      return value
 * 3. returnList = true:       return [type, value]
 */
CPythonVal* PropertyGetPropertyImpl(
    PyMOLGlobals* G, int prop_id, const char* prop_name, short returnList)
{
  CPropertyUnique* I = G->PropertyUnique;
  CPythonVal* retVal = nullptr;
  if (!prop_name) {
    retVal = PropertyGetNamesAsPyList(G, prop_id);
  } else {
    PropertyUniqueEntry* src_entry =
        PropertyFindPropertyUniqueEntry(G, prop_id, prop_name);
    if (src_entry) {
      //	int prop_id = src_entry->prop_id;
      auto prop_type = src_entry->type;
      CPythonVal* vVal = nullptr;
      switch (prop_type) {
      case PropertyType::Color:
      case PropertyType::Int:
        vVal = CPythonVal_New_Integer(std::get<int>(src_entry->value));
        break;
      case PropertyType::Boolean:
        vVal = CPythonVal_New_Boolean(std::get<int>(src_entry->value));
        break;
      case PropertyType::Float:
        vVal = CPythonVal_New_Float(std::get<double>(src_entry->value));
        break;
      case PropertyType::String: {
        const char* prop_value = OVLexicon_FetchCString(
            I->string_values, std::get<int>(src_entry->value));
        vVal = CPythonVal_New_String(prop_value, strlen(prop_value));
      } break;
      }
      if (returnList && vVal) {
        retVal = PyList_New(2);
        PyList_SetItem(retVal, 0, PyInt_FromLong(static_cast<int>(prop_type)));
        PyList_SetItem(retVal, 1, vVal); // steal vVal ref
      } else {
        retVal = vVal;
      }
    }
  }
  return retVal;
}

/**
 * Get a property as float. All numerical types will be casted to float,
 * string type will be parsed as float. Returns 0.0 on failure.
 */
float PropertyGetAsFloat(PyMOLGlobals* G, int prop_id, const char* prop_name)
{
  float fval;
  char* sval;
  PropertyUniqueEntry* entry =
      PropertyFindPropertyUniqueEntry(G, prop_id, prop_name);
  ok_assert(1, entry);
  switch (entry->type) {
  case PropertyType::Float:
    return std::get<double>(entry->value);
  case PropertyType::Int:
  case PropertyType::Color:
    return static_cast<float>(std::get<int>(entry->value));
  case PropertyType::Boolean:
    return std::get<int>(entry->value) ? 1.0 : 0.0;
  case PropertyType::String:
    sval = OVLexicon_FetchCString(
        G->PropertyUnique->string_values, std::get<int>(entry->value));
    ok_assert(1, sval && sscanf(sval, "%f", &fval));
    return fval;
  }
ok_except1:
  return 0.0;
}

/**
 * Get a property as string. For string type, return a borrowed reference to
 * the dictionary lookup. For numerical types, write the string representation
 * into the "buffer" pointer and return that. Return nullptr on failure.
 *
 * @param int atom property ID
 * @param prop_name name of property
 * @param buffer buffer to store property value if valid
 * @param fmt optional format for buffer output
 * @return buffer is valid; nullptr if fail.
 */
const char* PropertyGetAsString(PyMOLGlobals* G, int prop_id, const char* prop_name,
    char* buffer, PropertyFmtOptions fmt)
{
  PropertyUniqueEntry* entry =
      PropertyFindPropertyUniqueEntry(G, prop_id, prop_name);
  ok_assert(1, entry);
  switch (entry->type) {
  case PropertyType::Float:
    snprintf(buffer, 64, fmt.floatFmt.data(), std::get<double>(entry->value));
    return buffer;
  case PropertyType::Int:
  case PropertyType::Color:
    snprintf(buffer, 64, fmt.decFmt.data(), std::get<int>(entry->value));
    return buffer;
  case PropertyType::Boolean:
    sprintf(buffer, fmt.boolFmt.data(), std::get<int>(entry->value) ? 1 : 0);
    return buffer;
  case PropertyType::String:
    return OVLexicon_FetchCString(
        G->PropertyUnique->string_values, std::get<int>(entry->value));
  }
ok_except1:
  return nullptr;
}

/**
 * Gets the property as string.
 *
 * @param prop_id atom property id
 * @param prop_name property name
 * @param fmt optional format for buffer output
 * @return property as std::string
 */
pymol::Result<std::string> PropertyGetAsString(PyMOLGlobals* G, int prop_id,
    std::string_view prop_name, PropertyFmtOptions fmt)
{
  OrthoLineType buffer{};
  auto res = PropertyGetAsString(G, prop_id, prop_name.data(), buffer, fmt);
  if (!res) {
    return pymol::make_error("Invalid property");
  }
  return std::string(res);
}

/**
 * Get property entry handle
 *
 * @param[in] prop_id       Properties handle id (e.g. for an object or atom)
 * @param[in] propname_id   Property name id
 * @return  Pointer to property item or nullptr if not found
 */
PropertyUniqueEntry* PropertyFindPropertyUniqueEntry(
    PyMOLGlobals* G, int prop_id, int propname_id)
{
  auto I = G->PropertyUnique;
  PropertyUniqueEntry* item;

  auto offsetIt = I->id2offset.find(prop_id);
  if (offsetIt != I->id2offset.end()) {
    for (auto offset = offsetIt->second; offset; offset = item->next) {
      item = &I->entry[offset];
      if (item->prop_id == propname_id)
        return item;
    }
  }

  return nullptr;
}

/**
 * Get property entry handle
 */
static PropertyUniqueEntry* PropertyFindPropertyUniqueEntry(
    PyMOLGlobals* G, int prop_id, const char* propname)
{
  CPropertyUnique* I = G->PropertyUnique;
  OVreturn_word result;
  if (OVreturn_IS_OK(
          (result = OVLexicon_BorrowFromCString(I->propnames, propname)))) {
    int prop_id_tmp = result.word;
    return PropertyFindPropertyUniqueEntry(G, prop_id, prop_id_tmp);
  }
  return nullptr;
}

/**
 * Called from Python (pymol.properties.(set_property|set_atom_property)).
 * Types must be pre-handled with pymol.properties._typecast and are not
 * checked here.
 */
int PropertySetPropertyImpl(PyMOLGlobals* G, int prop_id, const char* propname,
    CPythonVal* value, PropertyType proptype)
{
#ifdef _PYMOL_NOPY
  return 0;
#else
  int keylen = strlen(propname);
  if (keylen > 1024) {
    PRINTFB(G, FB_Property, FB_Errors)
    " Property-Error: Name too long with %d characters (max 1024): "
    "'%.1023s...'\n",
        keylen, propname ENDFB(G);
    return 0;
  }

  union {
    int valueInt;
    const char* valueStr;
  };

  if (proptype == PropertyType::Auto && PConvPyStrToStrPtr(value, &valueStr)) {
    PropertySetromString(G, prop_id, propname, valueStr);
    return true;
  } else if (proptype == PropertyType::Color) {
    if (!PConvPyIntToInt(value, &valueInt)) {
      PRINTFB(G, FB_Property, FB_Errors)
      " Property-Error: wrong Python type for color\n" ENDFB(G);
      return false;
    }
    PropertyUniqueSet_color(G, prop_id, propname, valueInt);
    return true;
  }

  return PropertySet(G, prop_id, propname, value);
#endif
}

/**
 * For PSE loading
 */
int PropertyFromPyList(PyMOLGlobals* G, PyObject* list)
{
  int ok = true;
  int prop_id = 0;
  ov_size size;
  ov_size a;
  if (ok)
    ok = (list != nullptr);
  if (ok)
    ok = PyList_Check(list);
  if (ok) {
    size = PyList_Size(list);
    for (a = 0; a < size; a++) {
      CPythonVal* val = CPythonVal_PyList_GetItem(G, list, a);
      if (PyList_Check(val)) {
        int sizei = PyList_Size(val);
        if (sizei == 3) {
          char* prop_name;
          int prop_typeInt;
          CPythonVal* propval;
          int prop_name_len;
          CPythonVal* pVal = nullptr;
          if (!prop_id)
            prop_id = PropertyGetNewUniqueID(G);
          pVal = CPythonVal_PyList_GetItem(G, val, 0);
          prop_name_len = PyString_AsSomeString(pVal).length() + 1;
          prop_name = pymol::malloc<char>(prop_name_len);
          CPythonVal_PConvPyStrToStr(pVal, prop_name, prop_name_len);
          CPythonVal_Free(pVal);
          CPythonVal_PConvPyIntToInt_From_List(G, val, 1, &prop_typeInt);
          auto prop_type = static_cast<PropertyType>(prop_typeInt);
          propval = CPythonVal_PyList_GetItem(G, val, 2);
          switch (prop_type) {
          case PropertyType::Int:
          case PropertyType::Color: {
            int prop_value;
            CPythonVal_PConvPyIntToInt(propval, &prop_value);
            PropertySet(G, prop_id, prop_name, prop_value);
          } break;
          case PropertyType::Boolean: {
            int prop_value;
            CPythonVal_PConvPyBoolToInt(propval, &prop_value);
            PropertyUniqueSet_b(G, prop_id, prop_name, prop_value);
          } break;
          case PropertyType::Float: {
            double prop_value;
            CPythonVal_PConvPyFloatToDouble(propval, &prop_value);
            PropertySet(G, prop_id, prop_name, prop_value);
          } break;
          case PropertyType::String: {
            char* prop_value;
            int prop_value_len = PyString_AsSomeString(propval).length() + 1;
            prop_value = pymol::malloc<char>(prop_value_len);
            CPythonVal_PConvPyStrToStr(propval, prop_value, prop_value_len);
            PropertySet(G, prop_id, prop_name, prop_value);
            FreeP(prop_value);
          } break;
          }
          CPythonVal_Free(propval);
          FreeP(prop_name);
        }
      }
      CPythonVal_Free(val);
    }
  }
  return prop_id;
}

/**
 * For iterate: p.all
 */
int PropertyAddAllToDictItem(PyMOLGlobals* G, PyObject* propdict, int unique_id)
{
  int ok = true;
  CPropertyUnique* I = G->PropertyUnique;
  auto srcOffsetIt = I->id2offset.find(unique_id);
  if (srcOffsetIt != I->id2offset.end()) {
    int src_offset = srcOffsetIt->second;
    PropertyUniqueEntry* src_entry;
    while (ok && src_offset) {
      {
        src_entry = &I->entry[src_offset];
        {
          int prop_id = src_entry->prop_id;
          auto prop_type = src_entry->type;
          switch (prop_type) {
          case PropertyType::Int:
          case PropertyType::Color:
          case PropertyType::Boolean:
          case PropertyType::Float:
          case PropertyType::String:
            const char* prop_name =
                OVLexicon_FetchCString(I->propnames, prop_id);
            switch (prop_type) {
            case PropertyType::Int:
            case PropertyType::Color:
              PConvIntToPyDictItem(
                  propdict, prop_name, std::get<int>(src_entry->value));
              break;
            case PropertyType::Boolean:
              PyDict_SetItemString(propdict, prop_name,
                  std::get<int>(src_entry->value) ? Py_True : Py_False);
              break;
            case PropertyType::Float:
              PyDict_SetItemString(propdict, prop_name,
                  PyFloat_FromDouble(std::get<double>(src_entry->value)));
              break;
            case PropertyType::String: {
              const char* prop_value = OVLexicon_FetchCString(
                  I->string_values, std::get<int>(src_entry->value));
              PyDict_SetItemString(
                  propdict, prop_name, PyString_FromString(prop_value));
            } break;
            }
          }
        }
      }
      src_offset =
          I->entry[src_offset]
              .next; /* src_entry invalid, since I->entry may have changed */
    }
  }
  return ok;
}

/**
 * For PSE and ChemPy export
 */
PyObject* PropertyAsPyList(PyMOLGlobals* G, int unique_id, int include_type)
{
  int ok = true;
  CPropertyUnique* I = G->PropertyUnique;
  PyObject* retList = PyList_New(0);
  auto srcOffsetIt = I->id2offset.find(unique_id);
  if (srcOffsetIt != I->id2offset.end()) {
    int src_offset = srcOffsetIt->second;
    PropertyUniqueEntry* src_entry;
    while (ok && src_offset) {
      {
        src_entry = &I->entry[src_offset];
        {
          int prop_id = src_entry->prop_id;
          auto prop_type = src_entry->type;
          PyObject* pyvalue = nullptr;

          switch (prop_type) {
          case PropertyType::Int:
          case PropertyType::Color:
            pyvalue = PyInt_FromLong(std::get<int>(src_entry->value));
            break;
          case PropertyType::Boolean:
            pyvalue = std::holds_alternative<int>(src_entry->value) ? Py_True
                                                                    : Py_False;
            Py_INCREF(pyvalue);
            break;
          case PropertyType::Float:
            pyvalue = PyFloat_FromDouble(std::get<double>(src_entry->value));
            break;
          case PropertyType::String: {
            const char* prop_value = OVLexicon_FetchCString(
                I->string_values, std::get<int>(src_entry->value));
            pyvalue = PyString_FromString(prop_value);
          } break;
          default:
            printf(" %s-Error: unexpected type %d\n", __FUNCTION__,
                static_cast<int>(prop_type));
          }

          if (pyvalue) {
            int N = include_type ? 3 : 2;
            PyObject* propList = PyList_New(N);

            const char* prop_name =
                OVLexicon_FetchCString(I->propnames, prop_id);
            PyList_SetItem(propList, 0, PyString_FromString(prop_name));
            PyList_SetItem(propList, N - 1, pyvalue); // steal pyvalue ref

            if (include_type) {
              PyList_SetItem(
                  propList, 1, PyInt_FromLong(static_cast<int>(prop_type)));
            }

            PyList_Append(retList, propList);
            Py_DECREF(propList);
          }
        }
      }
      src_offset =
          I->entry[src_offset]
              .next; /* src_entry invalid, since I->entry may have changed */
    }
  }
  return retList;
}

#ifndef _PYMOL_NOPY
/**
 * Derive property type from Python type
 */
int PropertySet(
    PyMOLGlobals* G, int prop_id, const char* propname, PyObject* value)
{
  union {
    int valueInt;
    double valueFloat;
    const char* valueStr;
  };

  if (PConvPyFloatToDouble(value, &valueFloat)) {
    PropertySet(G, prop_id, propname, valueFloat);
  } else if (PConvPyBoolToInt(value, &valueInt)) {
    PropertyUniqueSet_b(G, prop_id, propname, valueInt);
  } else if (PConvPyIntToInt(value, &valueInt)) {
    PropertySet(G, prop_id, propname, valueInt);
  } else if (PConvPyStrToStrPtr(value, &valueStr)) {
    PropertySet(G, prop_id, propname, valueStr);
  } else if (!value || value == Py_None) {
    // PropertyDel
    PropertyUniqueSetTypedValue(
        G, prop_id, propname, PropertyType::Blank, nullptr);
  } else {
    PyErr_Format(PyExc_TypeError,
        "Unsupported type for property '%s', "
        "only supports bool, int, float, str",
        propname);
    return false;
  }

  return true;
}
#endif
