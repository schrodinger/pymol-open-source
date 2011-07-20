
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
#include"os_python.h"

#include"os_proprietary.h"

#include"os_predef.h"
#include"os_std.h"

#include"Base.h"
#include"OOMac.h"
#include"MemoryDebug.h"
#include"Ortho.h"
#include"Setting.h"
#include"Scene.h"
#include"ButMode.h"
#include"Executive.h"
#include"Editor.h"
#include"P.h"
#include"Util.h"
#include"main.h"
#include"PConv.h"
#include"Wizard.h"
#include"Seq.h"
#include"PyMOLOptions.h"
#include"OVContext.h"

static void *SettingPtr(CSetting * I, int index, ov_size size);

static CSetting *SettingCopyAll(PyMOLGlobals * G, CSetting * src, CSetting * dst)
{
  if(!dst) {
    dst = Calloc(CSetting, 1);
    if(dst) {
      SettingInit(G, dst);
    }
  }
  if(dst && src) {

    /* simply overwriting existing data (if any) ... in the future we
       may need to release references etc. before doing this */

    unsigned int size = VLAGetSize(src->info);
    VLACheck(dst->info, SettingRec, size);
    UtilCopyMem(dst->info, src->info, sizeof(SettingRec) * size);
    VLACheck(dst->data, char, src->size);
    dst->size = src->size;
    UtilCopyMem(dst->data, src->data, src->size);
  }
  return dst;
}

void SettingStoreDefault(PyMOLGlobals * G)
{
  G->Default = SettingCopyAll(G, G->Setting, G->Default);
}

void SettingPurgeDefault(PyMOLGlobals * G)
{
  if(G->Default) {
    SettingPurge(G->Default);
    FreeP(G->Default);
    G->Default = NULL;
  }
}

void SettingUniqueDetachChain(PyMOLGlobals * G, int unique_id)
{
  register CSettingUnique *I = G->SettingUnique;
  OVreturn_word result;
  if(OVreturn_IS_OK(result = OVOneToOne_GetForward(I->id2offset, unique_id))) {
    int offset = result.word;
    int next;

    OVOneToOne_DelForward(I->id2offset, unique_id);

    {
      SettingUniqueEntry *entry;
      while(offset) {
        entry = I->entry + offset;
        next = entry->next;
        entry->next = I->next_free;
        I->next_free = offset;
        offset = next;
      }
    }
  } else {
    /* uncaught error */
  }
}

static void SettingUniqueExpand(PyMOLGlobals * G)
{
  register CSettingUnique *I = G->SettingUnique;

  if(!I->next_free) {
    int new_n_alloc = (I->n_alloc * 3) / 2;
    int a;
    VLACheck(I->entry, SettingUniqueEntry, new_n_alloc);
    for(a = I->n_alloc; a < new_n_alloc; a++) {
      I->entry[a].next = I->next_free;
      I->next_free = a;
    }
    I->n_alloc = new_n_alloc;
  }
}

int SettingUniqueCheck(PyMOLGlobals * G, int unique_id, int setting_id)
{
  register CSettingUnique *I = G->SettingUnique;
  OVreturn_word result;
  if(OVreturn_IS_OK(result = OVOneToOne_GetForward(I->id2offset, unique_id))) {
    int offset = result.word;
    SettingUniqueEntry *entry;

    while(offset) {
      entry = I->entry + offset;
      if(entry->setting_id == setting_id) {
        return 1;
      }
      offset = entry->next;
    }
  }
  return 0;
}

static int SettingUniqueGetTypedValue(PyMOLGlobals * G, int unique_id, int setting_id,
                                      int setting_type, void *value)
{
  register CSettingUnique *I = G->SettingUnique;
  OVreturn_word result;
  if(OVreturn_IS_OK(result = OVOneToOne_GetForward(I->id2offset, unique_id))) {
    int offset = result.word;
    SettingUniqueEntry *entry;
    while(offset) {
      entry = I->entry + offset;
      if(entry->setting_id == setting_id) {
        if(entry->type == setting_type) {
          *(int *) value = entry->value.int_;
        } else
          switch (setting_type) {
          case cSetting_int:
          case cSetting_color:
          case cSetting_boolean:
            switch (entry->type) {
            case cSetting_float:
              *(int *) value = (int) (*(float *) &entry->value.float_);
              break;
            default:
              *(int *) value = entry->value.int_;
              break;
            }
            break;
          case cSetting_float:
            *(float *) value = (float) (*(int *) &entry->value.int_);
            break;
          }
        return 1;
      }
      offset = entry->next;
    }
  }
  return 0;
}

int SettingUniqueGet_b(PyMOLGlobals * G, int unique_id, int setting_id, int *value)
{
  return SettingUniqueGetTypedValue(G, unique_id, setting_id, cSetting_boolean, value);
}

int SettingUniqueGet_i(PyMOLGlobals * G, int unique_id, int setting_id, int *value)
{
  return SettingUniqueGetTypedValue(G, unique_id, setting_id, cSetting_int, value);
}

int SettingUniqueGet_f(PyMOLGlobals * G, int unique_id, int setting_id, float *value)
{
  return SettingUniqueGetTypedValue(G, unique_id, setting_id, cSetting_float, value);
}

int SettingUniqueGet_color(PyMOLGlobals * G, int unique_id, int setting_id, int *value)
{
  return SettingUniqueGetTypedValue(G, unique_id, setting_id, cSetting_color, value);
}

void SettingUniqueSetTypedValue(PyMOLGlobals * G, int unique_id, int setting_id,
                                int setting_type, void *value)

/* set value to NULL in order to delete setting */
{
  register CSettingUnique *I = G->SettingUnique;
  OVreturn_word result;

  if(OVreturn_IS_OK((result = OVOneToOne_GetForward(I->id2offset, unique_id)))) {       /* setting list exists for atom */
    int offset = result.word;
    int prev = 0;
    int found = false;
    while(offset) {
      SettingUniqueEntry *entry = I->entry + offset;
      if(entry->setting_id == setting_id) {
        found = true;           /* this setting is already defined */
        if(value) {             /* if redefining value */
          entry->value.int_ = *(int *) value;
          entry->type = setting_type;
        } else {                /* or NULL value means delete this setting */
          if(!prev) {           /* if first entry in list */
            OVOneToOne_DelForward(I->id2offset, unique_id);
            if(entry->next) {   /* set new list start */
              OVOneToOne_Set(I->id2offset, unique_id, entry->next);
            }
          } else {              /* otherwise excise from middle or end */
            I->entry[prev].next = entry->next;
          }
          entry->next = I->next_free;
          I->next_free = offset;
        }
        break;
      }
      prev = offset;
      offset = entry->next;
    }
    if((!found) && value) {     /* setting not found in existing list, so append new value */
      if(!I->next_free)
        SettingUniqueExpand(G);
      if(I->next_free) {
        offset = I->next_free;
        {
          SettingUniqueEntry *entry = I->entry + offset;
          I->next_free = entry->next;
          entry->next = 0;

          if(prev) {            /* append onto existing list */
            I->entry[prev].next = offset;
            entry->type = setting_type;
            entry->value.int_ = *(int *) value;
            entry->setting_id = setting_id;
          } else if(OVreturn_IS_OK(OVOneToOne_Set(I->id2offset, unique_id, offset))) {
            /* create new list */
            entry->type = setting_type;
            entry->value.int_ = *(int *) value;
            entry->setting_id = setting_id;
          }
        }
      }
    }
  } else if(value && (result.status == OVstatus_NOT_FOUND)) {   /* new setting list for atom */
    if(!I->next_free)
      SettingUniqueExpand(G);
    if(I->next_free) {
      int offset = I->next_free;
      SettingUniqueEntry *entry = I->entry + offset;

      if(OVreturn_IS_OK(OVOneToOne_Set(I->id2offset, unique_id, offset))) {
        I->next_free = entry->next;
        entry->type = setting_type;
        entry->value.int_ = *(int *) value;
        entry->setting_id = setting_id;
        entry->next = 0;
      }
    }
  } else {
    /* unhandled error */
  }
}

void SettingUniqueSet_i(PyMOLGlobals * G, int unique_id, int setting_id, int value)
{
  SettingUniqueSetTypedValue(G, unique_id, setting_id, cSetting_int, &value);
}

void SettingUniqueSet_f(PyMOLGlobals * G, int unique_id, int setting_id, float value)
{
  SettingUniqueSetTypedValue(G, unique_id, setting_id, cSetting_float, &value);
}

void SettingUniqueSet_b(PyMOLGlobals * G, int unique_id, int setting_id, int value)
{
  SettingUniqueSetTypedValue(G, unique_id, setting_id, cSetting_boolean, &value);
}

void SettingUniqueSet_color(PyMOLGlobals * G, int unique_id, int setting_id, int value)
{
  SettingUniqueSetTypedValue(G, unique_id, setting_id, cSetting_color, &value);
}

void SettingUniqueResetAll(PyMOLGlobals * G)
{
  register CSettingUnique *I = G->SettingUnique;

  OVOneToOne_Reset(I->id2offset);
  {
    int a;
    I->n_alloc = 10;
    VLAFreeP(I->entry);
    I->entry = VLACalloc(SettingUniqueEntry, I->n_alloc);
    /* note: intentially skip index 0  */
    for(a = 2; a < 10; a++) {
      I->entry[a].next = a - 1;
    }
    I->next_free = I->n_alloc - 1;
  }
}

int SettingUniqueCopyAll(PyMOLGlobals * G, int src_unique_id, int dst_unique_id)
{
  int ok = true;
  register CSettingUnique *I = G->SettingUnique;
  OVreturn_word dst_result;

  if(OVreturn_IS_OK((dst_result = OVOneToOne_GetForward(I->id2offset, dst_unique_id)))) {       /* setting list exists for atom */

    /* note, this code path not yet tested...doesn't occur */

    OVreturn_word src_result;
    if(OVreturn_IS_OK(src_result = OVOneToOne_GetForward(I->id2offset, src_unique_id))) {
      int src_offset = src_result.word;
      SettingUniqueEntry *src_entry;
      while(src_offset) {
        src_entry = I->entry + src_offset;

        {
          int setting_id = src_entry->setting_id;
          int setting_type = src_entry->type;
          int setting_value = src_entry->value.int_;
          int dst_offset = dst_result.word;
          int prev = 0;
          int found = false;
          while(dst_offset) {
            SettingUniqueEntry *dst_entry = I->entry + dst_offset;
            if(dst_entry->setting_id == setting_id) {
              found = true;     /* this setting is already defined */
              dst_entry->value.int_ = setting_value;
              dst_entry->type = setting_type;
              break;
            }
            prev = dst_offset;
            dst_offset = dst_entry->next;
          }
          if(!found) {          /* setting not found in existing list, so append new value */
            if(!I->next_free)
              SettingUniqueExpand(G);
            if(I->next_free) {
              dst_offset = I->next_free;
              {
                SettingUniqueEntry *dst_entry = I->entry + dst_offset;
                I->next_free = dst_entry->next;
                dst_entry->next = 0;
                if(prev) {      /* append onto existing list */
                  I->entry[prev].next = dst_offset;
                  dst_entry->type = setting_type;
                  dst_entry->value.int_ = setting_value;
                  dst_entry->setting_id = setting_id;
                } else
                  if(OVreturn_IS_OK
                     (OVOneToOne_Set(I->id2offset, dst_unique_id, dst_offset))) {
                  /* create new list */
                  dst_entry->type = setting_type;
                  dst_entry->value.int_ = setting_value;
                  dst_entry->setting_id = setting_id;
                }
              }
            }
          }
        }
        src_offset = I->entry[src_offset].next; /* src_entry invalid, since I->entry may have changed */
      }
    }
  } else if(dst_result.status == OVstatus_NOT_FOUND) {  /* new setting list for atom */
    OVreturn_word src_result;
    if(OVreturn_IS_OK(src_result = OVOneToOne_GetForward(I->id2offset, src_unique_id))) {
      int src_offset = src_result.word;
      int prev = 0;
      SettingUniqueEntry *src_entry;
      while(ok && src_offset) {
        if(!I->next_free)
          SettingUniqueExpand(G);
        {
          src_entry = I->entry + src_offset;
          {
            int setting_id = src_entry->setting_id;
            int setting_type = src_entry->type;
            int setting_value = src_entry->value.int_;
            if(I->next_free) {
              int dst_offset = I->next_free;
              SettingUniqueEntry *dst_entry = I->entry + dst_offset;
              I->next_free = dst_entry->next;

              if(!prev) {
                if(!OVreturn_IS_OK
                   (OVOneToOne_Set(I->id2offset, dst_unique_id, dst_offset))) {
                  ok = false;
                }
              } else {
                I->entry[prev].next = dst_offset;
              }

              if(ok) {
                dst_entry->type = setting_type;
                dst_entry->value.int_ = setting_value;
                dst_entry->setting_id = setting_id;
                dst_entry->next = 0;
              }
              prev = dst_offset;
            }
          }
        }
        src_offset = I->entry[src_offset].next; /* src_entry invalid, since I->entry may have changed */
      }
    }
  } else {
    ok = false;
    /* unhandled error */
  }

  return ok;
}

static void SettingUniqueInit(PyMOLGlobals * G)
{
  register CSettingUnique *I = G->SettingUnique;

  if((I = (G->SettingUnique = Calloc(CSettingUnique, 1)))) {
    I->id2offset = OVOneToOne_New(G->Context->heap);
    {
      int a;
      I->n_alloc = 10;
      I->entry = VLACalloc(SettingUniqueEntry, I->n_alloc);
      /* note: intentially skip index 0  */
      for(a = 2; a < 10; a++) {
        I->entry[a].next = a - 1;       /* 1-based linked list with 0 as sentinel */
      }
      I->next_free = I->n_alloc - 1;
    }
  }
}

static void SettingUniqueFree(PyMOLGlobals * G)
{
  register CSettingUnique *I = G->SettingUnique;
  VLAFreeP(I->entry);
  OVOneToOne_Del(I->id2offset);
  if(I->old2new)
    OVOneToOne_Del(I->old2new);
  FreeP(I);
}

int SettingUniqueConvertOldSessionID(PyMOLGlobals * G, int old_unique_id)
{
  register CSettingUnique *I = G->SettingUnique;
  int unique_id = old_unique_id;
  if(I->old2new) {
    OVreturn_word ret;
    if(OVreturn_IS_OK(ret = OVOneToOne_GetForward(I->old2new, old_unique_id))) {
      unique_id = ret.word;
    }
  }
  return unique_id;
}

int SettingUniqueFromPyList(PyMOLGlobals * G, PyObject * list, int partial_restore)
{
#ifdef _PYMOL_NOPY
  return 0;
#else
  int ok = true;
  register CSettingUnique *I = G->SettingUnique;
  if(!partial_restore) {
    SettingUniqueResetAll(G);
    if(I->old2new) {
      OVOneToOne_Del(I->old2new);
      I->old2new = NULL;
    }
  } else {
    if(!I->old2new) {
      I->old2new = OVOneToOne_New(G->Context->heap);
    } else {
      OVOneToOne_Reset(I->old2new);
    }
  }
  if(list)
    if(PyList_Check(list)) {
      ov_size n_id = PyList_Size(list);
      ov_size a;
      for(a = 0; a < n_id; a++) {
        PyObject *id_list = PyList_GetItem(list, a);
        int unique_id;
        if(ok)
          ok = PyList_Check(id_list);
        if(ok)
          ok = (PyList_Size(id_list) > 1);
        if(ok)
          ok = PConvPyIntToInt(PyList_GetItem(id_list, 0), &unique_id);
        if(ok && partial_restore) {
          if(AtomInfoIsUniqueIDActive(G, unique_id)) {
            /* if this ID is already active, then we need a substitute */
            int old_unique_id = unique_id;
            unique_id = AtomInfoGetNewUniqueID(G);
            OVOneToOne_Set(I->old2new, old_unique_id, unique_id);
          }
        }
        if(ok) {
          ov_size n_set = 0;

          PyObject *setting_list = PyList_GetItem(id_list, 1);
          if(ok)
            ok = PyList_Check(setting_list);
          if(ok)
            n_set = PyList_Size(setting_list);
          if(ok) {
            ov_size b;
            for(b = 0; b < n_set; b++) {
              PyObject *entry_list = PyList_GetItem(setting_list, b);
              if(ok)
                ok = PyList_Check(entry_list);
              if(ok)
                ok = (PyList_Size(entry_list) > 2);
              if(ok) {
                int setting_id;
                int setting_type;
                union {
                  int int_;
                  float float_;
                } value_store;
                if(ok)
                  ok = PConvPyIntToInt(PyList_GetItem(entry_list, 0), &setting_id);
                if(ok)
                  ok = PConvPyIntToInt(PyList_GetItem(entry_list, 1), &setting_type);
                if(ok)
                  switch (setting_type) {

                  case cSetting_int:
                  case cSetting_color:
                  case cSetting_boolean:
                    ok = PConvPyIntToInt(PyList_GetItem(entry_list, 2),
                                         &value_store.int_);
                    break;
                  case cSetting_float:
                    ok = PConvPyFloatToFloat(PyList_GetItem(entry_list, 2),
                                             &value_store.float_);
                    break;
                  }
                if(ok) {
                  SettingUniqueSetTypedValue(G, unique_id, setting_id,
                                             setting_type, &value_store.int_);
                }
              }
            }
          }
        }
      }
    }
  return ok;
#endif
}

PyObject *SettingUniqueAsPyList(PyMOLGlobals * G)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else
  PyObject *result = NULL;
  register CSettingUnique *I = G->SettingUnique;
  {
    ov_word hidden = 0;
    OVreturn_word ret;
    int n_entry = 0;
    while(1) {
      ret = OVOneToOne_IterateForward(I->id2offset, &hidden);
      if(ret.status != OVstatus_YES)
        break;
      n_entry++;
    }
    result = PyList_New(n_entry);
    if(result) {
      hidden = 0;
      n_entry = 0;
      while(1) {
        PyObject *setting_list = NULL;
        int save_offset, unique_id;
        ret = OVOneToOne_IterateForward(I->id2offset, &hidden);

        if(ret.status != OVstatus_YES)
          break;
        unique_id = ret.word;
        if(OVreturn_IS_OK(ret = OVOneToOne_GetForward(I->id2offset, unique_id))) {
          int offset = ret.word;
          int n_set = 0;

          /* count number of settings for this unique_id */

          SettingUniqueEntry *entry;
          save_offset = offset;
          while(offset) {
            entry = I->entry + offset;
            n_set++;
            offset = entry->next;
          }

          /* create and insert list for each setting */

          setting_list = PyList_New(n_set);
          n_set = 0;
          offset = save_offset;
          while(offset) {
            PyObject *setting_entry = PyList_New(3);
            entry = I->entry + offset;
            PyList_SetItem(setting_entry, 0, PyInt_FromLong(entry->setting_id));
            PyList_SetItem(setting_entry, 1, PyInt_FromLong(entry->type));
            switch (entry->type) {
            case cSetting_int:
            case cSetting_color:
            case cSetting_boolean:
              PyList_SetItem(setting_entry, 2, PyInt_FromLong(entry->value.int_));
              break;
            case cSetting_float:
              PyList_SetItem(setting_entry, 2,
                             PyFloat_FromDouble(*(float *) &entry->value.float_));
              break;
            }
            PyList_SetItem(setting_list, n_set, setting_entry);
            n_set++;
            offset = entry->next;
          }
        }

        /* add this unique_id set into the overall list */

        {
          PyObject *unique_list = PyList_New(2);
          PyList_SetItem(unique_list, 0, PyInt_FromLong(unique_id));
          PyList_SetItem(unique_list, 1, setting_list);
          PyList_SetItem(result, n_entry, unique_list);
        }
        n_entry++;
      }
    }
  }
  return (PConvAutoNone(result));
#endif
}

int SettingSetSmart_i(PyMOLGlobals * G, CSetting * set1, CSetting * set2, int index,
                      int value)
{
  int dummy;
  if(set1 && SettingGetIfDefined_i(G, set1, index, &dummy)) {
    return SettingSet_i(set1, index, value);
  }
  if(set2 && SettingGetIfDefined_i(G, set2, index, &dummy)) {
    return SettingSet_i(set2, index, value);
  }
  return SettingSetGlobal_i(G, index, value);
}

int SettingSetGlobalsFromPyList(PyMOLGlobals * G, PyObject * list)
{
#ifdef _PYMOL_NOPY
  return 0;
#else
  int ok = true;

  int session_migration = SettingGetGlobal_b(G, cSetting_session_migration);
  int session_version_check = SettingGetGlobal_b(G, cSetting_session_version_check);
  int full_screen = SettingGetGlobal_b(G, cSetting_full_screen);
  int internal_gui = SettingGetGlobal_b(G, cSetting_internal_gui);
  int internal_feedback = SettingGetGlobal_b(G, cSetting_internal_feedback);

  int stereo = SettingGetGlobal_b(G, cSetting_stereo);
  int text = SettingGetGlobal_b(G, cSetting_text);
  int texture_fonts = SettingGetGlobal_b(G, cSetting_texture_fonts);
  int use_display_lists = SettingGetGlobal_b(G, cSetting_use_display_lists);
  int max_threads = SettingGetGlobal_i(G, cSetting_max_threads);
  int nvidia_bugs = SettingGetGlobal_b(G, cSetting_nvidia_bugs);
  int ati_bugs = SettingGetGlobal_b(G, cSetting_ati_bugs);
  int stereo_mode = SettingGetGlobal_i(G, cSetting_stereo_mode);
  int stereo_double_pump_mono = SettingGetGlobal_b(G, cSetting_stereo_double_pump_mono);
  int show_progress = SettingGetGlobal_b(G, cSetting_show_progress);
  int defer_updates = SettingGetGlobal_b(G, cSetting_defer_updates);
  int suspend_updates = SettingGetGlobal_b(G, cSetting_suspend_updates);
  int cache_max = SettingGetGlobal_i(G, cSetting_cache_max);
  int logging = SettingGetGlobal_i(G, cSetting_logging);
  float no_idle = SettingGetGlobal_f(G, cSetting_no_idle);
  float slow_idle = SettingGetGlobal_f(G, cSetting_fast_idle);
  float fast_idle = SettingGetGlobal_f(G, cSetting_slow_idle);
  int mouse_grid = SettingGetGlobal_b(G, cSetting_mouse_grid);
  int mouse_z_scale = SettingGetGlobal_i(G,cSetting_mouse_scale);
  register CSetting *I = G->Setting;
  if(list)
    if(PyList_Check(list))
      ok = SettingFromPyList(I, list);

  SettingSet_i(I, cSetting_security, G->Security);      /* always override Security setting with global variable */
  SettingSet_b(I, cSetting_session_migration, session_migration);       /* preserve current migration info */
  SettingSet_b(I, cSetting_session_version_check, session_version_check);

  /* restore the following settings  */

  SettingSetGlobal_f(G, cSetting_no_idle, no_idle);
  SettingSetGlobal_f(G, cSetting_fast_idle, fast_idle);
  SettingSetGlobal_f(G, cSetting_slow_idle, slow_idle);

  SettingSet_b(I, cSetting_stereo, stereo);
  SettingSet_b(I, cSetting_text, text);
  SettingSet_b(I, cSetting_texture_fonts, texture_fonts);
  SettingSet_b(I, cSetting_use_display_lists, use_display_lists);
  SettingSet_i(I, cSetting_max_threads, max_threads);
  SettingSet_i(I, cSetting_nvidia_bugs, nvidia_bugs);
  SettingSet_i(I, cSetting_ati_bugs, ati_bugs);
  SettingSet_i(I, cSetting_cache_max, cache_max);
  SettingSet_i(I, cSetting_logging, logging);

  SettingSet_i(I, cSetting_stereo_mode, stereo_mode);
  SettingSet_b(I, cSetting_stereo_double_pump_mono, stereo_double_pump_mono);
  SettingSet_b(I, cSetting_full_screen, full_screen);
  SettingSet_b(I, cSetting_show_progress, show_progress);
  SettingSet_b(I, cSetting_defer_updates, defer_updates);
  SettingSet_b(I, cSetting_suspend_updates, suspend_updates);
  SettingSet_b(I, cSetting_session_changed, 0);

  SettingSet_b(I, cSetting_mouse_grid, mouse_grid);
  SettingSet_i(I, cSetting_mouse_z_scale, mouse_z_scale);
  if(G->Option->presentation) {
    SettingSet_b(I, cSetting_full_screen, full_screen);
    SettingSet_b(I, cSetting_presentation, 1);
    SettingSet_b(I, cSetting_internal_gui, internal_gui);
    SettingSet_b(I, cSetting_internal_feedback, internal_feedback);
  }
  if(G->Option->no_quit) {
    SettingSet_b(I, cSetting_presentation_auto_quit, 0);
  }

#ifdef _PYMOL_ACTIVEX
  SettingSet_i(I, cSetting_max_threads, 1);
  SettingSet_b(I, cSetting_async_builds, 0);
#endif
  return (ok);
#endif
}

PyObject *SettingGetGlobalsAsPyList(PyMOLGlobals * G)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else
  PyObject *result = NULL;
  register CSetting *I = G->Setting;
  result = SettingAsPyList(I);
  return (PConvAutoNone(result));
#endif
}

#ifndef _PYMOL_NOPY
static PyObject *get_list(CSetting * I, int index)
{
  PyObject *result = NULL;
  int setting_type = I->info[index].type;
  switch (setting_type) {

  case cSetting_boolean:
  case cSetting_int:
  case cSetting_color:
    result = PyList_New(3);
    PyList_SetItem(result, 0, PyInt_FromLong(index));
    PyList_SetItem(result, 1, PyInt_FromLong(setting_type));
    PyList_SetItem(result, 2,
                   PyInt_FromLong(*((int *) (I->data + I->info[index].offset))));
    break;
  case cSetting_float:
    result = PyList_New(3);
    PyList_SetItem(result, 0, PyInt_FromLong(index));
    PyList_SetItem(result, 1, PyInt_FromLong(setting_type));
    PyList_SetItem(result, 2,
                   PyFloat_FromDouble(*((float *) (I->data + I->info[index].offset))));
    break;
  case cSetting_float3:
    result = PyList_New(3);
    PyList_SetItem(result, 0, PyInt_FromLong(index));
    PyList_SetItem(result, 1, PyInt_FromLong(setting_type));
    PyList_SetItem(result, 2,
                   PConvFloatArrayToPyList(((float *) (I->data + I->info[index].offset)),
                                           3));
    break;
  case cSetting_string:
    result = PyList_New(3);
    PyList_SetItem(result, 0, PyInt_FromLong(index));
    PyList_SetItem(result, 1, PyInt_FromLong(setting_type));
    PyList_SetItem(result, 2,
                   PyString_FromString(((char *) (I->data + I->info[index].offset))));
    break;
  default:
    result = Py_None;
    break;
  }
  return (PConvAutoNone(result));
}
#endif

PyObject *SettingAsPyList(CSetting * I)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else

  PyObject *result = NULL;
  int cnt = 0;
  int a;

  if(I) {
    for(a = 0; a < cSetting_INIT; a++) {
      if(I->info[a].defined)
        cnt++;
    }
    result = PyList_New(cnt);
    cnt = 0;
    for(a = 0; a < cSetting_INIT; a++) {
      if(I->info[a].defined) {
        PyList_SetItem(result, cnt, get_list(I, a));
        cnt++;
      }
    }
  }
  return (PConvAutoNone(result));
#endif
}


/*========================================================================*/
#ifndef _PYMOL_NOPY
static int set_list(CSetting * I, PyObject * list)
{
  int ok = true;
  int index;
  int setting_type;
  char *str;
  if(list != Py_None) {
    if(ok)
      ok = (list != NULL);
    if(ok)
      ok = PyList_Check(list);
    if(ok)
      ok = PConvPyIntToInt(PyList_GetItem(list, 0), &index);
    if(ok)
      ok = PConvPyIntToInt(PyList_GetItem(list, 1), &setting_type);
    if(ok && (index < cSetting_INIT)) { /* ignore unknown settings */
      switch (index) {
        /* don't restore the folllowing settings,
           which are inherently system-dependent */
      case cSetting_stereo_double_pump_mono:
      case cSetting_max_threads:
      case cSetting_session_migration:
        break;
      default:
        if(ok)
          switch (setting_type) {
          case cSetting_boolean:
          case cSetting_int:
            ok = PConvPyIntToInt(PyList_GetItem(list, 2),
                                 (int *) SettingPtr(I, index, sizeof(int)));
            break;
          case cSetting_color:
            {
              int color = 0;
              ok = PConvPyIntToInt(PyList_GetItem(list, 2), &color);
              if(ok)
                color = ColorConvertOldSessionIndex(I->G, color);
              *((int *) SettingPtr(I, index, sizeof(int))) = color;
            }
            break;
          case cSetting_float:
            ok = PConvPyFloatToFloat(PyList_GetItem(list, 2),
                                     (float *) SettingPtr(I, index, sizeof(float)));
            break;
          case cSetting_float3:
            ok = PConvPyListToFloatArrayInPlaceAutoZero(PyList_GetItem(list, 2),
                                                        (float *) SettingPtr(I, index,
                                                                             3 *
                                                                             sizeof
                                                                             (float)), 3);
            break;
          case cSetting_string:
            ok = PConvPyStrToStrPtr(PyList_GetItem(list, 2), &str);
            if(ok) {
              strcpy(((char *) SettingPtr(I, index, strlen(str) + 1)), str);
            }
            break;
          }
      }
      if(ok)
        I->info[index].type = setting_type;
    }
  }
  return (ok);
}
#endif

/*========================================================================*/
CSetting *SettingNewFromPyList(PyMOLGlobals * G, PyObject * list)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else

  int ok = true;
  ov_size size;
  ov_size a;
  CSetting *I = NULL;
  if(ok)
    ok = (list != NULL);
  if(ok)
    ok = PyList_Check(list);
  if(ok) {
    I = SettingNew(G);
    size = PyList_Size(list);
    for(a = 0; a < size; a++) {
      if(ok)
        ok = set_list(I, PyList_GetItem(list, a));
    }
  }
  return (I);
#endif
}


/*========================================================================*/
int SettingFromPyList(CSetting * I, PyObject * list)
{
#ifdef _PYMOL_NOPY
  return 0;
#else

  int ok = true;
  ov_size size;
  ov_size a;
  if(ok)
    ok = (I != NULL);
  if(ok)
    ok = PyList_Check(list);
  if(ok) {
    size = PyList_Size(list);
    for(a = 0; a < size; a++) {
      if(!set_list(I, PyList_GetItem(list, a)))
        ok = false;
    }
  }
  return (ok);
#endif
}


/*========================================================================*/
#ifndef _PYMOL_NOPY
PyObject *SettingGetUpdateList(PyMOLGlobals * G, CSetting * I)
{                               /* assumes blocked interpreter */

  int a;
  int n;
  PyObject *result;

  if(!I)
    I = G->Setting;             /* fall back on global settings */

  n = VLAGetSize(I->info);
  result = PyList_New(0);
  for(a = 0; a < n; a++) {
    if(I->info[a].changed) {
      I->info[a].changed = false;
      PyList_Append(result, PyInt_FromLong(a));
    }
  }
  return (result);

}
#endif


/*========================================================================*/
void SettingCheckHandle(PyMOLGlobals * G, CSetting ** handle)
{
  if(!*handle)
    *handle = SettingNew(G);
}


/*========================================================================*/
int SettingGetTextValue(PyMOLGlobals * G, CSetting * set1, CSetting * set2, int index,
                        char *buffer)

/* not range checked */
{
  int type;
  int ok = true;
  float *ptr;
  type = SettingGetType(G, index);
  switch (type) {
  case cSetting_boolean:
    if(SettingGet_b(G, set1, set2, index))
      sprintf(buffer, "on");
    else
      sprintf(buffer, "off");
    break;
  case cSetting_int:
    sprintf(buffer, "%d", SettingGet_i(G, set1, set2, index));
    break;
  case cSetting_float:
    sprintf(buffer, "%1.5f", SettingGet_f(G, set1, set2, index));
    break;
  case cSetting_float3:
    ptr = SettingGet_3fv(G, set1, set2, index);
    sprintf(buffer, "[ %1.5f, %1.5f, %1.5f ]", ptr[0], ptr[1], ptr[2]);
    break;
  case cSetting_color:
    {
      int color = SettingGet_color(G, set1, set2, index);
      if(color < 0) {
        switch (color) {
        case cColorAtomic:
          strcpy(buffer, "atomic");
          break;
        case cColorObject:
          strcpy(buffer, "object");
          break;
        case cColorFront:
          strcpy(buffer, "front");
          break;
        case cColorBack:
          strcpy(buffer, "back");
          break;
        default:
          if(color > cColorExtCutoff) {
            strcpy(buffer, "default");
          } else {
            char *st = ColorGetName(G, color);
            if(st)
              strcpy(buffer, st);
            else
              strcpy(buffer, "invalid");
          }
          break;
        }
      } else {
        /* assuming valid color */
        strcpy(buffer, ColorGetName(G, color));
      }
    }
    break;
  case cSetting_string:
    strcpy(buffer, SettingGet_s(G, set1, set2, index));
    break;
  default:
    ok = false;
    break;
  }
  return (ok);
}


/*========================================================================*/
#ifndef _PYMOL_NOPY
int SettingSetFromTuple(PyMOLGlobals * G, CSetting * I, int index, PyObject * tuple)

/* must have interpret locked to make this call */
{
  PyObject *value;
  int type;
  int ok = true;
  if(!I)
    I = G->Setting;             /* fall back on global settings */

  /* this data structure has been pre-checked at the python level... */

  type = PyInt_AsLong(PyTuple_GetItem(tuple, 0));
  value = PyTuple_GetItem(tuple, 1);
  switch (type) {
  case cSetting_boolean:
    SettingSet_b(I, index, PyInt_AsLong(PyTuple_GetItem(value, 0)));
    break;
  case cSetting_int:
    SettingSet_i(I, index, PyInt_AsLong(PyTuple_GetItem(value, 0)));
    break;
  case cSetting_float:
    SettingSet_f(I, index, (float) PyFloat_AsDouble(PyTuple_GetItem(value, 0)));
    break;
  case cSetting_float3:
    SettingSet_3f(I, index,
                  (float) PyFloat_AsDouble(PyTuple_GetItem(value, 0)),
                  (float) PyFloat_AsDouble(PyTuple_GetItem(value, 1)),
                  (float) PyFloat_AsDouble(PyTuple_GetItem(value, 2)));
    break;
  case cSetting_color:
    SettingSet_color(I, index, PyString_AsString(PyTuple_GetItem(value, 0)));
    break;
  case cSetting_string:
    SettingSet_s(I, index, PyString_AsString(PyTuple_GetItem(value, 0)));
    break;
  default:
    ok = false;
    break;
  }
  return (ok);
}
#endif

/*========================================================================*/
int SettingStringToTypedValue(PyMOLGlobals * G, int index, char *st, int *type,
                              int *value)
{
  int ok = true;

  /* this data structure has been pre-checked at the python level... */

  *type = SettingGetType(G, index);

  switch (*type) {
  case cSetting_boolean:
    if((!*st) || (*st == '0') || (*st == 'F') || WordMatchExact(G, st, "on", true)
       || WordMatchExact(G, st, "false", true))
      *value = 0;
    else
      *value = 1;
    break;
  case cSetting_int:
    if(sscanf(st, "%d", value) != 1)
      ok = false;
    break;
  case cSetting_float:
    if(sscanf(st, "%f", (float *) value) != 1)
      ok = false;
    break;
  case cSetting_color:
    {
      int color_index = ColorGetIndex(G, st);
      if((color_index < 0) && (color_index > cColorExtCutoff))
        color_index = 0;
      *(value) = color_index;
    }
    break;
  default:
    ok = false;
    break;
  }
  return (ok);
}

int SettingSetFromString(PyMOLGlobals * G, CSetting * I, int index, char *st)
{
  int type;
  int ok = true;
  if(!I)
    I = G->Setting;             /* fall back on global settings */

  /* this data structure has been pre-checked at the python level... */

  type = SettingGetType(G, index);

  switch (type) {
  case cSetting_boolean:
    if((!*st) || (*st == '0') || (*st == 'F') || WordMatchExact(G, st, "on", true)
       || WordMatchExact(G, st, "false", true))
      SettingSet_b(I, index, 0);
    else
      SettingSet_b(I, index, 1);
    break;
  case cSetting_int:
    {
      int tmp;
      if(sscanf(st, "%d", &tmp) == 1)
        SettingSet_i(I, index, tmp);
      else
        ok = false;
    }
    break;
  case cSetting_float:
    {
      float tmp;
      if(sscanf(st, "%f", &tmp) == 1)
        SettingSet_f(I, index, tmp);
      else
        ok = false;
    }
    break;
  case cSetting_float3:
    {
      float tmp1, tmp2, tmp3;
      if(sscanf(st, "%f%f%f", &tmp1, &tmp2, &tmp3) == 3)
        SettingSet_3f(I, index, tmp1, tmp2, tmp3);
      else
        ok = false;
    }
    break;
  case cSetting_color:
    SettingSet_color(I, index, st);
    break;
  case cSetting_string:
    SettingSet_s(I, index, st);
    break;
  default:
    ok = false;
    break;
  }
  return (ok);
}


/*========================================================================*/
#ifndef _PYMOL_NOPY
PyObject *SettingGetTuple(PyMOLGlobals * G, CSetting * set1, CSetting * set2, int index)
{                               /* assumes blocked python interpeter */
  PyObject *result = NULL;
  float *ptr;
  int type = SettingGetType(G, index);

  switch (type) {
  case cSetting_boolean:
    result = Py_BuildValue("(i(i))", type, SettingGet_b(G, set1, set2, index));
    break;
  case cSetting_int:
    result = Py_BuildValue("(i(i))", type, SettingGet_i(G, set1, set2, index));
    break;
  case cSetting_float:
    result = Py_BuildValue("(i(f))", type, SettingGet_f(G, set1, set2, index));
    break;
  case cSetting_float3:
    ptr = SettingGet_3fv(G, set1, set2, index);
    result = Py_BuildValue("(i(fff))", type, ptr[0], ptr[1], ptr[2]);
    break;
  case cSetting_color:
    result = Py_BuildValue("(i(i))", type, SettingGet_color(G, set1, set2, index));
    break;
  case cSetting_string:
    result = Py_BuildValue("(i(s))", type, SettingGet_s(G, set1, set2, index));
    break;
  default:
    result = PConvAutoNone(Py_None);
    break;
  }
  return result;
}
#endif

/*========================================================================*/
#ifndef _PYMOL_NOPY
PyObject *SettingGetDefinedTuple(PyMOLGlobals * G, CSetting * set1, int index)
{                               /* Assumes blocked Python interpreter */
  PyObject *result = NULL;
  int defined = true;
  int type = SettingGetType(G, index);
  int int1;
  float float1, *vect1 = NULL;
  char *str1;
  switch (type) {
  case cSetting_boolean:
    defined = SettingGetIfDefined_b(G, set1, index, &int1);
    if(defined)
      result = Py_BuildValue("(i(i))", type, int1);
    break;
  case cSetting_int:
    defined = SettingGetIfDefined_i(G, set1, index, &int1);
    if(defined)
      result = Py_BuildValue("(i(i))", type, int1);
    break;
  case cSetting_float:
    defined = SettingGetIfDefined_f(G, set1, index, &float1);
    if(defined)
      result = Py_BuildValue("(i(f))", type, float1);
    break;
  case cSetting_float3:
    defined = SettingGetIfDefined_3fv(G, set1, index, &vect1);
    result = Py_BuildValue("(i(fff))", type, vect1[0], vect1[1], vect1[2]);
    break;
  case cSetting_color:
    defined = SettingGetIfDefined_color(G, set1, index, &int1);
    if(defined)
      result = Py_BuildValue("(i(i))", type, int1);
    break;
  case cSetting_string:
    defined = SettingGetIfDefined_s(G, set1, index, &str1);
    if(defined)
      result = Py_BuildValue("(i(s))", type, str1);
    break;
  default:
    break;
  }
  if(!defined) {
    result = Py_BuildValue("(i)", 0);
  }
  if(!result) {
    result = PConvAutoNone(Py_None);
  }
  return result;
}
#endif


/*========================================================================*/
CSetting *SettingNew(PyMOLGlobals * G)
{
  OOAlloc(G, CSetting);
  SettingInit(G, I);
  return (I);
}


/*========================================================================*/
void SettingPurge(CSetting * I)
{
  if(I) {
    VLAFreeP(I->data);
    VLAFreeP(I->info);
    I->size = 0;
  }
}


/*========================================================================*/
void SettingFreeP(CSetting * I)
{
  if(I)
    SettingPurge(I);
  OOFreeP(I);
}


/*========================================================================*/
void SettingInit(PyMOLGlobals * G, CSetting * I)
{
  I->G = G;
  I->size = sizeof(int);        /* insures offset is never zero, except when undef */
  I->data = VLAlloc(char, 10);
  I->info = VLAMalloc(cSetting_INIT, sizeof(SettingRec), 5, 1); /* auto-zero */
}


/*========================================================================*/
void SettingClear(CSetting * I, int index)
{
  if(I)
    I->info[index].defined = false;
}


/*========================================================================*/
static void *SettingPtr(CSetting * I, int index, ov_size size)
{
  /* note that this routine essentially leaks RAM in terms of not
     recovering space used for previous settings of smaller size */

  VLACheck(I->info, SettingRec, index);
  {
    SettingRec *sr = I->info + index;
    if(size < sizeof(int))
      size = sizeof(int);       /* make sure we're word aligned */
    while(size & (sizeof(int) - 1))
      size++;

    if((!sr->offset) || (sr->max_size < size)) {
      sr->offset = I->size;
      I->size += size;
      sr->max_size = size;
      VLACheck(I->data, char, I->size);
    }
    sr->defined = true;
    sr->changed = true;
    return (I->data + sr->offset);
  }
}


/*========================================================================*/
int SettingUnset(CSetting * I, int index)
{
  if(I) {
    SettingRec *sr = I->info + index;
    sr->defined = false;
    sr->changed = true;
  }
  return true;
}


/*========================================================================*/
int SettingGetType(PyMOLGlobals * G, int index)
{
  register CSetting *I = G->Setting;
  return (I->info[index].type);
}


/*========================================================================*/
static int get_i(CSetting * I, int index)
{
  PyMOLGlobals *G = I->G;
  int result;
  switch (I->info[index].type) {
  case cSetting_boolean:
  case cSetting_int:
  case cSetting_color:
    result = (*((int *) (I->data + I->info[index].offset)));
    break;
  case cSetting_float:
    result = (int) (*((float *) (I->data + I->info[index].offset)));
    break;
  default:
    PRINTFB(G, FB_Setting, FB_Errors)
      "Setting-Error: type read mismatch (int) %d\n", index ENDFB(G);
    result = 0;
    break;
  }
  return (result);
}


/*========================================================================*/
static int get_b(CSetting * I, int index)
{
  int result;
  PyMOLGlobals *G = I->G;
  switch (I->info[index].type) {
  case cSetting_boolean:
  case cSetting_int:
  case cSetting_color:
    result = (*((int *) (I->data + I->info[index].offset)));
    break;
  case cSetting_float:
    result = (int) (*((float *) (I->data + I->info[index].offset)));
    break;
  default:
    PRINTFB(G, FB_Setting, FB_Errors)
      "Setting-Error: type read mismatch (boolean) %d\n", index ENDFB(G);
    result = 0;
  }
  return (result);
}


/*========================================================================*/
static int get_color(CSetting * I, int index)
{
  int result;
  PyMOLGlobals *G = I->G;
  switch (I->info[index].type) {
  case cSetting_boolean:
  case cSetting_int:
  case cSetting_color:
    result = (*((int *) (I->data + I->info[index].offset)));
    break;
  case cSetting_float:
    result = (int) (*((float *) (I->data + I->info[index].offset)));
    break;
  default:
    PRINTFB(G, FB_Setting, FB_Errors)
      "Setting-Error: type read mismatch (color) %d\n", index ENDFB(G);
    result = 0;
  }
  return (result);
}


/*========================================================================*/
static float get_f(CSetting * I, int index)
{
  float result;
  PyMOLGlobals *G = I->G;
  switch (I->info[index].type) {
  case cSetting_boolean:
  case cSetting_int:
  case cSetting_color:
    result = (float) (*((int *) (I->data + I->info[index].offset)));
    break;
  case cSetting_float:
    result = (*((float *) (I->data + I->info[index].offset)));
    break;
  default:
    PRINTFB(G, FB_Setting, FB_Errors)
      "Setting-Error: type read mismatch (float) %d\n", index ENDFB(G);
    result = 0.0F;
  }
  return (result);
}


/*========================================================================*/
static char *get_s(CSetting * I, int index)
{
  char *result;
  PyMOLGlobals *G = I->G;
  switch (I->info[index].type) {
  case cSetting_string:
    result = ((char *) (I->data + I->info[index].offset));
    break;
  default:
    PRINTFB(G, FB_Setting, FB_Errors)
      "Setting-Error: type read mismatch (string) %d\n", index ENDFB(G);
    result = NULL;
  }
  return (result);
}


/*========================================================================*/
int SettingSet_b(CSetting * I, int index, int value)
{
  int ok = true;
  if(I) {
    VLACheck(I->info, SettingRec, index);
    {
      int setting_type = I->info[index].type;
      PyMOLGlobals *G = I->G;
      switch (setting_type) {
      case cSetting_blank:
      case cSetting_boolean:
      case cSetting_int:
      case cSetting_color:
	*((int *) SettingPtr(I, index, sizeof(int))) = value;
	if(setting_type == cSetting_blank)
	  I->info[index].type = cSetting_boolean;
	break;
      case cSetting_float:
	*((float *) SettingPtr(I, index, sizeof(float))) = (float) value;
	break;
      default:
	PRINTFB(G, FB_Setting, FB_Errors)
	  "Setting-Error: type set mismatch (boolean) %d\n", index ENDFB(G);
	ok = false;
      }
    }
  } else {
    ok = false;
  }
  return (ok);
}


/*========================================================================*/
int SettingSet_i(CSetting * I, int index, int value)
{
  int ok = true;
  if(I) {
    PyMOLGlobals *G = I->G;
    VLACheck(I->info, SettingRec, index);
    {
      int setting_type = I->info[index].type;
      switch (setting_type) {
      case cSetting_blank:
      case cSetting_boolean:
      case cSetting_int:
      case cSetting_color:
	*((int *) SettingPtr(I, index, sizeof(int))) = value;
	if(setting_type == cSetting_blank)
	  I->info[index].type = cSetting_int;
	break;
      case cSetting_float:
	*((float *) SettingPtr(I, index, sizeof(float))) = (float) value;
	break;
      default:
	PRINTFB(G, FB_Setting, FB_Errors)
	  "Setting-Error: type set mismatch (integer)\n" ENDFB(G);
	ok = false;
      }
    }
  } else {
    ok = false;
  }
  return (ok);
}


/*========================================================================*/
int SettingSet_color(CSetting * I, int index, char *value)
{
  int ok = true;
  int color_index;
  if(I) {
    PyMOLGlobals *G = I->G;
    color_index = ColorGetIndex(G, value);
    if((color_index == -1) && (strcmp(value, "-1") &&
                               strcmp(value, "-2") &&
                               strcmp(value, "-3") &&
                               strcmp(value, "-4") &&
                               strcmp(value, "-5") && strcmp(value, "default"))) {
      PRINTFB(G, FB_Setting, FB_Errors)
        "Setting-Error: unknown color '%s'\n", value ENDFB(G);
      ok = false;

    } else {
      VLACheck(I->info, SettingRec, index);
      {
	int setting_type = I->info[index].type;
	switch (setting_type) {
	case cSetting_blank:
	case cSetting_boolean:
	case cSetting_int:
	case cSetting_color:
	  *((int *) SettingPtr(I, index, sizeof(int))) = color_index;
	  if(setting_type == cSetting_blank)
	    I->info[index].type = cSetting_color;
	  break;
	case cSetting_float:
	  *((float *) SettingPtr(I, index, sizeof(float))) = (float) color_index;
	  break;
	default:
	  PRINTFB(G, FB_Setting, FB_Errors)
	    "Setting-Error: type set mismatch (color)\n" ENDFB(G);
	  ok = false;
	}
      }
    }
  }
  return (ok);
}


/*========================================================================*/
int SettingSet_f(CSetting * I, int index, float value)
{
  int ok = true;
  if(I) {
    PyMOLGlobals *G = I->G;
    VLACheck(I->info, SettingRec, index);
    {
      int setting_type = I->info[index].type;
      switch (setting_type) {
      case cSetting_boolean:
      case cSetting_int:
      case cSetting_color:
	*((int *) SettingPtr(I, index, sizeof(int))) = (int) value;
	break;
      case cSetting_blank:
      case cSetting_float:
	*((float *) SettingPtr(I, index, sizeof(float))) = value;
	if(setting_type == cSetting_blank)
	  I->info[index].type = cSetting_float;
	break;
      default:
	PRINTFB(G, FB_Setting, FB_Errors)
	  "Setting-Error: type set mismatch (float)\n" ENDFB(G);
	ok = false;
      }
    }
  } else {
    ok = false;
  }
  return (ok);
}


/*========================================================================*/
int SettingSet_s(CSetting * I, int index, char *value)
{
  int ok = true;
  if(I) {
    PyMOLGlobals *G = I->G;
    VLACheck(I->info, SettingRec, index);
    {
      int setting_type = I->info[index].type;
      switch (setting_type) {
      case cSetting_blank:
      case cSetting_string:
	strcpy(((char *) SettingPtr(I, index, strlen(value) + 1)), value);
	if(setting_type == cSetting_blank)
	  I->info[index].type = cSetting_string;
	break;
      default:
	PRINTFB(G, FB_Setting, FB_Errors)
	  "Setting-Error: type set mismatch (string)\n" ENDFB(G);
	ok = false;
      }
      if(setting_type == cSetting_blank)
	I->info[index].type = cSetting_string;
    }
  } else {
    ok = false;
  }
  return (ok);
}


/*========================================================================*/
int SettingSet_3f(CSetting * I, int index, float value1, float value2, float value3)
{
  int ok = false;
  float *ptr;
  if(I) {
    PyMOLGlobals *G = I->G;
    VLACheck(I->info, SettingRec, index);
    {
      int setting_type = I->info[index].type;
      switch (setting_type) {
      case cSetting_blank:
      case cSetting_float3:
	ptr = (float *) SettingPtr(I, index, sizeof(float) * 3);
	ptr[0] = value1;
	ptr[1] = value2;
	ptr[2] = value3;
	if(setting_type == cSetting_blank)
	  I->info[index].type = cSetting_float3;
	break;
      default:
	PRINTFB(G, FB_Setting, FB_Errors)
	  "Setting-Error: type set mismatch (float3)\n" ENDFB(G);
	ok = false;
      }
    }
  } else {
    ok = false;
  }
  return (ok);
}


/*========================================================================*/
int SettingSet_3fv(CSetting * I, int index, float *vector)
{
  float *ptr;
  ptr = (float *) SettingPtr(I, index, sizeof(float) * 3);
  copy3f(vector, ptr);
  I->info[index].type = cSetting_float3;
  return (true);
}


/*========================================================================*/
int SettingGetGlobal_b(PyMOLGlobals * G, int index)
{
  register CSetting *I = G->Setting;
  return (get_b(I, index));
}


/*========================================================================*/
int SettingGetGlobal_i(PyMOLGlobals * G, int index)
{
  register CSetting *I = G->Setting;
  return (get_i(I, index));
}


/*========================================================================*/
float SettingGetGlobal_f(PyMOLGlobals * G, int index)
{
  register CSetting *I = G->Setting;
  return (get_f(I, index));
}


/*========================================================================*/
char *SettingGetGlobal_s(PyMOLGlobals * G, int index)
{
  register CSetting *I = G->Setting;
  return (get_s(I, index));
}

int SettingGetGlobal_color(PyMOLGlobals * G, int index)
{
  register CSetting *I = G->Setting;
  return (get_color(I, index));
}


/*========================================================================*/
void SettingGetGlobal_3f(PyMOLGlobals * G, int index, float *value)
{
  register CSetting *I = G->Setting;
  float *ptr;
  ptr = (float *) (I->data + I->info[index].offset);
  copy3f(ptr, value);
}


/*========================================================================*/
float *SettingGetGlobal_3fv(PyMOLGlobals * G, int index)
{
  register CSetting *I = G->Setting;
  return (float *) (I->data + I->info[index].offset);
}


/*========================================================================*/
int SettingGet_b(PyMOLGlobals * G, CSetting * set1, CSetting * set2, int index)
{
  if(set1) {
    if(set1->info[index].defined) {
      return (get_b(set1, index));
    }
  }
  if(set2) {
    if(set2->info[index].defined) {
      return (get_b(set2, index));
    }
  }
  return (SettingGetGlobal_i(G, index));
}


/*========================================================================*/
int SettingGetIfDefined_b(PyMOLGlobals * G, CSetting * set1, int index, int *value)
{
  int result = false;
  if(set1) {
    if(set1->info[index].defined) {
      *value = get_b(set1, index);
      result = true;
    }
  }
  return (result);
}


/*========================================================================*/
int SettingGetIfDefined_i(PyMOLGlobals * G, CSetting * set1, int index, int *value)
{
  int result = false;
  if(set1) {
    if(set1->info[index].defined) {
      *value = get_i(set1, index);
      result = true;
    }
  }
  return (result);
}


/*========================================================================*/
int SettingGetIfDefined_color(PyMOLGlobals * G, CSetting * set1, int index, int *value)
{
  int result = false;
  if(set1) {
    if(set1->info[index].defined) {
      *value = get_color(set1, index);
      result = true;
    }
  }
  return (result);
}


/*========================================================================*/
int SettingGetIfDefined_f(PyMOLGlobals * G, CSetting * set1, int index, float *value)
{
  int result = false;
  if(set1) {
    if(set1->info[index].defined) {
      *value = get_f(set1, index);
      result = true;
    }
  }
  return (result);
}


/*========================================================================*/
int SettingGetIfDefined_3fv(PyMOLGlobals * G, CSetting * set1, int index, float **value)
{
  int result = false;
  if(set1) {
    if(set1->info[index].defined) {
      (*value) = (float *) (set1->data + set1->info[index].offset);
      result = true;
    }
  }
  return (result);
}


/*========================================================================*/
int SettingGetIfDefined_s(PyMOLGlobals * G, CSetting * set1, int index, char **value)
{
  int result = false;
  if(set1) {
    if(set1->info[index].defined) {
      *value = get_s(set1, index);
      result = true;
    }
  }
  return (result);
}


/*========================================================================*/
int SettingGet_i(PyMOLGlobals * G, CSetting * set1, CSetting * set2, int index)
{
  if(set1) {
    if(set1->info[index].defined) {
      return (get_i(set1, index));
    }
  }
  if(set2) {
    if(set2->info[index].defined) {
      return (get_i(set2, index));
    }
  }
  return (SettingGetGlobal_i(G, index));
}


/*========================================================================*/
int SettingGet_color(PyMOLGlobals * G, CSetting * set1, CSetting * set2, int index)
{
  if(set1) {
    if(set1->info[index].defined) {
      return (get_color(set1, index));
    }
  }
  if(set2) {
    if(set2->info[index].defined) {
      return (get_color(set2, index));
    }
  }
  return (SettingGetGlobal_i(G, index));
}


/*========================================================================*/
float SettingGet_f(PyMOLGlobals * G, CSetting * set1, CSetting * set2, int index)
{
  if(set1) {
    if(set1->info[index].defined) {
      return (get_f(set1, index));
    }
  }
  if(set2) {
    if(set2->info[index].defined) {
      return (get_f(set2, index));
    }
  }
  return (SettingGetGlobal_f(G, index));
}


/*========================================================================*/
char *SettingGet_s(PyMOLGlobals * G, CSetting * set1, CSetting * set2, int index)
{
  if(set1) {
    if(set1->info[index].defined) {
      return (get_s(set1, index));
    }
  }
  if(set2) {
    if(set2->info[index].defined) {
      return (get_s(set2, index));
    }
  }
  return (SettingGetGlobal_s(G, index));
}


/*========================================================================*/
void SettingGet_3f(PyMOLGlobals * G, CSetting * set1, CSetting * set2, int index,
                   float *value)
{
  float *ptr;
  if(set1) {
    if(set1->info[index].defined) {
      ptr = (float *) (set1->data + set1->info[index].offset);
      copy3f(ptr, value);
      return;
    }
  }
  if(set2) {
    if(set2->info[index].defined) {
      ptr = (float *) (set2->data + set2->info[index].offset);
      copy3f(ptr, value);
      return;
    }
  }
  SettingGetGlobal_3f(G, index, value);
}


/*========================================================================*/
float *SettingGet_3fv(PyMOLGlobals * G, CSetting * set1, CSetting * set2, int index)
{
  if(set1) {
    if(set1->info[index].defined) {
      return (float *) (set1->data + set1->info[index].offset);
    }
  }
  if(set2) {
    if(set2->info[index].defined) {
      return (float *) (set2->data + set2->info[index].offset);
    }
  }
  return (SettingGetGlobal_3fv(G, index));
}


/*========================================================================*/

/*========================================================================*/
int SettingGetIndex(PyMOLGlobals * G, char *name)
{                               /* can be called from any thread state */
#ifdef _PYMOL_NOPY
  /* we're going to need a C-based dictionary and settings name list for this situation */
  return 0;
#else
  PyObject *tmp;
  int unblock;
  int index = -1;

  unblock = PAutoBlock(G);
  if(P_setting) {
    tmp = PyObject_CallMethod(P_setting, "_get_index", "s", name);
    if(tmp) {
      if(PyInt_Check(tmp))
        index = PyInt_AsLong(tmp);
      Py_DECREF(tmp);
    }
  }
  PAutoUnblock(G, unblock);

  return (index);
#endif
}


/*========================================================================*/
int SettingGetName(PyMOLGlobals * G, int index, SettingName name)
{                               /* can be called from any thread state */
#ifdef _PYMOL_NOPY
  /* we're going to need a C-based dictionary and settings name list for this situation */
  name[0] = 0;
  return 0;
#else
  PyObject *tmp;
  int unblock;
  name[0] = 0;
  unblock = PAutoBlock(G);
  if(P_setting) {
    tmp = PyObject_CallMethod(P_setting, "_get_name", "i", index);
    if(tmp) {
      if(PyString_Check(tmp))
        UtilNCopy(name, PyString_AsString(tmp), sizeof(SettingName));
      Py_DECREF(tmp);
    }
  }
  PAutoUnblock(G, unblock);
  return (name[0] != 0);
#endif
}


/*========================================================================*/
void SettingGenerateSideEffects(PyMOLGlobals * G, int index, char *sele, int state)
{
  char all[] = "all";
  char *inv_sele;
  if(!sele) {
    inv_sele = all;
  } else if(sele[0] == 0) {
    inv_sele = all;
  } else {
    inv_sele = sele;
  }
  switch (index) {
  case cSetting_stereo:
    SceneUpdateStereo(G);
    break;
  case cSetting_pickable:
    ExecutiveInvalidateRep(G, inv_sele, cRepAll, cRepInvAll);
    SceneChanged(G);
    break;
  case cSetting_grid_mode:
  case cSetting_grid_slot:
    ExecutiveInvalidateGroups(G, false);
    SceneChanged(G);
    break;
  case cSetting_grid_max:
    SceneChanged(G);
    break;
  case cSetting_defer_builds_mode:
    ExecutiveRebuildAll(G);
    break;
  case cSetting_seq_view:
  case cSetting_seq_view_label_spacing:
  case cSetting_seq_view_label_mode:
  case cSetting_seq_view_label_start:
  case cSetting_seq_view_format:
  case cSetting_seq_view_color:
  case cSetting_seq_view_unaligned_mode:
    SeqChanged(G);
    break;
  case cSetting_seq_view_fill_color:
  case cSetting_seq_view_fill_char:
  case cSetting_seq_view_label_color:
    OrthoDirty(G);
    break;
  case cSetting_show_frame_rate:
    OrthoDirty(G);
    break;
  case cSetting_group_full_member_names:
  case cSetting_group_arrow_prefix:
    OrthoDirty(G);
    break;
  case cSetting_static_singletons:
    SeqChanged(G);
    break;
  case cSetting_seq_view_location:
  case cSetting_seq_view_overlay:
    PParse(G, "cmd.viewport(-1,-1)");
    break;
  case cSetting_stereo_mode:
    SceneUpdateStereoMode(G);
    break;
  case cSetting_dot_lighting:
  case cSetting_mesh_lighting:
  case cSetting_light:
  case cSetting_light2:
  case cSetting_light3:
  case cSetting_light_count:
  case cSetting_fog:
  case cSetting_field_of_view:
  case cSetting_fog_start:
  case cSetting_two_sided_lighting:
  case cSetting_transparency_mode:
  case cSetting_transparency_global_sort:
  case cSetting_dot_normals:
  case cSetting_mesh_normals:
  case cSetting_use_shaders:
    SceneInvalidate(G);
    break;
  case cSetting_stereo_shift:
  case cSetting_stereo_angle:
  case cSetting_stereo_dynamic_strength:
  case cSetting_texture_fonts:
    SceneInvalidate(G);
    break;
  case cSetting_scene_buttons:
  case cSetting_scene_buttons_mode:
    SceneInvalidate(G);
    break;
  case cSetting_dash_round_ends:
  case cSetting_dash_color:
  case cSetting_angle_color:
  case cSetting_dihedral_color:
    SceneInvalidate(G);
    break;
  case cSetting_mouse_selection_mode:
    OrthoDirty(G);
    break;
  case cSetting_internal_gui_control_size:
    WizardRefresh(G);
    OrthoDirty(G);
    break;
  case cSetting_hide_underscore_names:
    ExecutiveInvalidateGroups(G, false);
    OrthoDirty(G);
    break;
  case cSetting_gradient_spacing:
  case cSetting_gradient_max_length:
  case cSetting_gradient_min_length:
  case cSetting_gradient_normal_min_dot:
  case cSetting_gradient_step_size:
  case cSetting_gradient_min_slope:
  case cSetting_gradient_symmetry:
    ExecutiveInvalidateRep(G, inv_sele, cRepMesh, cRepInvRep);
    SceneChanged(G);
    break;
  case cSetting_min_mesh_spacing:
  case cSetting_mesh_grid_max:
  case cSetting_mesh_cutoff:
  case cSetting_mesh_carve_state:
  case cSetting_mesh_carve_cutoff:
  case cSetting_mesh_carve_selection:
  case cSetting_mesh_clear_state:
  case cSetting_mesh_clear_cutoff:
  case cSetting_mesh_clear_selection:
  case cSetting_mesh_mode:
  case cSetting_mesh_type:
  case cSetting_mesh_solvent:
  case cSetting_mesh_quality:
  case cSetting_mesh_skip:
    ExecutiveInvalidateRep(G, inv_sele, cRepMesh, cRepInvRep);
    SceneChanged(G);
    break;
  case cSetting_valence:
  case cSetting_valence_mode:
  case cSetting_valence_size:
  case cSetting_half_bonds:
  case cSetting_stick_transparency:
  case cSetting_line_stick_helper:
  case cSetting_hide_long_bonds:
  case cSetting_line_use_shader:
    ExecutiveInvalidateRep(G, inv_sele, cRepLine, cRepInvRep);
    ExecutiveInvalidateRep(G, inv_sele, cRepCyl, cRepInvRep);
    SceneChanged(G);
    break;
  case cSetting_slice_height_scale:
  case cSetting_slice_height_map:
  case cSetting_slice_grid:
  case cSetting_slice_dynamic_grid:
  case cSetting_slice_dynamic_grid_resolution:
    ExecutiveInvalidateRep(G, inv_sele, cRepSlice, cRepInvRep);
    SceneChanged(G);
    break;
  case cSetting_label_font_id:
  case cSetting_label_size:
    SceneChanged(G);
    break;
  case cSetting_retain_order:
  case cSetting_pdb_hetatm_sort:
  case cSetting_pdb_insertions_go_first:
    ExecutiveSort(G, inv_sele);
    break;
  case cSetting_roving_lines:
  case cSetting_roving_sticks:
  case cSetting_roving_spheres:
  case cSetting_roving_labels:
  case cSetting_roving_selection:
  case cSetting_roving_ribbon:
  case cSetting_roving_cartoon:
  case cSetting_roving_polar_contacts:
  case cSetting_roving_polar_cutoff:
  case cSetting_roving_nonbonded:
  case cSetting_roving_nb_spheres:
  case cSetting_roving_map1_level:
  case cSetting_roving_map2_level:
  case cSetting_roving_map3_level:
  case cSetting_roving_map1_name:
  case cSetting_roving_map2_name:
  case cSetting_roving_map3_name:
  case cSetting_roving_isosurface:
  case cSetting_roving_isomesh:
    SceneRovingChanged(G);
    break;
  case cSetting_roving_byres:
  case cSetting_roving_detail:
    SceneRovingDirty(G);
    break;
  case cSetting_dash_length:
  case cSetting_dash_gap:
  case cSetting_dash_radius:
  case cSetting_dash_width:
  case cSetting_angle_size:
  case cSetting_label_digits:
  case cSetting_label_distance_digits:
  case cSetting_label_angle_digits:
  case cSetting_label_dihedral_digits:
  case cSetting_angle_label_position:
  case cSetting_dihedral_size:
  case cSetting_dihedral_label_position:
    ExecutiveRebuildAllObjectDist(G);
    SceneChanged(G);
    break;
  case cSetting_button_mode:
    EditorMouseInvalid(G);
    OrthoDirty(G);
    break;
  case cSetting_stick_radius:
  case cSetting_stick_h_scale:
    ExecutiveInvalidateRep(G, inv_sele, cRepCyl, cRepInvRep);
    ExecutiveInvalidateRep(G, inv_sele, cRepCartoon, cRepInvRep);       /* base width */
    SceneChanged(G);
    break;
  case cSetting_stick_ball:
  case cSetting_stick_nub:
  case cSetting_stick_ball_ratio:
  case cSetting_stick_ball_color:
  case cSetting_stick_fixed_radius:
  case cSetting_stick_valence_scale:
  case cSetting_stick_quality:
  case cSetting_stick_overlap:
  case cSetting_stick_color:
  case cSetting_stick_use_shader:
    ExecutiveInvalidateRep(G, inv_sele, cRepCyl, cRepInvRep);
    SceneChanged(G);
    break;
  case cSetting_clamp_colors:
  case cSetting_ramp_blend_nearby_colors:
    ExecutiveInvalidateRep(G, inv_sele, cRepAll, cRepInvColor);
    SceneChanged(G);
    break;
  case cSetting_label_color:
  case cSetting_label_outline_color:
  case cSetting_label_position:
    ExecutiveRebuildAllObjectDist(G);
    ExecutiveInvalidateRep(G, inv_sele, cRepLabel, cRepInvRep);
    SceneChanged(G);
    break;
  case cSetting_cartoon_color:
    ExecutiveInvalidateRep(G, inv_sele, cRepCartoon, cRepInvRep);
    SceneChanged(G);
    break;
  case cSetting_ribbon_color:
  case cSetting_ribbon_nucleic_acid_mode:
    ExecutiveInvalidateRep(G, inv_sele, cRepRibbon, cRepInvRep);
    SceneChanged(G);
    break;
  case cSetting_all_states:
    SceneChanged(G);
    break;
  case cSetting_sel_counter:
    break;
  case cSetting_line_width:    /* auto-disable smooth lines if line width > 1 */
    /*    SettingSet(G,cSetting_line_smooth,0);  NO LONGER */
  case cSetting_line_color:
  case cSetting_line_radius:
    ExecutiveInvalidateRep(G, inv_sele, cRepLine, cRepInvRep);
    ExecutiveInvalidateRep(G, inv_sele, cRepNonbonded, cRepInvRep);
    SceneChanged(G);
    break;
  case cSetting_scene_current_name:
    SceneRestartFrameTimer(G);
    break;
  case cSetting_mesh_width:
    ExecutiveInvalidateRep(G, inv_sele, cRepMesh, cRepInvColor);
    SceneChanged(G);
    break;
  case cSetting_ellipsoid_probability:
  case cSetting_ellipsoid_scale:
  case cSetting_ellipsoid_color:
  case cSetting_ellipsoid_transparency:
    ExecutiveInvalidateRep(G, inv_sele, cRepEllipsoid, cRepInvRep);
    SceneChanged(G);
    break;
  case cSetting_ellipsoid_quality:
  case cSetting_cgo_ellipsoid_quality:
    ExecutiveInvalidateRep(G, inv_sele, cRepCGO, cRepInvRep);
    ExecutiveInvalidateRep(G, inv_sele, cRepEllipsoid, cRepInvRep);
    SceneChanged(G);
    break;
  case cSetting_mesh_color:
  case cSetting_mesh_negative_color:
    ExecutiveInvalidateRep(G, inv_sele, cRepMesh, cRepInvColor);
    SceneChanged(G);
    break;
  case cSetting_ray_color_ramps:
    ExecutiveInvalidateRep(G, inv_sele, cRepAll, cRepInvColor);
    SceneChanged(G);
    break;

  case cSetting_cull_spheres:
  case cSetting_sphere_scale:
  case cSetting_sphere_transparency:
  case cSetting_sphere_solvent:
  case cSetting_sphere_mode:
  case cSetting_sphere_point_max_size:
  case cSetting_sphere_point_size:
  case cSetting_sphere_use_shader:
    ExecutiveInvalidateRep(G, inv_sele, cRepSphere, cRepInvRep);
    SceneChanged(G);
    break;
  case cSetting_sphere_quality:
    ExecutiveInvalidateRep(G, inv_sele, cRepCyl, cRepInvRep);
    ExecutiveInvalidateRep(G, inv_sele, cRepNonbondedSphere, cRepInvRep);
    ExecutiveInvalidateRep(G, inv_sele, cRepSphere, cRepInvRep);
    SceneChanged(G);
    break;
  case cSetting_nonbonded_size:
  case cSetting_nonbonded_transparency:
    ExecutiveInvalidateRep(G, inv_sele, cRepNonbonded, cRepInvRep);
    ExecutiveInvalidateRep(G, inv_sele, cRepNonbondedSphere, cRepInvRep);
    SceneChanged(G);
    break;
  case cSetting_mesh_radius:
    ExecutiveInvalidateRep(G, inv_sele, cRepMesh, cRepInvColor);
    SceneChanged(G);
    break;
  case cSetting_surface_negative_color:
  case cSetting_surface_color:
  case cSetting_transparency:
  case cSetting_surface_ramp_above_mode:
    ExecutiveInvalidateRep(G, inv_sele, cRepSurface, cRepInvColor);
    SceneChanged(G);
    break;
  case cSetting_dot_color:
    ExecutiveInvalidateRep(G, inv_sele, cRepDot, cRepInvColor);
    SceneChanged(G);
    break;
  case cSetting_sphere_color:
    ExecutiveInvalidateRep(G, inv_sele, cRepSphere, cRepInvColor);
    SceneChanged(G);
    break;

  case cSetting_surface_quality:
  case cSetting_surface_mode:
  case cSetting_surface_normal:
  case cSetting_surface_type:
  case cSetting_surface_carve_state:
  case cSetting_surface_carve_cutoff:
  case cSetting_surface_carve_selection:
  case cSetting_surface_carve_normal_cutoff:
  case cSetting_surface_clear_state:
  case cSetting_surface_clear_cutoff:
  case cSetting_surface_clear_selection:
  case cSetting_surface_trim_cutoff:
  case cSetting_surface_trim_factor:
  case cSetting_surface_circumscribe:
  case cSetting_surface_solvent:
  case cSetting_surface_proximity:
  case cSetting_surface_cavity_mode:   
  case cSetting_surface_cavity_radius:
  case cSetting_surface_cavity_cutoff:   
  case cSetting_surface_use_shader:
  case cSetting_cavity_cull:
    ExecutiveInvalidateRep(G, inv_sele, cRepSurface, cRepInvRep);
    SceneChanged(G);
    break;
  case cSetting_surface_negative_visible:
    ExecutiveInvalidateRep(G, inv_sele, cRepSurface, cRepInvAll);
    SceneChanged(G);
    break;
  case cSetting_mesh_negative_visible:
    ExecutiveInvalidateRep(G, inv_sele, cRepMesh, cRepInvAll);
    SceneChanged(G);
    break;
  case cSetting_solvent_radius:
    ExecutiveInvalidateRep(G, inv_sele, cRepSurface, cRepInvRep);
    ExecutiveInvalidateRep(G, inv_sele, cRepMesh, cRepInvRep);
    ExecutiveInvalidateRep(G, inv_sele, cRepDot, cRepInvRep);
    SceneChanged(G);
    break;
  case cSetting_trace_atoms_mode:
    ExecutiveInvalidateRep(G, inv_sele, cRepRibbon, cRepInvRep);
    ExecutiveInvalidateRep(G, inv_sele, cRepCartoon, cRepInvRep);
    SceneChanged(G);
    break;
  case cSetting_ribbon_smooth:
  case cSetting_ribbon_power:
  case cSetting_ribbon_power_b:
  case cSetting_ribbon_sampling:
  case cSetting_ribbon_radius:
  case cSetting_ribbon_width:
  case cSetting_ribbon_throw:
  case cSetting_ribbon_trace_atoms:
  case cSetting_ribbon_transparency:
    ExecutiveInvalidateRep(G, inv_sele, cRepRibbon, cRepInvRep);
    SceneChanged(G);
    break;
  case cSetting_draw_mode:
    ExecutiveInvalidateRep(G, inv_sele, cRepCyl, cRepInvRep);
    ExecutiveInvalidateRep(G, inv_sele, cRepSphere, cRepInvRep);
    break;
  case cSetting_cartoon_side_chain_helper:
    ExecutiveInvalidateRep(G, inv_sele, cRepCartoon, cRepInvRep);
    ExecutiveInvalidateRep(G, inv_sele, cRepLine, cRepInvRep);
    ExecutiveInvalidateRep(G, inv_sele, cRepCyl, cRepInvRep);
    ExecutiveInvalidateRep(G, inv_sele, cRepSphere, cRepInvRep);
    SceneChanged(G);
    break;
  case cSetting_ribbon_side_chain_helper:
    ExecutiveInvalidateRep(G, inv_sele, cRepRibbon, cRepInvRep);
    ExecutiveInvalidateRep(G, inv_sele, cRepLine, cRepInvRep);
    ExecutiveInvalidateRep(G, inv_sele, cRepCyl, cRepInvRep);
    ExecutiveInvalidateRep(G, inv_sele, cRepSphere, cRepInvRep);
    SceneChanged(G);
    break;
  case cSetting_cartoon_transparency:
  case cSetting_cartoon_ring_transparency:
  case cSetting_cartoon_trace_atoms:
  case cSetting_cartoon_refine:
  case cSetting_cartoon_nucleic_acid_mode:
  case cSetting_cartoon_nucleic_acid_color:
  case cSetting_cartoon_ring_mode:
  case cSetting_cartoon_ring_finder:
  case cSetting_cartoon_ring_width:
  case cSetting_cartoon_ring_color:
  case cSetting_cartoon_ladder_mode:
  case cSetting_cartoon_ladder_radius:
  case cSetting_cartoon_ladder_color:
  case cSetting_cartoon_sampling:
  case cSetting_cartoon_loop_quality:
  case cSetting_ray_trace_mode:        /* affects loop quality */
  case cSetting_cartoon_loop_radius:
  case cSetting_cartoon_loop_cap:
  case cSetting_cartoon_tube_quality:
  case cSetting_cartoon_tube_radius:
  case cSetting_cartoon_tube_cap:
  case cSetting_cartoon_putty_quality:
  case cSetting_cartoon_putty_radius:
  case cSetting_cartoon_putty_range:
  case cSetting_cartoon_putty_scale_min:
  case cSetting_cartoon_putty_scale_max:
  case cSetting_cartoon_putty_scale_power:
  case cSetting_cartoon_putty_transform:
  case cSetting_cartoon_power:
  case cSetting_cartoon_power_b:
  case cSetting_cartoon_ring_radius:
  case cSetting_cartoon_rect_length:
  case cSetting_cartoon_rect_width:
  case cSetting_cartoon_oval_length:
  case cSetting_cartoon_oval_width:
  case cSetting_cartoon_oval_quality:
  case cSetting_cartoon_round_helices:
  case cSetting_cartoon_flat_sheets:
  case cSetting_cartoon_refine_normals:
  case cSetting_cartoon_smooth_loops:
  case cSetting_cartoon_dumbbell_width:
  case cSetting_cartoon_dumbbell_length:
  case cSetting_cartoon_dumbbell_radius:
  case cSetting_cartoon_fancy_helices:
  case cSetting_cartoon_fancy_sheets:
  case cSetting_cartoon_cylindrical_helices:
  case cSetting_cartoon_refine_tips:
  case cSetting_cartoon_helix_radius:
  case cSetting_cartoon_throw:
  case cSetting_cartoon_debug:
  case cSetting_cartoon_highlight_color:
  case cSetting_cartoon_discrete_colors:
  case cSetting_cartoon_smooth_first:
  case cSetting_cartoon_smooth_last:
  case cSetting_cartoon_smooth_cycles:
  case cSetting_cartoon_flat_cycles:
  case cSetting_cartoon_use_shader:
    ExecutiveInvalidateRep(G, inv_sele, cRepCartoon, cRepInvRep);
    SceneChanged(G);
    break;
  case cSetting_dot_width:
  case cSetting_dot_radius:
  case cSetting_dot_density:
  case cSetting_dot_mode:
  case cSetting_dot_solvent:
  case cSetting_dot_hydrogens:
  case cSetting_trim_dots:
    ExecutiveInvalidateRep(G, inv_sele, cRepDot, cRepInvRep);
    SceneChanged(G);
    break;
  case cSetting_bg_gradient:
    ExecutiveInvalidateRep(G, inv_sele, cRepAll, cRepInvColor);
    SceneChanged(G);
    break;
  case cSetting_bg_rgb_top:
    {
      /* clamp this value */
      float vv[3], *v = SettingGetfv(G, cSetting_bg_rgb_top);
      if((v[0] > 1.0F) || (v[1] > 1.0F) || (v[2] > 1.0F)) {
        vv[0] = v[0] / 255.0F;
        vv[1] = v[1] / 255.0F;
        vv[2] = v[2] / 255.0F;
        SettingSet_3fv(G->Setting, cSetting_bg_rgb_top, vv);
      }
      ColorUpdateFront(G, v);
    }
    ExecutiveInvalidateRep(G, inv_sele, cRepAll, cRepInvColor);
    SceneChanged(G);
    break;
  case cSetting_bg_rgb_bottom:
    {
      /* clamp this value */
      float vv[3], *v = SettingGetfv(G, cSetting_bg_rgb_bottom);
      if((v[0] > 1.0F) || (v[1] > 1.0F) || (v[2] > 1.0F)) {
        vv[0] = v[0] / 255.0F;
        vv[1] = v[1] / 255.0F;
        vv[2] = v[2] / 255.0F;
        SettingSet_3fv(G->Setting, cSetting_bg_rgb_bottom, vv);
      }
      ColorUpdateFront(G, v);
    }
    ExecutiveInvalidateRep(G, inv_sele, cRepAll, cRepInvColor);
    SceneChanged(G);
    break;
  case cSetting_bg_rgb:
    {
      /* clamp this value */
      float vv[3], *v = SettingGetfv(G, cSetting_bg_rgb);
      if((v[0] > 1.0F) || (v[1] > 1.0F) || (v[2] > 1.0F)) {
        vv[0] = v[0] / 255.0F;
        vv[1] = v[1] / 255.0F;
        vv[2] = v[2] / 255.0F;
        SettingSet_3fv(G->Setting, cSetting_bg_rgb, vv);
      }
      ColorUpdateFront(G, v);
    }
    ExecutiveInvalidateRep(G, inv_sele, cRepAll, cRepInvColor);
    SceneChanged(G);
    break;
  case cSetting_line_smooth:
  case cSetting_ortho:
  case cSetting_reflect:
  case cSetting_direct:
  case cSetting_ambient:
  case cSetting_gl_ambient:    /* deprecated */
  case cSetting_specular:
  case cSetting_specular_intensity:
  case cSetting_cgo_line_width:
  case cSetting_selection_width:
  case cSetting_selection_width_scale:
  case cSetting_selection_width_max:
    SceneInvalidate(G);
    break;
  case cSetting_depth_cue:
    SceneInvalidate(G);
    break;
  case cSetting_sculpting:
    OrthoDirty(G);
    break;
  case cSetting_auto_overlay:
    OrthoRemoveAutoOverlay(G); /* always start clean */
    break;
  case cSetting_overlay:
  case cSetting_overlay_lines:
  case cSetting_text:
    OrthoDirty(G);
    break;
  case cSetting_internal_gui_mode:
  case cSetting_internal_gui_width:
  case cSetting_internal_gui:
  case cSetting_internal_feedback:
  case cSetting_mouse_grid:
  case cSetting_movie_panel_row_height:
  case cSetting_movie_panel:
    OrthoCommandIn(G, "viewport");
    break;
  case cSetting_suspend_updates:
    if(!SettingGet(G, cSetting_suspend_updates)) {
      SceneChanged(G);          /* force big update upon resumption */
      OrthoDirty(G);
    }
    break;
  case cSetting_security:
    G->Security = (int) SettingGet(G, cSetting_security);
    break;
  case cSetting_state:
  case cSetting_frame:
    SceneChanged(G);
    break;
  case cSetting_dynamic_width:
  case cSetting_dynamic_width_factor:
  case cSetting_dynamic_width_min:
  case cSetting_dynamic_width_max:
    SceneChanged(G);
    break;
  case cSetting_rock:
  case cSetting_sweep_mode:
  case cSetting_sweep_phase:
  case cSetting_sweep_angle:
  case cSetting_sweep_speed:
    SceneRestartSweepTimer(G);
    break;
  case cSetting_motion_power:
  case cSetting_motion_bias:
  case cSetting_motion_simple:
  case cSetting_motion_linear:
  case cSetting_motion_hand:
  case cSetting_movie_loop:
    if(SettingGetGlobal_b(G, cSetting_movie_auto_interpolate)) {
      ExecutiveMotionReinterpolate(G);
    }
    break;
  case cSetting_volume_bit_depth:
    ExecutiveInvalidateRep(G, inv_sele, cRepVolume, cRepInvColor);
    SceneInvalidate(G);
    break;
  case cSetting_volume_data_range:
    ExecutiveInvalidateRep(G, inv_sele, cRepVolume, cRepInvAll);
    SceneInvalidate(G);
    break;
  case cSetting_volume_layers:
    ExecutiveInvalidateRep(G, inv_sele, cRepVolume, cRepInvAll);
    SceneInvalidate(G);
    break;
  case cSetting_atom_type_format:
    {
      char *setting;
      char *lsetting;
      int i;
      setting = SettingGetGlobal_s(G, cSetting_atom_type_format);
      lsetting = mmalloc(strlen(setting)+1);
      for (i=0; i<=strlen(setting); i++){
	lsetting[i] = tolower(setting[i]);
      }
      if (strcmp(lsetting, "mol2") && 
	  strcmp(lsetting, "sybyl") && 
	  strcmp(lsetting, "macromodel") &&
	  strcmp(lsetting, "mmd")){
	printf("lsetting='%s'\n", lsetting);
	PRINTFB(G, FB_Setting, FB_Warnings)
	  "Setting-Warning: atom_type_format needs to be either mol2/sybyl or macromodel/mmd setting back to default mol2\n"
	  ENDFB(G);
	SettingSet_s(G->Setting, cSetting_atom_type_format, "mol2");	
      } else if (strcmp(setting, lsetting)){
	SettingSet_s(G->Setting, cSetting_atom_type_format, lsetting);
      }
      mfree(lsetting);
      ExecutiveInvalidateRep(G, inv_sele, cRepAll, cRepInvAll);
    }
    break;
  default:
    break;
  }
}


/*========================================================================*/
int SettingSetfv(PyMOLGlobals * G, int index, float *v)
{
  /* Warren, are these side effects still relevant? */

  register CSetting *I = G->Setting;
  int ok = true;
  switch (index) {
  case cSetting_dot_mode:
    SettingSet_f(I, index, v[0]);
    /*I->Setting[index].Value[0]=v[0]; */
    break;
  case cSetting_bg_rgb:
    {
      float vv[3];

      if((v[0] > 1.0F) || (v[1] > 1.0F) || (v[2] > 1.0F)) {
        vv[0] = v[0] / 255.0F;
        vv[1] = v[1] / 255.0F;
        vv[2] = v[2] / 255.0F;
        SettingSet_3fv(I, index, vv);
      } else {
        SettingSet_3fv(I, index, v);
      }
      ColorUpdateFront(G, v);
    }
    ExecutiveInvalidateRep(G, "all", cRepAll, cRepInvColor);
    SceneChanged(G);
    break;
  case cSetting_light:
    SettingSet_3fv(I, index, v);
    SceneInvalidate(G);
    break;
  case cSetting_valence:
    ExecutiveInvalidateRep(G, "all", cRepLine, cRepInvRep);
    SettingSet_f(I, index, v[0]);
    SceneChanged(G);
    break;
  case cSetting_dash_length:
  case cSetting_dash_gap:
    ExecutiveInvalidateRep(G, "all", cRepDash, cRepInvRep);
    SettingSet_f(I, index, v[0]);
    SceneChanged(G);
    break;
  case cSetting_button_mode:
    SettingSet_f(I, index, v[0]);
    OrthoDirty(G);
    break;
  case cSetting_stick_radius:
  case cSetting_stick_quality:
  case cSetting_stick_overlap:
    ExecutiveInvalidateRep(G, "all", cRepCyl, cRepInvRep);
    SettingSet_f(I, index, v[0]);
    /*I->Setting[index].Value[0]=v[0];   */
    SceneChanged(G);
    break;
  case cSetting_label_color:
    ExecutiveInvalidateRep(G, "all", cRepLabel, cRepInvRep);
    SettingSet_f(I, index, v[0]);
    /* I->Setting[index].Value[0]=v[0]; */
    SceneChanged(G);
    break;
  case cSetting_all_states:
    SettingSet_f(I, index, v[0]);
    /* I->Setting[index].Value[0]=v[0];  */
    SceneChanged(G);
    break;
  case cSetting_dot_density:
    SettingSet_f(I, index, v[0]);
    /*I->Setting[index].Value[0]=v[0]; */
    break;
  case cSetting_sel_counter:
    SettingSet_f(I, index, v[0]);
    /*I->Setting[index].Value[0]=v[0]; */
    break;
  case cSetting_ortho:
  case cSetting_gl_ambient:
    SceneInvalidate(G);
    break;
  case cSetting_overlay:
  case cSetting_text:
    OrthoDirty(G);
  default:
    ok = SettingSet_f(I, index, v[0]);
    /*I->Setting[index].Value[0]=v[0]; */
    break;
  }
  return (ok);
}


/*========================================================================*/
int SettingSetGlobal_b(PyMOLGlobals * G, int index, int value)
{
  return (SettingSet_b(G->Setting, index, value));
}


/*========================================================================*/
int SettingSetGlobal_i(PyMOLGlobals * G, int index, int value)
{
  return (SettingSet_i(G->Setting, index, value));
}


/*========================================================================*/
int SettingSetGlobal_f(PyMOLGlobals * G, int index, float value)
{
  return (SettingSet_f(G->Setting, index, value));
}


/*========================================================================*/
int SettingSetGlobal_3f(PyMOLGlobals * G, int index, float value1, float value2,
                        float value3)
{
  return (SettingSet_3f(G->Setting, index, value1, value2, value3));
}


/*========================================================================*/
int SettingSet(PyMOLGlobals * G, int index, float v)
{
  return (SettingSetfv(G, index, &v));
}


/*========================================================================*/
int SettingSetNamed(PyMOLGlobals * G, char *name, char *value)
{
  int ok = true;
  int index = SettingGetIndex(G, name);
  float v, vv[3];
  SettingName realName;
  char buffer[1024] = "";
  if(index >= 0) {
    SettingGetName(G, index, realName);
    switch (index) {
    case cSetting_dot_mode:
      if(strcmp(value, "molecular") == 0) {
        v = 0.0;
        SettingSetfv(G, index, &v);
        sprintf(buffer, " Setting: %s set to %s\n", realName, value);
      } else if(strcmp(value, "solvent_accessible") == 0) {
        v = 1.0;
        SettingSetfv(G, index, &v);
        sprintf(buffer, " Setting: %s set to %s\n", realName, value);
      } else if(sscanf(value, "%f", &v) == 1) {
        SettingSetfv(G, index, &v);
        sprintf(buffer, " Setting: %s set to %s\n", realName, value);
      }
      break;
    case cSetting_bg_rgb:
    case cSetting_light:
      if(sscanf(value, "%f%f%f", vv, vv + 1, vv + 2) == 3) {
        SettingSetfv(G, index, vv);
        sprintf(buffer, " Setting: %s set to %5.3f %8.3f %8.3f\n", realName,
                *vv, *(vv + 1), *(vv + 2));
      }
      break;
    case cSetting_dot_density:
      sscanf(value, "%f", &v);
      SettingSetfv(G, index, &v);
      sprintf(buffer, " Setting: %s set to %d\n", realName, (int) v);
      break;
    case cSetting_text:
    case cSetting_overlay:
    case cSetting_sel_counter:
    case cSetting_dist_counter:
      sscanf(value, "%f", &v);
      SettingSetfv(G, index, &v);
      break;
    case cSetting_line_width:  /* auto-disable smooth lines if line width > 1 */
    case cSetting_mesh_width:
      sscanf(value, "%f", &v);
      SettingSetfv(G, index, &v);
      sprintf(buffer, " Setting: %s set to %5.3f\n", realName, v);
      SceneInvalidate(G);
      break;
    default:
      sscanf(value, "%f", &v);
      SettingSetfv(G, index, &v);
      sprintf(buffer, " Setting: %s set to %5.3f\n", realName, v);
      break;
    }
  } else {
    PRINTFB(G, FB_Setting, FB_Warnings)
      " Error: Non-Existent Settin\n" ENDFB(G);
    ok = false;
  }
  if(buffer[0]) {
    PRINTFB(G, FB_Setting, FB_Actions)
      "%s", buffer ENDFB(G);
  }
  return (ok);
}


/*========================================================================*/
float SettingGetNamed(PyMOLGlobals * G, char *name)
{
  return (SettingGet(G, SettingGetIndex(G, name)));
}


/*========================================================================*/
float SettingGet(PyMOLGlobals * G, int index)
{
  return (SettingGetGlobal_f(G, index));
}


/*========================================================================*/
float *SettingGetfv(PyMOLGlobals * G, int index)
{
  return (SettingGetGlobal_3fv(G, index));
}


/*========================================================================*/
void SettingFreeGlobal(PyMOLGlobals * G)
{
  register CSetting *I = G->Setting;
  SettingUniqueFree(G);
  SettingPurge(I);
  if(G->Default) {
    SettingPurge(G->Default);
    FreeP(G->Default);
  }
  FreeP(G->Setting);
}


/*========================================================================*/
void SettingInitGlobal(PyMOLGlobals * G, int alloc, int reset_gui, int use_default)
{
  register CSetting *I = G->Setting;

  /* use function pointers to prevent the compiler from inlining every
     call in this block (a waste of RAM and time) */

  int (*set_f) (CSetting * I, int index, float valueI) = SettingSet_f;
  int (*set_i) (CSetting * I, int index, int value) = SettingSet_i;
  int (*set_3f) (CSetting * I, int index, float value1, float value2, float value3) =
    SettingSet_3f;
  int (*set_b) (CSetting * I, int index, int value) = SettingSet_b;
  int (*set_color) (CSetting * I, int index, char *value) = SettingSet_color;
  int (*set_s) (CSetting * I, int index, char *value) = SettingSet_s;

  if(alloc || !I) {
    I = (G->Setting = Calloc(CSetting, 1));
    SettingUniqueInit(G);
    SettingInit(G, I);
  }

  if(G->Default && use_default) {

    SettingCopyAll(G, G->Default, G->Setting);

  } else {

    set_f(I, cSetting_bonding_vdw_cutoff, 0.2F);

    set_f(I, cSetting_min_mesh_spacing, 0.6F);

    set_i(I, cSetting_dot_density, 2);

    set_i(I, cSetting_dot_mode, 0);

    set_f(I, cSetting_solvent_radius, 1.4F);

    set_i(I, cSetting_sel_counter, 0);

    set_3f(I, cSetting_bg_rgb, 0.0F, 0.0F, 0.0F);

    set_f(I, cSetting_ambient, 0.14F);

    set_f(I, cSetting_direct, 0.45F);

    set_f(I, cSetting_reflect, 0.45F);

    set_3f(I, cSetting_light, -0.4F, -0.4F, -1.0F);

    set_i(I, cSetting_antialias, 1);

    set_i(I, cSetting_cavity_cull, 10);

    set_f(I, cSetting_gl_ambient, 0.12F);       /* no longer effective */

    set_b(I, cSetting_single_image, 0);

    set_f(I, cSetting_movie_delay, 30.0F);

    set_f(I, cSetting_ribbon_power, 2.0F);

    set_f(I, cSetting_ribbon_power_b, 0.5F);

    set_i(I, cSetting_ribbon_sampling, 1);

    set_f(I, cSetting_ribbon_radius, 0.0F);

    set_f(I, cSetting_stick_radius, 0.25F);

    set_i(I, cSetting_hash_max, 100);

    set_b(I, cSetting_ortho, 0);

    set_f(I, cSetting_power, 1.0F);

    set_f(I, cSetting_spec_reflect, -1.0F);

    set_f(I, cSetting_spec_power, -1.0F);

    set_f(I, cSetting_sweep_angle, 20.0F);

    set_f(I, cSetting_sweep_speed, 0.75F);

    set_b(I, cSetting_dot_hydrogens, 1);

    set_f(I, cSetting_dot_radius, 0.0F);

    set_b(I, cSetting_ray_trace_frames, 0);

    set_b(I, cSetting_cache_frames, 0);

    set_b(I, cSetting_trim_dots, 1);

    set_i(I, cSetting_cull_spheres, 0);

    set_f(I, cSetting_test1, 3.0F);

    set_f(I, cSetting_test2, -0.5F);

    set_f(I, cSetting_surface_best, 0.25F);

    set_f(I, cSetting_surface_normal, 0.5F);

    set_i(I, cSetting_surface_quality, 0);

    set_b(I, cSetting_surface_proximity, 1);

    set_f(I, cSetting_stereo_angle, 2.1F);

    set_f(I, cSetting_stereo_shift, 2.0F);

    set_b(I, cSetting_line_smooth, 1);

    set_f(I, cSetting_line_width, 1.49F);       /* under 1.5F to retain SGI antialiasing */

    set_b(I, cSetting_half_bonds, 0);

    set_i(I, cSetting_stick_quality, 8);

    set_f(I, cSetting_stick_overlap, 0.2F);

    set_f(I, cSetting_stick_nub, 0.7F);

    set_b(I, cSetting_all_states, 0);

    set_b(I, cSetting_pickable, 1);

    set_i(I, cSetting_sphere_quality, 1);

    set_b(I, cSetting_auto_show_lines, G->Option->sphere_mode < 0);

    set_f(I, cSetting_fast_idle, 10000.0F);     /* 1/100th of a sec. */

    set_f(I, cSetting_no_idle, 2000.0F);        /* 1/500th of a sec. */

    set_f(I, cSetting_slow_idle, 40000.0F);     /* 1/25th of a sec. */

    set_f(I, cSetting_idle_delay, 1.5F);

    set_f(I, cSetting_rock_delay, 30.0F);

    set_i(I, cSetting_dist_counter, 0);

    set_f(I, cSetting_dash_length, 0.15F);

    set_f(I, cSetting_dash_gap, 0.45F);

    set_i(I, cSetting_auto_zoom, G->Option->zoom_mode);

    set_i(I, cSetting_overlay, 0);

    set_b(I, cSetting_text, 0);

    set_i(I, cSetting_button_mode, 0);

    set_b(I, cSetting_valence, 0);

    set_f(I, cSetting_nonbonded_size, 0.25F);

    set_color(I, cSetting_label_color, "-6");

    set_f(I, cSetting_ray_trace_fog, -1.0F);

    set_f(I, cSetting_spheroid_scale, 1.0F);

    set_f(I, cSetting_ray_trace_fog_start, -1.0F);

    set_f(I, cSetting_spheroid_smooth, 1.1F);

    set_f(I, cSetting_spheroid_fill, 1.30F);

    set_b(I, cSetting_auto_show_nonbonded, G->Option->sphere_mode < 0);

    set_f(I, cSetting_mesh_radius, 0.000F);

#ifdef WIN32
    /* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
    set_b(I, cSetting_cache_display, 0);
    /* END PROPRIETARY CODE SEGMENT */
#else
    set_b(I, cSetting_cache_display, 1);
#endif

    set_b(I, cSetting_normal_workaround, 0);

    set_b(I, cSetting_backface_cull, 0); 

    set_f(I, cSetting_gamma, 1.0F);

    set_f(I, cSetting_dot_width, 2.0F);

    set_b(I, cSetting_auto_show_selections, 1);

    set_b(I, cSetting_auto_hide_selections, 1);

    set_f(I, cSetting_selection_width, 3.0F);

    set_f(I, cSetting_selection_overlay, 1.0F);

    set_b(I, cSetting_static_singletons, 1);

    set_i(I, cSetting_max_triangles, 1000000);  /* no longer used */

    set_b(I, cSetting_depth_cue, 1);

    set_f(I, cSetting_specular, 1.0F);

    set_f(I, cSetting_shininess, 55.0F);

    set_f(I, cSetting_fog, 1.0F);

    set_b(I, cSetting_isomesh_auto_state, 0);   /* no longer necessary? */

    set_f(I, cSetting_mesh_width, 1.0F);

    set_i(I, cSetting_cartoon_sampling, 7);

    set_f(I, cSetting_cartoon_loop_radius, 0.2F);

    set_f(I, cSetting_cartoon_loop_quality, 6.0F);

    set_f(I, cSetting_cartoon_power, 2.0F);

    set_f(I, cSetting_cartoon_power_b, 0.52F);

    set_f(I, cSetting_cartoon_rect_length, 1.40F);

    set_f(I, cSetting_cartoon_rect_width, 0.4F);

    if(reset_gui) {
      set_i(I, cSetting_internal_gui_width, cOrthoRightSceneMargin);

      set_b(I, cSetting_internal_gui, 1);
    }

    set_f(I, cSetting_cartoon_oval_length, 1.35F);

    set_f(I, cSetting_cartoon_oval_width, 0.25F);

    set_f(I, cSetting_cartoon_oval_quality, 10.0F);

    set_f(I, cSetting_cartoon_tube_radius, 0.5F);

    set_f(I, cSetting_cartoon_tube_quality, 9.0F);

    set_i(I, cSetting_cartoon_debug, 0);

    set_f(I, cSetting_ribbon_width, 3.0F);

    set_f(I, cSetting_dash_width, 2.5F);

    set_f(I, cSetting_dash_radius, 0.00F);

    set_f(I, cSetting_cgo_ray_width_scale, -0.15F);

    set_f(I, cSetting_line_radius, 0.0F);

    set_b(I, cSetting_cartoon_round_helices, 1);

    set_i(I, cSetting_cartoon_refine_normals, -1);

    set_b(I, cSetting_cartoon_flat_sheets, 1);

    set_b(I, cSetting_cartoon_smooth_loops, 0);

    set_f(I, cSetting_cartoon_dumbbell_length, 1.60F);

    set_f(I, cSetting_cartoon_dumbbell_width, 0.17F);

    set_f(I, cSetting_cartoon_dumbbell_radius, 0.16F);

    set_b(I, cSetting_cartoon_fancy_helices, 0);

    set_b(I, cSetting_cartoon_fancy_sheets, 1);

    set_b(I, cSetting_ignore_pdb_segi, 0);

    set_f(I, cSetting_ribbon_throw, 1.35F);

    set_f(I, cSetting_cartoon_throw, 1.35F);

    set_i(I, cSetting_cartoon_refine, 5);

    set_i(I, cSetting_cartoon_refine_tips, 10);

    set_b(I, cSetting_cartoon_discrete_colors, 0);

    set_b(I, cSetting_normalize_ccp4_maps, 1);

    set_f(I, cSetting_surface_poor, 0.85F);

    set_i(I, cSetting_internal_feedback, G->Option->internal_feedback);

    set_f(I, cSetting_cgo_line_width, 1.00F);

    set_f(I, cSetting_cgo_line_radius, -0.05F);

    set_i(I, cSetting_logging, 0);      /* 0 = off, 1 = regular (PML), 2 = python (PYM) */

    set_b(I, cSetting_robust_logs, 0);

    set_b(I, cSetting_log_box_selections, 1);

    set_b(I, cSetting_log_conformations, 1);

    set_f(I, cSetting_valence_size, 0.060F);

    set_f(I, cSetting_surface_miserable, 2.0F);

    set_i(I, cSetting_ray_opaque_background, -1);

    set_f(I, cSetting_transparency, 0.0F);

    set_i(I, cSetting_ray_texture, 0);

    set_3f(I, cSetting_ray_texture_settings, 0.1F, 5.0F, 1.0F);

    set_b(I, cSetting_suspend_updates, 0);

    set_b(I, cSetting_full_screen, 0);

    set_i(I, cSetting_surface_mode, 0); /* by flag is the default */

    set_color(I, cSetting_surface_color, "-1"); /* use atom colors by default */

    set_i(I, cSetting_mesh_mode, 0);    /* by flag is the default */

    set_color(I, cSetting_mesh_color, "-1");    /* use atom colors by default */

    set_b(I, cSetting_auto_indicate_flags, 0);

    set_i(I, cSetting_surface_debug, 0);

    set_f(I, cSetting_ray_improve_shadows, 0.1F);

    set_b(I, cSetting_smooth_color_triangle, 0);

    set_i(I, cSetting_ray_default_renderer, 0);

    set_f(I, cSetting_field_of_view, 20.0F);

    set_f(I, cSetting_reflect_power, 1.0F);

    set_b(I, cSetting_preserve_chempy_ids, 0);

    set_f(I, cSetting_sphere_scale, 1.0F);

    set_i(I, cSetting_two_sided_lighting, -1);

    set_f(I, cSetting_secondary_structure, 2.0F);       /* unused? */

    set_b(I, cSetting_auto_remove_hydrogens, 0);

    set_b(I, cSetting_raise_exceptions, 1);

    set_b(I, cSetting_stop_on_exceptions, 0);

    set_b(I, cSetting_sculpting, 0);

    set_b(I, cSetting_auto_sculpt, 0);

    set_f(I, cSetting_sculpt_vdw_scale, 0.97F);

    set_f(I, cSetting_sculpt_vdw_scale14, 0.90F);       /* 0.915 */

    set_f(I, cSetting_sculpt_vdw_weight, 1.0F);

    set_f(I, cSetting_sculpt_vdw_weight14, 0.2F);       /* 0.33 */

    set_f(I, cSetting_sculpt_bond_weight, 2.25F);

    set_f(I, cSetting_sculpt_angl_weight, 1.0F);

    set_f(I, cSetting_sculpt_pyra_weight, 1.0F);

    set_f(I, cSetting_sculpt_plan_weight, 1.0F);

    set_i(I, cSetting_sculpting_cycles, 10);

    set_f(I, cSetting_sphere_transparency, 0.0F);

    set_color(I, cSetting_sphere_color, "-1");  /* use atom colors by default */

    set_i(I, cSetting_sculpt_field_mask, 0x1FF);        /* all terms */

    set_f(I, cSetting_sculpt_hb_overlap, 1.0F);

    set_f(I, cSetting_sculpt_hb_overlap_base, 0.35F);

    set_b(I, cSetting_legacy_vdw_radii, 0);

    set_b(I, cSetting_sculpt_memory, 1);

    set_i(I, cSetting_connect_mode, 0);

    set_b(I, cSetting_cartoon_cylindrical_helices, 0);

    set_f(I, cSetting_cartoon_helix_radius, 2.25F);

    set_f(I, cSetting_connect_cutoff, 0.35F);

    set_b(I, cSetting_save_pdb_ss, 0);

    set_f(I, cSetting_sculpt_line_weight, 1.0F);

    set_i(I, cSetting_fit_iterations, 1000);

    set_f(I, cSetting_fit_tolerance, 0.0000001F);

    set_s(I, cSetting_batch_prefix, "tmp_pymol");

    if(!G->Option->stereo_mode) {
      if(G->StereoCapable || G->Option->blue_line) {
        set_i(I, cSetting_stereo_mode, 1);      /* quadbuffer if we can */
      } else {
        set_i(I, cSetting_stereo_mode, 2);      /* otherwise crosseye by default */
      }
    } else {
      set_i(I, cSetting_stereo_mode, G->Option->stereo_mode);
    }

    set_i(I, cSetting_cgo_sphere_quality, 1);

    set_b(I, cSetting_pdb_literal_names, 0);

    set_b(I, cSetting_wrap_output, 0);

    set_f(I, cSetting_fog_start, 0.45F);

    set_i(I, cSetting_frame, 1);

    set_i(I, cSetting_state, 1);

    set_b(I, cSetting_ray_shadows, 1);

    set_i(I, cSetting_ribbon_trace_atoms, 0);

    set_i(I, cSetting_security, 1);

    set_f(I, cSetting_stick_transparency, 0.0F);

    set_b(I, cSetting_ray_transparency_shadows, 1);

    set_i(I, cSetting_session_version_check, 0);

    set_f(I, cSetting_ray_transparency_specular, 0.6F);

    set_b(I, cSetting_stereo_double_pump_mono, 0);

    set_b(I, cSetting_sphere_solvent, 0);

    set_i(I, cSetting_mesh_quality, 2);

    set_i(I, cSetting_mesh_solvent, 0);

    set_b(I, cSetting_dot_solvent, 0);

    set_f(I, cSetting_ray_shadow_fudge, 0.001F);

    set_f(I, cSetting_ray_triangle_fudge, 0.0000001F);

    set_i(I, cSetting_debug_pick, 0);

    set_color(I, cSetting_dot_color, "-1");     /* use atom colors by default */

    set_f(I, cSetting_mouse_limit, 100.0F);

    set_f(I, cSetting_mouse_scale, 1.3F);

    set_i(I, cSetting_transparency_mode, 2);

    set_b(I, cSetting_clamp_colors, 1);

    set_f(I, cSetting_pymol_space_max_red, 0.90F);

    set_f(I, cSetting_pymol_space_max_green, 0.75F);

    set_f(I, cSetting_pymol_space_max_blue, 0.90F);

    set_f(I, cSetting_pymol_space_min_factor, 0.15F);

    set_b(I, cSetting_roving_origin, 1);

    set_f(I, cSetting_roving_sticks, 6.0F);

    set_f(I, cSetting_roving_lines, 10.0F);

    set_f(I, cSetting_roving_spheres, 0.0F);

    set_f(I, cSetting_roving_labels, 0.0F);

    set_f(I, cSetting_roving_delay, 0.2F);

    set_s(I, cSetting_roving_selection, "all");

    set_b(I, cSetting_roving_byres, 1);

    set_f(I, cSetting_roving_ribbon, -7.0F);

    set_f(I, cSetting_roving_cartoon, 0.0F);

    set_f(I, cSetting_roving_polar_contacts, 7.0F);

    set_f(I, cSetting_roving_polar_cutoff, 3.31F);

    set_f(I, cSetting_roving_nonbonded, 0.0F);

    set_i(I, cSetting_float_labels, 0);

    set_b(I, cSetting_roving_detail, 0);

    set_f(I, cSetting_roving_nb_spheres, 8.0F);

    set_color(I, cSetting_ribbon_color, "-1");  /* use atom colors by default */

    set_color(I, cSetting_cartoon_color, "-1"); /* use atom colors by default */

    set_i(I, cSetting_ribbon_smooth, 0);

    set_b(I, cSetting_auto_color, 1);

    set_i(I, cSetting_auto_color_next, 0);

    set_color(I, cSetting_ray_interior_color, "-1");    /* object color */

    set_color(I, cSetting_cartoon_highlight_color, "-1");       /* no color */

    /* In order to get electrostatic potentials in kT from the Coulomb equation... 

       PyMOL charges: Q, in e
       PyMOL distances: r, in Angstrom
       Coulomb Constant: K = 8.987552e9 ((J*m)/(C^2))
       Angstrom Convertor: 1 A = 1e-10 m
       Coulomb Convertor: 1 e = 1.60217733e-19 C
       Angstrom Convertor: 1 A = 10e-10 m
       Dielectric Constant: D (unitless)

       ePot = (KQ)/(Dr) = 

       8.987552e9 J*m     1.6021773e-19 C   1.6021773e-19 C      1 A        Q  
       ---------------- * --------------- * --------------- * ---------- * --- =
       C^2              1 e               1 e            1e-10 m     Dr

       2.3070795237e-18 J*A      Q
       = ---------------------- * ---
       e^2               Dr

       Boltzmann Constant: k = 1.380658e-23 (J/K)
       Temperature: 300 Kelvin

       kT = 1.380658e-23 * 300 = 4.141974e-21 J

       2.3070795237e-18 J*A         1 kT             Q
       ePot = --------------------- * ------------------ * ---
       e^2              4.141974e-21 J       Dr

       557.00000 kT*A    Q
       ePot = -------------- * --- which will give kT/e units when applied
       e^2           Dr
     */

    set_f(I, cSetting_coulomb_units_factor, 557.00000F);

    set_f(I, cSetting_coulomb_dielectric, 2.0F);

    set_b(I, cSetting_ray_interior_shadows, 0);

    set_i(I, cSetting_ray_interior_texture, -1);

    set_s(I, cSetting_roving_map1_name, "");

    set_s(I, cSetting_roving_map2_name, "");

    set_s(I, cSetting_roving_map3_name, "");

    set_f(I, cSetting_roving_map1_level, 1.0F);

    set_f(I, cSetting_roving_map2_level, 2.0F);

    set_f(I, cSetting_roving_map3_level, 3.0F);

    set_f(I, cSetting_roving_isomesh, 8.0F);

    set_f(I, cSetting_roving_isosurface, 0.0F);

    set_f(I, cSetting_scenes_changed, 1.0F);

    set_f(I, cSetting_gaussian_b_adjust, 0.0F);

    set_b(I, cSetting_pdb_standard_order, 1);

    set_i(I, cSetting_cartoon_smooth_first, 1);
    set_i(I, cSetting_cartoon_smooth_last, 1);
    set_i(I, cSetting_cartoon_smooth_cycles, 2);
    set_i(I, cSetting_cartoon_flat_cycles, 4);

#ifdef WIN32
#ifdef _PYMOL_ACTIVEX
    set_i(I, cSetting_max_threads, 1);
#else
    /* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
    {
      SYSTEM_INFO SysInfo;
      GetSystemInfo(&SysInfo);
      {
        DWORD count = SysInfo.dwNumberOfProcessors;
        if(count > 1) {
          set_i(I, cSetting_max_threads, count);
        } else {
          set_i(I, cSetting_max_threads, 1);
        }
      }
    }
    /* END PROPRIETARY CODE SEGMENT */
#endif
#else
    set_i(I, cSetting_max_threads, 1);
#endif

    set_i(I, cSetting_show_progress, 1);

    set_i(I, cSetting_use_display_lists, 0);    /* don't make this default
                                                   until we have a way of
                                                   reusing display list 
                                                   identifiers */

    set_i(I, cSetting_cache_memory, 0); /* doesn't seem to do any good :( */

    set_b(I, cSetting_simplify_display_lists, 0);

    set_i(I, cSetting_retain_order, 0);

    set_i(I, cSetting_pdb_hetatm_sort, 0);

    set_i(I, cSetting_pdb_use_ter_records, 1);

    set_i(I, cSetting_cartoon_trace_atoms, 0);

    set_i(I, cSetting_ray_oversample_cutoff, 120);

    /* note that this setting is ad-hoc and calibrated such that a
       gaussian_resolution of 2.0 returns maps with the straight atomic
       scattering factors (unblurred).  At resolution of 4.0, they are
       blurred 2X, 8.0:4X, and so forth.... */
    set_f(I, cSetting_gaussian_resolution, 2.0F);

    set_f(I, cSetting_gaussian_b_floor, 0.0F);

    set_i(I, cSetting_sculpt_nb_interval, 17);

    set_f(I, cSetting_sculpt_tors_weight, 0.05F);

    set_f(I, cSetting_sculpt_tors_tolerance, 0.05F);

    set_b(I, cSetting_stick_ball, false);

    set_f(I, cSetting_stick_ball_ratio, 1.0F);

    set_b(I, cSetting_stick_fixed_radius, false);

    set_f(I, cSetting_cartoon_transparency, 0.0F);

    set_b(I, cSetting_dash_round_ends, 1);

    set_f(I, cSetting_h_bond_max_angle, 63.0F);

    set_f(I, cSetting_h_bond_cutoff_center, 3.6F);

    set_f(I, cSetting_h_bond_cutoff_edge, 3.2F);

    set_f(I, cSetting_h_bond_power_a, 1.6F);

    set_f(I, cSetting_h_bond_power_b, 5.0F);

    set_f(I, cSetting_h_bond_cone, 180.0F);

    set_f(I, cSetting_ss_helix_psi_target, -48.0F);
    set_f(I, cSetting_ss_helix_psi_include, 55.0F);     /* 30 */
    set_f(I, cSetting_ss_helix_psi_exclude, 85.0F);

    set_f(I, cSetting_ss_helix_phi_target, -57.0F);
    set_f(I, cSetting_ss_helix_phi_include, 55.0F);
    set_f(I, cSetting_ss_helix_phi_exclude, 85.0F);

    set_f(I, cSetting_ss_strand_psi_target, 124.0F);
    set_f(I, cSetting_ss_strand_psi_include, 40.0F);
    set_f(I, cSetting_ss_strand_psi_exclude, 90.0F);    /* 80 */

    set_f(I, cSetting_ss_strand_phi_target, -129.0F);
    set_f(I, cSetting_ss_strand_phi_include, 40.0F);
    set_f(I, cSetting_ss_strand_phi_exclude, 100.0F);

    set_b(I, cSetting_movie_loop, 1);

    set_b(I, cSetting_pdb_retain_ids, 0);

    set_b(I, cSetting_pdb_no_end_record, 0);

    set_f(I, cSetting_cgo_dot_width, 2.0F);
    set_f(I, cSetting_cgo_dot_radius, -1.0F);
    set_b(I, cSetting_defer_updates, 0);
    set_b(I, cSetting_normalize_o_maps, 1);
    set_b(I, cSetting_swap_dsn6_bytes, 1);
    set_b(I, cSetting_pdb_insertions_go_first, 0);
    set_b(I, cSetting_roving_origin_z, 1);
    set_f(I, cSetting_roving_origin_z_cushion, 3.0F);
    set_f(I, cSetting_specular_intensity, 0.5F);
    set_i(I, cSetting_overlay_lines, 5);
    set_f(I, cSetting_ray_transparency_spec_cut, 0.9F);
    set_b(I, cSetting_internal_prompt, 1);
    set_b(I, cSetting_normalize_grd_maps, 0);

    set_b(I, cSetting_ray_blend_colors, 0);
    set_f(I, cSetting_ray_blend_red, 0.17F);
    set_f(I, cSetting_ray_blend_green, 0.25F);
    set_f(I, cSetting_ray_blend_blue, 0.14F);
    set_f(I, cSetting_png_screen_gamma, 2.4F);
    set_f(I, cSetting_png_file_gamma, 1.0F);
    set_b(I, cSetting_editor_label_fragments, 0);
    set_i(I, cSetting_internal_gui_control_size, 18);
    set_b(I, cSetting_auto_dss, 1);
    set_i(I, cSetting_transparency_picking_mode, 2);    /* auto */
    set_b(I, cSetting_virtual_trackball, 1);
    set_i(I, cSetting_transparency_picking_mode, 2);    /* auto */
    set_i(I, cSetting_pdb_reformat_names_mode, 0);      /*
                                                           0 = no reformatting, 
                                                           1 = pdb compliant,
                                                           2 = amber compliant,
                                                           3 = pdb I/O, but iupac inside
                                                         */
    set_f(I, cSetting_ray_pixel_scale, 1.30F);
#ifdef _PYMOL_FREETYPE
    set_i(I, cSetting_label_font_id, 5);
#else
    set_i(I, cSetting_label_font_id, 0);
#endif
    set_b(I, cSetting_pdb_conect_all, 0);
    set_s(I, cSetting_button_mode_name, "");
    set_i(I, cSetting_surface_type, 0);
    set_b(I, cSetting_dot_normals, 1);
    set_b(I, cSetting_session_migration, 1);
    set_b(I, cSetting_mesh_normals, 1);
    set_i(I, cSetting_mesh_type, 0);    /* 0 = lines, 1 = points, 2 = solid, 3 = gradient */

    set_b(I, cSetting_dot_lighting, 1);
    set_b(I, cSetting_mesh_lighting, 0);
    set_b(I, cSetting_surface_solvent, 0);
    set_i(I, cSetting_triangle_max_passes, 5);
    set_f(I, cSetting_ray_interior_reflect, 0.4F);
    set_i(I, cSetting_internal_gui_mode, 0);
    set_s(I, cSetting_surface_carve_selection, "");
    set_i(I, cSetting_surface_carve_state, 0);
    set_f(I, cSetting_surface_carve_cutoff, 0.0F);
    set_s(I, cSetting_surface_clear_selection, "");
    set_i(I, cSetting_surface_clear_state, 0);
    set_f(I, cSetting_surface_clear_cutoff, 0.0F);
    set_f(I, cSetting_surface_trim_cutoff, 0.2F);
    set_f(I, cSetting_surface_trim_factor, 2.0F);
    set_i(I, cSetting_ray_max_passes, 25);
    set_b(I, cSetting_active_selections, true);
    set_f(I, cSetting_ray_transparency_contrast, 1.0F);
    set_b(I, cSetting_seq_view, 0);
    set_i(I, cSetting_mouse_selection_mode, 1);
    set_i(I, cSetting_seq_view_label_spacing, 5);
    set_i(I, cSetting_seq_view_label_start, 1);
    set_i(I, cSetting_seq_view_format, 0);
    set_i(I, cSetting_seq_view_location, 0);
    set_b(I, cSetting_seq_view_overlay, 0);
    set_b(I, cSetting_auto_classify_atoms, 1);
    set_i(I, cSetting_cartoon_nucleic_acid_mode, 4);
    set_color(I, cSetting_seq_view_color, "-1");
    set_i(I, cSetting_seq_view_label_mode, 2);
    set_i(I, cSetting_surface_ramp_above_mode, 0);
    set_b(I, cSetting_stereo, 0);
    set_i(I, cSetting_wizard_prompt_mode, 1);
    set_f(I, cSetting_coulomb_cutoff, 10.0F);
    set_b(I, cSetting_slice_track_camera, 0);
    set_f(I, cSetting_slice_height_scale, 1.0F);
    set_b(I, cSetting_slice_height_map, 0);
    set_f(I, cSetting_slice_grid, 0.3F);
    set_b(I, cSetting_slice_dynamic_grid, 0);
    set_f(I, cSetting_slice_dynamic_grid_resolution, 3.0F);
    set_b(I, cSetting_pdb_insure_orthogonal, 1);
    set_f(I, cSetting_ray_direct_shade, 0.0F);  /* only meaningful with one light source */
    set_color(I, cSetting_stick_color, "-1");
    set_f(I, cSetting_cartoon_putty_radius, 0.40F);
    set_f(I, cSetting_cartoon_putty_quality, 11.0F);
    set_f(I, cSetting_cartoon_putty_scale_min, 0.6F);
    set_f(I, cSetting_cartoon_putty_scale_max, 4.0F);
    set_f(I, cSetting_cartoon_putty_scale_power, 1.5F);
    set_f(I, cSetting_cartoon_putty_range, 2.0F);
    set_b(I, cSetting_cartoon_side_chain_helper, 0);
    set_b(I, cSetting_surface_optimize_subsets, 1);
    set_i(I, cSetting_multiplex, -1);
    if(G->Option->multisample > 1)
      set_b(I, cSetting_texture_fonts, 1);
    else
      set_b(I, cSetting_texture_fonts, 0);
    set_b(I, cSetting_pqr_workarounds, 1);
    set_b(I, cSetting_animation, 1);
    set_f(I, cSetting_animation_duration, 0.75F);
    set_i(I, cSetting_scene_animation, -1);
    set_b(I, cSetting_line_stick_helper, 1);
    set_i(I, cSetting_ray_orthoscopic, -1);
    set_i(I, cSetting_ribbon_side_chain_helper, 0);
    set_f(I, cSetting_selection_width_max, 10.0F);
    set_f(I, cSetting_selection_width_scale, 2.0F);
    set_s(I, cSetting_scene_current_name, "");
    set_b(I, cSetting_presentation, G->Option->presentation);
    set_i(I, cSetting_presentation_mode, 1);
    set_b(I, cSetting_pdb_truncate_residue_name, false);
    set_b(I, cSetting_scene_loop, 0);
    set_i(I, cSetting_sweep_mode, 0);
    set_f(I, cSetting_sweep_phase, 0.0F);
    set_b(I, cSetting_scene_restart_movie_delay, 1);
    set_b(I, cSetting_mouse_restart_movie_delay, 0);
    set_f(I, cSetting_angle_size, 0.6666F);
    set_f(I, cSetting_angle_label_position, 0.5);
    set_f(I, cSetting_dihedral_size, 0.6666F);
    set_f(I, cSetting_dihedral_label_position, 1.2F);
    set_i(I, cSetting_defer_builds_mode, G->Option->defer_builds_mode);
    set_b(I, cSetting_seq_view_discrete_by_state, 1);
    set_f(I, cSetting_scene_animation_duration, 2.25F);
    set_s(I, cSetting_wildcard, "*");
    set_s(I, cSetting_atom_name_wildcard, "");
    set_b(I, cSetting_ignore_case, 1);
    set_b(I, cSetting_presentation_auto_quit, !G->Option->no_quit);
    set_b(I, cSetting_editor_auto_dihedral, 1);
    set_b(I, cSetting_presentation_auto_start, 1);
    set_b(I, cSetting_validate_object_names, 1);
    set_b(I, cSetting_unused_boolean_def_true, 1);
    set_b(I, cSetting_auto_show_spheres, G->Option->sphere_mode >= 0);
    set_i(I, cSetting_sphere_mode, G->Option->sphere_mode);
    set_f(I, cSetting_sphere_point_max_size, 18.0);
    set_f(I, cSetting_sphere_point_size, 1.0);
    set_b(I, cSetting_pdb_honor_model_number, false);
    set_b(I, cSetting_rank_assisted_sorts, true);
    set_i(I, cSetting_ribbon_nucleic_acid_mode, 0);
    set_i(I, cSetting_cartoon_ring_mode, 0);
    set_f(I, cSetting_cartoon_ring_width, 0.125F);
    set_color(I, cSetting_cartoon_ring_color, "-1");
    set_i(I, cSetting_cartoon_ring_finder, 1);
    set_i(I, cSetting_cartoon_tube_cap, 2);
    set_i(I, cSetting_cartoon_loop_cap, 1);
    set_i(I, cSetting_nvidia_bugs, 0);
    set_f(I, cSetting_image_dots_per_inch, 0.0F);
    /* default is to leave it unspecified in PNG file */
    set_b(I, cSetting_opaque_background, 1);
    set_b(I, cSetting_draw_frames, 0);
    set_b(I, cSetting_show_alpha_checker, 1);
    set_i(I, cSetting_matrix_mode, -1); /* -1: automatic behavior based on implied intent
                                            0: coordinates (pre-1.0 legacy default mode)
                                            1: per-object matrices (TTTs: version 1.0 default mode?)
                                            2: per-state matrices (partially implemented)
                                            3: per-group matrices (may come in the future) */
    set_b(I, cSetting_editor_auto_origin, 1);
    set_s(I, cSetting_session_file, "");
    set_f(I, cSetting_cgo_transparency, 0.0F);
    set_b(I, cSetting_legacy_mouse_zoom, 0);
    set_b(I, cSetting_auto_number_selections, 0);
    set_i(I, cSetting_sculpt_vdw_vis_mode, 0);
    set_f(I, cSetting_sculpt_vdw_vis_min, -0.1F);
    set_f(I, cSetting_sculpt_vdw_vis_mid, 0.1F);
    set_f(I, cSetting_sculpt_vdw_vis_max, 0.3F);
    set_i(I, cSetting_cartoon_ladder_mode, 1);
    set_f(I, cSetting_cartoon_ladder_radius, 0.25F);
    set_color(I, cSetting_cartoon_ladder_color, "-1");
    set_color(I, cSetting_cartoon_nucleic_acid_color, "-1");
    set_f(I, cSetting_cartoon_ring_transparency, -1.0F);
    set_f(I, cSetting_label_size, 14.0F);
    set_f(I, cSetting_spec_direct, 0.0F);
    set_i(I, cSetting_light_count, 2);
    set_3f(I, cSetting_light2, -0.55F, -0.7F, 0.15F);
    set_3f(I, cSetting_light3, 0.3F, -0.6F, -0.2F);
    set_b(I, cSetting_hide_underscore_names, 1);
    set_b(I, cSetting_selection_round_points, 0);
    set_i(I, cSetting_distance_exclusion, 5);
    set_i(I, cSetting_h_bond_exclusion, 3);
    set_i(I, cSetting_label_shadow_mode, 0);
    set_3f(I, cSetting_light4, -1.2F, 0.3F, -0.2F);
    set_3f(I, cSetting_light5, 0.3F, 0.6F, -0.75F);
    set_3f(I, cSetting_light6, -0.3F, 0.5F, 0.0F);
    set_3f(I, cSetting_light7, 0.9F, -0.1F, -0.15F);
    set_color(I, cSetting_label_outline_color, "-1");
    set_i(I, cSetting_ray_trace_mode, 0);
    set_f(I, cSetting_ray_trace_gain, 0.12F);
    set_b(I, cSetting_selection_visible_only, 0);
    set_3f(I, cSetting_label_position, 0.0F, 0.0F, 1.75F);
    set_f(I, cSetting_ray_trace_depth_factor, 0.1F);
    set_f(I, cSetting_ray_trace_slope_factor, 0.6F);
    set_f(I, cSetting_ray_trace_disco_factor, 0.05F);
    set_f(I, cSetting_ray_shadow_decay_factor, 0.0F);
    set_i(I, cSetting_ray_interior_mode, 0);
    set_f(I, cSetting_ray_legacy_lighting, 0.0F);
    set_b(I, cSetting_sculpt_auto_center, 0);
    set_i(I, cSetting_pdb_discrete_chains, -1);
    set_i(I, cSetting_pdb_unbond_cations, 1);
    set_f(I, cSetting_sculpt_tri_scale, 1.025F);        /* allow for some play here... */
    set_f(I, cSetting_sculpt_tri_weight, 1.0F);
    set_i(I, cSetting_sculpt_tri_min, 2);
    set_i(I, cSetting_sculpt_tri_max, 18);
    set_i(I, cSetting_sculpt_tri_mode, 0);
    set_s(I, cSetting_pdb_echo_tags, "HEADER, TITLE, COMPND");
    set_b(I, cSetting_connect_bonded, 0);
    set_f(I, cSetting_spec_direct_power, 55.0F);
    set_3f(I, cSetting_light8, 1.3F, 2.0F, 0.8F);
    set_3f(I, cSetting_light9, -1.7F, -0.5F, 1.2F);
    set_f(I, cSetting_ray_shadow_decay_range, 1.8F);
    set_i(I, cSetting_spec_count, -1);
    set_f(I, cSetting_sculpt_min_scale, 0.975F);
    set_f(I, cSetting_sculpt_min_weight, 0.75F);
    set_f(I, cSetting_sculpt_min_min, 4.0F);
    set_f(I, cSetting_sculpt_min_max, 12.0F);
    set_f(I, cSetting_sculpt_max_scale, 1.025F);
    set_f(I, cSetting_sculpt_max_weight, 0.75F);
    set_f(I, cSetting_sculpt_max_min, 4.0F);
    set_f(I, cSetting_sculpt_max_max, 12.0F);
    set_i(I, cSetting_surface_circumscribe, -1);
    set_f(I, cSetting_sculpt_avd_weight, 4.0F);
    set_f(I, cSetting_sculpt_avd_gap, -1.0F);
    set_f(I, cSetting_sculpt_avd_range, -1.0F);
    set_i(I, cSetting_sculpt_avd_excl, 7);
    set_b(I, cSetting_async_builds, 0);
    set_s(I, cSetting_fetch_path, ".");
    set_s(I, cSetting_fetch_host, "pdb");
    set_f(I, cSetting_cartoon_ring_radius, -1.0F);
    set_b(I, cSetting_ray_color_ramps, 0);
    set_f(I, cSetting_ray_hint_camera, 2.15F);
    set_f(I, cSetting_ray_hint_shadow, 0.65F);
    set_f(I, cSetting_stick_valence_scale, 1.0F);
    set_s(I, cSetting_seq_view_alignment, "");
    set_i(I, cSetting_seq_view_unaligned_mode, 0);
    set_color(I, cSetting_seq_view_unaligned_color, "-1");
    set_s(I, cSetting_seq_view_fill_char, "-");
    set_color(I, cSetting_seq_view_fill_color, "104");  /* grey50 */
    set_color(I, cSetting_seq_view_label_color, "white");       /* grey50 */
    set_f(I, cSetting_surface_carve_normal_cutoff, -1.0F);
    set_i(I, cSetting_trace_atoms_mode, 5);
    set_b(I, cSetting_session_changed, 0);
    set_b(I, cSetting_ray_clip_shadows, 0);
    set_f(I, cSetting_mouse_wheel_scale, 0.5F);
    set_f(I, cSetting_nonbonded_transparency, 0.0F);
    set_b(I, cSetting_ray_spec_local, 0);
    set_color(I, cSetting_line_color, "-1");
    set_f(I, cSetting_ray_label_specular, 1.0F);
    set_i(I, cSetting_mesh_skip, 0);
    set_i(I, cSetting_label_digits, 1);
    set_i(I, cSetting_label_distance_digits, -1);
    set_i(I, cSetting_label_angle_digits, -1);
    set_i(I, cSetting_label_dihedral_digits, -1);
    set_s(I, cSetting_label_anchor, "CA");
    set_b(I, cSetting_surface_negative_visible, 0);
    set_color(I, cSetting_surface_negative_color, "grey50");
    set_b(I, cSetting_mesh_negative_visible, 0);
    set_color(I, cSetting_mesh_negative_color, "grey30");
    set_i(I, cSetting_group_auto_mode, 1);
    set_i(I, cSetting_group_full_member_names, 0);
    set_f(I, cSetting_gradient_max_length, 100.0F);
    set_f(I, cSetting_gradient_min_length, 2.0F);
    set_f(I, cSetting_gradient_min_slope, 0.00001F);
    set_f(I, cSetting_gradient_normal_min_dot, 0.70F);
    set_f(I, cSetting_gradient_step_size, 0.25F);
    set_i(I, cSetting_gradient_spacing, 3);
    set_f(I, cSetting_gradient_symmetry, 0.0F);
    set_color(I, cSetting_ray_trace_color, "-6");
    set_b(I, cSetting_group_arrow_prefix, 0);
    set_b(I, cSetting_suppress_hidden, true);
    set_b(I, cSetting_session_compression, 0);
    set_f(I, cSetting_movie_fps, 30.0);
    set_f(I, cSetting_ray_transparency_oblique, 0.0F);
    set_f(I, cSetting_ray_trace_trans_cutoff, 0.05F);
    set_f(I, cSetting_ray_trace_persist_cutoff, 0.10F);
    set_f(I, cSetting_ray_transparency_oblique_power, 1.0F);
    set_f(I, cSetting_ray_scatter, 0.0F);
    set_b(I, cSetting_h_bond_from_proton, 1);
    set_b(I, cSetting_auto_copy_images, 0);
    set_i(I, cSetting_moe_separate_chains, -1);
    set_b(I, cSetting_transparency_global_sort, 0);
    set_b(I, cSetting_hide_long_bonds, 0);
    set_b(I, cSetting_auto_rename_duplicate_objects, 0);        /* to do */
    set_b(I, cSetting_pdb_hetatm_guess_valences, 1);
    set_i(I, cSetting_ellipsoid_quality, 1);
    set_i(I, cSetting_cgo_ellipsoid_quality, -1);
    set_b(I, cSetting_movie_animate_by_frame, 0);
    set_b(I, cSetting_ramp_blend_nearby_colors, 0);
    set_i(I, cSetting_auto_defer_builds, 500);  /* 500 or more states, then automatically defer builds */
    set_f(I, cSetting_ellipsoid_probability, 0.5F);
    set_f(I, cSetting_ellipsoid_scale, 1.0F);
    set_color(I, cSetting_ellipsoid_color, "-1");
    set_f(I, cSetting_ellipsoid_transparency, 0.0F);
    set_i(I, cSetting_movie_rock, -1);
    set_i(I, cSetting_cache_mode, 0);
    set_color(I, cSetting_dash_color, "-1");
    set_color(I, cSetting_angle_color, "-1");
    set_color(I, cSetting_dihedral_color, "-1");
    set_i(I, cSetting_grid_mode, 0);
    set_i(I, cSetting_cache_max, 25000000);     /* default: ~100 MB cache */
    set_i(I, cSetting_grid_slot, -1);
    set_i(I, cSetting_grid_max, -1);
    set_i(I, cSetting_cartoon_putty_transform, cPuttyTransformNormalizedNonlinear);
    set_b(I, cSetting_rock, 0);
    set_i(I, cSetting_cone_quality, 18);
    set_b(I, cSetting_pdb_formal_charges, 1);
    set_i(I, cSetting_ati_bugs, 0);
    set_i(I, cSetting_geometry_export_mode, 0);
    set_b(I, cSetting_mouse_grid, 1);
    set_s(I, cSetting_mesh_carve_selection, "");
    set_i(I, cSetting_mesh_carve_state, 0);
    set_f(I, cSetting_mesh_carve_cutoff, 0.0F);
    set_s(I, cSetting_mesh_clear_selection, "");
    set_i(I, cSetting_mesh_clear_state, 0);
    set_f(I, cSetting_mesh_clear_cutoff, 0.0F);
    set_f(I, cSetting_mesh_cutoff, 0.0F);
    set_i(I, cSetting_mesh_grid_max, 80);
    set_i(I, cSetting_session_cache_optimize, 0);
    set_f(I, cSetting_sdof_drag_scale, 0.5F);
    set_i(I, cSetting_scene_buttons_mode, 1);
    set_b(I, cSetting_scene_buttons, 0);
    set_b(I, cSetting_map_auto_expand_sym, 1);
    set_b(I, cSetting_image_copy_always, 0);
    set_i(I, cSetting_max_ups, 0);
    set_i(I, cSetting_auto_overlay, 0);
    set_color(I, cSetting_stick_ball_color, "-1");
    set_f(I, cSetting_stick_h_scale, 0.4F);
    set_f(I, cSetting_sculpt_pyra_inv_weight, 10.0F);
    set_b(I, cSetting_keep_alive, 0);
    set_i(I, cSetting_fit_kabsch, 0);
    set_f(I, cSetting_stereo_dynamic_strength, 0.5F);
    set_b(I, cSetting_dynamic_width, 1);
    set_f(I, cSetting_dynamic_width_factor, 0.06);
    set_f(I, cSetting_dynamic_width_min, 0.75);
    set_f(I, cSetting_dynamic_width_max, 2.5);
    set_i(I, cSetting_draw_mode, 0);
    set_i(I, cSetting_clean_electro_mode, 1);
    set_i(I, cSetting_valence_mode, 1);
    set_b(I, cSetting_show_frame_rate, 0);
    set_i(I, cSetting_movie_panel, 1);
    set_f(I, cSetting_mouse_z_scale,1.0);
    set_b(I, cSetting_movie_auto_store, -1);
    set_b(I, cSetting_movie_auto_interpolate, 1);
    set_i(I, cSetting_movie_panel_row_height, 15);
    set_i(I, cSetting_movie_quality, 60);
    set_i(I, cSetting_scene_frame_mode,-1);
    set_i(I, cSetting_surface_cavity_mode,0);
    set_f(I, cSetting_surface_cavity_radius, 7.0F);
    set_f(I, cSetting_surface_cavity_cutoff, -3.0F);
    set_f(I, cSetting_motion_power, 0.0F);
    set_f(I, cSetting_motion_bias, -1.0F);
    set_i(I, cSetting_motion_simple, 0);
    set_f(I, cSetting_motion_linear, 0.0F);
    set_i(I, cSetting_motion_hand, 1);
    set_b(I, cSetting_pdb_ignore_conect, 0);
    set_b(I, cSetting_editor_bond_cycle_mode, 1); /* >0 -> include aromatic */
    set_b(I, cSetting_dynamic_measures, 1);
    set_f(I, cSetting_neighbor_cutoff, 3.5F);
    set_f(I, cSetting_heavy_neighbor_cutoff, 3.5F);
    set_f(I, cSetting_polar_neighbor_cutoff, 3.5F);
    set_f(I, cSetting_surface_residue_cutoff, 2.5F);
    set_b(I, cSetting_surface_use_shader, 1);
    set_b(I, cSetting_cartoon_use_shader, 1);
    set_b(I, cSetting_stick_use_shader, 1);
    set_b(I, cSetting_line_use_shader, 1);
    set_b(I, cSetting_sphere_use_shader, 1);
    set_b(I, cSetting_use_shaders, 0);  /* disable by default until optimized shaders present; doesn't effect vol */
    set_s(I, cSetting_shader_path, "data/shaders");
    set_i(I, cSetting_volume_bit_depth, 8);
    set_color(I, cSetting_volume_color, "-1");
    set_f(I, cSetting_volume_layers, 256);
    set_f(I, cSetting_volume_data_range, 5.0);
    set_i(I, cSetting_auto_defer_atom_count, 0);
    set_s(I, cSetting_default_refmac_names, "FWT PHWT DELFWT PHDELWT");
    set_s(I, cSetting_default_phenix_names, "2FOFCWT PH2FOFCWT FOFCWT PHFOFCWT");
    set_s(I, cSetting_default_phenix_no_fill_names, "2FOFCWT_no_fil PH2FOFCWT_no_fill None None");
    set_s(I, cSetting_default_buster_names, "2FOFCWT PH2FOFCWT FOFCWT PHFOFCWT");
    set_s(I, cSetting_default_fofc_map_rep, "volume");
    set_s(I, cSetting_default_2fofc_map_rep, "volume");
    set_s(I, cSetting_atom_type_format, "mol2");

    set_b(I, cSetting_autoclose_dialogs, 1);

    set_b(I, cSetting_bg_gradient, 0);
    set_3f(I, cSetting_bg_rgb_top, 0.0F, 0.0F, 0.3F);
    set_3f(I, cSetting_bg_rgb_bottom, 0.2F, 0.2F, 0.5F);

    set_b(I, cSetting_ray_volume, 0);
    
    set_f(I, cSetting_ribbon_transparency, 0.0F);

    set_i(I, cSetting_state_counter_mode, -1);

  }
}
