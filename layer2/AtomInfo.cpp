
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
#include <utility>
#include <algorithm>

#include"os_python.h"

#include"os_predef.h"
#include"os_std.h"

#include"AtomInfo.h"
#include"Word.h"
#include"MemoryDebug.h"
#include"Err.h"
#include"Feedback.h"
#include"Util.h"
#include"Util2.h"
#include"Color.h"
#include"PConv.h"
#include"Ortho.h"
#include"OVOneToAny.h"
#include"OVContext.h"
#include"PyMOLObject.h"
#include"Setting.h"
#include"Executive.h"
#include "Lex.h"
#include "pymol/zstring_view.h"

#include <map>

struct _CAtomInfo {
  int NColor, CColor, DColor, HColor, OColor, SColor;
  int BrColor, ClColor, FColor, IColor;
  int PColor, MgColor, MnColor, NaColor, KColor, CaColor;
  int CuColor, FeColor, ZnColor;
  int SeColor;
  int DefaultColor;
  int NextUniqueID;
  OVOneToAny *ActiveIDs;
};

void AtomInfoCleanAtomName(char *name)
{
  char *p = name, *q = name;
  while(*p) {
    if((((*p) >= '0') && ((*p) <= '9')) ||
       (((*p) >= 'a') && ((*p) <= 'z')) ||
       (((*p) >= 'A') && ((*p) <= 'Z')) ||
       ((*p) == '.') ||
       ((*p) == '_') || ((*p) == '+') || ((*p) == '\'') || ((*p) == '*')) {
      *q++ = *p;
    }
    p++;
  }
  *q = 0;
}

#ifndef _PYMOL_NOPY
int AtomInfoSetSettingFromPyObject(PyMOLGlobals * G, AtomInfoType *ai, int setting_id, PyObject * val){
  if (val == Py_None)
    val = NULL;

  if (!val) {
    if (!ai->has_setting)
      return true;
  }

  AtomInfoCheckUniqueID(G, ai);
  ai->has_setting = true;

  return SettingUniqueSetPyObject(G, ai->unique_id, setting_id, val);
}
#endif

PyObject *SettingGetIfDefinedPyObject(PyMOLGlobals * G, AtomInfoType * ai, int setting_id) {
  if(ai->has_setting) {
    return SettingUniqueGetPyObject(G, ai->unique_id, setting_id);
  }
  return NULL;
}

static int AtomInfoPrimeUniqueIDs(PyMOLGlobals * G)
{
  CAtomInfo *I = G->AtomInfo;
  if(!I->ActiveIDs) {
    OVContext *C = G->Context;
    I->ActiveIDs = OVOneToAny_New(C->heap);
  }
  return (I->ActiveIDs != NULL);
}

int AtomInfoReserveUniqueID(PyMOLGlobals * G, int unique_id)
{
  CAtomInfo *I = G->AtomInfo;
  if(!I->ActiveIDs)
    AtomInfoPrimeUniqueIDs(G);
  if(I->ActiveIDs)
    return (OVreturn_IS_OK(OVOneToAny_SetKey(I->ActiveIDs, unique_id, 1)));
  return 0;
}

int AtomInfoIsUniqueIDActive(PyMOLGlobals * G, int unique_id)
{
  CAtomInfo *I = G->AtomInfo;
  if(!I->ActiveIDs)
    return 0;
  else
    return (OVreturn_IS_OK(OVOneToAny_GetKey(I->ActiveIDs, unique_id)));
  return 0;
}

int AtomInfoGetNewUniqueID(PyMOLGlobals * G)
{
  CAtomInfo *I = G->AtomInfo;
  int result = 0;
  AtomInfoPrimeUniqueIDs(G);
  if(I->ActiveIDs) {
    while(1) {
      result = I->NextUniqueID++;
      if(result) {              /* skip zero */
        if(OVOneToAny_GetKey(I->ActiveIDs, result).status == OVstatus_NOT_FOUND) {
          if(OVreturn_IS_ERROR(OVOneToAny_SetKey(I->ActiveIDs, result, 1)))
            result = 0;
          break;
        }
      }
    }
  }
  ExecutiveUniqueIDAtomDictInvalidate(G);
  return result;
}

int AtomInfoCheckUniqueID(PyMOLGlobals * G, AtomInfoType * ai)
{
  if(!ai->unique_id)
    ai->unique_id = AtomInfoGetNewUniqueID(G);
  return ai->unique_id;
}

int AtomInfoCheckUniqueBondID(PyMOLGlobals * G, BondType * bi)
{
  if(!bi->unique_id)
    bi->unique_id = AtomInfoGetNewUniqueID(G);
  return bi->unique_id;
}

void BondTypeInit(BondType *bt){
  bt->unique_id = 0;
  bt->has_setting = 0;
  bt->symop_2 = pymol::SymOp();
}

/**
 * Set atom indices and bond order, and invalidate all other fields.
 */
void BondTypeInit2(BondType *bond, int i1, int i2, int order)
{
  BondTypeInit(bond);
  bond->index[0] = i1;
  bond->index[1] = i2;
  bond->order = order;
}

int AtomInfoInit(PyMOLGlobals * G)
{
  CAtomInfo *I = NULL;
  if((I = (G->AtomInfo = pymol::calloc<CAtomInfo>(1)))) {
    AtomInfoPrimeColors(G);
    I->NextUniqueID = 1;
    return 1;
  } else
    return 0;
}

void AtomInfoFree(PyMOLGlobals * G)
{
  CAtomInfo *I = G->AtomInfo;
  OVOneToAny_DEL_AUTO_NULL(I->ActiveIDs);
  FreeP(G->AtomInfo);
}


/*========================================================================*/
void AtomInfoGetPDB3LetHydroName(PyMOLGlobals * G, const char *resn, const char *iname, char *oname)
{
  oname[0] = ' ';
  strcpy(oname + 1, iname);

  switch (resn[0]) {
  case 'A':
    switch (resn[1]) {
    case 'L':
      if(resn[2] == 'A')
        if((iname[0] == 'H') &&
           (iname[1] == 'B') && ((iname[2] >= '0') && (iname[2] <= '9'))) {
          oname[0] = iname[2];
          oname[1] = iname[0];
          oname[2] = iname[1];
          oname[3] = 0;
        }
      break;
    case 'R':
      if(resn[2] == 'G')
        if((iname[0] == 'H') &&
           ((iname[1] == 'B') || (iname[1] == 'G') || (iname[1] == 'D')) &&
           ((iname[2] >= '0') && (iname[2] <= '9'))) {
          oname[0] = iname[2];
          oname[1] = iname[0];
          oname[2] = iname[1];
          oname[3] = 0;
        }
      break;
    case 'S':
      switch (resn[2]) {
      case 'P':
        if((iname[0] == 'H') &&
           (iname[1] == 'B') && ((iname[2] >= '0') && (iname[2] <= '9'))) {
          oname[0] = iname[2];
          oname[1] = iname[0];
          oname[2] = iname[1];
          oname[3] = 0;
        }
        break;
      case 'N':
        if((iname[0] == 'H') &&
           (iname[1] == 'B') && ((iname[2] >= '0') && (iname[2] <= '9'))) {
          oname[0] = iname[2];
          oname[1] = iname[0];
          oname[2] = iname[1];
          oname[3] = 0;
        }
        break;
      }
      break;
    }
    break;
  case 'C':
    switch (resn[1]) {
    case 'Y':
      switch (resn[2]) {
      case 'S':
      case 'X':
        if((iname[0] == 'H') &&
           (iname[1] == 'B') && ((iname[2] >= '0') && (iname[2] <= '9'))) {
          oname[0] = iname[2];
          oname[1] = iname[0];
          oname[2] = iname[1];
          oname[3] = 0;
        }
        break;
      }
      break;
    }
    break;
  case 'G':
    switch (resn[1]) {
    case 'L':
      switch (resn[2]) {
      case 'N':
        if((iname[0] == 'H') &&
           ((iname[1] == 'B') || (iname[1] == 'G')) &&
           ((iname[2] >= '0') && (iname[2] <= '9'))) {
          oname[0] = iname[2];
          oname[1] = iname[0];
          oname[2] = iname[1];
          oname[3] = 0;
        }
        break;
      case 'U':
        if((iname[0] == 'H') &&
           ((iname[1] == 'B') || (iname[1] == 'G')) &&
           ((iname[2] >= '0') && (iname[2] <= '9'))) {
          oname[0] = iname[2];
          oname[1] = iname[0];
          oname[2] = iname[1];
          oname[3] = 0;
        }
        break;
      case 'Y':
        if((iname[0] == 'H') &&
           (iname[1] == 'A') && ((iname[2] >= '0') && (iname[2] <= '9'))) {
          oname[0] = iname[2];
          oname[1] = iname[0];
          oname[2] = iname[1];
          oname[3] = 0;
        }
        break;
      }
    }
    break;
  case 'H':
    switch (resn[1]) {
    case 'I':
      switch (resn[2]) {
      case 'S':
      case 'D':
      case 'E':
      case 'P':
        if((iname[0] == 'H') &&
           (iname[1] == 'B') && ((iname[2] >= '0') && (iname[2] <= '9'))) {
          oname[0] = iname[2];
          oname[1] = iname[0];
          oname[2] = iname[1];
          oname[3] = 0;
        }
        break;
      }
      break;
    case 'O':
      switch (resn[2]) {
      case 'H':
        break;
      }
      break;
    case '2':
      switch (resn[2]) {
      case 'O':
        break;
      }
      break;
    }
  case 'I':
    switch (resn[1]) {
    case 'L':
      switch (resn[2]) {
      case 'E':
        break;
      }
    }
    break;
  case 'L':
    switch (resn[1]) {
    case 'E':
      switch (resn[2]) {
      case 'U':
        if((iname[0] == 'H') &&
           (iname[1] == 'B') && ((iname[2] >= '0') && (iname[2] <= '9'))) {
          oname[0] = iname[2];
          oname[1] = iname[0];
          oname[2] = iname[1];
          oname[3] = 0;
        }
        break;
      }
      break;
    case 'Y':
      switch (resn[2]) {
      case 'S':
        if((iname[0] == 'H') &&
           ((iname[1] == 'B') || (iname[1] == 'G') || (iname[1] == 'D')
            || (iname[1] == 'E') || (iname[1] == 'Z')) && ((iname[2] >= '0')
                                                           && (iname[2] <= '9'))) {
          oname[0] = iname[2];
          oname[1] = iname[0];
          oname[2] = iname[1];
          oname[3] = 0;
        }
        break;
      }
      break;
    }
    break;
  case 'M':
    switch (resn[1]) {
    case 'E':
      switch (resn[2]) {
      case 'T':
        if((iname[0] == 'H') &&
           ((iname[1] == 'B') || (iname[1] == 'G') || (iname[1] == 'E')) &&
           ((iname[2] >= '0') && (iname[2] <= '9'))) {
          oname[0] = iname[2];
          oname[1] = iname[0];
          oname[2] = iname[1];
          oname[3] = 0;
        }
        break;
      }
    }
    break;
  case 'P':
    switch (resn[1]) {
    case 'H':
      switch (resn[2]) {
      case 'E':
        if((iname[0] == 'H') &&
           (iname[1] == 'B') && ((iname[2] >= '0') && (iname[2] <= '9'))) {
          oname[0] = iname[2];
          oname[1] = iname[0];
          oname[2] = iname[1];
          oname[3] = 0;
        }
        break;
      }
      break;
    case 'R':
      switch (resn[2]) {
      case 'O':
        if((iname[0] == 'H') &&
           ((iname[1] == 'B') || (iname[1] == 'G') || (iname[1] == 'D')) &&
           ((iname[2] >= '0') && (iname[2] <= '9'))) {
          oname[0] = iname[2];
          oname[1] = iname[0];
          oname[2] = iname[1];
          oname[3] = 0;
        }
        break;
      }
      break;
    }
    break;
  case 'S':
    switch (resn[1]) {
    case 'E':
      switch (resn[2]) {
      case 'R':
        if((iname[0] == 'H') &&
           (iname[1] == 'B') && ((iname[2] >= '0') && (iname[2] <= '9'))) {
          oname[0] = iname[2];
          oname[1] = iname[0];
          oname[2] = iname[1];
          oname[3] = 0;
        }
        break;
      }
      break;
    }
    break;
  case 'T':
    switch (resn[1]) {
    case 'H':
      switch (resn[2]) {
      case 'R':
        break;
      }
      break;
    case 'I':
      switch (resn[2]) {
      case 'P':
        break;
      }
      break;
    case 'R':
      switch (resn[2]) {
      case 'P':
        if((iname[0] == 'H') &&
           (iname[1] == 'B') && ((iname[2] >= '0') && (iname[2] <= '9'))) {
          oname[0] = iname[2];
          oname[1] = iname[0];
          oname[2] = iname[1];
          oname[3] = 0;
        }
        break;
      }
      break;
    case 'Y':
      switch (resn[2]) {
      case 'R':
        if((iname[0] == 'H') &&
           (iname[1] == 'B') && ((iname[2] >= '0') && (iname[2] <= '9'))) {
          oname[0] = iname[2];
          oname[1] = iname[0];
          oname[2] = iname[1];
          oname[3] = 0;
        }
        break;
      }
      break;
    }
    break;
  case 'V':
    switch (resn[1]) {
    case 'A':
      switch (resn[2]) {
      case 'L':
        break;
      }
      break;
    }
    break;
  case 'W':
    switch (resn[1]) {
    case 'A':
      switch (resn[2]) {
      case 'T':
        break;
      }
      break;
    }
    break;
  }
}

int AtomInfoKnownWaterResName(PyMOLGlobals * G, const char *resn)
{
  switch (resn[0]) {
  case 'H':
    switch (resn[1]) {
    case 'O':
      switch (resn[2]) {
      case 'H':
        return true;
        break;
      }
      break;
    case '2':
      switch (resn[2]) {
      case 'O':
        return true;
        break;
      }
      break;
    }
  case 'D':
    switch (resn[1]) {
    case 'O':
      switch (resn[2]) {
      case 'D':
        return true;
        break;
      }
      break;
    }
    break;
  case 'T':
    switch (resn[1]) {
    case 'I':
    case '3': // T3P
    case '4': // T4P
      switch (resn[2]) {
      case 'P':
        return true;
        break;
      }
      break;
    }
    break;
  case 'W':
    switch (resn[1]) {
    case 'A':
      switch (resn[2]) {
      case 'T':
        return true;
        break;
      }
      break;
    }
    break;
  case 'S':
    switch (resn[1]) {
    case 'O':
      switch (resn[2]) {
      case 'L':
        return true;
        break;
      }
      break;
    case 'P':
      switch (resn[2]) {
      case 'C':
        return true;
        break;
      }
      break;
    }
    break;
  }
  return false;
}

int AtomInfoKnownPolymerResName(const char *resn)
{
  return AtomInfoKnownProteinResName(resn) ||
         AtomInfoKnownNucleicResName(resn);
}

int AtomInfoKnownNucleicResName(const char *resn)
{
  if (resn[0] == 'D') {
    // Deoxy ribonucleotide
    ++resn;
  }

  switch (resn[0]) {
    case 'A':
    case 'C':
    case 'G':
    case 'I':
    case 'T':
    case 'U':
      if (!resn[1])
        return true;
  }

  return false;
}

int AtomInfoKnownProteinResName(const char *resn)
{
  switch (resn[0]) {
  case 'A':
    switch (resn[1]) {
    case 'L':
      switch (resn[2]) {
      case 'A': /* ALA */
        return true;
        break;
      }
      break;
    case 'R':
      switch (resn[2]) {
      case 'G': /* ARG */
        return true;
        break;
      }
      break;
    case 'S':
      switch (resn[2]) {
      case 'P': /* ASP */
        return true;
        break;
      case 'N': /* ASN */
        return true;
        break;
      }
      break;
    }
    break;
  case 'C':
    switch (resn[1]) {
    case 'Y':
      switch (resn[2]) {
      case 'S': /* CYS */
      case 'X': /* CYX */
        return true;
        break;
      }
      break;
    }
    break;
  case 'G':
    switch (resn[1]) {
    case 'L':
      switch (resn[2]) {
      case 'N': /* GLN */
        return true;
        break;
      case 'U': /* GLU */
        return true;
        break;
      case 'Y': /* GLY */
        return true;
        break;
      }
    }
    break;
  case 'H':
    switch (resn[1]) {
    case 'I':
      switch (resn[2]) {
      case 'S': /* HIS */
      case 'D': /* HID */
      case 'E': /* HIE */
      case 'P': /* HIP */
        return true;
        break;
      }
      break;
    }
    break;
  case 'I':
    switch (resn[1]) {
    case 'L':
      switch (resn[2]) {
      case 'E': /* ILE */
        return true;
        break; 
      }
    }
    break;
  case 'L':
    switch (resn[1]) {
    case 'E':
      switch (resn[2]) {
      case 'U': /* LEU */
        return true;
        break;
      }
      break;
    case 'Y':
      switch (resn[2]) {
      case 'S':
        return true;
        break; /* LYS */
      }
      break;
    }
    break;
  case 'M':
    switch (resn[1]) {
    case 'E':
      switch (resn[2]) {
      case 'T': /* MET */
        return true;
        break;
      }
    case 'S':
      switch (resn[2]) {
      case 'E': /* MSE */
        return true;
        break;
      }
    }
    break;
  case 'P':
    switch (resn[1]) {
    case 'H':
      switch (resn[2]) {
      case 'E': /* PHE */
        return true;
        break;
      }
      break;
    case 'R':
      switch (resn[2]) {
      case 'O': /* PRO */
        return true;
        break;
      }
      break;
    case 'T':
      switch (resn[2]) {
      case 'R': /* PTR */
        return true;
        break;
      }
      break;
    }
    break;
  case 'S':
    switch (resn[1]) {
    case 'E':
      switch (resn[2]) {
      case 'R': /* SER */
        return true;
        break;
      }
      break;
    }
    break;
  case 'T':
    switch (resn[1]) {
    case 'H':
      switch (resn[2]) {
      case 'R': /* THR */
        return true;
        break;
      }
      break;
    case 'R':
      switch (resn[2]) {
      case 'P': /* TRP */
        return true;
        break;
      }
      break;
    case 'Y':
      switch (resn[2]) {
      case 'R': /* TYR */
        return true;
        break;
      }
      break;
    }
    break;
  case 'V':
    switch (resn[1]) {
    case 'A':
      switch (resn[2]) {
      case 'L': /* VAL */
        return true;
        break;
      }
      break;
    }
    break;
  }
  return false;
}


/*========================================================================*/

static char getInscodeUpper(const AtomInfoType * ai) {
  char c = ai->inscode;
  if ('a' <= c && c <= 'z')
    return c - ('a' - 'A'); // 97 - 65
  return c;
}

// return false if buffer too small
bool AtomResiFromResv(char *resi, size_t size, int resv, char inscode) {
  if (inscode > ' ')
    return snprintf(resi, size, "%d%c", resv, inscode) < size;
  return snprintf(resi, size, "%d", resv) < size;
}

/*========================================================================*/
PyObject *AtomInfoAsPyList(PyMOLGlobals * G, const AtomInfoType * I)
{
  PyObject *result = NULL;

  result = PyList_New(48);
  
  int version = SettingGetGlobal_f(G, cSetting_pse_export_version) * 1000;
  char resi[8];

  // at some point, change this to !version || ...
  if (version >= 1810) {
    resi[0] = I->inscode;
    resi[1] = '\0';
  } else {
    AtomResiFromResv(resi, sizeof(resi), I);
  }

  PyList_SetItem(result, 0, PyInt_FromLong(I->resv));
  PyList_SetItem(result, 1, PyString_FromString(LexStr(G, I->chain)));
  PyList_SetItem(result, 2, PyString_FromString(I->alt));
  PyList_SetItem(result, 3, PyString_FromString(resi));
  PyList_SetItem(result, 4, PyString_FromString(LexStr(G, I->segi)));
  PyList_SetItem(result, 5, PyString_FromString(LexStr(G, I->resn)));
  PyList_SetItem(result, 6, PyString_FromString(LexStr(G, I->name)));
  PyList_SetItem(result, 7, PyString_FromString(I->elem));
  PyList_SetItem(result, 8, PyString_FromString(LexStr(G, I->textType)));
  PyList_SetItem(result, 9, PyString_FromString(LexStr(G, I->label)));
  PyList_SetItem(result, 10, PyString_FromString(I->ssType));
  PyList_SetItem(result, 11, PyInt_FromLong((int) I->isHydrogen())); // TODO redundant
  PyList_SetItem(result, 12, PyInt_FromLong(I->customType));
  PyList_SetItem(result, 13, PyInt_FromLong(I->priority));
  PyList_SetItem(result, 14, PyFloat_FromDouble(I->b));
  PyList_SetItem(result, 15, PyFloat_FromDouble(I->q));
  PyList_SetItem(result, 16, PyFloat_FromDouble(I->vdw));
  PyList_SetItem(result, 17, PyFloat_FromDouble(I->partialCharge));
  PyList_SetItem(result, 18, PyInt_FromLong(I->formalCharge));
  PyList_SetItem(result, 19, PyInt_FromLong((int) I->hetatm));
  PyList_SetItem(result, 20, PyInt_FromLong((int) I->visRep));
  PyList_SetItem(result, 21, PyInt_FromLong(I->color));
  PyList_SetItem(result, 22, PyInt_FromLong(I->id));
  PyList_SetItem(result, 23, PyInt_FromLong((char) I->cartoon));
  PyList_SetItem(result, 24, PyInt_FromLong(I->flags));
  PyList_SetItem(result, 25, PyInt_FromLong((int) I->bonded));
  PyList_SetItem(result, 26, PyInt_FromLong((int) I->chemFlag));
  PyList_SetItem(result, 27, PyInt_FromLong((int) I->geom));
  PyList_SetItem(result, 28, PyInt_FromLong((int) I->valence));
  PyList_SetItem(result, 29, PyInt_FromLong((int) I->masked));
  PyList_SetItem(result, 30, PyInt_FromLong((int) I->protekted));
  PyList_SetItem(result, 31, PyInt_FromLong((int) I->protons));
  PyList_SetItem(result, 32, PyInt_FromLong(I->unique_id));
  PyList_SetItem(result, 33, PyInt_FromLong((char) I->stereo));
  PyList_SetItem(result, 34, PyInt_FromLong(I->discrete_state));
  PyList_SetItem(result, 35, PyFloat_FromDouble(I->elec_radius));
  PyList_SetItem(result, 36, PyInt_FromLong(I->rank));
  PyList_SetItem(result, 37, PyInt_FromLong((int) I->hb_donor));
  PyList_SetItem(result, 38, PyInt_FromLong((int) I->hb_acceptor));
  PyList_SetItem(result, 39, PyInt_FromLong(0 /* atomic_color */));
  PyList_SetItem(result, 40, PyInt_FromLong((int) I->has_setting));

  const float anisou_stack[] {0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
  const float * anisou = I->anisou ? I->anisou : anisou_stack;
  for (int i = 0; i < 6; ++i) {
    PyList_SetItem(result, 41 + i, PyFloat_FromDouble(anisou[i]));
  }

  PyList_SetItem(result, 47, PyString_FromString(LexStr(G, I->custom)));

  return (PConvAutoNone(result));
}

int AtomInfoFromPyList(PyMOLGlobals * G, AtomInfoType * I, PyObject * list)
{
  int ok = true;
  int tmp_int;
  ov_size ll = 0;
  OrthoLineType temp = "";

#define PCONVPYSTRTOLEXIDX(i, n) { \
  ok = CPythonVal_PConvPyStrToStr_From_List(G, list, i, temp, sizeof(OrthoLineType)); \
  n = LexIdx(G, temp); }

  if(ok)
    ok = PyList_Check(list);
  if(ok)
    ll = PyList_Size(list);
  if(ok)
    ok = CPythonVal_PConvPyIntToInt_From_List(G, list, 0, &I->resv);
  if(ok)
    PCONVPYSTRTOLEXIDX(1, I->chain);
  if(ok)
    ok = CPythonVal_PConvPyStrToStr_From_List(G, list, 2, I->alt, sizeof(Chain));
  if(ok)
  {
    // get inscode from resi
    ok = CPythonVal_PConvPyStrToStr_From_List(G, list, 3, temp, sizeof(temp));
    int i = strlen(temp) - 1;
    if (i >= 0 && !isdigit(temp[i])) {
      I->setInscode(temp[i]);
    }
  }

  if(ok) PCONVPYSTRTOLEXIDX(4, I->segi);
  if(ok) PCONVPYSTRTOLEXIDX(5, I->resn);
  if(ok) PCONVPYSTRTOLEXIDX(6, I->name);
  if(ok)
    ok = CPythonVal_PConvPyStrToStr_From_List(G, list, 7, I->elem, sizeof(ElemName));
  if(ok) PCONVPYSTRTOLEXIDX(8, I->textType);
  if(ok) PCONVPYSTRTOLEXIDX(9, I->label);
  if(ok)
    ok = PConvPyStrToStr(PyList_GetItem(list, 10), I->ssType, sizeof(SSType));
  if(ok)
    ok = PConvPyIntToInt(PyList_GetItem(list, 12), &I->customType);
  if(ok)
    ok = PConvPyIntToInt(PyList_GetItem(list, 13), &I->priority);
  if(ok)
    ok = PConvPyFloatToFloat(PyList_GetItem(list, 14), &I->b);
  if(ok)
    ok = PConvPyFloatToFloat(PyList_GetItem(list, 15), &I->q);
  if(ok)
    ok = PConvPyFloatToFloat(PyList_GetItem(list, 16), &I->vdw);
  if(ok)
    ok = PConvPyFloatToFloat(PyList_GetItem(list, 17), &I->partialCharge);
  if(ok)
    if((ok = CPythonVal_PConvPyIntToInt_From_List(G, list, 18, &tmp_int)))
      I->formalCharge = tmp_int;
  if(ok)
    if((ok = CPythonVal_PConvPyIntToInt_From_List(G, list, 19, &tmp_int)))
      I->hetatm = tmp_int;
  if(ok){
    PyObject *val = PyList_GetItem(list, 20);
    if (PyList_Check(val)){
      ok = PConvPyListToBitmask(val, &I->visRep, cRepCnt);
    } else {
      ok = PConvPyIntToInt(val, &I->visRep);
    }
  }
  if(ok)
    ok = PConvPyIntToInt(PyList_GetItem(list, 21), &I->color);
  if(ok)
    I->color = ColorConvertOldSessionIndex(G, I->color);
  if(ok)
    ok = PConvPyIntToInt(PyList_GetItem(list, 22), &I->id);
  if(ok)
    if((ok = CPythonVal_PConvPyIntToInt_From_List(G, list, 23, &tmp_int)))
      I->cartoon = tmp_int;
  if(ok)
    ok = PConvPyIntToInt(PyList_GetItem(list, 24), (int *) &I->flags);
  if(ok)
    if((ok = CPythonVal_PConvPyIntToInt_From_List(G, list, 25, &tmp_int)))
      I->bonded = tmp_int;
  if(ok)
    if((ok = CPythonVal_PConvPyIntToInt_From_List(G, list, 26, &tmp_int)))
      I->chemFlag = tmp_int;
  if(ok)
    if((ok = CPythonVal_PConvPyIntToInt_From_List(G, list, 27, &tmp_int)))
      I->geom = tmp_int;
  if(ok)
    if((ok = CPythonVal_PConvPyIntToInt_From_List(G, list, 28, &tmp_int)))
      I->valence = tmp_int;
  if(ok)
    if((ok = CPythonVal_PConvPyIntToInt_From_List(G, list, 29, &tmp_int)))
      I->masked = tmp_int;
  if(ok)
    if((ok = CPythonVal_PConvPyIntToInt_From_List(G, list, 30, &tmp_int)))
      I->protekted = tmp_int;
  if(ok)
    ok = PConvPyIntToChar(PyList_GetItem(list, 31), (char *) &I->protons);
  if(ok)
    ok = PConvPyIntToInt(PyList_GetItem(list, 32), &I->unique_id);
  if(ok && I->unique_id) {      /* reserve existing IDs */
    I->unique_id = SettingUniqueConvertOldSessionID(G, I->unique_id);
  }
  if(ok)
    if((ok = CPythonVal_PConvPyIntToInt_From_List(G, list, 33, &tmp_int)))
      I->stereo = tmp_int;
  if(ok && (ll > 34))
    ok = PConvPyIntToInt(PyList_GetItem(list, 34), &I->discrete_state);
  if(ok && (ll > 35))
    ok = PConvPyFloatToFloat(PyList_GetItem(list, 35), &I->elec_radius);
  if(ok && (ll > 36))
    ok = PConvPyIntToInt(PyList_GetItem(list, 36), &I->rank);
  if(ok && (ll > 37))
    if ((ok = CPythonVal_PConvPyIntToInt_From_List(G, list, 37, &tmp_int)))
      I->hb_donor = tmp_int;
  if(ok && (ll > 38))
    if ((ok = CPythonVal_PConvPyIntToInt_From_List(G, list, 38, &tmp_int)))
      I->hb_acceptor = tmp_int;
  if(ok && (ll > 40))
    if ((ok = CPythonVal_PConvPyIntToInt_From_List(G, list, 40, &tmp_int)))
      I->has_setting = tmp_int;
  if(ok && (ll > 46)) {
    // only allocate if not all zero
    float u[6];
    for (int i = 0; ok && i < 6; ++i)
      ok = CPythonVal_PConvPyFloatToFloat_From_List(G, list, 41 + i, u + i);
    if(ok && (u[0] || u[1] || u[2] || u[3] || u[4] || u[5]))
      memcpy(I->get_anisou(), u, 6 * sizeof(float));
  }
  if(ok && (ll > 47)) {
    PCONVPYSTRTOLEXIDX(47, I->custom);
  }
  return (ok);
}

void AtomInfoCopy(PyMOLGlobals * G, const AtomInfoType * src, AtomInfoType * dst, int copy_properties)
{
  /* copy, handling resource management issues... */

  *dst = *src;
  dst->selEntry = 0;
  if(src->unique_id && src->has_setting) {
    dst->unique_id = AtomInfoGetNewUniqueID(G);
    if(!SettingUniqueCopyAll(G, src->unique_id, dst->unique_id))
      dst->has_setting = 0;
  } else {
    dst->unique_id = 0;
    dst->has_setting = 0;
  }
  LexInc(G, dst->label);
  LexInc(G, dst->textType);
  LexInc(G, dst->custom);
  LexInc(G, dst->chain);
  LexInc(G, dst->segi);
  LexInc(G, dst->resn);
  LexInc(G, dst->name);
#ifdef _PYMOL_IP_EXTRAS
#endif
  if (src->anisou) {
    dst->anisou = NULL;
    memcpy(dst->get_anisou(), src->anisou, 6 * sizeof(float));
  }
}

void AtomInfoBondCopy(PyMOLGlobals * G, const BondType * src, BondType * dst)
{
  *(dst) = *(src);

  if(src->unique_id && src->has_setting) {
    dst->unique_id = AtomInfoGetNewUniqueID(G);
    if(!SettingUniqueCopyAll(G, src->unique_id, dst->unique_id))
      dst->has_setting = 0;
  } else {
    dst->unique_id = 0;
    dst->has_setting = 0;
  }
}

void AtomInfoPurgeBond(PyMOLGlobals * G, BondType * bi)
{
  CAtomInfo *I = G->AtomInfo;
  if(bi->has_setting && bi->unique_id) {
    SettingUniqueDetachChain(G, bi->unique_id);
  }
  if(bi->unique_id && I->ActiveIDs) {
    OVOneToAny_DelKey(I->ActiveIDs, bi->unique_id);
    bi->unique_id = 0;
  }
}

void AtomInfoPurge(PyMOLGlobals * G, AtomInfoType * ai)
{
  CAtomInfo *I = G->AtomInfo;
  LexDec(G, ai->textType);
  LexDec(G, ai->custom);
  LexDec(G, ai->label);
  LexDec(G, ai->chain);
  ai->textType = 0;
  ai->custom = 0;
  ai->label = 0;
  ai->chain = 0;
  if(ai->has_setting && ai->unique_id) {
    SettingUniqueDetachChain(G, ai->unique_id);
  }
  if(ai->unique_id) {
    ExecutiveUniqueIDAtomDictInvalidate(G);

    if (I->ActiveIDs)
      OVOneToAny_DelKey(I->ActiveIDs, ai->unique_id);
  }
#ifdef _PYMOL_IP_EXTRAS
#endif
  DeleteAP(ai->anisou);
}


/*========================================================================*/
/**
 * Transfer `mask` selected atomic properties (such as b, q, text_type, etc.)
 * from `src` to `dst`, while keeping all atomic identifiers untouched.
 * Purges `src`.
 *
 * @param mask bitmask of `cAIC_*` bits
 */
void AtomInfoCombine(PyMOLGlobals * G, AtomInfoType * dst, AtomInfoType&& src_, int mask)
{
  AtomInfoType* src = &src_;
  if(mask & cAIC_tt) {
    std::swap(dst->textType, src->textType);
  }
  if(mask & cAIC_ct)
    dst->customType = src->customType;
  if(mask & cAIC_pc)
    dst->partialCharge = src->partialCharge;
  if(mask & cAIC_fc)
    dst->formalCharge = src->formalCharge;
  if(mask & cAIC_flags)
    dst->flags = src->flags;
  if(mask & cAIC_b)
    dst->b = src->b;
  if(mask & cAIC_q)
    dst->q = src->q;
  if(mask & cAIC_id)
    dst->id = src->id;
  if(mask & cAIC_state)
    dst->discrete_state = src->discrete_state;
  if(mask & cAIC_rank)
    dst->rank = src->rank;
  dst->temp1 = src->temp1;

  SWAP_NOREF(dst->has_setting, src->has_setting);
  std::swap(dst->unique_id, src->unique_id);
#ifdef _PYMOL_IP_EXTRAS
  std::swap(dst->prop_id, src->prop_id);
#endif

  /* keep all existing names, identifiers, etc. */
  /* also keep all existing selections,
     colors, masks, and visible representations */
    /* leaves dst->label untouched */

  AtomInfoPurge(G, src);
}


/*========================================================================*/
/**
 * Make atom names in `atInfo1` unique w.r.t.\ `atInfo0` (and to `atInfo1` itself).
 * @param atInfo0 List of reference atoms
 * @param n0 Size of atInfo0 list
 * @param atInfo1 List of atoms which need to be made unique
 * @param flag1 Optional whitelist mask for `atInfo1` or NULL
 * @param n1 Size of atInfo1 list
 * @param mol Optional reference molecule to limit to atoms with coordinates
 * @return Number of renamed atoms
 */
int AtomInfoUniquefyNames(PyMOLGlobals * G, const AtomInfoType * atInfo0, int n0,
                          AtomInfoType * atInfo1, int *flag1, int n1,
                          const ObjectMolecule * mol)
{
  /* makes sure all names in atInfo1 are unique WRT 0 and 1 */
  auto ignore_case = SettingGet<bool>(G, cSetting_ignore_case);

  /* tricky optimizations to avoid n^2 dependence in this operation */
  int result = 0;
  int a, b, c;
  const AtomInfoType *ai0, *lai1, *lai0;
  AtomInfoType *ai1;
  int st1, nd1, st0, nd0;       /* starts and ends */
  int matchFlag;
  int bracketFlag;
  WordType name;

  ai1 = atInfo1;
  lai0 = NULL;                  /* last atom compared against in each object */
  lai1 = NULL;
  st0 = 0;
  nd0 = 0;
  st1 = 0;
  nd1 = 0;
  c = 1;
  /* ai1->name is the atom we're currently on */

  b = 0;
  while(b < n1) {
    matchFlag = false;

    if(!ai1->name)
      matchFlag = true;

    if(!matchFlag) {
      /* check within object 1 */

      if(!lai1)
        bracketFlag = true;
      else if(!AtomInfoSameResidue(G, lai1, ai1))
        bracketFlag = true;
      else
        bracketFlag = false;
      if(bracketFlag) {
        c = 1;
        AtomInfoBracketResidue(G, atInfo1, n1, ai1, &st1, &nd1);
        lai1 = ai1;
      }

      ai0 = atInfo1 + st1;
      for(a = st1; a <= nd1; a++) {
        if(!WordMatchExact(G, ai1->name, ai0->name, ignore_case))
          ai0++;
        else if(!AtomInfoSameResidue(G, ai1, ai0))
          ai0++;
        else if(ai1 != ai0) {
          matchFlag = true;
          break;
        } else
          ai0++;
      }
    }

    if(!matchFlag) {
      if(atInfo0) {
        /* check within object 2 */

        if (!lai0 || !AtomInfoSameResidue(G, lai0, ai1)) {
          AtomInfoBracketResidue(G, atInfo0, n0, ai1, &st0, &nd0);
          lai0 = ai1;
        }

        for (a = st0; a <= nd0; ++a) {
          ai0 = atInfo0 + a;

          if (WordMatchExact(G, ai1->name, ai0->name, ignore_case) &&
              AtomInfoSameResidue(G, ai1, ai0) && ai1 != ai0 &&
              (!mol || mol->atomHasAnyCoordinates(a))) {
            matchFlag = true;
            break;
          }
        }
      }
    }

    if(matchFlag && ((!flag1) || flag1[ai1 - atInfo1])) {
      if(c < 100) {
        if((c < 10) && ai1->elem[1])    /* try to keep halogens 3 or under */
          sprintf(name, "%2s%1d", ai1->elem, c);
        else
          sprintf(name, "%1s%02d", ai1->elem, c);
      } else {
        sprintf(name, "%1d%1s%02d", c / 100, ai1->elem, c % 100);
      }
      LexAssign(G, ai1->name, name);
      result++;
      c = c + 1;
    } else {
      ai1++;
      b++;
    }
  }
  return result;
}

/**
 * Make atom names in `atoms` unique w.r.t.\ atoms-with-coordinates in `mol`.
 * @param mol Reference molecule
 * @param atoms List of atoms which need to be made unique
 * @param natoms Size of atoms list
 * @return Number of renamed atoms
 */
int AtomInfoUniquefyNames(
    const ObjectMolecule* mol, AtomInfoType* atoms, size_t natoms)
{
  return AtomInfoUniquefyNames(
      mol->G, mol->AtomInfo, mol->NAtom, atoms, nullptr, natoms, mol);
}

/*========================================================================*/
void AtomInfoBracketResidue(PyMOLGlobals * G, const AtomInfoType * ai0, int n0,
                            const AtomInfoType * ai, int *st, int *nd)
{
  /* inefficient but reliable way to find where residue atoms are located in an object 
   * for purpose of residue-based operations */
  int a;
  const AtomInfoType *ai1;

  *st = 0;
  *nd = n0 - 1;
  ai1 = ai0;
  for(a = 0; a < n0; a++) {
    if(!AtomInfoSameResidue(G, ai, ai1++))
      *st = a;
    else
      break;
  }
  ai1 = ai0 + n0 - 1;
  for(a = n0 - 1; a >= 0; a--) {
    if(!AtomInfoSameResidue(G, ai, ai1--)) {
      *nd = a;
    } else
      break;
  }
}


/*========================================================================*/
void AtomInfoBracketResidueFast(PyMOLGlobals * G, const AtomInfoType * ai0, int n0, int cur,
                                int *st, int *nd)
{
  /* efficient but unreliable way to find where residue atoms are located in an object 
   * for purpose of residue-based operations.
   * This function finds all atoms in the AtomInfoType array in both directions that 
   * belong to the same residue. It starts the search from the "cur" offset, and returns
   * the beginning offset (st) and the ending offset (nd).*/
  int a;
  const AtomInfoType *ai1;

  *st = cur;
  *nd = cur;
  ai0 = ai0 + cur;
  ai1 = ai0 - 1;
  for(a = cur - 1; a >= 0; a--) {
    if(!AtomInfoSameResidue(G, ai0, ai1--))
      break;
    *st = a;
  }
  ai1 = ai0 + 1;
  for(a = cur + 1; a < n0; a++) {
    if(!AtomInfoSameResidue(G, ai0, ai1++))
      break;
    *nd = a;
  }
}


/*========================================================================*/
int AtomInfoUpdateAutoColor(PyMOLGlobals * G)
{
  CAtomInfo *I = G->AtomInfo;
  if(SettingGetGlobal_b(G, cSetting_auto_color))
    I->CColor = ColorGetNext(G);
  else
    I->CColor = ColorGetIndex(G, "carbon");
  return I->CColor;
}


/*========================================================================*/
void AtomInfoPrimeColors(PyMOLGlobals * G)
{
  CAtomInfo *I = G->AtomInfo;

  I->NColor = ColorGetIndex(G, "nitrogen");
  I->CColor = ColorGetIndex(G, "carbon");

  I->HColor = ColorGetIndex(G, "hydrogen");
  I->OColor = ColorGetIndex(G, "oxygen");
  I->SColor = ColorGetIndex(G, "sulfur");

  I->ClColor = ColorGetIndex(G, "chlorine");
  I->BrColor = ColorGetIndex(G, "bromine");
  I->FColor = ColorGetIndex(G, "fluorine");
  I->IColor = ColorGetIndex(G, "iodine");

  I->PColor = ColorGetIndex(G, "phosphorus");

  I->MgColor = ColorGetIndex(G, "magnesium");
  I->MnColor = ColorGetIndex(G, "manganese");

  I->NaColor = ColorGetIndex(G, "sodium");
  I->KColor = ColorGetIndex(G, "potassium");
  I->CaColor = ColorGetIndex(G, "calcium");

  I->CuColor = ColorGetIndex(G, "copper");
  I->FeColor = ColorGetIndex(G, "iron");
  I->ZnColor = ColorGetIndex(G, "zinc");

  I->SeColor = ColorGetIndex(G, "selenium");
  I->DColor = ColorGetIndex(G, "deuterium");
}


/*========================================================================*/
float AtomInfoGetBondLength(PyMOLGlobals * G, const AtomInfoType * ai1, const AtomInfoType * ai2)
{
  float result = 1.6F;
  const AtomInfoType *a1, *a2;

  /* very simple for now ... 
     flush this out with pattern-based parameters later on */

  if(ai1->protons > ai2->protons) {
    a1 = ai2;
    a2 = ai1;
  } else {
    a1 = ai1;
    a2 = ai2;
  }
  switch (a1->protons) {


/* hydrogen =========================== */
  case cAN_H:
    switch (a2->protons) {
    case cAN_H:
      result = 0.74F;
      break;
    case cAN_C:
      result = 1.09F;
      break;
    case cAN_N:
      result = 1.01F;
      break;
    case cAN_S:
      result = 1.34F;
      break;
    case cAN_O:
      result = 0.96F;
      break;
    default:
      result = 1.09F;
      break;
    }
    break;

/* carbon =========================== */

  case cAN_C:                  /* carbon */
    switch (a1->geom) {

    case cAtomInfoLinear:      /* linear carbon ============= */
      switch (a2->geom) {
      case cAtomInfoLinear:
        switch (a2->protons) {
        case cAN_C:
          result = 1.20F;
          break;                /* C#^C */
        case cAN_N:
          result = 1.16F;
          break;                /* C#^N */
        default:
          result = 1.20F;
          break;
        }
        break;
      case cAtomInfoPlanar:
        switch (a2->protons) {
        case cAN_I:
          result = 2.14F;
          break;                /* normal single bond lengths */
        case cAN_Cl:
          result = 1.77F;
          break;
        case cAN_Br:
          result = 1.94F;
          break;
        case cAN_F:
          result = 1.35F;
          break;
        case cAN_S:
          result = 1.82F;
          break;
        case cAN_N:
          result = 1.47F;
          break;
        case cAN_O:
          result = 1.43F;
          break;
        case cAN_P:
          result = 1.84F;
          break;
        case cAN_C:
          result = 1.54F;
          break;
        default:
          result = 1.54F;
          break;
        }
        break;
      default:                 /* ?#C-^? */
        switch (a2->protons) {
        case cAN_I:
          result = 2.14F;
          break;                /* normal single bond lengths */
        case cAN_Cl:
          result = 1.77F;
          break;
        case cAN_Br:
          result = 1.94F;
          break;
        case cAN_F:
          result = 1.35F;
          break;
        case cAN_S:
          result = 1.82F;
          break;
        case cAN_N:
          result = 1.47F;
          break;
        case cAN_O:
          result = 1.43F;
          break;
        case cAN_P:
          result = 1.84F;
          break;
        case cAN_C:
          result = 1.54F;
          break;
        default:
          result = 1.54F;
          break;
        }
        break;
      }
      break;

    case cAtomInfoPlanar:      /* planer carbon ============= */
      switch (a2->geom) {
      case cAtomInfoLinear:
        switch (a2->protons) {
        case cAN_I:
          result = 2.14F;
          break;                /* normal single bond lengths */
        case cAN_Cl:
          result = 1.77F;
          break;
        case cAN_Br:
          result = 1.94F;
          break;
        case cAN_F:
          result = 1.35F;
          break;
        case cAN_S:
          result = 1.82F;
          break;
        case cAN_N:
          result = 1.47F;
          break;
        case cAN_O:
          result = 1.43F;
          break;
        case cAN_P:
          result = 1.84F;
          break;
        case cAN_C:
          result = 1.54F;
          break;
        default:
          result = 1.54F;
          break;
        }
        break;
      case cAtomInfoPlanar:
        switch (a2->protons) {
        case cAN_C:
          result = 1.34F;
          break;                /* C=^C or ?=C-^C=? */
        case cAN_O:
          result = 1.20F;
          break;                /* carbonyl */
        case cAN_N:
          result = 1.29F;
          break;                /* C=^N or ?=C-^N=? */
        case cAN_S:
          result = 1.60F;
          break;                /* C=^S or ?=C-^S-?=? */
        default:
          result = 1.34F;
          break;
        }
        break;
      default:                 /* ?#C-^? */
        switch (a2->protons) {
        case cAN_I:
          result = 2.14F;
          break;                /* normal single bond lengths */
        case cAN_Cl:
          result = 1.77F;
          break;
        case cAN_Br:
          result = 1.94F;
          break;
        case cAN_F:
          result = 1.35F;
          break;
        case cAN_S:
          result = 1.71F;
          break;                /* mmod */
        case cAN_N:
          result = 1.47F;
          break;
        case cAN_O:
          result = 1.43F;
          break;
        case cAN_P:
          result = 1.84F;
          break;
        case cAN_C:
          result = 1.54F;
          break;
        default:
          result = 1.54F;
          break;
        }
        break;
      }
      break;

    default:                   /* tetrahedral carbon */
      switch (a2->protons) {
      case cAN_I:
        result = 2.14F;
        break;                  /* normal single bond lengths */
      case cAN_Cl:
        result = 1.77F;
        break;
      case cAN_Br:
        result = 1.94F;
        break;
      case cAN_F:
        result = 1.35F;
        break;
      case cAN_S:
        result = 1.82F;
        break;
      case cAN_N:
        result = 1.47F;
        break;
      case cAN_O:
        result = 1.43F;
        break;
      case cAN_P:
        result = 1.84F;
        break;
      case cAN_C:
        result = 1.54F;
        break;
      default:
        result = 1.54F;
        break;
      }
      break;
    }
    break;


/* nitrogen ========================= */

  case cAN_N:
    if((a1->geom == cAtomInfoPlanar) && (a2->geom == cAtomInfoPlanar))
      switch (a2->protons) {
      case cAN_N:
        result = 1.25F;
        break;
      case cAN_O:
        result = 1.21F;
        break;
      case cAN_S:
        result = 1.53F;
        break;                  /* interpolated */
      default:
        result = 1.25F;
        break;
    } else {
      switch (a2->protons) {
      case cAN_N:
        result = 1.45F;
        break;
      case cAN_O:
        result = 1.40F;
        break;
      case cAN_S:
        result = 1.75F;
        break;                  /* interpolated */
      default:
        result = 1.45F;
        break;
      }
    }
    break;


/* oxygen =========================== */

  case cAN_O:
    if((a1->geom == cAtomInfoPlanar) && (a2->geom == cAtomInfoPlanar))
      switch (a2->protons) {
      case cAN_O:
        result = 1.35F;
        break;                  /* guess */
      case cAN_S:
        result = 1.44F;
        break;                  /* macromodel */
      default:
        result = 1.35F;
        break;
    } else if(a1->geom == cAtomInfoPlanar) {
      switch (a2->protons) {
      case cAN_O:
        result = 1.35F;
        break;                  /* guess */
      case cAN_S:
        result = 1.44F;
        break;                  /* macromodel */
      default:
        result = 1.35F;
        break;
      }
    } else {
      switch (a2->protons) {
      case cAN_O:
        result = 1.40F;
        break;
      case cAN_S:
        result = 1.75F;
        break;                  /* interpolated */
      default:
        result = 1.45F;
        break;
      }
    }
    break;


/* sulfur =========================== */

  case cAN_S:
    switch (a2->protons) {
    case cAN_S:
      result = 2.05F;
      break;                    /* interpolated */
    default:
      result = 1.82F;
      break;
    }
    break;

    /* fall-back to old method */
  default:

    result = 0.0;
    switch (a1->geom) {
    case cAtomInfoLinear:
      result += 1.20F;
      break;
    case cAtomInfoPlanar:
      result += 1.34F;
      break;
    default:
      result += 1.54F;
      break;
    }
    switch (a2->geom) {
    case cAtomInfoLinear:
      result += 1.20F;
      break;
    case cAtomInfoPlanar:
      result += 1.34F;
      break;
    default:
      result += 1.54F;
      break;
    }
    result /= 2.0F;
    break;
  }
  return (result);
}

/**
 * Periodic table of elements
 *
 * Note: Default VDW radius is 1.80
 */
const ElementTableItemType ElementTable[] = {
  {"lonepair",          "LP",   0.50,   0.000000},
  {"hydrogen",          "H",    1.20,   1.007940},
  {"helium",            "He",   1.40,   4.002602},
  {"lithium",           "Li",   1.82,   6.941000},
  {"beryllium",         "Be",   1.80,   9.012182},
  {"boron",             "B",    1.85,  10.811000},
  {"carbon",            "C",    1.70,  12.010700},
  {"nitrogen",          "N",    1.55,  14.006700},
  {"oxygen",            "O",    1.52,  15.999400},
  {"fluorine",          "F",    1.47,  18.998403},
  {"neon",              "Ne",   1.54,  20.179700},
  {"sodium",            "Na",   2.27,  22.989770},
  {"magnesium",         "Mg",   1.73,  24.305000},
  {"aluminum",          "Al",   2.00,  26.981538},
  {"silicon",           "Si",   2.10,  28.085500},
  {"phosphorus",        "P",    1.80,  30.973761},
  {"sulfur",            "S",    1.80,  32.065000},
  {"chlorine",          "Cl",   1.75,  35.453000},
  {"argon",             "Ar",   1.88,  39.948000},
  {"potassium",         "K",    2.75,  39.098300},
  {"calcium",           "Ca",   1.80,  40.078000},
  {"scandium",          "Sc",   1.80,  44.955910},
  {"titanium",          "Ti",   1.80,  47.867000},
  {"vanadium",          "V",    1.80,  50.941500},
  {"chromium",          "Cr",   1.80,  51.996100},
  {"manganese",         "Mn",   1.73,  54.938049},
  {"iron",              "Fe",   1.80,  55.845000},
  {"cobalt",            "Co",   1.80,  58.933200},
  {"nickel",            "Ni",   1.63,  58.693400},
  {"copper",            "Cu",   1.40,  63.546000},
  {"zinc",              "Zn",   1.39,  65.390000},
  {"gallium",           "Ga",   1.87,  69.723000},
  {"germanium",         "Ge",   1.80,  72.640000},
  {"arsenic",           "As",   1.85,  74.921600},
  {"selenium",          "Se",   1.90,  78.960000},
  {"bromine",           "Br",   1.85,  79.904000},
  {"krypton",           "Kr",   2.02,  83.800000},
  {"rubidium",          "Rb",   1.80,  85.467800},
  {"strontium",         "Sr",   1.80,  87.620000},
  {"yttrium",           "Y",    1.80,  88.905850},
  {"zirconium",         "Zr",   1.80,  91.224000},
  {"niobium",           "Nb",   1.80,  92.906380},
  {"molybdenum",        "Mo",   1.80,  95.940000},
  {"technetium",        "Tc",   1.80,  98.000000},
  {"ruthenium",         "Ru",   1.80, 101.070000},
  {"rhodium",           "Rh",   1.80, 102.905500},
  {"palladium",         "Pd",   1.63, 106.420000},
  {"silver",            "Ag",   1.72, 107.868200},
  {"cadmium",           "Cd",   1.58, 112.411000},
  {"indium",            "In",   1.93, 114.818000},
  {"tin",               "Sn",   2.17, 118.710000},
  {"antimony",          "Sb",   1.80, 121.760000},
  {"tellurium",         "Te",   2.06, 127.600000},
  {"iodine",            "I",    1.98, 126.904470},
  {"xenon",             "Xe",   2.16, 131.293000},
  {"cesium",            "Cs",   1.80, 132.905450},
  {"barium",            "Ba",   1.80, 137.327000},
  {"lanthanum",         "La",   1.80, 138.905500},
  {"cerium",            "Ce",   1.80, 140.116000},
  {"praseodymium",      "Pr",   1.80, 140.907650},
  {"neodymium",         "Nd",   1.80, 144.240000},
  {"promethium",        "Pm",   1.80, 145.000000},
  {"samarium",          "Sm",   1.80, 150.360000},
  {"europium",          "Eu",   1.80, 151.964000},
  {"gadolinium",        "Gd",   1.80, 157.250000},
  {"terbium",           "Tb",   1.80, 158.925340},
  {"dysprosium",        "Dy",   1.80, 162.500000},
  {"holmium",           "Ho",   1.80, 164.930320},
  {"erbium",            "Er",   1.80, 167.259000},
  {"thulium",           "Tm",   1.80, 168.934210},
  {"ytterbium",         "Yb",   1.80, 173.040000},
  {"lutetium",          "Lu",   1.80, 174.967000},
  {"hafnium",           "Hf",   1.80, 178.490000},
  {"tantalum",          "Ta",   1.80, 180.947900},
  {"tungsten",          "W",    1.80, 183.840000},
  {"rhenium",           "Re",   1.80, 186.207000},
  {"osmium",            "Os",   1.80, 190.230000},
  {"iridium",           "Ir",   1.80, 192.217000},
  {"platinum",          "Pt",   1.75, 195.078000},
  {"gold",              "Au",   1.66, 196.966550},
  {"mercury",           "Hg",   1.55, 200.590000},
  {"thallium",          "Tl",   1.96, 204.383300},
  {"lead",              "Pb",   2.02, 207.200000},
  {"bismuth",           "Bi",   1.80, 208.980380},
  {"polonium",          "Po",   1.80, 208.980000},
  {"astatine",          "At",   1.80, 209.990000},
  {"radon",             "Rn",   1.80, 222.020000},
  {"francium",          "Fr",   1.80, 223.020000},
  {"radium",            "Ra",   1.80, 226.030000},
  {"actinium",          "Ac",   1.80, 227.030000},
  {"thorium",           "Th",   1.80, 232.038100},
  {"protactinium",      "Pa",   1.80, 231.035880},
  {"uranium",           "U",    1.86, 238.028910},
  {"neptunium",         "Np",   1.80, 237.050000},
  {"plutonium",         "Pu",   1.80, 244.060000},
  {"americium",         "Am",   1.80, 243.060000},
  {"curium",            "Cm",   1.80, 247.070000},
  {"berkelium",         "Bk",   1.80, 247.070000},
  {"californium",       "Cf",   1.80, 251.080000},
  {"einsteinium",       "Es",   1.80, 252.080000},
  {"fermium",           "Fm",   1.80, 257.100000},
  {"mendelevium",       "Md",   1.80, 258.100000},
  {"nobelium",          "No",   1.80, 259.100000},
  {"lawrencium",        "Lr",   1.80, 262.110000},
  {"rutherfordium",     "Rf",   1.80, 261.110000},
  {"dubnium",           "Db",   1.80, 262.110000},
  {"seaborgium",        "Sg",   1.80, 266.120000},
  {"bohrium",           "Bh",   1.80, 264.120000},
  {"hassium",           "Hs",   1.80, 269.130000},
  {"meitnerium",        "Mt",   1.80, 268.140000},
  {"darmstadtium",      "Ds",   1.80, 281.000000},
  {"roentgenium",       "Rg",   1.80, 281.000000},
  {"copernicium",       "Cn",   1.80, 285.000000},
  {"nihonium",          "Nh",   1.80, 286.000000},
  {"flerovium",         "Fl",   1.80, 289.000000},
  {"moscovium",         "Mc",   1.80, 290.000000},
  {"livermorium",       "Lv",   1.80, 293.000000},
  {"tennessine",        "Ts",   1.80, 294.000000},
  {"oganesson",         "Og",   1.80, 294.000000},
  {NULL,                NULL,   0.00,   0.000000}
};

const int ElementTableSize = sizeof(ElementTable) / sizeof(ElementTable[0]) - 1;

/**
 * Assign atomic color, or G->AtomInfo->CColor in case of carbon.
 */
void AtomInfoAssignColors(PyMOLGlobals * G, AtomInfoType * at1)
{
  at1->color = AtomInfoGetColor(G, at1);
}

/**
 * Get atomic color, based on protons and elem
 * @return color index
 */
int AtomInfoGetColor(PyMOLGlobals * G, const AtomInfoType * at1)
{
  // fast lookup for most common elements
  switch (at1->protons) {
    case cAN_H:
      if (at1->elem[0] == 'D')
        return G->AtomInfo->DColor;
      return G->AtomInfo->HColor;
    case cAN_N: return G->AtomInfo->NColor;
    case cAN_C: return G->AtomInfo->CColor;
    case cAN_O: return G->AtomInfo->OColor;
    case cAN_P: return G->AtomInfo->PColor;
  }

  // general by-name lookup (exclude LP, PS)
  if (at1->protons > 0 && at1->protons < ElementTableSize)
    return ColorGetIndex(G, ElementTable[at1->protons].name);

  // special cases
  if (strcmp(at1->elem, "PS") == 0)
    return ColorGetIndex(G, "pseudoatom");
  if (strcmp(at1->elem, "LP") == 0)
    return ColorGetIndex(G, "lonepair");

  return G->AtomInfo->DefaultColor;
}

void AtomInfoFreeSortedIndexes(PyMOLGlobals * G, int **index, int **outdex)
{
  FreeP(*index);
  FreeP(*outdex);
}

static int AtomInfoNameCompare(PyMOLGlobals * G, const char *name1, const char *name2)
{
  const char *n1, *n2;
  int cmp;

  if((name1[0] >= '0') && (name1[0] <= '9'))
    n1 = name1 + 1;
  else
    n1 = name1;
  if((name2[0] >= '0') && (name2[0] <= '9'))
    n2 = name2 + 1;
  else
    n2 = name2;
  cmp = WordCompare(G, n1, n2, true);

  if(cmp)
    return cmp;
  return WordCompare(G, name1, name2, true);

}

static int AtomInfoNameCompare(PyMOLGlobals * G, const lexidx_t& name1, const lexidx_t& name2)
{
  if (name1 == name2)
    return 0;
  return AtomInfoNameCompare(G, LexStr(G, name1), LexStr(G, name2));
}

/**
 * Compares atoms based on all atom identifiers, discrete state, priority,
 * hetatm (optional) and rank (optional)
 *
 * Insertion code sorting depends on settings:
 * - pdb_insertions_go_first
 * - rank_assisted_sorts
 *
 * Returns:
 *   0: at1 == at2
 *  -1: at1 < at2
 *   1: at1 > at2
 */
static int AtomInfoCompare(PyMOLGlobals *G, const AtomInfoType *at1, const AtomInfoType *at2,
    bool ignore_hetatm, bool ignore_rank)
{
  int wc;

  if ((wc = WordCompare(G, at1->segi, at2->segi, false))) return wc;
  if ((wc = WordCompare(G, at1->chain, at2->chain, false))) return wc;
  if (!ignore_hetatm && at1->hetatm != at2->hetatm) return (at2->hetatm) ? -1 : 1;
  if (at1->resv != at2->resv) return (at1->resv < at2->resv) ? -1 : 1;

  if ((wc = getInscodeUpper(at1) - getInscodeUpper(at2))) {
    if (SettingGetGlobal_b(G, cSetting_pdb_insertions_go_first))
      return (!at1->inscode) ? 1 : (!at2->inscode) ? -1 : wc;

    if ((at1->rank != at2->rank) && SettingGetGlobal_b(G, cSetting_rank_assisted_sorts))
      return (at1->rank < at2->rank) ? -1 : 1;

    return wc;
  }

  if ((wc = WordCompare(G, at1->resn, at2->resn, true))) return wc;
  if (at1->discrete_state != at2->discrete_state) return (at1->discrete_state < at2->discrete_state) ? -1 : 1;

  // if this looks like a "bulk" het group with no residue number, then don't
  // compare by name/alt/priority (affects mmCIF with cif_use_auth=0)
  if (!ignore_rank && !at1->resv && at1->hetatm)
    goto rank_compare;

  if (at1->priority != at2->priority) return (at1->priority < at2->priority) ? -1 : 1;

  // Changed (PyMOL 2.1): name before alt
  if ((wc = AtomInfoNameCompare(G, at1->name, at2->name))) return wc;

  // Changed (PyMOL 2.1): empty alt goes first: '' < 'A' < 'B'
  if (at1->alt[0] != at2->alt[0]) {
    return (at1->alt[0] < at2->alt[0]) ? -1 : 1;
  }

rank_compare:
  if (!ignore_rank && at1->rank != at2->rank) return (at1->rank < at2->rank) ? -1 : 1;

  return 0;
}

int AtomInfoCompare(PyMOLGlobals * G, const AtomInfoType * at1, const AtomInfoType * at2)
{
  return AtomInfoCompare(G, at1, at2, false, false);
}

int AtomInfoCompareIgnoreRankHet(PyMOLGlobals * G, const AtomInfoType * at1, const AtomInfoType * at2)
{
  return AtomInfoCompare(G, at1, at2, true, true);
}

int AtomInfoCompareIgnoreRank(PyMOLGlobals * G, const AtomInfoType * at1, const AtomInfoType * at2)
{
  return AtomInfoCompare(G, at1, at2, false, true);
}

int AtomInfoCompareIgnoreHet(PyMOLGlobals * G, const AtomInfoType * at1, const AtomInfoType * at2)
{
  return AtomInfoCompare(G, at1, at2, true, false);
}

/**
 * Function only used for matching atoms of two aligned residues.
 *
 * Changed (PyMOL 2.1):
 * Consider alt codes to be equal if at least one is empty.
 * Otherwise, alt-code (C-alpha) atoms can't be aligned to non-alt atoms
 * (also affects morphing).
 */
int AtomInfoNameOrder(PyMOLGlobals * G, const AtomInfoType * at1, const AtomInfoType * at2)
{
  int result;
  if(at1->alt[0] == at2->alt[0] || !at1->alt[0] || !at2->alt[0]) {
    if(at1->priority == at2->priority) {
      result = AtomInfoNameCompare(G, at1->name, at2->name);
    } else if(at1->priority < at2->priority) {
      result = -1;
    } else {
      result = 1;
    }
  } else if(at1->alt[0] < at2->alt[0]) {
    result = -1;
  } else {
    result = 1;
  }
  return (result);
}

int AtomInfoSameResidue(PyMOLGlobals * G, const AtomInfoType * at1, const AtomInfoType * at2)
{
  return (
      at1->resv == at2->resv &&
      at1->chain == at2->chain &&
      at1->hetatm == at2->hetatm &&
      at1->discrete_state == at2->discrete_state &&
      at1->inscode == at2->inscode &&
      at1->segi == at2->segi &&
      WordMatchExact(G, at1->resn, at2->resn, SettingGet<bool>(G, cSetting_ignore_case)));
}

int AtomInfoSameResidueP(PyMOLGlobals * G, const AtomInfoType * at1, const AtomInfoType * at2)
{
  if(at1 && at2)
    return AtomInfoSameResidue(G, at1, at2);
  return 0;
}

int AtomInfoSameChainP(PyMOLGlobals * G, const AtomInfoType * at1, const AtomInfoType * at2)
{
  if(at1 && at2)
    if(at1->chain == at2->chain)
      if(at1->segi == at2->segi)
        return 1;
  return 0;
}

int AtomInfoSameSegmentP(PyMOLGlobals * G, const AtomInfoType * at1, const AtomInfoType * at2)
{
  if(at1 && at2)
    if(at1->segi == at2->segi)
      return 1;
  return 0;
}

int AtomInfoSequential(PyMOLGlobals * G, const AtomInfoType * at1, const AtomInfoType * at2, int mode)
{
  if(mode > 0) {
    if(at1->hetatm == at2->hetatm) {
      if(mode > 1) {
        if(at1->segi == at2->segi) {
          if(mode > 2) {
            if(at1->chain == at2->chain) {
              if(mode > 3) {
                if(at1->resv == at2->resv) {
                  if(mode > 4) {
                    if(at1->inscode == at2->inscode)
                      return 1;
                    if(at1->inscode + 1 == at2->inscode)
                      return 1;
                  } else {
                    return 1;   /* no resi check */
                  }
                } else if((at1->resv + 1) == at2->resv)
                  return 1;
              } else {
                return 1;       /* no resv check */
              }
            }
          } else {
            return 1;           /* no chain check */
          }
        }
      } else {
        return 1;               /* no segi check */
      }
    }
  } else {
    return 1;                   /* no hetatm check */
  }
  return 0;
}

/**
 * Used in "rms" and "update" for matching two selections
 * @return 1 if atoms match, 0 otherwise
 */
int AtomInfoMatch(PyMOLGlobals * G, const AtomInfoType * at1, const AtomInfoType * at2,
    bool ignore_case, bool ignore_case_chain)
{
  if(at1->resv == at2->resv)
  if(WordMatchExact(G, at1->chain, at2->chain, ignore_case_chain))
    if(WordMatchExact(G, at1->name, at2->name, ignore_case))
      if(WordMatchExact(G, at1->inscode, at2->inscode, ignore_case))
        if(WordMatchExact(G, at1->resn, at2->resn, ignore_case))
          if(WordMatchExact(G, at1->segi, at2->segi, ignore_case_chain))
            if(WordMatchExact(G, at1->alt[0], at2->alt[0], ignore_case))
              return 1;
  return 0;
}

int AtomInfoIsFreeCation(PyMOLGlobals * G, const AtomInfoType * I)
{
  switch (I->protons) {
  case cAN_Na:
  case cAN_K:
  case cAN_Ca:
  case cAN_Mg:
  case cAN_Mn:
  case cAN_Sr:
    return true;
  }
  return false;
}

int AtomInfoGetExpectedValence(PyMOLGlobals * G, const AtomInfoType * I)
{
  int result = -1;              /* negative indicates minimum expected valence (abs)
                                   but it could be higher  */

  if(I->formalCharge == 0) {
    switch (I->protons) {
    case cAN_H:
      result = 1;
      break;
    case cAN_C:
      result = 4;
      break;
    case cAN_N:
      result = 3;
      break;
    case cAN_O:
      result = 2;
      break;
    case cAN_F:
      result = 1;
      break;
    case cAN_Cl:
      result = 1;
      break;
    case cAN_Br:
      result = 1;
      break;
    case cAN_I:
      result = 1;
      break;
    case cAN_Na:
      result = 1;
      break;
    case cAN_Ca:
      result = 1;
      break;
    case cAN_K:
      result = 1;
      break;
    case cAN_Mg:
      result = 2;
      break;
    case cAN_Zn:
      result = -1;
      break;
    case cAN_S:
      result = -2;
      break;
    case cAN_P:
      result = -3;
      break;
    }
  } else if(I->formalCharge == 1) {
    switch (I->protons) {
    case cAN_N:
      result = 4;
      break;
    case cAN_O:
      result = 3;
      break;
    case cAN_Na:
      result = 0;
      break;
    case cAN_Ca:
      result = 0;
      break;
    case cAN_K:
      result = 0;
      break;
    case cAN_Mg:
      result = 1;
      break;
    case cAN_Zn:
      result = -1;
      break;
    case cAN_S:
      result = -2;
      break;
    case cAN_P:
      result = -3;
      break;
    }
  } else if(I->formalCharge == -1) {
    switch (I->protons) {
    case cAN_N:
      result = 2;
      break;
    case cAN_O:
      result = 1;
      break;
    case cAN_C:
      result = 3;
      break;
    case cAN_Zn:
      result = -1;
      break;
    case cAN_S:
      result = -2;
      break;
    case cAN_P:
      result = -3;
      break;
    }
  } else if(I->formalCharge == 2) {
    switch (I->protons) {
    case cAN_Mg:
      result = 0;
      break;
    case cAN_Zn:
      result = -1;
      break;
    case cAN_S:
      result = -2;
      break;
    case cAN_P:
      result = -3;
      break;
    }
  }
  return (result);
}

/**
 * Get number of protons for element symbol
 */
static int get_protons(const char * symbol)
{
  char titleized[4];
  static std::map<pymol::zstring_view, int> lookup;

  if (lookup.empty()) {
    for (int i = 0; i < ElementTableSize; i++)
      lookup[ElementTable[i].symbol] = i;

    lookup["Q"] = cAN_H;
    lookup["D"] = cAN_H;
  }

  // check second letter for lower case
  if (symbol[0] && isupper(symbol[1]) && strcmp(symbol, "LP") != 0) {
    UtilNCopy(titleized, symbol, 4);
    titleized[1] = tolower(symbol[1]);
    symbol = titleized;
  }

  // find in lookup dictionary
  auto it = lookup.find(symbol);

  if (it != lookup.end()) {
    return it->second;
  } else {
    // allow wacky names for C and H
    switch(symbol[0]) {
    case 'C':
      return cAN_C;
    case 'H':
      return cAN_H;
    }
  }

  return -1;
}

/**
 * Assign "protons" from "elem" or "name" property
 */
static void set_protons(PyMOLGlobals * G, AtomInfoType * I)
{
  int protons = get_protons(I->elem);

  if (protons < 0) {
    // try again with atom name, skip numbers
    const char * name = LexStr(G, I->name);
    while((*name >= '0') && (*name <= '9') && (*(name + 1)))
      name++;

    protons = get_protons(name);
  }

  I->protons = protons;
}

/**
 * Assign (based on name, elem, protons):
 *  - elem      (if empty string)
 *  - protons   (if < 1)
 *  - vdw       (if 0.0)
 *  - priority
 */
void AtomInfoAssignParameters(PyMOLGlobals * G, AtomInfoType * I)
{
  const char *n = NULL;
  char *e = NULL;
  int pri;

  e = I->elem;

  // elem from protons (except for LP, assume protons=0 if uninitialized)
  if(!*e && I->protons > 0) {
    atomicnumber2elem(e, I->protons);
  }

  // elem from name
  if(!*e) {                     /* try to guess atomic type from name */
    n = LexStr(G, I->name);
    while(((*n) >= '0') && ((*n) <= '9') && (*(n + 1)))
      n++;
    strncpy(e, n, cElemNameLen);
    switch (*e) {
    case '\0':
      break;
    case 'C':
      if(*(e + 1) == 'A') {
        if(!WordMatchExact(G, G->lex_const.CA, I->resn, true)
           && (!WordMatchExact(G, "CA+", LexStr(G, I->resn), true)))
          /* CA intpreted as carbon, not calcium */
          *(e + 1) = 0;
      } else if(!((*(e + 1) == 'a') ||  /* CA intpreted as carbon, not calcium */
                  (*(e + 1) == 'l') || (*(e + 1) == 'L') ||
                  (*(e + 1) == 'u') || (*(e + 1) == 'U') ||
                  (*(e + 1) == 'o') || (*(e + 1) == 'O') ||
                  (*(e + 1) == 's') || (*(e + 1) == 'S') ||
                  (*(e + 1) == 'r') || (*(e + 1) == 'R')
                ))
        *(e + 1) = 0;
      break;
    case 'H':
      if(!((*(e + 1) == 'e')
         ))
        *(e + 1) = 0;
      break;
    case 'D':                  /* take deuterium to hydrogen */
      *(e + 1) = 0;
      break;
    case 'N':
      if(!((*(e + 1) == 'i') || (*(e + 1) == 'I') ||
           (*(e + 1) == 'a') || (*(e + 1) == 'A') ||
           (*(e + 1) == 'b') || (*(e + 1) == 'B')
         ))
        *(e + 1) = 0;
      break;
    case 'S':
      if(!((*(e + 1) == 'e') || (*(e + 1) == 'E') ||
           (*(e + 1) == 'r') || (*(e + 1) == 'R') ||
           (*(e + 1) == 'c') || (*(e + 1) == 'C') ||
           (*(e + 1) == 'b') || (*(e + 1) == 'B')
         ))
        *(e + 1) = 0;

      break;
    case 'O':
      if(!((*(e + 1) == 's')))
        *(e + 1) = 0;
      break;
    case 'Q':
      *(e + 1) = 0;
      break;
    case 'p':
      if (p_strstartswith(n, "pseudo")) {
        strcpy(e, "PS");
      }
    }
    if(*(e + 1) &&
        (/* e != "LP" */ e[1] != 'P' || e[0] != 'L') &&
        (/* e != "PS" */ e[1] != 'S' || e[0] != 'P'))
      *(e + 1) = tolower(*(e + 1));
  }

  // priority
  n = LexStr(G, I->name);
  while((*n >= '0') && (*n <= '9') && (*(n + 1)))
    n++;
  if(toupper(*n) != I->elem[0]) {
    pri = 1000;                 /* unconventional atom name -- make no assignments */
  } else if(SettingGetGlobal_b(G, cSetting_pdb_standard_order)) {
    switch (*n) {

    case 'N':
    case 'C':
    case 'O':
    case 'S':
      switch (*(n + 1)) {
      case 0:
        switch (*n) {
        case 'N':
          pri = 1;
          break;
        case 'C':
          pri = 3;
          break;
        case 'O':
          pri = 4;
          break;
        default:
          pri = 1000;
          break;
        }
        break;
      case 'A':
        switch (*n) {
        case 'C':
          pri = 2;
          break;
        default:
          pri = 5;
          break;                /* generic alpha atom */
        }
        break;
      case 'B':
        pri = 6;
        break;                  /* generic beta atom, etc. */
      case 'G':
        pri = 7;
        break;
      case 'D':
        pri = 8;
        break;
      case 'E':
        pri = 9;
        break;
      case 'Z':
        pri = 10;
        break;
      case 'H':
        pri = 11;
        break;
      case 'I':
        pri = 12;
        break;
      case 'J':
        pri = 13;
        break;
      case 'K':
        pri = 14;
        break;
      case 'L':
        pri = 15;
        break;
      case 'M':
        pri = 16;
        break;
      case 'N':
        pri = 17;
        break;
      case 'X':
        switch (*(n + 2)) {
        case 'T':
          pri = 999;
          break;
        default:
        case 0:
          pri = 16;
          break;
        }
      case '0':
      case '1':
      case '2':
      case '3':
      case '4':
      case '5':
      case '6':
      case '7':
      case '8':
      case '9':
        pri = 0;
        n++;
        while(*n) {
          if(*n == 'P') {
            pri -= 200;
            break;
          } else if((*n == '*') || (*n == '\'')) {
            pri = (-100) - pri;
            break;
          } else if((*n < '0') || (*n > '9'))
            break;
          pri *= 10;
          pri += (*n - '0');
          n++;
        }
        pri += 300;
        break;
      default:
        pri = 500;
        break;
      }
      break;
    case 'P':                  /* this will place the phosphate before CNO numbered atoms */
      pri = 20;
      break;
    case 'D':
    case 'H':
      switch (*(n + 1)) {
      case 0:
        pri = 1001;
        break;
      case 'A':
      case 'B':
        pri = 1003;
        break;
      case 'G':
        pri = 1004;
        break;
      case 'D':
        pri = 1005;
        break;
      case 'E':
        pri = 1006;
        break;
      case 'Z':
        pri = 1007;
        break;
      case 'H':
        pri = 1008;
        break;
      case 'I':
        pri = 1009;
        break;
      case 'J':
        pri = 1010;
        break;
      case 'K':
        pri = 1011;
        break;
      case 'L':
        pri = 1012;
        break;
      case 'M':
        pri = 1013;
        break;
      case 'N':
        pri = 1002;
        break;
      case 'X':
        pri = 1999;
        break;
      case '0':
      case '1':
      case '2':
      case '3':
      case '4':
      case '5':
      case '6':
      case '7':
      case '8':
      case '9':
        pri = 1020;
        n++;
        while(*n) {
          pri *= 10;
          pri += (*n - '0');
          n++;
        }
        pri += 25;
        break;

      default:
        pri = 1500;
        break;
      }
      break;
    default:
      pri = 1000;
      break;
    }
  } else {
    switch (*n) {
    case 'N':
    case 'C':
    case 'O':
    case 'S':
      switch (*(n + 1)) {
      case 0:
        switch (*n) {
        case 'N':
          pri = 1;
          break;
        case 'C':
          pri = 997;
          break;
        case 'O':
          pri = 998;
          break;
        default:
          pri = 1000;
          break;
        }
        break;
      case 'A':
        pri = 3;
        break;                  /* generic alpha */
      case 'B':
        pri = 4;
        break;
      case 'G':
        pri = 5;
        break;
      case 'D':
        pri = 6;
        break;
      case 'E':
        pri = 7;
        break;
      case 'Z':
        pri = 8;
        break;
      case 'H':
        pri = 9;
        break;
      case 'I':
        pri = 10;
        break;
      case 'J':
        pri = 11;
        break;
      case 'K':
        pri = 12;
        break;
      case 'L':
        pri = 13;
        break;
      case 'M':
        pri = 14;
        break;
      case 'N':
        pri = 15;
        break;
      case 'X':
        switch (*(n + 2)) {
        case 'T':
          pri = 999;
          break;
        default:
        case 0:
          pri = 16;
          break;
        }
      case '0':
      case '1':
      case '2':
      case '3':
      case '4':
      case '5':
      case '6':
      case '7':
      case '8':
      case '9':
        pri = 0;
        n++;
        while(*n) {
          pri *= 10;
          pri += (*n - '0');
          n++;
        }
        pri += 25;
        break;
      default:
        pri = 500;
        break;
      }
      break;
    default:
      pri = 1000;
      break;
    }
  }

  I->priority = pri;

  // protons from elem or name
  if (I->protons < 1)
    set_protons(G, I);

  // vdw from protons
  if(I->vdw == 0.0) {
    if (I->protons > -1 && I->protons < ElementTableSize) {
      I->vdw = ElementTable[I->protons].vdw;
    } else {
      I->vdw = 1.80F;
    }
  }
}

/**
 * Get the element symbol. E.g. "He" for Helium.
 * @param[out] dst output buffer of length cElemNameLen
 * @param protons atomic number, e.g. 2 for Helium
 */
void atomicnumber2elem(char * dst, int protons) {
  if (protons > -1 && protons < ElementTableSize)
    strncpy(dst, ElementTable[protons].symbol, cElemNameLen);
}

/**
 * Get a string representation of the stereo configuration. Either
 * R/S chirality, if set, or the SDF odd/even parity, if set.
 */
const char * AtomInfoGetStereoAsStr(const AtomInfoType * ai) {
  switch (ai->mmstereo) {
    case 1 /* MMSTEREO_CHIRALITY_R */: return "R";
    case 2 /* MMSTEREO_CHIRALITY_S */: return "S";
  }

  switch (ai->stereo) {
    case SDF_CHIRALITY_ODD:  return "odd";
    case SDF_CHIRALITY_EVEN: return "even";
  }

  if (ai->mmstereo || ai->stereo) {
    return "?";
  }

  return "";
}

/**
 * Set stereochemistry. Valid are: R, S, N[one], E[ven], O[dd]
 */
void AtomInfoSetStereo(AtomInfoType * ai, const char * stereo) {
  switch (toupper(stereo[0])) {
    case 'R':  ai->mmstereo = 1; ai->stereo = 0; break; // MMSTEREO_CHIRALITY_R
    case 'S':  ai->mmstereo = 2; ai->stereo = 0; break; // MMSTEREO_CHIRALITY_S
    case 'E':  ai->mmstereo = 0; ai->stereo = SDF_CHIRALITY_EVEN;       break;
    case 'O':  ai->mmstereo = 0; ai->stereo = SDF_CHIRALITY_ODD;        break;
    case 'A': // ANS (s), ANR (r) pseudochirality
    case 'N': case 0: ai->mmstereo = ai->stereo = 0;    break;
    default:   ai->mmstereo = ai->stereo = 3;           break;
  }
}


/**
 * Get column aligned (left space padded) PDB residue name
 *
 * @param[out] resn output buffer
 */
void AtomInfoGetAlignedPDBResidueName(PyMOLGlobals * G,
    const AtomInfoType * ai,
    ResName & resn)
{
  sprintf(resn, "%3.4s", LexStr(G, ai->resn));
  if(SettingGetGlobal_b(G, cSetting_pdb_truncate_residue_name)) {
    resn[3] = 0;                /* enforce 3-letter residue name in PDB files */
  }
}


/**
 * Get column aligned (left space padded) PDB atom name
 *
 * @param resn space padded residue name
 * @param[out] name output buffer
 */
void AtomInfoGetAlignedPDBAtomName(PyMOLGlobals * G,
    const AtomInfoType * ai,
    const ResName & resn,
    AtomName & name)
{
  int literal = SettingGetGlobal_b(G, cSetting_pdb_literal_names);
  int reformat = SettingGetGlobal_i(G, cSetting_pdb_reformat_names_mode);

  // default is "literal=0" and "reformat=0" (no reformatting)

  const char * ai_name = LexStr(G, ai->name);
  auto ai_name_len = strlen(ai_name);
  bool start_column_1 = false;

  UtilNCopy(name, ai_name, 5);

  if(!ai->name) {
    if(!ai->elem[1])
      sprintf(name, " %s", ai->elem);
    else
      sprintf(name, "%s", ai->elem);
  } else if(!literal) {
    if(ai_name_len < 4) {  /* atom name less than length 4 */
      if(!isdigit(name[0])) {     /* doesn't start with a number */
        if((toupper(ai->elem[0]) == toupper(name[0])) && ((!ai->elem[1]) || /* symbol len = 1 */
              (toupper(ai->elem[1]) == toupper(name[1])))) {        /* matched len 2 */
          /* starts with corrent atomic symbol, so */
          if(!ai->elem[1]) { /* symbol len = 1 */
            switch (reformat) {
            case 1:            /* pdb with internal pdb */
            case 3:            /* pdb with internal iupac */
              if((ai->elem[0] == 'H') && ai_name_len > 2) {
                AtomInfoGetPDB3LetHydroName(G, resn, ai_name, name);
                break;
              }
            default:           /* otherwise, start in column 1 */
              start_column_1 = true;
              break;
            }
          }
        } else {                /* name doesn't start with atomic symbol */
          /* then just place it in column 1 as usual */
          start_column_1 = true;
        }
      } else {                  /* name starts with a number */
        switch (reformat) {
        case 2:                /* make Amber compliant */
          if((ai->elem[0] == name[1]) &&
             ((!ai->elem[1]) || (toupper(ai->elem[1]) == toupper(name[2])))) {
            /* rotate the name to place atom symbol in column 0 to comply with Amber PDB format */
            name[3] = name[0];
            name[0] = ' ';
          }
          break;
        }
      }                         /* just stick it in column 0 and hope for the best */
    } else {                    /* if name is length 4 */
      if((ai->elem[0] == name[0]) && ((!ai->elem[1]) ||     /* symbol len = 1 */
                                          (toupper(ai->elem[1]) == toupper(name[1])))) {    /* matched len 2 */
        /* name starts with the atomic symbol */
        if((!ai->elem[1]) && (ai->elem[0])) {   /* but if element is one letter... */
          switch (reformat) {
          case 1:              /* retaining PDB compliance throughout, or */
          case 3:              /* saving as PDB compliant, but use IUPAC within PyMOL */
            if(isdigit(name[3])) {  /* and last character is a number */
              /* rotate the name to place atom symbol in column 1 to comply with PDB format */
              name[0] = ai_name[3];
              name[1] = ai_name[0];
              name[2] = ai_name[1];
              name[3] = ai_name[2];
              name[4] = 0;
            }
            break;
          }
        }
      } else {                  /* name does not start with the symbol... */
        if(reformat == 2) {     /* AMBER compliance mode */
          if(isdigit(name[0])) {
            if((ai->elem[0] == name[1]) &&
               ((!(ai->elem[1])) || (toupper(ai->elem[1]) == toupper(name[2])))) {
              /* rotate the name to place atom symbol in column 0 to comply with Amber PDB format */
              name[0] = ai_name[1];
              name[1] = ai_name[2];
              name[2] = ai_name[3];
              name[3] = ai_name[0];
              name[4] = 0;
            }
          }
        }
      }
    }
  } else {                      /* LITERAL mode: preserve what was in the original PDB as best PyMOL can
                                   this should enable people to open and save amber pdb files without issues */
    if (ai_name_len < 4 && !(ai->elem[1] && /* elem len = 2 */
          toupper(ai->elem[0]) == toupper(name[0]) &&
          toupper(ai->elem[1]) == toupper(name[1]))) {
      start_column_1 = true;
    }
  }

  if (start_column_1) {
    name[0] = ' ';
    UtilNCopy(name + 1, ai_name, 4);
  }

  name[4] = 0;
}

bool BondType::hasSymOp() const
{
  return symop_2;
}
