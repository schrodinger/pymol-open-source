
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

#include"os_predef.h"
#include"os_std.h"

#include"AtomInfo.h"
#include"Word.h"
#include"MemoryDebug.h"
#include"Err.h"
#include"Feedback.h"
#include"Util.h"
#include"Color.h"
#include"PConv.h"
#include"Ortho.h"
#include"OVOneToAny.h"
#include"OVContext.h"
#include"PyMOLObject.h"

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
  int c = 0;
  while(*p) {
    if(c + 1 == sizeof(AtomName)) {
      break;
    }
    if((((*p) >= '0') && ((*p) <= '9')) ||
       (((*p) >= 'a') && ((*p) <= 'z')) ||
       (((*p) >= 'A') && ((*p) <= 'Z')) ||
       ((*p) == '.') ||
       ((*p) == '_') || ((*p) == '+') || ((*p) == '\'') || ((*p) == '*')) {
      *q++ = *p;
      c++;
    }
    p++;
  }
  *q = 0;
}

int AtomInfoCheckSetting(PyMOLGlobals * G, AtomInfoType * ai, int setting_id)
{
  if(!ai->has_setting) {
    return 0;
  } else {
    if(!SettingUniqueCheck(G, ai->unique_id, setting_id)) {
      return 0;
    } else {
      return 1;
    }
  }
}

int AtomInfoGetSetting_b(PyMOLGlobals * G, AtomInfoType * ai, int setting_id, int current,
                         int *effective)
{
  if(!ai->has_setting) {
    *effective = current;
    return 0;
  } else {
    if(!SettingUniqueGet_b(G, ai->unique_id, setting_id, effective)) {
      *effective = current;
      return 0;
    } else {
      return 1;
    }
  }
}

int AtomInfoGetSetting_i(PyMOLGlobals * G, AtomInfoType * ai, int setting_id, int current,
                         int *effective)
{
  if(!ai->has_setting) {
    *effective = current;
    return 0;
  } else {
    if(!SettingUniqueGet_i(G, ai->unique_id, setting_id, effective)) {
      *effective = current;
      return 0;
    } else {
      return 1;
    }
  }
}

int AtomInfoGetSetting_f(PyMOLGlobals * G, AtomInfoType * ai, int setting_id,
                         float current, float *effective)
{
  if(!ai->has_setting) {
    *effective = current;
    return 0;
  } else {
    if(!SettingUniqueGet_f(G, ai->unique_id, setting_id, effective)) {
      *effective = current;
      return 0;
    } else {
      return 1;
    }
  }
}

int AtomInfoGetSetting_color(PyMOLGlobals * G, AtomInfoType * ai, int setting_id,
                             int current, int *effective)
{
  if(!ai->has_setting) {
    *effective = current;
    return 0;
  } else {
    if(!SettingUniqueGet_color(G, ai->unique_id, setting_id, effective)) {
      *effective = current;
      return 0;
    } else {
      return 1;
    }
  }
}

int AtomInfoCheckBondSetting(PyMOLGlobals * G, BondType * bi, int setting_id)
{
  if(!bi->has_setting) {
    return 0;
  } else {
    if(!SettingUniqueCheck(G, bi->unique_id, setting_id)) {
      return 0;
    } else {
      return 1;
    }
  }
}

int AtomInfoGetBondSetting_b(PyMOLGlobals * G, BondType * bi, int setting_id, int current,
                             int *effective)
{
  if(!bi->has_setting) {
    *effective = current;
    return 0;
  } else {
    if(!SettingUniqueGet_b(G, bi->unique_id, setting_id, effective)) {
      *effective = current;
      return 0;
    } else {
      return 1;
    }
  }
}

int AtomInfoGetBondSetting_i(PyMOLGlobals * G, BondType * bi, int setting_id, int current,
                             int *effective)
{
  if(!bi->has_setting) {
    *effective = current;
    return 0;
  } else {
    if(!SettingUniqueGet_i(G, bi->unique_id, setting_id, effective)) {
      *effective = current;
      return 0;
    } else {
      return 1;
    }
  }
}

int AtomInfoGetBondSetting_f(PyMOLGlobals * G, BondType * bi, int setting_id,
                             float current, float *effective)
{
  if(!bi->has_setting) {
    *effective = current;
    return 0;
  } else {
    if(!SettingUniqueGet_f(G, bi->unique_id, setting_id, effective)) {
      *effective = current;
      return 0;
    } else {
      return 1;
    }
  }
}

int AtomInfoGetBondSetting_color(PyMOLGlobals * G, BondType * bi, int setting_id,
                                 int current, int *effective)
{
  if(!bi->has_setting) {
    *effective = current;
    return 0;
  } else {
    if(!SettingUniqueGet_color(G, bi->unique_id, setting_id, effective)) {
      *effective = current;
      return 0;
    } else {
      return 1;
    }
  }
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

int AtomInfoInit(PyMOLGlobals * G)
{
  register CAtomInfo *I = NULL;
  if((I = (G->AtomInfo = Calloc(CAtomInfo, 1)))) {
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
void AtomInfoGetPDB3LetHydroName(PyMOLGlobals * G, char *resn, char *iname, char *oname)
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

int AtomInfoKnownWaterResName(PyMOLGlobals * G, char *resn)
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

int AtomInfoKnownPolymerResName(char *resn)
{
  switch (resn[0]) {
  case 'A':
    switch (resn[1]) {
    case 0: /*  A*/
      return true;
      break;
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
    case 0: /* C */
      return true;
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
  case 'D':
    switch (resn[1]) {
    case 'G': /* DG */
      switch (resn[2]) {
      case 0:
        return true;
        break;
      }
    case 'C': /* DC */
      switch (resn[2]) {
      case 0:
        return true;
        break;
      }
    case 'T': /* DT */
      switch (resn[2]) {
      case 0:
        return true;
        break;
      }
    case 'A': /* DA */
      switch (resn[2]) {
      case 0:
        return true;
        break;
      }
    }
    break;
  case 'G':
    switch (resn[1]) {
    case 0: /* G */
      return true;
      break;
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
    case 0: /* T */
      return true;
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
  case 'U':
    switch (resn[1]) {
    case 0: /* U */
      return true;
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

static int AtomInfoInOrder(PyMOLGlobals * G, AtomInfoType * atom, int atom1, int atom2)
{
  return (AtomInfoCompare(G, atom + atom1, atom + atom2) <= 0);
}

static int AtomInfoInOrderIgnoreHet(PyMOLGlobals * G, AtomInfoType * atom, int atom1,
                                    int atom2)
{
  return (AtomInfoCompareIgnoreHet(G, atom + atom1, atom + atom2) <= 0);
}

static int AtomInfoInOrigOrder(PyMOLGlobals * G, AtomInfoType * atom, int atom1,
                               int atom2)
{
  if(atom[atom1].rank != atom[atom2].rank)
    return (atom[atom1].rank < atom[atom2].rank);
  else
    return (AtomInfoCompare(G, atom + atom1, atom + atom2) <= 0);
}

int AtomResvFromResi(char *resi)
{
  int result = 1;
  char *start = resi;
  while(*start) {
    if(sscanf(start, "%d", &result))
      break;
    else
      result = 1;
    start++;
  }
  return result;
}


/*========================================================================*/
PyObject *AtomInfoAsPyList(PyMOLGlobals * G, AtomInfoType * I)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else
  PyObject *result = NULL;

  result = PyList_New(47);
  PyList_SetItem(result, 0, PyInt_FromLong(I->resv));
  PyList_SetItem(result, 1, PyString_FromString(I->chain));
  PyList_SetItem(result, 2, PyString_FromString(I->alt));
  PyList_SetItem(result, 3, PyString_FromString(I->resi));
  PyList_SetItem(result, 4, PyString_FromString(I->segi));
  PyList_SetItem(result, 5, PyString_FromString(I->resn));
  PyList_SetItem(result, 6, PyString_FromString(I->name));
  PyList_SetItem(result, 7, PyString_FromString(I->elem));

  {
    char null_st[1] = "";
    char *st = null_st;

    if(I->textType)
      st = OVLexicon_FetchCString(G->Lexicon, I->textType);
    PyList_SetItem(result, 8, PyString_FromString(st));

    st = null_st;
    if(I->label)
      st = OVLexicon_FetchCString(G->Lexicon, I->label);
    PyList_SetItem(result, 9, PyString_FromString(st));
  }

  PyList_SetItem(result, 10, PyString_FromString(I->ssType));
  PyList_SetItem(result, 11, PyInt_FromLong((char) I->hydrogen));
  PyList_SetItem(result, 12, PyInt_FromLong(I->customType));
  PyList_SetItem(result, 13, PyInt_FromLong(I->priority));
  PyList_SetItem(result, 14, PyFloat_FromDouble(I->b));
  PyList_SetItem(result, 15, PyFloat_FromDouble(I->q));
  PyList_SetItem(result, 16, PyFloat_FromDouble(I->vdw));
  PyList_SetItem(result, 17, PyFloat_FromDouble(I->partialCharge));
  PyList_SetItem(result, 18, PyInt_FromLong(I->formalCharge));
  PyList_SetItem(result, 19, PyInt_FromLong((int) I->hetatm));
  PyList_SetItem(result, 20, PConvSCharArrayToPyList(I->visRep, cRepCnt));
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
  PyList_SetItem(result, 39, PyInt_FromLong((int) I->atomic_color));
  PyList_SetItem(result, 40, PyInt_FromLong((int) I->has_setting));
  PyList_SetItem(result, 41, PyFloat_FromDouble(I->U11));
  PyList_SetItem(result, 42, PyFloat_FromDouble(I->U22));
  PyList_SetItem(result, 43, PyFloat_FromDouble(I->U33));
  PyList_SetItem(result, 44, PyFloat_FromDouble(I->U12));
  PyList_SetItem(result, 45, PyFloat_FromDouble(I->U13));
  PyList_SetItem(result, 46, PyFloat_FromDouble(I->U23));

  return (PConvAutoNone(result));
#endif
}

int AtomInfoFromPyList(PyMOLGlobals * G, AtomInfoType * I, PyObject * list)
{
#ifdef _PYMOL_NOPY
  return 0;
#else
  int ok = true;
  int hetatm;
  ov_size ll = 0;
  if(ok)
    ok = PyList_Check(list);
  if(ok)
    ll = PyList_Size(list);
  if(ok)
    ok = PConvPyIntToInt(PyList_GetItem(list, 0), &I->resv);
  if(ok)
    ok = PConvPyStrToStr(PyList_GetItem(list, 1), I->chain, sizeof(Chain));
  if(ok)
    ok = PConvPyStrToStr(PyList_GetItem(list, 2), I->alt, sizeof(Chain));
  if(ok)
    ok = PConvPyStrToStr(PyList_GetItem(list, 3), I->resi, sizeof(ResIdent));
  if(ok)
    ok = PConvPyStrToStr(PyList_GetItem(list, 4), I->segi, sizeof(SegIdent));
  if(ok)
    ok = PConvPyStrToStr(PyList_GetItem(list, 5), I->resn, sizeof(ResName));
  if(ok)
    ok = PConvPyStrToStr(PyList_GetItem(list, 6), I->name, sizeof(AtomName));
  if(ok)
    ok = PConvPyStrToStr(PyList_GetItem(list, 7), I->elem, sizeof(ElemName));
  if(ok) {
    OrthoLineType temp;
    PConvPyStrToStr(PyList_GetItem(list, 8), temp, sizeof(OrthoLineType));
    I->textType = 0;
    if(temp[0]) {
      OVreturn_word result = OVLexicon_GetFromCString(G->Lexicon, temp);
      if(OVreturn_IS_OK(result)) {
        I->textType = result.word;
      }
    }
  }
  if(ok) {
    OrthoLineType temp;
    PConvPyStrToStr(PyList_GetItem(list, 9), temp, sizeof(OrthoLineType));
    I->textType = 0;
    if(temp[0]) {
      OVreturn_word result = OVLexicon_GetFromCString(G->Lexicon, temp);
      if(OVreturn_IS_OK(result)) {
        I->label = result.word;
      }
    }
  }
  if(ok)
    ok = PConvPyStrToStr(PyList_GetItem(list, 10), I->ssType, sizeof(SSType));
  if(ok)
    ok = PConvPyIntToChar(PyList_GetItem(list, 11), (char *) &I->hydrogen);
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
    ok = PConvPyIntToInt(PyList_GetItem(list, 18), &I->formalCharge);
  if(ok)
    ok = PConvPyIntToInt(PyList_GetItem(list, 19), &hetatm);
  if(ok)
    I->hetatm = hetatm;
  if(ok)
    ok =
      PConvPyListToSCharArrayInPlaceAutoZero(PyList_GetItem(list, 20), I->visRep,
                                             cRepCnt);
  if(ok)
    ok = PConvPyIntToInt(PyList_GetItem(list, 21), &I->color);
  if(ok)
    I->color = ColorConvertOldSessionIndex(G, I->color);
  if(ok)
    ok = PConvPyIntToInt(PyList_GetItem(list, 22), &I->id);
  if(ok)
    ok = PConvPyIntToChar(PyList_GetItem(list, 23), (char *) &I->cartoon);
  if(ok)
    ok = PConvPyIntToInt(PyList_GetItem(list, 24), (int *) &I->flags);
  if(ok)
    ok = PConvPyIntToChar(PyList_GetItem(list, 25), (char *) &I->bonded);
  if(ok)
    ok = PConvPyIntToChar(PyList_GetItem(list, 26), (char *) &I->chemFlag);
  if(ok)
    ok = PConvPyIntToChar(PyList_GetItem(list, 27), (char *) &I->geom);
  if(ok)
    ok = PConvPyIntToChar(PyList_GetItem(list, 28), (char *) &I->valence);
  if(ok)
    ok = PConvPyIntToChar(PyList_GetItem(list, 29), (char *) &I->masked);
  if(ok)
    ok = PConvPyIntToChar(PyList_GetItem(list, 30), (char *) &I->protekted);
  if(ok)
    ok = PConvPyIntToChar(PyList_GetItem(list, 31), (char *) &I->protons);
  if(ok)
    ok = PConvPyIntToInt(PyList_GetItem(list, 32), &I->unique_id);
  if(ok && I->unique_id) {      /* reserve existing IDs */
    CAtomInfo *II = G->AtomInfo;
    AtomInfoPrimeUniqueIDs(G);
    I->unique_id = SettingUniqueConvertOldSessionID(G, I->unique_id);
    OVOneToAny_SetKey(II->ActiveIDs, I->unique_id, 1);
  }
  if(ok)
    ok = PConvPyIntToChar(PyList_GetItem(list, 33), (char *) &I->stereo);
  if(ok && (ll > 34))
    ok = PConvPyIntToInt(PyList_GetItem(list, 34), &I->discrete_state);
  if(ok && (ll > 35))
    ok = PConvPyFloatToFloat(PyList_GetItem(list, 35), &I->elec_radius);
  if(ok && (ll > 36))
    ok = PConvPyIntToInt(PyList_GetItem(list, 36), &I->rank);
  if(ok && (ll > 37))
    ok = PConvPyIntToChar(PyList_GetItem(list, 37), (char *) &I->hb_donor);
  if(ok && (ll > 38))
    ok = PConvPyIntToChar(PyList_GetItem(list, 38), (char *) &I->hb_acceptor);
  if(ok && (ll > 39)) {
    ok = PConvPyIntToInt(PyList_GetItem(list, 39), &I->atomic_color);
    if(ok)
      I->atomic_color = ColorConvertOldSessionIndex(G, I->atomic_color);
  } else {
    I->atomic_color = AtomInfoGetColor(G, I);
  }
  if(ok && (ll > 40))
    ok = PConvPyIntToChar(PyList_GetItem(list, 40), (char *) &I->has_setting);
  if(ok && (ll > 46)) {
    if(ok)
      ok = PConvPyFloatToFloat(PyList_GetItem(list, 41), &I->U11);
    if(ok)
      ok = PConvPyFloatToFloat(PyList_GetItem(list, 42), &I->U22);
    if(ok)
      ok = PConvPyFloatToFloat(PyList_GetItem(list, 43), &I->U33);
    if(ok)
      ok = PConvPyFloatToFloat(PyList_GetItem(list, 44), &I->U12);
    if(ok)
      ok = PConvPyFloatToFloat(PyList_GetItem(list, 45), &I->U13);
    if(ok)
      ok = PConvPyFloatToFloat(PyList_GetItem(list, 46), &I->U23);
  }
  return (ok);
#endif
}

void AtomInfoCopy(PyMOLGlobals * G, AtomInfoType * src, AtomInfoType * dst)
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
  if(dst->label) {
    OVLexicon_IncRef(G->Lexicon, dst->label);
  }
  if(dst->textType) {
    OVLexicon_IncRef(G->Lexicon, dst->textType);
  }
}

void AtomInfoBondCopy(PyMOLGlobals * G, BondType * src, BondType * dst)
{
  *(dst) = *(src);

  if(src->unique_id && src->has_setting) {
    dst->unique_id = AtomInfoGetNewUniqueID(G);
    if(!SettingUniqueCopyAll(G, src->unique_id, dst->unique_id))
      dst->has_setting = 0;
  } else {
    src->unique_id = 0;
    src->has_setting = 0;
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
  if(ai->textType) {
    OVLexicon_DecRef(G->Lexicon, ai->textType);
  }
  if(ai->has_setting && ai->unique_id) {
    SettingUniqueDetachChain(G, ai->unique_id);
  }
  if(ai->unique_id && I->ActiveIDs) {
    OVOneToAny_DelKey(I->ActiveIDs, ai->unique_id);
  }
  if(ai->label) {
    /*    printf("purging %d [%s]\n", OVLexicon_GetNActive(G->Lexicon),
       OVLexicon_FetchCString(G->Lexicon,ai->label)); */
    OVLexicon_DecRef(G->Lexicon, ai->label);
  }
}


/*========================================================================*/
void AtomInfoCombine(PyMOLGlobals * G, AtomInfoType * dst, AtomInfoType * src, int mask)
{
  if(mask & cAIC_tt) {
    if(dst->textType) {
      OVLexicon_DecRef(G->Lexicon, dst->textType);
      dst->textType = 0;
    }
    dst->textType = src->textType;
  } else if(src->textType) {
    OVLexicon_DecRef(G->Lexicon, src->textType);
    src->textType = 0;
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
  dst->unique_id = src->unique_id;
  /* keep all existing names, identifiers, etc. */
  /* also keep all existing selections,
     colors, masks, and visible representations */
  {
    if(src->label) {            /* destroy src label if one exists */
      OVLexicon_DecRef(G->Lexicon, src->label);
      src->label = 0;
    }
    /* leaves dst->label untouched */
  }
}


/*========================================================================*/
int AtomInfoUniquefyNames(PyMOLGlobals * G, AtomInfoType * atInfo0, int n0,
                          AtomInfoType * atInfo1, int *flag1, int n1)
{
  /* makes sure all names in atInfo1 are unique WRT 0 and 1 */

  /* tricky optimizations to avoid n^2 dependence in this operation */
  int result = 0;
  int a, b, c;
  AtomInfoType *ai0, *ai1, *lai1, *lai0;
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

    if(!ai1->name[0])
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
        if(!WordMatchExact(G, ai1->name, ai0->name, true))
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

        if(!lai0)
          bracketFlag = true;
        else if(!AtomInfoSameResidue(G, lai0, ai1))
          bracketFlag = true;
        else
          bracketFlag = false;
        if(bracketFlag) {
          AtomInfoBracketResidue(G, atInfo0, n0, ai1, &st0, &nd0);
          lai0 = ai1;
        }
        ai0 = atInfo0 + st0;
        for(a = st0; a <= nd0; a++) {
          if(!WordMatchExact(G, ai1->name, ai0->name, true))
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
      name[4] = 0;              /* just is case we go over */
      strcpy(ai1->name, name);
      result++;
      c = c + 1;
    } else {
      ai1++;
      b++;
    }
  }
  return result;
}


/*========================================================================*/
void AtomInfoBracketResidue(PyMOLGlobals * G, AtomInfoType * ai0, int n0,
                            AtomInfoType * ai, int *st, int *nd)
{
  /* inefficient but reliable way to find where residue atoms are located in an object 
   * for purpose of residue-based operations */
  int a;
  AtomInfoType *ai1;

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
void AtomInfoBracketResidueFast(PyMOLGlobals * G, AtomInfoType * ai0, int n0, int cur,
                                int *st, int *nd)
{
  /* efficient but unreliable way to find where residue atoms are located in an object 
   * for purpose of residue-based operations */
  int a;
  AtomInfoType *ai1;

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
  if(SettingGet(G, cSetting_auto_color))
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
float AtomInfoGetBondLength(PyMOLGlobals * G, AtomInfoType * ai1, AtomInfoType * ai2)
{
  float result = 1.6F;
  AtomInfoType *a1, *a2;

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

void AtomInfoAssignColors(PyMOLGlobals * G, AtomInfoType * at1)
{
  at1->atomic_color = (at1->color = AtomInfoGetColor(G, at1));
}

int AtomInfoGetColor(PyMOLGlobals * G, AtomInfoType * at1)
{
  CAtomInfo *I = G->AtomInfo;
  char *n = at1->elem;
  int color = I->DefaultColor;

  if(!n[0]) {
    n = at1->name;
    while(((*n) >= '0') && ((*n) <= '9') && (*(n + 1)))
      n++;
  }
  {
    register char c1 = n[0], c2 = n[1];
    c2 = tolower(c2);
    if(c2 <= 32)
      c2 = 0;
    switch (c1) {
    case 'A':
      switch (c2) {
      case 'c':
        color = ColorGetIndex(G, "actinium");
        break;
      case 'g':
        color = ColorGetIndex(G, "silver");
        break;
      case 'l':
        color = ColorGetIndex(G, "aluminum");
        break;
      case 'm':
        color = ColorGetIndex(G, "americium");
        break;
      case 'r':
        color = ColorGetIndex(G, "argon");
        break;
      case 's':
        color = ColorGetIndex(G, "arsenic");
        break;
      case 't':
        color = ColorGetIndex(G, "astatine");
        break;
      case 'u':
        color = ColorGetIndex(G, "gold");
        break;
      }
      break;
    case 'B':
      switch (c2) {
      case 0:
        color = ColorGetIndex(G, "boron");
        break;
      case 'a':
        color = ColorGetIndex(G, "barium");
        break;
      case 'e':
        color = ColorGetIndex(G, "beryllium");
        break;
      case 'h':
        color = ColorGetIndex(G, "bohrium");
        break;
      case 'i':
        color = ColorGetIndex(G, "bismuth");
        break;
      case 'k':
        color = ColorGetIndex(G, "berkelium");
        break;
      case 'r':
        color = I->BrColor;
        break;
      }
      break;
    case 'C':
      switch (c2) {
      case 0:
        color = I->CColor;
        break;
      case 'a':
        color = ColorGetIndex(G, "calcium");
        break;
      case 'd':
        color = ColorGetIndex(G, "cadmium");
        break;
      case 'e':
        color = ColorGetIndex(G, "cerium");
        break;
      case 'f':
        color = ColorGetIndex(G, "californium");
        break;
      case 'l':
        color = I->ClColor;
        break;
      case 'm':
        color = ColorGetIndex(G, "curium");
        break;
      case 'o':
        color = ColorGetIndex(G, "cobalt");
        break;
      case 'r':
        color = ColorGetIndex(G, "chromium");
        break;
      case 's':
        color = ColorGetIndex(G, "cesium");
        break;
      case 'u':
        color = I->CuColor;
        break;
      }
      break;
    case 'D':
      switch (c2) {
      case 0:
        color = I->DColor;
        break;
      case 'b':
        color = ColorGetIndex(G, "dubnium");
        break;
      case 'y':
        color = ColorGetIndex(G, "dysprosium");
        break;
      }
      break;
    case 'E':
      switch (c2) {
      case 'r':
        color = ColorGetIndex(G, "erbium");
        break;
      case 's':
        color = ColorGetIndex(G, "einsteinium");
        break;
      case 'u':
        color = ColorGetIndex(G, "europium");
        break;
      }
      break;
    case 'F':
      switch (c2) {
      case 0:
        color = I->FColor;
        break;
      case 'e':
        color = I->FeColor;
        break;
      case 'm':
        color = ColorGetIndex(G, "fermium");
        break;
      case 'r':
        color = ColorGetIndex(G, "francium");
        break;
      }
      break;
    case 'G':
      switch (c2) {
      case 'a':
        color = ColorGetIndex(G, "gallium");
        break;
      case 'd':
        color = ColorGetIndex(G, "gadolinium");
        break;
      case 'e':
        color = ColorGetIndex(G, "germanium");
        break;
      }
      break;
    case 'H':
      switch (c2) {
      case 0:
        color = I->HColor;
        break;
      case 'e':
        color = ColorGetIndex(G, "helium");
        break;
      case 'f':
        color = ColorGetIndex(G, "hafnium");
        break;
      case 'g':
        color = ColorGetIndex(G, "mercury");
        break;
      case 'o':
        color = ColorGetIndex(G, "holmium");
        break;
      case 's':
        color = ColorGetIndex(G, "hassium");
        break;
      }
      break;
    case 'I':
      switch (c2) {
      case 0:
        color = I->IColor;
        break;
      case 'n':
        color = ColorGetIndex(G, "indium");
        break;
      case 'r':
        color = ColorGetIndex(G, "iridium");
        break;
      }
      break;
    case 'K':
      switch (c2) {
      case 0:
        color = ColorGetIndex(G, "potassium");
        break;
      case 'r':
        color = ColorGetIndex(G, "krypton");
        break;
      }
      break;
    case 'L':
      switch (c2) {
      case 'a':
        color = ColorGetIndex(G, "lanthanum");
        break;
      case 'i':
        color = ColorGetIndex(G, "lithium");
        break;
      case 'r':
        color = ColorGetIndex(G, "lawrencium");
        break;
      case 'p':
        color = ColorGetIndex(G, "lonepair");
        break;
      case 'u':
        color = ColorGetIndex(G, "lutetium");
        break;
      }
      break;
    case 'M':
      switch (c2) {
      case 'd':
        color = ColorGetIndex(G, "mendelevium");
        break;
      case 'g':
        color = I->MgColor;
        break;
      case 'n':
        color = I->MnColor;
        break;
      case 'o':
        color = ColorGetIndex(G, "molybdenum");
        break;
      case 't':
        color = ColorGetIndex(G, "meitnerium");
        break;
      }
      break;
    case 'N':
      switch (c2) {
      case 0:
        color = I->NColor;
        break;
      case 'a':
        color = I->NaColor;
        break;
      case 'b':
        color = ColorGetIndex(G, "niobium");
        break;
      case 'd':
        color = ColorGetIndex(G, "neodymium");
        break;
      case 'e':
        color = ColorGetIndex(G, "neon");
        break;
      case 'i':
        color = ColorGetIndex(G, "nickel");
        break;
      case 'o':
        color = ColorGetIndex(G, "nobelium");
        break;
      case 'p':
        color = ColorGetIndex(G, "neptunium");
        break;
      }
      break;
    case 'O':
      switch (c2) {
      case 0:
        color = I->OColor;
        break;
      case 's':
        color = ColorGetIndex(G, "osmium");
        break;
      }
      break;
    case 'P':
      switch (c2) {
      case 0:
        color = I->PColor;
        break;
      case 'a':
        color = ColorGetIndex(G, "protactinium");
        break;
      case 'b':
        color = ColorGetIndex(G, "lead");
        break;
      case 'd':
        color = ColorGetIndex(G, "palladium");
        break;
      case 'm':
        color = ColorGetIndex(G, "promethium");
        break;
      case 'o':
        color = ColorGetIndex(G, "polonium");
        break;
      case 'r':
        color = ColorGetIndex(G, "praseodymium");
        break;
      case 's':
        color = ColorGetIndex(G, "pseudoatom");
        break;
      case 't':
        color = ColorGetIndex(G, "platinum");
        break;
      case 'u':
        color = ColorGetIndex(G, "plutonium");
        break;
      }
      break;
    case 'R':
      switch (c2) {
      case 'a':
        color = ColorGetIndex(G, "radium");
        break;
      case 'b':
        color = ColorGetIndex(G, "rubidium");
        break;
      case 'e':
        color = ColorGetIndex(G, "rhenium");
        break;
      case 'f':
        color = ColorGetIndex(G, "rutherfordium");
        break;
      case 'h':
        color = ColorGetIndex(G, "rhodium");
        break;
      case 'n':
        color = ColorGetIndex(G, "radon");
        break;
      case 'u':
        color = ColorGetIndex(G, "ruthenium");
        break;
      }
      break;
    case 'S':
      switch (c2) {
      case 0:
        color = I->SColor;
        break;
      case 'b':
        color = ColorGetIndex(G, "antimony");
        break;
      case 'c':
        color = ColorGetIndex(G, "scandium");
        break;
      case 'e':
        color = I->SeColor;
        break;
      case 'g':
        color = ColorGetIndex(G, "seaborgium");
        break;
      case 'i':
        color = ColorGetIndex(G, "silicon");
        break;
      case 'm':
        color = ColorGetIndex(G, "samarium");
        break;
      case 'n':
        color = ColorGetIndex(G, "tin");
        break;
      case 'r':
        color = ColorGetIndex(G, "strontium");
        break;
      }
      break;
    case 'T':
      switch (c2) {
      case 'a':
        color = ColorGetIndex(G, "tantalum");
        break;
      case 'b':
        color = ColorGetIndex(G, "terbium");
        break;
      case 'c':
        color = ColorGetIndex(G, "technetium");
        break;
      case 'e':
        color = ColorGetIndex(G, "tellurium");
        break;
      case 'h':
        color = ColorGetIndex(G, "thorium");
        break;
      case 'i':
        color = ColorGetIndex(G, "titanium");
        break;
      case 'l':
        color = ColorGetIndex(G, "thallium");
        break;
      case 'm':
        color = ColorGetIndex(G, "thulium");
        break;
      }
      break;
    case 'U':
      switch (c2) {
      case 0:
        color = ColorGetIndex(G, "uranium");
        break;
      }
      break;
    case 'V':
      switch (c2) {
      case 0:
        color = ColorGetIndex(G, "vanadium");
        break;
      }
      break;
    case 'W':
      switch (c2) {
      case 0:
        color = ColorGetIndex(G, "tungsten");
        break;
      }
      break;
    case 'X':
      switch (c2) {
      case 'e':
        color = ColorGetIndex(G, "xenon");
        break;
      }
      break;
    case 'Y':
      switch (c2) {
      case 0:
        color = ColorGetIndex(G, "yttrium");
        break;
      case 'b':
        color = ColorGetIndex(G, "ytterbium");
        break;
      }
      break;
    case 'Z':
      switch (c2) {
      case 'n':
        color = I->ZnColor;
        break;
      case 'r':
        color = ColorGetIndex(G, "zirconium");
        break;
      }
      break;
    }
  }
  return (color);
}

int *AtomInfoGetSortedIndex(PyMOLGlobals * G, CObject * obj, AtomInfoType * rec, int n,
                            int **outdex)
{
  int *index;
  int a;
  CSetting *setting = NULL;
  index = Alloc(int, n + 1);
  ErrChkPtr(G, index);
  (*outdex) = Alloc(int, n + 1);
  ErrChkPtr(G, *outdex);
  if(obj)
    setting = obj->Setting;

  if(SettingGet_b(G, setting, NULL, cSetting_retain_order)) {
    UtilSortIndexGlobals(G, n, rec, index, (UtilOrderFnGlobals *) AtomInfoInOrigOrder);
  } else if(SettingGet_b(G, setting, NULL, cSetting_pdb_hetatm_sort)) {
    UtilSortIndexGlobals(G, n, rec, index, (UtilOrderFnGlobals *) AtomInfoInOrder);
  } else {
    UtilSortIndexGlobals(G, n, rec, index,
                         (UtilOrderFnGlobals *) AtomInfoInOrderIgnoreHet);
  }

  for(a = 0; a < n; a++)
    (*outdex)[index[a]] = a;
  return (index);
}

void AtomInfoFreeSortedIndexes(PyMOLGlobals * G, int *index, int *outdex)
{
  FreeP(index);
  FreeP(outdex);
}

static int AtomInfoNameCompare(PyMOLGlobals * G, char *name1, char *name2)
{
  char *n1, *n2;
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

int AtomInfoCompare(PyMOLGlobals * G, AtomInfoType * at1, AtomInfoType * at2)
{
  int result;
  int wc;
  /* order by segment, chain, residue value, residue id, residue name, priority,
     and lastly by name */

  /*  printf(":%s:%s:%d:%s:%s:%i:%s =? ",
     at1->segi,at1->chain,at1->resv,at1->resi,at1->resn,at1->priority,at1->name);
     printf(":%s:%s:%d:%s:%s:%i:%s\n",
     at2->segi,at2->chain,at2->resv,at2->resi,at2->resn,at2->priority,at2->name);
   */

  wc = WordCompare(G, at1->segi, at2->segi, true);
  if(!wc) {
    if(at1->chain[0] == at2->chain[0]) {
      if(at1->hetatm == at2->hetatm) {
        if(at1->resv == at2->resv) {
          wc = WordCompare(G, at1->resi, at2->resi, true);
          if(!wc) {
            wc = WordCompare(G, at1->resn, at2->resn, true);
            if(!wc) {
              if(at1->discrete_state == at2->discrete_state) {
                if(at1->priority == at2->priority) {
                  if(at1->alt[0] == at2->alt[0]) {
                    if((result = AtomInfoNameCompare(G, at1->name, at2->name)) == 0) {
                      if(at1->rank < at2->rank)
                        result = -1;
                      else if(at1->rank > at2->rank)
                        result = 1;
                      else
                        result = 0;
                    }
                  } else if((!at2->alt[0])
                            || (at1->alt[0] && ((at1->alt[0] < at2->alt[0])))) {
                    result = -1;
                  } else {
                    result = 1;
                  }
                } else if(at1->priority < at2->priority) {
                  result = -1;
                } else {
                  result = 1;
                }
              } else if(at1->discrete_state < at2->discrete_state) {
                result = -1;
              } else {
                result = 1;
              }
            } else {
              result = wc;
            }
          } else {
            /* NOTE: don't forget to synchronize with below */
            if(SettingGetGlobal_b(G, cSetting_pdb_insertions_go_first)) {
              ov_size sl1, sl2;
              sl1 = strlen(at1->resi);
              sl2 = strlen(at2->resi);
              if(sl1 == sl2)
                result = wc;
              else if(sl1 < sl2)        /* sort residue 188A before 188, etc. */
                result = 1;
              else
                result = -1;
            } else if((at1->rank != at2->rank) &&
                      SettingGetGlobal_b(G, cSetting_rank_assisted_sorts)) {
              /* use rank to resolve insertion code ambiguities */
              if(at1->rank < at2->rank)
                result = -1;
              else
                result = 1;
            } else {
              result = wc;
            }
          }
        } else if(at1->resv < at2->resv) {
          result = -1;
        } else {
          result = 1;
        }
      } else if(at2->hetatm) {
        result = -1;
      } else {
        result = 1;
      }
    } else if((!at2->chain[0]) || (at1->chain[0] && ((at1->chain[0] < at2->chain[0])))) {
      result = -1;
    } else {
      result = 1;
    }
  } else {
    result = wc;
  }
  return (result);
}

int AtomInfoCompareIgnoreRankHet(PyMOLGlobals * G, AtomInfoType * at1, AtomInfoType * at2)
{
  register int result;
  register int wc = 0;
  /* order by segment, chain, residue value, residue id, residue name, priority,
     and lastly by name */

  /*  printf(":%s:%s:%d:%s:%s:%i:%s =? ",
     at1->segi,at1->chain,at1->resv,at1->resi,at1->resn,at1->priority,at1->name);
     printf(":%s:%s:%d:%s:%s:%i:%s\n",
     at2->segi,at2->chain,at2->resv,at2->resi,at2->resn,at2->priority,at2->name);
   */
  {
    /* attempt to optimize the segi comparison for equality */
    register char *p1, *p2;
    p1 = at1->segi;
    p2 = at2->segi;
    if((p1[0] != p2[0]) ||
       (p1[0] && ((p1[1] != p2[1]) ||
                  (p1[1] && ((p1[2] != p2[2]) || (p1[2] && (p1[3] != p2[3])))))))
      wc = WordCompare(G, p1, p2, true);
  }
  if(!wc) {
    if(at1->chain[0] == at2->chain[0]) {
      if(at1->resv == at2->resv) {
        wc = WordCompare(G, at1->resi, at2->resi, true);
        if(!wc) {
          wc = WordCompare(G, at1->resn, at2->resn, true);
          if(!wc) {
            if(at1->discrete_state == at2->discrete_state) {
              if(at1->priority == at2->priority) {
                if(at1->alt[0] == at2->alt[0]) {
                  result = AtomInfoNameCompare(G, at1->name, at2->name);
                } else if((!at2->alt[0])
                          || (at1->alt[0] && ((at1->alt[0] < at2->alt[0])))) {
                  result = -1;
                } else {
                  result = 1;
                }
              } else if(at1->priority < at2->priority) {
                result = -1;
              } else {
                result = 1;
              }
            } else if(at1->discrete_state < at2->discrete_state) {
              result = -1;
            } else {
              result = 1;
            }
          } else {
            result = wc;
          }
        } else {
          /* NOTE: don't forget to synchronize with below */
          if(SettingGetGlobal_b(G, cSetting_pdb_insertions_go_first)) {
            ov_size sl1, sl2;
            sl1 = strlen(at1->resi);
            sl2 = strlen(at2->resi);
            if(sl1 == sl2)
              result = wc;
            else if(sl1 < sl2)  /* sort residue 188A before 188, etc. */
              result = 1;
            else
              result = -1;
          } else if((at1->rank != at2->rank) &&
                    SettingGetGlobal_b(G, cSetting_rank_assisted_sorts)) {
            /* use rank to resolve insertion code ambiguities */
            if(at1->rank < at2->rank)
              result = -1;
            else
              result = 1;
          } else {
            result = wc;
          }
        }
      } else if(at1->resv < at2->resv) {
        result = -1;
      } else {
        result = 1;
      }
    } else if((!at2->chain[0]) || (at1->chain[0] && ((at1->chain[0] < at2->chain[0])))) {
      result = -1;
    } else {
      result = 1;
    }
  } else {
    result = wc;
  }
  return (result);
}

int AtomInfoCompareIgnoreRank(PyMOLGlobals * G, AtomInfoType * at1, AtomInfoType * at2)
{
  register int result;
  register int wc = 0;
  /* order by segment, chain, residue value, residue id, residue name, priority,
     and lastly by name */

  /*  printf(":%s:%s:%d:%s:%s:%i:%s =? ",
     at1->segi,at1->chain,at1->resv,at1->resi,at1->resn,at1->priority,at1->name);
     printf(":%s:%s:%d:%s:%s:%i:%s\n",
     at2->segi,at2->chain,at2->resv,at2->resi,at2->resn,at2->priority,at2->name);
   */
  {
    /* attempt to optimize the segi comparison for equality */
    register char *p1, *p2;
    p1 = at1->segi;
    p2 = at2->segi;
    if((p1[0] != p2[0]) ||
       (p1[0] && ((p1[1] != p2[1]) ||
                  (p1[1] && ((p1[2] != p2[2]) || (p1[2] && (p1[3] != p2[3])))))))
      wc = WordCompare(G, p1, p2, true);
  }
  if(!wc) {
    if(at1->chain[0] == at2->chain[0]) {
      if(at1->hetatm == at2->hetatm) {
        if(at1->resv == at2->resv) {
          wc = WordCompare(G, at1->resi, at2->resi, true);
          if(!wc) {
            wc = WordCompare(G, at1->resn, at2->resn, true);
            if(!wc) {
              if(at1->discrete_state == at2->discrete_state) {
                if(at1->priority == at2->priority) {
                  if(at1->alt[0] == at2->alt[0]) {
                    result = AtomInfoNameCompare(G, at1->name, at2->name);
                  } else if((!at2->alt[0])
                            || (at1->alt[0] && ((at1->alt[0] < at2->alt[0])))) {
                    result = -1;
                  } else {
                    result = 1;
                  }
                } else if(at1->priority < at2->priority) {
                  result = -1;
                } else {
                  result = 1;
                }
              } else if(at1->discrete_state < at2->discrete_state) {
                result = -1;
              } else {
                result = 1;
              }
            } else {
              result = wc;
            }
          } else {
            /* NOTE: don't forget to synchronize with below */
            if(SettingGetGlobal_b(G, cSetting_pdb_insertions_go_first)) {
              ov_size sl1, sl2;
              sl1 = strlen(at1->resi);
              sl2 = strlen(at2->resi);
              if(sl1 == sl2)
                result = wc;
              else if(sl1 < sl2)        /* sort residue 188A before 188, etc. */
                result = 1;
              else
                result = -1;
            } else if((at1->rank != at2->rank) &&
                      SettingGetGlobal_b(G, cSetting_rank_assisted_sorts)) {
              /* use rank to resolve insertion code ambiguities */
              if(at1->rank < at2->rank)
                result = -1;
              else
                result = 1;
            } else {
              result = wc;
            }
          }
        } else if(at1->resv < at2->resv) {
          result = -1;
        } else {
          result = 1;
        }
      } else if(at2->hetatm) {
        result = -1;
      } else {
        result = 1;
      }
    } else if((!at2->chain[0]) || (at1->chain[0] && ((at1->chain[0] < at2->chain[0])))) {
      result = -1;
    } else {
      result = 1;
    }
  } else {
    result = wc;
  }
  return (result);
}

int AtomInfoCompareIgnoreHet(PyMOLGlobals * G, AtomInfoType * at1, AtomInfoType * at2)
{
  int result;
  int wc;
  /* order by segment, chain, residue value, residue id, residue name, priority,
     and lastly by name */

  /*  printf(":%s:%s:%d:%s:%s:%i:%s =? ",
     at1->segi,at1->chain,at1->resv,at1->resi,at1->resn,at1->priority,at1->name);
     printf(":%s:%s:%d:%s:%s:%i:%s\n",
     at2->segi,at2->chain,at2->resv,at2->resi,at2->resn,at2->priority,at2->name);
   */

  wc = WordCompare(G, at1->segi, at2->segi, true);
  if(!wc) {
    if(at1->chain[0] == at2->chain[0]) {
      if(at1->resv == at2->resv) {
        wc = WordCompare(G, at1->resi, at2->resi, true);
        if(!wc) {
          wc = WordCompare(G, at1->resn, at2->resn, true);
          if(!wc) {
            if(at1->discrete_state == at2->discrete_state) {
              if(at1->priority == at2->priority) {
                if(at1->alt[0] == at2->alt[0]) {
                  if((result = AtomInfoNameCompare(G, at1->name, at2->name)) == 0) {
                    if(at1->rank < at2->rank)
                      result = -1;
                    else if(at1->rank > at2->rank)
                      result = 1;
                    else
                      result = 0;
                  }
                } else if((!at2->alt[0])
                          || (at1->alt[0] && ((at1->alt[0] < at2->alt[0])))) {
                  result = -1;
                } else {
                  result = 1;
                }
              } else if(at1->priority < at2->priority) {
                result = -1;
              } else {
                result = 1;
              }
            } else if(at1->discrete_state < at2->discrete_state) {
              result = -1;
            } else {
              result = 1;
            }
          } else {
            result = wc;
          }
        } else {
          /* NOTE: don't forget to synchronize with above */

          if(SettingGetGlobal_b(G, cSetting_pdb_insertions_go_first)) {
            ov_size sl1, sl2;
            sl1 = strlen(at1->resi);
            sl2 = strlen(at2->resi);
            if(sl1 == sl2)
              result = wc;
            else if(sl1 < sl2)  /* sort residue 188A before 188, etc. */
              result = 1;
            else
              result = -1;
          } else if((at1->rank != at2->rank) &&
                    SettingGetGlobal_b(G, cSetting_rank_assisted_sorts)) {
            /* use rank to resolve insertion code ambiguities */
            if(at1->rank < at2->rank)
              result = -1;
            else
              result = 1;
          } else {
            result = wc;
          }
        }
      } else if(at1->resv < at2->resv) {
        result = -1;
      } else {
        result = 1;
      }
    } else if((!at2->chain[0]) || (at1->chain[0] && ((at1->chain[0] < at2->chain[0])))) {
      result = -1;
    } else {
      result = 1;
    }
  } else {
    result = wc;
  }
  return (result);
}

int AtomInfoNameOrder(PyMOLGlobals * G, AtomInfoType * at1, AtomInfoType * at2)
{
  int result;
  if(at1->alt[0] == at2->alt[0]) {
    if(at1->priority == at2->priority) {
      result = AtomInfoNameCompare(G, at1->name, at2->name);
    } else if(at1->priority < at2->priority) {
      result = -1;
    } else {
      result = 1;
    }
  } else if((!at2->alt[0]) || (at1->alt[0] && ((at1->alt[0] < at2->alt[0])))) {
    result = -1;
  } else {
    result = 1;
  }
  return (result);
}

int AtomInfoSameResidue(PyMOLGlobals * G, AtomInfoType * at1, AtomInfoType * at2)
{
  if(at1->hetatm == at2->hetatm)
    if(at1->chain[0] == at2->chain[0])
      if(at1->resv == at2->resv)

        if(at1->discrete_state == at2->discrete_state)
          if(WordMatch(G, at1->resi, at2->resi, true) < 0)
            if(WordMatch(G, at1->segi, at2->segi, true) < 0)
              if(WordMatch(G, at1->resn, at2->resn, true) < 0)
                return 1;
  return 0;
}

int AtomInfoSameResidueP(PyMOLGlobals * G, AtomInfoType * at1, AtomInfoType * at2)
{
  if(at1 && at2)
    if(at1->hetatm == at2->hetatm)
      if(at1->chain[0] == at2->chain[0])
        if(at1->resv == at2->resv)
          if(at1->discrete_state == at2->discrete_state)
            if(WordMatch(G, at1->resi, at2->resi, true) < 0)
              if(WordMatch(G, at1->segi, at2->segi, true) < 0)
                if(WordMatch(G, at1->resn, at2->resn, true) < 0)
                  return 1;
  return 0;
}

int AtomInfoSameChainP(PyMOLGlobals * G, AtomInfoType * at1, AtomInfoType * at2)
{
  if(at1 && at2)
    if(at1->chain[0] == at2->chain[0])
      if(WordMatch(G, at1->segi, at2->segi, true) < 0)
        return 1;
  return 0;
}

int AtomInfoSameSegmentP(PyMOLGlobals * G, AtomInfoType * at1, AtomInfoType * at2)
{
  if(at1 && at2)
    if(WordMatch(G, at1->segi, at2->segi, true) < 0)
      return 1;
  return 0;
}

int AtomInfoSequential(PyMOLGlobals * G, AtomInfoType * at1, AtomInfoType * at2, int mode)
{
  char last1 = 0, last2 = 0;
  char *p;
  if(mode > 0) {
    if(at1->hetatm == at2->hetatm) {
      if(mode > 1) {
        if(WordMatch(G, at1->segi, at2->segi, true) < 0) {
          if(mode > 2) {
            if(at1->chain[0] == at2->chain[0]) {
              if(mode > 3) {
                if(at1->resv == at2->resv) {
                  if(mode > 4) {
                    p = at1->resi;
                    while(*p) {
                      last1 = (*p++);
                    }
                    p = at2->resi;
                    while(*p) {
                      last2 = (*p++);
                    }
                    if(last1 == last2)
                      return 1;
                    if((last1 + 1) == last2)
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

int AtomInfoMatch(PyMOLGlobals * G, AtomInfoType * at1, AtomInfoType * at2)
{
  if((tolower(at1->chain[0])) == (tolower(at2->chain[0])))
    if(WordMatch(G, at1->name, at2->name, true) < 0)
      if(WordMatch(G, at1->resi, at2->resi, true) < 0)
        if(WordMatch(G, at1->resn, at2->resn, true) < 0)
          if(WordMatch(G, at1->segi, at2->segi, true) < 0)
            if((tolower(at1->alt[0])) == (tolower(at2->alt[0])))
              return 1;
  return 0;
}

int AtomInfoIsFreeCation(PyMOLGlobals * G, AtomInfoType * I)
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

int AtomInfoGetExpectedValence(PyMOLGlobals * G, AtomInfoType * I)
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

static void set_protons(AtomInfoType * I, char *elem)
{
  char *e = I->elem;

  if(elem)
    e = elem;
  while((*e >= '0') && (*e <= '9') && (*(e + 1)))
    e++;

  if(!e[1]) {                   /* single letter */
    switch (e[0]) {
    case 'H':
      I->protons = cAN_H;
      break;
    case 'D':
      I->protons = cAN_H;
      break;
    case 'Q':
      I->protons = cAN_H;
      break;                    /* for NMR structures */
    case 'C':
      I->protons = cAN_C;
      break;
    case 'N':
      I->protons = cAN_N;
      break;
    case 'O':
      I->protons = cAN_O;
      break;
    case 'F':
      I->protons = cAN_F;
      break;
    case 'S':
      I->protons = cAN_S;
      break;
    case 'P':
      I->protons = cAN_P;
      break;
    case 'K':
      I->protons = cAN_K;
      break;
    case 'I':
      I->protons = cAN_I;
      break;
    case 'U':
      I->protons = cAN_U;
      break;
    default:                   /* unrecognized element (possible garbage in this column?) */
      if(!elem)
        set_protons(I, I->name);
      break;
    }
  } else {
    switch (e[0]) {
    case 'A':
      switch (e[1]) {
      case 'G':
      case 'g':
        I->protons = cAN_Ag;
        break;
      case 'L':
      case 'l':
        I->protons = cAN_Al;
        break;
      case 'R':
      case 'r':
        I->protons = cAN_Ar;
        break;
      case 'S':
      case 's':
        I->protons = cAN_As;
        break;
      case 'U':
      case 'u':
        I->protons = cAN_Au;
        break;
      }
      break;
    case 'B':
      switch (e[1]) {
      case 'A':
      case 'a':
        I->protons = cAN_Ba;
        break;
      case 'E':
      case 'e':
        I->protons = cAN_Be;
        break;
      case 'R':
      case 'r':
        I->protons = cAN_Br;
        break;
      }
      break;
    case 'C':
      switch (e[1]) {
      case 'A':
      case 'a':
        I->protons = cAN_Ca;
        break;
      case 'D':
      case 'd':
        I->protons = cAN_Cd;
        break;
      case 'E':
      case 'e':
        I->protons = cAN_Ce;
        break;
      case 'L':
      case 'l':
        I->protons = cAN_Cl;
        break;
      case 'O':
      case 'o':
        I->protons = cAN_Co;
        break;
      case 'R':
      case 'r':
        I->protons = cAN_Cr;
        break;
      case 'S':
      case 's':
        I->protons = cAN_Cs;
        break;
      case 'U':
      case 'u':
        I->protons = cAN_Cu;
        break;
      default:
        I->protons = cAN_C;     /* failsafe for wacky carbons */
        break;
      }
      break;
    case 'F':
      switch (e[1]) {
      case 'E':
      case 'e':
        I->protons = cAN_Fe;
        break;
      }
      break;
    case 'G':
      switch (e[1]) {
      case 'A':
      case 'a':
        I->protons = cAN_Ga;
        break;
      case 'E':
      case 'e':
        I->protons = cAN_Ge;
        break;
      }
      break;
    case 'H':
      switch (e[1]) {
      case 'G':
      case 'g':
        I->protons = cAN_Hg;
        break;
      case 'E':
      case 'e':
        I->protons = cAN_He;
        break;
      default:
        I->protons = cAN_H;     /* if unrecognized but begins with H, then assume H */
        break;
      }
      break;
    case 'I':
      switch (e[1]) {
      case 'N':
      case 'n':
        I->protons = cAN_In;
        break;
      }
      break;
    case 'K':
      switch (e[1]) {
      case 'R':
      case 'r':
        I->protons = cAN_Kr;
        break;
        break;
      }
      break;
    case 'L':
      switch (e[1]) {
      case 'I':
      case 'i':
        I->protons = cAN_LP;
        break;
        break;
      case 'P':
      case 'p':
        I->protons = cAN_LP;
        break;
        break;
      }
      break;
    case 'M':
      switch (e[1]) {
      case 'G':
      case 'g':
        I->protons = cAN_Mg;
        break;
      case 'N':
      case 'n':
        I->protons = cAN_Mn;
        break;
      }
      break;
    case 'N':
      switch (e[1]) {
      case 'A':
      case 'a':
        I->protons = cAN_Na;
        break;
      case 'E':
      case 'e':
        I->protons = cAN_Ne;
        break;
      case 'I':
      case 'i':
        I->protons = cAN_Ni;
        break;
      }
      break;
    case 'P':
      switch (e[1]) {
      case 'B':
      case 'b':
        I->protons = cAN_Pb;
        break;
      case 'D':
      case 'd':
        I->protons = cAN_Pd;
        break;
      case 'T':
      case 't':
        I->protons = cAN_Pt;
        break;
      }
      break;
    case 'S':
      switch (e[1]) {
      case 'B':
      case 'b':
        I->protons = cAN_Sb;
        break;
      case 'E':
      case 'e':
        I->protons = cAN_Se;
        break;
      case 'I':
      case 'i':
        I->protons = cAN_Si;
        break;
      case 'N':
      case 'n':
        I->protons = cAN_Sn;
        break;
      case 'R':
      case 'r':
        I->protons = cAN_Sr;
        break;
      }
      break;
    case 'T':
      switch (e[1]) {
      case 'I':
      case 'i':
        I->protons = cAN_Ti;
        break;
      case 'E':
      case 'e':
        I->protons = cAN_Te;
        break;
      case 'L':
      case 'l':
        I->protons = cAN_Tl;
        break;
      }
      break;

    case 'Z':
      switch (e[1]) {
      case 'N':
      case 'n':
        I->protons = cAN_Zn;
        break;
      }
      break;
    default:                   /* unrecognized element (possible garbage?) */
      if(!elem)
        set_protons(I, I->name);
      break;
    }
  }
}

void AtomInfoAssignParameters(PyMOLGlobals * G, AtomInfoType * I)
{
  char *n = NULL, *e = NULL;
  int pri;
  float vdw;

  e = I->elem;
  if(!*e) {                     /* try to guess atomic type from name */
    n = I->name;
    while(((*n) >= '0') && ((*n) <= '9') && (*(n + 1)))
      n++;
    while((((*n) >= 'A') && ((*n) <= 'Z')) || (((*n) >= 'a') && ((*n) <= 'z'))) {
      *(e++) = *(n++);
    }
    *e = 0;
    e = I->elem;
    switch (*e) {
    case 'C':
      if(*(e + 1) == 'A') {
        if(!(WordMatch(G, "CA", I->resn, true) < 0)
           && (!(WordMatch(G, "CA+", I->resn, true) < 0)))
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
    default:
      break;
    }
    if(*(e + 1))
      *(e + 1) = tolower(*(e + 1));
  }
  I->hydrogen = ((((*I->elem) == 'H') || ((*I->elem) == 'D') || ((*I->elem) == 'Q'))
                 && (!(*(I->elem + 1))));
  n = I->name;
  while((*n >= '0') && (*n <= '9') && (*(n + 1)))
    n++;
  if(toupper(*n) != I->elem[0]) {
    pri = 1000;                 /* unconventional atom name -- make no assignments */
  } else if((int) SettingGet(G, cSetting_pdb_standard_order)) {
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

  e = I->elem;
  set_protons(I, NULL);

  /* vdw radii */

  switch (I->protons) {

  case cAN_H:
    vdw = 1.20F;
    break;
    /* NOTE: Rowland and Taylor J. Phys. Chem. 100 (18): 7384–91 suggest
       1.09 would be more consistent with Hydrogen crystal data */

  case cAN_He:
    vdw = 1.40F;
    break;

  case cAN_Li:
    vdw = 1.82F;
    break;

  case cAN_B:
    vdw = 1.85F;
    break;
  case cAN_C:
    vdw = 1.70F;
    break;
  case cAN_N:
    vdw = 1.55F;
    break;
  case cAN_O:
    vdw = 1.52F;
    break;
  case cAN_F:
    vdw = 1.47F;
    break;
  case cAN_Ne:
    vdw = 1.54F;
    break;

  case cAN_Na:
    vdw = 2.27F;
    break;
  case cAN_Mg:
    vdw = 1.73F;
    break;
  case cAN_Al:
    vdw = 2.00F;
    break;
  case cAN_Si:
    vdw = 2.10F;
    break;
  case cAN_P:
    vdw = 1.80F;
    break;
  case cAN_S:
    vdw = 1.80F;
    break;
  case cAN_Cl:
    vdw = 1.75F;
    break;
  case cAN_Ar:
    vdw = 1.88F;
    break;

  case cAN_K:
    vdw = 2.75F;
    break;

  case cAN_Mn:
    vdw = 1.73F;
    break;
  case cAN_Ni:
    vdw = 1.63F;
    break;
  case cAN_Cu:
    vdw = 1.40F;
    break;
  case cAN_Zn:
    vdw = 1.39F;
    break;
  case cAN_Ga:
    vdw = 1.87F;
    break;

  case cAN_As:
    vdw = 1.85F;
    break;
  case cAN_Se:
    vdw = 1.90F;
    break;
  case cAN_Br:
    vdw = 1.85F;
    break;
  case cAN_Kr:
    vdw = 2.02F;
    break;

  case cAN_Pd:
    vdw = 1.63F;
    break;
  case cAN_Ag:
    vdw = 1.72F;
    break;
  case cAN_Cd:
    vdw = 1.58F;
    break;
  case cAN_In:
    vdw = 1.93F;
    break;
  case cAN_Sn:
    vdw = 2.17F;
    break;

  case cAN_Te:
    vdw = 2.06F;
    break;
  case cAN_I:
    vdw = 1.98F;
    break;
  case cAN_Xe:
    vdw = 2.16F;
    break;

  case cAN_Pt:
    vdw = 1.75F;
    break;
  case cAN_Au:
    vdw = 1.66F;
    break;
  case cAN_Hg:
    vdw = 1.55F;
    break;
  case cAN_Tl:
    vdw = 1.96F;
    break;
  case cAN_Pb:
    vdw = 2.02F;
    break;

  case cAN_U:
    vdw = 1.86F;
    break;

  case cAN_LP:
    vdw = 0.5F;
    break;                      /* lone pairs @ 0.5 same as MOE? */

  default:
    vdw = 1.80F;
    break;                      /* default radius for known atoms with unknown radii */
  }

  if(SettingGet(G, cSetting_legacy_vdw_radii)) {        /* ver<0.75, old, incorrect VDW */
    if(!strcmp(e, "N"))
      vdw = 1.8F;               /* slow but compact */
    if(!strcmp(e, "C"))
      vdw = 1.8F;
    if(!strcmp(e, "Cl"))
      vdw = 1.8F;
    if(!strcmp(e, "O"))
      vdw = 1.5F;
    if(!strcmp(e, "Br"))
      vdw = 1.9F;
    if(!strcmp(e, "I"))
      vdw = 2.15F;
    if(!strcmp(e, "S"))
      vdw = 1.9F;
    if(!strcmp(e, "P"))
      vdw = 1.9F;
    if(!strcmp(e, "F"))
      vdw = 1.35F;
    if(!strcmp(e, "H"))
      vdw = 1.1F;
  }

  if(I->vdw == 0.0)             /* only assigned if not yet assigned */
    I->vdw = vdw;

#if 0
  if(!I->protons)
    I->protons = cAN_C;         /* default assumption */
#endif

  if(I->protons == cAN_H)
    I->hydrogen = true;
  /*  printf("I->name %s I->priority %d\n",I->name,I->priority); */
}
