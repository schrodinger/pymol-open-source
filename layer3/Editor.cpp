
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
#include"os_gl.h"
#include"Err.h"

#include"PConv.h"
#include"MemoryDebug.h"
#include"Vector.h"
#include"Matrix.h"
#include"ButMode.h"
#include"Scene.h"
#include"Editor.h"
#include"Selector.h"
#include"Ortho.h"
#include"main.h"
#include"Color.h"
#include"Setting.h"
#include"Util.h"
#include"Executive.h"
#include"P.h"
#include"CGO.h"
#include "Lex.h"

struct _CEditor {
  ObjectMolecule *DihedObject;
  WordType DragSeleName;
  int Active;
  int ActiveState;
  int DragIndex;
  int DragSelection;
  int DragHaveAxis, DragHaveBase, DragBondFlag, DragSlowFlag;
  int PickMode;                 /* 1 = atom, 2 = bond, 3 = multiatom */
  int NextPickSele;
  int BondMode;
  pymol::CObject *DragObject;
  int NFrag;
  float V0[3], V1[3], Axis[3], Center[3], DragBase[3];
  float *PosVLA;
  int ShowFrags;
  int DihedralInvalid;
  int MouseInvalid;
  int FavorOrigin;
  float FavoredOrigin[3];
  CGO *shaderCGO;
};

int EditorGetScheme(PyMOLGlobals * G)
{
  CEditor *I = G->Editor;
  int scheme = EDITOR_SCHEME_OBJ;

  if(EditorActive(G))
    scheme = EDITOR_SCHEME_FRAG;
  else if(I->DragObject) {
    if(I->DragIndex >= 0) {
      scheme = EDITOR_SCHEME_OBJ;
    } else {
      scheme = EDITOR_SCHEME_DRAG;
    }
  }
  return scheme;
}

void EditorFavorOrigin(PyMOLGlobals * G, float *v1)
{
  CEditor *I = G->Editor;
  if(v1) {
    I->FavorOrigin = true;
    copy3f(v1, I->FavoredOrigin);
  } else {
    I->FavorOrigin = false;
  }
}

static void EditorDrawDihedral(PyMOLGlobals * G)
{
  if(EditorActive(G) && EditorIsBondMode(G)
     && SettingGetGlobal_b(G, cSetting_editor_auto_dihedral)) {
    ObjectMolecule *obj1, *obj2;
    int at1, at2, at0, at3;
    int sele1 = SelectorIndexByName(G, cEditorSele1);
    int sele2 = SelectorIndexByName(G, cEditorSele2);
    if((sele1 >= 0) && (sele2 >= 0)) {
      obj1 = SelectorGetFastSingleAtomObjectIndex(G, sele1, &at1);
      obj2 = SelectorGetFastSingleAtomObjectIndex(G, sele2, &at2);
      if(obj1 && (obj1 == obj2)) {
        CEditor *I = G->Editor;

        I->DihedObject = obj1;
        at0 = ObjectMoleculeGetTopNeighbor(G, obj1, at1, at2);
        at3 = ObjectMoleculeGetTopNeighbor(G, obj1, at2, at1);

        if((at0 >= 0) && (at3 >= 0)) {
          /* find the highest priority atom attached to index1 */
          SelectorCreateOrderedFromObjectIndices(G, cEditorDihe1, obj1, &at0, 1);
          SelectorCreateOrderedFromObjectIndices(G, cEditorDihe2, obj2, &at3, 1);

          ExecutiveDihedral(G, cEditorDihedral, cEditorDihe1, cEditorSele1,
              cEditorSele2, cEditorDihe2, 0, true, true, false, true, -1);
          ExecutiveColor(G, cEditorDihedral, "white", 1, true);
          ExecutiveSetSettingFromString(G, cSetting_float_labels,
                                        "1", cEditorDihedral, 0, true, true);
#ifndef _PYMOL_FREETYPE
          ExecutiveSetSettingFromString(G, cSetting_label_font_id,
                                        "4", cEditorDihedral, 0, true, true);
#else
          ExecutiveSetSettingFromString(G, cSetting_label_font_id,
                                        "8", cEditorDihedral, 0, true, true);
          ExecutiveSetSettingFromString(G, cSetting_label_size,
                                        "20", cEditorDihedral, 0, true, true);
#endif
          ExecutiveSetSettingFromString(G, cSetting_label_color,
                                        "brightorange", cEditorDihedral, 0, true, true);
        }
      }
    }
  }
}

void EditorDihedralInvalid(PyMOLGlobals * G, ObjectMolecule * obj)
{
  CEditor *I = G->Editor;
  if(!obj)
    I->DihedralInvalid = true;
  else if(obj == I->DihedObject)
    I->DihedralInvalid = true;
}

void EditorMouseInvalid(PyMOLGlobals * G)
{
  CEditor *I = G->Editor;
  I->MouseInvalid = true;
}

static void EditorConfigMouse(PyMOLGlobals * G)
{

  int scheme = EditorGetScheme(G);
  const char *mouse_mode = SettingGetGlobal_s(G, cSetting_button_mode_name);

  if(mouse_mode && (!strcmp(mouse_mode, "3-Button Editing") ||
                    !strcmp(mouse_mode, "3-Button Motions"))) {
    /* WEAK! */
    int button;

    button = cButModeMiddleShft;

    {
      int action = ButModeGet(G, button);
      if((action == cButModeMovFrag) ||
         (action == cButModeMovObj) || (action == cButModeMovDrag)) {
        switch (scheme) {
        case EDITOR_SCHEME_OBJ:
          action = cButModeMovObj;
          break;
        case EDITOR_SCHEME_FRAG:
          action = cButModeMovFrag;
          break;
        case EDITOR_SCHEME_DRAG:
          action = cButModeMovDrag;
          break;
        }
        ButModeSet(G, button, action);
      }
    }

    button = cButModeLeftShft;
    {
      int action = ButModeGet(G, button);
      if((action == cButModeRotFrag) ||
         (action == cButModeRotObj) || (action == cButModeRotDrag)) {
        switch (scheme) {
        case EDITOR_SCHEME_OBJ:
          action = cButModeRotObj;
          break;
        case EDITOR_SCHEME_FRAG:
          action = cButModeRotFrag;
          break;
        case EDITOR_SCHEME_DRAG:
          action = cButModeRotDrag;
          break;
        }
        ButModeSet(G, button, action);
      }
    }

    button = cButModeRightShft;
    {
      int action = ButModeGet(G, button);
      if((action == cButModeMovFragZ) ||
         (action == cButModeMovObjZ) || (action == cButModeMovDragZ)) {
        switch (scheme) {
        case EDITOR_SCHEME_OBJ:
          action = cButModeMovObjZ;
          break;
        case EDITOR_SCHEME_FRAG:
          action = cButModeMovFragZ;
          break;
        case EDITOR_SCHEME_DRAG:
          action = cButModeMovDragZ;
          break;
        }
        ButModeSet(G, button, action);
      }
    }

    button = cButModeLeftCtrl;
    {
      int action = ButModeGet(G, button);
      if((action == cButModeMoveAtom) || (action == cButModeTorFrag)) {
        switch (scheme) {
        case EDITOR_SCHEME_OBJ:
          action = cButModeMoveAtom;
          break;
        case EDITOR_SCHEME_FRAG:
          action = cButModeTorFrag;
          break;
        case EDITOR_SCHEME_DRAG:
          action = cButModeMoveAtom;
          break;
        }
        ButModeSet(G, button, action);
      }
    }

    button = cButModeLeftDouble;
    {
      int action = ButModeGet(G, button);
      if((action == cButModeMoveAtom) || (action == cButModeTorFrag)) {
        switch (scheme) {
        case EDITOR_SCHEME_OBJ:
          action = cButModeMoveAtom;
          break;
        case EDITOR_SCHEME_FRAG:
          action = cButModeTorFrag;
          break;
        case EDITOR_SCHEME_DRAG:
          action = cButModeMoveAtom;
          break;
        }
        ButModeSet(G, button, action);
      }
    }

    button = cButModeLeftCtSh;
    {
      int action = ButModeGet(G, button);
      if((action == cButModeMoveAtom) || (action == cButModeMoveAtomZ)) {
        switch (scheme) {
        case EDITOR_SCHEME_OBJ:
          action = cButModeMoveAtomZ;
          break;
        case EDITOR_SCHEME_FRAG:
          action = cButModeMoveAtom;
          break;
        case EDITOR_SCHEME_DRAG:
          action = cButModeMoveAtomZ;
          break;
        }
        ButModeSet(G, button, action);
      }
    }
  }

}

void EditorUpdate(PyMOLGlobals * G)
{
  CEditor *I = G->Editor;
  if(I->DihedralInvalid) {
    EditorDrawDihedral(G);
    I->DihedralInvalid = false;
  }
  if(I->MouseInvalid) {
    EditorConfigMouse(G);
    I->MouseInvalid = false;
  }
}

static int EditorGetEffectiveState(PyMOLGlobals * G, pymol::CObject * obj, int state)
{
  if(obj && (obj->type == cObjectMolecule)) {
    ObjectMolecule *objMol = (ObjectMolecule*)(void*)obj;
    if(!objMol)
      objMol = SelectorGetFastSingleObjectMolecule(G, SelectorIndexByName(G, cEditorSele1));
    if(!objMol)
      objMol = SelectorGetFastSingleObjectMolecule(G, SelectorIndexByName(G, cEditorSele2));
    if(!objMol)
      objMol = SelectorGetFastSingleObjectMolecule(G, SelectorIndexByName(G, cEditorSele3));
    if(!objMol)
      objMol = SelectorGetFastSingleObjectMolecule(G, SelectorIndexByName(G, cEditorSele4));
    
    if(objMol) {
      if((objMol->NCSet == 1) && (state > 0))
        if(SettingGet_i(G, NULL, objMol->Setting.get(), cSetting_static_singletons))
          return 0;
    }
  }
  return state;
}

int EditorGetNFrag(PyMOLGlobals * G)
{
  CEditor *I = G->Editor;
  if(EditorActive(G)) {
    return I->NFrag;
  }
  return 0;
}

void EditorDefineExtraPks(PyMOLGlobals * G)
{
  WordType name;
  WordType buffer;

  if(EditorGetSinglePicked(G, name)) {
    sprintf(buffer, "(byres %s)", name);
    SelectorCreate(G, cEditorRes, buffer, NULL, true, NULL);
    sprintf(buffer, "(bychain %s)", name);
    SelectorCreate(G, cEditorChain, buffer, NULL, true, NULL);
    sprintf(buffer, "(byobject %s)", name);
    SelectorCreate(G, cEditorObject, buffer, NULL, true, NULL);

    if(SettingGetGlobal_b(G, cSetting_auto_hide_selections))
      ExecutiveHideSelections(G);
    EditorInvalidateShaderCGO(G);
  }
}

int EditorDeselectIfSelected(PyMOLGlobals * G, ObjectMolecule * obj, int index,
                             int update)
{
  CEditor *I = G->Editor;
  int result = false;
  int s, sele;
  if(obj) {
    if((index >= 0) && (index < obj->NAtom)) {
      s = obj->AtomInfo[index].selEntry;
      sele = SelectorIndexByName(G, cEditorSele1);
      if(SelectorIsMember(G, s, sele)) {
        ExecutiveDelete(G, cEditorSele1);
        result = true;
      }
      sele = SelectorIndexByName(G, cEditorSele2);
      if(SelectorIsMember(G, s, sele)) {
        ExecutiveDelete(G, cEditorSele2);
        result = true;
      }
      sele = SelectorIndexByName(G, cEditorSele3);
      if(SelectorIsMember(G, s, sele)) {
        ExecutiveDelete(G, cEditorSele3);
        result = true;
      }
      sele = SelectorIndexByName(G, cEditorSele4);
      if(SelectorIsMember(G, s, sele)) {
        ExecutiveDelete(G, cEditorSele4);
        result = true;
      }
      if(result && update)
        EditorActivate(G, I->ActiveState, I->BondMode);
    }
  }

  return result;
}

int EditorIsBondMode(PyMOLGlobals * G)
{
  CEditor *I = G->Editor;
  return (I->BondMode);
}

PyObject *EditorAsPyList(PyMOLGlobals * G)
{
  PyObject *result = NULL;
  CEditor *I = G->Editor;

  if(!EditorActive(G)) {
    result = PyList_New(0);     /* not editing? return null list */
  } else {
    result = PyList_New(3);
    PyList_SetItem(result, 0, PyString_FromString(""));
    PyList_SetItem(result, 1, PyInt_FromLong(I->ActiveState));
    PyList_SetItem(result, 2, PyInt_FromLong(I->BondMode));
  }
  return (PConvAutoNone(result));
}

int EditorFromPyList(PyMOLGlobals * G, PyObject * list)
{
  int ok = true;
  int active_flag = false;
  int active_state;
  WordType obj_name;
  int ll = 0;
  int bond_mode = true;

  if(ok)
    ok = (list != NULL);
  if(ok)
    ok = PyList_Check(list);
  if(ok)
    ll = PyList_Size(list);
  /* TO SUPPORT BACKWARDS COMPATIBILITY...
     Always check ll when adding new PyList_GetItem's */
  if(ok)
    active_flag = (PyList_Size(list) != 0);
  if(!active_flag) {
    EditorInactivate(G);
  } else {
    if(ok)
      ok = PConvPyStrToStr(PyList_GetItem(list, 0), obj_name, sizeof(WordType));
    if(ok)
      ok = PConvPyIntToInt(PyList_GetItem(list, 1), &active_state);
    if(ok && (ll > 2))
      ok = PConvPyIntToInt(PyList_GetItem(list, 2), &bond_mode);        /* newer session files */
    if(ok) {
      EditorActivate(G, active_state, bond_mode);
      EditorDefineExtraPks(G);
    } else {
      EditorInactivate(G);
    }
  }
  if(!ok) {
    EditorInactivate(G);
  }
  return (ok);
}

int EditorActive(PyMOLGlobals * G)
{
  CEditor *I = G->Editor;
  return (I->Active);
}

pymol::CObject *EditorDragObject(PyMOLGlobals * G)
{
  CEditor *I = G->Editor;
  return I->DragObject;
}

int EditorGetSinglePicked(PyMOLGlobals * G, char *name)
{
  int cnt = 0;
  int sele;
  if((sele = SelectorIndexByName(G, cEditorSele1)) >= 0) {
    cnt++;
    if(name)
      strcpy(name, cEditorSele1);
  }
  if((sele = SelectorIndexByName(G, cEditorSele2)) >= 0) {
    cnt++;
    if(name)
      strcpy(name, cEditorSele2);
  }
  if((sele = SelectorIndexByName(G, cEditorSele3)) >= 0) {
    cnt++;
    if(name)
      strcpy(name, cEditorSele3);
  }
  if((sele = SelectorIndexByName(G, cEditorSele4)) >= 0) {
    cnt++;
    if(name)
      strcpy(name, cEditorSele4);
  }
  return (cnt == 1);
}

void EditorGetNextMultiatom(PyMOLGlobals * G, char *name)
{
  CEditor *I = G->Editor;
  int sele;
  sele = SelectorIndexByName(G, cEditorSele1);
  if(sele < 0) {
    strcpy(name, cEditorSele1);
    I->NextPickSele = 0;
    return;
  }
  sele = SelectorIndexByName(G, cEditorSele2);
  if(sele < 0) {
    strcpy(name, cEditorSele2);
    I->NextPickSele = 1;
    return;
  }
  sele = SelectorIndexByName(G, cEditorSele3);
  if(sele < 0) {
    strcpy(name, cEditorSele3);
    I->NextPickSele = 2;
    return;
  }
  sele = SelectorIndexByName(G, cEditorSele4);
  if(sele < 0) {
    strcpy(name, cEditorSele4);
    I->NextPickSele = 3;
    return;
  }
  strcpy(name, cEditorSele4);
  I->NextPickSele = 3;
  return;
  /*
     I->NextPickSele = (++I->NextPickSele)&0x3;
     switch(I->NextPickSele) {
     case 0: strcpy(name,cEditorSele1); break;
     case 1: strcpy(name,cEditorSele2); break;
     case 2: strcpy(name,cEditorSele3); break;
     case 3: strcpy(name,cEditorSele4); break;
     }
     return;    
   */
}


/*========================================================================*/
int EditorLogState(PyMOLGlobals * G, int pkresi)
{
  CEditor *I = G->Editor;
  if(SettingGetGlobal_i(G, cSetting_logging)) {

    OrthoLineType buffer, buf1 = "None", buf2 = "None", buf3 = "None", buf4 = "None";
    int pkbond = 1;

    if(!EditorActive(G)) {
      PLog(G, "edit", cPLog_pml);
    } else {
      int sele1, sele2, sele3, sele4;
      ObjectMolecule *obj1 = NULL, *obj2 = NULL, *obj3 = NULL, *obj4 = NULL;
      int index1, index2, index3, index4;

      sele1 = SelectorIndexByName(G, cEditorSele1);
      sele2 = SelectorIndexByName(G, cEditorSele2);
      sele3 = SelectorIndexByName(G, cEditorSele3);
      sele4 = SelectorIndexByName(G, cEditorSele4);

      obj1 = SelectorGetFastSingleAtomObjectIndex(G, sele1, &index1);
      obj2 = SelectorGetFastSingleAtomObjectIndex(G, sele2, &index2);
      obj3 = SelectorGetFastSingleAtomObjectIndex(G, sele3, &index3);
      obj4 = SelectorGetFastSingleAtomObjectIndex(G, sele4, &index4);

      if((sele1 >= 0) && (sele2 >= 0) && I->BondMode && obj1 && obj2) {

        /* bond mode */
        ObjectMoleculeGetAtomSeleLog(obj1, index1, buf1, true);
        ObjectMoleculeGetAtomSeleLog(obj2, index2, buf2, true);

      } else {

        /* atom mode */
        pkbond = 0;

        if(obj1) {
          ObjectMoleculeGetAtomSeleLog(obj1, index1, buf1, true);
        }

        if(obj2) {
          ObjectMoleculeGetAtomSeleLog(obj2, index2, buf2, true);
        }

        if(obj3) {
          ObjectMoleculeGetAtomSeleLog(obj3, index3, buf3, true);
        }

        if(obj4) {
          ObjectMoleculeGetAtomSeleLog(obj4, index4, buf4, true);
        }
      }

      sprintf(buffer, "cmd.edit(%s,%s,%s,%s,pkresi=%d,pkbond=%d)",
              buf1, buf2, buf3, buf4, pkresi ? 1 : 0, pkbond ? 1 : 0);

      PLog(G, buffer, cPLog_pym);

    }
  }
  return 1;
}


/*========================================================================*/

pymol::Result<>
EditorInvert(PyMOLGlobals * G, int quiet)
{
  CEditor *I = G->Editor;
  int sele0, sele1, sele2;
  int i0, frg;
  int ia0 = -1;
  int ia1 = -1;
  float v[3], v0[3], v1[3];
  float n0[3], n1[3];
  float m[16];
  int state;
  int vf, vf0, vf1;
  int found = false;
  WordType name;
  ObjectMolecule *obj0, *obj1, *obj2;

  if(!EditorActive(G)) {
    return pymol::Error("Must pick an atom to invert");
  } else {
    sele0 = SelectorIndexByName(G, cEditorSele1);
    sele1 = SelectorIndexByName(G, cEditorSele2);
    sele2 = SelectorIndexByName(G, cEditorSele3);
    obj0 = SelectorGetFastSingleAtomObjectIndex(G, sele0, &i0);
    obj1 = SelectorGetFastSingleAtomObjectIndex(G, sele1, &ia0);
    obj2 = SelectorGetFastSingleAtomObjectIndex(G, sele2, &ia1);
    if(sele0 < 0) {
      return pymol::Error("Must pick atom to invert as pk1");
    } else if(sele1 < 0) {
      return pymol::Error("Must pick immobile atom in pk2");
    } else if(sele2 < 0) {
      return pymol::Error("Must pick immobile atom in pk3");
    } else if(!(obj0 && (obj0 == obj1) && (obj0 = obj2))) {
      return pymol::Error("Must pick three atoms in the same object");
    } else {
      state = SceneGetState(G);
      ObjectMoleculeSaveUndo(obj0, state, false);

      vf = ObjectMoleculeGetAtomVertex(obj0, state, i0, v);
      vf0 = ObjectMoleculeGetAtomVertex(obj0, state, ia0, v0);
      vf1 = ObjectMoleculeGetAtomVertex(obj0, state, ia1, v1);

      if(vf & vf0 & vf1) {
        subtract3f(v, v0, n0);
        subtract3f(v, v1, n1);
        normalize3f(n0);
        normalize3f(n1);

        add3f(n0, n1, n0);
        normalize3f(n0);

        get_rotation_about3f3fTTTf((float) cPI, n0, v, m);
        for(frg = 1; frg <= I->NFrag; frg++) {
          sprintf(name, "%s%1d", cEditorFragPref, frg);
          sele2 = SelectorIndexByName(G, name);

          if (ObjectMoleculeIsAtomBondedToSele(obj0, i0, sele2) &&
              (!ObjectMoleculeIsAtomBondedToSele(obj0, ia0, sele2)) &&
              (!ObjectMoleculeIsAtomBondedToSele(obj0, ia1, sele2))) {
            found = true;
              ObjectMoleculeTransformSelection(obj0, state, sele2, m, false, NULL, false,
                                               false);
          }
        }
        if(found) {
          if(!quiet) {
            PRINTFB(G, FB_Editor, FB_Actions)
              " Editor: Inverted atom.\n" ENDFB(G);
          }
        } else {
          return pymol::Error("No free fragments found for inversion");
        }

        SceneInvalidate(G);
        I->DragIndex = -1;
        I->DragSelection = -1;
        I->DragObject = NULL;
      }
    }
  }
  return {};
}


/*========================================================================*/
pymol::Result<> EditorTorsion(PyMOLGlobals * G, float angle)
{
  CEditor *I = G->Editor;
  int sele0, sele1, sele2;
  int i0, i1;
  float v0[3], v1[3];
  float d1[3], n0[3];
  float theta;
  float m[16];
  int state;
  int vf1, vf2;
  int ok = false;
  WordType sele;
  ObjectMolecule *obj0 = NULL, *obj1 = NULL, *obj2 = NULL;

  if(!EditorActive(G)) {
    ErrMessage(G, "Editor", "Must specify a bond first.");
  } else {
    sele0 = SelectorIndexByName(G, cEditorSele1);
    if(sele0 >= 0) {
      obj0 = SelectorGetFastSingleAtomObjectIndex(G, sele0, &i0);
      sele1 = SelectorIndexByName(G, cEditorSele2);
      obj1 = SelectorGetFastSingleAtomObjectIndex(G, sele1, &i1);
      strcpy(sele, cEditorFragPref);
      strcat(sele, "1");
      sele2 = SelectorIndexByName(G, sele);
      obj2 = SelectorGetFastSingleObjectMolecule(G, sele2);
      if(!((sele0 >= 0) && (sele1 >= 0) && (sele2 >= 0) && (obj0 == obj1))) {
        return pymol::Error("Must specify a bond first.");
      } else {
        if((i0 >= 0) && (i1 >= 0)) {
          state = SceneGetState(G);

          vf1 = ObjectMoleculeGetAtomVertex(obj0, state, i0, I->V0);
          vf2 = ObjectMoleculeGetAtomVertex(obj1, state, i1, I->V1);

          if(vf1 && vf2) {
            ObjectMoleculeSaveUndo(obj0, SceneGetState(G), false);

            subtract3f(I->V1, I->V0, I->Axis);
            average3f(I->V1, I->V0, I->Center);
            normalize3f(I->Axis);

            copy3f(I->V0, v1);
            copy3f(I->V1, v0);

            subtract3f(v1, v0, d1);
            copy3f(d1, n0);
            normalize3f(n0);

            theta = (float) (cPI * angle / 180.0);
            get_rotation_about3f3fTTTf(theta, n0, v1, m);
            ok =
              ObjectMoleculeTransformSelection(obj2, state, sele2, m, false, NULL, false,
                                               false);
            SceneInvalidate(G);

            I->DragIndex = -1;
            I->DragSelection = -1;
            I->DragObject = NULL;

            if(I->BondMode && SettingGetGlobal_b(G, cSetting_editor_auto_dihedral))
              EditorDihedralInvalid(G, NULL);
          }
        }
      }
    }
  }
  if(ok) {
    return {};
  } else {
    return pymol::Error("Error occurred.");
  }
}


/*========================================================================*/
pymol::Result<> EditorSelect(PyMOLGlobals* G, const char* str0,
    const char* str1, const char* str2, const char* str3, int pkresi,
    int pkbond, int quiet)
{
  SelectorTmp tmp0(G, str0);
  SelectorTmp tmp1(G, str1);
  SelectorTmp tmp2(G, str2);
  SelectorTmp tmp3(G, str3);
  auto s0 = tmp0.getName();
  auto s1 = tmp1.getName();
  auto s2 = tmp2.getName();
  auto s3 = tmp3.getName();
  int i0 = -1;
  int i1 = -1;
  int i2 = -1;
  int i3 = -1;
  int result = false;
  int ok = true;
  ObjectMolecule *obj0 = NULL, *obj1 = NULL, *obj2 = NULL, *obj3 = NULL;

  if(s0)
    if(!*s0)
      s0 = NULL;
  if(s1)
    if(!*s1)
      s1 = NULL;
  if(s2)
    if(!*s2)
      s2 = NULL;
  if(s3)
    if(!*s3)
      s3 = NULL;

  if(s0) {
    obj0 = SelectorGetFastSingleAtomObjectIndex(G, tmp0.getIndex(), &i0);
    ExecutiveDelete(G, cEditorSele1);
  }

  if(s1) {
    obj1 = SelectorGetFastSingleAtomObjectIndex(G, tmp1.getIndex(), &i1);
    ExecutiveDelete(G, cEditorSele2);
  }

  if(s2) {
    obj2 = SelectorGetFastSingleAtomObjectIndex(G, tmp2.getIndex(), &i2);
    ExecutiveDelete(G, cEditorSele3);
  }

  if(s3) {
    obj3 = SelectorGetFastSingleAtomObjectIndex(G, tmp3.getIndex(), &i3);
    ExecutiveDelete(G, cEditorSele4);
  }

  if(!(obj0 || obj1 || obj2 || obj3))
    ok = false;

  if(ok) {
    if(obj0)
      ObjectMoleculeVerifyChemistry(obj0, -1);
    if(obj1 && (obj1 != obj0))
      ObjectMoleculeVerifyChemistry(obj1, -1);
    if(obj2 && (obj2 != obj0) && (obj2 != obj1))
      ObjectMoleculeVerifyChemistry(obj2, -1);
    if(obj3 && (obj3 != obj0) && (obj3 != obj1) && (obj3 != obj2))
      ObjectMoleculeVerifyChemistry(obj3, -1);

    if(i0 >= 0)
      SelectorCreate(G, cEditorSele1, s0, NULL, quiet, NULL);
    if(i1 >= 0)
      SelectorCreate(G, cEditorSele2, s1, NULL, quiet, NULL);
    if(i2 >= 0)
      SelectorCreate(G, cEditorSele3, s2, NULL, quiet, NULL);
    if(i3 >= 0)
      SelectorCreate(G, cEditorSele4, s3, NULL, quiet, NULL);

    EditorActivate(G, SceneGetState(G), pkbond);

    if(pkresi)
      EditorDefineExtraPks(G);

    SceneInvalidate(G);
    result = true;

  } else {
    EditorInactivate(G);
    if(s0 && s0[0]) {
      return pymol::Error("Invalid input selection(s)");
    }
  }
  return {};
}


/*========================================================================*/
int EditorIsAnActiveObject(PyMOLGlobals * G, ObjectMolecule * obj)
{
  if(EditorActive(G)) {
    if(obj) {
      if(obj == SelectorGetFastSingleObjectMolecule(G,
                                                    SelectorIndexByName(G, cEditorSele1)))
        return true;
      if(obj == SelectorGetFastSingleObjectMolecule(G,
                                                    SelectorIndexByName(G, cEditorSele2)))
        return true;
      if(obj == SelectorGetFastSingleObjectMolecule(G,
                                                    SelectorIndexByName(G, cEditorSele3)))
        return true;
      if(obj == SelectorGetFastSingleObjectMolecule(G,
                                                    SelectorIndexByName(G, cEditorSele4)))
        return true;
    }
  }
  return false;
}


/*========================================================================*/
pymol::Result<> EditorCycleValence(PyMOLGlobals * G, int quiet)
{
  CEditor *I = G->Editor;

  if(EditorActive(G)) {
    for(const char* eSele : {cEditorSele3, cEditorSele4}) {
      if(SelectorIndexByName(G, eSele) >= 0) {
        return pymol::make_error("Only two picked selections allowed.");
      }
    }
    ObjectMolecule *obj0, *obj1;
    auto sele0 = SelectorIndexByName(G, cEditorSele1);
    if(sele0 >= 0) {
      auto sele1 = SelectorIndexByName(G, cEditorSele2);
      if(sele1 >= 0) {
        obj0 = SelectorGetFastSingleObjectMolecule(G, sele0);
        obj1 = SelectorGetFastSingleObjectMolecule(G, sele1);
        if(obj0 != obj1) {
          return pymol::make_error(
              "Both pk selections must belong to the same molecule.");
        }
        if(I->BondMode) {
          ObjectMoleculeVerifyChemistry(obj0, -1);
          ObjectMoleculeAdjustBonds(obj0, sele0, sele1, 0, 0);
        } else {
          return pymol::make_error("Invalid bond.");
        }
      } else {
        return pymol::make_error("No valid pk2 selection.");
      }
    } else {
      return pymol::make_error("No valid pk1 selection.");
    }
  }
  return {};
}


/*========================================================================*/
pymol::Result<> EditorAttach(PyMOLGlobals * G, const char *elem, int geom, int valence,
                  const char *name, int quiet)
{

  if(EditorActive(G)) {

    for(const char* eSele : {cEditorSele3, cEditorSele4}) {
      if(SelectorIndexByName(G, eSele) >= 0) {
        return pymol::make_error("Only 1 or 2 picked selections allowed.");
      }
    }

    auto sele0 = SelectorIndexByName(G, cEditorSele1);
    if(sele0 >= 0) {
      auto sele1 = SelectorIndexByName(G, cEditorSele2);
      auto obj0 = SelectorGetFastSingleObjectMolecule(G, sele0);
      auto obj1 = SelectorGetFastSingleObjectMolecule(G, sele1);
      if(obj0) {
        if(obj0->DiscreteFlag) {
          return pymol::make_error("Can't attach atoms onto discrete objects.");
        } else {
          ObjectMoleculeVerifyChemistry(obj0, -1);      /* remember chemistry for later */
          if(obj1) {
            if(obj0 == obj1) {
              /* bond mode - behave like replace */
              EditorReplace(G, elem, geom, valence, name, quiet);
            } else {
               return pymol::make_error("Picked atoms must belong to the same object.");
            }
          } else {
            pymol::vla<AtomInfoType> atInfo(1);
            auto ai = &atInfo[0];
            /* atom mode */
            auto i0 = ObjectMoleculeGetAtomIndex(obj0, sele0);       /* slow */
            if(i0 >= 0) {
              UtilNCopy(ai->elem, elem, sizeof(ElemName));
              ai->geom = geom;
              ai->valence = valence;
              if(name[0])
                LexAssign(G, ai->name, name);
              if(!ObjectMoleculeAttach(obj0, i0, std::move(atInfo))) {
                return pymol::make_error("Could not attach atom.");
              }
            }
          }
        }
      } else {
        return pymol::make_error("Invalid object.");
      }
    } else {
      return pymol::make_error("Invalid pk1 selection.");
    }
  }
  return {};
}


/*========================================================================*/
pymol::Result<> EditorRemove(PyMOLGlobals * G, int hydrogen, int quiet)
{
#define cEditorRemoveSele "_EditorRemove"

  if(!EditorActive(G)) {
    return pymol::make_error("Editor not active");
  }

  CEditor *I = G->Editor;
  int sele0 = SelectorIndexByName(G, cEditorSele1);
  ObjectMolecule *obj0 = SelectorGetFastSingleObjectMolecule(G, sele0);
  ObjectMoleculeVerifyChemistry(obj0, -1);    /* remember chemistry for later */
  if(!((sele0 >= 0) && obj0)) {
    return pymol::make_error("Invalid pk selection");
  }
  int sele1 = SelectorIndexByName(G, cEditorSele2);
  ObjectMolecule *obj1 = SelectorGetFastSingleObjectMolecule(G, sele1);
  if((sele1 >= 0) && (obj0 == obj1) && I->BondMode) {
    /* bond mode */
    ObjectMoleculeRemoveBonds(obj0, sele0, sele1);
    EditorInactivate(G);
  } else {
    int h_flag = false;

    if(SelectorIndexByName(G, cEditorSet) < 0) {
      int i0 = 0;
      /* only one atom picked */

      if(hydrogen) {
        auto buf = pymol::string_format("((neighbor %s) and hydro)", cEditorSele1);
        h_flag = SelectorCreate(G, cEditorRemoveSele, buf.c_str(), NULL, false, NULL).result();
      }

      if(SelectorGetFastSingleAtomObjectIndex(G, sele0, &i0)) {
        /* atom mode */
        if(i0 >= 0) {
          ExecutiveRemoveAtoms(G, cEditorSele1, quiet);
        }
      }
    } else {                /* multiple atoms picked */

      if(hydrogen) {
        auto buf = pymol::string_format("((neighbor %s) and hydro)", cEditorSet);
        h_flag = SelectorCreate(G, cEditorRemoveSele, buf.c_str(), NULL, false, NULL).result();
      }
      ExecutiveRemoveAtoms(G, cEditorSet, quiet);
    }

    EditorInactivate(G);
    if(h_flag) {
      ExecutiveRemoveAtoms(G, cEditorRemoveSele, quiet);
      SelectorDelete(G, cEditorRemoveSele);
    }
  }
#undef cEditorRemoveSele
  return {};
}


/*========================================================================*/
pymol::Result<> EditorHFill(PyMOLGlobals * G, int quiet)
{
  int sele0, sele1;
  int i0;
  ObjectMolecule *obj0 = NULL, *obj1 = NULL;

  if(EditorActive(G)) {
    sele0 = SelectorIndexByName(G, cEditorSele1);
    obj0 = SelectorGetFastSingleObjectMolecule(G, sele0);
    ObjectMoleculeVerifyChemistry(obj0, -1);    /* remember chemistry for later */
    if(sele0 >= 0) {
      
      sele1 = SelectorIndexByName(G, cEditorSele2);
      if(sele0 >= 0) {
        std::string s1, s2;
	if(sele1 >= 0){
          s2 = pymol::string_format("%s|%s", cEditorSele1, cEditorSele2);
          s1 = pymol::string_format("(neighbor (%s)) & hydro & !(%s)", s2, s2);
	} else {
          s2 = cEditorSele1;
          s1 = pymol::string_format("(neighbor (%s)) & hydro", s2);
	}
	
	ExecutiveRemoveAtoms(G, s1.c_str(), quiet);
	i0 = ObjectMoleculeGetAtomIndex(obj0, sele0);
	obj0->AtomInfo[i0].chemFlag = false;
        ExecutiveAddHydrogens(G, cEditorSele1, quiet);

	if(sele1 >= 0) {
	  obj1 = SelectorGetFastSingleObjectMolecule(G, sele1);
	  i0 = ObjectMoleculeGetAtomIndex(obj1, sele1);
	  obj1->AtomInfo[i0].chemFlag = false;
          ExecutiveAddHydrogens(G, cEditorSele2, quiet);
	}
      }
    } else {
      return pymol::Error("Nothing picked.");
    }
  } else {
    return pymol::Error("Editor not active.");
  }
  return {};
}


/*========================================================================*/
pymol::Result<> EditorHFix(PyMOLGlobals * G, const char *sele, int quiet)
{
  int sele0, sele1;
  ObjectMolecule *obj0, *obj1;
  if((!sele) || (!sele[0])) {   /* if selection is empty, then apply to picked atoms */
    if(EditorActive(G)) {
      sele0 = SelectorIndexByName(G, cEditorSele1);
      if(sele0 >= 0) {
        obj0 = SelectorGetFastSingleObjectMolecule(G, sele0);
        ObjectMoleculeVerifyChemistry(obj0, -1);
        ExecutiveFixHydrogens(G, cEditorSele1, quiet);
      }
      sele1 = SelectorIndexByName(G, cEditorSele2);
      if(sele1 >= 0) {
        obj1 = SelectorGetFastSingleObjectMolecule(G, sele1);
        ObjectMoleculeVerifyChemistry(obj1, -1);
        ExecutiveFixHydrogens(G, cEditorSele2, quiet);
      }
    } else {
      return pymol::Error("No valid selection and active editor.");
    }
  } else {
    ExecutiveFixHydrogens(G, sele, quiet);
  }
  return {};
}


/*========================================================================*/
pymol::Result<> EditorReplace(PyMOLGlobals* G, const char* elem, int geom,
    int valence, const char* name, int quiet)
{
  int i0;
  int sele0;
  AtomInfoType ai;
  ObjectMolecule *obj0 = NULL;
  int ok = true;
  UtilZeroMem(&ai, sizeof(AtomInfoType));
  if(EditorActive(G)) {
    for (const char* eSele : {cEditorSele2, cEditorSele3, cEditorSele4}) {
      if (SelectorIndexByName(G, eSele) >= 0) {
        return pymol::make_error("Only one picked selection allowed.");
      }
    }

    sele0 = SelectorIndexByName(G, cEditorSele1);
    obj0 = SelectorGetFastSingleObjectMolecule(G, sele0);
    if(obj0->DiscreteFlag) {
      return pymol::make_error("Can't attach atoms onto discrete objects.");
    } else {
      ObjectMoleculeVerifyChemistry(obj0, -1);  /* remember chemistry for later */
      if(sele0 >= 0) {
        i0 = ObjectMoleculeGetAtomIndex(obj0, sele0);   /* slow */
        if(i0 >= 0) {
          UtilNCopy(ai.elem, elem, sizeof(ElemName));
          if(name[0])
            LexAssign(G, ai.name, name);
          ai.geom = geom;
          ai.valence = valence;
          if (ok)
	    ok &= ObjectMoleculePrepareAtom(obj0, i0, &ai);
          if (ok)
	    ok &= ObjectMoleculePreposReplAtom(obj0, i0, &ai);
          ObjectMoleculeReplaceAtom(obj0, i0, std::move(ai));
          ObjectMoleculeVerifyChemistry(obj0, -1);
          ObjectMoleculeFillOpenValences(obj0, i0);
          if (ok)
	    ok &= ObjectMoleculeSort(obj0);
          ObjectMoleculeUpdateIDNumbers(obj0);
          EditorInactivate(G);
        }
      }
    }
  }
  if(ok) {
  return {};
  } else {
    return pymol::make_error("Could not replace atom.");
  }
}

static void draw_bond(PyMOLGlobals * G, float *v0, float *v1, CGO *shaderCGO)
{

  float v[3], v2[3], v3[3];
  float d0[3], n0[3], n1[3], n2[3];
  float x[50], y[50];
  int nEdge;
  int c, a;
  float tube_size1 = 0.5F;
  float tube_size3 = 0.45F;

  nEdge = SettingGetGlobal_i(G, cSetting_stick_quality) * 2;
  if(nEdge > 50)
    nEdge = 50;
  if(nEdge < 3)
    nEdge = 3;

  subdivide(nEdge, x, y);

  subtract3f(v1, v0, d0);
  average3f(v1, v0, v2);
  average3f(v0, v2, v3);
  average3f(v2, v3, v2);
  copy3f(d0, n0);
  get_system1f3f(n0, n1, n2);
  if (shaderCGO){
    CGOColorv(shaderCGO, ColorGet(G, 0));
    CGOBegin(shaderCGO, GL_TRIANGLE_STRIP);
    for(a = 0; a <= nEdge; a++) {
      c = a % nEdge;
      v[0] = n1[0] * x[c] + n2[0] * y[c];
      v[1] = n1[1] * x[c] + n2[1] * y[c];
      v[2] = n1[2] * x[c] + n2[2] * y[c];
      normalize3f(v);
      CGONormalv(shaderCGO, v);
      v[0] = v2[0] + n1[0] * tube_size1 * x[c] + n2[0] * tube_size1 * y[c];
      v[1] = v2[1] + n1[1] * tube_size1 * x[c] + n2[1] * tube_size1 * y[c];
      v[2] = v2[2] + n1[2] * tube_size1 * x[c] + n2[2] * tube_size1 * y[c];
      CGOVertexv(shaderCGO, v);
      
      v[0] = v3[0] + n1[0] * tube_size1 * x[c] + n2[0] * tube_size1 * y[c];
      v[1] = v3[1] + n1[1] * tube_size1 * x[c] + n2[1] * tube_size1 * y[c];
      v[2] = v3[2] + n1[2] * tube_size1 * x[c] + n2[2] * tube_size1 * y[c];
      CGOVertexv(shaderCGO, v);
    }
    CGOEnd(shaderCGO);

    CGOBegin(shaderCGO, GL_TRIANGLE_STRIP);
    CGONormalv(shaderCGO, n0);
    for(a = 0; a <= nEdge; a++) {
      c = a % nEdge;
      v[0] = v2[0] + n1[0] * tube_size3 * x[c] + n2[0] * tube_size3 * y[c];
      v[1] = v2[1] + n1[1] * tube_size3 * x[c] + n2[1] * tube_size3 * y[c];
      v[2] = v2[2] + n1[2] * tube_size3 * x[c] + n2[2] * tube_size3 * y[c];
      CGOVertexv(shaderCGO, v);
      v[0] = v2[0] + n1[0] * tube_size1 * x[c] + n2[0] * tube_size1 * y[c];
      v[1] = v2[1] + n1[1] * tube_size1 * x[c] + n2[1] * tube_size1 * y[c];
      v[2] = v2[2] + n1[2] * tube_size1 * x[c] + n2[2] * tube_size1 * y[c];
      CGOVertexv(shaderCGO, v);
    }
    CGOEnd(shaderCGO);

    CGOBegin(shaderCGO, GL_TRIANGLE_STRIP);
    scale3f(n0, -1.0F, v);
    CGONormalv(shaderCGO, v);
    for(a = 0; a <= nEdge; a++) {
      c = a % nEdge;
      v[0] = v3[0] + n1[0] * tube_size1 * x[c] + n2[0] * tube_size1 * y[c];
      v[1] = v3[1] + n1[1] * tube_size1 * x[c] + n2[1] * tube_size1 * y[c];
      v[2] = v3[2] + n1[2] * tube_size1 * x[c] + n2[2] * tube_size1 * y[c];
      CGOVertexv(shaderCGO, v);
      v[0] = v3[0] + n1[0] * tube_size3 * x[c] + n2[0] * tube_size3 * y[c];
      v[1] = v3[1] + n1[1] * tube_size3 * x[c] + n2[1] * tube_size3 * y[c];
      v[2] = v3[2] + n1[2] * tube_size3 * x[c] + n2[2] * tube_size3 * y[c];
      CGOVertexv(shaderCGO, v);
    }
    CGOEnd(shaderCGO);
  } else {
#ifdef PURE_OPENGL_ES_2
    /* TODO */
#else
  glColor3fv(ColorGet(G, 0));
  glBegin(GL_TRIANGLE_STRIP);
  for(a = 0; a <= nEdge; a++) {
    c = a % nEdge;
    v[0] = n1[0] * x[c] + n2[0] * y[c];
    v[1] = n1[1] * x[c] + n2[1] * y[c];
    v[2] = n1[2] * x[c] + n2[2] * y[c];
    normalize3f(v);
    glNormal3fv(v);
    v[0] = v2[0] + n1[0] * tube_size1 * x[c] + n2[0] * tube_size1 * y[c];
    v[1] = v2[1] + n1[1] * tube_size1 * x[c] + n2[1] * tube_size1 * y[c];
    v[2] = v2[2] + n1[2] * tube_size1 * x[c] + n2[2] * tube_size1 * y[c];
    glVertex3fv(v);

    v[0] = v3[0] + n1[0] * tube_size1 * x[c] + n2[0] * tube_size1 * y[c];
    v[1] = v3[1] + n1[1] * tube_size1 * x[c] + n2[1] * tube_size1 * y[c];
    v[2] = v3[2] + n1[2] * tube_size1 * x[c] + n2[2] * tube_size1 * y[c];
    glVertex3fv(v);
  }
  glEnd();
#endif

#ifdef PURE_OPENGL_ES_2
    /* TODO */
#else
  glBegin(GL_TRIANGLE_STRIP);
  glNormal3fv(n0);
  for(a = 0; a <= nEdge; a++) {
    c = a % nEdge;
    v[0] = v2[0] + n1[0] * tube_size3 * x[c] + n2[0] * tube_size3 * y[c];
    v[1] = v2[1] + n1[1] * tube_size3 * x[c] + n2[1] * tube_size3 * y[c];
    v[2] = v2[2] + n1[2] * tube_size3 * x[c] + n2[2] * tube_size3 * y[c];
    glVertex3fv(v);
    v[0] = v2[0] + n1[0] * tube_size1 * x[c] + n2[0] * tube_size1 * y[c];
    v[1] = v2[1] + n1[1] * tube_size1 * x[c] + n2[1] * tube_size1 * y[c];
    v[2] = v2[2] + n1[2] * tube_size1 * x[c] + n2[2] * tube_size1 * y[c];
    glVertex3fv(v);
  }
  glEnd();
#endif

#ifdef PURE_OPENGL_ES_2
    /* TODO */
#else
  glBegin(GL_TRIANGLE_STRIP);
  scale3f(n0, -1.0F, v);
  glNormal3fv(v);
  for(a = 0; a <= nEdge; a++) {
    c = a % nEdge;
    v[0] = v3[0] + n1[0] * tube_size1 * x[c] + n2[0] * tube_size1 * y[c];
    v[1] = v3[1] + n1[1] * tube_size1 * x[c] + n2[1] * tube_size1 * y[c];
    v[2] = v3[2] + n1[2] * tube_size1 * x[c] + n2[2] * tube_size1 * y[c];
    glVertex3fv(v);
    v[0] = v3[0] + n1[0] * tube_size3 * x[c] + n2[0] * tube_size3 * y[c];
    v[1] = v3[1] + n1[1] * tube_size3 * x[c] + n2[1] * tube_size3 * y[c];
    v[2] = v3[2] + n1[2] * tube_size3 * x[c] + n2[2] * tube_size3 * y[c];
    glVertex3fv(v);
  }
  glEnd();
#endif
  }
}

static void draw_globe(PyMOLGlobals * G, float *v2, int number, CGO *shaderCGO)
{
  float v[3];
  float n0[3], n1[3], n2[3];
  float x[50], y[50];
  int nEdge;
  int a, c;
  float radius = 0.5F;
  float width_base = 0.10F;
  float width = 0.0F;
  float offset = 0.0F;
  int cycle_counter;

  nEdge = SettingGetGlobal_i(G, cSetting_stick_quality) * 2;
  if(nEdge > 50)
    nEdge = 50;
  if(nEdge < 3)
    nEdge = 3;

  subdivide(nEdge, x, y);

  n0[0] = 1.0;
  n0[1] = 0.0;
  n0[2] = 0.0;
  get_system1f3f(n0, n1, n2);

#ifndef PURE_OPENGL_ES_2
  glColor3fv(ColorGet(G, 0));
#endif
  if (shaderCGO)
    CGOColorv(shaderCGO, ColorGet(G, 0));
  cycle_counter = number;
  while(cycle_counter) {

    switch (number) {
    case 1:
      width = width_base;
      offset = 0.0F;
      break;

    case 2:
      switch (cycle_counter) {
      case 2:
        width = width_base / 2;
        offset = width_base;
        break;
      case 1:
        offset = -width_base;
        break;
      }
      break;

    case 3:
      switch (cycle_counter) {
      case 3:
        width = width_base / 2.8F;
        offset = 1.33F * width_base;
        break;
      case 2:
        offset = 0.0F;
        break;
      case 1:
        offset = -1.33F * width_base;
        break;
      }
      break;

    case 4:
      switch (cycle_counter) {
      case 4:
        width = width_base / 3.2F;
        offset = 2 * width_base;
        break;
      case 3:
        offset = 0.66F * width_base;
        break;
      case 2:
        offset = -0.66F * width_base;
        break;
      case 1:
        offset = -2 * width_base;
        break;
      }
    }
    if (shaderCGO){
      CGOBegin(shaderCGO, GL_TRIANGLE_STRIP);
      for(a = 0; a <= nEdge; a++) {
	c = a % nEdge;
	v[0] = n1[0] * x[c] + n2[0] * y[c];
	v[1] = n1[1] * x[c] + n2[1] * y[c];
	v[2] = n1[2] * x[c] + n2[2] * y[c];
	normalize3f(v);
	CGONormalv(shaderCGO, v);
	v[0] =
	  v2[0] + n1[0] * radius * x[c] + n2[0] * radius * y[c] + n0[0] * (offset + width);
	v[1] =
	  v2[1] + n1[1] * radius * x[c] + n2[1] * radius * y[c] + n0[1] * (offset + width);
	v[2] =
	  v2[2] + n1[2] * radius * x[c] + n2[2] * radius * y[c] + n0[2] * (offset + width);
	CGOVertexv(shaderCGO, v);
	v[0] =
	  v2[0] + n1[0] * radius * x[c] + n2[0] * radius * y[c] + n0[0] * (offset - width);
	v[1] =
	  v2[1] + n1[1] * radius * x[c] + n2[1] * radius * y[c] + n0[1] * (offset - width);
	v[2] =
	  v2[2] + n1[2] * radius * x[c] + n2[2] * radius * y[c] + n0[2] * (offset - width);
	CGOVertexv(shaderCGO, v);
      }
      CGOEnd(shaderCGO);

      CGOBegin(shaderCGO, GL_TRIANGLE_STRIP);
      for(a = 0; a <= nEdge; a++) {
	c = a % nEdge;
	v[0] = n2[0] * x[c] + n0[0] * y[c];
	v[1] = n2[1] * x[c] + n0[1] * y[c];
	v[2] = n2[2] * x[c] + n0[2] * y[c];
	normalize3f(v);
	CGONormalv(shaderCGO, v);
	v[0] =
	  v2[0] + n2[0] * radius * x[c] + n0[0] * radius * y[c] + n1[0] * (offset + width);
	v[1] =
	  v2[1] + n2[1] * radius * x[c] + n0[1] * radius * y[c] + n1[1] * (offset + width);
	v[2] =
	  v2[2] + n2[2] * radius * x[c] + n0[2] * radius * y[c] + n1[2] * (offset + width);
	CGOVertexv(shaderCGO, v);
	v[0] =
	  v2[0] + n2[0] * radius * x[c] + n0[0] * radius * y[c] + n1[0] * (offset - width);
	v[1] =
        v2[1] + n2[1] * radius * x[c] + n0[1] * radius * y[c] + n1[1] * (offset - width);
	v[2] =
	  v2[2] + n2[2] * radius * x[c] + n0[2] * radius * y[c] + n1[2] * (offset - width);
	CGOVertexv(shaderCGO, v);
      }
      CGOEnd(shaderCGO);
      
      CGOBegin(shaderCGO, GL_TRIANGLE_STRIP);
      for(a = 0; a <= nEdge; a++) {
	c = a % nEdge;
	v[0] = n0[0] * x[c] + n1[0] * y[c];
	v[1] = n0[1] * x[c] + n1[1] * y[c];
	v[2] = n0[2] * x[c] + n1[2] * y[c];
	normalize3f(v);
	CGONormalv(shaderCGO, v);
	v[0] =
	  v2[0] + n0[0] * radius * x[c] + n1[0] * radius * y[c] + n2[0] * (offset + width);
	v[1] =
	  v2[1] + n0[1] * radius * x[c] + n1[1] * radius * y[c] + n2[1] * (offset + width);
	v[2] =
	  v2[2] + n0[2] * radius * x[c] + n1[2] * radius * y[c] + n2[2] * (offset + width);
	CGOVertexv(shaderCGO, v);
	v[0] =
	  v2[0] + n0[0] * radius * x[c] + n1[0] * radius * y[c] + n2[0] * (offset - width);
	v[1] =
	  v2[1] + n0[1] * radius * x[c] + n1[1] * radius * y[c] + n2[1] * (offset - width);
	v[2] =
	  v2[2] + n0[2] * radius * x[c] + n1[2] * radius * y[c] + n2[2] * (offset - width);
	CGOVertexv(shaderCGO, v);
      }
      CGOEnd(shaderCGO);
    } else {
#ifdef PURE_OPENGL_ES_2
    /* TODO */
#else
    glBegin(GL_TRIANGLE_STRIP);
    for(a = 0; a <= nEdge; a++) {
      c = a % nEdge;
      v[0] = n1[0] * x[c] + n2[0] * y[c];
      v[1] = n1[1] * x[c] + n2[1] * y[c];
      v[2] = n1[2] * x[c] + n2[2] * y[c];
      normalize3f(v);
      glNormal3fv(v);
      v[0] =
        v2[0] + n1[0] * radius * x[c] + n2[0] * radius * y[c] + n0[0] * (offset + width);
      v[1] =
        v2[1] + n1[1] * radius * x[c] + n2[1] * radius * y[c] + n0[1] * (offset + width);
      v[2] =
        v2[2] + n1[2] * radius * x[c] + n2[2] * radius * y[c] + n0[2] * (offset + width);
      glVertex3fv(v);
      v[0] =
        v2[0] + n1[0] * radius * x[c] + n2[0] * radius * y[c] + n0[0] * (offset - width);
      v[1] =
        v2[1] + n1[1] * radius * x[c] + n2[1] * radius * y[c] + n0[1] * (offset - width);
      v[2] =
        v2[2] + n1[2] * radius * x[c] + n2[2] * radius * y[c] + n0[2] * (offset - width);
      glVertex3fv(v);
    }
    glEnd();
#endif

#ifdef PURE_OPENGL_ES_2
    /* TODO */
#else
    glBegin(GL_TRIANGLE_STRIP);
    for(a = 0; a <= nEdge; a++) {
      c = a % nEdge;
      v[0] = n2[0] * x[c] + n0[0] * y[c];
      v[1] = n2[1] * x[c] + n0[1] * y[c];
      v[2] = n2[2] * x[c] + n0[2] * y[c];
      normalize3f(v);
      glNormal3fv(v);
      v[0] =
        v2[0] + n2[0] * radius * x[c] + n0[0] * radius * y[c] + n1[0] * (offset + width);
      v[1] =
        v2[1] + n2[1] * radius * x[c] + n0[1] * radius * y[c] + n1[1] * (offset + width);
      v[2] =
        v2[2] + n2[2] * radius * x[c] + n0[2] * radius * y[c] + n1[2] * (offset + width);
      glVertex3fv(v);
      v[0] =
        v2[0] + n2[0] * radius * x[c] + n0[0] * radius * y[c] + n1[0] * (offset - width);
      v[1] =
        v2[1] + n2[1] * radius * x[c] + n0[1] * radius * y[c] + n1[1] * (offset - width);
      v[2] =
        v2[2] + n2[2] * radius * x[c] + n0[2] * radius * y[c] + n1[2] * (offset - width);
      glVertex3fv(v);
    }
    glEnd();
#endif

#ifdef PURE_OPENGL_ES_2
    /* TODO */
#else
    glBegin(GL_TRIANGLE_STRIP);
    for(a = 0; a <= nEdge; a++) {
      c = a % nEdge;
      v[0] = n0[0] * x[c] + n1[0] * y[c];
      v[1] = n0[1] * x[c] + n1[1] * y[c];
      v[2] = n0[2] * x[c] + n1[2] * y[c];
      normalize3f(v);
      glNormal3fv(v);
      v[0] =
        v2[0] + n0[0] * radius * x[c] + n1[0] * radius * y[c] + n2[0] * (offset + width);
      v[1] =
        v2[1] + n0[1] * radius * x[c] + n1[1] * radius * y[c] + n2[1] * (offset + width);
      v[2] =
        v2[2] + n0[2] * radius * x[c] + n1[2] * radius * y[c] + n2[2] * (offset + width);
      glVertex3fv(v);
      v[0] =
        v2[0] + n0[0] * radius * x[c] + n1[0] * radius * y[c] + n2[0] * (offset - width);
      v[1] =
        v2[1] + n0[1] * radius * x[c] + n1[1] * radius * y[c] + n2[1] * (offset - width);
      v[2] =
        v2[2] + n0[2] * radius * x[c] + n1[2] * radius * y[c] + n2[2] * (offset - width);
      glVertex3fv(v);
    }
    glEnd();
#endif
    }
    cycle_counter--;
  }

}


/*
static void draw_string(float *v,char *l)
{
  glDisable(GL_DEPTH_TEST);	 
  glDisable(GL_LIGHTING);
  if(*l) {
    glColor3f(1.0,0.0,0.5);
    glRasterPos4f(v[0],v[1],v[2],1.0);
  }
  
  while(*l) {
    p_g lutBi tmapChar acter(P_G LUT_BITMAP_8_BY_13,*(l++));
  }

  glEnable(GL_LIGHTING);
  glEnable(GL_DEPTH_TEST);	 
}
*/


/*========================================================================*/
void EditorRender(PyMOLGlobals * G, int state)
{
  CEditor *I = G->Editor;
  int sele1, sele2, sele3, sele4;
  float v0[3], v1[3];
  float vp[12], *vv;
  /*  int v_cnt; */
  ObjectMolecule *obj1 = NULL, *obj2 = NULL, *obj3 = NULL, *obj4 = NULL;
  int index1, index2, index3, index4;
  int st, frozen;
  CGO *shaderCGO = NULL;

  if(EditorActive(G)) {
    int use_shader = SettingGetGlobal_b(G, cSetting_use_shaders);
    if (use_shader){
      if (!I->shaderCGO){
	shaderCGO = CGONew(G);
      } else {
	CGORenderGL(I->shaderCGO, NULL, NULL, NULL, NULL, NULL);
	return;
      }
    } else {
      CGOFree(I->shaderCGO);
    }

    PRINTFD(G, FB_Editor)
      " EditorRender-Debug: rendering...\n" ENDFD;

    if(G->HaveGUI && G->ValidContext) {

      sele1 = SelectorIndexByName(G, cEditorSele1);
      sele2 = SelectorIndexByName(G, cEditorSele2);
      sele3 = SelectorIndexByName(G, cEditorSele3);
      sele4 = SelectorIndexByName(G, cEditorSele4);

      obj1 = SelectorGetFastSingleAtomObjectIndex(G, sele1, &index1);
      obj2 = SelectorGetFastSingleAtomObjectIndex(G, sele2, &index2);
      obj3 = SelectorGetFastSingleAtomObjectIndex(G, sele3, &index3);
      obj4 = SelectorGetFastSingleAtomObjectIndex(G, sele4, &index4);

      /*      printf("%d %d %d %d\n",sele1,sele2,sele3,sele4);
         printf("%p %p %p %p\n",obj1,obj2,obj3,obj4);
         printf("%d %d %d %d\n",index1,index2,index3,index4); */

      if((sele1 >= 0) && (sele2 >= 0) && I->BondMode && obj1 && obj2) {
        /* bond mode */

        ObjectMoleculeGetAtomTxfVertex(obj1, state, index1, v0);
        ObjectMoleculeGetAtomTxfVertex(obj2, state, index2, v1);
        draw_bond(G, v0, v1, shaderCGO);

      } else {
        /* atom mode */

        vv = vp;

        if(obj1) {
	  /* if the user froze a state, use it instead of the global */
	  if((frozen = SettingGetIfDefined_i(obj1->G, obj1->Setting.get(), cSetting_state, &st))) {
	    state = st-1;
	  }
          if(ObjectMoleculeGetAtomTxfVertex(obj1, state, index1, vv)) {
            draw_globe(G, vv, 1, shaderCGO);
            vv += 3;
          }
        }

        if(obj2) {
	  if((frozen = SettingGetIfDefined_i(obj2->G, obj2->Setting.get(), cSetting_state, &st))) {
	    state = st-1;
	  }
          if(ObjectMoleculeGetAtomTxfVertex(obj2, state, index2, vv)) {
            draw_globe(G, vv, 2, shaderCGO);
            vv += 3;
          }
        }

        if(obj3) {
	  if((frozen = SettingGetIfDefined_i(obj3->G, obj3->Setting.get(), cSetting_state, &st))) {
	    state = st-1;
	  }
          if(ObjectMoleculeGetAtomTxfVertex(obj3, state, index3, vv)) {
            draw_globe(G, vv, 3, shaderCGO);
            vv += 3;
          }
        }

        if(obj4) {
	  if((frozen = SettingGetIfDefined_i(obj4->G, obj4->Setting.get(), cSetting_state, &st))) {
	    state = st-1;
	  }
          if(ObjectMoleculeGetAtomTxfVertex(obj4, state, index4, vv)) {
            draw_globe(G, vv, 4, shaderCGO);
            vv += 3;
          }
        }
      }
    }
    if (shaderCGO){
      int ok = true;
      CGOStop(shaderCGO);
      CHECKOK(ok, shaderCGO);
      if (ok){
        CGO* tmpCGO = CGONew(G);
        if (ok) ok &= CGOEnable(tmpCGO, GL_DEFAULT_SHADER);
	if (ok) ok &= CGODisable(tmpCGO, GL_TWO_SIDED_LIGHTING);
        tmpCGO->free_append(CGOOptimizeToVBONotIndexedNoShader(shaderCGO));
        if (ok) ok &= CGODisable(tmpCGO, GL_DEFAULT_SHADER);
	if (ok) ok &= CGOStop(tmpCGO);
	I->shaderCGO = tmpCGO;
	I->shaderCGO->use_shader = true;
      }
      CGOFree(shaderCGO);
      if (ok){
	CGORenderGL(I->shaderCGO, NULL, NULL, NULL, NULL, NULL);
      }
    }
  }
}


/*========================================================================*/
void EditorInactivate(PyMOLGlobals * G)
{
  CEditor *I = G->Editor;

  PRINTFD(G, FB_Editor)
    " EditorInactivate-Debug: callend.\n" ENDFD;

  if (I->Active) {
    // force refresh of the object menu panel (PYMOL-3411)
    OrthoInvalidateDoDraw(G);
  }

  I->DihedObject = NULL;
  I->DragObject = NULL;
  I->BondMode = false;
  I->ShowFrags = false;
  I->NFrag = 0;
  I->Active = false;
  SelectorDeletePrefixSet(G, cEditorFragPref);
  SelectorDeletePrefixSet(G, cEditorBasePref);
  ExecutiveDelete(G, cEditorSele1);
  ExecutiveDelete(G, cEditorSele2);
  ExecutiveDelete(G, cEditorSele3);
  ExecutiveDelete(G, cEditorSele4);
  ExecutiveDelete(G, cEditorSet);
  ExecutiveDelete(G, cEditorBond);
  ExecutiveDelete(G, cEditorRes);
  ExecutiveDelete(G, cEditorChain);
  ExecutiveDelete(G, cEditorObject);
  ExecutiveDelete(G, cEditorComp);
  ExecutiveDelete(G, cEditorLink);
  ExecutiveDelete(G, cEditorDihedral);
  ExecutiveDelete(G, cEditorDihe1);
  ExecutiveDelete(G, cEditorDihe2);
  ExecutiveDelete(G, cEditorMeasure);
  EditorMouseInvalid(G);
  EditorInvalidateShaderCGO(G);
  SceneInvalidate(G);
}


/*========================================================================*/
/**
 * Create a transient distance, angle, or dihedral measurement between
 * the pk1 - pk4 atoms.
 *
 * Assumes that the cEditorMeasure object does not exist yet.
 */
static
void EditorAutoMeasure(PyMOLGlobals * G,
    int sele1, int sele2, int sele3, int sele4, int state)
{
  if (sele1 < 0 || sele2 < 0)
    return;

  if (sele3 < 0) {
    ExecutiveDistance(G, cEditorMeasure, cEditorSele1, cEditorSele2, 0, -1.f, true,
        true, false /* reset */, state, false);
  } else if (sele4 < 0) {
    ExecutiveAngle(G, cEditorMeasure, cEditorSele1, cEditorSele2, cEditorSele3,
        0, true, false /* reset */, false, true, state);
  } else {
    ExecutiveDihedral(G, cEditorMeasure, cEditorSele1, cEditorSele2,
        cEditorSele3, cEditorSele4, 0, true, false /* reset */, false, true,
        state);
  }

  ExecutiveColor(G, cEditorMeasure, "gray", 0x1, true);
}


/*========================================================================*/
void EditorActivate(PyMOLGlobals * G, int state, int enable_bond)
{
  int sele1, sele2, sele3, sele4;

  CEditor *I = G->Editor;

  sele1 = SelectorIndexByName(G, cEditorSele1);
  sele2 = SelectorIndexByName(G, cEditorSele2);
  sele3 = SelectorIndexByName(G, cEditorSele3);
  sele4 = SelectorIndexByName(G, cEditorSele4);

  if((sele1 >= 0) || (sele2 >= 0) || (sele3 >= 0) || (sele4 >= 0)) {

    I->Active = true;
    ExecutiveDelete(G, cEditorComp);
    ExecutiveDelete(G, cEditorRes);
    ExecutiveDelete(G, cEditorChain);
    ExecutiveDelete(G, cEditorObject);
    ExecutiveDelete(G, cEditorBond);
    ExecutiveDelete(G, cEditorDihedral);
    ExecutiveDelete(G, cEditorDihe1);
    ExecutiveDelete(G, cEditorDihe2);
    ExecutiveDelete(G, cEditorMeasure);

    I->BondMode = enable_bond;
    I->NFrag = SelectorSubdivide(G, cEditorFragPref,
                                 sele1, sele2,
                                 sele3, sele4,
                                 cEditorBasePref, cEditorComp, &I->BondMode);
    /* just returns 'state' */
    state = EditorGetEffectiveState(G, NULL, state);
    I->ActiveState = state;

    I->ShowFrags = false;
    if(SettingGetGlobal_b(G, cSetting_auto_hide_selections))
      ExecutiveHideSelections(G);

    if(I->BondMode && SettingGetGlobal_b(G, cSetting_editor_auto_dihedral))
      EditorDihedralInvalid(G, NULL);

    if (!I->BondMode && SettingGetGlobal_b(G, cSetting_editor_auto_measure))
      EditorAutoMeasure(G, sele1, sele2, sele3, sele4, state);
  } else {
    EditorInactivate(G);
  }
  EditorMouseInvalid(G);
  EditorInvalidateShaderCGO(G);
}


/*========================================================================*/
void EditorSetDrag(PyMOLGlobals * G, pymol::CObject * obj, int sele, int quiet, int state)
{
  EditorInactivate(G);
  state = EditorGetEffectiveState(G, obj, state);
  if(obj->type == cObjectMolecule) {
    ObjectMolecule *objMol = (ObjectMolecule*)(void*)obj;
    if(ObjectMoleculeCheckFullStateSelection(objMol, sele, state)) {
      int matrix_mode = SettingGet_i(G, obj->Setting.get(), NULL, cSetting_matrix_mode);
      if(matrix_mode>=1) {
        /* force / coerce object matrix drags? */
        sele = -1;
      }
    }
  }
  EditorPrepareDrag(G, obj, sele, -1, state, 0);
}

void EditorReadyDrag(PyMOLGlobals * G, int state)
{
  CEditor *I = G->Editor;
  if(I->DragObject && (I->DragIndex == -1)) {
    EditorPrepareDrag(G, I->DragObject, I->DragSelection, -1, state, 0);
  }
}


/*========================================================================*/
void EditorPrepareDrag(PyMOLGlobals * G, pymol::CObject * obj,
                       int sele, int index, int state, int mode)
{
  int frg;
  int sele0 = -1, sele1 = -1, sele2 = -1, sele3 = -1;
  int s;
  WordType name;
  int seleFlag = false;
  int i0, i1, i2, i3;
  CEditor *I = G->Editor;
  int log_trans = SettingGetGlobal_b(G, cSetting_log_conformations);
  int drag_sele = -1;
  ObjectMolecule *objMol = NULL;

  PRINTFD(G, FB_Editor)
    " EditorPrepareDrag-Debug: entered. obj %p index %d\n", (void *) obj, index ENDFD;

  if(obj->type == cObjectMolecule)
    objMol = (ObjectMolecule*)(void*)obj;

  state = EditorGetEffectiveState(G, obj, state);

  /* TODO: if user is drags label, then the editor must be deactivated */

  if((!EditorActive(G))||(!objMol)) {
    /* non-anchored dragging of objects and now selections */

    float mn[3], mx[3];

    I->DragObject = obj;
    I->DragIndex = index;       /* set to -1 when in "mouse drag" mode */
    I->DragSelection = sele;
    I->DragHaveBase = false;
    if(sele >= 0) {
      auto sele_name = SelectorGetNameFromIndex(G, sele);
      if(sele_name) {
        strcpy(I->DragSeleName, sele_name);
        if(SettingGetGlobal_b(G, cSetting_editor_auto_origin)) {
          if(I->FavorOrigin) {
            I->DragHaveBase = true;
            copy3f(I->FavoredOrigin, I->DragBase);
          } else {
            if(ExecutiveGetExtent(G, sele_name, mn, mx, true, state, true)) {
              average3f(mn, mx, I->DragBase);
              I->DragHaveBase = true;
            }
          }
        }
      } else {
        I->DragSeleName[0] = 0;
      }
    } else {
      if(SettingGetGlobal_b(G, cSetting_editor_auto_origin)) {
        if(I->FavorOrigin) {
          I->DragHaveBase = true;
          copy3f(I->FavoredOrigin, I->DragBase);
        } else {
          if(ExecutiveGetExtent(G, obj->Name, mn, mx, true, state, true)) {
            average3f(mn, mx, I->DragBase);
            I->DragHaveBase = true;
          }
        }
      }
    }
  } else {
    
    /* anchored / fragment dragging  */
    for(frg = 1; frg <= I->NFrag; frg++) {
      sprintf(name, "%s%1d", cEditorFragPref, frg);
      drag_sele = SelectorIndexByName(G, name);
      if(drag_sele >= 0) {
        s = objMol->AtomInfo[index].selEntry;
        seleFlag = SelectorIsMember(G, s, drag_sele);
      }
      if(seleFlag) {
        strcpy(I->DragSeleName, name);
        break;
      }
    }
    if(seleFlag) {              /* normal selection */

      PRINTFB(G, FB_Editor, FB_Blather)
        " Editor: grabbing (%s).", name ENDFB(G);

      I->DragIndex = index;
      I->DragSelection = drag_sele;
      I->DragObject = obj;
      I->DragHaveAxis = false;
      I->DragHaveBase = false;
      I->DragBondFlag = false;
      I->DragSlowFlag = false;

      sprintf(name, "%s%1d", cEditorBasePref, frg);     /* get relevant base vertex of bond */
      sele1 = SelectorIndexByName(G, name);
      if(sele1 >= 0) {
        i1 = ObjectMoleculeGetAtomIndex(objMol, sele1);
        if(i1 >= 0) {
          ObjectMoleculeGetAtomTxfVertex(objMol, state, i1, I->DragBase);
          I->DragHaveBase = true;
          /*printf("base %s\n",name); */
        }
      }

      /* get axis or base atom */

      {
        int cnt = 0;

        if((sele0 = SelectorIndexByName(G, cEditorSele1)) >= 0) {
          if(SelectorIsAtomBondedToSele(G, objMol, sele0, drag_sele))
            cnt++;
          else
            sele0 = -1;
        }
        if((sele1 = SelectorIndexByName(G, cEditorSele2)) >= 0) {
          if(SelectorIsAtomBondedToSele(G, objMol, sele1, drag_sele))
            cnt++;
          else
            sele1 = -1;
        }
        if((sele2 = SelectorIndexByName(G, cEditorSele3)) >= 0) {
          if(SelectorIsAtomBondedToSele(G, objMol, sele2, drag_sele))
            cnt++;
          else
            sele2 = -1;
        }
        if((sele3 = SelectorIndexByName(G, cEditorSele4)) >= 0) {
          if(SelectorIsAtomBondedToSele(G, objMol, sele3, drag_sele))
            cnt++;
          else
            sele3 = -1;
        }

        i0 = ObjectMoleculeGetAtomIndex(objMol, sele0);
        i1 = ObjectMoleculeGetAtomIndex(objMol, sele1);
        i2 = ObjectMoleculeGetAtomIndex(objMol, sele2);
        i3 = ObjectMoleculeGetAtomIndex(objMol, sele3);

        if(cnt > 1) {           /* bond/multiatom mode */

          I->DragBondFlag = I->BondMode;

          zero3f(I->Center);
          if(i0 >= 0) {
            ObjectMoleculeGetAtomTxfVertex(objMol, state, i0, I->V0);
          } else if(i1 >= 0) {
            ObjectMoleculeGetAtomTxfVertex(objMol, state, i1, I->V0);
          } else if(i2 >= 0) {
            ObjectMoleculeGetAtomTxfVertex(objMol, state, i2, I->V0);
          } else if(i3 >= 0) {
            ObjectMoleculeGetAtomTxfVertex(objMol, state, i3, I->V0);
          }

          if(i0 >= 0) {
            ObjectMoleculeGetAtomTxfVertex(objMol, state, i0, I->V1);
            add3f(I->V1, I->Center, I->Center);
          }
          if(i1 >= 0) {
            ObjectMoleculeGetAtomTxfVertex(objMol, state, i1, I->V1);
            add3f(I->V1, I->Center, I->Center);
          }
          if(i2 >= 0) {
            ObjectMoleculeGetAtomTxfVertex(objMol, state, i2, I->V1);
            add3f(I->V1, I->Center, I->Center);
          }
          if(i3 >= 0) {
            ObjectMoleculeGetAtomTxfVertex(objMol, state, i3, I->V1);
            add3f(I->V1, I->Center, I->Center);
          }

          {
            float div = 1.0F / cnt;
            scale3f(I->Center, div, I->Center);
          }

          subtract3f(I->Center, I->V0, I->Axis);

          normalize3f(I->Axis);
          I->DragHaveAxis = true;

          if(SettingGetGlobal_b(G, cSetting_editor_auto_origin)) {
            if(I->FavorOrigin) {
              I->DragHaveBase = true;
              copy3f(I->FavoredOrigin, I->DragBase);
            } else {
              copy3f(I->Center, I->DragBase);
              I->DragHaveBase = true;
            }
          }

        } else {                /* atom mode */

          if(i0 >= 0) {
            ObjectMoleculeGetAtomTxfVertex(objMol, state, i0, I->V0);
          } else if(i1 >= 0) {
            ObjectMoleculeGetAtomTxfVertex(objMol, state, i1, I->V0);
          } else if(i2 >= 0) {
            ObjectMoleculeGetAtomTxfVertex(objMol, state, i2, I->V0);
          } else if(i3 >= 0) {
            ObjectMoleculeGetAtomTxfVertex(objMol, state, i3, I->V0);
          }
          if(I->DragHaveBase) {

            copy3f(I->DragBase, I->V1);
            subtract3f(I->V1, I->V0, I->Axis);
            average3f(I->V1, I->V0, I->Center);
            normalize3f(I->Axis);
            I->DragHaveAxis = true;
            if(mode == cButModeRotFrag) {
              copy3f(I->V0, I->DragBase);
            }

          }
        }
      }
    } else {                    /* clicked directly on an anchor atom */

      sele0 = SelectorIndexByName(G, cEditorSele1);
      if(sele0 < 0)
        sele0 = SelectorIndexByName(G, cEditorSele2);
      if(sele0 < 0)
        sele0 = SelectorIndexByName(G, cEditorSele3);
      if(sele0 < 0)
        sele0 = SelectorIndexByName(G, cEditorSele4);
      if(sele0 >= 0) {
        s = objMol->AtomInfo[index].selEntry;
        seleFlag = SelectorIsMember(G, s, sele0);
      }

      PRINTFB(G, FB_Editor, FB_Actions)
        " Editor: grabbing all fragments." ENDFB(G);
      I->DragIndex = index;
      I->DragSelection = SelectorIndexByName(G, cEditorComp);
      strcpy(I->DragSeleName, cEditorComp);
      I->DragObject = obj;
      I->DragHaveAxis = false;
      I->DragHaveBase = false;
      I->DragBondFlag = false;

      I->DragSlowFlag = true;

      if(sele0 >= 0) {          /* just provide a base vector, no valid axis exists */
        i1 = ObjectMoleculeGetAtomIndex(objMol, sele0);
        if(i1 >= 0) {
          ObjectMoleculeGetAtomTxfVertex(objMol, state, i1, I->DragBase);
          I->DragHaveBase = true;
          I->DragBondFlag = true;
        }
      }
    }
    if(!seleFlag) {
      I->DragIndex = -1;
      I->DragSelection = -1;
      I->DragObject = NULL;
    }
  }
  if(I->DragObject) {
    I->ShowFrags = false;
    if(objMol) {
      ObjectMoleculeSaveUndo(objMol, state, log_trans);
      if(SettingGetGlobal_b(G, cSetting_auto_sculpt)) {
        SettingSetGlobal_b(G, cSetting_sculpting, 1);
        if(!objMol->Sculpt)
          ObjectMoleculeSculptImprint(objMol, state, -1, 0);
      }
    }
  }
  if(log_trans)
    PLogFlush(G);

  PRINTFD(G, FB_Editor)
    " EditorPrepDrag-Debug: leaving Index %d Sele %d Object %p\n Axis %d Base %d BondFlag %d SlowFlag %d seleFlag %d\n",
    I->DragIndex, I->DragSelection, (void *) I->DragObject,
    I->DragHaveAxis, I->DragHaveBase, I->DragBondFlag, I->DragSlowFlag, seleFlag ENDFD;
}

int EditorDraggingObjectMatrix(PyMOLGlobals *G)
{
  CEditor *I = G->Editor;
  if(I->DragObject && (I->DragSelection < 0) && (I->DragIndex == -1)) {
      return true;
  }
  return false;
}

void EditorDrag(PyMOLGlobals * G, pymol::CObject * obj, int index, int mode, int state,
                float *pt, float *mov, float *z_dir)
{
  CEditor *I = G->Editor;
  float v0[3], v1[3], v2[3], v3[3], v4[4], cp[3];
  float d0[3], d1[3], d2[3], n0[3], n1[3], n2[3];
  float opp, adj, theta;
  float m[16];
  int log_trans = SettingGetGlobal_b(G, cSetting_log_conformations);

  PRINTFD(G, FB_Editor)
    " EditorDrag-Debug: entered. obj %p state %d index %d mode %d \nIndex %d Sele %d Object %p\n Axis %d Base %d BondFlag %d SlowFlag %d\n",
    (void *) obj, state, index, mode,
    I->DragIndex, I->DragSelection, (void *) I->DragObject,
    I->DragHaveAxis, I->DragHaveBase, I->DragBondFlag, I->DragSlowFlag ENDFD;

  if((index < 0) && (!obj))
    obj = I->DragObject;

  if(obj) {
    ObjectMolecule *objMol = NULL;
    if(obj->type == cObjectMolecule)
      objMol = (ObjectMolecule*)(void*)obj;

    state = EditorGetEffectiveState(G, obj, state);

    if((index == I->DragIndex) && (obj == I->DragObject)) {
      if(!EditorActive(G)) {
        int matrix_mode = SettingGet_i(G, I->DragObject->Setting.get(),
                                       NULL, cSetting_matrix_mode);
        if(matrix_mode<0)
          matrix_mode = EditorDraggingObjectMatrix(G) ? 1 : 0;
        
        /* always force use of matrix mode for non-molecular objects */
        if((!objMol)&&(!matrix_mode))
          matrix_mode = 1; 

        /* non-achored actions */
        switch (mode) {
        case cButModeRotDrag:
          if(I->DragHaveBase) {
            copy3f(I->DragBase, v3);
          } else {
            SceneOriginGet(G, v3);
          }

          get_rotation_about3f3fTTTf(pt[0], mov, v3, m);
          if(matrix_mode && (I->DragSelection < 0)) {
            switch (matrix_mode) {
            case 1:
              ObjectCombineTTT(obj, m, false, SettingGetGlobal_b(G,cSetting_movie_auto_store));
              break;
            case 2:
              if(objMol)
                ObjectMoleculeTransformState44f(objMol, state, m, log_trans, false, true);
              break;
            }
          } else {
            if(objMol)
              ObjectMoleculeTransformSelection(objMol, state, I->DragSelection,
                                               m, log_trans, I->DragSeleName, false, true);
          }
          SceneInvalidate(G);
          break;
        case cButModeRotFrag:
        case cButModeRotObj:
        case cButModeRotView:
          if(I->DragHaveBase) {
            copy3f(I->DragBase, v3);
          } else {
            SceneOriginGet(G, v3);
          }
          subtract3f(pt, v3, n0);
          add3f(pt, mov, n1);
          subtract3f(n1, v3, n1);
          normalize3f(n0);
          normalize3f(n1);
          cross_product3f(n0, n1, cp);
          theta = (float) asin(length3f(cp));
          normalize23f(cp, n2);
          get_rotation_about3f3fTTTf(theta, n2, v3, m);
          /* matrix m now contains a valid TTT rotation in global
             coordinate space that could be applied directly to the
             coordinates to effect the desired rotation */
          if(mode == cButModeRotView) {
            /* modify the object's TTT */
            ObjectCombineTTT(obj, m, false, 
                             SettingGetGlobal_b(G,cSetting_movie_auto_store));
          } else {
            if(matrix_mode) {
              switch (matrix_mode) {
              case 1:
                ObjectCombineTTT(obj, m, false, 
                                 SettingGetGlobal_b(G,cSetting_movie_auto_store));
                break;
              case 2:
                if(objMol)
                  ObjectMoleculeTransformState44f(objMol, state, m, log_trans, false, true);
                break;
              }
            } else {
              if(objMol)
                ObjectMoleculeTransformSelection(objMol, state, I->DragSelection,
                                                 m, log_trans, I->DragSeleName, false,
                                                 true);
            }
          }
          SceneInvalidate(G);
          break;
        case cButModeTorFrag:
          if(objMol) {
            ObjectMoleculeMoveAtom(objMol, state, index, mov, 1, log_trans);
            SceneInvalidate(G);
          }
          break;
        case cButModeMovView:
        case cButModeMovViewZ:
          ObjectTranslateTTT(obj, mov, SettingGetGlobal_b(G,cSetting_movie_auto_store));
          break;
        case cButModeMovObj:
        case cButModeMovObjZ:
        case cButModeMovFrag:
        case cButModeMovFragZ:
        case cButModeMovDrag:
        case cButModeMovDragZ:
          if(matrix_mode && (I->DragSelection < 0)) {
            identity44f(m);
            m[3] = mov[0];
            m[7] = mov[1];
            m[11] = mov[2];
            switch (matrix_mode) {
            case 1:
              ObjectCombineTTT(obj, m, false, SettingGetGlobal_b(G,cSetting_movie_auto_store));
              break;
            case 2:
              if(objMol)
                ObjectMoleculeTransformState44f(objMol, state, m, log_trans, true, true);
              break;
            }
          } else {
            identity44f(m);
            copy3f(mov, m + 12);        /* questionable... */
            if(objMol) 
              ObjectMoleculeTransformSelection(objMol, state, I->DragSelection, m,
                                               log_trans, I->DragSeleName, false, true);
          }
          SceneInvalidate(G);
          break;
        }
      } else {
        switch (mode) {
        case cButModeRotFrag:
        case cButModeRotObj:
          if(I->DragHaveBase) {
            copy3f(I->DragBase, v3);
          } else {
            copy3f(I->V0, v3);
          }
          if(I->DragSlowFlag) {
            SceneGetViewNormal(G, v4);
            scale3f(v4, -1.0F, v4);
            add3f(v3, v4, v4);
              subtract3f(pt, v4, n0);
            add3f(pt, mov, n1);
            subtract3f(n1, v4, n1);
          } else {
            subtract3f(pt, v3, n0);
            add3f(pt, mov, n1);
            subtract3f(n1, v3, n1);
          }
          normalize3f(n0);
          normalize3f(n1);
          cross_product3f(n0, n1, cp);
          theta = (float) asin(length3f(cp));
          normalize23f(cp, n2);

          get_rotation_about3f3fTTTf(theta, n2, v3, m);
          if(objMol) {
            ObjectMoleculeTransformSelection(objMol, state, I->DragSelection,
                                             m, log_trans, I->DragSeleName, false, true);
            SceneInvalidate(G);
          }
          break;
        case cButModeTorFrag:
        case cButModePkTorBnd:
          if(I->DragHaveAxis) {
            subtract3f(pt, I->Center, d0);
            if(dot_product3f(d0, I->Axis) < 0.0) {
              copy3f(I->V0, v1);
              copy3f(I->V1, v0);
            } else {
              copy3f(I->V0, v0);
              copy3f(I->V1, v1);
            }
            subtract3f(v1, v0, d1);
            copy3f(d1, n0);
            normalize3f(n0);
            cross_product3f(n0, d0, n1);
            normalize3f(n1);

            project3f(d0, n0, v2);
            add3f(I->Center, v2, v2);   /* v2 is the perpendicular point on the axis */
            subtract3f(pt, v2, d2);
            opp = (float) length3f(mov);
            adj = (float) length3f(d2);
            if(adj > R_SMALL4) {
              theta = (float) atan(opp / adj);
              if(dot_product3f(n1, mov) < 0.0)
                theta = -theta;
              get_rotation_about3f3fTTTf(theta, n0, v1, m);
              if(objMol)
                ObjectMoleculeTransformSelection(objMol, state, I->DragSelection, m,
                                                 log_trans, I->DragSeleName, false, true);
            } else {

              if(z_dir) {       /* NULL-safety */
                cross_product3f(I->Axis, z_dir, d0);
                theta = -dot_product3f(d0, mov);
                get_rotation_about3f3fTTTf(theta, n0, v1, m);
                if(objMol)
                  ObjectMoleculeTransformSelection(objMol, state, I->DragSelection, m,
                                                   log_trans, I->DragSeleName, false, true);
              }

            }
            if(I->BondMode && SettingGetGlobal_b(G, cSetting_editor_auto_dihedral))
              EditorDihedralInvalid(G, NULL);
          }

          SceneInvalidate(G);
          break;
        case cButModeMovFrag:
        case cButModeMovFragZ:
          identity44f(m);
          copy3f(mov, m + 12);  /* questionable */
          if(objMol)
            ObjectMoleculeTransformSelection(objMol, state, I->DragSelection,
                                           m, log_trans, I->DragSeleName, false, true);
          SceneInvalidate(G);
          break;
        }
      }
    }
    ExecutiveInvalidateSelectionIndicatorsCGO(G);
    EditorInvalidateShaderCGO(G);    
  }
  PRINTFD(G, FB_Editor)
    " EditorDrag-Debug: leaving...\n" ENDFD;

}

void EditorRemoveStale(PyMOLGlobals* G)
{
  if (!EditorActive(G)) {
    return;
  }

  for (const auto& sele :
      {cEditorSele1, cEditorSele2, cEditorSele3, cEditorSele4}) {
    auto tarSele = SelectorIndexByName(G, sele);
    if (tarSele > 0) {
      int at;
      auto obj = SelectorGetFastSingleAtomObjectIndex(G, tarSele, &at);
      if (!obj) {
        ExecutiveDelete(G, sele);
      }
    }
  }

  EditorActivate(G, -1, true);
}

/*========================================================================*/
int EditorInit(PyMOLGlobals * G)
{
  CEditor *I = NULL;
  if((I = (G->Editor = pymol::calloc<CEditor>(1)))) {

    I->DihedObject = NULL;
    I->NFrag = 0;
    I->Active = false;
    I->DragObject = NULL;
    I->DragIndex = -1;
    I->DragSelection = -1;
    I->NextPickSele = 0;
    I->BondMode = false;
    I->PosVLA = VLAlloc(float, 30);
    I->DihedralInvalid = false;
    I->MouseInvalid = false;
    I->FavorOrigin = false;
    I->shaderCGO = NULL;
    return 1;
  } else
    return 0;
}


/*========================================================================*/
void EditorFree(PyMOLGlobals * G)
{
  CEditor *I = G->Editor;
  VLAFreeP(I->PosVLA);
  FreeP(G->Editor);
}

void EditorInvalidateShaderCGO(PyMOLGlobals * G){
  CEditor *I = G->Editor;
  CGOFree(I->shaderCGO);
}
