
/* 
   A* -------------------------------------------------------------------
   B* This file contains source code for the PyMOL computer program
   C* Copyright (c) Schrodinger, LLC. 
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
#include "os_std.h"
#include "os_gl.h"

#include <clocale>

#include "MemoryDebug.h"

#include "Base.h"

#include "OVContext.h"

#include "MemoryDebug.h"
#include "Err.h"
#include "Util.h"
#include "Selector.h"
#include "Color.h"
#include "Ortho.h"
#include "Scene.h"
#include "PyMOLObject.h"
#include "Executive.h"
#include "Word.h"
#include "RepMesh.h"
#include "ObjectMolecule.h"
#include "Control.h"
#include "Sphere.h"
#include "Setting.h"
#include "Ray.h"
#include "Util.h"
#include "Movie.h"
#include "P.h"
#include "Editor.h"
#include "SculptCache.h"
#include "Isosurf.h"
#include "Tetsurf.h"
#include "PConv.h"
#include "VFont.h"
#include "Wizard.h"
#include "Text.h"
#include "Character.h"
#include "Seq.h"
#include "Seeker.h"
#include "Texture.h"
#include "TestPyMOL.h"
#include "TypeFace.h"
#include "PlugIOManager.h"
#include "MovieScene.h"
#include "Lex.h"
#include "SelectorDef.h"

#ifdef _PYMOL_OPENVR
#include "OpenVRMode.h"
#endif

#include "PyMOL.h"
#include "PyMOLGlobals.h"
#include "PyMOLOptions.h"
#include "Feedback.h"
#include "GraphicsUtil.h"
#include "pymol/zstring_view.h"

#include "ShaderMgr.h"
#include "Version.h"

#ifndef _PYMOL_NOPY
PyMOLGlobals *SingletonPyMOLGlobals = NULL;
#endif

#ifdef _PYMOL_LIB_HAS_PYTHON
#define PYMOL_API_LOCK if(I->PythonInitStage && (!I->ModalDraw)) { PLockAPIAndUnblock(I->G); {
#define PYMOL_API_LOCK_MODAL if(I->PythonInitStage) { PLockAPIAndUnblock(I->G); {
#define PYMOL_API_TRYLOCK if(I->PythonInitStage && (!I->ModalDraw)) { if(PTryLockAPIAndUnblock(I->G)) {
#define PYMOL_API_UNLOCK PBlockAndUnlockAPI(I->G); }}
#define PYMOL_API_UNLOCK_NO_FLUSH PBlockAndUnlockAPI(I->G); }}
#else
#define PYMOL_API_LOCK if(!I->ModalDraw) {
#define PYMOL_API_LOCK_MODAL {
#define PYMOL_API_TRYLOCK if(!I->ModalDraw) {
#define PYMOL_API_UNLOCK }
#define PYMOL_API_UNLOCK_NO_FLUSH }
#endif
#define IDLE_AND_READY 3

typedef struct _CPyMOL {
  PyMOLGlobals *G;
  int FakeDragFlag;
  int RedisplayFlag;
  int PassiveFlag;
  int SwapFlag;
  int BusyFlag;
  int InterruptFlag;
  int ReshapeFlag;
  int ClickReadyFlag;
  int DrawnFlag;
  ObjectNameType ClickedObject;
  int ClickedIndex, ClickedButton, ClickedModifiers, ClickedX, ClickedY, ClickedHavePos, ClickedPosState;
  int ClickedBondIndex;
  float ClickedPos[3];
  int ImageRequestedFlag, ImageReadyFlag;
  int DraggedFlag;
  int Reshape[PYMOL_RESHAPE_SIZE];
  int Progress[PYMOL_PROGRESS_SIZE];
  int ProgressChanged;
  int IdleAndReady;
  int ExpireCount;
  bool done_ConfigureShaders;

  PyMOLModalDrawFn *ModalDraw;

  PyMOLSwapBuffersFn *SwapFn;


/* Python stuff */
#ifndef _PYMOL_NOPY
  int PythonInitStage;
#endif
  /* dynamically mapped string constants */

  OVLexicon *Lex;
  OVOneToOne *Rep;
  ov_word lex_everything, lex_sticks, lex_spheres, lex_surface;
  ov_word lex_labels, lex_nb_spheres, lex_cartoon, lex_ribbon;
  ov_word lex_lines, lex_mesh, lex_dots, lex_dashes, lex_nonbonded;
  ov_word lex_cell, lex_cgo, lex_callback, lex_extent, lex_slice;

  OVOneToOne *Clip;
  ov_word lex_near, lex_far, lex_move, lex_slab, lex_atoms;

  OVOneToOne *Reinit;
  ov_word lex_settings;

  OVOneToOne *SelectList;
  ov_word lex_index, lex_id, lex_rank;

  OVOneToOne *Setting;

#ifdef _PYMOL_LIB
  OVOneToOne *MouseButtonCodeLexicon;
  ov_word lex_left, lex_middle, lex_right;
  ov_word lex_wheel;
  ov_word lex_double_left, lex_double_middle, lex_double_right;
  ov_word lex_single_left, lex_single_middle, lex_single_right;

  OVOneToOne *MouseButtonModCodeLexicon;
  ov_word lex_none, lex_shft, lex_ctrl, lex_ctsh;
  ov_word lex_alt, lex_alsh, lex_ctal, lex_ctas;

  OVOneToOne *MouseButtonActionCodeLexicon;
  ov_word lex_but_rota, lex_but_move, lex_but_movz, lex_but_clip, lex_but_rotz;
  ov_word lex_but_clpn, lex_but_clpf, lex_but_lb, lex_but_mb, lex_but_rb;
  ov_word lex_but_plus_lb, lex_but_plus_mb, lex_but_plus_rb, lex_but_pkat, lex_but_pkbd;
  ov_word lex_but_rotf, lex_but_torf, lex_but_movf, lex_but_orig, lex_but_plus_lbx;
  ov_word lex_but_minus_lbx, lex_but_lbbx, lex_but_none, lex_but_cent, lex_but_pktb;
  ov_word lex_but_slab, lex_but_movs, lex_but_pk1;
  ov_word lex_but_mova, lex_but_menu, lex_but_sele, lex_but_plus_minus;
  ov_word lex_but_plus_box, lex_but_minus_box, lex_but_mvsz, lex_but_dgrt, lex_but_dgmv;
  ov_word lex_but_dgmz, lex_but_roto, lex_but_movo, lex_but_mvoz, lex_but_mvfz;
  ov_word lex_but_mvaz, lex_but_drgm, lex_but_rotv, lex_but_movv, lex_but_mvvz;
  ov_word lex_but_drgo, lex_but_imsz, lex_but_imvz, lex_but_box, lex_but_irtz;

  OVOneToOne *MouseModeLexicon;
#include "buttonmodes_lex_def.h"

  OVOneToOne *PaletteLexicon;
#include "palettes_lex_def.h"

#endif

  AtomPropertyInfo AtomPropertyInfos[NUM_ATOM_PROPERTIES];
  OVOneToOne *AtomPropertyLexicon;
  ov_word lex_atom_prop_model, lex_atom_prop_index, lex_atom_prop_type,
    lex_atom_prop_name, lex_atom_prop_resn, lex_atom_prop_resi,
    lex_atom_prop_resv, lex_atom_prop_chain, lex_atom_prop_alt,
    lex_atom_prop_segi, lex_atom_prop_elem,
    lex_atom_prop_ss, lex_atom_prop_text_type, 
    lex_atom_prop_custom, lex_atom_prop_label, lex_atom_prop_numeric_type,
    lex_atom_prop_q, lex_atom_prop_b, lex_atom_prop_vdw,
    lex_atom_prop_elec_radius, lex_atom_prop_partial_charge, lex_atom_prop_formal_charge,
    lex_atom_prop_stereo, lex_atom_prop_cartoon, lex_atom_prop_color,
    lex_atom_prop_ID, lex_atom_prop_rank, lex_atom_prop_flags,
    lex_atom_prop_geom, lex_atom_prop_valence,
    lex_atom_prop_x, lex_atom_prop_y, lex_atom_prop_z,
    lex_atom_prop_settings, lex_atom_prop_properties,
    lex_atom_prop_reps,
    lex_atom_prop_protons,
    lex_atom_prop_oneletter,
    lex_atom_prop_s, lex_atom_prop_p, lex_atom_prop_state; 
  /*
    lex_atom_prop_, lex_atom_prop_, lex_atom_prop_,
    lex_atom_prop_, lex_atom_prop_, lex_atom_prop_,*/

} _CPyMOL;


/* convenience functions -- inline */

inline PyMOLstatus get_status_ok(int ok)
{
  if(ok)
    return PyMOLstatus_SUCCESS;
  else
    return PyMOLstatus_FAILURE;
}

inline PyMOLreturn_status return_status_ok(int ok)
{
  PyMOLreturn_status result;
  result.status = get_status_ok(ok);
  return result;
}

inline PyMOLreturn_status return_status(int status)
{
  PyMOLreturn_status result;
  result.status = status;
  return result;
}

inline PyMOLreturn_float return_result(const pymol::Result<float>& res)
{
  PyMOLreturn_float result = {PyMOLstatus_FAILURE};
  if (res) {
    result.status = PyMOLstatus_SUCCESS;
    result.value = res.result();
  }
  return result;
}

static PyMOLreturn_string_array return_result(
    const pymol::Result<std::vector<const char*>>& res)
{
  PyMOLreturn_string_array result = {
      PyMOLstatus_SUCCESS, 0, nullptr};
  if (!res) {
    result.status = PyMOLstatus_FAILURE;
  } else if (!res.result().empty()) {
    const auto& vec = res.result();
    result.size = vec.size();
    result.array = VLAlloc(char*, result.size);

    // allocate space for concatenated string in first array element
    size_t reslen = 0;
    for (const char* s : vec) {
      reslen += strlen(s) + 1;
    }
    result.array[0] = VLAlloc(char, reslen);

    // copy elements
    for (size_t pl = 0, i = 0; i != vec.size(); ++i) {
      result.array[i] = result.array[0] + pl;
      strcpy(result.array[i], vec[i]);
      pl += strlen(vec[i]) + 1;
    }
  }
  return result;
}

#if defined(__cplusplus) && !defined(_WEBGL)
extern "C" {
#endif

#ifdef _PYMOL_LIB
int initial_button_modes[cButModeInputCount];
#endif

static OVstatus PyMOL_InitAPI(CPyMOL * I)
{
  OVContext *C = I->G->Context;
  OVreturn_word result;
  I->Lex = OVLexicon_New(C->heap);
  if(!I->Lex)
    return_OVstatus_FAILURE;

  /* the following preprocessor macros may require GNU's cpp or VC++
     we'll see... */

#define LEX(ARG)  \
  if(!OVreturn_IS_OK( (result= OVLexicon_GetFromCString(I->Lex,#ARG))))  \
    return_OVstatus_FAILURE \
    else \
      I -> lex_ ## ARG = result.word;

  /* string constants that are accepted on input */

#define LEX_REP(NAME,CODE) LEX(NAME) \
    if(!OVreturn_IS_OK( OVOneToOne_Set(I->Rep,I->lex_ ## NAME, CODE)))  \
      return_OVstatus_FAILURE;

  I->Rep = OVOneToOne_New(C->heap);
  if(!I->Rep)
    return_OVstatus_FAILURE;

  LEX_REP(everything, -1);
  LEX_REP(sticks, 0);
  LEX_REP(spheres, 1);
  LEX_REP(surface, 2);
  LEX_REP(labels, 3);
  LEX_REP(nb_spheres, 4);
  LEX_REP(cartoon, 5);
  LEX_REP(ribbon, 6);
  LEX_REP(lines, 7);
  LEX_REP(mesh, 8);
  LEX_REP(dots, 9);
  LEX_REP(dashes, 10);
  LEX_REP(nonbonded, 11);
  LEX_REP(cell, 12);
  LEX_REP(cgo, 13);
  LEX_REP(callback, 14);
  LEX_REP(extent, 15);
  LEX_REP(slice, 16);

  /* workaround for unexplained bug with nested macro on VC6 */

#define LEX_CLIP(NAME,CODE) {if(!OVreturn_IS_OK( (result= OVLexicon_GetFromCString(I->Lex,#NAME))))  \
    return_OVstatus_FAILURE \
    else \
    I -> lex_ ## NAME = result.word;} \
    if(!OVreturn_IS_OK( OVOneToOne_Set(I->Clip,I->lex_ ## NAME, CODE)))  \
      return_OVstatus_FAILURE;

  I->Clip = OVOneToOne_New(C->heap);
  if(!I->Clip)
    return_OVstatus_FAILURE;

  LEX_CLIP(near, 0);
  LEX_CLIP(far, 1);
  LEX_CLIP(move, 2);
  LEX_CLIP(slab, 3);
  LEX_CLIP(atoms, 4);

#define LEX_REINIT(NAME,CODE) {if(!OVreturn_IS_OK( (result= OVLexicon_GetFromCString(I->Lex,#NAME))))  \
    return_OVstatus_FAILURE \
    else \
    I -> lex_ ## NAME = result.word;} \
    if(!OVreturn_IS_OK( OVOneToOne_Set(I->Reinit,I->lex_ ## NAME, CODE)))  \
      return_OVstatus_FAILURE;

  I->Reinit = OVOneToOne_New(C->heap);
  if(!I->Reinit)
    return_OVstatus_FAILURE;

  LEX_REINIT(everything, 0);
  LEX_REINIT(settings, 1);

#define LEX_SELLIST(NAME,CODE) {if(!OVreturn_IS_OK( (result= OVLexicon_GetFromCString(I->Lex,#NAME))))  \
    return_OVstatus_FAILURE \
    else \
    I -> lex_ ## NAME = result.word;} \
    if(!OVreturn_IS_OK( OVOneToOne_Set(I->SelectList,I->lex_ ## NAME, CODE)))  \
      return_OVstatus_FAILURE;

  I->SelectList = OVOneToOne_New(C->heap);
  if(!I->SelectList)
    return_OVstatus_FAILURE;

  LEX_SELLIST(index, 0);
  LEX_SELLIST(id, 1);
  LEX_SELLIST(rank, 2);

  I->Setting = OVOneToOne_New(C->heap);
  if(!I->Setting)
    return_OVstatus_FAILURE;

  if(!CPyMOLInitSetting(I->Lex, I->Setting))
    return_OVstatus_FAILURE;

#ifdef _PYMOL_LIB

  I->MouseButtonCodeLexicon = OVOneToOne_New(C->heap);
  if(!I->MouseButtonCodeLexicon)
    return_OVstatus_FAILURE;

#define LEX_MOUSECODE(NAME,CODE) LEX(NAME) \
    if(!OVreturn_IS_OK( OVOneToOne_Set(I->MouseButtonCodeLexicon,I->lex_ ## NAME, CODE)))  \
      return_OVstatus_FAILURE;

  LEX_MOUSECODE(left, 0);
  LEX_MOUSECODE(middle, 1);
  LEX_MOUSECODE(right, 2);
  LEX_MOUSECODE(wheel, 3);
  LEX_MOUSECODE(double_left, 4);
  LEX_MOUSECODE(double_middle, 5);
  LEX_MOUSECODE(double_right, 6);
  LEX_MOUSECODE(single_left, 7);
  LEX_MOUSECODE(single_middle, 8);
  LEX_MOUSECODE(single_right, 9);

  I->MouseButtonModCodeLexicon = OVOneToOne_New(C->heap);
  if(!I->MouseButtonModCodeLexicon)
    return_OVstatus_FAILURE;

#define LEX_BUTTONMODCODE(NAME,CODE) LEX(NAME) \
    if(!OVreturn_IS_OK( OVOneToOne_Set(I->MouseButtonModCodeLexicon,I->lex_ ## NAME, CODE)))  \
      return_OVstatus_FAILURE;

  LEX_BUTTONMODCODE(none, 0);
  LEX_BUTTONMODCODE(shft, 1);
  LEX_BUTTONMODCODE(ctrl, 2);
  LEX_BUTTONMODCODE(ctsh, 3);
  LEX_BUTTONMODCODE(alt, 4);
  LEX_BUTTONMODCODE(alsh, 5);
  LEX_BUTTONMODCODE(ctal, 6);
  LEX_BUTTONMODCODE(ctas, 7);

  I->MouseButtonActionCodeLexicon = OVOneToOne_New(C->heap);
  if(!I->MouseButtonActionCodeLexicon)
    return_OVstatus_FAILURE;

#define LEX_BUT(ARG)  \
  if(!OVreturn_IS_OK( (result= OVLexicon_GetFromCString(I->Lex,#ARG))))  \
    return_OVstatus_FAILURE \
    else \
      I -> lex_but_ ## ARG = result.word;

#define LEX_BUTTONACTIONCODE(NAME,CODE) LEX_BUT(NAME) \
    if(!OVreturn_IS_OK( OVOneToOne_Set(I->MouseButtonActionCodeLexicon,I->lex_but_ ## NAME, CODE)))  \
      return_OVstatus_FAILURE;
#define LEX_BUTTONACTIONCODEWITHSTRING(NAME,STRARG,CODE) \
  if(!OVreturn_IS_OK( (result= OVLexicon_GetFromCString(I->Lex,STRARG))))  \
    return_OVstatus_FAILURE \
    else \
      I -> lex_but_ ## NAME = result.word; \
  if(!OVreturn_IS_OK( OVOneToOne_Set(I->MouseButtonActionCodeLexicon,I->lex_but_ ## NAME, CODE)))  \
    return_OVstatus_FAILURE;


  LEX_BUTTONACTIONCODE(rota, 0);
  LEX_BUTTONACTIONCODE(move, 1);
  LEX_BUTTONACTIONCODE(movz, 2);
  LEX_BUTTONACTIONCODE(clip, 3);
  LEX_BUTTONACTIONCODE(rotz, 4);
  LEX_BUTTONACTIONCODE(clpn, 5);
  LEX_BUTTONACTIONCODE(clpf, 6);
  LEX_BUTTONACTIONCODE(lb, 7);
  LEX_BUTTONACTIONCODE(mb, 8);
  LEX_BUTTONACTIONCODE(rb, 9);
  LEX_BUTTONACTIONCODEWITHSTRING(plus_lb, "+lb", 10);
  LEX_BUTTONACTIONCODEWITHSTRING(plus_mb, "+mb", 11);
  LEX_BUTTONACTIONCODEWITHSTRING(plus_rb, "+rb", 12);
  LEX_BUTTONACTIONCODE(pkat, 13);
  LEX_BUTTONACTIONCODE(pkbd, 14);
  LEX_BUTTONACTIONCODE(rotf, 15);
  LEX_BUTTONACTIONCODE(torf, 16);
  LEX_BUTTONACTIONCODE(movf, 17);
  LEX_BUTTONACTIONCODE(orig, 18);
  LEX_BUTTONACTIONCODEWITHSTRING(plus_lbx, "+lbx", 19);
  LEX_BUTTONACTIONCODEWITHSTRING(minus_lbx, "-lbx", 20);
  LEX_BUTTONACTIONCODE(lbbx, 21);
  LEX_BUTTONACTIONCODE(none, 22);
  LEX_BUTTONACTIONCODE(cent, 23);
  LEX_BUTTONACTIONCODE(pktb, 24);
  LEX_BUTTONACTIONCODE(slab, 25);
  LEX_BUTTONACTIONCODE(movs, 26);
  LEX_BUTTONACTIONCODE(pk1, 27);
  LEX_BUTTONACTIONCODE(mova, 28);
  LEX_BUTTONACTIONCODE(menu, 29);
  LEX_BUTTONACTIONCODE(sele, 30);
  LEX_BUTTONACTIONCODEWITHSTRING(plus_minus,"+/-", 31);
  LEX_BUTTONACTIONCODEWITHSTRING(plus_box, "+box", 32);
  LEX_BUTTONACTIONCODEWITHSTRING(minus_box, "-box", 33);
  LEX_BUTTONACTIONCODE(mvsz, 34);
  LEX_BUTTONACTIONCODE(dgrt, 36);
  LEX_BUTTONACTIONCODE(dgmv, 37);
  LEX_BUTTONACTIONCODE(dgmz, 38);
  LEX_BUTTONACTIONCODE(roto, 39);
  LEX_BUTTONACTIONCODE(movo, 40);
  LEX_BUTTONACTIONCODE(mvoz, 41);
  LEX_BUTTONACTIONCODE(mvfz, 42);
  LEX_BUTTONACTIONCODE(mvaz, 43);
  LEX_BUTTONACTIONCODE(drgm, 44);
  LEX_BUTTONACTIONCODE(rotv, 45);
  LEX_BUTTONACTIONCODE(movv, 46);
  LEX_BUTTONACTIONCODE(mvvz, 47);
  LEX_BUTTONACTIONCODE(drgo, 49);
  LEX_BUTTONACTIONCODE(imsz, 50);
  LEX_BUTTONACTIONCODE(imvz, 51);
  LEX_BUTTONACTIONCODE(box, 52);
  LEX_BUTTONACTIONCODE(irtz, 53);

  I->MouseModeLexicon = OVOneToOne_New(C->heap);
  if(!I->MouseModeLexicon)
    return_OVstatus_FAILURE;

#define LEX_MOUSEMODECODE(NAME,CODE) LEX(NAME) \
    if(!OVreturn_IS_OK( OVOneToOne_Set(I->MouseModeLexicon,I->lex_ ## NAME, CODE)))  \
      return_OVstatus_FAILURE;

#include "buttonmodes_lex_init.h"

  {
    int a;
    /* These are set by default for the modes, basically, any single or double
       click is a simple click (i.e., cButModeSimpleClick), mouse button
       actions are initialized to a potential click, wheel actions are set to none.
       This is very similar to what is done in PyMOL_SetMouseButtonMode(), 
       and it makes it easier to specify new modes without needing to set
       every mouse function */
    for(a = cButModeLeftDouble /* 16 */; a <= cButModeRightCtrlAltShftSingle /* 63 */; a++) {
      /* all single and double clicks */
      initial_button_modes[a] = cButModeSimpleClick;
    }
    for(a = cButModeLeftAlt /* 68 */; a <= cButModeRightCtrlAltShft /* 79 */; a++) {
      /* all button modes with Alt */
      initial_button_modes[a] = cButModePotentialClick;
    }
    for(a = cButModeLeftNone /* 0 */; a <= cButModeRightCtSh /* 11 */; a++) {
      /* all button modes without Alt */
      initial_button_modes[a] = cButModePotentialClick;
    }
    for(a = cButModeWheelNone /* 12 */; a <= cButModeWheelCtSh /* 15 */; a++) {
      initial_button_modes[a] = cButModeNone;
    }  
    for(a = cButModeWheelAlt /* 64 */; a <= cButModeWheelCtrlAltShft /* 67 */; a++) {
      initial_button_modes[a] = cButModeNone;
    }
  }

  I->PaletteLexicon = OVOneToOne_New(C->heap);
  if(!I->PaletteLexicon)
    return_OVstatus_FAILURE;

#define LEX_PALETTE(NAME,CODE) LEX(NAME) \
    if(!OVreturn_IS_OK( OVOneToOne_Set(I->PaletteLexicon,I->lex_ ## NAME, CODE)))  \
      return_OVstatus_FAILURE;

#include "palettes_lex_init.h"

#endif


  I->AtomPropertyLexicon = OVOneToOne_New(C->heap);
  if(!I->AtomPropertyLexicon)
    return_OVstatus_FAILURE;

#define LEX_ATM_PROP(ARG)  \
  if(!OVreturn_IS_OK( (result= OVLexicon_GetFromCString(I->Lex,#ARG))))  \
    return_OVstatus_FAILURE \
    else \
      I -> lex_atom_prop_ ## ARG = result.word;
#define LEX_ATOM_PROP(NAME,CODE,TYPE,OFFSET) LEX_ATM_PROP(NAME)		\
    if(!OVreturn_IS_OK( OVOneToOne_Set(I->AtomPropertyLexicon,I->lex_atom_prop_ ## NAME, CODE)))  \
      return_OVstatus_FAILURE;     \
    I->AtomPropertyInfos[CODE].id = CODE;    \
    I->AtomPropertyInfos[CODE].Ptype = TYPE;    \
    I->AtomPropertyInfos[CODE].offset = OFFSET;  \
    I->AtomPropertyInfos[CODE].maxlen = 0;

#define LEX_ATOM_PROP_S(NAME,CODE,TYPE,OFFSET,MAXLEN) LEX_ATM_PROP(NAME)	\
    if(!OVreturn_IS_OK( OVOneToOne_Set(I->AtomPropertyLexicon,I->lex_atom_prop_ ## NAME, CODE)))  \
      return_OVstatus_FAILURE;     \
    I->AtomPropertyInfos[CODE].id = CODE;    \
    I->AtomPropertyInfos[CODE].Ptype = TYPE;    \
    I->AtomPropertyInfos[CODE].offset = OFFSET;  \
    I->AtomPropertyInfos[CODE].maxlen = MAXLEN;

  /*TEMP*/
  LEX_ATOM_PROP(model, 0, cPType_model, 0);
  LEX_ATOM_PROP(index, 1, cPType_index, 0);
  LEX_ATOM_PROP(type, 2, cPType_char_as_type, 0);
  LEX_ATOM_PROP(name, 3, cPType_int_as_string, offsetof(AtomInfoType,name));
  LEX_ATOM_PROP(resn, 4, cPType_int_as_string, offsetof(AtomInfoType,resn));
  LEX_ATOM_PROP(resi, 5, 0, 0);
  LEX_ATOM_PROP(resv, 6, cPType_int, offsetof(AtomInfoType,resv));
  LEX_ATOM_PROP(chain, 7, cPType_int_as_string, offsetof(AtomInfoType,chain));
  LEX_ATOM_PROP_S(alt, 8, cPType_string, offsetof(AtomInfoType,alt), 1);
  LEX_ATOM_PROP(segi, 9, cPType_int_as_string, offsetof(AtomInfoType,segi));
  LEX_ATOM_PROP_S(elem, 10, cPType_string, offsetof(AtomInfoType,elem), cElemNameLen);
  LEX_ATOM_PROP_S(ss, 11, cPType_string, offsetof(AtomInfoType,ssType), 1);
  LEX_ATOM_PROP(text_type, 12, cPType_int_as_string, offsetof(AtomInfoType,textType));
  LEX_ATOM_PROP(custom, 13, cPType_int_as_string, offsetof(AtomInfoType,custom));
  LEX_ATOM_PROP(label, 14, cPType_int_as_string, offsetof(AtomInfoType,label));
  LEX_ATOM_PROP(numeric_type, 15, cPType_int_custom_type, offsetof(AtomInfoType,customType));
  LEX_ATOM_PROP(q, 16, cPType_float, offsetof(AtomInfoType,q));
  LEX_ATOM_PROP(b, 17, cPType_float, offsetof(AtomInfoType,b));
  LEX_ATOM_PROP(vdw, 18, cPType_float, offsetof(AtomInfoType,vdw));
  LEX_ATOM_PROP(elec_radius, 19, cPType_float, offsetof(AtomInfoType,elec_radius));
  LEX_ATOM_PROP(partial_charge, 20, cPType_float, offsetof(AtomInfoType,partialCharge));
  LEX_ATOM_PROP(formal_charge, 21, cPType_schar, offsetof(AtomInfoType,formalCharge));
  LEX_ATOM_PROP(stereo, 22, 0, 0);
  LEX_ATOM_PROP(cartoon, 23, cPType_schar, offsetof(AtomInfoType,cartoon));
  LEX_ATOM_PROP(color, 24, cPType_int, offsetof(AtomInfoType,color));
  LEX_ATOM_PROP(ID, 25, cPType_int, offsetof(AtomInfoType,id));
  LEX_ATOM_PROP(rank, 26, cPType_int, offsetof(AtomInfoType,rank));
  LEX_ATOM_PROP(flags, 27, cPType_uint32, offsetof(AtomInfoType,flags));
  LEX_ATOM_PROP(geom, 28, cPType_schar, offsetof(AtomInfoType,geom));
  LEX_ATOM_PROP(valence, 29, cPType_schar, offsetof(AtomInfoType,valence));
  LEX_ATOM_PROP(x, 30, cPType_xyz_float, 0);
  LEX_ATOM_PROP(y, 31, cPType_xyz_float, 1);
  LEX_ATOM_PROP(z, 32, cPType_xyz_float, 2);
  LEX_ATOM_PROP(settings, 33, cPType_settings, 0);
  LEX_ATOM_PROP(properties, 34, cPType_properties, 0);
  LEX_ATOM_PROP(s, 35, cPType_settings, 0);
  LEX_ATOM_PROP(p, 36, cPType_properties, 0);
  LEX_ATOM_PROP(state, 37, cPType_state, 0);
  LEX_ATOM_PROP(reps, 38, cPType_int, offsetof(AtomInfoType, visRep));
  LEX_ATOM_PROP(protons, 39, cPType_schar, offsetof(AtomInfoType, protons));
  LEX_ATOM_PROP(oneletter, 40, 0, 0);
  //  LEX_ATOM_PROP(, );

  return_OVstatus_SUCCESS;
}

int PyMOL_NewG3DStream(CPyMOL * I, int **array_ptr)
{
  int *return_vla = ExecutiveGetG3d(I->G);
  int result = OVstatus_FAILURE;
  if(return_vla) {
    result = VLAGetSize(return_vla) * (sizeof(G3dPrimitive) / sizeof(int));
  }
  if(array_ptr)
    *array_ptr = return_vla;
  return result;
}

int PyMOL_DelG3DStream(CPyMOL * I, int *array_ptr)
{
  VLAFreeP(array_ptr);
  return OVstatus_SUCCESS;
}

static OVstatus PyMOL_PurgeAPI(CPyMOL * I)
{
  OVOneToOne_DEL_AUTO_NULL(I->Setting);
  OVOneToOne_DEL_AUTO_NULL(I->Clip);
  OVOneToOne_DEL_AUTO_NULL(I->SelectList);
  OVOneToOne_DEL_AUTO_NULL(I->Reinit);
  OVOneToOne_DEL_AUTO_NULL(I->Rep);
#ifdef _PYMOL_LIB
  OVOneToOne_DEL_AUTO_NULL(I->MouseButtonCodeLexicon);
  OVOneToOne_DEL_AUTO_NULL(I->MouseButtonModCodeLexicon);
  OVOneToOne_DEL_AUTO_NULL(I->MouseButtonActionCodeLexicon);
  OVOneToOne_DEL_AUTO_NULL(I->MouseModeLexicon);
  OVOneToOne_DEL_AUTO_NULL(I->PaletteLexicon);
  OVOneToOne_DEL_AUTO_NULL(I->CartoonLexicon);
  OVOneToOne_DEL_AUTO_NULL(I->FlagLexicon);
  OVOneToOne_DEL_AUTO_NULL(I->FlagActionLexicon);
#endif

  OVOneToOne_DEL_AUTO_NULL(I->AtomPropertyLexicon);

  OVLexicon_DEL_AUTO_NULL(I->Lex);
  return_OVstatus_SUCCESS;
}

int PyMOL_FreeResultArray(CPyMOL * I, void *array)
{
  if(array) {
    VLAFreeP(array);
    return PyMOLstatus_SUCCESS;
  } else {
    return PyMOLstatus_FAILURE;
  }
}

PyMOLreturn_status PyMOL_CmdDraw(CPyMOL * I, int width, int height,
                                 int antialias, int quiet)
{
  PyMOLreturn_status result = { PyMOLstatus_FAILURE };
  PYMOL_API_LOCK
    result.status =
    get_status_ok(ExecutiveDrawCmd(I->G, width, height, antialias, false, quiet));
  I->ImageRequestedFlag = true;
  I->ImageReadyFlag = false;
  PYMOL_API_UNLOCK return result;
}

PyMOLreturn_status PyMOL_CmdCapture(CPyMOL * I, int quiet)
{
  PyMOLreturn_status result = { PyMOLstatus_FAILURE };
  PYMOL_API_LOCK
    result.status = get_status_ok(ExecutiveDrawCmd(I->G, -1, -1, 0, true, quiet));
  I->ImageRequestedFlag = true;
  I->ImageReadyFlag = false;
  PYMOL_API_UNLOCK return result;
}

PyMOLreturn_status PyMOL_CmdRay(CPyMOL * I, int width, int height, int antialias,
                                float angle, float shift, int renderer, int defer,
                                int quiet)
{
  PyMOLreturn_status result = { PyMOLstatus_FAILURE };
  PYMOL_API_LOCK if(renderer < 0)
    renderer = SettingGetGlobal_i(I->G, cSetting_ray_default_renderer);
  SceneInvalidateCopy(I->G, true);
  result.status =
    get_status_ok(ExecutiveRay
                  (I->G, width, height, renderer, angle, shift, quiet, defer, antialias));
  if(defer) {
    I->ImageRequestedFlag = true;
    I->ImageReadyFlag = false;
  } else {
    I->ImageRequestedFlag = false;
    if(SceneHasImage(I->G)) {
      I->ImageReadyFlag = true;
    } else {
      I->ImageReadyFlag = false;
    }
  }
  PYMOL_API_UNLOCK return result;
}

PyMOLreturn_status PyMOL_CmdSetView(CPyMOL * I, float *view, int view_len,
                                    float animate, int quiet)
{
  PyMOLreturn_status result = { PyMOLstatus_FAILURE };
  SceneViewType tmp;
  PYMOL_API_LOCK if(view_len >= 18) {
    int a;
    UtilZeroMem(tmp, sizeof(tmp));
    tmp[15] = 1.0F;
    for(a = 0; a < 3; a++) {
      tmp[a] = view[a];
      tmp[a + 4] = view[a + 3];
      tmp[a + 8] = view[a + 6];
      tmp[a + 16] = view[a + 9];
      tmp[a + 19] = view[a + 12];
      tmp[a + 22] = view[a + 15];
    }
    SceneSetView(I->G, tmp, quiet, animate, 0); /* TO DO -- add hand to the API */
    result.status = get_status_ok(true);
  } else {
    result.status = get_status_ok(false);
  }
  PYMOL_API_UNLOCK return result;
}

PyMOLreturn_float_array PyMOL_CmdGetView(CPyMOL * I, int quiet)
{
  PyMOLreturn_float_array result = { PyMOLstatus_FAILURE };
  SceneViewType tmp;
  PYMOL_API_LOCK result.size = 18;
  result.array = VLAlloc(float, result.size);
  if(result.array) {
    int a;
    SceneGetView(I->G, tmp);
    for(a = 0; a < 3; a++) {
      result.array[a] = tmp[a];
      result.array[a + 3] = tmp[a + 4];
      result.array[a + 6] = tmp[a + 8];
      result.array[a + 9] = tmp[a + 16];
      result.array[a + 12] = tmp[a + 19];
      result.array[a + 15] = tmp[a + 22];
    }
    result.status = get_status_ok(true);
  } else {
    result.status = get_status_ok(false);
  }
  PYMOL_API_UNLOCK return result;
}


PyMOLreturn_float_array PyMOL_CmdAlign(CPyMOL * I, const char *source, const char *target,
                                       float cutoff, int cycles, float gap, float extend,
                                       int max_gap, const char *object, const char *matrix,
                                       int source_state, int target_state, int quiet,
                                       int max_skip, int transform, int reset)
{
  PyMOLreturn_float_array result = { PyMOLstatus_FAILURE };
  PYMOL_API_LOCK OrthoLineType s2 = "", s3 = "";
  int ok = false;
  ExecutiveRMSInfo rms_info;
  result.size = 7;
  result.array = VLAlloc(float, result.size);
  if(!result.array) {
    ok = false;
  } else {
    ok = ((SelectorGetTmp(I->G, source, s2) >= 0) &&
          (SelectorGetTmp(I->G, target, s3) >= 0));
    if(ok) {
      const float _0 = 0.0F;    /* GCC compiler bug workaround */
      const float _m1 = -1.0F;
      ok = ExecutiveAlign(I->G, s2, s3, matrix, gap, extend, max_gap,
                          max_skip, cutoff, cycles, quiet, object,
                          source_state - 1, target_state - 1,
                          &rms_info, transform, reset, _m1, _0, _0, _0, _0, _0, 0, _0);
      if(ok) {
        result.array[0] = rms_info.final_rms;
        result.array[1] = rms_info.final_n_atom;
        result.array[2] = rms_info.n_cycles_run;
        result.array[3] = rms_info.initial_rms;
        result.array[4] = rms_info.initial_n_atom;
        result.array[5] = rms_info.raw_alignment_score;
        result.array[6] = rms_info.n_residues_aligned;
      }
    }
  }
  SelectorFreeTmp(I->G, s2);
  SelectorFreeTmp(I->G, s3);
  if(!ok) {
    VLAFreeP(result.array);
  }
  result.status = get_status_ok(ok);

  PYMOL_API_UNLOCK return result;

}

PyMOLreturn_status PyMOL_CmdDelete(CPyMOL * I, const char *name, int quiet)
{
  PYMOL_API_LOCK ExecutiveDelete(I->G, name);
  PyMOL_NeedRedisplay(I);  /* this should really only get called if ExecutiveDelete deletes something */
  PYMOL_API_UNLOCK return return_status_ok(true);       /* TO DO: return a real result */
}

PyMOLreturn_status PyMOL_CmdZoom(CPyMOL * I, const char *selection, float buffer,
                                 int state, int complete, float animate, int quiet)
{
  int ok = false;
  PYMOL_API_LOCK
    auto result = ExecutiveWindowZoom(I->G, selection, buffer, state - 1,
                             complete, animate, quiet);
    ok = static_cast<bool>(result);
  PYMOL_API_UNLOCK return return_status_ok(ok);
}

PyMOLreturn_status PyMOL_CmdOrient(CPyMOL * I, const char *selection, float buffer,
                                   int state, int complete, float animate, int quiet)
{
  int ok = true;
  PYMOL_API_LOCK
  auto res = ExecutiveOrient(
      I->G, selection, state - 1, animate, complete, buffer, quiet);
  ok = static_cast<bool>(res);
  PYMOL_API_UNLOCK return return_status_ok(ok);
}

PyMOLreturn_status PyMOL_CmdCenter(CPyMOL * I, const char *selection, int state, int origin,
                                   float animate, int quiet)
{
  int ok = false;
  PYMOL_API_LOCK
  auto result = ExecutiveCenter(I->G, selection,
      state - 1, origin, animate, nullptr, quiet);
  ok = static_cast<bool>(result);
  PYMOL_API_UNLOCK return return_status_ok(ok);
}

PyMOLreturn_status PyMOL_CmdOrigin(CPyMOL * I, const char *selection, int state, int quiet)
{
  int ok = true;
  PYMOL_API_LOCK
  float v[3] = { 0.0F, 0.0F, 0.0F };
  auto result = ExecutiveOrigin(I->G, selection, true, "", v, state - 1);
  ok = static_cast<bool>(result);
  PYMOL_API_UNLOCK return return_status_ok(ok);
}

PyMOLreturn_status PyMOL_CmdOriginAt(CPyMOL * I, float x, float y, float z, int quiet)
{
  int ok = true;
  PYMOL_API_LOCK float v[3];
  v[0] = x;
  v[1] = y;
  v[2] = z;
  auto result = ExecutiveOrigin(I->G, "", true, "", v, 0);
  ok = static_cast<bool>(ok);
  PYMOL_API_UNLOCK return return_status_ok(ok);
}

static OVreturn_word get_rep_id(CPyMOL * I, const char *representation)
{
  OVreturn_word result;

  if(!OVreturn_IS_OK((result = OVLexicon_BorrowFromCString(I->Lex, representation))))
    return result;
  return OVOneToOne_GetForward(I->Rep, result.word);
}

OVreturn_word get_setting_id(CPyMOL * I, const char *setting)
{
  OVreturn_word result;
  if(!OVreturn_IS_OK((result = OVLexicon_BorrowFromCString(I->Lex, setting))))
    return result;
  return OVOneToOne_GetForward(I->Setting, result.word);
}

static OVreturn_word get_reinit_id(CPyMOL * I, const char *reinit)
{
  OVreturn_word result;
  if(!OVreturn_IS_OK((result = OVLexicon_BorrowFromCString(I->Lex, reinit))))
    return result;
  return OVOneToOne_GetForward(I->Reinit, result.word);
}

static OVreturn_word get_select_list_mode(CPyMOL * I, const char *mode)
{
  OVreturn_word result;
  if(!OVreturn_IS_OK((result = OVLexicon_BorrowFromCString(I->Lex, mode))))
    return result;
  return OVOneToOne_GetForward(I->SelectList, result.word);
}

PyMOLreturn_status PyMOL_CmdClip(CPyMOL * I,
                                 const char *mode, float amount,
                                 const char *selection,
                                 int state, int quiet)
{
  int ok = true;
  PYMOL_API_LOCK
    SelectorTmp2 s1(I->G, selection);
    SceneClipFromMode(I->G, mode, amount, s1.getName(), state - 1);
  PYMOL_API_UNLOCK return return_status_ok(ok);
}

PyMOLreturn_status PyMOL_CmdLabel(CPyMOL * I, const char *selection, const char *text, int quiet)
{
  int ok = true;
  PYMOL_API_LOCK
  auto result = ExecutiveLabel(I->G, selection, text, quiet, cExecutiveLabelEvalAlt);
  ok = static_cast<bool>(result);
  PYMOL_API_UNLOCK return return_status_ok(ok);
}

PyMOLreturn_status PyMOL_CmdSelect(CPyMOL * I, const char *name, const char *selection, int quiet)
{
  int ret = -1;
  PYMOL_API_LOCK

  auto res = SelectorCreate(I->G, name, selection, NULL, quiet, NULL);
  ret = res ? res.result() : -1;

  PYMOL_API_UNLOCK return return_status_ok(ret >= 0); // if ret is negative it should fail
}

PyMOLreturn_status PyMOL_CmdSelectList(CPyMOL * I, const char *name, const char *object, int *list,
                                       int list_len, int state, const char *mode, int quiet)
{
  PyMOLreturn_status result = { PyMOLstatus_FAILURE };
  PYMOL_API_LOCK OVreturn_word mode_id;
  if(OVreturn_IS_OK((mode_id = get_select_list_mode(I, mode)))) {
    auto res =
      ExecutiveSelectList(I->G, name, object, list, list_len, state - 1, mode_id.word,
                          quiet);
    result = return_status_ok(bool(res));
  }
  PYMOL_API_UNLOCK return result;
}

PyMOLreturn_status PyMOL_CmdShow(CPyMOL * I,
                                 const char *representation,
                                 const char *selection,
                                 int quiet)
{
  int ok = true;
  PYMOL_API_LOCK OrthoLineType s1;
  OVreturn_word rep_id;
  if(OVreturn_IS_OK((rep_id = get_rep_id(I, representation)))) {
    SelectorGetTmp2(I->G, selection, s1);
    if (!s1[0]){  /* This doesn't catch patterns that don't match, but everything else */
      ok = false;
    } else {
      ExecutiveSetRepVisib(I->G, s1, rep_id.word, true);
      PyMOL_NeedRedisplay(I);  /* this should really only get called if ExecutiveSetRepVisib changes something */
      SelectorFreeTmp(I->G, s1);
    }
  } else {
    ok = false;
  }
  PYMOL_API_UNLOCK return return_status_ok(ok);
}

PyMOLreturn_status PyMOL_CmdHide(CPyMOL * I,
                                 const char *representation,
                                 const char *selection,
                                 int quiet)
{
  int ok = true;
  PYMOL_API_LOCK OrthoLineType s1;
  OVreturn_word rep_id;
  if(OVreturn_IS_OK((rep_id = get_rep_id(I, representation)))) {
    SelectorGetTmp2(I->G, selection, s1);
    if (!s1[0]){  /* This doesn't catch patterns that don't match, but everything else */
      ok = false;
    } else {
      ExecutiveSetRepVisib(I->G, s1, rep_id.word, false);
      SelectorFreeTmp(I->G, s1);
    }
  } else {
    ok = false;
  }
  PYMOL_API_UNLOCK return return_status_ok(ok);
}

PyMOLreturn_status PyMOL_CmdEnable(CPyMOL * I, const char *name, int quiet)
{
  int ok = false;
  PYMOL_API_LOCK if(name[0] == '(') {
    auto result1 = ExecutiveSetOnOffBySele(I->G, name, true);
    ok = static_cast<bool>(result1);
  } else {
    auto result2 =
        ExecutiveSetObjVisib(I->G, name, true, false); /* TO DO: parents */
    ok = static_cast<bool>(result2);
  }
  PYMOL_API_UNLOCK return return_status_ok(ok);
}

PyMOLreturn_status PyMOL_CmdDisable(CPyMOL * I, const char *name, int quiet)
{
  int ok = false;
  PYMOL_API_LOCK if(name[0] == '(') {
    auto result = ExecutiveSetOnOffBySele(I->G, name, false);
    ok = static_cast<bool>(result);
  } else {
    auto result = ExecutiveSetObjVisib(I->G, name, false, false);
    ok = static_cast<bool>(result);
  }
  PYMOL_API_UNLOCK return return_status_ok(ok);
}

PyMOLreturn_status PyMOL_CmdSetBond(CPyMOL * I, const char *setting, const char *value,
                                    const char *selection1, const char *selection2,
                                    int state, int quiet, int side_effects)
{
  int ok = true;
  PYMOL_API_LOCK {
    OVreturn_word setting_id;
    OrthoLineType s1 = "";
    OrthoLineType s2 = "";
    if(ok) ok = OVreturn_IS_OK((setting_id = get_setting_id(I, setting)));
    if(ok) ok = (SelectorGetTmp(I->G, selection1, s1) >= 0);
    if(ok) {
      if(selection2 && selection2[0]) {
        ok = (SelectorGetTmp(I->G, selection2, s2) >= 0);
      } else {
        ok = (SelectorGetTmp(I->G, selection1, s2) >= 0);
      }
    }
    if(ok) {
      ok = ExecutiveSetBondSettingFromString(I->G, setting_id.word, value,
                                             s1, s2,
                                             state - 1, quiet, side_effects);
    }
    SelectorFreeTmp(I->G, s1);
    SelectorFreeTmp(I->G, s2);
  } PYMOL_API_UNLOCK 
      return return_status_ok(ok);
}

PyMOLreturn_status PyMOL_CmdUnsetBond(CPyMOL * I, const char *setting,
                                      const char *selection1, const char *selection2,
                                      int state, int quiet, int side_effects)
{
  int ok = true;
  PYMOL_API_LOCK {
    OVreturn_word setting_id;
    OrthoLineType s1 = "";
    OrthoLineType s2 = "";
    if(ok) ok = OVreturn_IS_OK((setting_id = get_setting_id(I, setting)));
    if(ok) ok = (SelectorGetTmp(I->G, selection1, s1) >= 0);
    if(ok) {
      if(selection2 && selection2[0]) {
        ok = (SelectorGetTmp(I->G, selection2, s2) >= 0);
      } else {
        ok = (SelectorGetTmp(I->G, selection1, s2) >= 0);
      }
    }
    if(ok) {
      ok = ExecutiveUnsetBondSetting(I->G, setting_id.word, 
                                     s1, s2,
                                     state - 1, quiet, side_effects);
    }
    SelectorFreeTmp(I->G, s1);
    SelectorFreeTmp(I->G, s2);
  } PYMOL_API_UNLOCK 
      return return_status_ok(ok);
}

PyMOLreturn_status PyMOL_CmdSet(CPyMOL * I,
                                const char *setting,
                                const char *value,
                                const char *selection,
                                int state, int quiet, int side_effects)
{
  int ok = true;
  PYMOL_API_LOCK {
    OVreturn_word setting_id;
    OrthoLineType s1 = "";
    if(ok) ok = OVreturn_IS_OK((setting_id = get_setting_id(I, setting)));
    if(ok) ok = (SelectorGetTmp2(I->G, selection, s1) >= 0);
    
    if(ok) {
      ExecutiveSetSettingFromString(I->G, setting_id.word, value, s1,
                                    state - 1, quiet, side_effects);
    }
    SelectorFreeTmp(I->G, s1);
  }
  PYMOL_API_UNLOCK return return_status_ok(ok);
}

PyMOLreturn_value PyMOL_CmdGet(CPyMOL * I,
                                const char *setting,
                                const char *selection,
                                int state, int quiet){
  int ok = true;
  PyMOLreturn_value result = { PyMOLstatus_SUCCESS };

  PYMOL_API_LOCK {
    OVreturn_word setting_id;
    OrthoLineType s1 = "";
    if(ok) ok = OVreturn_IS_OK((setting_id = get_setting_id(I, setting)));
    if(ok) ok = (SelectorGetTmp2(I->G, selection, s1) >= 0);
    
    if(ok) {
      ExecutiveGetSettingFromString(I->G, &result, setting_id.word, s1,
                                    state - 1, quiet);
    }
    SelectorFreeTmp(I->G, s1);
  }
  PYMOL_API_UNLOCK return result;
}


PyMOLreturn_status PyMOL_CmdUnset(CPyMOL * I, const char *setting, const char *selection,
                                  int state, int quiet, int side_effects)
{
  int ok = true;
  PYMOL_API_LOCK {
    OVreturn_word setting_id;
    OrthoLineType s1 = "";
    if(ok) ok = OVreturn_IS_OK((setting_id = get_setting_id(I, setting)));
    if(ok) ok = (SelectorGetTmp2(I->G, selection, s1) >= 0);
    if(ok) {
      ExecutiveUnsetSetting(I->G, setting_id.word, s1,
                            state - 1, quiet, side_effects);
    }
    SelectorFreeTmp(I->G, s1);
  }
  PYMOL_API_UNLOCK return return_status_ok(ok);
}

PyMOLreturn_status PyMOL_CmdColor(CPyMOL * I, const char *color, const char *selection, int flags,
                                  int quiet)
{
  int ok = true;
  PYMOL_API_LOCK
  auto result = ExecutiveColorFromSele(I->G, selection, color, flags, quiet);
  ok = static_cast<bool>(result);
  PYMOL_API_UNLOCK return return_status_ok(ok);
}

/* -- JV */
PyMOLreturn_status PyMOL_CmdBackgroundColor(CPyMOL * I, const char *value) {
  int ok = true;
  PYMOL_API_LOCK

  int idx = ColorGetIndex(I->G, value);
  if(idx >= 0){
    SettingSetGlobal_i(I->G, cSetting_bg_rgb, idx);
  } else {
    ErrMessage(I->G, "Color", "Bad color name.");
  }
  PYMOL_API_UNLOCK return return_status_ok(ok);
}

PyMOLreturn_status PyMOL_CmdReinitialize(CPyMOL * I,
    const char *what,
    const char *object_name)
{
  int ok = true;
  OVreturn_word what_id;
  PYMOL_API_LOCK if(OVreturn_IS_OK((what_id = get_reinit_id(I, what)))) {
    auto res = ExecutiveReinitialize(I->G, what_id.word, object_name);
    ok = static_cast<bool>(res);
  }
  PYMOL_API_UNLOCK return return_status_ok(ok);
}

PyMOLreturn_int PyMOL_CmdGetMovieLength(CPyMOL * I,int quiet)
{
  int ok = true;
  PyMOLreturn_int result;
  result.status = PyMOLstatus_FAILURE;
  result.value = 0;

  PYMOL_API_LOCK
  if(ok) {
    result.value = MovieGetLength(I->G);
    result.status = get_status_ok(ok);
  };
  PYMOL_API_UNLOCK return result;
}

PyMOLreturn_float PyMOL_CmdGetDistance(CPyMOL * I,
                                       const char *selection1,
                                       const char *selection2, int state, int quiet)
{
  PyMOLreturn_float result;
  PYMOL_API_LOCK {
    result = return_result(ExecutiveGetDistance(I->G,
        selection1,
        selection2,
        state));
  }
  PYMOL_API_UNLOCK return result;
}

PyMOLreturn_float PyMOL_CmdDistance(CPyMOL * I,
                                    const char *name,
                                    const char *selection1,
                                    const char *selection2,
                                    int mode,
                                    float cutoff,
                                    int label, int reset, int zoom, int state, int quiet)
{
  PyMOLreturn_float result;
  PYMOL_API_LOCK {
    int defState1 = -4, defState2 = -4;
    auto res = ExecutiveDistance(I->G, name,
        selection1, selection2, mode, cutoff, label, quiet, reset, state, zoom, defState1, defState2);
    result = return_result(res);
  }
  PYMOL_API_UNLOCK return result;
}

PyMOLreturn_float PyMOL_CmdGetAngle(CPyMOL * I,
                                    const char *selection1,
                                    const char *selection2,
                                    const char *selection3, int state, int quiet)
{
  PyMOLreturn_float result;
  PYMOL_API_LOCK {
    result = return_result(ExecutiveGetAngle(I->G,
        selection1,
        selection2,
        selection3,
        state));
  }
  PYMOL_API_UNLOCK return result;
}

PyMOLreturn_float PyMOL_CmdAngle(CPyMOL * I,
                                 const char *name,
                                 const char *selection1,
                                 const char *selection2,
                                 const char *selection3,
                                 int mode,
                                 int label, int reset, int zoom, int state, int quiet)
{
  PyMOLreturn_float result;
  PYMOL_API_LOCK {
    int defState1 = -4, defState2 = -4, defState3 = -3;
    auto res = ExecutiveAngle(I->G, name,
        selection1, selection2, selection3, mode, label, reset, zoom, quiet,
        state, defState1, defState2, defState3);
    result = return_result(res);
  }
  PYMOL_API_UNLOCK return result;
}

PyMOLreturn_float PyMOL_CmdGetDihedral(CPyMOL * I,
                                       const char *selection1,
                                       const char *selection2,
                                       const char *selection3,
                                       const char *selection4, int state, int quiet)
{
  PyMOLreturn_float result;
  PYMOL_API_LOCK {
    result = return_result(ExecutiveGetDihe(I->G,
        selection1,
        selection2,
        selection3,
        selection4,
        state));
  }
  PYMOL_API_UNLOCK return result;
}

PyMOLreturn_float PyMOL_CmdDihedral(CPyMOL * I,
                                    const char *name,
                                    const char *selection1,
                                    const char *selection2,
                                    const char *selection3,
                                    const char *selection4,
                                    int mode,
                                    int label, int reset, int zoom, int state, int quiet)
{
  PyMOLreturn_float result;
  PYMOL_API_LOCK {
    auto res = ExecutiveDihedral(I->G,
        name, selection1, selection2, selection3, selection4, mode, label,
        reset, zoom, quiet, state);
    result = return_result(res);
  }
  PYMOL_API_UNLOCK return result;
}

PyMOLreturn_status PyMOL_CmdIsodot(CPyMOL * I, const char *name, const char *map_name, float level,
                                   const char *selection, float buffer, int state, float carve,
                                   int source_state, int quiet)
{
  PyMOLreturn_status result = { PyMOLstatus_FAILURE };
  PYMOL_API_LOCK
  auto res = ExecutiveIsomeshEtc(I->G, name, map_name, level, selection, buffer,
      state - 1, carve, source_state - 1, quiet, 1, level);
  result.status = get_status_ok(bool(res));
  PYMOL_API_UNLOCK return result;

}

PyMOLreturn_status PyMOL_CmdIsomesh(CPyMOL * I, const char *name, const char *map_name, float level,
                                    const char *selection, float buffer, int state, float carve,
                                    int source_state, int quiet)
{
  PyMOLreturn_status result = { PyMOLstatus_FAILURE };
  PYMOL_API_LOCK
  auto res = ExecutiveIsomeshEtc(I->G, name, map_name, level, selection, buffer,
      state - 1, carve, source_state - 1, quiet, 0, level);
  result.status = get_status_ok(bool(res));
  PYMOL_API_UNLOCK return result;
}

PyMOLreturn_status PyMOL_CmdIsosurface(CPyMOL * I, const char *name, const char *map_name,
                                       float level, const char *selection, float buffer,
                                       int state, float carve, int source_state, int side,
                                       int mode, int quiet)
{
  PyMOLreturn_status result = { PyMOLstatus_FAILURE };
  PYMOL_API_LOCK
  auto res = ExecutiveIsosurfaceEtc(I->G, name, map_name, level, selection,
      buffer, state - 1, carve, source_state - 1, side, quiet, mode);
  result.status = get_status_ok(bool(res));
  PYMOL_API_UNLOCK return result;
}

PyMOLreturn_status PyMOL_CmdGradient(CPyMOL * I, const char *name, const char *map_name,
                                     float minimum, float maximum, const char *selection,
                                     float buffer, int state, float carve,
                                     int source_state, int quiet)
{
  PyMOLreturn_status result = { PyMOLstatus_FAILURE };
  PYMOL_API_LOCK
  auto res = ExecutiveIsomeshEtc(I->G, name, map_name, minimum, selection,
      buffer, state - 1, carve, source_state - 1, quiet, 3, maximum);
  result.status = get_status_ok(bool(res));
  PYMOL_API_UNLOCK return result;
}

PyMOLreturn_float PyMOL_CmdIsolevel(CPyMOL * I, const char *name, float level, int state,
                                    int query, int quiet)
{
  PyMOLreturn_float result;
  PYMOL_API_LOCK
    if(query) {
      auto res = ExecutiveGetIsolevel(I->G, name, state - 1);
      result = return_result(res);
  } else {
      auto res = ExecutiveIsolevel(I->G, name, level, state - 1, quiet);
      result.status = get_status_ok(static_cast<bool>(res));
      result.value = level;
  }
  PYMOL_API_UNLOCK return result;
}

static int word_count(const char *src)
{                               /* only works for ascii */
  int cnt = 0;
  while((*src) && ((*src) < 33))        /* skip leading whitespace */
    src++;
  while(*src) {
    if((*src) > 32) {
      cnt++;
      while((*src) && ((*src) > 32))
        src++;
    }
    while((*src) && ((*src) < 33))
      src++;
  }
  return cnt;
}

static const char *next_word(const char *src, char *dst, int buf_size)
{                               /* only works for ascii */
  while((*src) && ((*src) < 33))        /* skip leading whitespace */
    src++;
  while(*src) {
    if((*src) > 32) {
      while((*src) && ((*src) > 32) && (buf_size > 1)) {
        *(dst++) = *(src++);
        buf_size--;
      }
      break;
    }
  }
  dst[0] = 0;
  return src;
}

PyMOLreturn_status PyMOL_CmdRampNew(CPyMOL * I, const char *name, const char *map, float *range,
                                    int n_level, const char *color, int state, const char *selection,
                                    float beyond, float within, float sigma,
                                    int zero, int calc_mode, int quiet)
{
  int ok = true;
  PyMOLreturn_status result = { PyMOLstatus_FAILURE };
  OrthoLineType s1 = "";
  float *color_vla = NULL;
  float *range_vla = NULL;
  PYMOL_API_LOCK if(selection && selection[0]) {
    if(ok)
      ok = (SelectorGetTmp(I->G, selection, s1) >= 0);
  }
  if(ok) {
    if(range && n_level) {
      range_vla = VLAlloc(float, n_level);
      UtilCopyMem(range_vla, range, sizeof(float) * n_level);
    }
  }

  if(ok && color) {
    int n_color = word_count(color);
    /* to do */
    if(color && n_color) {
      color_vla = VLAlloc(float, n_color * 3);
      if(color_vla) {
        WordType colorName;
        int a;
        for(a = 0; a < n_color; a++) {
          color = next_word(color, colorName, sizeof(colorName));
          {
            const float *src = ColorGetNamed(I->G, colorName);
            float *dst = color_vla + 3 * a;
            copy3f(src, dst);
          }
        }
      }
    }
  }
  if(ok) {
    auto res =
        ExecutiveRampNew(I->G, name, map, pymol::vla_take_ownership(range_vla),
            pymol::vla_take_ownership(color_vla), state, s1, beyond, within,
            sigma, zero, calc_mode, quiet);
    ok = static_cast<bool>(res);
    result.status = get_status_ok(ok);
  } else {
    result.status = PyMOLstatus_FAILURE;
  }
  SelectorFreeTmp(I->G, s1);
  PYMOL_API_UNLOCK return result;
}

/**
 * Supported file formats and their internal codes
 */
struct {
  const char * name;
  cLoadType_t code_buffer;
  cLoadType_t code_filename;
} const ContentTypeTable[] = {
  // molecules
  {"pdb",           cLoadTypePDBStr,    cLoadTypePDB},
  {"vdb",           cLoadTypeVDBStr,    cLoadTypeUnknown},
  {"cif",           cLoadTypeCIFStr,    cLoadTypeCIF},
  {"mmtf",          cLoadTypeMMTFStr,   cLoadTypeMMTF},
  {"mae",           cLoadTypeMAEStr,    cLoadTypeMAE},
  {"sdf",           cLoadTypeSDF2Str,   cLoadTypeSDF2},
  {"mol",           cLoadTypeMOLStr,    cLoadTypeMOL},
  {"mol2",          cLoadTypeMOL2Str,   cLoadTypeMOL2},
  {"xyz",           cLoadTypeXYZStr,    cLoadTypeXYZ},
  {"pqr",           cLoadTypeUnknown,   cLoadTypePQR},
  {"macromodel",    cLoadTypeMMDStr,    cLoadTypeMMD},
  // maps
  {"ccp4",          cLoadTypeCCP4Str,   cLoadTypeCCP4Map},
  {"mrc",           cLoadTypeMRCStr,    cLoadTypeMRC},
  {"map",           cLoadTypeCCP4UnspecifiedStr, cLoadTypeCCP4Unspecified},
  {"xplor",         cLoadTypeXPLORStr,  cLoadTypeXPLORMap},
  {"phi",           cLoadTypePHIStr,    cLoadTypePHIMap},
  {"dx",            cLoadTypeDXStr,     cLoadTypeDXMap},
  // special
  {"cgo",           cLoadTypeCGO,       cLoadTypeUnknown},
  {NULL,            cLoadTypeUnknown,   cLoadTypeUnknown}
};

/**
 * Proxy for "ExecutiveLoad" with string "content_format" (and "content_type")
 * argument.
 *
 * content:     Either file name or file contents, depending on "content_type"
 * content_type:        "filename", "string", "raw", or "cgo"
 * content_length:      Length of "content", if it's not a file name or a
 *                      null-terminated string (pass -1).
 * content_format:      The file format, e.g. "pdb", "sdf", "mol2", ...
 * object_name:         New object name. Can be empty if "content_type" is
 *                      "filename".
 */
static PyMOLreturn_status Loader(CPyMOL * I, const char *content, const char *content_type,
                                 int content_length, const char *content_format,
                                 const char *object_name, int state,
                                 int discrete, int finish,
                                 int quiet, int multiplex, int zoom)
{
  PyMOLGlobals * G = I->G;
  bool content_is_filename = false;
  int ok = true;
  WordType obj_name;

  // `content` can be a file name, or the file contents
  if (strcmp(content_type, "filename") == 0) {
    content_is_filename = true;
  } else if (strcmp(content_type, "string") == 0) {
    if (content_length < 0)
      content_length = strlen(content);
  } else if (
      strcmp(content_type, "raw") != 0 &&
      strcmp(content_type, "cgo") != 0) {
    PRINTFB(G, FB_Executive, FB_Errors)
      " Error: Unknown content type '%s'\n", content_type ENDFB(G);
    ok = false;
  }

  if(ok) {
    {                           /* if object_name is blank and content is a filename, then 
                                   compute the object_name from the file prefix */
      if((!object_name[0]) && content_is_filename) {
        const char *start, *stop;
        stop = start = content + strlen(content) - 1;
        while(start > content) {        /* known path separators */
          if((start[-1] == ':') || (start[-1] == '\'') || (start[-1] == '/'))
            break;
          start--;
        }
        while(stop > start) {
          if(*stop == '.')
            break;
          stop--;
        }
        if(stop == start)
          stop = content + strlen(content);
        if((stop - start) >= sizeof(WordType))
          stop = start + sizeof(WordType) - 1;
        {
          char *q;
          const char *p = start;
          q = obj_name;
          while(p < stop) {
            *(q++) = *(p++);
          }
          *q = 0;
          object_name = obj_name;
        }
      }
    }
    {
      cLoadType_t pymol_content_type = cLoadTypeUnknown;

      for (auto it = ContentTypeTable; it->name; ++it) {
        if (strcmp(it->name, content_format) == 0) {
          pymol_content_type = content_is_filename ?
            it->code_filename : it->code_buffer;
          break;
        }
      }

      if (pymol_content_type == cLoadTypeUnknown) {
        PRINTFB(G, FB_Executive, FB_Errors)
          " Error: Unknown content format '%s' with type '%s'\n",
          content_format, content_type ENDFB(G);
        ok = false;
      }

      if(ok) {
        auto result = ExecutiveLoad(I->G,
            content_is_filename ? content : nullptr,
            content_is_filename ? nullptr : content,
                           content_length,
                           pymol_content_type,
                           object_name,
                           state - 1, zoom, discrete, finish, multiplex, quiet, NULL, 0, NULL);
        ok = static_cast<bool>(result);
      }
    }
  }
  if (ok)
    PyMOL_NeedRedisplay(I);
  return return_status_ok(ok);
}

PyMOLreturn_status PyMOL_CmdLoad(CPyMOL * I,
                                 const char *content,
                                 const char *content_type,
                                 const char *content_format,
                                 const char *object_name, int state,
                                 int discrete, int finish,
                                 int quiet, int multiplex, int zoom)
{
  PyMOLreturn_status status = { PyMOLstatus_FAILURE };
  PYMOL_API_LOCK
    status = Loader(I, content, content_type, -1, content_format, object_name,
                    state, discrete, finish, quiet, multiplex, zoom);
  PYMOL_API_UNLOCK return status;
}

PyMOLreturn_status PyMOL_CmdLoadRaw(CPyMOL * I,
                                    const char *content,
                                    int content_length,
                                    const char *content_format,
                                    const char *object_name, int state,
                                    int discrete, int finish,
                                    int quiet, int multiplex, int zoom)
{
  PyMOLreturn_status status = { PyMOLstatus_FAILURE };
  PYMOL_API_LOCK
    status = Loader(I, content, "raw", content_length, content_format,
                    object_name, state, discrete, finish, quiet, multiplex, zoom);
  PYMOL_API_UNLOCK return status;
}

PyMOLreturn_status PyMOL_CmdLoadCGO(CPyMOL * I,
                                    const float *content,
                                    int content_length,
                                    const char *object_name, int state, int quiet, int zoom)
{
  PyMOLreturn_status status = { PyMOLstatus_FAILURE };
  PYMOL_API_LOCK
    status = Loader(I, (char *) content, "cgo", content_length, "cgo",
                    object_name, state, 0, 1, quiet, 0, zoom);
  PYMOL_API_UNLOCK return status;
}

PyMOLreturn_status PyMOL_CmdCreate(CPyMOL * I,  
                                   const char *name,
                                   const char *selection, int source_state,
                                   int target_state, int discrete,
                                   int zoom, int quiet, int singletons,
                                   const char *extract, int copy_properties)
{
  int ok = true;
  PYMOL_API_LOCK
  auto result = ExecutiveSeleToObject(I->G,
      name, selection, source_state, target_state, discrete, zoom, quiet,
      singletons, copy_properties);
  ok = static_cast<bool>(result);
  PYMOL_API_UNLOCK return return_status_ok(ok);
}

PyMOLreturn_status PyMOL_CmdPseudoatom(CPyMOL * I, const char *object_name, const char *selection,
				       const char *name, const char *resn, const char *resi, const char *chain,
				       const char *segi, const char *elem, float vdw, int hetatm,
				       float b, float q, const char *color, const char *label, 
				       int use_xyz, float x, float y, float z,
				       int state, int mode, int quiet)
{
  int ok = true;
  PYMOL_API_LOCK
  if(ok) {
    int color_index = ColorGetIndex(I->G, color);
    if(ok) {
      float pos_tmp[3], *pos = pos_tmp;
      if(use_xyz) {
	pos[0] = x;
	pos[1] = y;
	pos[2] = z;
      } else {
	pos = NULL;
      }
      auto pseudoatom_name = ExecutivePreparePseudoatomName(I->G, object_name);
      auto res = ExecutivePseudoatom(I->G, pseudoatom_name, selection, name,
          resn, resi, chain, segi, elem, vdw, hetatm, b, q, label, pos,
          color_index, state - 1, mode, quiet);
      ok = static_cast<bool>(res);
    }
  }
  PYMOL_API_UNLOCK return return_status_ok(ok);
}

PyMOLreturn_status PyMOL_CmdTurn(CPyMOL * I, char axis, float angle){
  PyMOLreturn_status result = { PyMOLstatus_SUCCESS };
  PYMOL_API_LOCK PyMOLGlobals * G = I->G;
  switch (axis){
  case 'x':
    SceneRotate(G, angle, 1.0, 0.0, 0.0);
    break;
  case 'y':
    SceneRotate(G, angle, 0.0, 1.0, 0.0);
    break;
  case 'z':
    SceneRotate(G, angle, 0.0, 0.0, 1.0);
    break;
  default:
    result.status = PyMOLstatus_FAILURE;
    break;
  }
  PYMOL_API_UNLOCK return result;
}

PyMOLreturn_status PyMOL_CmdMPlay(CPyMOL * I, int cmd){
  PyMOLreturn_status result = { PyMOLstatus_SUCCESS };
  PYMOL_API_LOCK PyMOLGlobals * G = I->G;
  MoviePlay(G, cmd);
  PYMOL_API_UNLOCK return result;
}

PyMOLreturn_status PyMOL_CmdSetFeedbackMask(CPyMOL * I, int action, int module, int mask){
  PyMOLreturn_status result = { PyMOLstatus_SUCCESS };
  PYMOL_API_LOCK PyMOLGlobals * G = I->G;
  switch (action){
  case 0:
    G->Feedback->setMask(module, (uchar) mask);
    break;
  case 1:
    G->Feedback->enable(module, (uchar) mask);
    break;
  case 2:
    G->Feedback->disable(module, (uchar) mask);
    break;
  case 3:
    G->Feedback->push();
    break;
  case 4:
    G->Feedback->pop();
    break;
  }
  PYMOL_API_UNLOCK return result;
}

static const CPyMOLOptions Defaults = {
  true,                         /* pmgui */
#ifndef _PYMOL_NOPY
  true,                         /* internal_gui */
#else
  false,
#endif
#ifndef _PYMOL_NOPY
  true,                         /* show_splash */
#else
  false,
#endif
#ifndef _PYMOL_NOPY
  1,                            /* internal_feedback */
#else
  0,
#endif
  true,                         /* security */
  false,                        /* game mode */
  0,                            /* force_stereo */
  640,                          /* winX */
  480,                          /* winY */
  false,                        /* blue_line */
  0,                            /* winPX */
  175,                          /* winPY */
  true,                         /* external_gui */
  true,                         /* siginthand */
  false,                        /* reuse helper */
  false,                        /* auto reinitialize */
  false,                        /* keep thread alive */
  false,                        /* quiet */
  false,                        /* incentive product */
  "",                           /* after_load_script */
  0,                            /* multisample */
  1,                            /* window_visible */
  0,                            /* read_stdin */
  0,                            /* presentation */
  0,                            /* defer builds mode */
  0,                            /* full screen mode */
  -1,                           /* sphere mode */
  0,                            /* stereo capable */
  0,                            /* stereo mode */
  -1,                           /* zoom mode */
  0,                            /* launch_status */
  0,                            /* no quit */
  0,                            /* gldebug */
  false,                        /* no openvr stub */
};

CPyMOLOptions *PyMOLOptions_New(void)
{
  CPyMOLOptions *result = NULL;
  result = pymol::calloc<CPyMOLOptions>(1);
  if(result)
    *result = Defaults;
  return result;
}

void PyMOLOptions_Free(CPyMOLOptions * options)
{
  FreeP(options);
}

void PyMOL_ResetProgress(CPyMOL * I)
{
  I->ProgressChanged = true;
  UtilZeroMem(I->Progress, sizeof(int) * 6);
}

void PyMOL_SetProgress(CPyMOL * I, int offset, int current, int range)
{
  switch (offset) {
  case PYMOL_PROGRESS_SLOW:
  case PYMOL_PROGRESS_MED:
  case PYMOL_PROGRESS_FAST:
    if(current != I->Progress[offset]) {
      I->Progress[offset] = current;
      I->ProgressChanged = true;
    }
    if(range != I->Progress[offset + 1]) {
      I->Progress[offset + 1] = range;
      I->ProgressChanged = true;
    }
  }
}

int PyMOL_GetProgress(CPyMOL * I, int *progress, int reset)
{
  int a;
  int result = I->ProgressChanged;
  for(a = 0; a < PYMOL_PROGRESS_SIZE; a++) {
    progress[a] = I->Progress[a];
  }
  if(reset)
    I->ProgressChanged = false;
  return result;
}

int PyMOL_GetProgressChanged(CPyMOL * I, int reset)
{
  int result = I->ProgressChanged;
  if(reset)
    I->ProgressChanged = false;
  return result;
}

CPyMOL *PyMOL_New(void)
{
  assert(!Defaults.stereo_capable);
  return PyMOL_NewWithOptions(&Defaults);
}

CPyMOL *PyMOL_NewWithOptions(const CPyMOLOptions * option)
{
  auto result = pymol::calloc<CPyMOL>(1);
  assert(result);

  auto G = pymol::calloc<PyMOLGlobals>(1);
  assert(G);

  result->G = G;
  G->PyMOL = result;

  PyMOL_ResetProgress(result);

  G->Option = pymol::calloc<CPyMOLOptions>(1);
  assert(G->Option);

  if (!option) {
    assert(!Defaults.stereo_capable);
    option = &Defaults;
  }

  *(G->Option) = *option;

  G->Security = option->security;
  G->StereoCapable = option->stereo_capable;

  return result;
}

void PyMOL_Start(CPyMOL * I)
{
  PyMOLGlobals *G = I->G;

  // It's possible to change this from Python, functions which rely on
  // C locale should reset it before doing printf, atof, etc.
  std::setlocale(LC_NUMERIC, "C");

  G->Context = OVContext_New();
  G->Lexicon = OVLexicon_New(G->Context->heap);

  if(OVreturn_IS_ERROR(PyMOL_InitAPI(I))) {
    printf("ERROR: PyMOL internal C API initialization failed.\n");
  }

  // global lexicon "constants"
#define LEX_CONSTANTS_IMPL
#include "lex_constants.h"

  G->Feedback = new CFeedback(G, G->Option->quiet);
  WordInit(G);
  UtilInit(G);
  ColorInit(G);
  CGORendererInit(G);
  ShaderMgrInit(G);
  SettingInitGlobal(G, true, true, false);
  SettingSetGlobal_i(G, cSetting_internal_gui, G->Option->internal_gui);
  SettingSetGlobal_i(G, cSetting_internal_feedback, G->Option->internal_feedback);
  TextureInit(G);
  TypeInit(G);
  TextInit(G);
  CharacterInit(G);
  PlugIOManagerInit(G);
  SphereInit(G);
  // OpenVRInit() called in ExecutiveStereo
  OrthoInit(G, G->Option->show_splash);
  SceneInit(G);
  MovieScenesInit(G);
  WizardInit(G);                /* must come after ortho & scene */
  G->Movie = new CMovie(G);
  G->SelectorMgr = new CSelectorManager();
  G->Selector = new CSelector(G, G->SelectorMgr);
  SeqInit(G);
  SeekerInit(G);
  ButModeInit(G);
  ControlInit(G);
  AtomInfoInit(G);
  SculptCacheInit(G);
  VFontInit(G);
  ExecutiveInit(G);
  IsosurfInit(G);
  TetsurfInit(G);
  EditorInit(G);
#ifdef TRACKER_UNIT_TEST
  TrackerUnitTest(G);
#endif

  I->DrawnFlag = false;
  I->RedisplayFlag = true;
  G->Ready = true;
}

/* This function is necessary to be called from PyMOL_StartWithPython
   which is called before the PYMOL_API is instantiated, thus
   it is not necessary (and you can't) lock the API */
void PyMOL_ConfigureShadersGL_WithoutLock(CPyMOL * I){
    I->done_ConfigureShaders = false;
    // ShaderMgr->Config() moved to PyMOL_DrawWithoutLock
}

/* This function is called from CMol and needs to lock 
   the PYMOL_API */
void PyMOL_ConfigureShadersGL(CPyMOL * I){
  PYMOL_API_LOCK 
    PyMOL_ConfigureShadersGL_WithoutLock(I);
  PYMOL_API_UNLOCK
}

#ifndef _PYMOL_NOPY

void PyMOL_StartWithPython(CPyMOL * I)
{
  PyMOL_Start(I);

  /* now locate all the C to Python function hooks and objects we need */

  PInit(I->G, false);

  /* and begin the initialization sequence */

  I->PythonInitStage = 1;
}

#endif

void PyMOL_Stop(CPyMOL * I)
{
  PyMOLGlobals *G = I->G;
  G->Terminating = true;
  TetsurfFree(G);
  IsosurfFree(G);
  WizardFree(G);
  EditorFree(G);
  ExecutiveFree(G);
  VFontFree(G);
  SculptCacheFree(G);
  AtomInfoFree(G);
  ButModeFree(G);
  ControlFree(G);
  SeekerFree(G);
  SeqFree(G);
  DeleteP(G->Selector);
  DeleteP(G->SelectorMgr);
  DeleteP(G->Movie);
  SceneFree(G);
  MovieScenesFree(G);
  OrthoFree(G);
#ifdef _PYMOL_OPENVR
  OpenVRFree(G);
#endif
  DeleteP(G->ShaderMgr);
  SettingFreeGlobal(G);
  CharacterFree(G);
  TextFree(G);
  TypeFree(G);
  TextureFree(G);
  SphereFree(G);
  PlugIOManagerFree(G);
  PFree(G);
  CGORendererFree(G);
  ColorFree(G);
  UtilFree(G);
  WordFree(G);
  DeleteP(G->Feedback);

  PyMOL_PurgeAPI(I);
  /*    printf("%d \n", OVLexicon_GetNActive(G->Lexicon)); */
  OVLexicon_Del(G->Lexicon);
  OVContext_Del(G->Context);
}

void PyMOL_Free(CPyMOL * I)
{
#if !defined(_PYMOL_ACTIVEX) && !defined(_MACPYMOL_XCODE)
  PYMOL_API_LOCK
#endif
    /* take PyMOL down gracefully */
    PyMOLOptions_Free(I->G->Option);

#ifndef _PYMOL_NOPY
  FreeP(I->G->P_inst);
  if(I->G == SingletonPyMOLGlobals)
    SingletonPyMOLGlobals = NULL;
#endif

  FreeP(I->G);
  FreeP(I);
  return;
#if !defined(_PYMOL_ACTIVEX) && !defined(_MACPYMOL_XCODE)
  PYMOL_API_UNLOCK;
#endif
}

struct PyMOLGlobals* PyMOL_GetGlobals(CPyMOL * I)
{
  return I->G;
}

struct PyMOLGlobals** PyMOL_GetGlobalsHandle(CPyMOL * I)
{
  return &(I->G);
}

void PyMOL_LockAPIAndUnblock(CPyMOL * I)
{
  PyMOLGlobals *G = I->G;
  (void)G;
  PLockAPIAndUnblock(G);
}

void PyMOL_BlockAndUnlockAPI(CPyMOL * I)
{
  PyMOLGlobals *G = I->G;
  (void)G;
  PBlockAndUnlockAPI(G);
}

void PyMOL_AdaptToHardware(CPyMOL * I)
{
  PYMOL_API_LOCK PyMOLGlobals * G = I->G;
  if(G->HaveGUI) {
    PyMOL_PushValidContext(I);
    {
      char *vendor = (char *) glGetString(GL_VENDOR);
      char *renderer = (char *) glGetString(GL_RENDERER);
      char *version = (char *) glGetString(GL_VERSION);
      if(vendor && version) {
        /* work around broken lighting under Windows GDI Generic */
        if((strcmp(vendor, "Microsoft Corporation") == 0) &&
           (strcmp(renderer, "GDI Generic") == 0)) {
          ExecutiveSetSettingFromString(I->G, cSetting_light_count, "1", "", 0, 1, 0);
          ExecutiveSetSettingFromString(I->G, cSetting_spec_direct, "0.7", "", 0, 1, 0);
        }
      }
    }
    PyMOL_PopValidContext(I);
  }
PYMOL_API_UNLOCK}

static void setup_gl_state(void)
{

  /* get us into a well defined GL state */

  /*glMatrixMode(GL_PROJECTION);
     glLoadIdentity();
     glMatrixMode(GL_MODELVIEW);
     glLoadIdentity(); */

#ifndef PURE_OPENGL_ES_2
  glDisable(GL_ALPHA_TEST);
  glDisable(GL_COLOR_LOGIC_OP);
  glDisable(GL_COLOR_MATERIAL);
  glDisable(GL_FOG);
  glDisable(GL_LIGHTING);
  glDisable(GL_LIGHT0);
  glDisable(GL_LIGHT1);
  glDisable(GL_LINE_SMOOTH);
  glDisable(GL_NORMALIZE);
#endif

  glDisable(GL_BLEND);
  glDisable(GL_CULL_FACE);
  glDisable(GL_DEPTH_TEST);
  glDisable(GL_DITHER);
#ifndef PURE_OPENGL_ES_2
  glDisable(GL_POLYGON_SMOOTH);
#endif
}
void PyMOL_DrawWithoutLock(CPyMOL * I);

void PyMOL_Draw(CPyMOL * I){
  PYMOL_API_LOCK_MODAL
  PyMOL_DrawWithoutLock(I);
  PYMOL_API_UNLOCK
}

static void PyMOL_LaunchStatus_Feedback(PyMOLGlobals * G)
{
  G->LaunchStatus |= G->Option->launch_status;

  if(G->StereoCapable) {
    OrthoAddOutput(G,
        " OpenGL quad-buffer stereo 3D detected and enabled.\n");;
  } else {
    if(G->LaunchStatus & cPyMOLGlobals_LaunchStatus_StereoFailed) {
      G->Feedback->addColored(
          "Error: The requested stereo 3D visualization mode is not available.\n",
          FB_Errors);
    }
  }

  if(G->LaunchStatus & cPyMOLGlobals_LaunchStatus_MultisampleFailed) {
    G->Feedback->addColored(
        "Error: The requested multisampling mode is not available.\n",
        FB_Errors);
  }
}

#ifndef PURE_OPENGL_ES_2
static void check_gl_stereo_capable(PyMOLGlobals * G)
{
  // quad buffer stereo available?
  GLboolean state;
  glGetBooleanv(GL_STEREO, &state);
  G->StereoCapable = state || G->Option->force_stereo > 0;
  if (!state && G->Option->force_stereo > 0) {
    printf("Warning: forcing stereo despite GL_STEREO=0\n");
  }

  // stereo request feedback
  if (state && G->Option->stereo_mode == cStereo_default) {
    SettingSetGlobal_i(G, cSetting_stereo_mode, cStereo_quadbuffer);
  } else if (!state && G->Option->stereo_mode == cStereo_quadbuffer) {
    G->LaunchStatus |= cPyMOLGlobals_LaunchStatus_StereoFailed;
  }

  // multisample request feedback
  if (G->Option->multisample) {
    GLint samplebuffers = 0;
    glGetIntegerv(GL_SAMPLE_BUFFERS, &samplebuffers);
    if (!samplebuffers) {
      G->LaunchStatus |= cPyMOLGlobals_LaunchStatus_MultisampleFailed;
    }
  }

  // GL_BACK if GL_DOUBLEBUFFER else GL_FRONT
  // With QOpenGLWidget -> framebuffer object
  GLint buf;
  glGetIntegerv(GL_DRAW_BUFFER0, &buf);
  if (!buf) {
    printf("Warning: GL_DRAW_BUFFER0=0 -> using GL_BACK\n");
    buf = GL_BACK;
  }
  G->DRAW_BUFFER0 = buf;

  // double buffer check
  glGetBooleanv(GL_DOUBLEBUFFER, &state);
  if (!state && buf <= GL_BACK) {
    printf("Warning: GL_DOUBLEBUFFER=0\n");
  }

  // default framebuffer
  glGetIntegerv(GL_FRAMEBUFFER_BINDING, &buf);
  G->ShaderMgr->default_framebuffer_id = buf;
}
#endif

void PyMOL_DrawWithoutLock(CPyMOL * I)
{
  if (!I->done_ConfigureShaders) {
    I->done_ConfigureShaders = true;

    I->G->HaveGUI = I->G->Option->pmgui;

#ifndef PURE_OPENGL_ES_2
    // stereo test with PyQt5 on Linux is broken (QTBUG-59636), so we
    // test for stereo here
    if (I->G->HaveGUI)
    {
      check_gl_stereo_capable(I->G);
    }
#endif

    PyMOL_LaunchStatus_Feedback(I->G);

    I->G->ShaderMgr->Config();

    // OpenGL debugging (glewInit must be called first)
    if (I->G->Option->gldebug) {
#ifdef GL_DEBUG_OUTPUT
      if (!glDebugMessageCallback) {
        printf("glDebugMessageCallback not available\n");
      } else {
        glDebugMessageCallback(gl_debug_proc, NULL);
        glEnable(GL_DEBUG_OUTPUT);
      }
#else
      printf("GL_DEBUG_OUTPUT not available\n");
#endif
    }
  }

  PyMOLGlobals * G = I->G;
  if(I->ModalDraw) {
    if(G->HaveGUI) {
      PyMOL_PushValidContext(I);
      setup_gl_state();
    }
    {
      PyMOLModalDrawFn *fn = I->ModalDraw;
      I->ModalDraw = NULL;      /* always resets to NULL! */
      fn(G);
    }

    if(G->HaveGUI) {
      PyMOL_PopValidContext(I);
    }
  } else {

    if(I->DraggedFlag) {
      if(ControlIdling(I->G)) {
        ExecutiveSculptIterateAll(I->G);
      }
      I->DraggedFlag = false;
    }

    if(G->HaveGUI) {
      PyMOL_PushValidContext(I);

      setup_gl_state();

      if(!I->DrawnFlag) {
        SceneSetCardInfo(G, (char *) glGetString(GL_VENDOR),
                         (char *) glGetString(GL_RENDERER),
                         (char *) glGetString(GL_VERSION));
        if(G->Option->show_splash && !G->Option->quiet) {

          PRINTFB(G, FB_OpenGL, FB_Results)
            " OpenGL graphics engine:\n"
            "  GL_VENDOR:   %s\n"
            "  GL_RENDERER: %s\n"
            "  GL_VERSION:  %s\n",
            (char *) glGetString(GL_VENDOR),
            (char *) glGetString(GL_RENDERER),
            (char *) glGetString(GL_VERSION) ENDFB(G);
          if(Feedback(G, FB_OpenGL, FB_Blather)) {
            printf("  GL_EXTENSIONS: %s\n", (char *) glGetString(GL_EXTENSIONS));
          }
        }
        I->DrawnFlag = true;
      }
    } else {
      I->DrawnFlag = true;
    }

    I->RedisplayFlag = false;

    OrthoBusyPrime(G);
    ExecutiveDrawNow(G);

    if(I->ImageRequestedFlag) {
      if(SceneHasImage(G)) {
        I->ImageReadyFlag = true;
        I->ImageRequestedFlag = false;
      } else {
        I->ImageReadyFlag = false;
      }
    } else if(I->ImageReadyFlag) {
      if(!SceneHasImage(G))
        I->ImageReadyFlag = false;
    }

    if(G->HaveGUI)
      PyMOL_PopValidContext(I);
  }
}

void PyMOL_Key(CPyMOL * I, unsigned char k, int x, int y, int modifiers)
{
  PYMOL_API_LOCK PyMOLGlobals * G = I->G;
  if(!WizardDoKey(G, k, x, y, modifiers))
    OrthoKey(G, k, x, y, modifiers);
  PyMOL_NeedRedisplay(G->PyMOL);
PYMOL_API_UNLOCK}

void PyMOL_Special(CPyMOL * I, int k, int x, int y, int modifiers)
{
  PYMOL_API_LOCK PyMOLGlobals * G = I->G;

  int grabbed = false;
  char buffer[255];
  (void)buffer;
  if(!grabbed)
    grabbed = WizardDoSpecial(G, (unsigned char) k, x, y, modifiers);

  switch (k) {
  case P_GLUT_KEY_UP:
  case P_GLUT_KEY_DOWN:
    grabbed = 1;
    OrthoSpecial(G, k, x, y, modifiers);
    break;
  case P_GLUT_KEY_LEFT:
  case P_GLUT_KEY_RIGHT:
    if(OrthoArrowsGrabbed(G)) {
      grabbed = 1;
      OrthoSpecial(G, k, x, y, modifiers);
    }
    break;
  }

#ifndef _PYMOL_NOPY
  if(!grabbed) {
    sprintf(buffer, "_special %d,%d,%d,%d", k, x, y, modifiers);
    PLog(G, buffer, cPLog_pml);
    PParse(G, buffer);
    PFlush(G);
  }
#endif

PYMOL_API_UNLOCK}

void PyMOL_Reshape(CPyMOL * I, int width, int height, int force)
{
  PYMOL_API_LOCK PyMOLGlobals * G = I->G;

  G->Option->winX = width;
  G->Option->winY = height;
  
  OrthoReshape(G, width, height, force);
PYMOL_API_UNLOCK}

int PyMOL_Idle(CPyMOL * I)
{
  int did_work = false;
  PYMOL_API_TRYLOCK PyMOLGlobals * G = I->G;

  I->DraggedFlag = false;
  if(I->IdleAndReady < IDLE_AND_READY) {
    if(I->DrawnFlag)
      I->IdleAndReady++;
  }
  if(I->FakeDragFlag == 1) {
    I->FakeDragFlag = false;
    OrthoFakeDrag(G);
    did_work = true;
  }

  if(ControlIdling(G)) {
    ExecutiveSculptIterateAll(G);
    ControlSdofIterate(G);
    did_work = true;
  }

  SceneIdle(G);

  if(SceneRovingCheckDirty(G)) {
    SceneRovingUpdate(G);
    did_work = true;
  }
#ifndef _PYMOL_NOPY
  if(PFlush(G)) {
    did_work = true;
  }

  if(I->PythonInitStage > 0) {
    if(I->PythonInitStage < 2) {
      I->PythonInitStage++;
    } else {
      I->PythonInitStage = -1;
      PBlock(G);

      PXDecRef(PYOBJECT_CALLMETHOD
               (G->P_inst->obj, "adapt_to_hardware", "O", G->P_inst->obj));
      
      if(PyErr_Occurred())
        PyErr_Print();

      PXDecRef(PYOBJECT_CALLMETHOD(G->P_inst->obj, "exec_deferred", "O", G->P_inst->obj));

      if(PyErr_Occurred())
        PyErr_Print();

      PUnblock(G);
      PFlush(G);
    }
  }
#endif

  /* reset the interrupt flag if we're not doing anything */

  if(!(did_work || I->ModalDraw))
    if(PyMOL_GetInterrupt(I, false))
      PyMOL_SetInterrupt(I, false);

  PYMOL_API_UNLOCK_NO_FLUSH return (did_work || I->ModalDraw);
}

void PyMOL_ExpireIfIdle(CPyMOL * I)
{
  PYMOL_API_LOCK PyMOLGlobals * G = I->G;
  int final_init_done = true;
#ifndef _PYMOL_NOPY
  final_init_done = (I->PythonInitStage == -1);
#endif

  if(!G->HaveGUI) {
    if(final_init_done) {
      if(!OrthoCommandWaiting(G)) {
        if((!G->Option->keep_thread_alive) && (!G->Option->read_stdin)) {
          I->ExpireCount++;
          if(I->ExpireCount == 10) {
            PParse(G, "_quit");
          }
        }
      }
    }
  }
  PYMOL_API_UNLOCK;
}

void PyMOL_NeedFakeDrag(CPyMOL * I)
{
  I->FakeDragFlag = true;
}

void PyMOL_NeedRedisplay(CPyMOL * I)
{
  I->RedisplayFlag = true;
}

void PyMOL_NeedSwap(CPyMOL * I)
{
  I->SwapFlag = true;
}

void PyMOL_NeedReshape(CPyMOL * I, int mode, int x, int y, int width, int height)
{ 
  PyMOLGlobals *G = I->G;
  if(width < 0) {
    if (!G->HaveGUI)
      return;

    width = SceneGetBlock(G)->getWidth();
    if(SettingGetGlobal_b(G, cSetting_internal_gui))
      width += DIP2PIXEL(SettingGetGlobal_i( G, cSetting_internal_gui_width));
  }

  /* if height is negative, force a reshape based on the current height */

  if(height < 0) {
    int internal_feedback;
    height = SceneGetBlock(G)->getHeight();
    internal_feedback = SettingGetGlobal_i(G, cSetting_internal_feedback);
    if(internal_feedback)
      height += (internal_feedback - 1) * cOrthoLineHeight + cOrthoBottomSceneMargin;
    if(SettingGetGlobal_b(G, cSetting_seq_view)
       && !SettingGetGlobal_b(G, cSetting_seq_view_overlay)) {
      height += SeqGetHeight(G);
    }
    height += MovieGetPanelHeight(G);
  }

  if(G->HaveGUI) {
    float sf = DIP2PIXEL(1);

    I->Reshape[1] = x / sf;
    I->Reshape[2] = y / sf;
    I->Reshape[3] = width / sf;
    I->Reshape[4] = height / sf;

    I->ReshapeFlag = true;
    I->Reshape[0] = mode;
    PyMOL_NeedRedisplay(I);
  } else {
    /* if no gui, then force immediate reshape */
    PyMOLGlobals *G = I->G;

    G->Option->winX = width;
    G->Option->winY = height;

    OrthoReshape(G, width, height, true);
  }
}

int PyMOL_GetIdleAndReady(CPyMOL * I)
{
  return (I->IdleAndReady == IDLE_AND_READY);
}

int PyMOL_GetReshape(CPyMOL * I)
{
  return I->ReshapeFlag;
}

PyMOLreturn_int_array PyMOL_GetReshapeInfo(CPyMOL * I, int reset)
{
  PyMOLreturn_int_array result = { PyMOLstatus_SUCCESS, PYMOL_RESHAPE_SIZE, NULL };
  PYMOL_API_LOCK
  if(reset)
    I->ReshapeFlag = false;
  result.array = VLAlloc(int, PYMOL_RESHAPE_SIZE);
  if(!result.array) {
    result.status = PyMOLstatus_FAILURE;
  } else {
    int a;
    for(a = 0; a < PYMOL_RESHAPE_SIZE; a++)
      result.array[a] = I->Reshape[a];
  }

 PYMOL_API_UNLOCK
   return result;
}

void PyMOL_SetPassive(CPyMOL * I, int onOff)
{
  I->PassiveFlag = onOff;
}

void PyMOL_SetClickReady(CPyMOL * I, const char *name, int index, int button, int mod, int x,
                         int y, const float *pos, int state, int bond)
{
  I->ClickReadyFlag = true;
  I->ClickedIndex = index;
  I->ClickedBondIndex = bond;
  I->ClickedButton = button;
  I->ClickedModifiers = mod;
  I->ClickedX = x;
  I->ClickedY = y;
  I->ClickedPosState = state;

  strcpy(I->ClickedObject, name ? name : "");

  if ((I->ClickedHavePos = bool(pos))) {
    copy3f(pos, I->ClickedPos);
  } else {
    zero3f(I->ClickedPos);
  }
}

int PyMOL_GetClickReady(CPyMOL * I, int reset)
{
  int result = I->ClickReadyFlag;
  if(reset) {
    I->ClickReadyFlag = false;
  }
  return result;
}

char *PyMOL_GetClickString(CPyMOL * I, int reset)
{
  char *result = NULL;
  PYMOL_API_LOCK int ready = I->ClickReadyFlag;
  if(reset)
    I->ClickReadyFlag = false;
  if(ready) {
    size_t result_size = OrthoLineLength + 1;
    result = pymol::malloc<char>(result_size);
    if(result) {
      const char* butstr = "left";
      switch (I->ClickedButton) {
      case P_GLUT_SINGLE_LEFT:
        butstr = "single_left";
        break;
      case P_GLUT_SINGLE_MIDDLE:
        butstr = "single_middle";
        break;
      case P_GLUT_SINGLE_RIGHT:
        butstr = "single_right";
        break;
      case P_GLUT_DOUBLE_LEFT:
        butstr = "double_left";
        break;
      case P_GLUT_DOUBLE_MIDDLE:
        butstr = "double_middle";
        break;
      case P_GLUT_DOUBLE_RIGHT:
        butstr = "double_right";
        break;
      }

      WordType modstr = "";
      if(cOrthoCTRL & I->ClickedModifiers) {
        strcat(modstr, " ctrl");
      }
      if(cOrthoALT & I->ClickedModifiers) {
        strcat(modstr, " alt");
      }
      if(cOrthoSHIFT & I->ClickedModifiers) {
        strcat(modstr, " shift");
      }

      result[0] = 0;
      if(!I->ClickedObject[0]) {
        strcat(result, "type=none\n");
      } else if (auto cobj =
                     ExecutiveFindObjectByName(I->G, I->ClickedObject)) {
        switch (cobj->type) {
        case cObjectMolecule:
          strcat(result, "type=object:molecule\n");
          break;
        case cObjectCGO:
          strcat(result, "type=object:cgo\n");
          break;
        default:
          strcat(result, "type=object\n");
          break;
        }

        snprintf(result + strlen(result), result_size - strlen(result),
            "object=%s\n"
            "index=%d\n"
            "bond=%d\n",
            I->ClickedObject, //
            I->ClickedIndex + 1 /* 1-based */,
            I->ClickedBondIndex /* 0-based or cPickable_t */);

        auto obj = dynamic_cast<ObjectMolecule const*>(cobj);
        if(obj && (I->ClickedIndex < obj->NAtom)) {
          const AtomInfoType* ai = obj->AtomInfo + I->ClickedIndex;
          char inscode_str[2] = { ai->inscode, '\0' };
          snprintf(result + strlen(result), result_size - strlen(result),
              "rank=%d\n"
              "id=%d\n"
              "segi=%s\n"
              "chain=%s\n"
              "resn=%s\n"
              "resi=%d%s\n"
              "name=%s\n"
              "alt=%s\n",
                  ai->rank,
                  ai->id,
                  LexStr(I->G, ai->segi),
                  LexStr(I->G, ai->chain),
                  LexStr(I->G, ai->resn),
                  ai->resv, inscode_str,
                  LexStr(I->G, ai->name),
                  ai->alt);
        }
      }

      snprintf(result + strlen(result), result_size - strlen(result),
          "click=%s\n"
          "mod_keys=%s\n"
          "x=%d\n"
          "y=%d\n",
          butstr, modstr + (modstr[0] == ' ' ? 1 : 0), I->ClickedX,
          I->ClickedY);

      if (I->ClickedHavePos) {
        snprintf(result + strlen(result), result_size - strlen(result),
            "px=%.7g\n"
            "py=%.7g\n"
            "pz=%.7g\n"
            "state=%d\n",
            I->ClickedPos[0], I->ClickedPos[1], I->ClickedPos[2],
            I->ClickedPosState);
      }

      // strip trailing newline
      assert(pymol::zstring_view(result).ends_with('\n'));
      result[strlen(result) - 1] = '\0';
    }
  }
  PYMOL_API_UNLOCK return (result);
}

int PyMOL_GetImageReady(CPyMOL * I, int reset)
{
  int result = I->ImageReadyFlag;
  if(reset) {
    I->ImageReadyFlag = false;
  }
  return result;
}

PyMOLreturn_int_array PyMOL_GetImageInfo(CPyMOL * I)
{
  PyMOLreturn_int_array result = { PyMOLstatus_SUCCESS, 2, NULL };
  PYMOL_API_LOCK result.array = VLAlloc(int, 2);
  if(!result.array) {
    result.status = PyMOLstatus_FAILURE;
  } else {
    std::tie(result.array[0], result.array[1]) = SceneGetImageSize(I->G);
  }
  PYMOL_API_UNLOCK return result;
}

int PyMOL_GetImageData(CPyMOL * I,
                       int width, int height,
                       int row_bytes, void *buffer, int mode, int reset)
{
  int ok = true;
  PYMOL_API_LOCK if(reset)
    I->ImageReadyFlag = false;
  ok = SceneCopyExternal(I->G, width, height, row_bytes, (unsigned char *) buffer, mode);
  PYMOL_API_UNLOCK return get_status_ok(ok);
}

PyMOLreturn_int_array PyMOL_GetImageDataReturned(CPyMOL * I,
                       int width, int height,
                       int row_bytes, int mode, int reset)
{
  PyMOLreturn_int_array result = { PyMOLstatus_SUCCESS, 0, NULL };
  int ok = true;
  int size;
  void *buffer;
  PYMOL_API_LOCK
    if(reset){
      I->ImageReadyFlag = false;
    }
  size = width*height;
  buffer = VLAlloc(int, size);
  ((int*)buffer)[0] = ('A'<<24)|('B'<<16)|('G'<<8)|'R';
  ok = SceneCopyExternal(I->G, width, height, row_bytes, (unsigned char *) buffer, mode);
  if(ok) {
    result.array = (int*) buffer;
    result.size = size;
  } else {
    result.status = PyMOLstatus_FAILURE;
  }
  
  PYMOL_API_UNLOCK return result;
}

int PyMOL_FreeResultString(CPyMOL * I, char *st)
{
  PYMOL_API_LOCK FreeP(st);
  PYMOL_API_UNLOCK return get_status_ok((st != NULL));
}

int PyMOL_GetRedisplay(CPyMOL * I, int reset)
{
  int result = false;

  PYMOL_API_TRYLOCK PyMOLGlobals * G = I->G;
  result = I->RedisplayFlag;

  if(result) {
    if(SettingGet_b(G, NULL, NULL, cSetting_defer_updates)) {
      result = false;
    } else {
      if(reset)
        I->RedisplayFlag = false;
    }
  }
  PYMOL_API_UNLOCK_NO_FLUSH return (result || I->ModalDraw);    /* always true when ModalDraw is set */
}

int PyMOL_GetPassive(CPyMOL * I, int reset)
{                               /* lock intentionally omitted */
  int result = I->PassiveFlag;
  if(reset)
    I->PassiveFlag = false;
  return result;
}

int PyMOL_GetModalDraw(CPyMOL * I)
{
  if(I)
    return (I->ModalDraw != NULL);
  return false;
}

void PyMOL_SetModalDraw(CPyMOL * I, PyMOLModalDrawFn * fn)
{
  I->ModalDraw = fn;
}

int PyMOL_GetSwap(CPyMOL * I, int reset)
{                               /* lock intentionally omitted */
  int result = I->SwapFlag;
  if(reset)
    I->SwapFlag = false;
  return result;
}

int PyMOL_GetBusy(CPyMOL * I, int reset)
{                               /* lock intentionally omitted */
  int result = I->BusyFlag;
  if(reset)
    PyMOL_SetBusy(I, false);
  return result;
}

void PyMOL_SetBusy(CPyMOL * I, int value)
{                               /* lock intentionally omitted */
  if(!I->BusyFlag)
    /* if we weren't busy before, then reset the progress indicators */
    PyMOL_ResetProgress(I);

  I->BusyFlag = value;
}

int PyMOL_GetInterrupt(CPyMOL * I, int reset)
{                               /* lock intentionally omitted */
  if(I) {
    int result = I->InterruptFlag;

    if(reset)
      PyMOL_SetInterrupt(I, false);
    return result;
  } else
    return false;
}

void PyMOL_SetInterrupt(CPyMOL * I, int value)
{                               /* lock intentionally omitted */
  if(I) {
    I->InterruptFlag = value;
    if(I->G)
      I->G->Interrupt = value;
  }
}

void PyMOL_Drag(CPyMOL * I, int x, int y, int modifiers)
{
  PYMOL_API_LOCK OrthoDrag(I->G, x, y, modifiers);
  I->DraggedFlag = true;
PYMOL_API_UNLOCK}

/**
 * Mouse button and keyboard press handler
 *
 * button: mouse button or key code
 * state:
 *   -2 = key press with GLUT_KEY_* special code
 *   -1 = key press with ascii code
 *    0 = mouse down
 *    1 = mouse up
 * x, y: mouse pointer position
 * modifiers: SHIFT/CTRL/ALT bitmask
 */
void PyMOL_Button(CPyMOL * I, int button, int state, int x, int y, int modifiers)
{
  PYMOL_API_LOCK
    if (state == -1) {
      PyMOL_Key(I, (unsigned char)button, x, y, modifiers);
    } else if (state == -2) {
      PyMOL_Special(I, button, x, y, modifiers);
    } else {
      OrthoButton(I->G, button, state, x, y, modifiers);
    }
PYMOL_API_UNLOCK}

void PyMOL_SetSwapBuffersFn(CPyMOL * I, PyMOLSwapBuffersFn * fn)
{
  I->SwapFn = fn;
}

void PyMOL_SwapBuffers(CPyMOL * I)
{
  if(I->SwapFn && I->G->ValidContext) {
    I->SwapFn();
    I->SwapFlag = false;
  } else {
    I->SwapFlag = true;
  }
}

void PyMOL_RunTest(CPyMOL * I, int group, int test)
{
  PYMOL_API_LOCK TestPyMOLRun(I->G, group, test);
PYMOL_API_UNLOCK}

void PyMOL_PushValidContext(CPyMOL * I)
{
  if(I && I->G)
    I->G->ValidContext++;
}

void PyMOL_PopValidContext(CPyMOL * I)
{
  if(I && I->G && (I->G->ValidContext > 0))
    I->G->ValidContext--;
}

void PyMOL_SetStereoCapable(CPyMOL * I, int stereoCapable){
  PYMOL_API_LOCK PyMOLGlobals * G = I->G;
  G->StereoCapable = stereoCapable;
  if (SettingGetGlobal_b(I->G, cSetting_stereo_mode)==0){
    /* if users haven't set stereo_mode, then set it to default */
    if (G->StereoCapable){
      SettingSetGlobal_i(I->G, cSetting_stereo_mode, cStereo_quadbuffer);      /* quadbuffer if we can */
    } else {
      SettingSetGlobal_i(I->G, cSetting_stereo_mode, cStereo_crosseye);      /* otherwise crosseye by default */
    }
  } else if (G->StereoCapable && SettingGetGlobal_b(G, cSetting_stereo)){
    SettingSetGlobal_i(I->G, cSetting_stereo_mode, SettingGetGlobal_b(I->G, cSetting_stereo_mode));
  }
  SceneUpdateStereo(I->G);
  PYMOL_API_UNLOCK
}

void PyMOL_InitializeCMol(CPyMOL * I){
  PYMOL_API_LOCK PyMOLGlobals * G = I->G;
  /* Set stereo_mode to 0, so that PyMOL_SetStereoCapable()
     can determine whether user has changed the stereo_mode.  If
     users have not changed it, then set to quadbuffer if we can,
     or crosseye by default (see above) */
  SettingSetGlobal_i(G, cSetting_stereo_mode, 0);
  PYMOL_API_UNLOCK
}

void PyMOL_SetDefaultMouse(CPyMOL * I)
{
  PYMOL_API_LOCK PyMOLGlobals * G = I->G;

  ButModeSet(G, cButModeLeftNone, cButModeRotXYZ);
  ButModeSet(G, cButModeMiddleNone, cButModeTransXY);
  ButModeSet(G, cButModeRightNone, cButModeTransZ);

  ButModeSet(G, cButModeLeftShft, cButModePotentialClick);
  ButModeSet(G, cButModeMiddleShft, cButModePotentialClick);
  ButModeSet(G, cButModeRightShft, cButModeClipNF);

  ButModeSet(G, cButModeLeftCtrl, cButModePotentialClick);
  ButModeSet(G, cButModeMiddleCtrl, cButModePotentialClick);
  ButModeSet(G, cButModeRightCtrl, cButModePotentialClick);

  ButModeSet(G, cButModeLeftCtSh, cButModePotentialClick);
  ButModeSet(G, cButModeMiddleCtSh, cButModePotentialClick);
  ButModeSet(G, cButModeRightCtSh, cButModePotentialClick);

  ButModeSet(G, cButModeWheelNone, cButModeScaleSlab);
  ButModeSet(G, cButModeWheelShft, cButModeMoveSlab);
  ButModeSet(G, cButModeWheelCtrl, cButModeMoveSlabAndZoom);
  ButModeSet(G, cButModeWheelCtSh, cButModeTransZ);

  ButModeSet(G, cButModeMiddleCtSh, cButModeOrigAt); /* SET TWICE?!? */

  ButModeSet(G, cButModeLeftSingle, cButModeSimpleClick);
  ButModeSet(G, cButModeMiddleSingle, cButModeCent);
  ButModeSet(G, cButModeRightSingle, cButModeSimpleClick);

  ButModeSet(G, cButModeLeftDouble, cButModeSimpleClick);
  ButModeSet(G, cButModeRightDouble, cButModeSimpleClick);

  {
    int a;
    for(a = cButModeLeftShftDouble; a <= cButModeRightCtrlAltShftSingle; a++) {
      ButModeSet(G, a, cButModeSimpleClick);
    }
    for(a = cButModeLeftAlt; a <= cButModeRightCtrlAltShft; a++) {
      ButModeSet(G, a, cButModePotentialClick);
    }

  }
  G->Feedback->currentMask(FB_Scene) &= ~(FB_Results); /* suppress click messages */
PYMOL_API_UNLOCK}

PyMOLreturn_status PyMOL_CmdRock(CPyMOL * I, int mode){
  PyMOLreturn_status result = { PyMOLstatus_SUCCESS };
  PYMOL_API_LOCK PyMOLGlobals * G = I->G;
  ControlRock(G, mode);
  PYMOL_API_UNLOCK return result;
}


PyMOLreturn_string_array PyMOL_CmdGetNames(CPyMOL * I, int mode, const char *s0, int enabled_only){
  PyMOLreturn_string_array result = { PyMOLstatus_FAILURE };
  PYMOL_API_LOCK PyMOLGlobals * G = I->G;
  result = return_result(ExecutiveGetNames(G, mode, enabled_only, s0));
  PYMOL_API_UNLOCK 
  return (result);
}

PyMOLreturn_status PyMOL_CmdMapNew(CPyMOL * I, const char *name, int type, float grid_spacing, 
				   const char *selection, int state, int normalize,
				   int zoom, int quiet){
  PyMOLreturn_status result = { PyMOLstatus_FAILURE };
  float minCorner[3], maxCorner[3];
  PYMOL_API_LOCK 

  minCorner[0] = minCorner[1] = minCorner[2] = 0.;
  maxCorner[0] = maxCorner[1] = maxCorner[2] = 1.;
  auto res = ExecutiveMapNew(I->G, name, type, grid_spacing,
      selection, -1., minCorner, maxCorner, state, 0, quiet, 0, normalize, 1.,
      -1., 0.);
  result.status = get_status_ok(static_cast<bool>(res));

  PYMOL_API_UNLOCK return result;
}

#ifdef _PYMOL_LIB
PyMOLreturn_status PyMOL_SetIsEnabledCallback(CPyMOL * I, void *CallbackObject, void (*enabledCallback)(void *, const char *, int )){
  PyMOLreturn_status result = { PyMOLstatus_FAILURE };
  PYMOL_API_LOCK 
    I->G->CallbackObject = CallbackObject;
    I->G->enabledCallback = enabledCallback;
  result.status =  PyMOLstatus_SUCCESS;
  PYMOL_API_UNLOCK 
  return result;
}

PyMOLreturn_int_array PyMOL_GetRepsInSceneForObject(CPyMOL * I, const char *name){
  PyMOLreturn_int_array result = { PyMOLstatus_SUCCESS, 0, NULL };
  int *retarr = 0;
  PYMOL_API_LOCK
    retarr = ExecutiveGetRepsInSceneForObject(I->G, name);
  if(!retarr) {
    result.status = PyMOLstatus_FAILURE;
  } else {
    result.size = VLAGetSize(retarr);
    result.array = retarr;
  }
  PYMOL_API_UNLOCK return result;
}

PyMOLreturn_int_array PyMOL_GetRepsForObject(CPyMOL * I, const char *name){
  PyMOLreturn_int_array result = { PyMOLstatus_SUCCESS, 0, NULL };
  int *retarr = 0;
  PYMOL_API_LOCK
    retarr = ExecutiveGetRepsForObject(I->G, name);
  if(!retarr) {
    result.status = PyMOLstatus_FAILURE;
  } else {
    result.size = VLAGetSize(retarr);
    result.array = retarr;
  }
  PYMOL_API_UNLOCK return result;
}

static OVreturn_word get_button_code(CPyMOL * I, char *code);
static OVreturn_word get_button_mod_code(CPyMOL * I, char *modcode);
static OVreturn_word get_button_action_code(CPyMOL * I, char *actioncode);
static OVreturn_word get_mouse_mode(CPyMOL * I, char *mousemode);

PyMOLreturn_status PyMOL_SetButton(CPyMOL * I, const char *buttonarg, const char *modifierarg, const char *actionarg){
  int ok = true;
  OVreturn_word button_num, but_mod_num, act_code;
  int but_code;
  PyMOLreturn_status result = { PyMOLstatus_FAILURE };
  OrthoLineType button, modifier, action;
  PYMOL_API_LOCK
    UtilNCopyToLower(button, buttonarg, strlen(buttonarg)+1);
    UtilNCopyToLower(modifier, modifierarg, strlen(modifierarg)+1);
    UtilNCopyToLower(action, actionarg, strlen(actionarg)+1);
    ok = OVreturn_IS_OK(button_num = get_button_code(I, button));
  if (ok) ok = OVreturn_IS_OK((but_mod_num = get_button_mod_code(I, modifier)));
  if (ok) ok = OVreturn_IS_OK((act_code = get_button_action_code(I, action)));

  if (ok){
    /* This is directly from the button() function in controlling.py */
    if (button_num.word < 3){ // normal button (L,M,R)
      if (but_mod_num.word < 4){ 
	// none, shft, ctrl, ctsh
	but_code = button_num.word + 3*but_mod_num.word;
      } else {
	// alt, alsh, alct, alcs
	but_code = button_num.word + 68 + 3*(but_mod_num.word-4);
      }
    } else if (button_num.word < 4){ // wheel
      if (but_mod_num.word < 4){
	// none, shft, ctrl, ctsh
	but_code = 12 + but_mod_num.word;
      } else {
	but_code = 64 + but_mod_num.word - 4;
      }
    } else {
      // single and double clicks
      but_code = (16 + button_num.word -4) + but_mod_num.word * 6;
    }
    ButModeSet(I->G, but_code, act_code.word);
    result.status =  PyMOLstatus_SUCCESS;
  }
  PYMOL_API_UNLOCK return result;
}
#include "buttonmodes.h"

PyMOLreturn_status PyMOL_SetMouseButtonMode(CPyMOL * I, const char *modename){
  int ok = true, i, start;
  OVreturn_word mode;
  PyMOLreturn_status result = { PyMOLstatus_FAILURE };

  PYMOL_API_LOCK
    {
        char *nmodename = (char*)malloc(strlen(modename)+1);
	UtilNCopyToLower((char*)nmodename, modename, strlen(modename)+1);
        ok = OVreturn_IS_OK(mode = get_mouse_mode(I, (char*)nmodename));
        free(nmodename);
    }
  if (ok){
    result.status =  PyMOLstatus_SUCCESS;
    {
      int a;
      /* for Button modes, first initialize all buttons, so that 
	 previous functionality does not linger */
      for(a = 0; a < cButModeInputCount; a++) {
	ButModeSet(I->G, a, initial_button_modes[a]);
      }
    }
    start = button_mode_start[mode.word];
    for (i=0; i<n_button_mode[mode.word]; i++, start+=2){
      ButModeSet(I->G, all_buttons[start], all_buttons[start+1]);
    }
  }
  PYMOL_API_UNLOCK return result;  
}

static OVreturn_word get_button_code(CPyMOL * I, char *code)
{
  OVreturn_word result;
  if(!OVreturn_IS_OK((result = OVLexicon_BorrowFromCString(I->Lex, code))))
    return result;
  return OVOneToOne_GetForward(I->MouseButtonCodeLexicon, result.word);
}
static OVreturn_word get_button_mod_code(CPyMOL * I, char *modcode)
{
  OVreturn_word result;
  if(!OVreturn_IS_OK((result = OVLexicon_BorrowFromCString(I->Lex, modcode))))
    return result;
  return OVOneToOne_GetForward(I->MouseButtonModCodeLexicon, result.word);
}
static OVreturn_word get_button_action_code(CPyMOL * I, char *actioncode)
{
  OVreturn_word result;
  if(!OVreturn_IS_OK((result = OVLexicon_BorrowFromCString(I->Lex, actioncode))))
    return result;
  return OVOneToOne_GetForward(I->MouseButtonActionCodeLexicon, result.word);
}
static OVreturn_word get_mouse_mode(CPyMOL * I, char *mousemode)
{
  OVreturn_word result;
  if(!OVreturn_IS_OK((result = OVLexicon_BorrowFromCString(I->Lex, mousemode))))
    return result;
  return OVOneToOne_GetForward(I->MouseModeLexicon, result.word);
}

PyMOLreturn_status PyMOL_ZoomScene(CPyMOL * I, float scale){
  PyMOLreturn_status result = { PyMOLstatus_FAILURE };
  PYMOL_API_LOCK  
    SceneZoom(I->G, scale);
    PyMOL_NeedRedisplay(I);
  result.status =  PyMOLstatus_SUCCESS;
  PYMOL_API_UNLOCK return result;
}

PyMOLreturn_status PyMOL_TranslateScene(CPyMOL * I, float x, float y, float z){
  PyMOLreturn_status result = { PyMOLstatus_FAILURE };
  PYMOL_API_LOCK  
    SceneTranslate(I->G, x, y, z);
  result.status =  PyMOLstatus_SUCCESS;
  PYMOL_API_UNLOCK return result;
}

#include "palettes.h"

static OVreturn_word get_palette(CPyMOL * I, char *palette)
{
  OVreturn_word result;
  if(!OVreturn_IS_OK((result = OVLexicon_BorrowFromCString(I->Lex, palette))))
    return result;
  return OVOneToOne_GetForward(I->PaletteLexicon, result.word);
}

PyMOLreturn_float_array PyMOL_Spectrum(CPyMOL * I, const char *expression, const char *pal, const char *selection, float minimum, float maximum, int byres, int quiet){
  PyMOLreturn_float_array result = { PyMOLstatus_FAILURE };
  PYMOL_API_LOCK
  int ok = true;
  int digits, first, last, array_pl, ret;
  float min_ret, max_ret;
  char prefix[2];
  OVreturn_word pal_word;
  char *palette = (char*)malloc(strlen(pal)+1);
    UtilNCopyToLower((char*)palette, pal, strlen(pal)+1);
  
  if (ok)
    ok = OVreturn_IS_OK(pal_word = get_palette(I, (char*)palette));  
  free(palette);
  prefix[0] = palette_prefix[pal_word.word];
  prefix[1] = 0;
  array_pl = pal_word.word * 3;
  digits = palette_data[array_pl++];
  first = palette_data[array_pl++];
  last = palette_data[array_pl++];

  ret = ExecutiveSpectrum(I->G, selection, expression, minimum, maximum,
      first, last, prefix, digits, byres, quiet, &min_ret, &max_ret);

  if (ret){
    result.size = 2;
    result.array = VLAlloc(float, 2);
    result.array[0] = min_ret;
    result.array[1] = max_ret;
    result.status = PyMOLstatus_SUCCESS;
  } else {
    result.status = PyMOLstatus_FAILURE;
  }
  PYMOL_API_UNLOCK return result;
}

#endif

PyMOLreturn_value PyMOL_GetVersion(CPyMOL * I){
  int ok = true;
  PyMOLreturn_value result;
  result.status = PyMOLstatus_FAILURE;

  PYMOL_API_LOCK
  if(ok) {
    result.type = PYMOL_RETURN_VALUE_IS_STRING;
    result.string = mstrdup(_PyMOL_VERSION);
    result.status = PyMOLstatus_SUCCESS;
  };
  PYMOL_API_UNLOCK return result;
}

AtomPropertyInfo *PyMOL_GetAtomPropertyInfo(CPyMOL * I, const char *atompropname)
{
  OVreturn_word result;
  if(!OVreturn_IS_OK((result = OVLexicon_BorrowFromCString(I->Lex, atompropname))))
    return NULL;
  result = OVOneToOne_GetForward(I->AtomPropertyLexicon, result.word);
  if(!OVreturn_IS_OK(result))
    return NULL;
  return &I->AtomPropertyInfos[result.word];
}

#ifdef __cplusplus
}
#endif
