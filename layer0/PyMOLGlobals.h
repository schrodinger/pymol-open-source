

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
#ifndef _H_PyMOLGlobals
#define _H_PyMOLGlobals

namespace pymol
{
class cif_file;
class cif_data;
}; // namespace pymol

/* retina scale factor for ortho gui */
extern int _gScaleFactor;
inline int DIP2PIXEL(int v) { return v * _gScaleFactor; }
inline float DIP2PIXEL(float v) { return v * _gScaleFactor; }

/* all of the private singleton classes associated with a PyMOL instance */


/* this gets included in virtually every PyMOL source file, so keep it
   short and sweet */

typedef int lexidx_t;
typedef int lexborrow_t;

using SelectorID_t = int;
using SelectorMemberOffset_t = int;

using StateIndex_t = int; ///< 0-based state index (C/C++/JavaScript)
using StateIndexPython_t = int; ///< 1-based state index (Python/Settings)

constexpr StateIndex_t cStateAll = -1;
constexpr StateIndex_t cStateCurrent = -2;

typedef struct _CMemoryCache CMemoryCache;
struct CIsosurf;
typedef struct _CTetsurf CTetsurf;
typedef struct _CSphere CSphere;
class CFeedback;
typedef struct _CUtil CUtil;
struct CColor;
struct CMovie;
struct CControl;
struct CButMode;
class COrtho;
typedef struct _CWord CWord;
typedef struct _CCGORenderer CCGORenderer;
typedef struct _CCharacter CCharacter;
struct CPop;
class CScene;
struct CSeq;
struct CSetting;
struct CSettingUnique;
struct CText;
struct CWizard;
typedef struct _CAtomInfo CAtomInfo;
typedef struct _CSculptCache CSculptCache;
typedef struct _CVFont CVFont;
typedef struct _CEditor CEditor;
struct CExecutive;
struct CSeeker;
struct CSelector;
struct CSelectorManager;
typedef struct _CTexture CTexture;
typedef struct _CType CType;
typedef struct _CMain CMain;
typedef struct _CPlugIOManager CPlugIOManager;
struct COpenVR;
struct ObjectMolecule;

class CShaderMgr;
class CMovieScenes;

#ifndef _PYMOL_NOPY
typedef struct _CP_inst CP_inst;
#endif

#ifndef OVLexicon_DEFINED
typedef struct _OVLexicon OVLexicon;
#define OVLexicon_DEFINED
#endif

#ifndef CPyMOLOptions_DEFINED
typedef struct _CPyMOLOptions CPyMOLOptions;
#define CPyMOLOptions_DEFINED
#endif

#ifndef OVCONTEXT_DEFINED
typedef struct _OVContext OVContext;
#define OVCONTEXT_DEFINED
#endif

#ifndef OVONETOONE_DEFINED
typedef struct _OVOneToOne OVOneToOne;
#define OVONETOONE_DEFINED
#endif

#ifndef CPyMOL_DEFINED
typedef struct _CPyMOL CPyMOL;
#define CPyMOL_DEFINED
#endif

#ifndef CGO_DEFINED
class CGO;
#define CGO_DEFINED
#endif

#define cPyMOLGlobals_LaunchStatus_StereoFailed 0x1
#define cPyMOLGlobals_LaunchStatus_MultisampleFailed 0x2

struct PyMOLGlobals {

  /* singleton objects */

  CMemoryCache *MemoryCache;    /* could probably eliminate this... */
  CIsosurf *Isosurf;
  CTetsurf *Tetsurf;
  CSphere *Sphere;
  CFeedback *Feedback;
  CUtil *Util;
  CColor *Color;
  CMovie *Movie;
  CControl *Control;
  CButMode *ButMode;
  COrtho *Ortho;
  CWord *Word;
  CCGORenderer *CGORenderer;
  CCharacter *Character;
  CPop *Pop;
  CScene *Scene;
  CGO *DebugCGO;                /* for debugging rendering */
  CSeq *Seq;
  CSetting *Setting, *Default;
  CSettingUnique *SettingUnique;
  CText *Text;
  CWizard *Wizard;
  CAtomInfo *AtomInfo;
  CSculptCache *SculptCache;
  CVFont *VFont;
  CEditor *Editor;
  CExecutive *Executive;
  CSeeker *Seeker;
  CSelectorManager* SelectorMgr;
  CSelector *Selector;
  CTexture *Texture;
  CType *Type;
  OVContext *Context;
  CMain *Main;                  /* host/platform-specific "main" code */
  CPyMOLOptions *Option;
  CPyMOL *PyMOL;                /* the instance */
  OVLexicon *Lexicon;           /* lexicon for data (e.g. label) strings */
  CPlugIOManager *PlugIOManager;
  CShaderMgr* ShaderMgr;
  COpenVR* OpenVR;

#ifndef _PYMOL_NOPY
  CP_inst *P_inst;
#endif

  /* global variables */

  int HaveGUI;                  /* do we have an OpenGL graphics window or are we
                                 * command-line only? */

  int ValidContext;             /* are we guaranteed to have a valid OpenGL
                                 * context at this exact moment? */

  int Ready;                    /* is the program fully initialized and ready to receive
                                 * messages? */

  int Interrupt;                /* set when we are attempting to abort time-consuming calculations */

  int Terminating;              /* is the program shutting down? */

  /* note that the following four options are also contained in
   * PyMOLOption global -- they exist here as independent globals only
   * because changes haven't yet been made throughout code */

  int StereoCapable;            /* the current graphics context quad buffered? */

  int LaunchStatus;             /* to enable deferred output regarding launch status */

  int Security;                 /* do we warn before potentially executing any
                                 * Python code and ask for their informed consent? */

  int DragDirtyFlag;            /* do we need an extra callback to handle a mouse drag? */

#ifdef _PYMOL_LIB
  void *CallbackObject;
  void (*enabledCallback)(void *, const char *, int );
#endif

  // user defined scenes
  CMovieScenes * scenes;

  // for glDrawBuffer (e.g. GL_BACK, unless we're using QOpenGLWidget)
  int DRAW_BUFFER0;

  struct { lexidx_t
#include "lex_constants.h"
    _; } lex_const;
};


/* if we're running PyMOL as a global singleton (old way / backward
   compatible) then this global variable will contain a pointer to
   PyMOL global state variables */

#ifndef _PYMOL_NOPY
extern PyMOLGlobals *SingletonPyMOLGlobals;
#endif

#define ORTHOCGOARG , CGO *orthoCGO
#define ORTHOCGOARGB  orthoCGO
#define ORTHOCGOARGVAR , orthoCGO
#define ORTHOFASTCGOARGVAR , orthoFastCGO
#define ORTHOCGOARGNULL , NULL

#endif
