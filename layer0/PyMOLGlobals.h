/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* copyright 1998-2000 by Warrn Lyford Delano of DeLano Scientific. 
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

/* all of the private singleton classes associated with a PyMOL instance */

/* this gets included in virtually every PyMOL source file, so keep it
   short and sweet */

typedef struct _CMemoryCache CMemoryCache;
typedef struct _CIsosurf CIsosurf;
typedef struct _CTetsurf CTetsurf;
typedef struct _CSphere CSphere;
typedef struct _CFeedback CFeedback;
typedef struct _CUtil CUtil;
typedef struct _CColor CColor;
typedef struct _CMovie CMovie;
typedef struct _CControl CControl;
typedef struct _CButMode CButMode;
typedef struct _COrtho  COrtho;
typedef struct _CWord  CWord;
typedef struct _CCGORenderer CCGORenderer;
typedef struct _CCharacter CCharacter;
typedef struct _CPop CPop;
typedef struct _CScene CScene;
typedef struct _CSeq CSeq;
typedef struct _CSetting CSetting;
typedef struct _CText CText;
typedef struct _CWizard CWizard;
typedef struct _CAtomInfo CAtomInfo;
typedef struct _CSculptCache CSculptCache;
typedef struct _CVFont CVFont;
typedef struct _CEditor CEditor;
typedef struct _CExecutive CExecutive;
typedef struct _CSeeker CSeeker;
typedef struct _CSelector CSelector;
typedef struct _CTexture CTexture;
typedef struct _CMain CMain;
typedef struct _CGO CGO;

#ifndef CPyMOLOptions_DEFINED
typedef struct _CPyMOLOptions CPyMOLOptions;
#define CPyMOLOptions_DEFINED
#endif

#ifndef OVCONTEXT_DEFINED
typedef struct _OVContext OVContext;
#define OVCONTEXT_DEFINED
#endif

#ifndef CPyMOL_DEFINED
typedef struct _CPyMOL CPyMOL;
#define CPyMOL_DEFINED
#endif

typedef struct _PyMOLGlobals PyMOLGlobals;
struct _PyMOLGlobals {

  /* singleton objects */

  CMemoryCache *MemoryCache; /* could probably eliminate this... */
  CIsosurf  *Isosurf;
  CTetsurf  *Tetsurf;
  CSphere   *Sphere;
  CFeedback *Feedback;
  CUtil     *Util;
  CColor    *Color;
  CMovie    *Movie;
  CControl  *Control;
  CButMode  *ButMode;
  COrtho    *Ortho;
  CWord     *Word;
  CCGORenderer *CGORenderer;
  CCharacter   *Character;
  CPop      *Pop;
  CScene    *Scene;
  CGO       *DebugCGO; /* for debugging rendering */
  CSeq      *Seq;
  CSetting  *Setting;
  CText     *Text;
  CWizard   *Wizard;
  CAtomInfo *AtomInfo;
  CSculptCache *SculptCache;
  CVFont    *VFont;
  CEditor   *Editor;
  CExecutive *Executive;
  CSeeker   *Seeker;
  CSelector *Selector;
  CTexture  *Texture;
  OVContext *Context;
  CMain     *Main; /* host/platform-specific "main" code */
  CPyMOLOptions *Option; 
  CPyMOL    *PyMOL; /* the instance */
  /* global variables */

  int HaveGUI; /* do we have an OpenGL graphics window or are we
                * command-line only? */

  int ValidContext; /* are we guaranteed to have a valid OpenGL
                     * context at this exact moment? */

  int Ready; /* is the program fully initialized and ready to receive
                   * messages? */

  int Terminating; /* is the program shutting down? */
  
  /* note that the following four options are also contained in
   * PyMOLOption global -- they exist here as independent globals only
   * because changes haven't yet been made throughout code */

  int StereoCapable; /* the current graphics context quad buffered? */
  
  int Security; /* do we warn the use before potentially executing any
                 * Python code and ask for their informed consent? */

  int DragDirtyFlag; /* do we need an extra callback to handle a mouse drag? */


};

/* the following transitional global will disappear once we've
   completely removed global state from PyMOL's code */

#ifndef _PYMOL_NOPY
extern PyMOLGlobals *TempPyMOLGlobals;
#endif

/* not a global, but CRay widely used and Ray.h definitely isn't a
 * lightweight include... */

typedef struct _CRay               CRay;


#endif

