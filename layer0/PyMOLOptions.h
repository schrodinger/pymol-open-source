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

#ifndef _H_PyMOLOptions
#define _H_PyMOLOptions

#define PYMOL_MAX_OPT_STR 1025

struct _COption {
  int pmgui, 
    internal_gui, 
    show_splash,
    internal_feedback, 
    security, 
    game_mode,
    force_stereo, /* 1 = force stereo (if possible); -1 = force mono; 0 = autodetect */
    winX, 
    winY, 
    blue_line,
    winPX, 
    winPY, 
    external_gui, 
    siginthand,
    reuse_helper, 
    auto_reinitialize, 
    keep_thread_alive, 
    quiet, 
    incentive_product;

  char after_load_script[PYMOL_MAX_OPT_STR];

  int multisample,
    window_visible, 
    read_stdin;
  
  /* WARNING: for the sake of forward compability, never delete or
     more any fields in the above */

  /* WARNING: don't delete items or change order unless you also update
     ClassPyMOL.c where this global structure is initialized */
};

#ifndef COption_DEFINED
typedef struct _COption COption;
#define COption_DEFINED
#endif

#endif
