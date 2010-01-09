

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

#ifndef _H_PyMOLOptions
#define _H_PyMOLOptions

#define PYMOL_MAX_OPT_STR 1025

struct _CPyMOLOptions {
  int pmgui, internal_gui, show_splash, internal_feedback, security, game_mode, force_stereo,   /* 1 = force stereo (if possible); -1 = force mono; 0 = autodetect */
   
    winX,
    winY,
    blue_line,
    winPX,
    winPY,
    external_gui,
    siginthand,
    reuse_helper, auto_reinitialize, keep_thread_alive, quiet, incentive_product;

  char after_load_script[PYMOL_MAX_OPT_STR];

  int multisample, window_visible, read_stdin, presentation, defer_builds_mode, full_screen, sphere_mode, stereo_capable,       /* for informing PyMOL as to the capabilities of the context */
    stereo_mode, zoom_mode, no_quit;    /* prevent any action from quitting or killing PyMOL */
  
  /* WARNING: for the sake of forward compability, never delete or
     move any fields in the above ...initialization struct in PyMOL.c */

  /* WARNING: don't add, delete, or change item order unless you also update
     PyMOL.c "CPyMOLOptions Defaults"where this global structure is initialized */
};

#ifndef CPyMOLOptions_DEFINED
typedef struct _CPyMOLOptions CPyMOLOptions;
#define CPyMOLOptions_DEFINED
#endif

#endif
