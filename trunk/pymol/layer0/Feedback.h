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

#ifndef _H_Feedback
#define _H_Feedback

#include"PyMOLGlobals.h"

struct _CFeedback {
  char *Mask;
  char *Stack;
  int Depth;
};

/* 

IMPORTANT DEVELOPER NOTICE:

All non-debugging output should pass through the PRINTF and ENDF
macros currently defined below, or through the FeedbackAdd or
FeedbackAutoAdd routines.

Feedback bits are:

Results -- DEFAULT: ON

output from a definite action which gives a result, such as an RMS
fit, a measured surface area, etc.

Errors -- DEFAULT: ON

complaints which will cause failure at some level.

Actions -- DEFAULT: ON

Output regarding actions in progress or completed, but which
don't return a particular result.  Example: loading an object or
creating a selection.

Warnings -- DEFAULT: ON

Questionable situations which will not necessarily result in task
failure.  Examples: creation of atom selection which includes no
atoms.  RMS fitting with <4 atoms. 

Details -- DEFAULT: ON

Verbose output reflecting details about what is going on, such as the
number of primitives in a raytracing scene. DEFAULT: ON

Blather -- DEFAULT: OFF

Output which doesn't fit into the above catogories, and is not likely
to be required except in extreme cases, but doesn't fall into the
category of debugging.

Debugging -- DEFAULT: OFF

Text output while would only be of interest to a developer. 

NOTE: Debugging output is the only kind of output which should be sent
directly to standard output (actually, standard error).

NOTE: Debugging output should always be preceeded b the enclosing
function name.

*/



/* WARNING: The following constants are replicated in Python for the purpose
 * of minimize program startup time */

/* Discrete Systems and/or Code Modules */

#define FB_All               0 /* only used for setting */

/* NOTE, the following don't have to be packed, or in order -- we just
   need to record what the maximum index is.  Rember that that
   feedback architecture is purely a performance hack, so expect some
   inconvenience... 
 */

/* layer 0 */


#define FB_Isomesh                   1
#define FB_Map                       2
#define FB_Matrix                    3
#define FB_MyPNG                     4
#define FB_Triangle                  5
#define FB_Match                     6
#define FB_Raw                       7
#define FB_Isosurface                8
#define FB_OpenGL                    9

/* layer 1 */

#define FB_Color                     10
#define FB_CGO                       11
#define FB_Feedback                  12
#define FB_Scene                     13
#define FB_Threads                   14  /* part of P.c */
#define FB_Symmetry                  15
#define FB_Ray                       16
#define FB_Setting                   17
#define FB_Object                    18
#define FB_Ortho                     19
#define FB_Movie                     20
#define FB_Python                    21 /* part of P.c */
#define FB_Extrude                   22
#define FB_Rep                       23
#define FB_Shaker                    24

/* layer 2 */

#define FB_CoordSet                  25
#define FB_DistSet                   26
#define FB_GadgetSet                 27

#define FB_ObjectMolecule            30
#define FB_ObjectMap                 31
#define FB_ObjectMesh                32
#define FB_ObjectDist                33 
#define FB_ObjectCGO                 34
#define FB_ObjectCallback            35
#define FB_ObjectSurface             36
#define FB_ObjectGadget              37
#define FB_ObjectSlice               38

#define FB_RepAngle                  43
#define FB_RepDihedral               44
#define FB_RepWireBond               45
#define FB_RepCylBond                46

#define FB_RepLabel                  48
#define FB_RepSphere                 49
#define FB_RepSurface                50
#define FB_RepMesh                   51
#define FB_RepDot                    52
#define FB_RepNonbonded              53
#define FB_RepNonbondedSphere        54
#define FB_RepDistDash               55
#define FB_RepDistLabel              56
#define FB_RepRibbon                 57
#define FB_RepCartoon                58
#define FB_Sculpt                    59
#define FB_VFont                     60

/* layer 3 */

#define FB_Executive                 70
#define FB_Selector                  71
#define FB_Editor                    72

/* layer 4 */

#define FB_Export                    75
#define FB_CCmd                      76  /* "cmd" is just the python version */
#define FB_API                       77  /* APIEntry/Exit */

/* layer 5 */

#define FB_Main                      80  

#define FB_Total                     81 /* highest index + 1 */

/* Feedback level bit masks */

#define FB_None            0x00

#define FB_Output          0x01
/* python/text output */
#define FB_Results         0x02
/* limited to actual results of an operation...requested measurements, etc. */
#define FB_Errors          0x04
#define FB_Actions         0x08
/* advisories regarding the completion of a */
#define FB_Warnings        0x10
#define FB_Details         0x20
#define FB_Blather         0x40
#define FB_Debugging       0x80

#define FB_Everything      0xFF 

int FeedbackInit(PyMOLGlobals *G,int quiet);
void FeedbackFree(PyMOLGlobals *G);
void FeedbackPush(PyMOLGlobals *G);
void FeedbackPop(PyMOLGlobals *G);

void FeedbackAutoAdd(PyMOLGlobals *G,unsigned int sysmod,unsigned char mask,char *str);
void FeedbackAdd(PyMOLGlobals *G,char *str);

void FeedbackSetMask(PyMOLGlobals *G,unsigned int sysmod,unsigned char mask);
void FeedbackDisable(PyMOLGlobals *G,unsigned int sysmod,unsigned char mask);
void FeedbackEnable(PyMOLGlobals *G,unsigned int sysmod,unsigned char mask);

/* Mechanism: a high-speed bit test, with no range checking 
 * in order to avoid penalizing performance-senstive code
 * modules which may contain live debugging code.  
 */

#define Feedback(G,sysmod,mask) (G->Feedback->Mask[sysmod]&mask) 

/* FEEDBACK_MAX_OUTPUT should be as small as is reasonable
 * since this much space gets consumed on the stack
 * every time we have a PRINTF macro.  One might consider
 * rewriting these macros to consume heap space instead.
*/

#define FEEDBACK_MAX_OUTPUT 255
typedef char FeedbackLineType[FEEDBACK_MAX_OUTPUT];

/* Print Feedback Macros -- this the most flexible and cross-OS
 * portable solution I've come up with for sending output with
 * variable arguments.
*/

#define PRINTFB(G,sysmod,mask) { FeedbackLineType _FBstr; if(Feedback(G,sysmod,mask)) {sprintf( _FBstr,
#define ENDFB(G) );  FeedbackAdd(G,_FBstr);}}

#define PRINTF { FeedbackLineType _FBstr; sprintf( _FBstr,
#define ENDF(G)  ); FeedbackAdd(G,_FBstr);}

/* debugging: goes to stderr */

#define PRINTFD(G,sysmod) {if(Feedback(G,sysmod,FB_Debugging)) { fprintf(stderr,
#define ENDFD   );fflush(stderr);}}

/* convenient vector dumping routine */

#define ENDFD3f(v) );fprintf(stderr,": %8.3f %8.3f %8.3f\n",v[0],v[1],v[2]);fflush(stderr);}}

#endif
