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

#define FB_Isosurf                   1
#define FB_Map                       2
#define FB_Matrix                    3
#define FB_MyPNG                     4

/* layer 1 */

#define FB_Feedback                  12
#define FB_Scene                     13
#define FB_Threads                   14  /* part of P.c */
#define FB_Symmetry                  15
#define FB_Ray                       16
#define FB_Setting                   17
#define FB_Object                    18
#define FB_Ortho                     19

/* layer 2 */

#define FB_CoordSet                  25
#define FB_DistSet                   26

#define FB_ObjectMolecule            30
#define FB_ObjectMap                 31
#define FB_ObjectMesh                32
#define FB_ObjectDist                33 
#define FB_ObjectCGO                 34
#define FB_ObjectCallback            35

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

#define FB_Results         0x01
#define FB_Errors          0x02
#define FB_Actions         0x04
#define FB_Warnings        0x08
#define FB_Details         0x10
#define FB_Blather         0x20
#define FB_Debugging       0x80

#define FB_Everything      0xFF 

extern char *FeedbackMask;

void FeedbackInit(void);
void FeedbackFree(void);
void FeedbackPush(void);
void FeedbackPop(void);

void FeedbackAutoAdd(unsigned int sysmod,unsigned char mask,char *str);
void FeedbackAdd(char *str);

void FeedbackSetMask(unsigned int sysmod,unsigned char mask);
void FeedbackDisable(unsigned int sysmod,unsigned char mask);
void FeedbackEnable(unsigned int sysmod,unsigned char mask);

/* Mechanism: a high-speed bit test, with no range checking 
 * in order to avoid penalizing performance-senstive code
 * modules which may contain live debugging code.  
 */

#define Feedback(sysmod,mask) (FeedbackMask[sysmod]&mask) 

#define FEEDBACK_MAX_OUTPUT 65535  
typedef char FeedbackLineType[FEEDBACK_MAX_OUTPUT];

/* Print Feedback Macros -- this the most flexible and cross-OS
 * portable solution I've come up with for sending output with
 * variable arguments.
*/

#define PRINTFB(sysmod,mask) { FeedbackLineType _FBstr; if(FeedbackMask[sysmod]&mask) {sprintf( _FBstr,
#define ENDFB );  FeedbackAdd(_FBstr);}}

#define PRINTF { FeedbackLineType _FBstr; sprintf( _FBstr,
#define ENDF   ); FeedbackAdd(_FBstr);}

/* debugging: goes to stderr */

#define PRINTFD(sysmod) {if(Feedback(sysmod,FB_Debugging)) fprintf(stderr,
#define ENDFD   );}

#endif
