
#ifndef _H_Feedback
#define _H_Feedback

/* 

IMPORTANT DEVELOPER NOTICE:

All non-debugging output should pass through the PRINTF and ENDF
macros currently defined below, or through the FeedbackAdd or
FeedbackAutoAdd routines.

Feedback bits are:

results -- DEFAULT: ON

output from a definite action which gives a result, such as an RMS
fit, a measured surface area, etc.

errors -- DEFAULT: ON

complaints which will cause failure at some level.

actions -- DEFAULT: ON

Output regarding actions in progress or completed, but which
don't return a particular result.  Example: loading an object or
creating a selection.

warnings -- DEFAULT: ON

Questionable situations which will not necessarily result in task
failure.  

details -- DEFAULT: ON

Verbose output reflecting details about what is going on.

blather -- DEFAULT: OFF

Output which doesn't fit into the above catogories, and is not likely
to be required except in extreme cases, but doesn't fall into the
category of debugging.

debugging -- DEFAULT: OFF

Text output while would only be of interest to a developer. 

NOTE: Debugging output is the only kind of output which should be sent
directly to standard output (actually, standard error).

NOTE: Debugging output should always be preceeded by the enclosing
function name.

*/



/* WARNING: The following constants are replicated in Python for the purpose
 * of minimize program startup time */

/* Discrete Systems and/or Code Modules */

#define FB_all               0 /* only used for setting */

#define FB_feedback_                 1
#define FB_smiles_parsing            2
#define FB_smiles_creation           3


#define FB_total                     20 /* highest index + 1 */

/* Feedback level bit masks */

#define FB_none            0x00


#define FB_results         0x01
#define FB_errors          0x02
#define FB_actions         0x04
#define FB_warnings        0x08
#define FB_details         0x10
#define FB_blather         0x20
#define FB_debugging       0x80

#define FB_everything      0xFF 

extern char *feedback_Mask;

void feedback_Init(void);
void feedback_Free(void);
void feedback_Push(void);
void feedback_Pop(void);

void feedback_SetMask(unsigned int sysmod,unsigned char mask);
void feedback_Disable(unsigned int sysmod,unsigned char mask);
void feedback_Enable(unsigned int sysmod,unsigned char mask);

/* Mechanism: a high-speed bit test, with no range checking 
 * in order to avoid penalizing performance-senstive code
 * modules which may contain live debugging code.  
 */

#define feedback_(sysmod,mask) (feedback_Mask[sysmod]&mask) 

#define FEEDBACK_MAX_OUTPUT 1024
typedef char feedback_LineType[FEEDBACK_MAX_OUTPUT];

/* Print feedback_ Macros -- this the most flexible and cross-OS
 * portable solution I've come up with for sending output with
 * variable arguments.
*/

#define PRINTFB(sysmod,mask) { if(feedback_(sysmod,mask)) { printf(
#define ENDFB );}}

#define PRINTF { printf(
#define ENDF   );}

/* debugging: goes to stderr */

#define PRINTFD(sysmod) {if(feedback_(sysmod,FB_debugging)) fprintf(stderr,
#define ENDFD   );}

#endif
