/* 
 * COPYRIGHT NOTICE: This file contains original source code from the
 * Jenarix (TM) Library, Copyright (C) 2007-8 by Warren L. Delano of
 * DeLano Scientific LLC, Palo Alto, California, United States.
 * Please see the accompanying LICENSE file for further information.
 * All rights not explicitly granted in that LICENSE file are
 * reserved.  It is unlawful to modify or remove this notice.
 * TRADEMARK NOTICE: Jenarix is a Trademark of DeLano Scientific LLC.
*/
#ifndef _H_ov_status
#define _H_ov_status

/* 
   Maintenance of status is a major headache in C due to the lack of
   exception handling.  Basically, every function call lacking default
   failsafe behavior must check for errors and handle them
   appropriately.  

   To minimize this headache, the following patterns are being
   adopted:

   - A fallible C function must either return a status code or exhibit
     failsafe behavior.  It is never acceptable leak memory or crash.

   - Negative return status codes indicate an error condition whereas
     zero or positive return status codes imply success and may convey
     additional information.

   - An example of failsafe behavior would be a function that returns
     a pointer on success or NULL on failure, or a function which
     performs an integer computation but simply returns zero when an
     error occurs.  In both cases, if the caller cares about knowing
     when an error occurs, then the caller should use an alternate
     function with explicit error checking.  

   - Internally, implementations of fallible functions typically
     maintain a variable named 'status' that is initialized to
     OV_STATUS_SUCCESS (or OV_STATUS_FAILURE) and are then updated as
     appropriate throughout the function to keep track of progress.

   - OV_OK and OV_ERR macros provide a convient way of testing status
     returned from calls (but watch out for short-circuit behavior!)

     if(OV_OK( status = fn_call(...) )) { ... }

   - OV_PTR macro provides a convenient way of generating status based
     on confirming that one or more pointers are not NULL pointers:

     if(OV_OK( status = OV_PTR(ptr) )) {  ... } 

     if(OV_OK( status = OV_PTR(ptr1 && ptr2 && ptr3) )) { ... }

*/

/* macros for interpreting status results: */

#define OV_OK(s) ((s)>=0)
#define OV_ERR(s) ((s)<0)

/* asserting non-NULL pointer */

#define OV_PTR(p) ( (p) ? OV_SUCCESS : OV_STATUS_NULL_PTR )

/* successful returns are always >=0, and function or method results
   (if any) must be valid */

#define OV_STATUS_SUCCESS 0

/* default success status is also passively negative */

#define OV_STATUS_YES     1
#define OV_STATUS_NO      OV_STATUS_SUCCESS

/* error returns are always <0 */

#define OV_STATUS_FAILURE                     -1
#define OV_STATUS_NULL_PTR                    -2
#define OV_STATUS_OUT_OF_MEMORY               -3
#define OV_STATUS_NOT_FOUND                   -4
#define OV_STATUS_DUPLICATE                   -5
#define OV_STATUS_NO_EFFECT                   -6
#define OV_STATUS_MISMATCH                    -7
#define OV_STATUS_INVALID_INC_REF             -8
#define OV_STATUS_INVALID_DEC_REF             -9
#define OV_STATUS_BAD_PTR_OR_CORRUPT_HEAP    -10
#define OV_STATUS_IMPOSSIBLE_CAST            -11
#define OV_STATUS_THREAD_CREATION_FAILURE    -12
#define OV_STATUS_THREADING_ERROR            -13
#define OV_STATUS_MUTEX_ERROR                -14
#define OV_STATUS_COND_ERROR                 -15
#define OV_STATUS_RLOCK_ERROR                -16
#define OV_STATUS_OS_ERROR                   -17
#define OV_STATUS_INVALID_PIPE               -18
#define OV_STATUS_BROKEN_PIPE                -19
#define OV_STATUS_PIPE_TRANSPORT_ERROR       -20
#define OV_STATUS_INVALID_HANDLE             -21
#define OV_STATUS_ACTIVE_HANDLE_FREED_ERROR  -22
#define OV_STATUS_INVALID_OBJECT             -23
#define OV_STATUS_BUFFER_EXHAUSTED           -24
#define OV_STATUS_INVALID_IDENTIFIER         -25
#define OV_STATUS_TLS_ERROR                  -26
#define OV_STATUS_IMPOSSIBLE_CONVERSION      -27

/* status codes below this number are dynamic exception identifiers
   which consume resources owned by the local environment
   (e.g. ov_node) */

#define OV_STATUS_EXCEPTION_START         -1024

/* convenience aliases */

#define OV_SUCCESS OV_STATUS_SUCCESS
#define OV_YES  OV_STATUS_YES
#define OV_NO   OV_STATUS_NO

#endif
