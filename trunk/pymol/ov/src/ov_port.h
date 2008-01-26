/* MACHINE PROCESSED SOURCE CODE -- DO NOT EDIT */

#ifndef _H_ov_port
#define _H_ov_port

#define OV_JENARIX

#ifdef OV_JENARIX
#include "ov_defines.h"
#include "ov_os.h"
#else
/* headers */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#ifndef WIN32
#include <malloc.h>
#endif

#ifndef OV_NULL
#define OV_NULL ((void*)0)
#endif

#ifndef OV_FALSE
#define OV_FALSE 0
#endif

#ifndef OV_TRUE
#define OV_TRUE 1
#endif

/* how do we inline functions in header files? */

#define OV_INLINE __inline__ static

/* memory management */

#define ov_os_malloc malloc
#define ov_os_calloc calloc
#define ov_os_realloc realloc
#define ov_os_free free

/* termination */

#define ov_os_abort abort
#define ov_os_exit exit

#endif
#endif


