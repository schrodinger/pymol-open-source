/* MACHINE PROCESSED SOURCE CODE -- DO NOT EDIT */

#ifndef _H_ov_port
#define _H_ov_port

#ifdef OV_JENARIX

#include "ov_public.h"

#else
/* headers */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#ifdef __linux__
#include <malloc.h>
#else
#include <stddef.h>
#endif

#include "ov_defines.h"
#include "ov_status.h"

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


