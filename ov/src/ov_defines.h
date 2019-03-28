
/* 
 * COPYRIGHT NOTICE: This file contains original source code from the
 * Jenarix (TM) Library, Copyright (C) 2007-8 by Warren L. Delano of
 * DeLano Scientific LLC, Palo Alto, California, United States.
 * Please see the accompanying LICENSE file for further information.
 * All rights not explicitly granted in that LICENSE file are
 * reserved.  It is unlawful to modify or remove this notice.
 * TRADEMARK NOTICE: Jenarix is a Trademark of DeLano Scientific LLC.
*/
#ifndef _H_ov_defines
#define _H_ov_defines


/* defines */

#define OV_INLINE inline
#define OV_STATIC static

#ifndef OV_FALSE
#define OV_FALSE 0
#endif

#ifndef OV_TRUE
#define OV_TRUE 1
#endif


/* NULL pointer */

#ifdef NULL
#define OV_NULL NULL
#else
#ifdef __cplusplus
#define OV_NULL 0
#else
#define OV_NULL ((void*)0)
#endif
#endif


/* if heap tracker is on AND we're running multithreading then we need
   a global mutex for the tracker itself */

#ifdef OV_HEAP_TRACKER
#ifndef OV_OS_FAKE_THREADS
#define OV_HEAP_TRACKER_MUTEX
#endif
#endif


/* workaround for compilers which disallow [0] size arrays */

#define OV_ZERO_ARRAY_SIZE 1


/*
  OV__FILE__ -> NULL     for maximum efficiency
  OV__FILE__ -> __FILE__ for debugging
*/

#define OV__FILE__ __FILE__


/*
  OV__LINE__ -> 0        for maximum efficiency
  OV__LINE__ -> __LINE__ for debugging 
*/
#define OV__LINE__ __LINE__

#endif
