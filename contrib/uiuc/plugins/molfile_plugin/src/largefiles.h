/***************************************************************************
 *cr
 *cr            (C) Copyright 1995-2006 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: largefiles.h,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.1 $       $Date: 2006/01/20 16:38:17 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *   Platform dependent defines for enabling 64-bit file I/O on 32-bit machines
 *
 ***************************************************************************/
 
#if defined(_AIX)
/* Define to enable large file extensions on AIX */
#define _LARGE_FILE
#define _LARGE_FILES
#else
/* Defines which enable LFS I/O interfaces for large (>2GB) file support
 * on 32-bit machines.  These must be defined before inclusion of any
 * system headers.
 */
#define _LARGEFILE_SOURCE
#define _FILE_OFFSET_BITS 64
#endif

