/***************************************************************************
 *cr
 *cr            (C) Copyright 2009 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: vmdconio.h,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.4 $       $Date: 2015/10/11 22:36:27 $
 *
 ***************************************************************************/

/** @file
 * APIs for console output management.  The calling application may 
 * optionally provide callback routines for console output that direct 
 * output to GUI consoles and other places besides stdout.
 */

#ifndef VMDCON_PLUGIN_H
#define VMDCON_PLUGIN_H

/* this has to correspond to vmdconsole.h */
#define VMDCON_ALL       0      /**< "print all messages" log level   */
#define VMDCON_INFO      1      /**< informational messages log level */
#define VMDCON_WARN      2      /**< warning messages" log level      */
#define VMDCON_ERROR     3      /**< error messages log level         */
#define VMDCON_ALWAYS    4      /**< print always log level           */
#define VMDCON_LOG       5      /**< store only in syslog log level   */

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* set default */
#if !defined(THISPLUGIN)
#define THISPLUGIN plugin
#endif

/* forward declaration */
static molfile_plugin_t THISPLUGIN;

/* 
 * Emulate printf. unfortunately, we cannot rely on 
 * snprintf being available, so we have to write to
 * a very large buffer and then free it. :-( 
 */
static int vmdcon_printf(int lvl, const char *fmt, ...) {
  va_list ap;
  char *buf;
  int len;

  /* expand formatted output into a single string */
  buf = (char *)malloc(MOLFILE_BIGBUFSIZ);
  va_start(ap, fmt);
  len = vsprintf(buf, fmt, ap);

  /* 
   * Check result. we may get a segfault, but if not
   * let the user know that he/she is in trouble. 
   */
  if (len >= MOLFILE_BIGBUFSIZ) {
    fprintf(stderr,"WARNING! buffer overflow in vmdcon_printf. %d vs %d.\n",
            len, MOLFILE_BIGBUFSIZ);
    free(buf);
#if 0
    errno=ERANGE; /* this is inherently thread-unsafe */
#endif
    return -1;
  }

  /* 
   * write to registered console output function.
   * fall back to stdout, if vmdcon not available. 
   */
  if (THISPLUGIN.cons_fputs) {
    (*THISPLUGIN.cons_fputs)(lvl, buf);
  } else {
    fputs(buf, stdout);
  }
  free(buf);
  return 0;    
}

#ifdef __cplusplus
}
#endif

#endif /* VMDCON_PLUGIN_H */

