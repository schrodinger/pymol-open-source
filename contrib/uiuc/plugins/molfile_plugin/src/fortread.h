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
 *      $RCSfile: fortread.h,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.5 $       $Date: 2006/01/05 00:05:52 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *   Formatted fortran file reading routines used in various plugins.
 *   Each function reads the next record from the file (provided it contains
 *   no more than n elements), optionally swapping its contents before
 *   writing it into dest. 
 *   Returns the number of elements on success, 0 on failure.
 *
 *   TODO: These should perhaps rewind the file to the beginning of the
 *   record on failure.
 *
 ***************************************************************************/

#ifndef FORTREAD_H
#define FORTREAD_H

#include <stdlib.h>
#include <stdio.h>

#include "endianswap.h"

/* Only works with aligned four-byte quantities, swap4_aligned() will */
/* cause a bus error on sume platforms if used on unaligned data      */
static int fortread_4(void *dest, int n, int swap, FILE *fd) {
  int dataBegin, dataEnd, count;

  if (fread(&dataBegin, sizeof(int), 1, fd) != 1) return 0;
  if (swap) swap4_aligned(&dataBegin, 1);
  if ((dataBegin <= 0) || (n < dataBegin/4)) return 0;

  count = fread(dest, 4, dataBegin/4, fd);
  if (count != dataBegin/4) return 0;
  if (swap) swap4_aligned(dest, dataBegin/4);

  if (fread(&dataEnd, sizeof(int), 1, fd) != 1) return 0;
  if (swap) swap4_aligned(&dataBegin, 1);
  if (dataEnd != dataBegin) return 0;

  return count;
}

#endif

