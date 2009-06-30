/***************************************************************************
 *cr
 *cr            (C) Copyright 1995-2009 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: fortread.h,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.9 $       $Date: 2009/04/29 15:45:30 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *   Unformatted and formatted fortran file reading routines for
 *   use in various plugins to simplify dealing with fortran i/o quirks.
 ***************************************************************************/

#ifndef FORTREAD_H
#define FORTREAD_H

#include <stdlib.h>
#include <stdio.h>

#include "endianswap.h"

#define FORT_RECLEN_32BIT 4
#define FORT_RECLEN_64BIT 8

/*
 * consume a record length marker
 */
static int fort_eat_recmark(FILE *ifp, int recmarklen) {
  int tmp;
  if (recmarklen == FORT_RECLEN_64BIT) {
    if (fread(&tmp, 4, 1, ifp) != 1) 
      return -1;
    if (fread(&tmp, 4, 1, ifp) != 1) 
      return -1;
  } else {
    if (fread(&tmp, 4, 1, ifp) != 1) 
      return -1;
  }

  return 0;
}

/*
 * Determine endianness and length of Fortran 
 * record length markers within unformatted binary files.
 */
static int fort_get_endian_reclen(int reclen, int word0, int word1, 
                                  int *swap, int *recmarklen) {
  /* check for 32-bit record length markers */
  if (reclen == word0) {
    *swap=0;
    *recmarklen=FORT_RECLEN_32BIT;
    return 0;
  } else {
    int tmp = word0;
    swap4_aligned(&tmp, 1);   
    if (tmp == reclen) {
      *swap=0;
      *recmarklen=FORT_RECLEN_32BIT;
      return 0;
    }
  }
  
  /* check for 64-bit record length markers */ 
  if (reclen == (word0+word1)) {
    *swap=0;
    *recmarklen=FORT_RECLEN_64BIT;
  } else {
    int tmp0=word0;
    int tmp1=word1;
    swap4_aligned(&tmp0, 1);   
    swap4_aligned(&tmp1, 1);   
    *swap=1;
    *recmarklen=FORT_RECLEN_64BIT;
  }
    
  return -1; /* completely unrecognized record length marker */
}



/*  Unformatted reads.
 *
 *   Each function reads the next record from the file (provided it contains
 *   no more than n elements), optionally swapping its contents before
 *   writing it into dest. 
 *   Returns the number of elements on success, 0 on failure.
 *
 *   TODO: These should perhaps rewind the file to the beginning of the
 *   record on failure.
 */

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

/* Formatted reads:
 *
 * copy at most 'len' characters from source to target. 
 *
 * leading white space is skipped over but counts towards 'len'.
 * the copy stops at first whitspace or a '\0'.
 * unlike strncpy(3) the result will always be \0 terminated.
 *
 * intended for copying (short) strings from formatted fortran  
 * i/o files that must not contain whitespace (e.g. residue names, 
 * atom name/types etc. in .pdb, .psf and alike.). 
 */
static void strnwscpy(char *target, const char *source, const int len) {
  int i, c;

  for (i=0, c=0; i<len; ++i) {
    if (*source == '\0' || (c > 0 && *source == ' ')) {
      break;
    }

    if (*source == ' ') { 
      source++;
    } else {
      *target++ = *source++;
      c++;
    }
  }
  *target = '\0';
}

#endif

