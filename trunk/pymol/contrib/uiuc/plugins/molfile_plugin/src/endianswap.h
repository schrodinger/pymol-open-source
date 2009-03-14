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
 *      $RCSfile: endianswap.h,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.5 $       $Date: 2006/01/05 00:05:52 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *   Byte swapping routines used in various plugins
 *   There are two versions of each routine, one that's safe to use in
 *   all cases (but is slow) and one that is only safe to use on memory 
 *   addresses that are aligned to the word size that's being byte-swapped
 *   but are much much much faster.  Use the aligned versions of these
 *   routines whenever possible.  The 'ndata' length count parameters and
 *   internal loops should be safe to use on huge memory arrays on 64-bit
 *   machines.
 *
 ***************************************************************************/

#ifndef ENDIAN_SWAP_H
#define ENDIAN_SWAP_H

#ifdef __GNUC__
#if __GNUC__ > 3
typedef int my_aliased_int __attribute__((may_alias));
#else
typedef int my_aliased_int;
#endif
#else
typedef int my_aliased_int;
#endif

/* works on unaligned 2-byte quantities */
static void swap2_unaligned(void *v, long ndata) {
  long i;
  char * dataptr = (char *) v;
  char tmp;

  for (i = 0; i < ndata-1; i += 2) {
    tmp = dataptr[i];
    dataptr[i] = dataptr[i+1];
    dataptr[i+1] = tmp;
  }
}


/* works on unaligned 4-byte quantities */
static void swap4_unaligned(void *v, long ndata) {
  long i;
  char *dataptr;
  char tmp;

  dataptr = (char *) v; 
  for (i=0; i<ndata; i++) {
    tmp = dataptr[0];
    dataptr[0] = dataptr[3];
    dataptr[3] = tmp;
    tmp = dataptr[1];
    dataptr[1] = dataptr[2];
    dataptr[2] = tmp;
    dataptr += 4;
  }
}


/* works on unaligned 8-byte quantities */
static void swap8_unaligned(void *v, long ndata) {
  char *data = (char *) v;
  long i;
  char byteArray[8];
  char *bytePointer;

  for (i=0; i<ndata; i++) {
    bytePointer = data + (i<<3);
    byteArray[0]  =  *bytePointer;
    byteArray[1]  =  *(bytePointer+1);
    byteArray[2]  =  *(bytePointer+2);
    byteArray[3]  =  *(bytePointer+3);
    byteArray[4]  =  *(bytePointer+4);
    byteArray[5]  =  *(bytePointer+5);
    byteArray[6]  =  *(bytePointer+6);
    byteArray[7]  =  *(bytePointer+7);

    *bytePointer     = byteArray[7];
    *(bytePointer+1) = byteArray[6];
    *(bytePointer+2) = byteArray[5];
    *(bytePointer+3) = byteArray[4];
    *(bytePointer+4) = byteArray[3];
    *(bytePointer+5) = byteArray[2];
    *(bytePointer+6) = byteArray[1];
    *(bytePointer+7) = byteArray[0];
  }
}


/* Only works with aligned 2-byte quantities, will cause a bus error */
/* on some platforms if used on unaligned data.                      */
static void swap2_aligned(void *v, long ndata) {
  short *data = (short *) v;
  long i;
  short *N; 

  for (i=0; i<ndata; i++) {
    N = data + i;
    *N=(((*N>>8)&0xff) | ((*N&0xff)<<8));  
  }
}


/* Only works with aligned 4-byte quantities, will cause a bus error */
/* on some platforms if used on unaligned data.                      */
static void swap4_aligned(void *v, long ndata) {
  my_aliased_int *data = (my_aliased_int *) v;
  long i;
  my_aliased_int *N;
  for (i=0; i<ndata; i++) {
    N = data + i;
    *N=(((*N>>24)&0xff) | ((*N&0xff)<<24) | 
        ((*N>>8)&0xff00) | ((*N&0xff00)<<8));
  }
}


/* Only works with aligned 8-byte quantities, will cause a bus error */
/* on some platforms if used on unaligned data.                      */
static void swap8_aligned(void *v, long ndata) {
  /* Use int* internally to prevent bugs caused by some compilers */
  /* and hardware that would potentially load data into an FP reg */
  /* and hose everything, such as the old "jmemcpy()" bug in NAMD */
  my_aliased_int *data = (my_aliased_int *) v;  
  long i;
  my_aliased_int *N; 
  int t0, t1;

  for (i=0; i<ndata; i++) {
    N = data + (i<<1);
    t0 = N[0];
    t0=(((t0>>24)&0xff) | ((t0&0xff)<<24) | 
        ((t0>>8)&0xff00) | ((t0&0xff00)<<8));

    t1 = N[1];
    t1=(((t1>>24)&0xff) | ((t1&0xff)<<24) | 
        ((t1>>8)&0xff00) | ((t1&0xff00)<<8));

    N[0] = t1; 
    N[1] = t0; 
  }
}

#if 0
/* Other implementations that might be faster in some cases */

/* swaps the endianism of an eight byte word. */
void mdio_swap8(double *i) {
        char c;
        char *n;
        n = (char *) i;
        c = n[0];
        n[0] = n[7];
        n[7] = c;
        c = n[1];
        n[1] = n[6];
        n[6] = c;
        c = n[2];
        n[2] = n[5];
        n[5] = c;
        c = n[3];
        n[3] = n[4];
        n[4] = c;
}

#endif

#endif

