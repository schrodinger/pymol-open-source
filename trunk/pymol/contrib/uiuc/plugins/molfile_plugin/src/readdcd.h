/***************************************************************************
 *cr
 *cr            (C) Copyright 1995-2005 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: readdcd.h,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.33 $       $Date: 2005/01/20 16:48:19 $
 *
 ***************************************************************************/

#ifndef READ_DCD_H
#define READ_DCD_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "endianswap.h"
#include "fastio.h"

/*  DEFINE ERROR CODES THAT MAY BE RETURNED BY DCD ROUTINES */
#define DCD_SUCCESS      0  /* No problems                     */
#define DCD_EOF         -1  /* Normal EOF                      */
#define DCD_DNE         -2  /* DCD file does not exist         */
#define DCD_OPENFAILED  -3  /* Open of DCD file failed         */
#define DCD_BADREAD     -4  /* read call on DCD file failed    */
#define DCD_BADEOF      -5  /* premature EOF found in DCD file */
#define DCD_BADFORMAT   -6  /* format of DCD file is wrong     */
#define DCD_FILEEXISTS  -7  /* output file already exists      */
#define DCD_BADMALLOC   -8  /* malloc failed                   */

/*
 * Read the header information from a dcd file.
 * Input: fd - a file struct opened for binary reading.
 * Output: 0 on success, negative error code on failure.
 * Side effects: *natoms set to number of atoms per frame
 *               *nsets set to number of frames in dcd file
 *               *istart set to starting timestep of dcd file
 *               *nsavc set to timesteps between dcd saves
 *               *delta set to value of trajectory timestep
 *               *nfixed set to number of fixed atoms 
 *               *freeind may be set to heap-allocated space
 *               *reverse set to one if reverse-endian, zero if not.
 *               *charmm set to internal code for handling charmm data.
 */
static int read_dcdheader(fio_fd fd, int *natoms, int *nsets, int *istart, int *nsavc, 
                   double *delta, int *nfixed, int **freeind, 
                   float **fixedcoords, int *reverse, int *charmm);

/* 
 * Read a dcd timestep from a dcd file
 * Input: fd - a file struct opened for binary reading, from which the 
 *             header information has already been read.
 *        natoms, nfixed, first, *freeind, reverse, charmm - the corresponding 
 *             items as set by read_dcdheader
 *        first - true if this is the first frame we are reading.
 *        x, y, z: space for natoms each of floats.
 *        unitcell - space for six floats to hold the unit cell data.  
 *                   Not set if no unit cell data is present.
 * Output: 0 on success, negative error code on failure.
 * Side effects: x, y, z contain the coordinates for the timestep read.
 *               unitcell holds unit cell data if present.
 */
static int read_dcdstep(fio_fd fd, int natoms, float *x, float *y, float *z, 
                        float *unitcell, int nfixed, int first, int *freeind, 
                        float *fixedcoords, int reverse, int charmm);

/* 
 * Skip past a timestep.  If there are fixed atoms, this cannot be used with
 * the first timestep.  
 * Input: fd - a file struct from which the header has already been read
 *        natoms - number of atoms per timestep
 *        nfixed - number of fixed atoms
 *        charmm - charmm flags as returned by read_dcdheader
 * Output: 0 on success, negative error code on failure.
 * Side effects: One timestep will be skipped; fd will be positioned at the
 *               next timestep.
 */
static int skip_dcdstep(fio_fd fd, int natoms, int nfixed, int charmm);

/*
 * clean up dcd data
 * Input: nfixed, freeind - elements as returned by read_dcdheader
 * Output: None
 * Side effects: Space pointed to by freeind is freed if necessary.
 */
static void close_dcd_read(int *freeind, float *fixedcoords);

/*
 * Write a header for a new dcd file
 * Input: fd - file struct opened for binary writing
 *        remarks - string to be put in the remarks section of the header.  
 *                  The string will be truncated to 70 characters.
 *        natoms, istart, nsavc, delta - see comments in read_dcdheader
 * Output: 0 on success, negative error code on failure.
 * Side effects: Header information is written to the dcd file.
 */
static int write_dcdheader(fio_fd fd, const char *remarks, int natoms, 
                    int istart, int nsavc, double delta, int with_unitcell, 
                    int charmm);

/* 
 * Write a timestep to a dcd file
 * Input: fd - a file struct for which a dcd header has already been written
 *       curframe: Count of frames written to this file, starting with 1.
 *       curstep: Count of timesteps elapsed = istart + curframe * nsavc.
 *        natoms - number of elements in x, y, z arrays
 *        x, y, z: pointers to atom coordinates
 * Output: 0 on success, negative error code on failure.
 * Side effects: coordinates are written to the dcd file.
 */
static int write_dcdstep(fio_fd fd, int curstep, int curframe, 
                 int natoms, const float *x, const float *y, const float *z,
                 const double *unitcell, int charmm);



#define DCD_IS_CHARMM       0x01
#define DCD_HAS_4DIMS       0x02
#define DCD_HAS_EXTRA_BLOCK 0x04

/* READ Macro to make porting easier */
#define READ(fd, buf, size) \
        fio_fread(((void *) buf), (size), 1, (fd))


/* WRITE Macro to make porting easier */
#define WRITE(fd, buf, size) \
        fio_fwrite(((void *) buf), (size), 1, (fd))

/* XXX This is broken - fread never returns -1 */
#define CHECK_FREAD(X, msg)  if (X==-1) \
			     { \
				return(DCD_BADREAD); \
			     }

#define CHECK_FEOF(X, msg)  if (X==0) \
			     { \
				return(DCD_BADEOF); \
			     }

static int read_dcdheader(fio_fd fd, int *N, int *NSET, int *ISTART, 
                   int *NSAVC, double *DELTA, int *NAMNF, 
                   int **FREEINDEXES, float **fixedcoords, int *reverseEndian, 
                   int *charmm)
{
  int input_integer;  /* buffer space */
  int i, ret_val;
  char hdrbuf[84];    /* char buffer used to store header */
  int NTITLE;

  /*  First thing in the file should be an 84 */
  ret_val = READ(fd, &input_integer, sizeof(int));
  CHECK_FREAD(ret_val, "reading first int from dcd file");
  CHECK_FEOF(ret_val, "reading first int from dcd file");

  /* Check magic number in file header and determine byte order*/
  if (input_integer != 84) {
    /* check to see if its merely reversed endianism     */
    /* rather than a totally incorrect file magic number */
    swap4_aligned(&input_integer, 1);

    if (input_integer == 84) {
      *reverseEndian=1;
    } else {
      /* not simply reversed endianism, but something rather more evil */
      return DCD_BADFORMAT;
    }
  } else {
    *reverseEndian=0;    
  }

  /* Buffer the entire header for random access */
  ret_val = READ(fd, hdrbuf, 84);
  CHECK_FREAD(ret_val, "buffering header");
  CHECK_FEOF(ret_val, "buffering header");

  /* Check for the ID string "COORD" */
  if (hdrbuf[0] != 'C' || hdrbuf[1] != 'O' ||
      hdrbuf[2] != 'R' || hdrbuf[3] != 'D') {
    return DCD_BADFORMAT;
  }

  /* CHARMm-genereate DCD files set the last integer in the     */
  /* header, which is unused by X-PLOR, to its version number.  */
  /* Checking if this is nonzero tells us this is a CHARMm file */
  /* and to look for other CHARMm flags.                        */
  if (*((int *) (hdrbuf + 80)) != 0) {
    (*charmm) = DCD_IS_CHARMM;
    if (*((int *) (hdrbuf + 44)) != 0)
      (*charmm) |= DCD_HAS_EXTRA_BLOCK;

    if (*((int *) (hdrbuf + 48)) == 1)
      (*charmm) |= DCD_HAS_4DIMS;
  } else {
    (*charmm) = 0;
  }

  /* Store the number of sets of coordinates (NSET) */
  (*NSET) = *((int *) (hdrbuf + 4));
  if (*reverseEndian) swap4_unaligned(NSET, 1);

  /* Store ISTART, the starting timestep */
  (*ISTART) = *((int *) (hdrbuf + 8));
  if (*reverseEndian) swap4_unaligned(ISTART, 1);

  /* Store NSAVC, the number of timesteps between dcd saves */
  (*NSAVC) = *((int *) (hdrbuf + 12));
  if (*reverseEndian) swap4_unaligned(NSAVC, 1);

  /* Store NAMNF, the number of fixed atoms */
  (*NAMNF) = *((int *) (hdrbuf + 36));
  if (*reverseEndian) swap4_unaligned(NAMNF, 1);

  /* Read in the timestep, DELTA */
  /* Note: DELTA is stored as a double with X-PLOR but as a float with CHARMm */
  if ((*charmm) & DCD_IS_CHARMM) {
    float ftmp;
    ftmp = *((float *)(hdrbuf+40)); /* is this safe on Alpha? */
    if (*reverseEndian)
      swap4_aligned(&ftmp, 1);

    *DELTA = (double)ftmp;
  } else {
    (*DELTA) = *((double *)(hdrbuf + 40));
    if (*reverseEndian) swap8_unaligned(DELTA, 1);
  }

  /* Get the end size of the first block */
  ret_val = READ(fd, &input_integer, sizeof(int));
  CHECK_FREAD(ret_val, "reading second 84 from dcd file");
  CHECK_FEOF(ret_val, "reading second 84 from dcd file");
  if (*reverseEndian) swap4_aligned(&input_integer, 1);

  if (input_integer != 84) {
    return DCD_BADFORMAT;
  }

  /* Read in the size of the next block */
  ret_val = READ(fd, &input_integer, sizeof(int));
  CHECK_FREAD(ret_val, "reading size of title block");
  CHECK_FEOF(ret_val, "reading size of title block");
  if (*reverseEndian) swap4_aligned(&input_integer, 1);

  if (((input_integer-4) % 80) == 0) {
    /* Read NTITLE, the number of 80 character title strings there are */
    ret_val = READ(fd, &NTITLE, sizeof(int));
    CHECK_FREAD(ret_val, "reading NTITLE");
    CHECK_FEOF(ret_val, "reading NTITLE");
    if (*reverseEndian) swap4_aligned(&NTITLE, 1);

    for (i=0; i<NTITLE; i++) {
      fio_fseek(fd, 80, FIO_SEEK_CUR);
      CHECK_FEOF(ret_val, "reading TITLE");
    }

    /* Get the ending size for this block */
    ret_val = READ(fd, &input_integer, sizeof(int));

    CHECK_FREAD(ret_val, "reading size of title block");
    CHECK_FEOF(ret_val, "reading size of title block");
  } else {
    return DCD_BADFORMAT;
  }

  /* Read in an integer '4' */
  ret_val = READ(fd, &input_integer, sizeof(int));
  CHECK_FREAD(ret_val, "reading a '4'");
  CHECK_FEOF(ret_val, "reading a '4'");
  if (*reverseEndian) swap4_aligned(&input_integer, 1);

  if (input_integer != 4) {
    return DCD_BADFORMAT;
  }

  /* Read in the number of atoms */
  ret_val = READ(fd, N, sizeof(int));
  CHECK_FREAD(ret_val, "reading number of atoms");
  CHECK_FEOF(ret_val, "reading number of atoms");
  if (*reverseEndian) swap4_aligned(N, 1);

  /* Read in an integer '4' */
  ret_val = READ(fd, &input_integer, sizeof(int));
  CHECK_FREAD(ret_val, "reading a '4'");
  CHECK_FEOF(ret_val, "reading a '4'");
  if (*reverseEndian) swap4_aligned(&input_integer, 1);

  if (input_integer != 4) {
    return DCD_BADFORMAT;
  }

  *FREEINDEXES = NULL;
  *fixedcoords = NULL;
  if (*NAMNF != 0) {
    (*FREEINDEXES) = (int *) calloc(((*N)-(*NAMNF)), sizeof(int));
    if (*FREEINDEXES == NULL)
      return DCD_BADMALLOC;

    *fixedcoords = (float *) calloc((*N)*4 - (*NAMNF), sizeof(float));
    if (*fixedcoords == NULL)
      return DCD_BADMALLOC;

    /* Read in index array size */
    ret_val = READ(fd, &input_integer, sizeof(int));
    CHECK_FREAD(ret_val, "reading size of index array");
    CHECK_FEOF(ret_val, "reading size of index array");
    if (*reverseEndian) swap4_aligned(&input_integer, 1);

    if (input_integer != ((*N)-(*NAMNF))*4) {
      return DCD_BADFORMAT;
    }

    ret_val = READ(fd, (*FREEINDEXES), ((*N)-(*NAMNF))*sizeof(int));
    CHECK_FREAD(ret_val, "reading size of index array");
    CHECK_FEOF(ret_val, "reading size of index array");

    if (*reverseEndian)
      swap4_aligned((*FREEINDEXES), ((*N)-(*NAMNF)));

    ret_val = READ(fd, &input_integer, sizeof(int));
    CHECK_FREAD(ret_val, "reading size of index array");
    CHECK_FEOF(ret_val, "reading size of index array");
    if (*reverseEndian) swap4_aligned(&input_integer, 1);

    if (input_integer != ((*N)-(*NAMNF))*4) {
      return DCD_BADFORMAT;
    }
  }

  return DCD_SUCCESS;
}

static int read_charmm_extrablock(fio_fd fd, int charmm, int reverseEndian,
                                  float *unitcell) {
  int i, input_integer;

  if ((charmm & DCD_IS_CHARMM) && (charmm & DCD_HAS_EXTRA_BLOCK)) {
    /* Leading integer must be 48 */
    if (fio_fread(&input_integer, sizeof(int), 1, fd) != 1)
      return DCD_BADREAD; 
    if (reverseEndian) swap4_aligned(&input_integer, 1);
    if (input_integer == 48) {
      double tmp[6];
      if (fio_fread(tmp, 48, 1, fd) != 1) return DCD_BADREAD;
      if (reverseEndian) 
        swap8_aligned(tmp, 6);
      for (i=0; i<6; i++) unitcell[i] = (float)tmp[i];
    } else {
      /* unrecognized block, just skip it */
      if (fio_fseek(fd, input_integer, FIO_SEEK_CUR)) return DCD_BADREAD;
    }
    if (fio_fread(&input_integer, sizeof(int), 1, fd) != 1) return DCD_BADREAD; 
  } 

  return DCD_SUCCESS;
}

static int read_fixed_atoms(fio_fd fd, int N, int num_free, const int *indexes,
                            int reverseEndian, const float *fixedcoords, 
                            float *freeatoms, float *pos) {
  int i, input_integer;
  
  /* Read leading integer */
  if (fio_fread(&input_integer, sizeof(int), 1, fd) != 1) return DCD_BADREAD;
  if (reverseEndian) swap4_aligned(&input_integer, 1);
  if (input_integer != 4*num_free) return DCD_BADFORMAT;
  
  /* Read free atom coordinates */
  if (fio_fread(freeatoms, 4*num_free, 1, fd) != 1) return DCD_BADREAD;
  if (reverseEndian)
    swap4_aligned(freeatoms, num_free);

  /* Copy fixed and free atom coordinates into position buffer */
  memcpy(pos, fixedcoords, 4*N);
  for (i=0; i<num_free; i++)
    pos[indexes[i]-1] = freeatoms[i];

  /* Read trailing integer */ 
  if (fio_fread(&input_integer, sizeof(int), 1, fd) != 1) return DCD_BADREAD;
  if (reverseEndian) swap4_aligned(&input_integer, 1);
  if (input_integer != 4*num_free) return DCD_BADFORMAT;

  return DCD_SUCCESS;
}
  
static int read_charmm_4dim(fio_fd fd, int charmm, int reverseEndian) {
  int input_integer;

  /* If this is a CHARMm file and contains a 4th dimension block, */
  /* we must skip past it to avoid problems                       */
  if ((charmm & DCD_IS_CHARMM) && (charmm & DCD_HAS_4DIMS)) {
    if (fio_fread(&input_integer, sizeof(int), 1, fd) != 1) return DCD_BADREAD;  
    if (reverseEndian) swap4_aligned(&input_integer, 1);
    if (fio_fseek(fd, input_integer, FIO_SEEK_CUR)) return DCD_BADREAD;
    if (fio_fread(&input_integer, sizeof(int), 1, fd) != 1) return DCD_BADREAD;  
  }

  return DCD_SUCCESS;
}

static int read_dcdstep(fio_fd fd, int N, float *X, float *Y, float *Z, 
                        float *unitcell, int num_fixed,
                        int first, int *indexes, float *fixedcoords, 
                        int reverseEndian, int charmm) {
  int ret_val;   /* Return value from read */

  if ((num_fixed==0) || first) {
    int tmpbuf[6];      /* temp storage for reading formatting info */
    fio_iovec iov[7];   /* I/O vector for fio_readv() call          */
    fio_size_t readlen; /* number of bytes actually read            */
    int i;

    /* if there are no fixed atoms or this is the first timestep read */
    /* then we read all coordinates normally.                         */

    /* read the charmm periodic cell information */
    /* XXX this too should be read together with the other items in a */
    /*     single fio_readv() call in order to prevent lots of extra  */
    /*     kernel/user context switches.                              */
    ret_val = read_charmm_extrablock(fd, charmm, reverseEndian, unitcell);
    if (ret_val) return ret_val;

    /* setup the I/O vector for the call to fio_readv() */
    iov[0].iov_base = (fio_caddr_t) &tmpbuf[0]; /* read format integer    */
    iov[0].iov_len  = sizeof(int);

    iov[1].iov_base = (fio_caddr_t) X;          /* read X coordinates     */
    iov[1].iov_len  = sizeof(float)*N;

    iov[2].iov_base = (fio_caddr_t) &tmpbuf[1]; /* read 2 format integers */
    iov[2].iov_len  = sizeof(int) * 2;

    iov[3].iov_base = (fio_caddr_t) Y;          /* read Y coordinates     */
    iov[3].iov_len  = sizeof(float)*N;

    iov[4].iov_base = (fio_caddr_t) &tmpbuf[3]; /* read 2 format integers */
    iov[4].iov_len  = sizeof(int) * 2;

    iov[5].iov_base = (fio_caddr_t) Z;          /* read Y coordinates     */
    iov[5].iov_len  = sizeof(float)*N;

    iov[6].iov_base = (fio_caddr_t) &tmpbuf[5]; /* read format integer    */
    iov[6].iov_len  = sizeof(int);

    readlen = fio_readv(fd, &iov[0], 7);

    if (readlen != (6*sizeof(int) + 3*N*sizeof(float)))
      return DCD_BADREAD;

    /* convert endianism if necessary */
    if (reverseEndian) {
      swap4_aligned(&tmpbuf[0], 6);
      swap4_aligned(X, N);
      swap4_aligned(Y, N);
      swap4_aligned(Z, N);
    }

    /* double-check the fortran format size values for safety */
    for (i=0; i<6; i++) {
      if (tmpbuf[i] != sizeof(float)*N) return DCD_BADFORMAT;
    }

    /* copy fixed atom coordinates into fixedcoords array if this was the */
    /* first timestep, to be used from now on.  We just copy all atoms.   */
    if (num_fixed && first) {
      memcpy(fixedcoords, X, N*sizeof(float));
      memcpy(fixedcoords+N, Y, N*sizeof(float));
      memcpy(fixedcoords+2*N, Z, N*sizeof(float));
    }

    /* read in the optional charmm 4th array */
    /* XXX this too should be read together with the other items in a */
    /*     single fio_readv() call in order to prevent lots of extra  */
    /*     kernel/user context switches.                              */
    ret_val = read_charmm_4dim(fd, charmm, reverseEndian);
    if (ret_val) return ret_val;
  } else {
    /* if there are fixed atoms, and this isn't the first frame, then we */
    /* only read in the non-fixed atoms for all subsequent timesteps.    */
    ret_val = read_charmm_extrablock(fd, charmm, reverseEndian, unitcell);
    if (ret_val) return ret_val;
    ret_val = read_fixed_atoms(fd, N, N-num_fixed, indexes, reverseEndian,
                               fixedcoords, fixedcoords+3*N, X);
    if (ret_val) return ret_val;
    ret_val = read_fixed_atoms(fd, N, N-num_fixed, indexes, reverseEndian,
                               fixedcoords+N, fixedcoords+3*N, Y);
    if (ret_val) return ret_val;
    ret_val = read_fixed_atoms(fd, N, N-num_fixed, indexes, reverseEndian,
                               fixedcoords+2*N, fixedcoords+3*N, Z);
    if (ret_val) return ret_val;
    ret_val = read_charmm_4dim(fd, charmm, reverseEndian);
    if (ret_val) return ret_val;
  }

  return DCD_SUCCESS;
}

static int skip_dcdstep(fio_fd fd, int natoms, int nfixed, int charmm) {
  int seekoffset = 0;

  /* Skip charmm extra block */
  if ((charmm & DCD_IS_CHARMM) && (charmm & DCD_HAS_EXTRA_BLOCK)) {
    seekoffset += 4 + 48 + 4;
  }

  /* For each atom set, seek past an int, the free atoms, and another int. */
  seekoffset += 3 * (2 + natoms - nfixed) * 4;

  /* Assume that charmm 4th dim is the same size as the other three. */
  if ((charmm & DCD_IS_CHARMM) && (charmm & DCD_HAS_4DIMS)) {
    seekoffset += (2 + natoms - nfixed) * 4;
  }
 
  if (fio_fseek(fd, seekoffset, FIO_SEEK_CUR)) return DCD_BADEOF;

  return DCD_SUCCESS;
}

#define NFILE_POS 8L
#define NSTEP_POS 20L

static int write_dcdstep(fio_fd fd, int curframe, int curstep, int N, 
                  const float *X, const float *Y, const float *Z, 
                  const double *unitcell, int charmm) {
  int out_integer;

  if (charmm) {
    /* write out optional unit cell */
    if (unitcell != NULL) {
      out_integer = 48; /* 48 bytes (6 doubles) */
      fio_write_int32(fd, out_integer);
      WRITE(fd, unitcell, out_integer);
      fio_write_int32(fd, out_integer);
    }
  }

  /* write out coordinates */
  out_integer = N*4; /* N*4 bytes per X/Y/Z array (N floats per array) */
  fio_write_int32(fd, out_integer);
  WRITE(fd, X, out_integer);
  fio_write_int32(fd, out_integer);
  fio_write_int32(fd, out_integer);
  WRITE(fd, Y, out_integer);
  fio_write_int32(fd, out_integer);
  fio_write_int32(fd, out_integer);
  WRITE(fd, Z, out_integer);
  fio_write_int32(fd, out_integer);

  /* update the DCD header information */
  fio_fseek(fd, NFILE_POS, FIO_SEEK_SET);
  fio_write_int32(fd, curframe);
  fio_fseek(fd, NSTEP_POS, FIO_SEEK_SET);
  fio_write_int32(fd, curstep);
  fio_fseek(fd, 0, FIO_SEEK_END);

  return DCD_SUCCESS;
}

static int write_dcdheader(fio_fd fd, const char *remarks, int N, 
                    int ISTART, int NSAVC, double DELTA, int with_unitcell,
                    int charmm) {
  int out_integer;
  float out_float;
  char title_string[200];
  time_t cur_time;
  struct tm *tmbuf;
  char time_str[81];

  out_integer = 84;
  WRITE(fd, (char *) & out_integer, sizeof(int));
  strcpy(title_string, "CORD");
  WRITE(fd, title_string, 4);
  fio_write_int32(fd, 0);      /* Number of frames in file, none written yet   */
  fio_write_int32(fd, ISTART); /* Starting timestep                            */
  fio_write_int32(fd, NSAVC);  /* Timesteps between frames written to the file */
  fio_write_int32(fd, 0);      /* Number of timesteps in simulation            */
  fio_write_int32(fd, 0);      /* NAMD writes NSTEP or ISTART - NSAVC here?    */
  fio_write_int32(fd, 0);
  fio_write_int32(fd, 0);
  fio_write_int32(fd, 0);
  fio_write_int32(fd, 0);
  if (charmm) {
    out_float = DELTA;
    WRITE(fd, (char *) &out_float, sizeof(float));
    if (with_unitcell) {
      fio_write_int32(fd, 1);
    } else {
      fio_write_int32(fd, 0);
    }
  } else {
    WRITE(fd, (char *) &DELTA, sizeof(double));
  }
  fio_write_int32(fd, 0);
  fio_write_int32(fd, 0);
  fio_write_int32(fd, 0);
  fio_write_int32(fd, 0);
  fio_write_int32(fd, 0);
  fio_write_int32(fd, 0);
  fio_write_int32(fd, 0);
  fio_write_int32(fd, 0);
  if (charmm) {
    fio_write_int32(fd, 24); /* Pretend to be Charmm version 24 */
  } else {
    fio_write_int32(fd, 0);
  }
  fio_write_int32(fd, 84);
  fio_write_int32(fd, 164);
  fio_write_int32(fd, 2);

  strncpy(title_string, remarks, 80);
  title_string[79] = '\0';
  WRITE(fd, title_string, 80);

  cur_time=time(NULL);
  tmbuf=localtime(&cur_time);
  strftime(time_str, 80, "REMARKS Created %d %B, %Y at %R", tmbuf);
  WRITE(fd, time_str, 80);

  fio_write_int32(fd, 164);
  fio_write_int32(fd, 4);
  fio_write_int32(fd, N);
  fio_write_int32(fd, 4);

  return DCD_SUCCESS;
}

static void close_dcd_read(int *indexes, float *fixedcoords) {
  free(indexes);
  free(fixedcoords);
}

#endif

