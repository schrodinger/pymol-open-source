/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_dcdplugin
#define STATIC_PLUGIN 1

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
 *      $RCSfile: dcdplugin.c,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.72 $       $Date: 2009/01/14 22:16:52 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *   Code for reading and writing CHARMM, NAMD, and X-PLOR format 
 *   molecular dynamic trajectory files.
 *
 * TODO:
 *   Integrate improvements from the NAMD source tree
 *    - NAMD's writer code has better type-correctness for the sizes
 *      of "int".  NAMD uses "int32" explicitly, which is required on
 *      machines like the T3E.  VMD's version of the code doesn't do that
 *      presently.
 *
 *  Try various alternative I/O API options:
 *   - use mmap(), with read-once flags
 *   - use O_DIRECT open mode on new revs of Linux kernel 
 *   - use directio() call on a file descriptor to enable on Solaris
 *   - use aio_open()/read()/write()
 *   - use readv/writev() etc.
 *
 *  Standalone test binary compilation flags:
 *  cc -fast -xarch=v9a -I../../include -DTEST_DCDPLUGIN dcdplugin.c \
 *    -o ~/bin/readdcd -lm
 *
 *  Profiling flags:
 *  cc -xpg -fast -xarch=v9a -g -I../../include -DTEST_DCDPLUGIN dcdplugin.c \
 *    -o ~/bin/readdcd -lm
 *
 ***************************************************************************/

#include "largefiles.h"   /* platform dependent 64-bit file I/O defines */

#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "endianswap.h"
#include "fastio.h"
#include "molfile_plugin.h"

#ifndef M_PI_2
#define M_PI_2 1.57079632679489661922
#endif

#define RECSCALE32BIT 1
#define RECSCALE64BIT 2
#define RECSCALEMAX   2

typedef struct {
  fio_fd fd;
  int natoms;
  int nsets;
  int setsread;
  int istart;
  int nsavc;
  double delta;
  int nfixed;
  float *x, *y, *z;
  int *freeind;
  float *fixedcoords;
  int reverse;
  int charmm;  
  int first;
  int with_unitcell;
} dcdhandle;

/* Define error codes that may be returned by the DCD routines */
#define DCD_SUCCESS      0  /* No problems                     */
#define DCD_EOF         -1  /* Normal EOF                      */
#define DCD_DNE         -2  /* DCD file does not exist         */
#define DCD_OPENFAILED  -3  /* Open of DCD file failed         */
#define DCD_BADREAD     -4  /* read call on DCD file failed    */
#define DCD_BADEOF      -5  /* premature EOF found in DCD file */
#define DCD_BADFORMAT   -6  /* format of DCD file is wrong     */
#define DCD_FILEEXISTS  -7  /* output file already exists      */
#define DCD_BADMALLOC   -8  /* malloc failed                   */
#define DCD_BADWRITE    -9  /* write call on DCD file failed   */

/* Define feature flags for this DCD file */
#define DCD_IS_XPLOR        0x00
#define DCD_IS_CHARMM       0x01
#define DCD_HAS_4DIMS       0x02
#define DCD_HAS_EXTRA_BLOCK 0x04
#define DCD_HAS_64BIT_REC   0x08

/* defines used by write_dcdstep */
#define NFILE_POS 8L
#define NSTEP_POS 20L

/* READ Macro to make porting easier */
#define READ(fd, buf, size)  fio_fread(((void *) buf), (size), 1, (fd))

/* WRITE Macro to make porting easier */
#define WRITE(fd, buf, size) fio_fwrite(((void *) buf), (size), 1, (fd))

/* XXX This is broken - fread never returns -1 */
#define CHECK_FREAD(X, msg) if (X==-1) { return(DCD_BADREAD); }
#define CHECK_FEOF(X, msg)  if (X==0)  { return(DCD_BADEOF); }


/* print DCD error in a human readable way */
static void print_dcderror(const char *func, int errcode) {
  const char *errstr;

  switch (errcode) {
    case DCD_EOF:         errstr = "end of file"; break;
    case DCD_DNE:         errstr = "file not found"; break;
    case DCD_OPENFAILED:  errstr = "file open failed"; break;
    case DCD_BADREAD:     errstr = "error during read"; break;
    case DCD_BADEOF:      errstr = "premature end of file"; break;
    case DCD_BADFORMAT:   errstr = "corruption or unrecognized file structure"; break;
    case DCD_FILEEXISTS:  errstr = "output file already exists"; break;
    case DCD_BADMALLOC:   errstr = "memory allocation failed"; break;
    case DCD_BADWRITE:    errstr = "error during write"; break;
    case DCD_SUCCESS:     
    default:
      errstr = "no error";
      break;
  } 
  printf("dcdplugin) %s: %s\n", func, errstr); 
}


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
static int read_dcdheader(fio_fd fd, int *N, int *NSET, int *ISTART, 
                   int *NSAVC, double *DELTA, int *NAMNF, 
                   int **FREEINDEXES, float **fixedcoords, int *reverseEndian, 
                   int *charmm)
{
  unsigned int input_integer[2];  /* buffer space */
  int i, ret_val, rec_scale;
  union hdrbuf_union {
  char charvalue[84];    /* char buffer used to store header */
  int intvalue;
  };
  union hdrbuf_union hdrbuf;
  int NTITLE;
  int dcdcordmagic;
  char *corp = (char *) &dcdcordmagic;

  /* coordinate dcd file magic string 'CORD' */
  corp[0] = 'C';
  corp[1] = 'O';
  corp[2] = 'R';
  corp[3] = 'D';

  /* First thing in the file should be an 84.
   * some 64-bit compiles have a 64-bit record length indicator,
   * so we have to read two ints and check in a more complicated 
   * way. :-( */
  ret_val = READ(fd, input_integer, 2*sizeof(unsigned int));
  CHECK_FREAD(ret_val, "reading first int from dcd file");
  CHECK_FEOF(ret_val, "reading first int from dcd file");

  /* Check magic number in file header and determine byte order*/
  if ((input_integer[0]+input_integer[1]) == 84) {
    *reverseEndian=0;
    rec_scale=RECSCALE64BIT;
    printf("dcdplugin) detected CHARMM -i8 64-bit DCD file of native endianness\n");
  } else if (input_integer[0] == 84 && input_integer[1] == dcdcordmagic) {
    *reverseEndian=0;
    rec_scale=RECSCALE32BIT;
    printf("dcdplugin) detected standard 32-bit DCD file of native endianness\n");
  } else {
    /* now try reverse endian */
    swap4_aligned(input_integer, 2); /* will have to unswap magic if 32-bit */
    if ((input_integer[0]+input_integer[1]) == 84) {
      *reverseEndian=1;
      rec_scale=RECSCALE64BIT;
      printf("dcdplugin) detected CHARMM -i8 64-bit DCD file of opposite endianness\n");
    } else {
      swap4_aligned(&input_integer[1], 1); /* unswap magic (see above) */
      if (input_integer[0] == 84 && input_integer[1] == dcdcordmagic) {
        *reverseEndian=1;
        rec_scale=RECSCALE32BIT;
        printf("dcdplugin) detected standard 32-bit DCD file of opposite endianness\n");
      } else {
        /* not simply reversed endianism or -i8, something rather more evil */
        printf("dcdplugin) unrecognized DCD header:\n");
        printf("dcdplugin)   [0]: %10d  [1]: %10d\n", input_integer[0], input_integer[1]);
        printf("dcdplugin)   [0]: 0x%08x  [1]: 0x%08x\n", input_integer[0], input_integer[1]);
        return DCD_BADFORMAT;

      }
    }
  }

  /* check for magic string, in case of long record markers */
  if (rec_scale == RECSCALE64BIT) { 
    ret_val = READ(fd, input_integer, sizeof(unsigned int));
    if (input_integer[0] != dcdcordmagic) {
      printf("dcdplugin) failed to find CORD magic in CHARMM -i8 64-bit DCD file\n");
      return DCD_BADFORMAT;
    }
  }

  /* Buffer the entire header for random access */
  ret_val = READ(fd, hdrbuf.charvalue, 80);
  CHECK_FREAD(ret_val, "buffering header");
  CHECK_FEOF(ret_val, "buffering header");

  /* CHARMm-genereate DCD files set the last integer in the     */
  /* header, which is unused by X-PLOR, to its version number.  */
  /* Checking if this is nonzero tells us this is a CHARMm file */
  /* and to look for other CHARMm flags.                        */
  if (*((int *) (hdrbuf.intvalue + 76)) != 0) {
    (*charmm) = DCD_IS_CHARMM;
    if (*((int *) (hdrbuf.intvalue + 40)) != 0)
      (*charmm) |= DCD_HAS_EXTRA_BLOCK;

    if (*((int *) (hdrbuf.intvalue + 44)) == 1)
      (*charmm) |= DCD_HAS_4DIMS;

    if (rec_scale == RECSCALE64BIT)
      (*charmm) |= DCD_HAS_64BIT_REC;
  
  } else {
    (*charmm) = DCD_IS_XPLOR; /* must be an X-PLOR format DCD file */
  }

  if (*charmm & DCD_IS_CHARMM) {
    /* CHARMM and NAMD versions 2.1b1 and later */
    printf("dcdplugin) CHARMM format DCD file (also NAMD 2.1 and later)\n");
  } else {
    /* CHARMM and NAMD versions prior to 2.1b1  */
    printf("dcdplugin) X-PLOR format DCD file (also NAMD 2.0 and earlier)\n");
  }

  /* Store the number of sets of coordinates (NSET) */
  (*NSET) = *((int *)(hdrbuf.intvalue)); 
  if (*reverseEndian) swap4_unaligned(NSET, 1);

  /* Store ISTART, the starting timestep */
  (*ISTART) = *((int *) (hdrbuf.intvalue + 4));
  if (*reverseEndian) swap4_unaligned(ISTART, 1);

  /* Store NSAVC, the number of timesteps between dcd saves */
  (*NSAVC) = *((int *) (hdrbuf.intvalue + 8));
  if (*reverseEndian) swap4_unaligned(NSAVC, 1);

  /* Store NAMNF, the number of fixed atoms */
  (*NAMNF) = *((int *) (hdrbuf.intvalue + 32));
  if (*reverseEndian) swap4_unaligned(NAMNF, 1);

  /* Read in the timestep, DELTA */
  /* Note: DELTA is stored as a double with X-PLOR but as a float with CHARMm */
  if ((*charmm) & DCD_IS_CHARMM) {
    float ftmp;
    ftmp = *((float *)(hdrbuf.intvalue+36)); /* is this safe on Alpha? */
    if (*reverseEndian)
      swap4_aligned(&ftmp, 1);

    *DELTA = (double)ftmp;
  } else {
    (*DELTA) = *((double *)(hdrbuf.intvalue + 36));
    if (*reverseEndian) swap8_unaligned(DELTA, 1);
  }

  /* Get the end size of the first block */
  ret_val = READ(fd, input_integer, rec_scale*sizeof(int));
  CHECK_FREAD(ret_val, "reading second 84 from dcd file");
  CHECK_FEOF(ret_val, "reading second 84 from dcd file");
  if (*reverseEndian) swap4_aligned(input_integer, rec_scale);

  if (rec_scale == RECSCALE64BIT) {
    if ((input_integer[0]+input_integer[1]) != 84) {
      return DCD_BADFORMAT;
    }
  } else {
    if (input_integer[0] != 84) {
      return DCD_BADFORMAT;
    }
  }
  
  /* Read in the size of the next block */
  input_integer[1] = 0;
  ret_val = READ(fd, input_integer, rec_scale*sizeof(int));
  CHECK_FREAD(ret_val, "reading size of title block");
  CHECK_FEOF(ret_val, "reading size of title block");
  if (*reverseEndian) swap4_aligned(input_integer, rec_scale);

  if ((((input_integer[0]+input_integer[1])-4) % 80) == 0) {
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
    ret_val = READ(fd, input_integer, rec_scale*sizeof(int));
    CHECK_FREAD(ret_val, "reading size of title block");
    CHECK_FEOF(ret_val, "reading size of title block");
  } else {
    return DCD_BADFORMAT;
  }

  /* Read in an integer '4' */
  input_integer[1] = 0;
  ret_val = READ(fd, input_integer, rec_scale*sizeof(int));
  
  CHECK_FREAD(ret_val, "reading a '4'");
  CHECK_FEOF(ret_val, "reading a '4'");
  if (*reverseEndian) swap4_aligned(input_integer, rec_scale);

  if ((input_integer[0]+input_integer[1]) != 4) {
    return DCD_BADFORMAT;
  }

  /* Read in the number of atoms */
  ret_val = READ(fd, N, sizeof(int));
  CHECK_FREAD(ret_val, "reading number of atoms");
  CHECK_FEOF(ret_val, "reading number of atoms");
  if (*reverseEndian) swap4_aligned(N, 1);

  /* Read in an integer '4' */
  input_integer[1] = 0;
  ret_val = READ(fd, input_integer, rec_scale*sizeof(int));
  CHECK_FREAD(ret_val, "reading a '4'");
  CHECK_FEOF(ret_val, "reading a '4'");
  if (*reverseEndian) swap4_aligned(input_integer, rec_scale);

  if ((input_integer[0]+input_integer[1]) != 4) {
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
    input_integer[1]=0;
    ret_val = READ(fd, input_integer, rec_scale*sizeof(int));
    CHECK_FREAD(ret_val, "reading size of index array");
    CHECK_FEOF(ret_val, "reading size of index array");
    if (*reverseEndian) swap4_aligned(input_integer, rec_scale);

    if ((input_integer[0]+input_integer[1]) != ((*N)-(*NAMNF))*4) {
      return DCD_BADFORMAT;
    }

    ret_val = READ(fd, (*FREEINDEXES), ((*N)-(*NAMNF))*sizeof(int));
    CHECK_FREAD(ret_val, "reading size of index array");
    CHECK_FEOF(ret_val, "reading size of index array");

    if (*reverseEndian)
      swap4_aligned((*FREEINDEXES), ((*N)-(*NAMNF)));

    input_integer[1]=0;
    ret_val = READ(fd, input_integer, rec_scale*sizeof(int));
    CHECK_FREAD(ret_val, "reading size of index array");
    CHECK_FEOF(ret_val, "reading size of index array");
    if (*reverseEndian) swap4_aligned(input_integer, rec_scale);

    if ((input_integer[0]+input_integer[1]) != ((*N)-(*NAMNF))*4) {
      return DCD_BADFORMAT;
    }
  }

  return DCD_SUCCESS;
}

static int read_charmm_extrablock(fio_fd fd, int charmm, int reverseEndian,
                                  float *unitcell) {
  int i, input_integer[2], rec_scale;

  if (charmm & DCD_HAS_64BIT_REC) {
    rec_scale = RECSCALE64BIT;
  } else {
    rec_scale = RECSCALE32BIT;
  }

  if ((charmm & DCD_IS_CHARMM) && (charmm & DCD_HAS_EXTRA_BLOCK)) {
    /* Leading integer must be 48 */
    input_integer[1] = 0;
    if (fio_fread(input_integer, sizeof(int), rec_scale, fd) != rec_scale)
      return DCD_BADREAD; 
    if (reverseEndian) swap4_aligned(input_integer, rec_scale);
    if ((input_integer[0]+input_integer[1]) == 48) {
      double tmp[6];
      if (fio_fread(tmp, 48, 1, fd) != 1) return DCD_BADREAD;
      if (reverseEndian) 
        swap8_aligned(tmp, 6);
      for (i=0; i<6; i++) unitcell[i] = (float)tmp[i];
    } else {
      /* unrecognized block, just skip it */
      if (fio_fseek(fd, (input_integer[0]+input_integer[1]), FIO_SEEK_CUR)) return DCD_BADREAD;
    }
    if (fio_fread(input_integer, sizeof(int), rec_scale, fd) != rec_scale) return DCD_BADREAD; 
  } 

  return DCD_SUCCESS;
}

static int read_fixed_atoms(fio_fd fd, int N, int num_free, const int *indexes,
                            int reverseEndian, const float *fixedcoords, 
                            float *freeatoms, float *pos, int charmm) {
  int i, input_integer[2], rec_scale;
  
  if(charmm & DCD_HAS_64BIT_REC) {
    rec_scale=RECSCALE64BIT;
  } else {
    rec_scale=RECSCALE32BIT;
  }
  
  /* Read leading integer */
  input_integer[1]=0;
  if (fio_fread(input_integer, sizeof(int), rec_scale, fd) != rec_scale) return DCD_BADREAD;
  if (reverseEndian) swap4_aligned(input_integer, rec_scale);
  if ((input_integer[0]+input_integer[1]) != 4*num_free) return DCD_BADFORMAT;
  
  /* Read free atom coordinates */
  if (fio_fread(freeatoms, 4*num_free, 1, fd) != 1) return DCD_BADREAD;
  if (reverseEndian)
    swap4_aligned(freeatoms, num_free);

  /* Copy fixed and free atom coordinates into position buffer */
  memcpy(pos, fixedcoords, 4*N);
  for (i=0; i<num_free; i++)
    pos[indexes[i]-1] = freeatoms[i];

  /* Read trailing integer */ 
  input_integer[1]=0;
  if (fio_fread(input_integer, sizeof(int), rec_scale, fd) != rec_scale) return DCD_BADREAD;
  if (reverseEndian) swap4_aligned(input_integer, rec_scale);
  if ((input_integer[0]+input_integer[1]) != 4*num_free) return DCD_BADFORMAT;

  return DCD_SUCCESS;
}
  
static int read_charmm_4dim(fio_fd fd, int charmm, int reverseEndian) {
  int input_integer[2],rec_scale;

  if (charmm & DCD_HAS_64BIT_REC) {
    rec_scale=RECSCALE64BIT;
  } else {
    rec_scale=RECSCALE32BIT;
  }
    
  /* If this is a CHARMm file and contains a 4th dimension block, */
  /* we must skip past it to avoid problems                       */
  if ((charmm & DCD_IS_CHARMM) && (charmm & DCD_HAS_4DIMS)) {
    input_integer[1]=0;
    if (fio_fread(input_integer, sizeof(int), rec_scale, fd) != rec_scale) return DCD_BADREAD;  
    if (reverseEndian) swap4_aligned(input_integer, rec_scale);
    if (fio_fseek(fd, (input_integer[0]+input_integer[1]), FIO_SEEK_CUR)) return DCD_BADREAD;
    if (fio_fread(input_integer, sizeof(int), rec_scale, fd) != rec_scale) return DCD_BADREAD;  
  }

  return DCD_SUCCESS;
}

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
static int read_dcdstep(fio_fd fd, int N, float *X, float *Y, float *Z, 
                        float *unitcell, int num_fixed,
                        int first, int *indexes, float *fixedcoords, 
                        int reverseEndian, int charmm) {
  int ret_val, rec_scale;   /* Return value from read */
  
  if (charmm & DCD_HAS_64BIT_REC) {
    rec_scale=RECSCALE64BIT;
  } else {
    rec_scale=RECSCALE32BIT;
  }
  
  if ((num_fixed==0) || first) {
    /* temp storage for reading formatting info */
    /* note: has to be max size we'll ever use  */
    int tmpbuf[6*RECSCALEMAX]; 

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
    iov[0].iov_len  = rec_scale*sizeof(int);

    iov[1].iov_base = (fio_caddr_t) X;          /* read X coordinates     */
    iov[1].iov_len  = sizeof(float)*N;

    iov[2].iov_base = (fio_caddr_t) &tmpbuf[1*rec_scale]; /* read 2 format integers */
    iov[2].iov_len  = rec_scale*sizeof(int) * 2;

    iov[3].iov_base = (fio_caddr_t) Y;          /* read Y coordinates     */
    iov[3].iov_len  = sizeof(float)*N;

    iov[4].iov_base = (fio_caddr_t) &tmpbuf[3*rec_scale]; /* read 2 format integers */
    iov[4].iov_len  = rec_scale*sizeof(int) * 2;

    iov[5].iov_base = (fio_caddr_t) Z;          /* read Y coordinates     */
    iov[5].iov_len  = sizeof(float)*N;

    iov[6].iov_base = (fio_caddr_t) &tmpbuf[5*rec_scale]; /* read format integer    */
    iov[6].iov_len  = rec_scale*sizeof(int);

    readlen = fio_readv(fd, &iov[0], 7);

    if (readlen != (rec_scale*6*sizeof(int) + 3*N*sizeof(float)))
      return DCD_BADREAD;

    /* convert endianism if necessary */
    if (reverseEndian) {
      swap4_aligned(&tmpbuf[0], rec_scale*6);
      swap4_aligned(X, N);
      swap4_aligned(Y, N);
      swap4_aligned(Z, N);
    }

    /* double-check the fortran format size values for safety */
    if(rec_scale == 1) {
      for (i=0; i<6; i++) {
        if (tmpbuf[i] != sizeof(float)*N) return DCD_BADFORMAT;
      }
    } else {
      for (i=0; i<6; i++) {
          if ((tmpbuf[2*i]+tmpbuf[2*i+1]) != sizeof(float)*N) return DCD_BADFORMAT;
      }
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
                               fixedcoords, fixedcoords+3*N, X, charmm);
    if (ret_val) return ret_val;
    ret_val = read_fixed_atoms(fd, N, N-num_fixed, indexes, reverseEndian,
                               fixedcoords+N, fixedcoords+3*N, Y, charmm);
    if (ret_val) return ret_val;
    ret_val = read_fixed_atoms(fd, N, N-num_fixed, indexes, reverseEndian,
                               fixedcoords+2*N, fixedcoords+3*N, Z, charmm);
    if (ret_val) return ret_val;
    ret_val = read_charmm_4dim(fd, charmm, reverseEndian);
    if (ret_val) return ret_val;
  }

  return DCD_SUCCESS;
}


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
static int skip_dcdstep(fio_fd fd, int natoms, int nfixed, int charmm) {
  
  int seekoffset = 0;
  int rec_scale;

  if (charmm & DCD_HAS_64BIT_REC) {
    rec_scale=RECSCALE64BIT;
  } else {
    rec_scale=RECSCALE32BIT;
  }

  /* Skip charmm extra block */
  if ((charmm & DCD_IS_CHARMM) && (charmm & DCD_HAS_EXTRA_BLOCK)) {
    seekoffset += 4*rec_scale + 48 + 4*rec_scale;
  }

  /* For each atom set, seek past an int, the free atoms, and another int. */
  seekoffset += 3 * (2*rec_scale + natoms - nfixed) * 4;

  /* Assume that charmm 4th dim is the same size as the other three. */
  if ((charmm & DCD_IS_CHARMM) && (charmm & DCD_HAS_4DIMS)) {
    seekoffset += (2*rec_scale + natoms - nfixed) * 4;
  }
 
  if (fio_fseek(fd, seekoffset, FIO_SEEK_CUR)) return DCD_BADEOF;

  return DCD_SUCCESS;
}


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
static int write_dcdstep(fio_fd fd, int curframe, int curstep, int N, 
                  const float *X, const float *Y, const float *Z, 
                  const double *unitcell, int charmm) {
  int out_integer;

  if (charmm) {
    /* write out optional unit cell */
    if (unitcell != NULL) {
      out_integer = 48; /* 48 bytes (6 floats) */
      fio_write_int32(fd, out_integer);
      WRITE(fd, unitcell, out_integer);
      fio_write_int32(fd, out_integer);
    }
  }

  /* write out coordinates */
  out_integer = N*4; /* N*4 bytes per X/Y/Z array (N floats per array) */
  fio_write_int32(fd, out_integer);
  if (fio_fwrite((void *) X, out_integer, 1, fd) != 1) return DCD_BADWRITE;
  fio_write_int32(fd, out_integer);
  fio_write_int32(fd, out_integer);
  if (fio_fwrite((void *) Y, out_integer, 1, fd) != 1) return DCD_BADWRITE;
  fio_write_int32(fd, out_integer);
  fio_write_int32(fd, out_integer);
  if (fio_fwrite((void *) Z, out_integer, 1, fd) != 1) return DCD_BADWRITE;
  fio_write_int32(fd, out_integer);

  /* update the DCD header information */
  fio_fseek(fd, NFILE_POS, FIO_SEEK_SET);
  fio_write_int32(fd, curframe);
  fio_fseek(fd, NSTEP_POS, FIO_SEEK_SET);
  fio_write_int32(fd, curstep);
  fio_fseek(fd, 0, FIO_SEEK_END);

  return DCD_SUCCESS;
}

/*
 * Write a header for a new dcd file
 * Input: fd - file struct opened for binary writing
 *        remarks - string to be put in the remarks section of the header.  
 *                  The string will be truncated to 70 characters.
 *        natoms, istart, nsavc, delta - see comments in read_dcdheader
 * Output: 0 on success, negative error code on failure.
 * Side effects: Header information is written to the dcd file.
 */
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
    fio_write_int32(fd, 24); /* Pretend to be CHARMM version 24 */
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


/*
 * clean up dcd data
 * Input: nfixed, freeind - elements as returned by read_dcdheader
 * Output: None
 * Side effects: Space pointed to by freeind is freed if necessary.
 */
static void close_dcd_read(int *indexes, float *fixedcoords) {
  free(indexes);
  free(fixedcoords);
}





static void *open_dcd_read(const char *path, const char *filetype, 
    int *natoms) {
  dcdhandle *dcd;
  fio_fd fd;
  int rc;
  struct stat stbuf;

  if (!path) return NULL;

  /* See if the file exists, and get its size */
  memset(&stbuf, 0, sizeof(struct stat));
  if (stat(path, &stbuf)) {
    printf("dcdplugin) Could not access file '%s'.\n", path);
    return NULL;
  }

  if (fio_open(path, FIO_READ, &fd) < 0) {
    printf("dcdplugin) Could not open file '%s' for reading.\n", path);
    return NULL;
  }

  dcd = (dcdhandle *)malloc(sizeof(dcdhandle));
  memset(dcd, 0, sizeof(dcdhandle));
  dcd->fd = fd;

  if ((rc = read_dcdheader(dcd->fd, &dcd->natoms, &dcd->nsets, &dcd->istart, 
         &dcd->nsavc, &dcd->delta, &dcd->nfixed, &dcd->freeind, 
         &dcd->fixedcoords, &dcd->reverse, &dcd->charmm))) {
    print_dcderror("read_dcdheader", rc);
    fio_fclose(dcd->fd);
    free(dcd);
    return NULL;
  }

  /*
   * Check that the file is big enough to really hold the number of sets
   * it claims to have.  Then we'll use nsets to keep track of where EOF
   * should be.
   */
  {
    fio_size_t ndims, firstframesize, framesize, extrablocksize;
    fio_size_t trjsize, filesize, curpos;
    int newnsets;

    extrablocksize = dcd->charmm & DCD_HAS_EXTRA_BLOCK ? 48 + 8 : 0;
    ndims = dcd->charmm & DCD_HAS_4DIMS ? 4 : 3;
    firstframesize = (dcd->natoms+2) * ndims * sizeof(float) + extrablocksize;
    framesize = (dcd->natoms-dcd->nfixed+2) * ndims * sizeof(float) 
      + extrablocksize;

    /* 
     * It's safe to use ftell, even though ftell returns a long, because the 
     * header size is < 4GB.
     */

    curpos = fio_ftell(dcd->fd); /* save current offset (end of header) */

#if defined(_MSC_VER) && defined(FASTIO_NATIVEWIN32)
    /* the stat() call is not 64-bit savvy on Windows             */
    /* so we have to use the fastio fseek/ftell routines for this */
    /* until we add a portable filesize routine for this purpose  */
    fio_fseek(dcd->fd, 0, FIO_SEEK_END);       /* seek to end of file */
    filesize = fio_ftell(dcd->fd);
    fio_fseek(dcd->fd, curpos, FIO_SEEK_SET);  /* return to end of header */
#else
    filesize = stbuf.st_size; /* this works ok on Unix machines */
#endif
    trjsize = filesize - curpos - firstframesize;
    if (trjsize < 0) {
      printf("dcdplugin) file '%s' appears to contain no timesteps.\n", path);
      fio_fclose(dcd->fd);
      free(dcd);
      return NULL;
    }

    newnsets = trjsize / framesize + 1;

    if (dcd->nsets > 0 && newnsets != dcd->nsets) {
      printf("dcdplugin) Warning: DCD header claims %d frames, file size indicates there are actually %d frames\n", dcd->nsets, newnsets);
    }

    dcd->nsets = newnsets; 
    dcd->setsread = 0;
  }

  dcd->first = 1;
  dcd->x = (float *)malloc(dcd->natoms * sizeof(float));
  dcd->y = (float *)malloc(dcd->natoms * sizeof(float));
  dcd->z = (float *)malloc(dcd->natoms * sizeof(float));
  if (!dcd->x || !dcd->y || !dcd->z) {
    printf("dcdplugin) Unable to allocate space for %d atoms.\n", dcd->natoms);
    if (dcd->x)
      free(dcd->x);
    if (dcd->y)
      free(dcd->y);
    if (dcd->z)
      free(dcd->z);
    fio_fclose(dcd->fd);
    free(dcd);
    return NULL;
  }
  *natoms = dcd->natoms;
  return dcd;
}


static int read_next_timestep(void *v, int natoms, molfile_timestep_t *ts) {
  dcdhandle *dcd;
  int i, j, rc;
  float unitcell[6];
  unitcell[0] = unitcell[2] = unitcell[5] = 1.0f;
  unitcell[1] = unitcell[3] = unitcell[4] = 90.0f;
  dcd = (dcdhandle *)v;

  /* Check for EOF here; that way all EOF's encountered later must be errors */
  if (dcd->setsread == dcd->nsets) return MOLFILE_EOF;
  dcd->setsread++;
  if (!ts) {
    if (dcd->first && dcd->nfixed) {
      /* We can't just skip it because we need the fixed atom coordinates */
      rc = read_dcdstep(dcd->fd, dcd->natoms, dcd->x, dcd->y, dcd->z, 
          unitcell, dcd->nfixed, dcd->first, dcd->freeind, dcd->fixedcoords, 
             dcd->reverse, dcd->charmm);
      dcd->first = 0;
      return rc; /* XXX this needs to be updated */
    }
    dcd->first = 0;
    /* XXX this needs to be changed */
    return skip_dcdstep(dcd->fd, dcd->natoms, dcd->nfixed, dcd->charmm);
  }
  rc = read_dcdstep(dcd->fd, dcd->natoms, dcd->x, dcd->y, dcd->z, unitcell,
             dcd->nfixed, dcd->first, dcd->freeind, dcd->fixedcoords, 
             dcd->reverse, dcd->charmm);
  dcd->first = 0;
  if (rc < 0) {  
    print_dcderror("read_dcdstep", rc);
    return MOLFILE_ERROR;
  }

  /* copy timestep data from plugin-local buffers to VMD's buffer */
  /* XXX 
   *   This code is still the root of all evil.  Just doing this extra copy
   *   cuts the I/O rate of the DCD reader from 728 MB/sec down to
   *   394 MB/sec when reading from a ram filesystem.  
   *   For a physical disk filesystem, the I/O rate goes from 
   *   187 MB/sec down to 122 MB/sec.  Clearly this extra copy has to go.
   */
  {
    int natoms = dcd->natoms;
    float *nts = ts->coords;
    const float *bufx = dcd->x;
    const float *bufy = dcd->y;
    const float *bufz = dcd->z;

    for (i=0, j=0; i<natoms; i++, j+=3) {
      nts[j    ] = bufx[i];
      nts[j + 1] = bufy[i];
      nts[j + 2] = bufz[i];
    }
  }

  ts->A = unitcell[0];
  ts->B = unitcell[2];
  ts->C = unitcell[5];

  if (unitcell[1] >= -1.0 && unitcell[1] <= 1.0 &&
      unitcell[3] >= -1.0 && unitcell[3] <= 1.0 &&
      unitcell[4] >= -1.0 && unitcell[4] <= 1.0) {
    /* This file was generated by CHARMM, or by NAMD > 2.5, with the angle */
    /* cosines of the periodic cell angles written to the DCD file.        */ 
    /* This formulation improves rounding behavior for orthogonal cells    */
    /* so that the angles end up at precisely 90 degrees, unlike acos().   */
    ts->alpha = 90.0 - asin(unitcell[4]) * 90.0 / M_PI_2; /* cosBC */
    ts->beta  = 90.0 - asin(unitcell[3]) * 90.0 / M_PI_2; /* cosAC */
    ts->gamma = 90.0 - asin(unitcell[1]) * 90.0 / M_PI_2; /* cosAB */
  } else {
    /* This file was likely generated by NAMD 2.5 and the periodic cell    */
    /* angles are specified in degrees rather than angle cosines.          */
    ts->alpha = unitcell[4]; /* angle between B and C */
    ts->beta  = unitcell[3]; /* angle between A and C */
    ts->gamma = unitcell[1]; /* angle between A and B */
  }
 
  return MOLFILE_SUCCESS;
}
 

static void close_file_read(void *v) {
  dcdhandle *dcd = (dcdhandle *)v;
  close_dcd_read(dcd->freeind, dcd->fixedcoords);
  fio_fclose(dcd->fd);
  free(dcd->x);
  free(dcd->y);
  free(dcd->z);
  free(dcd); 
}


static void *open_dcd_write(const char *path, const char *filetype, 
    int natoms) {
  dcdhandle *dcd;
  fio_fd fd;
  int rc;
  int istart, nsavc;
  double delta;
  int with_unitcell;
  int charmm;

  if (fio_open(path, FIO_WRITE, &fd) < 0) {
    printf("dcdplugin) Could not open file '%s' for writing\n", path);
    return NULL;
  }

  dcd = (dcdhandle *)malloc(sizeof(dcdhandle));
  memset(dcd, 0, sizeof(dcdhandle));
  dcd->fd = fd;

  istart = 0;             /* starting timestep of DCD file                  */
  nsavc = 1;              /* number of timesteps between written DCD frames */
  delta = 1.0;            /* length of a timestep                           */

  if (getenv("VMDDCDWRITEXPLORFORMAT") != NULL) {
    with_unitcell = 0;      /* no unit cell info */
    charmm = DCD_IS_XPLOR;  /* X-PLOR format */
    printf("dcdplugin) WARNING: Writing DCD file in X-PLOR format, \n");
    printf("dcdplugin) WARNING: unit cell information will be lost!\n");
  } else {
    with_unitcell = 1;      /* contains unit cell infor (Charmm format) */
    charmm = DCD_IS_CHARMM; /* charmm-formatted DCD file                */ 
    if (with_unitcell) 
      charmm |= DCD_HAS_EXTRA_BLOCK;
  }
 
  rc = write_dcdheader(dcd->fd, "Created by DCD plugin", natoms, 
                       istart, nsavc, delta, with_unitcell, charmm);

  if (rc < 0) {
    print_dcderror("write_dcdheader", rc);
    fio_fclose(dcd->fd);
    free(dcd);
    return NULL;
  }

  dcd->natoms = natoms;
  dcd->nsets = 0;
  dcd->istart = istart;
  dcd->nsavc = nsavc;
  dcd->with_unitcell = with_unitcell;
  dcd->charmm = charmm;
  dcd->x = (float *)malloc(natoms * sizeof(float));
  dcd->y = (float *)malloc(natoms * sizeof(float));
  dcd->z = (float *)malloc(natoms * sizeof(float));
  return dcd;
}


static int write_timestep(void *v, const molfile_timestep_t *ts) { 
  dcdhandle *dcd = (dcdhandle *)v;
  int i, rc, curstep;
  float *pos = ts->coords;
  double unitcell[6];
  unitcell[0] = unitcell[2] = unitcell[5] = 1.0f;
  unitcell[1] = unitcell[3] = unitcell[4] = 90.0f;

  /* copy atom coords into separate X/Y/Z arrays for writing */
  for (i=0; i<dcd->natoms; i++) {
    dcd->x[i] = *(pos++); 
    dcd->y[i] = *(pos++); 
    dcd->z[i] = *(pos++); 
  }
  dcd->nsets++;
  curstep = dcd->istart + dcd->nsets * dcd->nsavc;

  unitcell[0] = ts->A;
  unitcell[2] = ts->B;
  unitcell[5] = ts->C;
  unitcell[1] = sin((M_PI_2 / 90.0) * (90.0 - ts->gamma)); /* cosAB */
  unitcell[3] = sin((M_PI_2 / 90.0) * (90.0 - ts->beta));  /* cosAC */
  unitcell[4] = sin((M_PI_2 / 90.0) * (90.0 - ts->alpha)); /* cosBC */

  rc = write_dcdstep(dcd->fd, dcd->nsets, curstep, dcd->natoms, 
                     dcd->x, dcd->y, dcd->z,
                     dcd->with_unitcell ? unitcell : NULL,
                     dcd->charmm);
  if (rc < 0) {
    print_dcderror("write_dcdstep", rc);
    return MOLFILE_ERROR;
  }

  return MOLFILE_SUCCESS;
}

static void close_file_write(void *v) {
  dcdhandle *dcd = (dcdhandle *)v;
  fio_fclose(dcd->fd);
  free(dcd->x);
  free(dcd->y);
  free(dcd->z);
  free(dcd);
}


/*
 * Initialization stuff here
 */
static molfile_plugin_t plugin;

VMDPLUGIN_API int VMDPLUGIN_init() {
  memset(&plugin, 0, sizeof(molfile_plugin_t));
  plugin.abiversion = vmdplugin_ABIVERSION;
  plugin.type = MOLFILE_PLUGIN_TYPE;
  plugin.name = "dcd";
  plugin.prettyname = "CHARMM,NAMD,XPLOR DCD Trajectory";
  plugin.author = "Justin Gullingsrud, John Stone";
  plugin.majorv = 1;
  plugin.minorv = 10;
  plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
  plugin.filename_extension = "dcd";
  plugin.open_file_read = open_dcd_read;
  plugin.read_next_timestep = read_next_timestep;
  plugin.close_file_read = close_file_read;
  plugin.open_file_write = open_dcd_write;
  plugin.write_timestep = write_timestep;
  plugin.close_file_write = close_file_write;
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)(void *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_fini() {
  return VMDPLUGIN_SUCCESS;
}

  
#ifdef TEST_DCDPLUGIN

#include <sys/time.h>

/* get the time of day from the system clock, and store it (in seconds) */
double time_of_day(void) {
#if defined(_MSC_VER)
  double t;

  t = GetTickCount();
  t = t / 1000.0;

  return t;
#else
  struct timeval tm;
  struct timezone tz;

  gettimeofday(&tm, &tz);
  return((double)(tm.tv_sec) + (double)(tm.tv_usec)/1000000.0);
#endif
}

int main(int argc, char *argv[]) {
  molfile_timestep_t timestep;
  void *v;
  dcdhandle *dcd;
  int i, natoms;
  float sizeMB =0.0, totalMB = 0.0;
  double starttime, endtime, totaltime = 0.0;

  while (--argc) {
    ++argv; 
    natoms = 0;
    v = open_dcd_read(*argv, "dcd", &natoms);
    if (!v) {
      fprintf(stderr, "main) open_dcd_read failed for file %s\n", *argv);
      return 1;
    }
    dcd = (dcdhandle *)v;
    sizeMB = ((natoms * 3.0) * dcd->nsets * 4.0) / (1024.0 * 1024.0);
    totalMB += sizeMB; 
    printf("main) file: %s\n", *argv);
    printf("  %d atoms, %d frames, size: %6.1fMB\n", natoms, dcd->nsets, sizeMB);

    starttime = time_of_day();
    timestep.coords = (float *)malloc(3*sizeof(float)*natoms);
    for (i=0; i<dcd->nsets; i++) {
      int rc = read_next_timestep(v, natoms, &timestep);
      if (rc) {
        fprintf(stderr, "error in read_next_timestep on frame %d\n", i);
        return 1;
      }
    }
    endtime = time_of_day();
    close_file_read(v);
    totaltime += endtime - starttime;
    printf("  Time: %5.1f seconds\n", endtime - starttime);
    printf("  Speed: %5.1f MB/sec, %5.1f timesteps/sec\n", sizeMB / (endtime - starttime), (dcd->nsets / (endtime - starttime)));
  }
  printf("Overall Size: %6.1f MB\n", totalMB);
  printf("Overall Time: %6.1f seconds\n", totaltime);
  printf("Overall Speed: %5.1f MB/sec\n", totalMB / totaltime);
  return 0;
}
      
#endif

