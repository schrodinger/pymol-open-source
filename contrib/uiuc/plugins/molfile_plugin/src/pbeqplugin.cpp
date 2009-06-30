/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_pbeqplugin
#define STATIC_PLUGIN 1

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
 *      $RCSfile: pbeqplugin.C,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.5 $       $Date: 2009/04/29 15:45:32 $
 *
 ***************************************************************************/

/* 
 * "unformatted" binary potential map, created by CHARMM PBEQ module
 *
 * Format (fortran): 
 *      INTEGER UNIT
 *      INTEGER NCLX,NCLY,NCLZ
 *      REAL*8  DCEL
 *      REAL*8  XBCEN,YBCEN,ZBCEN,EPSW,EPSP,CONC,TMEMB,ZMEMB,EPSM
 *      REAL*4  PHI(*)
 *C
 *C     Local variable
 *      INTEGER I
 *
 *      WRITE(UNIT) NCLX,NCLY,NCLZ,DCEL,XBCEN,YBCEN,ZBCEN
 *      WRITE(UNIT) EPSW,EPSP,CONC,TMEMB,ZMEMB,EPSM
 *      WRITE(UNIT)(PHI(I),I=1,NCLX*NCLY*NCLZ)
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <string.h>

#if defined(_AIX)
#include <strings.h>
#endif

#include "molfile_plugin.h"
#include "endianswap.h"

typedef struct {
  FILE *fd;
  int nsets;
  int ndata;
  int nclx;
  int ncly;
  int nclz;
  int swap;
  molfile_volumetric_t *vol;
} pbeq_t;


static void *open_pbeq_read(const char *filepath, const char *filetype,
    int *natoms) {
  FILE *fd;
  pbeq_t *pbeq;
  int nclx, ncly, nclz;
  int trash, length;
  double dcel; 
  double xbcen, ybcen, zbcen;
  double epsw, epsp, conc, tmemb, zmemb, epsm;
  int swap=0; 

  fd = fopen(filepath, "rb");
  if (!fd) {
    printf("pbeqplugin) Error opening file %s.\n", filepath);
    return NULL;
  }

  // skip first Fortran length record for
  // WRITE(UNIT) NCLX,NCLY,NCLZ,DCEL,XBCEN,YBCEN,ZBCEN
  if (fread(&length, 4, 1, fd) != 1)
    return NULL;

  if (fread(&nclx, 4, 1, fd) != 1)
    return NULL;
  if (fread(&ncly, 4, 1, fd) != 1)
    return NULL;
  if (fread(&nclz, 4, 1, fd) != 1)
    return NULL;

  // test endianness first
  if (length != 44) {
    swap = 1;
    swap4_aligned(&length, 1);
    if (length != 44) {
      printf("pbeqplugin) length record != 44, unrecognized format (length: %d)\n", length);
      return NULL;
    }
  }

  if (swap) {
    swap4_aligned(&nclx, 1);
    swap4_aligned(&ncly, 1);
    swap4_aligned(&nclz, 1);
  } 

  // this is a risky strategy for detecting endianness, 
  // but it works so far, and the charmm potential maps don't
  // have version numbers, magic numbers, or anything else that we
  // might otherwise use for this purpose.
  if ((nclx > 4000 && ncly > 4000 && nclz > 4000) ||
      (nclx * ncly * nclz < 0)) {
    printf("pbeqplugin) inconclusive byte ordering, bailing out\n");
    return NULL;
  }

  // read the rest of the header
  if (fread(&dcel, 8, 1, fd) != 1) 
    return NULL;
  if (fread(&xbcen, 8, 1, fd) != 1) 
    return NULL;
  if (fread(&ybcen, 8, 1, fd) != 1) 
    return NULL;
  if (fread(&zbcen, 8, 1, fd) != 1) 
    return NULL;

  // skip second Fortran length record for
  // WRITE(UNIT) NCLX,NCLY,NCLZ,DCEL,XBCEN,YBCEN,ZBCEN
  if (fread(&trash, 4, 1, fd) != 1)
    return NULL;

  // skip first Fortran length record for 
  // WRITE(UNIT) EPSW,EPSP,CONC,TMEMB,ZMEMB,EPSM
  if (fread(&trash, 4, 1, fd) != 1)
    return NULL;

  if (fread(&epsw, 8, 1, fd) != 1) 
    return NULL;
  if (fread(&epsp, 8, 1, fd) != 1) 
    return NULL;
  if (fread(&conc, 8, 1, fd) != 1) 
    return NULL;
  if (fread(&tmemb, 8, 1, fd) != 1) 
    return NULL;
  if (fread(&zmemb, 8, 1, fd) != 1) 
    return NULL;
  if (fread(&epsm, 8, 1, fd) != 1) 
    return NULL;

  // skip second Fortran length record for 
  // WRITE(UNIT) EPSW,EPSP,CONC,TMEMB,ZMEMB,EPSM
  if (fread(&trash, 4, 1, fd) != 1)
    return NULL;

  // byte swap the header data if necessary
  if (swap) {
    swap8_aligned(&dcel, 1);
    swap8_aligned(&xbcen, 1);
    swap8_aligned(&ybcen, 1);
    swap8_aligned(&zbcen, 1);
    swap8_aligned(&epsw, 1);
    swap8_aligned(&epsp, 1);
    swap8_aligned(&conc, 1);
    swap8_aligned(&tmemb, 1);
    swap8_aligned(&zmemb, 1);
    swap8_aligned(&epsm, 1);
  }

#if 0
  // print header info for debugging of early versions
  printf("pbeqplugin) nclx:%d nxly:%d nclz:%d\n", nclx, ncly, nclz);
  printf("pbeqplugin) dcel: %f\n", dcel); 
  printf("pbeqplugin) x/y/zbcen: %g %g %g\n", xbcen, ybcen, zbcen); 
  printf("pbeqplugin) epsw/p: %g %g  conc: %g  tmemb: %g zmemb: %g  epsm: %g\n",
         epsw, epsp, conc, tmemb, zmemb, epsm);
#endif

  /* Allocate and initialize the pbeq structure */
  pbeq = new pbeq_t;
  pbeq->fd = fd;
  pbeq->vol = NULL;
  *natoms = MOLFILE_NUMATOMS_NONE;
  pbeq->nsets = 1; /* this file contains only one data set */
  pbeq->ndata = nclx*ncly*nclz;
  pbeq->nclx = nclx;
  pbeq->ncly = ncly;
  pbeq->nclz = nclz;
  pbeq->swap = swap;

  pbeq->vol = new molfile_volumetric_t[1];
  strcpy(pbeq->vol[0].dataname, "CHARMM PBEQ Potential Map");

  pbeq->vol[0].origin[0] = -0.5*((nclx-1) * dcel) + xbcen;
  pbeq->vol[0].origin[1] = -0.5*((ncly-1) * dcel) + ybcen;
  pbeq->vol[0].origin[2] = -0.5*((nclz-1) * dcel) + zbcen;

  // print origin info, for debuggin of early versions
  printf("pbeqplugin) box LL origin: %g %g %g\n", 
         pbeq->vol[0].origin[0],
         pbeq->vol[0].origin[1],
         pbeq->vol[0].origin[2]);

  pbeq->vol[0].xaxis[0] = (nclx-1) * dcel;
  pbeq->vol[0].xaxis[1] = 0;
  pbeq->vol[0].xaxis[2] = 0;

  pbeq->vol[0].yaxis[0] = 0;
  pbeq->vol[0].yaxis[1] = (ncly-1) * dcel;
  pbeq->vol[0].yaxis[2] = 0;

  pbeq->vol[0].zaxis[0] = 0;
  pbeq->vol[0].zaxis[1] = 0;
  pbeq->vol[0].zaxis[2] = (nclz-1) * dcel;

  pbeq->vol[0].xsize = nclx;
  pbeq->vol[0].ysize = ncly;
  pbeq->vol[0].zsize = nclz;

  pbeq->vol[0].has_color = 0;

  return pbeq;
}

static int read_pbeq_metadata(void *v, int *nsets, 
  molfile_volumetric_t **metadata) {
  pbeq_t *pbeq = (pbeq_t *)v;
  *nsets = pbeq->nsets; 
  *metadata = pbeq->vol;  

  return MOLFILE_SUCCESS;
}

static int read_pbeq_data(void *v, int set, float *datablock,
                         float *colorblock) {
  pbeq_t *pbeq = (pbeq_t *)v;
  int ndata = pbeq->ndata;
  int nclx = pbeq->nclx;
  int ncly = pbeq->ncly;
  int nclz = pbeq->nclz;
  FILE *fd = pbeq->fd;
  int trash;

  // skip first Fortran length record for 
  // WRITE(UNIT)(PHI(I),I=1,NCLX*NCLY*NCLZ)
  if (fread(&trash, 4, 1, fd) != 1)
    return MOLFILE_ERROR;

  /* Read the densities. Order for file is z fast, y medium, x slow */
  int x, y, z;
  int count=0;
  for (x=0; x<nclx; x++) {
    for (y=0; y<ncly; y++) {
      for (z=0; z<nclz; z++) {
        int addr = z*nclx*ncly + y*nclx + x;
        if (fread(datablock + addr, 4, 1, fd) != 1) {
          printf("pbeqplugin) Error reading potential map cell: %d,%d,%d\n", x, y, z);
          printf("pbeqplugin) offset: %d\n", ftell(fd));
          return MOLFILE_ERROR;
        }
        count++;
      }
    }
  }


#if 0
  // skip last Fortran length record for 
  // WRITE(UNIT)(PHI(I),I=1,NCLX*NCLY*NCLZ)
  if (fread(&trash, 4, 1, fd) != 1)
    return NULL;

  // print debugging info for early versions
  printf("pbeqplugin) read %d phi values from block\n", count);
  for (x=0; x<1000000; x++) {
    int trash;
    if (fread(&trash, 4, 1, fd) != 1) {
      printf("pbeqplugin) read %d extra phi values past the end of the block\n", x);
      break;
    }
  }   
#endif

  if (pbeq->swap) {
    swap4_aligned(datablock, ndata);
  }

  return MOLFILE_SUCCESS;
}

static void close_pbeq_read(void *v) {
  pbeq_t *pbeq = (pbeq_t *)v;

  fclose(pbeq->fd);
  if (pbeq->vol != NULL)
    delete [] pbeq->vol; 
  delete pbeq;
}

/*
 * Initialization stuff here
 */
static molfile_plugin_t plugin;

VMDPLUGIN_EXTERN int VMDPLUGIN_init(void) { 
  memset(&plugin, 0, sizeof(molfile_plugin_t));
  plugin.abiversion = vmdplugin_ABIVERSION;
  plugin.type = MOLFILE_PLUGIN_TYPE;
  plugin.name = "pbeq";
  plugin.prettyname = "CHARMM PBEQ Binary Potential Map";
  plugin.author = "John Stone";
  plugin.majorv = 0;
  plugin.minorv = 3;
  plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
  plugin.filename_extension = "pbeq, phi80";
  plugin.open_file_read = open_pbeq_read;
  plugin.read_volumetric_metadata = read_pbeq_metadata;
  plugin.read_volumetric_data = read_pbeq_data;
  plugin.close_file_read = close_pbeq_read;
  return VMDPLUGIN_SUCCESS; 
}

VMDPLUGIN_EXTERN int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_EXTERN int VMDPLUGIN_fini(void) { return VMDPLUGIN_SUCCESS; }

