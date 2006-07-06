/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_grdplugin
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
 *      $RCSfile: grdplugin.C,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.13 $       $Date: 2006/02/23 19:36:44 $
 *
 ***************************************************************************/

/* 
 * "unformatted" binary potential map, as used by Grasp and DelPhi
 *
 * Format (fortran): 
 * character*20 uplbl
 * character*10 nxtlbl,character*60 toplbl
 * real*4 phi(n,n,n)
 * character*16 botlbl
 * real*4 scale,oldmid(3)
 *
 * Where n is the length in grid units of each edge of the grid.
 *
 * More information can be found at:
 * <http://honiglab.cpmc.columbia.edu/grasp/grasp_contents.html#A.2>
 * <http://trantor.bioc.columbia.edu/delphi/doc/file_format.html>
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
  int swap;
  molfile_volumetric_t *vol;
} grd_t;


static void *open_grd_read(const char *filepath, const char *filetype,
    int *natoms) {
  FILE *fd;
  grd_t *grd;
  char uplbl[21], nxtlbl[11], toplbl[61];
  int swap, recordSize, gridSize, iGrid;
  float scale, midX, midY, midZ;
  
  fd = fopen(filepath, "rb");
  if (!fd) {
    fprintf(stderr, "Error opening file.\n");
    return NULL;
  }

  /* Check byte order. The first four bytes of the file always make up the
   * integer 20.
   */
  if (fread(&recordSize, 4, 1, fd) != 1) {
    fprintf(stderr, "Error reading file header: uplbl.\n");
    return NULL;
  }
  if (recordSize == 20) {
    swap = 0;
  }
  else {
    swap4_aligned(&recordSize, 1);
    if (recordSize == 20) {
      swap = 1;
    }
    else {
      fprintf(stderr, "Improperly formatted file header: uplbl.\n");
      return NULL;
    }
  }
  
  /* Check for a valid phimap 
   * XXX - Some programs write gibberish for this record, don't worry about
   * its contents
   * character*20 uplbl
   */
  if ( (fread(uplbl, 1, 20, fd) != 20) ||
       (fread(&recordSize, 4, 1, fd) != 1) ) {
    fprintf(stderr, "Error: uplbl does not match.\n");
    return NULL;
  }

  /* Read in the next record:
   * character*10 nxtlbl, character*60 toplbl 
   * The labels themselves are currently ignored, but they may be useful in
   * the future.
   */ 
  if (fread(&recordSize, 4, 1, fd) != 1) {
    fprintf(stderr, "Error reading file header: nxtlbl.\n");
    return NULL;
  }
  if (swap) {
    swap4_aligned(&recordSize, 1);
  }
  if (recordSize != 70) {
    fprintf(stderr, "Improperly formatted file header: nxtlbl.\n");
    return NULL;
  }
  if ( (fread(nxtlbl, 1, 10, fd) != 10) ||
       (fread(toplbl, 1, 60, fd) != 60) ||
       (fread(&recordSize, 4, 1, fd) != 1) ) {
    fprintf(stderr, "Error reading nxtlbl.\n");
    return NULL;
  }
  
  /* Find the number of data points in the file
   * The next integer gives the number of bytes used to store the data.
   */
  if (fread(&recordSize, 4, 1, fd) != 1) {
    fprintf(stderr, "Error reading file header: grid.\n");
    return NULL;
  }
  if (swap) {
    swap4_aligned(&recordSize, 1);
  }
  iGrid = recordSize / 4;

  /* Find the length in grid units of the edge of the cube, make sure it's
   * an integer 
   */
  gridSize = (int) (pow((double) iGrid, (double) 1.0/3.0) + 0.5);
  if ((gridSize*gridSize*gridSize) != iGrid) {
    fprintf(stderr, "Error: non-cube grid.\n");
    return NULL;
  }

  /* Read the scale and midpoint coordinates from the end of the file.
   */
  if ( (fseek(fd, -20, SEEK_END) != 0) ||
       (fread(&scale, sizeof(float), 1, fd) != 1) ||
       (fread(&midX, sizeof(float), 1, fd) != 1) ||
       (fread(&midY, sizeof(float), 1, fd) != 1) ||
       (fread(&midZ, sizeof(float), 1, fd) != 1) ) {
    fprintf(stderr, "Error reading scale and midpoint.\n");
    return NULL;
  }
  if (swap) {
    swap4_aligned(&scale, 1);
    swap4_aligned(&midX, 1);
    swap4_aligned(&midY, 1);
    swap4_aligned(&midZ, 1);
  }

  /* Allocate and initialize the grd structure */
  grd = new grd_t;
  grd->fd = fd;
  grd->vol = NULL;
  *natoms = MOLFILE_NUMATOMS_NONE;
  grd->nsets = 1; /* this file contains only one data set */
  grd->ndata = iGrid;
  grd->swap = swap;

  grd->vol = new molfile_volumetric_t[1];
  strcpy(grd->vol[0].dataname, "PHIMAP Electron Density Map");

  /* <midX, midY, midZ> is the middle point of the grid. */
  grd->vol[0].origin[0] = -0.5*(gridSize+1.0) / scale + midX;
  grd->vol[0].origin[1] = -0.5*(gridSize+1.0) / scale + midY;
  grd->vol[0].origin[2] = -0.5*(gridSize+1.0) / scale + midZ;

  grd->vol[0].xaxis[0] = gridSize / scale;
  grd->vol[0].xaxis[1] = 0;
  grd->vol[0].xaxis[2] = 0;

  grd->vol[0].yaxis[0] = 0;
  grd->vol[0].yaxis[1] = gridSize / scale;
  grd->vol[0].yaxis[2] = 0;

  grd->vol[0].zaxis[0] = 0;
  grd->vol[0].zaxis[1] = 0;
  grd->vol[0].zaxis[2] = gridSize / scale;

  grd->vol[0].xsize = gridSize;
  grd->vol[0].ysize = gridSize;
  grd->vol[0].zsize = gridSize;

  grd->vol[0].has_color = 0;

  return grd;
}

static int read_grd_metadata(void *v, int *nsets, 
  molfile_volumetric_t **metadata) {
  grd_t *grd = (grd_t *)v;
  *nsets = grd->nsets; 
  *metadata = grd->vol;  

  return MOLFILE_SUCCESS;
}

static int read_grd_data(void *v, int set, float *datablock,
                         float *colorblock) {
  grd_t *grd = (grd_t *)v;
  int ndata = grd->ndata;
  FILE *fd = grd->fd;

  /* Skip the header */
  fseek(fd, 110, SEEK_SET);

  /* Read the densities. Order for file is x fast, y medium, z slow */
  if (fread(datablock, sizeof(float), ndata, fd) != (unsigned int) ndata) {
    fprintf(stderr, "Error reading grid data.\n");
    return MOLFILE_ERROR;
  }

  if (grd->swap) {
    swap4_aligned(datablock, ndata);
  }

  return MOLFILE_SUCCESS;
}

static void close_grd_read(void *v) {
  grd_t *grd = (grd_t *)v;

  fclose(grd->fd);
  if (grd->vol != NULL)
    delete [] grd->vol; 
  delete grd;
}

/*
 * Initialization stuff here
 */
static molfile_plugin_t plugin = {
  vmdplugin_ABIVERSION,                /* ABI version */
  MOLFILE_PLUGIN_TYPE, 	               /* plugin type */
  "grd",                               /* file format description */
  "GRASP,Delphi Binary Potential Map", /* file format description */
  "Eamon Caddigan",                    /* author(s) */
  0,                                   /* major version */
  5,                                   /* minor version */
  VMDPLUGIN_THREADSAFE,                /* is reentrant */
  "phi,grd"                            /* filename extension */
};

VMDPLUGIN_EXTERN int VMDPLUGIN_init(void) { return VMDPLUGIN_SUCCESS; }
VMDPLUGIN_EXTERN int VMDPLUGIN_fini(void) { return VMDPLUGIN_SUCCESS; }
VMDPLUGIN_EXTERN int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  plugin.open_file_read = open_grd_read;
  plugin.read_volumetric_metadata = read_grd_metadata;
  plugin.read_volumetric_data = read_grd_data;
  plugin.close_file_read = close_grd_read;
  (*cb)(v, (vmdplugin_t *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

