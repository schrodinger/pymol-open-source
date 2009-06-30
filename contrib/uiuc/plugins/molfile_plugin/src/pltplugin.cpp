/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_pltplugin
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
 *      $RCSfile: pltplugin.C,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.13 $       $Date: 2009/04/29 15:45:32 $
 *
 ***************************************************************************/

/* 
 * plt format electron density maps from gOpenMol.
 *
 * More info for format can be found at 
 * <http://www.csc.fi/gopenmol/developers/plt_format.phtml>
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
  int swap;
  molfile_volumetric_t *vol;
} plt_t;


static void *open_plt_read(const char *filepath, const char *filetype,
    int *natoms) {
  FILE *fd;
  plt_t *plt;
  int swap=0;
  // File header data:
  int intHeader[5];
  float floatHeader[6];
  
  fd = fopen(filepath, "rb");
  if (!fd) {
    fprintf(stderr, "pltplugin) Error opening file.\n");
    return NULL;
  }

  // Integer header info: rank (always 3), surface type, z length, y length, 
  // x length.
  fread(intHeader, sizeof(int), 5, fd);
  if (intHeader[0] != 3) {
    // check if the bytes need to be swapped
    swap4_aligned(intHeader, 5);
    if (intHeader[0] == 3) 
      swap = 1;
    else {
      fprintf(stderr, "pltplugin) Incorrect header.\n");
      return NULL;
    }
  }

  // Float header info: z min, z max, y min, y max, xmin, x max.
  fread(floatHeader, sizeof(float), 6, fd);
  if (swap)
    swap4_aligned(floatHeader, 6);

  // Allocate and initialize the plt structure
  plt = new plt_t;
  plt->fd = fd;
  plt->vol = NULL;
  *natoms = MOLFILE_NUMATOMS_NONE;
  plt->nsets = 1; // this file contains only one data set
  plt->swap = swap;

  plt->vol = new molfile_volumetric_t[1];
  strcpy(plt->vol[0].dataname, "PLT Electron Density Map");

  // Best guesses for unit cell information, as none is included in the plt
  // file format.
  plt->vol[0].origin[0] = floatHeader[4];
  plt->vol[0].origin[1] = floatHeader[2];
  plt->vol[0].origin[2] = floatHeader[0];

  plt->vol[0].xaxis[0] = floatHeader[5] - floatHeader[4];
  plt->vol[0].xaxis[1] = 0;
  plt->vol[0].xaxis[2] = 0;

  plt->vol[0].yaxis[0] = 0;
  plt->vol[0].yaxis[1] = floatHeader[3] - floatHeader[2];
  plt->vol[0].yaxis[2] = 0;
  
  plt->vol[0].zaxis[0] = 0;
  plt->vol[0].zaxis[1] = 0;
  plt->vol[0].zaxis[2] = floatHeader[1] - floatHeader[0];

  plt->vol[0].xsize = intHeader[4];
  plt->vol[0].ysize = intHeader[3];
  plt->vol[0].zsize = intHeader[2];

  plt->vol[0].has_color = 0;

  return plt;
}

static int read_plt_metadata(void *v, int *nsets, 
  molfile_volumetric_t **metadata) {
  plt_t *plt = (plt_t *)v;
  *nsets = plt->nsets; 
  *metadata = plt->vol;  

  return MOLFILE_SUCCESS;
}

static int read_plt_data(void *v, int set, float *datablock,
                         float *colorblock) {
  plt_t *plt = (plt_t *)v;
  int swap, ndata;
  FILE *fd = plt->fd;

  swap = plt->swap;
  ndata = plt->vol[0].xsize * plt->vol[0].ysize * plt->vol[0].zsize;

  // Read the densities. Order for file is x fast, y medium, z slow
  if ( fread(datablock, sizeof(float), ndata, fd) != ndata ) {
    fprintf(stderr, "pltplugin) Error reading data, not enough values read.\n");
    return MOLFILE_ERROR;
  }

  if (swap) 
    swap4_aligned(datablock, ndata);

  return MOLFILE_SUCCESS;
}

static void close_plt_read(void *v) {
  plt_t *plt = (plt_t *)v;

  fclose(plt->fd);
  if (plt->vol != NULL)
    delete [] plt->vol; 
  delete plt;
}

/*
 * Initialization stuff here
 */
static molfile_plugin_t plugin;

VMDPLUGIN_EXTERN int VMDPLUGIN_init(void) {
  memset(&plugin, 0, sizeof(molfile_plugin_t));
  plugin.abiversion = vmdplugin_ABIVERSION;
  plugin.type = MOLFILE_PLUGIN_TYPE;
  plugin.name = "plt";
  plugin.prettyname = "gOpenmol plt";
  plugin.author = "Eamon Caddigan";
  plugin.majorv = 0;
  plugin.minorv = 4;
  plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
  plugin.filename_extension = "plt";
  plugin.open_file_read = open_plt_read;
  plugin.read_volumetric_metadata = read_plt_metadata;
  plugin.read_volumetric_data = read_plt_data;
  plugin.close_file_read = close_plt_read;
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_EXTERN int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_EXTERN int VMDPLUGIN_fini(void) { return VMDPLUGIN_SUCCESS; }

