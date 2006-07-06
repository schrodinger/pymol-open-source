/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_situsplugin
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
 *      $RCSfile: situsplugin.C,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.6 $       $Date: 2006/03/30 02:34:33 $
 *
 ***************************************************************************/

/* 
 * Situs EM map file reader
 *   map is a cubic lattice
 *   all data is stored in plain ASCII for convenience
 *
 * Format of the file is:
 *   voxel size in Angstroms
 *   coordinates of first voxel
 *   integer X/Y/Z counts
 *   voxels follow in X fastest, Y next fastest, Z slowest
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#if defined(_AIX)
#include <strings.h>
#endif

#include "molfile_plugin.h"

typedef struct {
  FILE *fd;
  int nsets;
  molfile_volumetric_t *vol;
} situs_t;


static void *open_situs_read(const char *filepath, const char *filetype,
    int *natoms) {
  FILE *fd;
  situs_t *situs;
  float scale;
  int xsize, ysize, zsize;
  float orig[3]; 
  
  fd = fopen(filepath, "r");
  if (!fd) {
    printf("situsplugin) Error opening file.\n");
    return NULL;
  }

  /* get the scale */
  if (fscanf(fd, "%f", &scale) != 1) {;
    printf("situsplugin) Error reading voxel scale.\n");
    return NULL;
  }

  if (fscanf(fd, "%f %f %f", orig, orig+1, orig+2) != 3) {
    printf("situsplugin) Error reading grid origin.\n");
    return NULL;
  }

  /* get the number of grid points */
  if (fscanf(fd, "%d %d %d", &xsize, &ysize, &zsize) != 3) {
    printf("situsplugin) Error reading grid dimensions.\n");
    return NULL;
  }

  /* allocate and initialize the situs structure */
  situs = new situs_t;
  situs->fd = fd;
  situs->vol = NULL;
  *natoms = MOLFILE_NUMATOMS_NONE;
  situs->nsets = 1; /* this file contains only one data set */

  situs->vol = new molfile_volumetric_t[1];
  strcpy(situs->vol[0].dataname, "Situs map");

  /* Set the unit cell origin and basis vectors */
  for (int i=0; i<3; i++) {
    situs->vol[0].origin[i] = orig[i];
    situs->vol[0].xaxis[i] = 0.0;
    situs->vol[0].yaxis[i] = 0.0;
    situs->vol[0].zaxis[i] = 0.0;
  }
  situs->vol[0].xaxis[0] = scale * (xsize-1);
  situs->vol[0].yaxis[1] = scale * (ysize-1);
  situs->vol[0].zaxis[2] = scale * (zsize-1);

  situs->vol[0].xsize = xsize;
  situs->vol[0].ysize = ysize;
  situs->vol[0].zsize = zsize;

  situs->vol[0].has_color = 0; /* Situs maps contain no color information. */

  return situs;
}

static int read_situs_metadata(void *v, int *nsets, 
  molfile_volumetric_t **metadata) {
  situs_t *situs = (situs_t *)v;
  *nsets = situs->nsets; 
  *metadata = situs->vol;  

  return MOLFILE_SUCCESS;
}

static int read_situs_data(void *v, int set, float *datablock,
                         float *colorblock) {
  situs_t *situs = (situs_t *)v;
  FILE *fd = situs->fd;
  int xsize, ysize, zsize, xysize, count, total;

  xsize = situs->vol[0].xsize;
  ysize = situs->vol[0].ysize;
  zsize = situs->vol[0].zsize;
  xysize = xsize * ysize;
  total = xysize * zsize;

  /* Read the values from the file */
  for (count=0; count<total; count++) {
    if (fscanf(fd, "%f", datablock + count) != 1) {
      printf("situsplugin) Failed reading situs map data\n");
      return MOLFILE_ERROR;
    }
  }   

  return MOLFILE_SUCCESS;
}

static void close_situs_read(void *v) {
  situs_t *situs = (situs_t *)v;
  
  fclose(situs->fd);
  if (situs->vol != NULL)
    delete [] situs->vol; 
  delete situs;
}

/*
 * Initialization stuff here
 */
static molfile_plugin_t plugin = {
  vmdplugin_ABIVERSION,   /* ABI version */
  MOLFILE_PLUGIN_TYPE, 	  /* plugin type */
  "situs",                /* file format description */
  "Situs Density Map",    /* file format description */
  "John Stone",           /* author(s) */
  0,                      /* major version */
  1,                      /* minor version */
  VMDPLUGIN_THREADSAFE,   /* is reentrant */
  "situs"                 /* filename extension */
};

VMDPLUGIN_EXTERN int VMDPLUGIN_init(void) { return VMDPLUGIN_SUCCESS; }
VMDPLUGIN_EXTERN int VMDPLUGIN_fini(void) { return VMDPLUGIN_SUCCESS; }
VMDPLUGIN_EXTERN int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  plugin.open_file_read = open_situs_read;
  plugin.read_volumetric_metadata = read_situs_metadata;
  plugin.read_volumetric_data = read_situs_data;
  plugin.close_file_read = close_situs_read;
  (*cb)(v, (vmdplugin_t *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

