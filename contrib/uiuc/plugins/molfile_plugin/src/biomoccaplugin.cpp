/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_biomoccaplugin
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
 *      $RCSfile: biomoccaplugin.C,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.5 $       $Date: 2009/04/29 15:45:28 $
 *
 ***************************************************************************/

/* 
 * Biomocca volumetric map file reader
 *   Biomocca is written by the CEG at UIUC:
 *     http://www.ceg.uiuc.edu/
 *
 * File format (simple ASCII text): 
 * Xcenter Ycenter Zcenter (in Angstroms)
 * Nx(number of cells on the x axis)  Ny  Nz
 * d (cell spacing, in Angstroms)
 * Voxel values (-1, 0, 1, ...) stored in Z/Y/X fortran style order
 *
 * Meaning of voxel values: 
 * -1 for lipid
 *  0 for channel or solvent baths
 *  1 stands for the protein
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
} biomocca_t;


static void *open_biomocca_read(const char *filepath, const char *filetype,
    int *natoms) {
  FILE *fd;
  biomocca_t *biomocca;
  float scale;
  int xsize, ysize, zsize;
  float orig[3];
  
  fd = fopen(filepath, "r");
  if (!fd) {
    printf("biomoccaplugin) Error opening file.\n");
    return NULL;
  }

  if (fscanf(fd, "%f %f %f", orig, orig+1, orig+2) != 3) {
    printf("biomoccaplugin) Error reading grid origin.\n");
    return NULL;
  }

  /* get the number of grid points */
  if (fscanf(fd, "%d %d %d", &xsize, &ysize, &zsize) != 3) {
    printf("biomoccaplugin) Error reading grid dimensions.\n");
    return NULL;
  }

  /* get the voxel scale */
  if (fscanf(fd, "%f", &scale) != 1) {;
    printf("biomoccaplugin) Error reading voxel scale.\n");
    return NULL;
  }

  /* allocate and initialize the biomocca structure */
  biomocca = new biomocca_t;
  biomocca->fd = fd;
  biomocca->vol = NULL;
  *natoms = MOLFILE_NUMATOMS_NONE;
  biomocca->nsets = 1; /* this file contains only one data set */

  biomocca->vol = new molfile_volumetric_t[1];
  strcpy(biomocca->vol[0].dataname, "BioMocca map");

  /* Set the unit cell origin and basis vectors */
  for (int i=0; i<3; i++) {
    biomocca->vol[0].origin[i] = orig[i];
    biomocca->vol[0].xaxis[i] = 0.0;
    biomocca->vol[0].yaxis[i] = 0.0;
    biomocca->vol[0].zaxis[i] = 0.0;
  }

  biomocca->vol[0].xaxis[0] = scale * (xsize-1);
  biomocca->vol[0].yaxis[1] = scale * (ysize-1);
  biomocca->vol[0].zaxis[2] = scale * (zsize-1);

  biomocca->vol[0].origin[0] -= 0.5 * biomocca->vol[0].xaxis[0];
  biomocca->vol[0].origin[1] -= 0.5 * biomocca->vol[0].yaxis[1];
  biomocca->vol[0].origin[2] -= 0.5 * biomocca->vol[0].zaxis[2];

  biomocca->vol[0].xsize = xsize;
  biomocca->vol[0].ysize = ysize;
  biomocca->vol[0].zsize = zsize;

  biomocca->vol[0].has_color = 0; /* BioMocca maps contain no color info */

  return biomocca;
}

static int read_biomocca_metadata(void *v, int *nsets, 
  molfile_volumetric_t **metadata) {
  biomocca_t *biomocca = (biomocca_t *)v;
  *nsets = biomocca->nsets; 
  *metadata = biomocca->vol;  

  return MOLFILE_SUCCESS;
}

static int read_biomocca_data(void *v, int set, float *datablock,
                         float *colorblock) {
  biomocca_t *biomocca = (biomocca_t *)v;
  FILE *fd = biomocca->fd;
  int x, y, z, xsize, ysize, zsize, xysize;

  xsize = biomocca->vol[0].xsize;
  ysize = biomocca->vol[0].ysize;
  zsize = biomocca->vol[0].zsize;
  xysize = xsize * ysize;

  for (x=0; x<xsize; x++) {
    for (y=0; y<ysize; y++) {
      for (z=0; z<zsize; z++) {
        if (fscanf(fd, "%f", datablock + z*xysize + y*xsize + x) != 1) {
          printf("biomoccaplugin) Failed reading biomocca map data\n");
          return MOLFILE_ERROR;
        }
  
      }
    }
  }

  return MOLFILE_SUCCESS;
}

static void close_biomocca_read(void *v) {
  biomocca_t *biomocca = (biomocca_t *)v;
  
  fclose(biomocca->fd);
  if (biomocca->vol != NULL)
    delete [] biomocca->vol; 
  delete biomocca;
}

/*
 * Initialization stuff here
 */
static molfile_plugin_t plugin;

VMDPLUGIN_EXTERN int VMDPLUGIN_init(void) { 
  memset(&plugin, 0, sizeof(molfile_plugin_t)); 
  plugin.abiversion = vmdplugin_ABIVERSION;
  plugin.type = MOLFILE_PLUGIN_TYPE;
  plugin.name = "biomocca";
  plugin.prettyname = "Biomocca Volumetric Map";
  plugin.author = "John Stone";
  plugin.majorv = 0;
  plugin.minorv = 2;
  plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
  plugin.filename_extension = "bmcg";
  plugin.open_file_read = open_biomocca_read;
  plugin.read_volumetric_metadata = read_biomocca_metadata;
  plugin.read_volumetric_data = read_biomocca_data;
  plugin.close_file_read = close_biomocca_read;
 
  return VMDPLUGIN_SUCCESS; 
}

VMDPLUGIN_EXTERN int VMDPLUGIN_fini(void) { return VMDPLUGIN_SUCCESS; }

VMDPLUGIN_EXTERN int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

