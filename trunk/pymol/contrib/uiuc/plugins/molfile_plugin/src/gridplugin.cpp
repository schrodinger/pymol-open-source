/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_gridplugin
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
 *      $RCSfile: gridplugin.C,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.11 $       $Date: 2006/02/23 19:36:44 $
 *
 ***************************************************************************/

/* 
 * Binary potential map format used by Molecular Discovery GRID, 
 * UHBD, and other packages.
 *
 * Files begin with a 160-byte formatted fortran header.
 * For each plane in the grid, there is a 12-byte formatted fortran record
 * giving the grid coordinates, followed by a nx*ny*4 byte formatted fortran
 * record with the grid data.
 *
 * The header has the following format (addresses in bytes):
 * header[0..71]:     char * 72, grid title
 * header[72..79]:    float * 2, grid clearance and cutoff energy (unused by
 *                    VMD)
 * header[80..99]:    unknown
 * header[100..111]:  int * 3, number of planes in each direction
 * header[112..127]:  float * 4, grid spacing in angstroms, followed by grid
 *                    origin
 * header[128..151]:  float * 6, VDW radius, neff, alph, q, emin, rmin
 *                    (unused by VMD)
 * header[152..159]:  int * 2, jd, ja (unused by VMD)
 *
 * XXX - Not sure if slicing order is the same in every file. Also, other
 * values from plane metadata seem to have no use -- I'm probably missing
 * something. 
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#if defined(_AIX)
#include <strings.h>
#endif

#include "molfile_plugin.h"
#include "endianswap.h"
#include "fortread.h"

typedef struct {
  FILE *fd;
  int swap;
  molfile_volumetric_t *vol;
} grid_t;


static void *open_grid_read(const char *filepath, const char *filetype,
    int *natoms) {
  FILE *fd;
  grid_t *grid;
  float header[64], ra, rx, ry, rz;
  int dataBegin, swap, blocksize, nnx, nny, nnz;
  
  fd = fopen(filepath, "rb");
  if (!fd) {
    fprintf(stderr, "Error opening file.\n");
    return NULL;
  }

  // Use the first four-byte integer in the file to determine the file's
  // byte-order
  fread(&dataBegin, sizeof(int), 1, fd);
  if ( (dataBegin > 255) || (dataBegin < 0) ) {
    // check if the bytes need to be swapped
    swap4_aligned(&dataBegin, 1);
    if (dataBegin <= 255) {
      swap = 1;
    } else {
      fprintf(stderr, "Cannot read file: header block is too large.\n");
      return NULL;
    }
  }
  else {
    swap = 0;
  }

  // Read the header
  rewind(fd);
  blocksize = fortread_4(header, 64, swap, fd);

  if (blocksize != 40) {
    fprintf(stderr, "Incorrect header size.\n");
    return NULL;
  }

  // number of planes in each dimension
  nnx = ((int *)header)[25];
  nny = ((int *)header)[26];
  nnz = ((int *)header)[27];

  // grid spacing in angstroms
  ra = header[28];

  // reference point for grid position
  rx = header[29];
  ry = header[30];
  rz = header[31];

  // Allocate and initialize the structure 
  grid = new grid_t;
  grid->fd = fd;
  grid->vol = NULL;
  *natoms = MOLFILE_NUMATOMS_NONE;
  grid->swap = swap;

  grid->vol = new molfile_volumetric_t[1];
  strcpy(grid->vol[0].dataname, "GRID Electron Density Map");

  // XXX - origin seems to be shifted by one grid point in each direction
  // from the reference point. 
  grid->vol[0].origin[0] = rx + ra;
  grid->vol[0].origin[1] = ry + ra;
  grid->vol[0].origin[2] = rz + ra;

  grid->vol[0].xaxis[0] = nnx * ra;
  grid->vol[0].xaxis[1] = 0;
  grid->vol[0].xaxis[2] = 0;

  grid->vol[0].yaxis[0] = 0;
  grid->vol[0].yaxis[1] = nny * ra;
  grid->vol[0].yaxis[2] = 0;

  grid->vol[0].zaxis[0] = 0;
  grid->vol[0].zaxis[1] = 0;
  grid->vol[0].zaxis[2] = nnz * ra;

  grid->vol[0].xsize = nnx;
  grid->vol[0].ysize = nny;
  grid->vol[0].zsize = nnz;

  grid->vol[0].has_color = 0;   // This file has no color

  return grid;
}

static int read_grid_metadata(void *v, int *nsets, 
  molfile_volumetric_t **metadata) {
  grid_t *grid = (grid_t *)v;
  *nsets = 1;                   // This file contains only one data set.
  *metadata = grid->vol;  

  return MOLFILE_SUCCESS;
}

static int read_grid_data(void *v, int set, float *datablock,
                         float *colorblock) {
  grid_t *grid = (grid_t *)v;
  int planeHeader[3], planeSize, i, z;
  float *planeData;

  planeSize = grid->vol[0].xsize * grid->vol[0].ysize;

  planeData = new float[planeSize];

  for (i = 0; i < grid->vol[0].zsize; i++) {
    // read the plane metadata
    if (fortread_4(planeHeader, 3, grid->swap, grid->fd) != 3) {
      fprintf(stderr, "Error reading plane metadata.\n");
      delete [] planeData;
      return MOLFILE_ERROR;
    }
    z = planeHeader[0] - 1;

    // read the plane data
    if (fortread_4(planeData, planeSize, grid->swap, grid->fd) != planeSize) {
      fprintf(stderr, "Error reading plane data.\n");
      delete [] planeData;
      return MOLFILE_ERROR;
    }

    // copy the plane data to the datablock
    memcpy(datablock + z*planeSize, planeData, planeSize * sizeof(float));
  }

  delete [] planeData;
  return MOLFILE_SUCCESS;
}

static void close_grid_read(void *v) {
  grid_t *grid = (grid_t *)v;

  fclose(grid->fd);
  if (grid->vol != NULL)
    delete [] grid->vol; 
  delete grid;
}

/*
 * Initialization stuff here
 */
static molfile_plugin_t plugin = {
  vmdplugin_ABIVERSION,              // ABI version
  MOLFILE_PLUGIN_TYPE,               // plugin type
  "grid",                            // file format description
  "GRID,UHBD Binary Potential Map",  // file format description
  "Eamon Caddigan",                  // author(s)
  0,                                 // major version
  2,                                 // minor version
  VMDPLUGIN_THREADSAFE,              // is reentrant
  "grid"                             // filename extension
};

VMDPLUGIN_EXTERN int VMDPLUGIN_init(void) { return VMDPLUGIN_SUCCESS; }
VMDPLUGIN_EXTERN int VMDPLUGIN_fini(void) { return VMDPLUGIN_SUCCESS; }
VMDPLUGIN_EXTERN int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  plugin.open_file_read = open_grid_read;
  plugin.read_volumetric_metadata = read_grid_metadata;
  plugin.read_volumetric_data = read_grid_data;
  plugin.close_file_read = close_grid_read;
  (*cb)(v, (vmdplugin_t *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

