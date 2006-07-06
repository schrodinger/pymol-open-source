/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_fs4plugin
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
 *      $RCSfile: fs4plugin.C,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.18 $       $Date: 2006/02/23 19:36:44 $
 *
 ***************************************************************************/

/* 
 * fsfour format density maps
 *
 * More info for this format can be found at 
 * <http://www.csb.yale.edu/userguides/datamanip/phases/FSFOUR.html>
 *
 * Old versions of the cns2fsfour and ccp2fsfour utilities produced a
 * slightly broken variant of this format, omitting the data that isn't used
 * by XTalView. This plugin currently reads these files, but the cell
 * dimensions will be incorrect.
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
#include "fortread.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct {
  FILE *fd;
  int nsets;
  int swap;
  int f, m, s, x, y, z; // indecies mapping fast-medium-slow to x-y-z
  float scale;
  molfile_volumetric_t *vol;
} fs4_t;


static void *open_fs4_read(const char *filepath, const char *filetype,
    int *natoms) {
  fs4_t *fs4;
  FILE *fd;
  float header[32], scale, fmsCellSize[3], alpha, beta, gamma, z1, z2, z3;
  int dataBegin, blocksize, geom[16], fmsGridSize[3], norn, swap=0;

  fd = fopen(filepath, "rb");
  if (!fd) {
    fprintf(stderr, "Error opening file.\n");
    return NULL;
  }

  // Use the first four-byte integer in the file to determine the file's
  // byte-order
  fread(&dataBegin, sizeof(int), 1, fd);
  if (dataBegin > 255) {
    // check if the bytes need to be swapped
    swap4_aligned(&dataBegin, 1);
    if (dataBegin <= 255) {
      swap = 1;
    } else {
      fprintf(stderr, "Cannot read file: header block is too large.\n");
      return NULL;
    }
  }

  // Read the header
  rewind(fd);
  blocksize = fortread_4(header, 32, swap, fd);

  // Handle files produced by old versions of (cns|ccp)2fsfour
  if (blocksize == 28) {
    printf("Recognized %s cns2fsfour map.\n", 
        swap ? "opposite-endian" : "same-endian");

    // Read the geometry block
    blocksize = fortread_4(geom, 16, swap, fd);
    if (blocksize != 7) {
      fprintf(stderr, "Incorrect size for geometry block.\n");
      return NULL;
    }

    fmsGridSize[0] = geom[0];
    fmsGridSize[1] = geom[1];
    fmsGridSize[2] = geom[2];
    norn = geom[4];

    // Warn about assumptions
    printf("Warning: file does not contain unit cell lengths or angles.\n");

    scale = 50.0;
    fmsCellSize[0] = 1.0;
    fmsCellSize[1] = 1.0;
    fmsCellSize[2] = 1.0;
    alpha = 90.0;
    beta = 90.0;
    gamma = 90.0;
  }

  // Handle standard fsfour files
  else if (blocksize == 31) {
    printf("Recognize standard fsfour map.\n");

    fmsCellSize[0] = header[21];
    fmsCellSize[1] = header[22];
    fmsCellSize[2] = header[23];
    alpha = header[24];
    beta = header[25];
    gamma = header[26];
    
    // Skip the symmetry block if one present
    blocksize = fortread_4(geom, 16, swap, fd);
    if (blocksize == 9) {
      printf("Skipping symmetry block.\n");
      blocksize = fortread_4(geom, 16, swap, fd);
    }

    // Read the geometry block
    if (blocksize != 13) {
      fprintf(stderr, "Incorrect size for geometry block.\n");
      return NULL;
    }

    fmsGridSize[0] = geom[0];
    fmsGridSize[1] = geom[1];
    fmsGridSize[2] = geom[2];

    scale = *((float *) geom + 3);
    if (scale == 0) {
      scale = 50;
    }

    norn = geom[4];
    if ((norn < 0) || (norn > 2)) {
      fprintf(stderr, "norn out of range.\n");
      return NULL;
    }
  }

  // Unrecognized format
  else {
    fprintf(stderr, "Unrecognized map format.\n");
    return NULL;
  }

  // Convert degrees to radians
  alpha *= M_PI / 180.0;
  beta  *= M_PI / 180.0;
  gamma *= M_PI / 180.0;

  // Warn about assumptions
  printf("Warning: file does not contain molecule center.\nCentering at <0, 0, 0>\n");

  // Allocate and initialize the fs4 structure
  fs4 = new fs4_t;
  fs4->fd = fd;
  fs4->vol = NULL;
  *natoms = MOLFILE_NUMATOMS_NONE;
  fs4->nsets = 1; // this file contains only one data set
  fs4->swap = swap;
  fs4->scale = scale;
  if (norn == 0) {
    // x fast, z medium, y slow
    fs4->x = 0;
    fs4->y = 2;
    fs4->z = 1;
    fs4->f = 0;
    fs4->m = 2;
    fs4->s = 1;
  }
  else if (norn == 1) {
    // y fast, z medium, x slow
    fs4->x = 2;
    fs4->y = 0;
    fs4->z = 1;
    fs4->f = 1;
    fs4->m = 2;
    fs4->s = 0;
  }
  else { // norn ==2
    // x fast, y medium, z slow
    fs4->x = 0;
    fs4->y = 1;
    fs4->z = 2;
    fs4->f = 0;
    fs4->m = 1;
    fs4->s = 2;
  }

  fs4->vol = new molfile_volumetric_t[1];
  strcpy(fs4->vol[0].dataname, "Fsfour Electron Density Map");

  // Place the origin at {0 0 0}
  fs4->vol[0].origin[0] = 0.0;
  fs4->vol[0].origin[1] = 0.0;
  fs4->vol[0].origin[2] = 0.0;

  fs4->vol[0].xaxis[0] = fmsCellSize[0];
  fs4->vol[0].xaxis[1] = 0.0;
  fs4->vol[0].xaxis[2] = 0.0;

  fs4->vol[0].yaxis[0] = cos(gamma) * fmsCellSize[1];
  fs4->vol[0].yaxis[1] = sin(gamma) * fmsCellSize[1];
  fs4->vol[0].yaxis[2] = 0.0;
 
  z1 = cos(beta);
  z2 = (cos(alpha) - cos(beta)*cos(gamma)) / sin(gamma);
  z3 = sqrt(1.0 - z1*z1 - z2*z2);
  fs4->vol[0].zaxis[0] = z1 * fmsCellSize[2];
  fs4->vol[0].zaxis[1] = z2 * fmsCellSize[2];
  fs4->vol[0].zaxis[2] = z3 * fmsCellSize[2];

  fs4->vol[0].xsize = fmsGridSize[fs4->x];
  fs4->vol[0].ysize = fmsGridSize[fs4->y];
  fs4->vol[0].zsize = fmsGridSize[fs4->z];

  fs4->vol[0].has_color = 0;

  return fs4;
}


static int read_fs4_metadata(void *v, int *nsets, 
  molfile_volumetric_t **metadata) {
  fs4_t *fs4 = (fs4_t *)v;
  *nsets = fs4->nsets; 
  *metadata = fs4->vol;  

  return MOLFILE_SUCCESS;
}


static int read_fs4_data(void *v, int set, float *dstBlock, 
                         float *colorblock) {
  fs4_t *fs4 = (fs4_t *) v;
  int *srcBlock, index;
  int col, row, plane, colSize, rowSize, planeSize;
  int xyzGridSize[3], xyzIndexIncrement[3];

  xyzGridSize[0] = fs4->vol[0].xsize;
  xyzGridSize[1] = fs4->vol[0].ysize;
  xyzGridSize[2] = fs4->vol[0].zsize;
  xyzIndexIncrement[0] = 1;
  xyzIndexIncrement[1] = xyzGridSize[0];
  xyzIndexIncrement[2] = xyzGridSize[0] * xyzGridSize[1];

  colSize = xyzGridSize[fs4->f];
  rowSize = xyzGridSize[fs4->m];
  planeSize = xyzGridSize[fs4->s];

  srcBlock = new int[colSize];

  index = 0;
  for (plane = 0; plane < planeSize; plane++) {
    for (row = 0; row < rowSize; row++) {

      // Read one row of data
      if (fortread_4(srcBlock, colSize, fs4->swap, fs4->fd) != colSize) {
        fprintf(stderr, "Error reading data: incorrect record size.\n");
        delete [] srcBlock;
        return MOLFILE_ERROR;
      }

      for (col = 0; col < colSize; col++) {
        dstBlock[index] = (float) srcBlock[col] / fs4->scale;
        index += xyzIndexIncrement[fs4->f];
      }

      index += xyzIndexIncrement[fs4->m] - colSize*xyzIndexIncrement[fs4->f];
    } // end for (row)

    index += xyzIndexIncrement[fs4->s] - rowSize*xyzIndexIncrement[fs4->m];
  } // end for (plane)

  delete [] srcBlock;
  return MOLFILE_SUCCESS;
}

static void close_fs4_read(void *v) {
  fs4_t *fs4 = (fs4_t *)v;

  fclose(fs4->fd);
  if (fs4->vol != NULL)
    delete [] fs4->vol; 
  delete fs4;
}

/*
 * Initialization stuff here
 */
static molfile_plugin_t plugin = {
  vmdplugin_ABIVERSION,   // ABI version
  MOLFILE_PLUGIN_TYPE, 	  // plugin type
  "fs",                   // short file format description
  "FS4 Density Map",      // pretty file format description
  "Eamon Caddigan",       // author(s)
  0,                      // major version
  5,                      // minor version
  VMDPLUGIN_THREADSAFE,   // is reentrant
  "fs,fs4"                // filename extension
};

VMDPLUGIN_EXTERN int VMDPLUGIN_init(void) { return VMDPLUGIN_SUCCESS; }
VMDPLUGIN_EXTERN int VMDPLUGIN_fini(void) { return VMDPLUGIN_SUCCESS; }
VMDPLUGIN_EXTERN int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  plugin.open_file_read = open_fs4_read;
  plugin.read_volumetric_metadata = read_fs4_metadata;
  plugin.read_volumetric_data = read_fs4_data;
  plugin.close_file_read = close_fs4_read;
  (*cb)(v, (vmdplugin_t *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

#ifdef TEST_FS4_PLUGIN

int main(int argc, char **argv) {
  fs4_t *fs4;
  int natoms;
  char *filetype;
  float *datablock;
  
  while (--argc) {
    ++argv;

    fs4 = (fs4_t*) open_fs4_read(*argv, "fs", &natoms);
    if (!fs4) {
      fprintf(stderr, "open_fs4_read failed for file %s\n", *argv);
      return 1;
    }

    printf("a:\t%f\nb:\t%f\nc:\t%f\n",
           fs4->vol[0].xaxis[0], fs4->vol[0].yaxis[1], fs4->vol[0].zaxis[2]);
    printf("ncol:\t%d\nnrow:\t%d\nnplane:\t%d\n", 
           fs4->vol[0].xsize, fs4->vol[0].ysize, fs4->vol[0].zsize);

    datablock = new float[fs4->vol[0].xsize * fs4->vol[0].ysize * 
                          fs4->vol[0].zsize];

    if (read_fs4_data(fs4, 0, datablock, NULL) != 0) {
      fprintf(stderr, "read_fs4_data failed for file %s\n", *argv);
      return 1;
    }

    delete [] datablock;
  }

  return 0;
}

# endif

