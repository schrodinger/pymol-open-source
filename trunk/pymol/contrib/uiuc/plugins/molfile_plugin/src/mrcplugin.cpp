/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_mrcplugin
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
 *      $RCSfile: mrcplugin.C,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.7 $       $Date: 2006/02/23 19:36:45 $
 *
 ***************************************************************************/

/*
 * MRC EM map file format: 
 *   http://iims.ebi.ac.uk/3dem-mrc-maps/distribution/mrc_maps.txt
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "molfile_plugin.h"
#include "endianswap.h"

typedef struct {
  FILE *fd;
  int nsets;
  int swap;
  int xyz2crs[3];
  long dataOffset;
  molfile_volumetric_t *vol;
} mrc_t;

static void *open_mrc_read(const char *filepath, const char *filetype,
    int *natoms) {
  FILE *fd;
  mrc_t *mrc;
  char mapString[5], symData[82];
  int origin[3], extent[3], grid[3], crs2xyz[3], mode, symBytes;
  int swap, i, xIndex, yIndex, zIndex;
  long dataOffset;
  float cellDimensions[3], cellAngles[3], xaxis[3], yaxis[3], zaxis[3];
  float alpha, beta, gamma, xScale, yScale, zScale, z1, z2, z3;
  
  fd = fopen(filepath, "rb");
  if (!fd) {
    fprintf(stderr, "Error opening file %s\n", filepath);
    return NULL;
  }

  if ( (fread(extent, sizeof(int), 3, fd) != 3) ||
       (fread(&mode, sizeof(int), 1, fd) != 1) ||
       (fread(origin, sizeof(int), 3, fd) != 3) ||
       (fread(grid, sizeof(int), 3, fd) != 3) ||
       (fread(cellDimensions, sizeof(float), 3, fd) != 3) ||
       (fread(cellAngles, sizeof(float), 3, fd) != 3) ||
       (fread(crs2xyz, sizeof(int), 3, fd) != 3) ) {
    fprintf(stderr, "Improperly formatted line.\n");
    return NULL;
  }

  // Check the number of bytes used for storing symmetry operators
  fseek(fd, 92, SEEK_SET);
  if ( (fread(&symBytes, sizeof(int), 1, fd) != 1) ) {
    fprintf(stderr, "Problem reading the file.\n");
    return NULL;
  }

  // Check for the string "MAP" at byte 208, indicating a CCP4 or MRC file.
  fseek(fd, 208, SEEK_SET);
  if ((fread(mapString, 1, 3, fd) != 3) || (strcmp(mapString, "MAP") != 0) ) {
    fprintf(stderr, "File not in MRC or CCP4 format.\n");
    return NULL;
  }

// XXX
printf("MRC file mode: %d\n", mode);

  swap = 0;
  // Check the data type of the file.
  if (mode != 2) {
    // Check if the byte-order is flipped
    swap4_aligned(&mode, 1);
    if (mode != 2) {
      fprintf(stderr, "Non-real data types are currently not supported.\n");
      return NULL;
    }
    else
      swap = 1;
  }

  // Swap all the information obtained from the header
  if (swap == 1) {
    swap4_aligned(extent, 3);
    swap4_aligned(origin, 3);
    swap4_aligned(grid, 3);
    swap4_aligned(cellDimensions, 3);
    swap4_aligned(cellAngles, 3);
    swap4_aligned(crs2xyz, 3);
    swap4_aligned(&symBytes, 1);
  }

// XXX
printf("MRC file mode: %d\n", mode);
printf("MRC stats:\n");
printf("  extent: %d %d %d\n", extent[0], extent[1], extent[2]);
printf("  origin: %d %d %d\n", origin[0], origin[1], origin[2]);
printf("    grid: %d %d %d\n",   grid[0],   grid[1],   grid[2]);
printf(" crs2xyz: %d %d %d\n", crs2xyz[0], crs2xyz[1], crs2xyz[2]);
printf(" celldim: %f %f %f\n", cellDimensions[0], cellDimensions[1], cellDimensions[2]);
printf(" cellang: %f %f %f\n", cellAngles[0], cellAngles[1], cellAngles[2]);
printf(" symbytes: %d\n", symBytes);

#if 1
  // Check the dataOffset: this fixes the problem caused by files claiming
  // to have symmetry records when they do not.
  fseek(fd, 0, SEEK_END);
  int filesize = ftell(fd);
  dataOffset = filesize - 4*(extent[0]*extent[1]*extent[2]);

printf("filesize: %d\n", filesize);
printf("expected file size = (%d - %d) = %d\n", filesize, 4*(extent[0]*extent[1]*extent[2]), filesize - 4*(extent[0]*extent[1]*extent[2]));
printf("dataOffset: %ld\n", dataOffset);

  if (dataOffset != (1024 + symBytes)) {
    if (dataOffset == 1024) {
      // Bogus symmetry record information
      fprintf(stdout, "Warning: file contains bogus symmetry record.\n");
      symBytes = 0;
    }
    else {
      fprintf(stderr, "File size does not match header.\n");
      return NULL;
    }
  }
#endif

  // Read symmetry records -- organized as 80-byte lines of text.
  if (symBytes != 0) {
    fprintf(stdout, "Symmetry records found:\n");
    fseek(fd, 1024, SEEK_SET);
    for (i = 0; i < symBytes/80; i++) {
      fgets(symData, 81, fd);
      fprintf(stdout, "%s\n", symData);
    }
  }

  // check extent and grid interval counts
  if (grid[0] == 0 && extent[0] != 0) {
    grid[0] = extent[0] - 1;
    printf("mrcplugin) Warning: Fixed X interval count\n");
  }
  if (grid[1] == 0 && extent[1] != 0) {
    grid[1] = extent[1] - 1;
    printf("mrcplugin) Warning: Fixed Y interval count\n");
  }
  if (grid[2] == 0 && extent[2] != 0) {
    grid[2] = extent[2] - 1;
    printf("mrcplugin) Warning: Fixed Z interval count\n");
  } 
  xScale = cellDimensions[0] / ((float) grid[0]);
  yScale = cellDimensions[1] / ((float) grid[1]);
  zScale = cellDimensions[2] / ((float) grid[2]);
printf(" scale: %f %f %f\n", xScale, yScale, zScale);

  // Allocate and initialize the mrc structure
  mrc = new mrc_t;
  mrc->fd = fd;
  mrc->vol = NULL;
  *natoms = MOLFILE_NUMATOMS_NONE;
  mrc->nsets = 1; // this EDM file contains only one data set
  mrc->swap = swap;
  mrc->dataOffset = dataOffset;

  mrc->vol = new molfile_volumetric_t[1];
  strcpy(mrc->vol[0].dataname, "MRC EM Map");

  // Mapping between MRC column, row, section and VMD x, y, z.
  mrc->xyz2crs[crs2xyz[0]-1] = 0;
  mrc->xyz2crs[crs2xyz[1]-1] = 1;
  mrc->xyz2crs[crs2xyz[2]-1] = 2;
  xIndex = mrc->xyz2crs[0];
  yIndex = mrc->xyz2crs[1];
  zIndex = mrc->xyz2crs[2];

  // calculate non-orthogonal unit cell coordinates
  alpha = (M_PI / 180.0) * cellAngles[0];
  beta = (M_PI / 180.0) * cellAngles[1];
  gamma = (M_PI / 180.0) * cellAngles[2];

  xaxis[0] = xScale;
  xaxis[1] = 0;
  xaxis[2] = 0;

  yaxis[0] = cos(gamma) * yScale;
  yaxis[1] = sin(gamma) * yScale;
  yaxis[2] = 0;

  z1 = cos(beta);
  z2 = (cos(alpha) - cos(beta)*cos(gamma)) / sin(gamma);
  z3 = sqrt(1.0 - z1*z1 - z2*z2);
  zaxis[0] = z1 * zScale;
  zaxis[1] = z2 * zScale;
  zaxis[2] = z3 * zScale;

  mrc->vol[0].origin[0] = xaxis[0] * origin[xIndex] + 
                           yaxis[0] * origin[yIndex] +
                           zaxis[0] * origin[zIndex];
  mrc->vol[0].origin[1] = yaxis[1] * origin[yIndex] +
                           zaxis[1] * origin[zIndex];
  mrc->vol[0].origin[2] = zaxis[2] * origin[zIndex];

  mrc->vol[0].xaxis[0] = xaxis[0] * (extent[xIndex]-1);
  mrc->vol[0].xaxis[1] = 0;
  mrc->vol[0].xaxis[2] = 0;

  mrc->vol[0].yaxis[0] = yaxis[0] * (extent[yIndex]-1);
  mrc->vol[0].yaxis[1] = yaxis[1] * (extent[yIndex]-1);
  mrc->vol[0].yaxis[2] = 0;

  mrc->vol[0].zaxis[0] = zaxis[0] * (extent[zIndex]-1);
  mrc->vol[0].zaxis[1] = zaxis[1] * (extent[zIndex]-1);
  mrc->vol[0].zaxis[2] = zaxis[2] * (extent[zIndex]-1);

  mrc->vol[0].xsize = extent[xIndex];
  mrc->vol[0].ysize = extent[yIndex];
  mrc->vol[0].zsize = extent[zIndex];

  mrc->vol[0].has_color = 0;

printf("origin: %f %f %f\n", origin[0], origin[1], origin[2]);
printf(" xaxis: %f %f %f\n", xaxis[0], xaxis[1], xaxis[2]);
printf(" yaxis: %f %f %f\n", yaxis[0], yaxis[1], yaxis[2]);
printf(" zaxis: %f %f %f\n", zaxis[0], zaxis[1], zaxis[2]);

  return mrc;
}

static int read_mrc_metadata(void *v, int *nsets, 
  molfile_volumetric_t **metadata) {
  mrc_t *mrc = (mrc_t *)v;
  *nsets = mrc->nsets; 
  *metadata = mrc->vol;  

  return MOLFILE_SUCCESS;
}

static int read_mrc_data(void *v, int set, float *datablock,
                         float *colorblock) {
  mrc_t *mrc = (mrc_t *)v;
  float *rowdata;
  int x, y, z, xSize, ySize, zSize, xySize, extent[3], coord[3];
  FILE *fd = mrc->fd;

  xSize = mrc->vol[0].xsize;
  ySize = mrc->vol[0].ysize;
  zSize = mrc->vol[0].zsize;
  xySize = xSize * ySize;

  // coord = <col, row, sec>
  // extent = <colSize, rowSize, secSize>
  extent[mrc->xyz2crs[0]] = xSize;
  extent[mrc->xyz2crs[1]] = ySize;
  extent[mrc->xyz2crs[2]] = zSize;

  rowdata = new float[extent[0]];

  fseek(fd, mrc->dataOffset, SEEK_SET);

  for (coord[2] = 0; coord[2] < extent[2]; coord[2]++) {
    for (coord[1] = 0; coord[1] < extent[1]; coord[1]++) {
      // Read an entire row of data from the file, then write it into the
      // datablock with the correct slice ordering.
      if (feof(fd)) {
        fprintf(stderr, "Unexpected end-of-file.\n");
        return MOLFILE_ERROR;
      }
      if (ferror(fd)) {
        fprintf(stderr, "Problem reading the file.\n");
        return MOLFILE_ERROR;
      }
      if ( fread(rowdata, sizeof(float), extent[0], fd) != extent[0] ) {
        fprintf(stderr, "Error reading data row %d/%d\n", coord[1], extent[1]);
        return MOLFILE_ERROR;
      }

      for (coord[0] = 0; coord[0] < extent[0]; coord[0]++) {
        x = coord[mrc->xyz2crs[0]];
        y = coord[mrc->xyz2crs[1]];
        z = coord[mrc->xyz2crs[2]];
        datablock[x + y*xSize + z*xySize] = rowdata[coord[0]];
      }
    }
  }

  if (mrc->swap == 1)
    swap4_aligned(datablock, xySize * zSize);

  delete [] rowdata;

  return MOLFILE_SUCCESS;
}

static void close_mrc_read(void *v) {
  mrc_t *mrc = (mrc_t *)v;

  fclose(mrc->fd);
  delete [] mrc->vol; 
  delete mrc;
}

/*
 * Initialization stuff here
 */
static molfile_plugin_t plugin = {
  vmdplugin_ABIVERSION,   // ABI version
  MOLFILE_PLUGIN_TYPE, 	  // plugin type
  "mrc",                  // file format description
  "MRC",                  // file format description
  "Eamon Caddigan",       // author(s)
  0,                      // major version
  5,                      // minor version
  VMDPLUGIN_THREADSAFE,   // is reentrant
  "mrc"                   // filename extension
};

VMDPLUGIN_EXTERN int VMDPLUGIN_init(void) { return VMDPLUGIN_SUCCESS; }
VMDPLUGIN_EXTERN int VMDPLUGIN_fini(void) { return VMDPLUGIN_SUCCESS; }
VMDPLUGIN_EXTERN int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  plugin.open_file_read = open_mrc_read;
  plugin.read_volumetric_metadata = read_mrc_metadata;
  plugin.read_volumetric_data = read_mrc_data;
  plugin.close_file_read = close_mrc_read;
  (*cb)(v, (vmdplugin_t *)(void *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

