/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_ccp4plugin
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
 *      $RCSfile: ccp4plugin.C,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.29 $       $Date: 2007/08/09 19:58:45 $
 *
 ***************************************************************************/

/*
 * CCP4 electron density map file format description:
 *   http://www.ccp4.ac.uk/html/maplib.html
 *   http://iims.ebi.ac.uk/3dem-mrc-maps/distribution/mrc_maps.txt
 *
 * TODO: Fix translation/scaling problems found when using non-orthogonal
 *       unit cells.
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

#define CCP4HDSIZE 1024

typedef struct {
  FILE *fd;
  int nsets;
  int swap;
  int xyz2crs[3];
  long dataOffset;
  molfile_volumetric_t *vol;
} ccp4_t;

static void *open_ccp4_read(const char *filepath, const char *filetype,
    int *natoms) {
  FILE *fd;
  ccp4_t *ccp4;
  char mapString[4], symData[81];
  int origin[3], extent[3], grid[3], crs2xyz[3], mode, symBytes;
  int swap, i, xIndex, yIndex, zIndex;
  long dataOffset, filesize;
  float cellDimensions[3], cellAngles[3], xaxis[3], yaxis[3], zaxis[3];
  float alpha, beta, gamma, xScale, yScale, zScale, z1, z2, z3;
  
  fd = fopen(filepath, "rb");
  if (!fd) {
    printf("ccp4plugin) Error opening file %s\n", filepath);
    return NULL;
  }

  if ( (fread(extent, sizeof(int), 3, fd) != 3) ||
       (fread(&mode, sizeof(int), 1, fd) != 1) ||
       (fread(origin, sizeof(int), 3, fd) != 3) ||
       (fread(grid, sizeof(int), 3, fd) != 3) ||
       (fread(cellDimensions, sizeof(float), 3, fd) != 3) ||
       (fread(cellAngles, sizeof(float), 3, fd) != 3) ||
       (fread(crs2xyz, sizeof(int), 3, fd) != 3) ) {
    printf("ccp4plugin) Error: Improperly formatted line.\n");
    return NULL;
  }

  // Check the number of bytes used for storing symmetry operators
  fseek(fd, 92, SEEK_SET);
  if ((fread(&symBytes, sizeof(int), 1, fd) != 1) ) {
    printf("ccp4plugin) Error: Failed reading symmetry bytes record.\n");
    return NULL;
  }

  // Check for the string "MAP" at byte 208, indicating a CCP4 file.
  fseek(fd, 208, SEEK_SET);
  if ( (fgets(mapString, 4, fd) == NULL) ||
       (strcmp(mapString, "MAP") != 0) ) {
    printf("ccp4plugin) Error: 'MAP' string missing, not a valid CCP4 file.\n");
    return NULL;
  }

  swap = 0;
  // Check the data type of the file.
  if (mode != 2) {
    // Check if the byte-order is flipped
    swap4_aligned(&mode, 1);
    if (mode != 2) {
      printf("ccp4plugin) Error: Non-real (32-bit float) data types are unsupported.\n");
      return NULL;
    } else {
      swap = 1; // enable byte swapping
    }
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


#if 1
  printf("ccp4plugin) extent: %d x %d x %d\n", extent[0], extent[1], extent[2]);
  printf("ccp4plugin) origin: %d x %d x %d\n", origin[0], origin[1], origin[2]);
  printf("ccp4plugin)   grid: %d x %d x %d\n", grid[0], grid[1], grid[2]);
  printf("ccp4plugin) celldim: %f x %f x %f\n", 
         cellDimensions[0], cellDimensions[1], cellDimensions[2]);
  printf("cpp4plugin)cellangles: %f, %f, %f\n", 
         cellAngles[0], cellAngles[1], cellAngles[2]);
  printf("ccp4plugin) crs2xyz: %d %d %d\n", 
         crs2xyz[0], crs2xyz[1], crs2xyz[2]);
  printf("ccp4plugin) symBytes: %d\n", symBytes);
#endif


  // Check the dataOffset: this fixes the problem caused by files claiming
  // to have symmetry records when they do not.
  fseek(fd, 0, SEEK_END);
  filesize = ftell(fd);
  dataOffset = filesize - 4*(extent[0]*extent[1]*extent[2]);
  if (dataOffset != (CCP4HDSIZE + symBytes)) {
    if (dataOffset == CCP4HDSIZE) {
      // Bogus symmetry record information
      printf("ccp4plugin) Warning: file contains bogus symmetry record.\n");
      symBytes = 0;
    } else if (dataOffset < CCP4HDSIZE) {
      printf("ccp4plugin) Error: File appears truncated and doesn't match header.\n");
      return NULL;
    } else if ((dataOffset > CCP4HDSIZE) && (dataOffset < (1024*1024))) {
      // Fix for loading SPIDER files which are larger than usual
      // In this specific case, we must absolutely trust the symBytes record
      dataOffset = CCP4HDSIZE + symBytes; 
      printf("ccp4plugin) Warning: File is larger than expected and doesn't match header.\n");
      printf("ccp4plugin) Warning: Continuing file load, good luck!\n");
    } else {
      printf("ccp4plugin) Error: File is MUCH larger than expected and doesn't match header.\n");
      return NULL;
    }
  }

  // Read symmetry records -- organized as 80-byte lines of text.
  if (symBytes != 0) {
    printf("ccp4plugin) Symmetry records found:\n");
    fseek(fd, CCP4HDSIZE, SEEK_SET);
    for (i = 0; i < symBytes/80; i++) {
      fgets(symData, 81, fd);
      printf("ccp4plugin) %s\n", symData);
    }
  }

  // check extent and grid interval counts
  if (grid[0] == 0 && extent[0] > 0) {
    grid[0] = extent[0] - 1;
    printf("ccp4plugin) Warning: Fixed X interval count\n");
  }
  if (grid[1] == 0 && extent[1] > 0) {
    grid[1] = extent[1] - 1;
    printf("ccp4plugin) Warning: Fixed Y interval count\n");
  }
  if (grid[2] == 0 && extent[2] > 0) {
    grid[2] = extent[2] - 1;
    printf("ccp4plugin) Warning: Fixed Z interval count\n");
  }
 

  // Allocate and initialize the ccp4 structure
  ccp4 = new ccp4_t;
  ccp4->fd = fd;
  ccp4->vol = NULL;
  *natoms = MOLFILE_NUMATOMS_NONE;
  ccp4->nsets = 1; // this EDM file contains only one data set
  ccp4->swap = swap;
  ccp4->dataOffset = dataOffset;

  ccp4->vol = new molfile_volumetric_t[1];
  strcpy(ccp4->vol[0].dataname, "CCP4 Electron Density Map");

  // Mapping between CCP4 column, row, section and VMD x, y, z.
  if (crs2xyz[0] == 0 && crs2xyz[1] == 0 && crs2xyz[2] == 0) {
    printf("ccp4plugin) Warning: All crs2xyz records are zero.\n");
    printf("ccp4plugin) Warning: Setting crs2xyz to 1, 2, 3\n");
    crs2xyz[0] = 1;
    crs2xyz[0] = 2;
    crs2xyz[0] = 3;
  }

  ccp4->xyz2crs[crs2xyz[0]-1] = 0;
  ccp4->xyz2crs[crs2xyz[1]-1] = 1;
  ccp4->xyz2crs[crs2xyz[2]-1] = 2;
  xIndex = ccp4->xyz2crs[0];
  yIndex = ccp4->xyz2crs[1];
  zIndex = ccp4->xyz2crs[2];

  // calculate non-orthogonal unit cell coordinates
  alpha = (M_PI / 180.0) * cellAngles[0];
  beta = (M_PI / 180.0) * cellAngles[1];
  gamma = (M_PI / 180.0) * cellAngles[2];

  if (cellDimensions[0] == 0.0 && 
      cellDimensions[1] == 0.0 &&
      cellDimensions[2] == 0.0) {
    printf("ccp4plugin) Warning: Cell dimensions are all zero.\n");
    printf("ccp4plugin) Warning: Setting to 1.0, 1.0, 1.0 for viewing.\n");
    printf("ccp4plugin) Warning: Map file will not align with other structures.\n");
    cellDimensions[0] = 1.0;
    cellDimensions[1] = 1.0;
    cellDimensions[2] = 1.0;
  } 


  xScale = cellDimensions[0] / grid[0];
  yScale = cellDimensions[1] / grid[1];
  zScale = cellDimensions[2] / grid[2];

  // calculate non-orthogonal unit cell coordinates
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

  ccp4->vol[0].origin[0] = xaxis[0] * origin[xIndex] + 
                           yaxis[0] * origin[yIndex] +
                           zaxis[0] * origin[zIndex];
  ccp4->vol[0].origin[1] = yaxis[1] * origin[yIndex] +
                           zaxis[1] * origin[zIndex];
  ccp4->vol[0].origin[2] = zaxis[2] * origin[zIndex];

  ccp4->vol[0].xaxis[0] = xaxis[0] * (extent[xIndex]-1);
  ccp4->vol[0].xaxis[1] = 0;
  ccp4->vol[0].xaxis[2] = 0;

  ccp4->vol[0].yaxis[0] = yaxis[0] * (extent[yIndex]-1);
  ccp4->vol[0].yaxis[1] = yaxis[1] * (extent[yIndex]-1);
  ccp4->vol[0].yaxis[2] = 0;

  ccp4->vol[0].zaxis[0] = zaxis[0] * (extent[zIndex]-1);
  ccp4->vol[0].zaxis[1] = zaxis[1] * (extent[zIndex]-1);
  ccp4->vol[0].zaxis[2] = zaxis[2] * (extent[zIndex]-1);

  ccp4->vol[0].xsize = extent[xIndex];
  ccp4->vol[0].ysize = extent[yIndex];
  ccp4->vol[0].zsize = extent[zIndex];

  ccp4->vol[0].has_color = 0;

  return ccp4;
}

static int read_ccp4_metadata(void *v, int *nsets, 
  molfile_volumetric_t **metadata) {
  ccp4_t *ccp4 = (ccp4_t *)v;
  *nsets = ccp4->nsets; 
  *metadata = ccp4->vol;  

  return MOLFILE_SUCCESS;
}

static int read_ccp4_data(void *v, int set, float *datablock,
                         float *colorblock) {
  ccp4_t *ccp4 = (ccp4_t *)v;
  float *rowdata;
  int x, y, z, xSize, ySize, zSize, xySize, extent[3], coord[3];
  FILE *fd = ccp4->fd;

  xSize = ccp4->vol[0].xsize;
  ySize = ccp4->vol[0].ysize;
  zSize = ccp4->vol[0].zsize;
  xySize = xSize * ySize;

  // coord = <col, row, sec>
  // extent = <colSize, rowSize, secSize>
  extent[ccp4->xyz2crs[0]] = xSize;
  extent[ccp4->xyz2crs[1]] = ySize;
  extent[ccp4->xyz2crs[2]] = zSize;

  rowdata = new float[extent[0]];

  fseek(fd, ccp4->dataOffset, SEEK_SET);

  for (coord[2] = 0; coord[2] < extent[2]; coord[2]++) {
    for (coord[1] = 0; coord[1] < extent[1]; coord[1]++) {
      // Read an entire row of data from the file, then write it into the
      // datablock with the correct slice ordering.
      if (feof(fd)) {
        printf("ccp4plugin) Unexpected end-of-file.\n");
        return MOLFILE_ERROR;
      }
      if (ferror(fd)) {
        printf("ccp4plugin) Problem reading the file.\n");
        return MOLFILE_ERROR;
      }
      if ( fread(rowdata, sizeof(float), extent[0], fd) != extent[0] ) {
        printf("ccp4plugin) Error reading data row.\n");
        return MOLFILE_ERROR;
      }

      for (coord[0] = 0; coord[0] < extent[0]; coord[0]++) {
        x = coord[ccp4->xyz2crs[0]];
        y = coord[ccp4->xyz2crs[1]];
        z = coord[ccp4->xyz2crs[2]];
        datablock[x + y*xSize + z*xySize] = rowdata[coord[0]];
      }
    }
  }

  if (ccp4->swap == 1)
    swap4_aligned(datablock, xySize * zSize);

  delete [] rowdata;

  return MOLFILE_SUCCESS;
}

static void close_ccp4_read(void *v) {
  ccp4_t *ccp4 = (ccp4_t *)v;

  fclose(ccp4->fd);
  delete [] ccp4->vol; 
  delete ccp4;
}

/*
 * Initialization stuff here
 */
static molfile_plugin_t plugin;

VMDPLUGIN_EXTERN int VMDPLUGIN_init(void) { 
  memset(&plugin, 0, sizeof(molfile_plugin_t));
  plugin.abiversion = vmdplugin_ABIVERSION;
  plugin.type = MOLFILE_PLUGIN_TYPE;
  plugin.name = "ccp4";
  plugin.prettyname = "CCP4, MRC Density Map";
  plugin.author = "Eamon Caddigan, John Stone";
  plugin.majorv = 1;
  plugin.minorv = 2;
  plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
  plugin.filename_extension = "ccp4,mrc,map";
  plugin.open_file_read = open_ccp4_read;
  plugin.read_volumetric_metadata = read_ccp4_metadata;
  plugin.read_volumetric_data = read_ccp4_data;
  plugin.close_file_read = close_ccp4_read;
  return VMDPLUGIN_SUCCESS; 
}


VMDPLUGIN_EXTERN int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)(void *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_EXTERN int VMDPLUGIN_fini(void) { return VMDPLUGIN_SUCCESS; }
