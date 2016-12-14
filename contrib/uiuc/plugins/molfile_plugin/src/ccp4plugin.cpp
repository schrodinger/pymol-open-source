/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_ccp4plugin
#define STATIC_PLUGIN 1

/***************************************************************************
 *cr
 *cr            (C) Copyright 1995-2016 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: ccp4plugin.C,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.38 $       $Date: 2016/11/28 05:01:53 $
 *
 ***************************************************************************/

//
// CCP4 electron density map file format description:
//   http://www2.mrc-lmb.cam.ac.uk/image2000.html
//   http://www.ccp4.ac.uk/html/maplib.html
//   http://iims.ebi.ac.uk/3dem-mrc-maps/distribution/mrc_maps.txt
//
// TODO: Fix translation/scaling problems found when using 
//       non-orthogonal unit cells.
//
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

// CCP4 and IMOD MRC image types
// based on notes in:
//   http://bio3d.colorado.edu/imod/doc/mrc_format.txt
#define MRC_TYPE_BYTE   0
#define MRC_TYPE_SHORT  1
#define MRC_TYPE_FLOAT  2
#define MRC_TYPE_SHORT2 3
#define MRC_TYPE_FLOAT2 4
#define MRC_TYPE_USHORT 6  /* non-standard */
#define MRC_TYPE_UCHAR3 16 /* non-standard */

#define IMOD_FILE_ID            0x494d4f44
#define IMOD_MAGIC_STAMP        1146047817
#define IMOD_FLAG_SIGNED               0x1
#define IMOD_FLAG_HEADER_SPACING       0x2
#define IMOD_FLAG_ORIGIN_INVERTED_SIGN 0x4

typedef struct {
  FILE *fd;
  int voxtype; 
  int imodstamp;
  int imodflags;
  int nsets;
  int swap;
  int xyz2crs[3];
  long dataOffset;
  molfile_volumetric_t *vol;
} ccp4_t;

typedef struct {
  unsigned char red;
  unsigned char blue;
  unsigned char green;
} uchar3;

static void *open_ccp4_read(const char *filepath, const char *filetype,
    int *natoms) {
  FILE *fd;
  ccp4_t *ccp4;
  char mapString[4], symData[81];
  int nxyzstart[3], extent[3], grid[3], crs2xyz[3], voxtype, symBytes;
  float origin2k[3];
  int swap, i, xIndex, yIndex, zIndex;
  long dataOffset, filesize;
  float cellDimensions[3], cellAngles[3], xaxis[3], yaxis[3], zaxis[3];
  float alpha, beta, gamma, xScale, yScale, zScale, z1, z2, z3;
  int imodstamp=0,imodflags=0;
  
  fd = fopen(filepath, "rb");
  if (!fd) {
    printf("ccp4plugin) Error opening file %s\n", filepath);
    return NULL;
  }

  if ((fread(extent, sizeof(int), 3, fd) != 3) ||
      (fread(&voxtype, sizeof(int), 1, fd) != 1) ||
      (fread(nxyzstart, sizeof(int), 3, fd) != 3) ||
      (fread(grid, sizeof(int), 3, fd) != 3) ||
      (fread(cellDimensions, sizeof(float), 3, fd) != 3) ||
      (fread(cellAngles, sizeof(float), 3, fd) != 3) ||
      (fread(crs2xyz, sizeof(int), 3, fd) != 3)) {
    printf("ccp4plugin) Error: Improperly formatted line.\n");
    return NULL;
  }

  // Check the number of bytes used for storing symmetry operators
  // (word 23, byte 92)
  fseek(fd, 23 * 4, SEEK_SET);
  if ((fread(&symBytes, sizeof(int), 1, fd) != 1) ) {
    printf("ccp4plugin) Error: Failed reading symmetry bytes record.\n");
    return NULL;
  }

  // read MRC2000 Origin record at word 49, byte 196, and use it if necessary
  // http://www2.mrc-lmb.cam.ac.uk/image2000.html
  fseek(fd, 49 * 4, SEEK_SET);
  if (fread(origin2k, sizeof(float), 3, fd) != 3) {
    printf("ccp4plugin) Error: unable to read ORIGIN records at offset 196.\n");
  }

  // Check for IMOD stamp at offset 152, indicating an IMOD file 
  fseek(fd, 152, SEEK_SET);
  if (fread(&imodstamp, sizeof(int), 1, fd) != 1) {
    printf("ccp4plugin) Error: failed to read IMOD stamp from MRC file.\n");
    return NULL;
  }
  if (fread(&imodflags, sizeof(int), 1, fd) != 1) {
    printf("ccp4plugin) Error: failed to read IMOD flags from MRC file.\n");
    return NULL;
  }

  // Check file endianism using some heuristics
  swap = 0;
  int tmp[3];
  memcpy(tmp, extent, 3*sizeof(int));
  if (tmp[0] > 65536 || tmp[1] > 65536 || tmp[2] > 65536) {
    swap4_aligned(tmp, 3);
    if (tmp[0] > 65536 || tmp[1] > 65536 || tmp[2] > 65536) {
      swap = 0;
      printf("ccp4plugin) Guessing file endianism: native\n");
    } else {
      swap = 1;
      printf("ccp4plugin) Guessing file endianism: swapped\n");
    }
  }
  if (voxtype > 16 && swap == 0) {
    tmp[0] = voxtype;
    swap4_aligned(tmp, 1);
    if (tmp[0] <= 16) {
      swap = 1;
      printf("ccp4plugin) Guessing file endianism: swapped\n");
    }
  }

  // Byte-swap header information if needed
  if (swap == 1) {
    swap4_aligned(extent, 3);
    swap4_aligned(&voxtype, 1);
    swap4_aligned(nxyzstart, 3);
    swap4_aligned(origin2k, 3);
    swap4_aligned(grid, 3);
    swap4_aligned(cellDimensions, 3);
    swap4_aligned(cellAngles, 3);
    swap4_aligned(crs2xyz, 3);
    swap4_aligned(&symBytes, 1);
    swap4_aligned(&imodstamp, 1);
    swap4_aligned(&imodflags, 1);
  }

  // Check for the string "MAP" at word 52 byte 208, indicating a CCP4 file.
  fseek(fd, 52 * 4, SEEK_SET);
  if (fgets(mapString, 4, fd) == NULL) {
    printf("ccp4plugin) Error: unable to read 'MAP' string, not a valid CCP4/IMOD MRC file.\n");
    return NULL;
  }
  if ((strcmp(mapString, "MAP") != 0) && (imodstamp != IMOD_MAGIC_STAMP)) {
    //Older versions of IMOD (2.6.19 or below) do not have the "MAP " string.
    //If the IMOD stamp is there its probably a valid mrc file
    printf("ccp4plugin) Warning: 'MAP' string missing which usually indicates that this is\n\
not a valid IMOD file. Some older versions of IMOD did not include the 'MAP'\n\
string so file loading will continue but may fail.\n");
  }

  // Check if we found an IMOD file or not
  if (imodstamp == IMOD_MAGIC_STAMP) {
    printf("ccp4plugin) MRC file generated by IMOD-compatible program.\n");
    if (imodflags & IMOD_FLAG_SIGNED)
      printf("ccp4plugin) IMOD flag: data uses signed-bytes\n");
    else 
      printf("ccp4plugin) IMOD flag: data uses unsigned-bytes\n");

    if (imodflags & IMOD_FLAG_HEADER_SPACING)
      printf("ccp4plugin) IMOD flag: pixel spacing set in extended header\n");

    if (imodflags & IMOD_FLAG_ORIGIN_INVERTED_SIGN)
      printf("ccp4plugin) IMOD flag: origin sign is inverted.\n");
  } else {
    printf("ccp4plugin) No IMOD stamp found.\n");
    imodflags = 0;
  }

  // Check the data type of the file.
  switch (voxtype) {
    case MRC_TYPE_BYTE:
      printf("ccp4plugin) voxel type: byte\n");
      break;

    case MRC_TYPE_SHORT:
      printf("ccp4plugin) voxel type: short (16-bit signed int)\n");
      break;

    case MRC_TYPE_FLOAT:
      printf("ccp4plugin) voxel type: float (32-bit real)\n");
      break;

    case MRC_TYPE_SHORT2:
      printf("ccp4plugin) voxel type: short2 (2x 16-bit signed int)\n");
      printf("ccp4plugin) Error: unimplemented voxel format\n");
      return NULL;

    case MRC_TYPE_FLOAT2:
      printf("ccp4plugin) voxel type: float2 (2x 32-bit real)\n");
      printf("ccp4plugin) Error: unimplemented voxel format\n");
      return NULL;

    case MRC_TYPE_USHORT:
      printf("ccp4plugin) voxel type: ushort (16-bit unsigned int)\n");
      break;

    case MRC_TYPE_UCHAR3:
      printf("ccp4plugin) voxel type: uchar3 (3x unsigned char)\n");
      break;

    default:
      printf("ccp4plugin) Error: Only byte, short (16-bit integer) or float (32-bit real) data types are supported.\n");
      return NULL;
  }

#if 1
  printf("ccp4plugin)    extent: %d x %d x %d\n",
         extent[0], extent[1], extent[2]);
  printf("ccp4plugin) nxyzstart: %d x %d x %d\n", 
         nxyzstart[0], nxyzstart[1], nxyzstart[2]);
  printf("ccp4plugin)  origin2k: %f x %f x %f\n", 
         origin2k[0], origin2k[1], origin2k[2]);
  printf("ccp4plugin)      grid: %d x %d x %d\n", grid[0], grid[1], grid[2]);
  printf("ccp4plugin)   celldim: %f x %f x %f\n", 
         cellDimensions[0], cellDimensions[1], cellDimensions[2]);
  printf("cpp4plugin)cellangles: %f, %f, %f\n", 
         cellAngles[0], cellAngles[1], cellAngles[2]);
  printf("ccp4plugin)   crs2xyz: %d %d %d\n", 
         crs2xyz[0], crs2xyz[1], crs2xyz[2]);
  printf("ccp4plugin)  symBytes: %d\n", symBytes);
#endif

  // Check the dataOffset: this fixes the problem caused by files claiming
  // to have symmetry records when they do not.
  fseek(fd, 0, SEEK_END);
  filesize = ftell(fd);

  // compute data offset using file size and voxel type info
  if (voxtype == MRC_TYPE_BYTE) {
    dataOffset = filesize - sizeof(char)*(extent[0]*extent[1]*extent[2]);
  } else if (voxtype == MRC_TYPE_FLOAT) {
    dataOffset = filesize - sizeof(float)*(extent[0]*extent[1]*extent[2]);
  } else if (voxtype == MRC_TYPE_SHORT || voxtype == MRC_TYPE_USHORT) {
    dataOffset = filesize - sizeof(short)*(extent[0]*extent[1]*extent[2]);
  } else if (voxtype == MRC_TYPE_UCHAR3) {
    dataOffset = filesize - sizeof(uchar3)*(extent[0]*extent[1]*extent[2]);
  } else {
    printf("ccp4plugin) unimplemented voxel type!\n");
  }

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
  ccp4->voxtype = voxtype;
  ccp4->imodstamp = imodstamp;
  ccp4->imodflags = imodflags;
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

#if 1
  // Handle both MRC-2000 and older format maps
  if (origin2k[0] == 0.0f && origin2k[1] == 0.0f && origin2k[2] == 0.0f) {
    printf("ccp4plugin) using CCP4 n[xyz]start origin\n");
    ccp4->vol[0].origin[0] = xaxis[0] * nxyzstart[xIndex] + 
                             yaxis[0] * nxyzstart[yIndex] +
                             zaxis[0] * nxyzstart[zIndex];
    ccp4->vol[0].origin[1] = yaxis[1] * nxyzstart[yIndex] +
                             zaxis[1] * nxyzstart[zIndex];
    ccp4->vol[0].origin[2] = zaxis[2] * nxyzstart[zIndex];
  } else {
    // Use ORIGIN records rather than old n[xyz]start records
    //   http://www2.mrc-lmb.cam.ac.uk/image2000.html
    // XXX the ORIGIN field is only used by the EM community, and
    //     has undefined meaning for non-orthogonal maps and/or
    //     non-cubic voxels, etc.
    printf("ccp4plugin) using MRC2000 origin\n");
    ccp4->vol[0].origin[0] = origin2k[xIndex];
    ccp4->vol[0].origin[1] = origin2k[yIndex];
    ccp4->vol[0].origin[2] = origin2k[zIndex];
  }
#else
  // old code that only pays attention to old MRC nxstart/nystart/nzstart
  ccp4->vol[0].origin[0] = xaxis[0] * nxyzstart[xIndex] + 
                           yaxis[0] * nxyzstart[yIndex] +
                           zaxis[0] * nxyzstart[zIndex];
  ccp4->vol[0].origin[1] = yaxis[1] * nxyzstart[yIndex] +
                           zaxis[1] * nxyzstart[zIndex];
  ccp4->vol[0].origin[2] = zaxis[2] * nxyzstart[zIndex];
#endif

#if 0
  printf("ccp4plugin) origin: %.3f %.3f %.3f\n",
         ccp4->vol[0].origin[0],
         ccp4->vol[0].origin[1],
         ccp4->vol[0].origin[2]);
#endif

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

  fseek(fd, ccp4->dataOffset, SEEK_SET);

  // Read entire rows of data from the file, then write into the
  // datablock with the correct slice ordering.
  if ((ccp4->voxtype == MRC_TYPE_BYTE) && (ccp4->imodflags & IMOD_FLAG_SIGNED)) {
    printf("ccp4plugin) reading signed-byte voxel data\n");
    char *rowdata = new char[extent[0]];
    for (coord[2] = 0; coord[2] < extent[2]; coord[2]++) {
      for (coord[1] = 0; coord[1] < extent[1]; coord[1]++) {
        if (feof(fd)) {
          printf("ccp4plugin) Unexpected end-of-file.\n");
          delete [] rowdata;
          return MOLFILE_ERROR;
        }
        if (ferror(fd)) {
          printf("ccp4plugin) Problem reading the file.\n");
          delete [] rowdata;
          return MOLFILE_ERROR;
        }
        if ( fread(rowdata, sizeof(char), extent[0], fd) != extent[0] ) {
          printf("ccp4plugin) Error reading data row.\n");
          delete [] rowdata;
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
    delete [] rowdata;
  } else if ((ccp4->voxtype == MRC_TYPE_BYTE) && !(ccp4->imodflags & IMOD_FLAG_SIGNED)) {
    printf("ccp4plugin) reading unsigned-byte voxel data\n");
    unsigned char *rowdata = new unsigned char[extent[0]];
    for (coord[2] = 0; coord[2] < extent[2]; coord[2]++) {
      for (coord[1] = 0; coord[1] < extent[1]; coord[1]++) {
        if (feof(fd)) {
          printf("ccp4plugin) Unexpected end-of-file.\n");
          delete [] rowdata;
          return MOLFILE_ERROR;
        }
        if (ferror(fd)) {
          printf("ccp4plugin) Problem reading the file.\n");
          delete [] rowdata;
          return MOLFILE_ERROR;
        }
        if ( fread(rowdata, sizeof(unsigned char), extent[0], fd) != extent[0] ) {
          printf("ccp4plugin) Error reading data row.\n");
          delete [] rowdata;
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
    delete [] rowdata;
  } else if (ccp4->voxtype == MRC_TYPE_FLOAT) {
    printf("ccp4plugin) reading float (32-bit real) voxel data\n");
    float *rowdata = new float[extent[0]];
    int x, y, z;
    for (coord[2] = 0; coord[2] < extent[2]; coord[2]++) {
      for (coord[1] = 0; coord[1] < extent[1]; coord[1]++) {
        if (feof(fd)) {
          printf("ccp4plugin) Unexpected end-of-file.\n");
          delete [] rowdata;
          return MOLFILE_ERROR;
        }
        if (ferror(fd)) {
          printf("ccp4plugin) Problem reading the file.\n");
          delete [] rowdata;
          return MOLFILE_ERROR;
        }
        if (fread(rowdata, sizeof(float), extent[0], fd) != extent[0] ) {
          printf("ccp4plugin) Error reading data row.\n");
          delete [] rowdata;
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
    delete [] rowdata;
    if (ccp4->swap == 1)
      swap4_aligned(datablock, xySize * zSize);
  } else if (ccp4->voxtype == MRC_TYPE_SHORT) {
    printf("ccp4plugin) reading short (16-bit int) voxel data\n");
    short *rowdata = new short[extent[0]];
    for (coord[2] = 0; coord[2] < extent[2]; coord[2]++) {
      for (coord[1] = 0; coord[1] < extent[1]; coord[1]++) {
        if (feof(fd)) {
          printf("ccp4plugin) Unexpected end-of-file.\n");
          delete [] rowdata;
          return MOLFILE_ERROR;
        }
        if (ferror(fd)) {
          printf("ccp4plugin) Problem reading the file.\n");
          delete [] rowdata;
          return MOLFILE_ERROR;
        }
        if (fread(rowdata, sizeof(short), extent[0], fd) != extent[0] ) {
          printf("ccp4plugin) Error reading data row.\n");
          delete [] rowdata;
          return MOLFILE_ERROR;
        }
        if (ccp4->swap == 1)
          swap2_aligned(rowdata, extent[0]);
        for (coord[0] = 0; coord[0] < extent[0]; coord[0]++) {
          x = coord[ccp4->xyz2crs[0]];
          y = coord[ccp4->xyz2crs[1]];
          z = coord[ccp4->xyz2crs[2]];
          datablock[x + y*xSize + z*xySize] = rowdata[coord[0]];
        }
      }
    }
    delete [] rowdata;
  } else if (ccp4->voxtype == MRC_TYPE_SHORT2) {
    /* IMOD developers said that this is not used anymore and not worth our time to implement */
    printf("TYPE_SHORT2 not implemented yet...\n");
    return MOLFILE_ERROR;
  } else if (ccp4->voxtype == MRC_TYPE_USHORT) {
    printf("ccp4plugin) reading unsigned short (16-bit int) voxel data\n");
    unsigned short *rowdata = new unsigned short[extent[0]];
    for (coord[2] = 0; coord[2] < extent[2]; coord[2]++) {
      for (coord[1] = 0; coord[1] < extent[1]; coord[1]++) {
        if (feof(fd)) {
          printf("ccp4plugin) Unexpected end-of-file.\n");
          delete [] rowdata;
          return MOLFILE_ERROR;
        }
        if (ferror(fd)) {
          printf("ccp4plugin) Problem reading the file.\n");
          delete [] rowdata;
          return MOLFILE_ERROR;
        }
        if (fread(rowdata, sizeof(unsigned short), extent[0], fd) != extent[0] ) {
          printf("ccp4plugin) Error reading data row.\n");
          delete [] rowdata;
          return MOLFILE_ERROR;
        }
        if (ccp4->swap == 1)
          swap2_aligned(rowdata, extent[0]);
        for (coord[0] = 0; coord[0] < extent[0]; coord[0]++) {
          x = coord[ccp4->xyz2crs[0]];
          y = coord[ccp4->xyz2crs[1]];
          z = coord[ccp4->xyz2crs[2]];
          datablock[x + y*xSize + z*xySize] = rowdata[coord[0]];
        }
      }
    }
    delete [] rowdata;
  } else if (ccp4->voxtype == MRC_TYPE_UCHAR3) {
    printf("ccp4plugin) reading unsigned char * 3 (8-bit uchar * 3) voxel data\n");
    uchar3 *rowdata = new uchar3[extent[0]];
    float grayscale;
    for (coord[2] = 0; coord[2] < extent[2]; coord[2]++) {
      for (coord[1] = 0; coord[1] < extent[1]; coord[1]++) {
        if (feof(fd)) {
          printf("ccp4plugin) Unexpected end-of-file.\n");
          delete [] rowdata;
          return MOLFILE_ERROR;
        }
        if (ferror(fd)) {
          printf("ccp4plugin) Problem reading the file.\n");
          delete [] rowdata;
          return MOLFILE_ERROR;
        }
        if ( fread(rowdata, sizeof(uchar3), extent[0], fd) != extent[0] ) {
          printf("ccp4plugin) Error reading data row.\n");
          delete [] rowdata;
          return MOLFILE_ERROR;
        }
        for (coord[0] = 0; coord[0] < extent[0]; coord[0]++) {
          x = coord[ccp4->xyz2crs[0]];
          y = coord[ccp4->xyz2crs[1]];
          z = coord[ccp4->xyz2crs[2]];
          grayscale = rowdata[coord[0]].red + rowdata[coord[0]].blue + rowdata[coord[0]].green;
          datablock[x + y*xSize + z*xySize] = grayscale/3.0;
        }
      }
    }
    delete [] rowdata;
  }
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

VMDPLUGIN_API int VMDPLUGIN_init(void) { 
  memset(&plugin, 0, sizeof(molfile_plugin_t));
  plugin.abiversion = vmdplugin_ABIVERSION;
  plugin.type = MOLFILE_PLUGIN_TYPE;
  plugin.name = "ccp4";
  plugin.prettyname = "CCP4, MRC Density Map";
  plugin.author = "Eamon Caddigan, Brendan McMorrow, John Stone";
  plugin.majorv = 1;
  plugin.minorv = 7;
  plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
  plugin.filename_extension = "ccp4,mrc,map";
  plugin.open_file_read = open_ccp4_read;
  plugin.read_volumetric_metadata = read_ccp4_metadata;
  plugin.read_volumetric_data = read_ccp4_data;
  plugin.close_file_read = close_ccp4_read;
  return VMDPLUGIN_SUCCESS; 
}


VMDPLUGIN_API int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_fini(void) { return VMDPLUGIN_SUCCESS; }
