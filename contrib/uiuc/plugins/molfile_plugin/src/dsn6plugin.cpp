/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_dsn6plugin
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
 *      $RCSfile: dsn6plugin.C,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.24 $       $Date: 2009/04/29 15:45:29 $
 *
 ***************************************************************************/

/* 
 * DSN6 format electron density maps.
 *
 * More info for format can be found at 
 * <http://www.uoxray.uoregon.edu/tnt/manual/node104.html>
 * TODO: Check byte-swapping and alignment issues in read_data()
 *
 * apparently there are some gotchas there, and mapman does some things to 
 * fix up some variants of DSN6:
 * Gerard "DVD" Kleywegt  gerard@xray.bmc.uu.se: 
 *   "dale's description is largely correct, but in elements 13 to 15
 *    the actual cell angles are stored (multiplied by a scale factor)
 *    rather than their cosines. also, turbo-frodo-style dsn6 maps
 *    differ slightly from o-style dsn6 maps. shown below is a chunk
 *    of code (yes, fortran) from mapman that fills the header record"
 * See original email for mapman code snippet:
 *      http://o-info.bioxray.dk/pipermail/o-info/2002-June/005993.html
 * 
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

#if defined(WIN32) || defined(WIN64)
#define strcasecmp stricmp
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "molfile_plugin.h"
#include "endianswap.h"

typedef struct {
  FILE *fd;
  int nsets;
  float prod, plus;
  molfile_volumetric_t *vol;
} dsn6_t;


static void *open_dsn6_read(const char *filepath, const char *filetype,
    int *natoms) {
  FILE *fd;
  dsn6_t *dsn6;
  short fileHeader[19];
  int start_x, start_y, start_z, extent_x, extent_y, extent_z;
  float scale, A, B, C, alpha, beta, gamma, 
        xaxis[3], yaxis[3], zaxis[3], z1, z2, z3;
  
  fd = fopen(filepath, "rb");
  if (!fd) {
    fprintf(stderr, "Error opening file.\n");
    return NULL;
  }

  // Read the header into a 19-element int array. The integers are stored
  // according to the endianness of the machine used to write the file, so
  // swapping may be necessary.
  fread(fileHeader, sizeof(short), 19, fd);

  // Check byte-order, swapping if necessary
  if (fileHeader[18] == 25600)
    swap2_aligned(fileHeader, 19);
  else if (fileHeader[18] != 100) {
    fprintf(stderr, "Error reading file header.\n");
    return NULL;
  }
  // else fileHeader[18] is 100, byte-order is fine 
  // (this value is hard-coded into the file format)

  // Unit cell origin, in grid coordinates
  start_x = fileHeader[0];
  start_y = fileHeader[1];
  start_z = fileHeader[2];

  // Unit cell size
  extent_x = fileHeader[3];
  extent_y = fileHeader[4];
  extent_z = fileHeader[5];

  // Constant cell scaling factor
  scale = 1.0 / fileHeader[17];

  // Unit cell edge, in angstroms / sample
  A = scale * fileHeader[9] / fileHeader[6];
  B = scale * fileHeader[10] / fileHeader[7];
  C = scale * fileHeader[11] / fileHeader[8];

  // Cell angles, in radians
  alpha = scale * fileHeader[12] * M_PI / 180.0;
  beta = scale * fileHeader[13] * M_PI / 180.0;
  gamma = scale * fileHeader[14] * M_PI / 180.0;

  // Allocate and initialize the dsn6 structure
  dsn6 = new dsn6_t;
  dsn6->fd = fd;
  dsn6->vol = NULL;
  *natoms = MOLFILE_NUMATOMS_NONE;
  dsn6->nsets = 1; // this file contains only one data set

  // Grid data scaling constants
  dsn6->prod = (float) fileHeader[15] / fileHeader[18];
  dsn6->plus = fileHeader[16];

  dsn6->vol = new molfile_volumetric_t[1];
  strcpy(dsn6->vol[0].dataname, "DSN6 Electron Density Map");

  // Calculate non-orthogonal unit cell coordinates
  xaxis[0] = A;
  xaxis[1] = 0;
  xaxis[2] = 0;

  yaxis[0] = cos(gamma) * B;
  yaxis[1] = sin(gamma) * B;
  yaxis[2] = 0;

  z1 = cos(beta);
  z2 = (cos(alpha) - cos(beta)*cos(gamma)) / sin(gamma);
  z3 = sqrt(1.0 - z1*z1 - z2*z2);
  zaxis[0] = z1 * C;
  zaxis[1] = z2 * C;
  zaxis[2] = z3 * C;

  // Convert the origin from grid space to cartesian coordinates
  dsn6->vol[0].origin[0] = xaxis[0] * start_x + yaxis[0] * start_y + 
                           zaxis[0] * start_z;
  dsn6->vol[0].origin[1] = yaxis[1] * start_y + zaxis[1] * start_z;
  dsn6->vol[0].origin[2] = zaxis[2] * start_z;

  dsn6->vol[0].xaxis[0] = xaxis[0] * (extent_x-1);
  dsn6->vol[0].xaxis[1] = 0;
  dsn6->vol[0].xaxis[2] = 0;

  dsn6->vol[0].yaxis[0] = yaxis[0] * (extent_y-1);
  dsn6->vol[0].yaxis[1] = yaxis[1] * (extent_y-1);
  dsn6->vol[0].yaxis[2] = 0;

  dsn6->vol[0].zaxis[0] = zaxis[0] * (extent_z-1);
  dsn6->vol[0].zaxis[1] = zaxis[1] * (extent_z-1);
  dsn6->vol[0].zaxis[2] = zaxis[2] * (extent_z-1);

  dsn6->vol[0].xsize = extent_x;
  dsn6->vol[0].ysize = extent_y;
  dsn6->vol[0].zsize = extent_z;

  dsn6->vol[0].has_color = 0;

  return dsn6;
}

static int read_dsn6_metadata(void *v, int *nsets, 
  molfile_volumetric_t **metadata) {
  dsn6_t *dsn6 = (dsn6_t *)v;
  *nsets = dsn6->nsets; 
  *metadata = dsn6->vol;  

  return MOLFILE_SUCCESS;
}

static int read_dsn6_data(void *v, int set, float *datablock,
                         float *colorblock) {
  dsn6_t *dsn6 = (dsn6_t *)v;
  float * cell = datablock;
  unsigned char brick[512];
  unsigned char* brickPtr = NULL;
  int xsize, ysize, zsize, xysize, xbrix, ybrix, zbrix, cellIndex;
  int x, y, z, xbrik, ybrik, zbrik;
  FILE * fd = dsn6->fd;
  float div, plus;

  // Read 512-byte "bricks" of data. Each brick contains data for 8*8*8 
  // gridpoints.
  fseek(fd, 512, SEEK_SET);

  div = 1.0 / dsn6->prod;
  plus = dsn6->plus;

  xsize = dsn6->vol[0].xsize;
  ysize = dsn6->vol[0].ysize;
  zsize = dsn6->vol[0].zsize;
  xysize = xsize * ysize;

  xbrix = (int) ceil((float) xsize / 8.0);
  ybrix = (int) ceil((float) ysize / 8.0);
  zbrix = (int) ceil((float) zsize / 8.0);

  cellIndex = 0;
  for (zbrik = 0; zbrik < zbrix; zbrik++) {
    for (ybrik = 0; ybrik < ybrix; ybrik++) {
      for (xbrik = 0; xbrik < xbrix; xbrik++) {
        // Read the next brick into the buffer and swap its bytes.
        if (feof(fd)) {
          fprintf(stderr, "Unexpected end-of-file.\n");
          return MOLFILE_ERROR;
        }
        if (ferror(fd)) {
          fprintf(stderr, "Error reading file.\n");
          return MOLFILE_ERROR;
        }

        fread(brick, sizeof(char), 512, fd);
        swap2_unaligned(brick, 512*sizeof(char));
        brickPtr = brick;

        for (z = 0; z < 8; z++) {
          if ((z + zbrik*8) >= zsize) {
            cellIndex += (8 - z) * xysize;
            break;
          }

          for (y = 0; y < 8; y++) {
            if ((y + ybrik*8) >= ysize) {
              cellIndex += (8 - y) * xsize;
              brickPtr += (8 - y) * 8;
              break;
            }

            for (x = 0; x < 8; x++) {
              if ((x + xbrik*8) >= xsize) {
                cellIndex += 8 - x;
                brickPtr += 8 - x;
                break;
              }

              // cell[(x+xbrik*8) + (y+ybrik*8)*xsize + (z+zbrik*8)*xysize] =
              // div * ((float) brick[x + 8*y + 64*z] - plus)
              cell[cellIndex] = div * ((float) *brickPtr - plus);

              brickPtr++;
              cellIndex++;
            } // end for(x)
           
            cellIndex += xsize - 8;
          } // end for(y)
         
          cellIndex += xysize - 8*xsize;
        } // end for(z)
      
        cellIndex += 8 - 8*xysize; 
      } // end for(xbrik)

      cellIndex += 8 * (xsize - xbrix);
    } // end for(ybrik)

    cellIndex += 8 * (xysize - xsize*ybrik);
  } // end for(zbrik)
 
  return MOLFILE_SUCCESS;
}

static void close_dsn6_read(void *v) {
  dsn6_t *dsn6 = (dsn6_t *)v;

  fclose(dsn6->fd);
  if (dsn6->vol != NULL)
    delete [] dsn6->vol; 
  delete dsn6;
}

/*
 * Initialization stuff here
 */
static molfile_plugin_t plugin;

VMDPLUGIN_EXTERN int VMDPLUGIN_init(void) { 
  memset(&plugin, 0, sizeof(molfile_plugin_t));
  plugin.abiversion = vmdplugin_ABIVERSION;
  plugin.type = MOLFILE_PLUGIN_TYPE;
  plugin.name = "DSN6";
  plugin.prettyname = "DSN6";
  plugin.author = "Eamon Caddigan";
  plugin.majorv = 0;
  plugin.minorv = 6;
  plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
  plugin.filename_extension = "ds6,dsn6,omap";
  plugin.open_file_read = open_dsn6_read;
  plugin.read_volumetric_metadata = read_dsn6_metadata;
  plugin.read_volumetric_data = read_dsn6_data;
  plugin.close_file_read = close_dsn6_read;
  return VMDPLUGIN_SUCCESS; 
}

VMDPLUGIN_EXTERN int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_EXTERN int VMDPLUGIN_fini(void) { return VMDPLUGIN_SUCCESS; }

