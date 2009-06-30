/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_brixplugin
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
 *      $RCSfile: brixplugin.C,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.19 $       $Date: 2009/04/29 15:45:28 $
 *
 ***************************************************************************/

/* 
 * BRIX format electron density maps.
 *
 * More info for format can be found at 
 * <http://zombie.imsb.au.dk/~mok/brix/brix-1.html#ss1.2>
 * 
 * TODO: Speed up inner loop of read_data, perhaps reading one full brick
 *       from the file at a time.
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

typedef struct {
  FILE *fd;
  int nsets;
  float prod, plus;
  molfile_volumetric_t *vol;
} brix_t;


static void *open_brix_read(const char *filepath, const char *filetype,
    int *natoms) {
  FILE *fd;
  brix_t *brix;
  char keyWord[81];
  // File header data:
  int org_x, org_y, org_z, ext_x, ext_y, ext_z;
  float grid_x, grid_y, grid_z, cell_x, cell_y, cell_z, 
        cell_alpha, cell_beta, cell_gamma, prod, plus, sigma, 
        xaxis[3], yaxis[3], zaxis[3], z1, z2, z3;
  
  fd = fopen(filepath, "rb");
  if (!fd) {
    fprintf(stderr, "brixplugin) Error opening file.\n");
    return NULL;
  }

  // read in BRIX file header information -- stored as a 512-element char array
  fscanf(fd, "%3s", keyWord);
  if (strcmp(keyWord, ":-)") != 0) {
    fprintf(stderr, "brixplugin) Error improperly formatted header.\n");
    return NULL;
  }

  // "Origin: the origin of the map in grid units, along X, Y, Z."
  fscanf(fd, " %s %d %d %d", keyWord, &org_x, &org_y, &org_z);
  if (strcasecmp(keyWord, "origin") != 0) {
    fprintf(stderr, "brixplugin) Error reading origin.\n");
    return NULL;
  }

  // "Extent: the extent (size) of the map, along X, Y, Z, in grid units"
  fscanf(fd, " %s %d %d %d", keyWord, &ext_x, &ext_y, &ext_z);
  if (strcasecmp(keyWord, "extent") != 0) {
    fprintf(stderr, "brixplugin) Error reading extent.\n");
    return NULL;
  }

  // "Grid: number of grid points along the whole unit cell, along X, Y, Z."
  fscanf(fd, " %s %f %f %f", keyWord, &grid_x, &grid_y, &grid_z);
  if (strcasecmp(keyWord, "grid") != 0) {
    fprintf(stderr, "brixplugin) Error reading grid.\n");
    return NULL;
  }

  // "Cell: Unit cell parameters"
  // cell x, y, and z are the length of the unit cell along each axis,
  // cell alpha, beta, and cell_gamma are the angles between axes.
  fscanf(fd, " %s %f %f %f %f %f %f", keyWord,
         &cell_x, &cell_y, &cell_z, &cell_alpha, &cell_beta, &cell_gamma);
  if (strcasecmp(keyWord, "cell") != 0) {
    fprintf(stderr, "brixplugin) Error reading cell.\n");
    return NULL;
  }

  // Convert angles to radians
  cell_alpha *= M_PI / 180.0;
  cell_beta *= M_PI / 180.0;
  cell_gamma *= M_PI / 180.0;

  // "Prod, plus: Constants that bring the electron density from byte to 
  // normal scale."
  fscanf(fd, " %s %f", keyWord, &prod);
  if (strcasecmp(keyWord, "prod") != 0) {
    fprintf(stderr, "brixplugin) Error reading prod.\n");
    return NULL;
  }
  fscanf(fd, " %s %f", keyWord, &plus);
  if (strcasecmp(keyWord, "plus") != 0) {
    fprintf(stderr, "brixplugin) Error reading plus.\n");
    return NULL;
  }

  // "Sigma: Rms deviation of electron density map."
  fscanf(fd, " %s %f", keyWord, &sigma);
  if (strcasecmp(keyWord, "sigma") != 0) {
    fprintf(stderr, "brixplugin) Error reading sigma.\n");
    return NULL;
  }

  // Allocate and initialize the brix structure
  brix = new brix_t;
  brix->fd = fd;
  brix->vol = NULL;
  *natoms = MOLFILE_NUMATOMS_NONE;
  brix->nsets = 1; // this file contains only one data set

  brix->prod = prod;
  brix->plus = plus;

  brix->vol = new molfile_volumetric_t[1];
  strcpy(brix->vol[0].dataname, "BRIX Electron Density Map");

  // Calculate non-orthogonal unit cell coordinates
  xaxis[0] = cell_x / grid_x;
  xaxis[1] = 0;
  xaxis[2] = 0;

  yaxis[0] = cos(cell_gamma) * cell_y / grid_y;
  yaxis[1] = sin(cell_gamma) * cell_y / grid_y;
  yaxis[2] = 0;

  z1 = cos(cell_beta);
  z2 = (cos(cell_alpha) - cos(cell_beta)*cos(cell_gamma)) / sin(cell_gamma);
  z3 = sqrt(1.0 - z1*z1 - z2*z2);
  zaxis[0] = z1 * cell_z / grid_z;
  zaxis[1] = z2 * cell_z / grid_z;
  zaxis[2] = z3 * cell_z / grid_z;

  brix->vol[0].origin[0] = xaxis[0] * org_x + yaxis[0] * org_y + 
                           zaxis[0] * org_z;
  brix->vol[0].origin[1] = yaxis[1] * org_y + zaxis[1] * org_z;
  brix->vol[0].origin[2] = zaxis[2] * org_z;

  brix->vol[0].xaxis[0] = xaxis[0] * (ext_x - 1);
  brix->vol[0].xaxis[1] = 0;
  brix->vol[0].xaxis[2] = 0;

  brix->vol[0].yaxis[0] = yaxis[0] * (ext_y - 1);
  brix->vol[0].yaxis[1] = yaxis[1] * (ext_y - 1);
  brix->vol[0].yaxis[2] = 0;

  brix->vol[0].zaxis[0] = zaxis[0] * (ext_z - 1);
  brix->vol[0].zaxis[1] = zaxis[1] * (ext_z - 1);
  brix->vol[0].zaxis[2] = zaxis[2] * (ext_z - 1);

  brix->vol[0].xsize = ext_x;
  brix->vol[0].ysize = ext_y;
  brix->vol[0].zsize = ext_z;

  brix->vol[0].has_color = 0;

  return brix;
}

static int read_brix_metadata(void *v, int *nsets, 
  molfile_volumetric_t **metadata) {
  brix_t *brix = (brix_t *)v;
  *nsets = brix->nsets; 
  *metadata = brix->vol;  

  return MOLFILE_SUCCESS;
}

static int read_brix_data(void *v, int set, float *datablock,
                         float *colorblock) {
  brix_t *brix = (brix_t *)v;
  float * cell = datablock;
  int xsize, ysize, zsize, xysize, xbrix, ybrix, zbrix, cellIndex;
  int x, y, z, xbrik, ybrik, zbrik;
  unsigned char brick[512];
  FILE * fd = brix->fd;
  float div, plus;

  // Read 512-byte "bricks" of data. Each brick contains data for 8*8*8 
  // gridpoints.
  fseek(fd, 512, SEEK_SET);

  div = 1.0 / brix->prod;
  plus = brix->plus;

  xsize = brix->vol[0].xsize;
  ysize = brix->vol[0].ysize;
  zsize = brix->vol[0].zsize;
  xysize = xsize * ysize;

  xbrix = (int) ceil((float) xsize / 8.0);
  ybrix = (int) ceil((float) ysize / 8.0);
  zbrix = (int) ceil((float) zsize / 8.0);

  // For every density value, 
  // cellIndex = (x + xbrik*8) + (y + ybrik*8)*xsize + (z + zbrik*8) * xysize
  cellIndex = 0;
  for (zbrik = 0; zbrik < zbrix; zbrik++) {
    for (ybrik = 0; ybrik < ybrix; ybrik++) {
      for (xbrik = 0; xbrik < xbrix; xbrik++) {
        if (feof(fd)) {
          fprintf(stderr, "brixplugin) Unexpected end-of-file.\n");
          return MOLFILE_ERROR;
        }
        if (ferror(fd)) {
          fprintf(stderr, "brixplugin) Error reading file.\n");
          return MOLFILE_ERROR;
        }

        // Read a data brick into the buffer.
        fread(brick, sizeof(char), 512, fd);

        // Copy the brick values into the cell
        for (z = 0; z < 8; z++) {
          for (y = 0; y < 8; y++) {
            for (x = 0; x < 8; x++) {

              if ( ((x + xbrik*8) < xsize) && ((y + ybrik*8) < ysize) &&
                   ((z + zbrik*8) < zsize) )
                cell[cellIndex] = 
                  div * ((float) brick[x+y*8+z*64] - plus);

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

static void close_brix_read(void *v) {
  brix_t *brix = (brix_t *)v;

  fclose(brix->fd);
  if (brix->vol != NULL)
    delete [] brix->vol; 
  delete brix;
}

/*
 * Initialization stuff here
 */
static molfile_plugin_t plugin;

VMDPLUGIN_EXTERN int VMDPLUGIN_init(void) { 
  memset(&plugin, 0, sizeof(molfile_plugin_t));
  plugin.abiversion = vmdplugin_ABIVERSION;
  plugin.type = MOLFILE_PLUGIN_TYPE;
  plugin.name = "brix";
  plugin.prettyname = "BRIX Density Map";
  plugin.author = "Eamon Caddigan";
  plugin.majorv = 0;
  plugin.minorv = 8;
  plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
  plugin.filename_extension = "brix,brx";
  plugin.open_file_read = open_brix_read;
  plugin.read_volumetric_metadata = read_brix_metadata;
  plugin.read_volumetric_data = read_brix_data;
  plugin.close_file_read = close_brix_read;

  return VMDPLUGIN_SUCCESS; 
}

VMDPLUGIN_EXTERN int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_EXTERN int VMDPLUGIN_fini(void) { return VMDPLUGIN_SUCCESS; }

