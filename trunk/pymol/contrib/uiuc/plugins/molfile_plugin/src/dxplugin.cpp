/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_dxplugin
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
 *      $RCSfile: dxplugin.C,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.21 $       $Date: 2006/02/23 19:36:44 $
 *
 ***************************************************************************/

/* 
 * DX potential maps
 *
 * Format of the file is:
 * # Comments
 * .
 * .
 * .
 * object 1 class gridpositions counts xn yn zn
 * origin xorg yorg zorg
 * delta xdel 0 0
 * delta 0 ydel 0
 * delta 0 0 zdel
 * object 2 class gridconnections counts xn yn zn
 * object 3 class array type double rank 0 items [ xn*yn*zn ] data follows
 * f1 f2 f3
 * f4 f5 f6
 * .
 * .
 * .
 * object "Dataset name" class field
 
 * Where xn, yn, and zn are the number of data points along each axis;
 * xorg, yorg, and zorg is the origin of the grid, in angstroms;
 * xdel, ydel, and zdel are the scaling factors to convert grid units to
 * angstroms.
 * Grid data follows, with three values per line, ordered z fast, y medium,
 * and x slow.
 *
 * Note: this plugin was written to read the OpenDX files created by the
 * APBS program, and may not support all of the features in the file format.
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
#define strcasecmp  stricmp
#define strncasecmp strnicmp
#endif

#include "molfile_plugin.h"

#define LINESIZE 85

typedef struct {
  FILE *fd;
  int nsets;
  molfile_volumetric_t *vol;
} dx_t;


// Get a string from a stream, printing any errors that occur
static char *dxgets(char *s, int n, FILE *stream) {
  char *returnVal;

  if (feof(stream)) {
    fprintf(stderr, "Unexpected end-of-file.\n");
    return NULL;
  } else if (ferror(stream)) {
    fprintf(stderr, "Error reading file.\n");
    return NULL;
  } else {
    returnVal = fgets(s, n, stream);
    if (returnVal == NULL) {
      fprintf(stderr, "Error reading line.\n");
    }
  }

  return returnVal;
}


static void *open_dx_read(const char *filepath, const char *filetype,
    int *natoms) {
  FILE *fd;
  dx_t *dx;
  char inbuf[LINESIZE];
  int xsize, ysize, zsize;
  float orig[3], xdelta[3], ydelta[3], zdelta[3];
  
  fd = fopen(filepath, "rb");
  if (!fd) {
    fprintf(stderr, "Error opening file.\n");
    return NULL;
  }

  /* skip comments */
  do {
    if (dxgets(inbuf, LINESIZE, fd) == NULL) 
      return NULL;
  } while (inbuf[0] == '#');

  /* get the number of grid points */
  if (sscanf(inbuf, "object 1 class gridpositions counts %d %d %d", &xsize, &ysize, &zsize) != 3) {
    fprintf(stderr, "Error reading grid dimensions.\n");
    return NULL;
  }

  /* get the cell origin */
  if (dxgets(inbuf, LINESIZE, fd) == NULL) {
    return NULL;
  }
  if (sscanf(inbuf, "origin %e %e %e", orig, orig+1, orig+2) != 3) {
    fprintf(stderr, "Error reading grid origin.\n");
    return NULL;
  }

  /* get the cell dimensions */
  if (dxgets(inbuf, LINESIZE, fd) == NULL) {
    return NULL;
  }
  if (sscanf(inbuf, "delta %e %e %e", xdelta, xdelta+1, xdelta+2) != 3) {
    fprintf(stderr, "Error reading cell x-dimension.\n");
    return NULL;
  }

  if (dxgets(inbuf, LINESIZE, fd) == NULL) {
    return NULL;
  }
  if (sscanf(inbuf, "delta %e %e %e", ydelta, ydelta+1, ydelta+2) != 3) {
    fprintf(stderr, "Error reading cell y-dimension.\n");
    return NULL;
  }

  if (dxgets(inbuf, LINESIZE, fd) == NULL) {
    return NULL;
  }
  if (sscanf(inbuf, "delta %e %e %e", zdelta, zdelta+1, zdelta+2) != 3) {
    fprintf(stderr, "Error reading cell z-dimension.\n");
    return NULL;
  }

  /* skip the last two lines of the header (described at the beginning of
   * the code), which aren't utilized by APBS-produced DX files.  */
  if (dxgets(inbuf, LINESIZE, fd) == NULL) 
    return NULL;
  if (dxgets(inbuf, LINESIZE, fd) == NULL)
    return NULL;
 
  /* allocate and initialize the dx structure */
  dx = new dx_t;
  dx->fd = fd;
  dx->vol = NULL;
  *natoms = MOLFILE_NUMATOMS_NONE;
  dx->nsets = 1; /* this file contains only one data set */

  dx->vol = new molfile_volumetric_t[1];
  strcpy(dx->vol[0].dataname, "DX map");

  /* Set the unit cell origin and basis vectors */
  for (int i=0; i<3; i++) {
    dx->vol[0].origin[i] = orig[i];
    dx->vol[0].xaxis[i] = xdelta[i] * (xsize-1);
    dx->vol[0].yaxis[i] = ydelta[i] * (ysize-1);
    dx->vol[0].zaxis[i] = zdelta[i] * (zsize-1);
  }

  dx->vol[0].xsize = xsize;
  dx->vol[0].ysize = ysize;
  dx->vol[0].zsize = zsize;

  /* DX files contain no color information. Taken from edmplugin.C */
  dx->vol[0].has_color = 0;

  return dx;
}

static int read_dx_metadata(void *v, int *nsets, 
  molfile_volumetric_t **metadata) {
  dx_t *dx = (dx_t *)v;
  *nsets = dx->nsets; 
  *metadata = dx->vol;  

  return MOLFILE_SUCCESS;
}

static int read_dx_data(void *v, int set, float *datablock,
                         float *colorblock) {
  dx_t *dx = (dx_t *)v;
  FILE *fd = dx->fd;
  char inbuf[LINESIZE];
  float grid[3];
  int x, y, z, xsize, ysize, zsize, xysize, count, total, i;

  xsize = dx->vol[0].xsize;
  ysize = dx->vol[0].ysize;
  zsize = dx->vol[0].zsize;
  xysize = xsize * ysize;
  total = xysize * zsize;

  /* Read the values from the file */
  x = y = z = 0;
  for (count = 0; count < total/3; count++) {
    if (dxgets(inbuf, LINESIZE, fd) == NULL )
      return MOLFILE_ERROR;
    if (sscanf(inbuf, "%e %e %e", &grid[0], &grid[1], &grid[2]) != 3) {
      fprintf(stderr, "Error reading grid data.\n");
      return MOLFILE_ERROR;
    }

    for (i = 0; i < 3; i++) { 
      datablock[x + y*xsize + z*xysize] = grid[i];
      z++;
      if (z >= zsize) {
        z = 0;
        y++;
        if (y >= ysize) {
          y = 0;
          x++;
        }
      }
    }
  }

  if ((total%3) != 0) {
    if (dxgets(inbuf, LINESIZE, fd) == NULL )
      return MOLFILE_ERROR;

    count = sscanf(inbuf, "%e %e %e", &grid[0], &grid[1], &grid[2]);
    if (count != (total%3)) {
      fprintf(stderr, "Error: incorrect number of data points.\n");
      return MOLFILE_ERROR;
    }

    for (i = 0; i < count; i++) {
      datablock[x + y*xsize + z*xysize] = grid[i];
      z++;
    }
  }
  
  char dxname[256];
  while (dxgets(inbuf, LINESIZE, dx->fd)) {
    if (sscanf(inbuf, "object \"%[^\"]\" class field", dxname) == 1) {
      // a dataset name has been found; override the default
      strcpy(dx->vol[0].dataname, dxname);
      break;
    }
  }

  return MOLFILE_SUCCESS;
}

static void close_dx_read(void *v) {
  dx_t *dx = (dx_t *)v;
  
  fclose(dx->fd);
  if (dx->vol != NULL)
    delete [] dx->vol; 
  delete dx;
}

/*
 * Initialization stuff here
 */
static molfile_plugin_t plugin = {
  vmdplugin_ABIVERSION,   /* ABI version */
  MOLFILE_PLUGIN_TYPE, 	  /* plugin type */
  "dx",                   /* short file format description */
  "DX",                   /* pretty file format description */
  "Eamon Caddigan",       /* author(s) */
  0,                      /* major version */
  6,                      /* minor version */
  VMDPLUGIN_THREADSAFE,   /* is reentrant */
  "dx"                    /* filename extension */
};

VMDPLUGIN_EXTERN int VMDPLUGIN_init(void) { return VMDPLUGIN_SUCCESS; }
VMDPLUGIN_EXTERN int VMDPLUGIN_fini(void) { return VMDPLUGIN_SUCCESS; }
VMDPLUGIN_EXTERN int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  plugin.open_file_read = open_dx_read;
  plugin.read_volumetric_metadata = read_dx_metadata;
  plugin.read_volumetric_data = read_dx_data;
  plugin.close_file_read = close_dx_read;
  (*cb)(v, (vmdplugin_t *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

