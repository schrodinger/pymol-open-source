/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_mapplugin
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
 *      $RCSfile: mapplugin.C,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.14 $       $Date: 2009/04/29 15:45:31 $
 *
 ***************************************************************************/

/* 
 * Autodock Grid Map File format plugin
 *
 * More info for this format can be found at
 * <http://www.scripps.edu/pub/olson-web/gmm/autodock/ad305/
 *  Using_AutoDock_305.21.html#pgfId=75765>
 * 
 */

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <string.h>

#include "molfile_plugin.h"

#define LINESIZE 85

typedef struct {
  FILE *fd;
  int nsets;
  molfile_volumetric_t *vol;
} gridmap_t;


// Get a string from a stream, printing any errors that occur
static char *mapgets(char *s, int n, FILE *stream) {
  char *returnVal;

  if (feof(stream)) {
    fprintf(stderr, "mapplugin) Unexpected end-of-file.\n");
    returnVal = NULL;
  }
  else if (ferror(stream)) {
    fprintf(stderr, "mapplugin) Error reading file.\n");
    return NULL;
  }
  else {
    returnVal = fgets(s, n, stream);
    if (returnVal == NULL) {
      fprintf(stderr, "mapplugin) Error reading line.\n");
    }
  }

  return returnVal;
}


static void *open_map_read(const char *filepath, const char *filetype,
    int *natoms) {
  FILE *fd;
  gridmap_t *map;
  char inbuf[LINESIZE];

  float spacing, midX, midY, midZ;
  int xsize, ysize, zsize;
  
  fd = fopen(filepath, "rb");
  if (!fd) {
    fprintf(stderr, "mapplugin) Error opening file.\n");
    return NULL;
  }

  /* Skip the header */
  if (mapgets(inbuf, LINESIZE, fd) == NULL) 
    return NULL;
  if (mapgets(inbuf, LINESIZE, fd) == NULL) 
    return NULL;
  if (mapgets(inbuf, LINESIZE, fd) == NULL) 
    return NULL;

  /* Space between grid points */
  if (mapgets(inbuf, LINESIZE, fd) == NULL) 
    return NULL;
  if (sscanf(inbuf, "SPACING %f", &spacing) != 1)
    return NULL;

  /* Grid size in grid units */
  if (mapgets(inbuf, LINESIZE, fd) == NULL) 
    return NULL;
  if (sscanf(inbuf, "NELEMENTS %d %d %d", &xsize, &ysize, &zsize) != 3) {
    fprintf(stderr, "mapplugin) Cannot read NELEMENTS.\n");
    return NULL;
  }

  /* XXX - I don't know why this is necessary */
  xsize++;
  ysize++;
  zsize++;

  /* Center of the cell */
  if (mapgets(inbuf, LINESIZE, fd) == NULL) 
    return NULL;
  if (sscanf(inbuf, "CENTER %f %f %f", &midX, &midY, &midZ) != 3)
    return NULL;

  /* Allocate and initialize the map structure */
  map = new gridmap_t;
  map->fd = fd;
  map->vol = NULL;
  *natoms = MOLFILE_NUMATOMS_NONE;
  map->nsets = 1; /* this file contains only one data set */

  map->vol = new molfile_volumetric_t[1];
  strcpy(map->vol[0].dataname, "Grid Map File");

  /* <midX, midY, midZ> is the middle point of the grid. */
  map->vol[0].origin[0] = -0.5*(xsize+1.0)* spacing  + midX;
  map->vol[0].origin[1] = -0.5*(ysize+1.0)* spacing  + midY;
  map->vol[0].origin[2] = -0.5*(zsize+1.0)* spacing  + midZ;

  map->vol[0].xaxis[0] = xsize * spacing;
  map->vol[0].xaxis[1] = 0;
  map->vol[0].xaxis[2] = 0;

  map->vol[0].yaxis[0] = 0;
  map->vol[0].yaxis[1] = ysize * spacing;
  map->vol[0].yaxis[2] = 0;
  
  map->vol[0].zaxis[0] = 0;
  map->vol[0].zaxis[1] = 0;
  map->vol[0].zaxis[2] = zsize * spacing;

  map->vol[0].xsize = xsize;
  map->vol[0].ysize = ysize;
  map->vol[0].zsize = zsize;

  map->vol[0].has_color = 0;

  return map;
}

static int read_map_metadata(void *v, int *nsets, 
  molfile_volumetric_t **metadata) {
  gridmap_t *map = (gridmap_t *)v;
  *nsets = map->nsets; 
  *metadata = map->vol;  

  return MOLFILE_SUCCESS;
}

static int read_map_data(void *v, int set, float *datablock,
                         float *colorblock) {
  gridmap_t *map = (gridmap_t *)v;
  FILE *fd = map->fd;
  float *cellIndex;
  char inbuf[LINESIZE];
  int count, ndata;

  cellIndex = datablock;
  count = 0;
  ndata = map->vol[0].xsize * map->vol[0].ysize * map->vol[0].zsize;

  /* Read the densities. Order for file is x fast, y medium, z slow */
  while (count < ndata) {
    if (mapgets(inbuf, LINESIZE, fd) == NULL) {
      return MOLFILE_ERROR;
    }

    *cellIndex = atof(inbuf);

    cellIndex++;
    count++;
  }

  return MOLFILE_SUCCESS;
}

static void close_map_read(void *v) {
  gridmap_t *map = (gridmap_t *)v;

  fclose(map->fd);
  if (map->vol != NULL)
    delete [] map->vol; 
  delete map;
}

/*
 * Initialization stuff here
 */
static molfile_plugin_t plugin;

VMDPLUGIN_EXTERN int VMDPLUGIN_init(void) {
  memset(&plugin, 0, sizeof(molfile_plugin_t));
  plugin.abiversion = vmdplugin_ABIVERSION;
  plugin.type = MOLFILE_PLUGIN_TYPE;
  plugin.name = "map";
  plugin.prettyname = "Autodock Grid Map";
  plugin.author = "Eamon Caddigan";
  plugin.majorv = 0;
  plugin.minorv = 6;
  plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
  plugin.filename_extension = "map";
  plugin.open_file_read = open_map_read;
  plugin.read_volumetric_metadata = read_map_metadata;
  plugin.read_volumetric_data = read_map_data;
  plugin.close_file_read = close_map_read;
  return VMDPLUGIN_SUCCESS; 
}

VMDPLUGIN_EXTERN int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_EXTERN int VMDPLUGIN_fini(void) { return VMDPLUGIN_SUCCESS; }

