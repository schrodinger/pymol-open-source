/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_phiplugin
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
 *      $RCSfile: phiplugin.C,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.26 $       $Date: 2008/01/09 19:36:51 $
 *
 ***************************************************************************/

/* 
 * "Formatted ASCII '.big'" potential maps from Delphi
 *   This format is created by the 'ASCIIPHI' program which was available with
 *   Delphi V3:
 *     http://www.csb.yale.edu/userguides/datamanip/delphi/manual.html#ASCIIPHI
 *
 *   More info for this format can be found at:
 *     http://www.msg.ucsf.edu/local/programs/insightII/doc/life/insight2000.1/delphi/B_Utilities.html
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
  int ndata;
  molfile_volumetric_t *vol;
} phi_t;


// Get a string from a stream, printing any errors that occur
static char *phigets(char *s, int n, FILE *stream) {
  char *returnVal;

  if (feof(stream)) {
    fprintf(stderr, "phiplugin) Unexpected end-of-file.\n");
    returnVal = NULL; 
  }
  else if (ferror(stream)) {
    fprintf(stderr, "phiplugin) Error reading file.\n");
    return NULL;
  }
  else {
    returnVal = fgets(s, n, stream);
    if (returnVal == NULL) {
      fprintf(stderr, "phiplugin) Error reading line.\n");
    }
  }

  return returnVal;
}


static void *open_phi_read(const char *filepath, const char *filetype,
    int *natoms) {
  FILE *fd;
  phi_t *phi;
  char inbuf[LINESIZE];

  float scale, midX, midY, midZ;
  float cellSize, iGrid = 0.0;
  
  fd = fopen(filepath, "rb");
  if (!fd) {
    fprintf(stderr, "phiplugin) Error opening file.\n");
    return NULL;
  }

  /* Skip the header */
  if (phigets(inbuf, LINESIZE, fd) == NULL) {
    return NULL;
  }
  if (phigets(inbuf, LINESIZE, fd) == NULL) {
    return NULL;
  }
  if (phigets(inbuf, LINESIZE, fd) == NULL) {
    return NULL;
  }

  /* Unit cell information is located at the *end* of the file. */
  if (phigets(inbuf, LINESIZE, fd) == NULL) {
    return NULL;
  }
  while (strncasecmp(inbuf, " end of phimap", 14) != 0) {
    /* use integer division so trailing whitespace isn't included in the
     * count */
    iGrid += strlen(inbuf) / 4; 

    if (phigets(inbuf, LINESIZE, fd) == NULL) {
      return NULL;
    }
  } 

  /* Find the cube-root of the number of datapoints (this will give the
   * number of datapoints in each direction) and make sure it's an integer
   */
  cellSize = pow((double) iGrid, (double) 1.0/3.0);
  if (fabs((double)(cellSize - floor(cellSize))) > 1e-8) {
    return NULL;
  }

  /* Read the unit cell information */
  if (phigets(inbuf, LINESIZE, fd) == NULL) {
    return NULL;
  }
  sscanf(inbuf, " %f %f %f %f", &scale, &midX, &midY, &midZ);

  /* Allocate and initialize the phi structure */
  phi = new phi_t;
  phi->fd = fd;
  phi->vol = NULL;
  phi->ndata = (int) iGrid;
  *natoms = MOLFILE_NUMATOMS_NONE;
  phi->nsets = 1; /* this file contains only one data set */

  phi->vol = new molfile_volumetric_t[1];
  strcpy(phi->vol[0].dataname, "PHIMAP Electron Density Map");

  /* <midX, midY, midZ> is the middle point of the grid. */
  phi->vol[0].origin[0] = -0.5*(cellSize+1.0) / scale + midX;
  phi->vol[0].origin[1] = -0.5*(cellSize+1.0) / scale + midY;
  phi->vol[0].origin[2] = -0.5*(cellSize+1.0) / scale + midZ;

  phi->vol[0].xaxis[0] = cellSize / scale;
  phi->vol[0].xaxis[1] = 0;
  phi->vol[0].xaxis[2] = 0;

  phi->vol[0].yaxis[0] = 0;
  phi->vol[0].yaxis[1] = cellSize / scale;
  phi->vol[0].yaxis[2] = 0;
  
  phi->vol[0].zaxis[0] = 0;
  phi->vol[0].zaxis[1] = 0;
  phi->vol[0].zaxis[2] = cellSize / scale;

  phi->vol[0].xsize = (int) cellSize;
  phi->vol[0].ysize = (int) cellSize;
  phi->vol[0].zsize = (int) cellSize;

  phi->vol[0].has_color = 0;

  return phi;
}

static int read_phi_metadata(void *v, int *nsets, 
  molfile_volumetric_t **metadata) {
  phi_t *phi = (phi_t *)v;
  *nsets = phi->nsets; 
  *metadata = phi->vol;  

  return MOLFILE_SUCCESS;
}

static int read_phi_data(void *v, int set, float *datablock,
                         float *colorblock) {
  phi_t *phi = (phi_t *)v;
  float *cellIndex;
  int value, ndata, count = 0;
  FILE *fd = phi->fd;
  char inbuf[LINESIZE], currNum[5], *currChar;

  cellIndex = datablock;
  ndata = phi->ndata;
  memset(currNum, 0, 5);

  /* Skip the header */
  rewind(fd);
  if (phigets(inbuf, LINESIZE, fd) == NULL) {
    return MOLFILE_ERROR;
  }
  if (phigets(inbuf, LINESIZE, fd) == NULL) {
    return MOLFILE_ERROR;
  }
  if (phigets(inbuf, LINESIZE, fd) == NULL) {
    return MOLFILE_ERROR;
  }

  /* Read the densities. Order for file is x fast, y medium, z slow */
  while (count < ndata) {
    if (phigets(inbuf, LINESIZE, fd) == NULL) {
      return MOLFILE_ERROR;
    }

    for (currChar = inbuf; (*currChar != '\n') && (*currChar != '\0'); 
         currChar += 4) {
      strncpy(currNum, currChar, 4);
      value = atoi(currNum);
      /* Scale -- units are in kT/e (25.6mV, 0.593 kcal/mole at 25°C). */
      *cellIndex = 0.01 * (value - 5000);
      cellIndex++;
      count++;
    }
  }

  return MOLFILE_SUCCESS;
}

static void close_phi_read(void *v) {
  phi_t *phi = (phi_t *)v;

  fclose(phi->fd);
  if (phi->vol != NULL)
    delete [] phi->vol; 
  delete phi;
}

/*
 * Initialization stuff here
 */
static molfile_plugin_t plugin;

VMDPLUGIN_EXTERN int VMDPLUGIN_init(void) {
  memset(&plugin, 0, sizeof(molfile_plugin_t));
  plugin.abiversion = vmdplugin_ABIVERSION;
  plugin.type = MOLFILE_PLUGIN_TYPE;
  plugin.name = "delphibig";
  plugin.prettyname = "Delphi 'Big' Formatted Potential Map";
  plugin.author = "Eamon Caddigan";
  plugin.majorv = 0;
  plugin.minorv = 7;
  plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
  plugin.filename_extension = "big";
  plugin.open_file_read = open_phi_read;
  plugin.read_volumetric_metadata = read_phi_metadata;
  plugin.read_volumetric_data = read_phi_data;
  plugin.close_file_read = close_phi_read;
  return VMDPLUGIN_SUCCESS; 
}

VMDPLUGIN_EXTERN int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)(void *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_EXTERN int VMDPLUGIN_fini(void) { return VMDPLUGIN_SUCCESS; }

