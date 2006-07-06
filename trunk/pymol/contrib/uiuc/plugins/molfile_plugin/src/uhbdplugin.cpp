/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_uhbdplugin
#define STATIC_PLUGIN 1

/***************************************************************************
 *cr
 *cr            (C) Copyright 1995-2005 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: uhbdplugin.C,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.8 $       $Date: 2006/03/14 08:16:28 $
 *
 ***************************************************************************/

/*
 * Plugin by Alexander Spaar for reading UHBD grid files
 * UHBD related docs:
 *   http://adrik.bchs.uh.edu/uhbd.html
 *   http://adrik.bchs.uh.edu/uhbd/
 *   http://mccammon.ucsd.edu/uhbd.html
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
} uhbd_t;


// Get a string from a stream, printing any errors that occur
static char *uhbdgets(char *s, int n, FILE *stream, const char *msg) {
  char *returnVal;

  if (feof(stream)) {
    printf(msg);
    printf("uhbdplugin) Unexpected end-of-file.\n");
    return NULL;
  } else if (ferror(stream)) {
    printf(msg);
    printf("uhbdplugin) Error reading file.\n");
    return NULL;
  } else {
    returnVal = fgets(s, n, stream);
    if (returnVal == NULL) {
      printf(msg);
      printf("uhbdplugin) Encountered EOF or error reading line.\n");
    }
  }

  return returnVal;
}


static void *open_uhbd_read(const char *filepath, const char *filetype,
    int *natoms) {
  FILE *fd;
  uhbd_t *uhbd;
  char inbuf[LINESIZE];
  int xsize, ysize, zsize;
  float orig[3], ra, o[3], s[3]; //xdelta[3], ydelta[3], zdelta[3];
  
  if ((fd = fopen(filepath, "rb")) == NULL) {
    printf("uhbdplugin) Error opening file.\n");
    return NULL;
  }

  // Read the header
  if (uhbdgets(inbuf, LINESIZE, fd, 
      "uhbdplugin) error while skipping header lines\n") == NULL) 
    return NULL;
  if (uhbdgets(inbuf, LINESIZE, fd,
      "uhbdplugin) error while skipping header lines\n") == NULL) 
    return NULL;

  /* get grid dimensions, spacing and origin */
  if (uhbdgets(inbuf, LINESIZE, fd,
      "uhbdplugin) error while getting grid dimensions\n") == NULL) {
    return NULL;
  }
  if (sscanf(inbuf, "%d %d %d %e %e %e %e", &xsize, &ysize, &zsize, &ra, orig, orig+1, orig+2)  != 7) {
    printf("uhbdplugin) Error reading grid dimensions, spacing and origin.\n");
    return NULL;
  }
  if (uhbdgets(inbuf, LINESIZE, fd,
      "uhbdplugin) error while skipping header lines\n") == NULL) 
    return NULL;
  if (uhbdgets(inbuf, LINESIZE, fd,
      "uhbdplugin) error while skipping header lines\n") == NULL) 
    return NULL;

  /* allocate and initialize the uhbd structure */
  uhbd = new uhbd_t;
  uhbd->fd = fd;
  uhbd->vol = NULL;
  *natoms = MOLFILE_NUMATOMS_NONE;
  uhbd->nsets = 1; /* this file contains only one data set */

  uhbd->vol = new molfile_volumetric_t[1];
  strcpy(uhbd->vol[0].dataname, "UHBD ascii Electron Density Map");

  /* Set the unit cell origin and basis vectors */
  for (int i=0; i<3; i++) {
    uhbd->vol[0].origin[i] = orig[i] + ra;
    o[i] = uhbd->vol[0].origin[i];
  }

  uhbd->vol[0].xaxis[0] = ra * (xsize-1);
  uhbd->vol[0].xaxis[1] = 0;
  uhbd->vol[0].xaxis[2] = 0;

  uhbd->vol[0].yaxis[0] = 0;
  uhbd->vol[0].yaxis[1] = ra * (ysize-1);
  uhbd->vol[0].yaxis[2] = 0;

  uhbd->vol[0].zaxis[0] = 0;
  uhbd->vol[0].zaxis[1] = 0;
  uhbd->vol[0].zaxis[2] = ra * (zsize-1);

  s[0] = uhbd->vol[0].xaxis[0];
  s[1] = uhbd->vol[0].yaxis[1];
  s[2] = uhbd->vol[0].zaxis[2];

  uhbd->vol[0].xsize = xsize;
  uhbd->vol[0].ysize = ysize;
  uhbd->vol[0].zsize = zsize;

  /* UHBD files contain no color information. */
  uhbd->vol[0].has_color = 0;

  return uhbd;
}

static int read_uhbd_metadata(void *v, int *nsets, 
  molfile_volumetric_t **metadata) {
  uhbd_t *uhbd = (uhbd_t *)v;
  *nsets = uhbd->nsets; 
  *metadata = uhbd->vol;  

  return MOLFILE_SUCCESS;
}

static int read_uhbd_data(void *v, int set, float *datablock,
                         float *colorblock) {
  uhbd_t *uhbd = (uhbd_t *)v;
  FILE *fd = uhbd->fd;
  char inbuf[LINESIZE];
  float grid[6];
  int z, xsize, ysize, zsize, xysize, count, count2, total, i;

  xsize = uhbd->vol[0].xsize;
  ysize = uhbd->vol[0].ysize;
  zsize = uhbd->vol[0].zsize;
  xysize = xsize * ysize;
  total = xysize * zsize;

  /* Read the values from the file */
  for (z = 0; z < zsize; z++) {
    // read header
    if (uhbdgets(inbuf, LINESIZE, fd, 
        "uhbdplugin) error while getting density plane indices\n") == NULL)
      return MOLFILE_ERROR;

    // read data
    for (count = 0; count < xysize/6; count++) {
      if (uhbdgets(inbuf, LINESIZE, fd,
          "uhbdplugin) error while getting density values\n") == NULL)
        return MOLFILE_ERROR;

      if (sscanf(inbuf, "%e %e %e %e %e %e", &grid[0], &grid[1], &grid[2], &grid[3], &grid[4], &grid[5]) != 6) {
        printf("uhbdplugin) Error reading grid data.\n");
        return MOLFILE_ERROR;
      }
    
      for (i = 0; i < 6; i++) { 
        datablock[i + count*6 + z*xysize] = grid[i];
      }
    }

    if ((xysize%6) != 0) {
      if (uhbdgets(inbuf, LINESIZE, fd, 
          "uhbdplugin) error reading data elements modulo 6\n") == NULL)
        return MOLFILE_ERROR;

      count2 = sscanf(inbuf, "%e %e %e %e %e %e", &grid[0], &grid[1], &grid[2], &grid[3], &grid[4], &grid[5]);
      if (count2 != (xysize%6)) {
        printf("uhbdplugin) Error: incorrect number of data points.\n");
        return MOLFILE_ERROR;
      }

      for (i = 0; i < count2; i++) {
        datablock[i + (count+1)*6 + z*xysize] = grid[i];
      }
    }
  }

  return MOLFILE_SUCCESS;
}


static void close_uhbd_read(void *v) {
  uhbd_t *uhbd = (uhbd_t *)v;
  
  fclose(uhbd->fd);
  if (uhbd->vol != NULL)
    delete [] uhbd->vol; 
  delete uhbd;
}


/*
 * Initialization stuff here
 */
static molfile_plugin_t plugin = {
  vmdplugin_ABIVERSION,   /* ABI version */
  MOLFILE_PLUGIN_TYPE,    /* plugin type */
  "uhbd",                 /* short file format description */
  "UHBD Grid",            /* pretty file format description */
  "Alexander Spaar",      /* author(s) */
  0,                      /* major version */
  2,                      /* minor version */
  VMDPLUGIN_THREADSAFE,   /* is reentrant */
  "uhbdgrd,grd"           /* filename extension */
};


VMDPLUGIN_EXTERN int VMDPLUGIN_init(void) { return VMDPLUGIN_SUCCESS; }
VMDPLUGIN_EXTERN int VMDPLUGIN_fini(void) { return VMDPLUGIN_SUCCESS; }
VMDPLUGIN_EXTERN int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  plugin.open_file_read = open_uhbd_read;
  plugin.read_volumetric_metadata = read_uhbd_metadata;
  plugin.read_volumetric_data = read_uhbd_data;
  plugin.close_file_read = close_uhbd_read;
  (*cb)(v, (vmdplugin_t *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

