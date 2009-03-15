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
 *      $Revision: 1.11 $       $Date: 2008/01/09 19:40:53 $
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
#include "endianswap.h"

#define LINESIZE 85

typedef struct {
  FILE *fd;
  int nsets;
  molfile_volumetric_t *vol;
  float scale;  // 0 for ascii files, otherwise signals binary and provides
                // scale factor for written data.
  int doswap;   // true if byteswapping needed
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
  int headersize = 0; // for binary format
  int doswap = 0; // true if byteswapping needed
  float scale=0;    // scale parameter in binary uhbd files
  
  if ((fd = fopen(filepath, "rb")) == NULL) {
    printf("uhbdplugin) Error opening file.\n");
    return NULL;
  }

  // check if the first part of file is a binary 160, which would indicate
  // that this is a binary rather than ascii UHBD file.
  fread(&headersize, sizeof(int), 1, fd);
  if (headersize == 160) {
    printf("uhbdplugin) Detected binary .grd file in native endian\n");
    doswap = 0;
  } else {
    swap4_unaligned(&headersize, 1);
    if (headersize == 160) {
      printf("uhbdplugin) Detected binary .grd file in opposite endian\n");
      doswap = 1;
    } else {
      headersize = 0;
    }
  }
  if (headersize == 160) {
    int iparams[8];
    float fparams[4];
    // we're doing binary
    // read in the entire header, plus the trailing header size
    char buf[164];
    if (fread(buf, 1, 160, fd) != 160) {
      fprintf(stderr, "uhbdplugin) Error: incomplete header in .grd file.\n");
      fclose(fd);
      return NULL;
    }
    // format of header is 72 character title, followed by:
    // scale dum2 grdflag, idum2 km one km im jm km h ox oy oz
    // The first two and last four parameters are floats, the rest ints.
    memcpy(&scale,  buf + 72, sizeof(float));
    memcpy(iparams, buf + 72 + 8, 32);
    memcpy(fparams, buf + 72 + 40, 16);
    if (doswap) {
      swap4_unaligned(&scale, 1);
      swap4_unaligned(iparams, 8);
      swap4_unaligned(fparams, 4);
    }
    xsize = iparams[5];
    ysize = iparams[6];
    zsize = iparams[7];
    ra = fparams[0];
    orig[0] = fparams[1];
    orig[1] = fparams[2];
    orig[2] = fparams[3];

  } else {
    rewind(fd);
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
  }

  /* allocate and initialize the uhbd structure */
  uhbd = new uhbd_t;
  uhbd->fd = fd;
  uhbd->vol = NULL;
  *natoms = MOLFILE_NUMATOMS_NONE;
  uhbd->nsets = 1; /* this file contains only one data set */
  uhbd->scale = scale; // set by binary files only
  uhbd->doswap = doswap;

  uhbd->vol = new molfile_volumetric_t[1];
  strcpy(uhbd->vol[0].dataname, 
      headersize ? "UHBD binary Electron Density Map"
                 : "UHBD ascii Electron Density Map");

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

  if (uhbd->scale) {
    int headerblock[6];
    // The binary data is written one z-slice at a time.  First, a 
    // block of three ints indicating which slice is next, of the form
    // kk, im, jm, where kk is the on-based slice index, and im and jm
    // are the x and y dimensions.  We just read past them and also read the
    // start of the next block, which is the actual data.
    for (z=0; z<zsize; z++) {
      if (fread(headerblock, sizeof(int), 6, fd) != 6) {
        fprintf(stderr, "uhbdplugin) Error reading header block in binary uhbd file\n");
        return MOLFILE_ERROR;
      }
      // TODO: some sanity checks on the header block.
      if (fread(datablock + z*xysize, sizeof(float), xysize, fd) != xysize) {
        fprintf(stderr, "uhbdplugin) Error reading data block in binary uhbd file\n");
        return MOLFILE_ERROR;
      }
      // read the trailing block delimiter
      fseek(fd, 4, SEEK_CUR);
    }
    if (uhbd->doswap) {
      swap4_aligned(datablock, total);
    }
    return MOLFILE_SUCCESS;
  }

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
static molfile_plugin_t plugin;

VMDPLUGIN_EXTERN int VMDPLUGIN_init(void) {
  memset(&plugin, 0, sizeof(molfile_plugin_t));
  plugin.abiversion = vmdplugin_ABIVERSION;
  plugin.type = MOLFILE_PLUGIN_TYPE;
  plugin.name = "uhbd";
  plugin.prettyname = "UHBD Grid";
  plugin.author = "Alexander Spaar, Justin Gullingsrud";
  plugin.majorv = 0;
  plugin.minorv = 4;
  plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
  plugin.filename_extension = "uhbdgrd,grd";
  plugin.open_file_read = open_uhbd_read;
  plugin.read_volumetric_metadata = read_uhbd_metadata;
  plugin.read_volumetric_data = read_uhbd_data;
  plugin.close_file_read = close_uhbd_read;
  return VMDPLUGIN_SUCCESS;
}


VMDPLUGIN_EXTERN int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)(void *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_EXTERN int VMDPLUGIN_fini(void) { return VMDPLUGIN_SUCCESS; }

