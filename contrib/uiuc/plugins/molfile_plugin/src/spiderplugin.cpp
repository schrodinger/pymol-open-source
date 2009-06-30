/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_spiderplugin
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
 *      $RCSfile: spiderplugin.C,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.15 $       $Date: 2009/04/29 15:45:34 $
 *
 ***************************************************************************/

/* 
 * SPIDER volumetric image datasets
 *   http://www.wadsworth.org/spider_doc/spider/docs/image_doc.html
 *
 * TODO:
 *  - Add code to determine axis scaling factors, axis angles, offsets, etc
 */

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include "endianswap.h"

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
  int byteswap;
  int nslice;
  int nrow;
  int nlabel;
  int iform;
  int imami;
  float fmax;
  float fmin;
  float av;
  float sig;
  int nsam;
  int headrec;
  int iangle;
  float phi;
  float theta;
  float gamma;
  float xoffset;
  float yoffset;
  float zoffset;
  float scale;
  int headbyt;
  int reclen;
  int nstack;
  int inuse;
  int maxim; 
} spider_t;


static void *open_spider_read(const char *filepath, const char *filetype,
    int *natoms) {
  FILE *fd;
  spider_t *vol;
  int total;

  /* file header buffer union */ 
  union buffer {
    float fbuf[256];  
    char  cbuf[1024];
  } h;
 
  fd = fopen(filepath, "rb");
  if (!fd) {
    fprintf(stderr, "spiderplugin) Error opening file.\n");
    return NULL;
  }

  /* allocate and initialize the spider structure */
  vol = new spider_t;
  vol->fd = fd;
  vol->vol = NULL;
  vol->byteswap = 0;
  *natoms = MOLFILE_NUMATOMS_NONE;
  vol->nsets = 1; /* this file contains only one data set */

  vol->vol = new molfile_volumetric_t[1];
  strcpy(vol->vol[0].dataname, "SPIDER map");

  // read SPIDER file header
  if (fread(&h.cbuf, 1024, 1, fd) < 1) {
    printf("spiderplugin) failed to read file header\n");
    return NULL; 
  } 

  // perform sanity checks on header values to see if we 
  // need to do byte swapping, or abort.
  vol->nslice   = (h.fbuf[0] < 0) ? -h.fbuf[0] : h.fbuf[0];
  vol->nrow     = h.fbuf[1];
  vol->nsam     = h.fbuf[11];
  total = vol->nslice * vol->nrow * vol->nsam;

  if (total <= 0 || 
      vol->nsam   <= 0 || vol->nsam   > 100000 ||
      vol->nrow   <= 0 || vol->nrow   > 100000 ||
      vol->nslice <= 0 || vol->nslice > 100000) { 

    printf("spiderplugin) Non-native endianness or unusual file detected\n");

    // byte swap the entire header in hopes of making sense of this gibberish
    vol->byteswap = 1;
    swap4_aligned(&h.fbuf, 256);

    vol->nslice   = (h.fbuf[0] < 0) ? -h.fbuf[0] : h.fbuf[0];
    vol->nrow     = h.fbuf[1];
    vol->nsam     = h.fbuf[11];
    total = vol->nslice * vol->nrow * vol->nsam;

    // check to see if we still have gibberish or not, bail out if we do.
    if (total <= 0 || 
      vol->nsam   <= 0 || vol->nsam   > 100000 ||
      vol->nrow   <= 0 || vol->nrow   > 100000 ||
      vol->nslice <= 0 || vol->nslice > 100000) { 
      printf("spiderplugin) bad header values in file fail sanity checks\n");
      delete [] vol->vol;
      delete vol;
      return NULL;
    }
  }
  if (vol->byteswap) {
    printf("spiderplugin) Enabling byte swapping\n");
  }

  vol->nlabel   = h.fbuf[3];
  vol->iform    = h.fbuf[4];
  vol->imami    = h.fbuf[5];
  vol->fmax     = h.fbuf[6];
  vol->fmin     = h.fbuf[7];
  vol->av       = h.fbuf[8];
  vol->sig      = h.fbuf[9];
  vol->headrec  = h.fbuf[12];
  vol->iangle   = h.fbuf[13];
  vol->phi      = h.fbuf[14];
  vol->theta    = h.fbuf[15];
  vol->gamma    = h.fbuf[16];
  vol->xoffset  = h.fbuf[17];
  vol->yoffset  = h.fbuf[18];
  vol->zoffset  = h.fbuf[19];
  vol->scale    = h.fbuf[20];
  vol->headbyt  = h.fbuf[21];
  vol->reclen   = h.fbuf[22];
  vol->nstack   = h.fbuf[23];
  vol->inuse    = h.fbuf[24];
  vol->maxim    = h.fbuf[25];

printf("spider  nslice: %d\n", vol->nslice);
printf("spider    nrow: %d\n", vol->nrow);
printf("spider    nsam: %d\n", vol->nsam);
printf("spider   iform: %d\n", vol->iform);
printf("spider   scale: %f\n", vol->scale);
printf("spider xoffset: %f\n", vol->xoffset);
printf("spider yoffset: %f\n", vol->yoffset);
printf("spider zoffset: %f\n", vol->zoffset);
printf("spider     phi: %f\n", vol->phi);
printf("spider   theta: %f\n", vol->theta);
printf("spider   gamma: %f\n", vol->gamma);

  /* correct bad headbyt and reclen SPIDER files */
  if (vol->iform < 4 && (vol->reclen < (vol->nsam * 4)))
    vol->reclen = vol->nsam * 4; 

  int headrec = 1024 / vol->reclen;
  if (vol->reclen < 1024 && (1024 % (vol->reclen)) != 0)
     headrec++;
  int headbyt = headrec * vol->reclen;
 
  if (vol->iform < 4 && (vol->headbyt < headbyt))
    vol->headbyt = headbyt;

printf("spider headbyt: %d\n", vol->headbyt);

  /* seek to data offset */
  fseek(fd, vol->headbyt, SEEK_SET);

  /* SPIDER files contain no color information */
  vol->vol[0].has_color = 0;

  vol->vol[0].xsize = vol->nsam;
  vol->vol[0].ysize = vol->nrow;
  vol->vol[0].zsize = vol->nslice;

  /* Set the unit cell origin and basis vectors */
  float vz[3] = {0.0, 0.0, 0.0};
  memcpy(vol->vol[0].xaxis, &vz, sizeof(vz));
  memcpy(vol->vol[0].yaxis, &vz, sizeof(vz));
  memcpy(vol->vol[0].zaxis, &vz, sizeof(vz));

  /* the scale value may be zero, if so, just reset to 1.0 */
  float vscale = vol->scale;
  if (vscale == 0.0) 
    vscale = 1.0;

  /* the data is stored in y/x/-z order and coordinate handedness  */
  /* we should probably rewrite the loader loop to shuffle x/y     */
  /* so that future conversions to other formats don't leave it in */
  /* an unusual packing order.  For now this works however         */

  float xlen = vscale * (vol->vol[0].ysize-1);
  float ylen = vscale * (vol->vol[0].xsize-1);
  float zlen = vscale * (vol->vol[0].zsize-1);

  vol->vol[0].xaxis[1] =  xlen;
  vol->vol[0].yaxis[0] =  ylen;
  vol->vol[0].zaxis[2] = -zlen;

  vol->vol[0].origin[0] = vol->yoffset - (0.5 * ylen);
  vol->vol[0].origin[1] = vol->xoffset - (0.5 * xlen);
  vol->vol[0].origin[2] = vol->zoffset + (0.5 * zlen);

printf("spider final offset: (%f, %f, %f)\n",
  vol->vol[0].origin[0], vol->vol[0].origin[1], vol->vol[0].origin[2]);

printf("spider final axes:\n");
printf("  X (%f, %f, %f)\n",
  vol->vol[0].xaxis[0], 
  vol->vol[0].xaxis[1], 
  vol->vol[0].xaxis[2]);

printf("  Y (%f, %f, %f)\n",
  vol->vol[0].yaxis[0], 
  vol->vol[0].yaxis[1], 
  vol->vol[0].yaxis[2]);

printf("  Z (%f, %f, %f)\n",
  vol->vol[0].zaxis[0], 
  vol->vol[0].zaxis[1], 
  vol->vol[0].zaxis[2]);

  return vol;
}

static int read_spider_metadata(void *v, int *nsets, 
  molfile_volumetric_t **metadata) {
  spider_t *vol = (spider_t *)v;
  *nsets = vol->nsets;
  *metadata = vol->vol;  

  return MOLFILE_SUCCESS;
}

static int read_spider_data(void *v, int set, float *datablock,
                         float *colorblock) {
  spider_t *vol = (spider_t *)v;
  FILE *fd = vol->fd;
  int x, y, z, xsize, ysize, zsize, xysize, total;

  xsize = vol->vol[0].xsize;
  ysize = vol->vol[0].ysize;
  zsize = vol->vol[0].zsize;
  xysize = xsize * ysize;
  total = xysize * zsize;

  // Read the values from the file
  fread(datablock, total * sizeof(float), 1, fd);

  // perform byte swapping if necessary
  if (vol->byteswap) 
    swap4_aligned(datablock, total);

  return MOLFILE_SUCCESS;
}

static void close_spider_read(void *v) {
  spider_t *vol = (spider_t *)v;
  
  fclose(vol->fd);
  if (vol->vol != NULL)
    delete [] vol->vol; 
  delete vol;
}

/*
 * Initialization stuff here
 */
static molfile_plugin_t plugin;

VMDPLUGIN_EXTERN int VMDPLUGIN_init(void) {
  memset(&plugin, 0, sizeof(molfile_plugin_t));
  plugin.abiversion = vmdplugin_ABIVERSION;
  plugin.type = MOLFILE_PLUGIN_TYPE;
  plugin.name = "spider";
  plugin.prettyname = "SPIDER Density Map";
  plugin.author = "John Stone";
  plugin.majorv = 0;
  plugin.minorv = 6;
  plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
  plugin.filename_extension = "spi,spider";
  plugin.open_file_read = open_spider_read;
  plugin.read_volumetric_metadata = read_spider_metadata;
  plugin.read_volumetric_data = read_spider_data;
  plugin.close_file_read = close_spider_read;
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_EXTERN int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_EXTERN int VMDPLUGIN_fini(void) { return VMDPLUGIN_SUCCESS; }

