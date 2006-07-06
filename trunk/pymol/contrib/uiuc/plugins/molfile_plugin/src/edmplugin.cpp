/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_edmplugin
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
 *      $RCSfile: edmplugin.C,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.25 $       $Date: 2006/02/23 19:36:44 $
 *
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "molfile_plugin.h"

typedef struct {
  FILE *fd;
  int nsets;
  molfile_volumetric_t *vol;
} edm_t;

static void eatline(FILE * fd) {
  char readbuf[1025];
  fgets(readbuf, 1024, fd);    // go on to next line
}  

static void *open_edm_read(const char *filepath, const char *filetype,
    int *natoms) {
  FILE *fd;
  edm_t *edm;
  int ntitle, na, nb, nc, xsize, ysize, zsize;
  int amin, amax, bmin, bmax, cmin, cmax;
  float a, b, c, alpha, beta, gamma;
  float xdelta, ydelta, zdelta;
  float alpha1, beta1, gamma1;
  float xaxis[3], yaxis[3], zaxis[3], z1, z2, z3;
  int i, convcnt;
  
  fd = fopen(filepath, "rb");
  if (!fd) 
    return NULL;

  edm = new edm_t;
  edm->fd = fd;
  edm->vol = NULL;
  *natoms = MOLFILE_NUMATOMS_NONE;

  edm->vol = new molfile_volumetric_t[1];
 
  edm->nsets = 1; // this EDM file contains only one data set

  // read in EDM file header information
  eatline(edm->fd);               // skip first header line 

  convcnt = fscanf(edm->fd, "%d", &ntitle); // read number of title lines
  if (convcnt != 1) {
    printf("edmplugin: failed to read in title line count\n");
    fclose(edm->fd);
    delete [] edm->vol;
    delete edm;
    return NULL;
  }
    
  eatline(edm->fd);               // go on to next line

  // skip past title and comment lines in header
  for (i=0; i<ntitle; i++) {
    eatline(edm->fd);             // skip a line
  }

  // read in the box dimensions and grid spacing deltas
  convcnt = fscanf(edm->fd, "%d %d %d %d %d %d %d %d %d",
         &na, &amin, &amax, &nb, &bmin, &bmax, &nc, &cmin, &cmax);
  if (convcnt != 9) {
    printf("edmplugin: failed to read in box dimensions\n");
    fclose(edm->fd);
    delete [] edm->vol;
    delete edm;
    return NULL;
  }

  eatline(edm->fd);               // go on to next line
  
  // calculate number of samples in each dimension
  xsize = amax - amin + 1;    
  ysize = bmax - bmin + 1;    
  zsize = cmax - cmin + 1;    
  edm->vol[0].xsize = xsize;
  edm->vol[0].ysize = ysize;
  edm->vol[0].zsize = zsize;
  edm->vol[0].has_color = 0;

  // read in 6 values for unit cell box orientation 
  convcnt = fscanf(edm->fd, "%f %f %f %f %f %f", 
                   &a, &b, &c, &alpha, &beta, &gamma);
  if (convcnt != 6) {
    printf("edmplugin: failed to read in box lengths and angles\n");
    fclose(edm->fd);
    delete [] edm->vol;
    delete edm;
    return NULL;
  }
  eatline(edm->fd);               // go on to next line

  // find box coordinates 
  xdelta = a / (float) na;
  ydelta = b / (float) nb;
  zdelta = c / (float) nc;

  strcpy(edm->vol[0].dataname, "X-PLOR Electron Density Map");

  // convert degrees to radians
  alpha1 = 3.14159265358979323846 * alpha / 180.0;
  beta1  = 3.14159265358979323846 *  beta / 180.0;
  gamma1 = 3.14159265358979323846 * gamma / 180.0;

  // calculate non-orthogonal unit cell coordinates
  xaxis[0] = xdelta;
  xaxis[1] = 0;
  xaxis[2] = 0;

  yaxis[0] = cos(gamma1) * ydelta;
  yaxis[1] = sin(gamma1) * ydelta;
  yaxis[2] = 0;

  z1 = cos(beta1);
  z2 = (cos(alpha1) - cos(beta1)*cos(gamma1)) / sin(gamma1);
  z3 = sqrt(1.0 - z1*z1 - z2*z2);
  zaxis[0] = z1 * zdelta;
  zaxis[1] = z2 * zdelta;
  zaxis[2] = z3 * zdelta;

  edm->vol[0].origin[0] = xaxis[0] * amin + yaxis[0] * bmin + zaxis[0] * cmin;
  edm->vol[0].origin[1] = yaxis[1] * bmin + zaxis[1] * cmin;
  edm->vol[0].origin[2] = zaxis[2] * cmin;

  edm->vol[0].xaxis[0] = xaxis[0] * xsize;
  edm->vol[0].xaxis[1] = 0;
  edm->vol[0].xaxis[2] = 0;

  edm->vol[0].yaxis[0] = yaxis[0] * ysize;
  edm->vol[0].yaxis[1] = yaxis[1] * ysize;
  edm->vol[0].yaxis[2] = 0;

  edm->vol[0].zaxis[0] = zaxis[0] * zsize;
  edm->vol[0].zaxis[1] = zaxis[1] * zsize;
  edm->vol[0].zaxis[2] = zaxis[2] * zsize;

  // Check that the EDM file is stored in the "ZYX" format we expect,
  // and return NULL if it is not a supported file type.
  char planeorder[4];
  memset(planeorder, 0, sizeof(planeorder));
  convcnt = fscanf(edm->fd, "%3s", planeorder);
  if (convcnt != 1) {
    printf("edmplugin: failed to read in plane order\n");
    fclose(edm->fd);
    delete [] edm->vol;
    delete edm;
    return NULL;
  }

  if (strcmp(planeorder, "ZYX")) { 
    printf("edmplugin: unsupported plane ordering %s\n", planeorder);
    fclose(edm->fd);
    delete [] edm->vol;
    delete edm;
    return NULL;
  }
  eatline(edm->fd);               // go on to next line

  return edm;
}

static int read_edm_metadata(void *v, int *nsets, 
  molfile_volumetric_t **metadata) {
  edm_t *edm = (edm_t *)v;
  *nsets = edm->nsets; 
  *metadata = edm->vol;  

  return MOLFILE_SUCCESS;
}

static int read_edm_data(void *v, int set, float *datablock,
                         float *colorblock) {
  edm_t *edm = (edm_t *)v;
  float * cell = datablock;
  int z, l, sentinel, convcnt;
  char readbuf[16];
 
  int xsize = edm->vol[0].xsize;
  int ysize = edm->vol[0].ysize;
  int zsize = edm->vol[0].zsize;

  // number of lines of text per slice
  int nperslice = xsize * ysize;
  int nlines = (int)(nperslice / 6.0);
  if ((nlines * 6) < nperslice) {
    nlines++;
  }
  int leftover = (nperslice - (nlines - 1) * 6);

  for (z=0; z<zsize; z++) {
    int c;
    eatline(edm->fd);                // Eat the Z-plane index and throw away

#if 1
    // read one plane of data, once cell at a time
    for (c=0; c<nperslice; c++) {
      convcnt = fscanf(edm->fd, "%f", cell);
      if (convcnt != 1) {
        printf("edmplugin: failed reading cell data\n");
        printf("edmplugin: cell %d of %d, slice %d\n", c, nperslice, z);
        return MOLFILE_ERROR; // bad file format encountered
      }
      cell++;
    }
    eatline(edm->fd);                // go on to next line
#else
    for (l=1; l<nlines; l++) {
      for (c=0; c<6; c++) {
        fgets(readbuf, 13, edm->fd); // read in 12 chars (fgets reads N-1)
        convcnt = sscanf(readbuf, "%f", cell); // convert ascii to float
        if (convcnt != 1) {
          printf("edmplugin: failed reading cell data\n");
          printf("edmplugin: cell on line %d cell %d, of %d lines\n", l, c, nlines);
          return MOLFILE_ERROR; // bad file format encountered
        }
        cell++;
      }
      eatline(edm->fd);              // go on to next line
    }

    // read any remaining partial line of text for this plane.
    if (leftover > 0) {
      for (c=0; c<leftover; c++) {
        fgets(readbuf, 13, edm->fd); 
        convcnt = sscanf(readbuf, "%f", cell);
        if (convcnt != 1) {
          printf("edmplugin: failed reading partial line cell data\n");
          return MOLFILE_ERROR; // bad file format encountered
        }
        cell++;
      }
    } 
#endif
  }

  // read the -9999 end-of-file sentinel record
  sentinel = 0;
  fgets(readbuf, 13, edm->fd);
  sscanf(readbuf, "%d", &sentinel);  
  if (sentinel != -9999) {
    printf("edmplugin: EOF sentinel != -9999\n");
    // return MOLFILE_ERROR; // bad file format encountered, no sentinel record
  }
 
  return MOLFILE_SUCCESS;
}

static void close_edm_read(void *v) {
  edm_t *edm = (edm_t *)v;

  fclose(edm->fd);
  delete [] edm->vol; 
  delete edm;
}

/*
 * Initialization stuff here
 */
static molfile_plugin_t plugin = {
  vmdplugin_ABIVERSION,         // ABI version
  MOLFILE_PLUGIN_TYPE,          // plugin type
  "edm",                        // file format description
  "XPLOR Electron Density Map", // file format description
  "John E. Stone",              // author(s)
  0,                            // major version
  4,                            // minor version
  VMDPLUGIN_THREADSAFE,         // is reentrant
  "edm"                         // filename extension
};

VMDPLUGIN_EXTERN int VMDPLUGIN_init(void) { return VMDPLUGIN_SUCCESS; }
VMDPLUGIN_EXTERN int VMDPLUGIN_fini(void) { return VMDPLUGIN_SUCCESS; }
VMDPLUGIN_EXTERN int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  plugin.open_file_read = open_edm_read;
  plugin.read_volumetric_metadata = read_edm_metadata;
  plugin.read_volumetric_data = read_edm_data;
  plugin.close_file_read = close_edm_read;
  (*cb)(v, (vmdplugin_t *)&plugin);
  return VMDPLUGIN_SUCCESS;
}



#ifdef TEST_EDMPLUGIN

int main(int argc, char *argv[]) {
  int natoms;
  void *v;
  int i, nsets, set;
  molfile_volumetric_t * meta;

  while (--argc) {
    ++argv;
    v = open_edm_read(*argv, "edm", &natoms);
    if (!v) {
      fprintf(stderr, "open_edm_read failed for file %s\n", *argv);
      return 1;
    }

    // try loading the EDM metadata now
    if (read_edm_metadata(v, &nsets, &meta)) {
      return 1; // failed to load edm file
    }

    for (set=0; set<nsets; set++) {
      printf("Loading volume set: %d\n", set);   
      
      int elements = meta[set].xsize * meta[set].ysize * meta[set].zsize;
      printf("   Grid Elements: %d\n", elements);
      printf(" Grid dimensions: X: %d Y: %d Z: %d\n", 
             meta[set].xsize, meta[set].ysize, meta[set].zsize);

      float * voldata = (float *) malloc(sizeof(float) * elements);
      float * coldata = NULL;

      if (meta[set].has_color) {
        coldata = (float *) malloc(sizeof(float) * elements * 3);
      }

      // try loading the EDM data sets now
      if (read_edm_data(v, set, voldata, coldata)) {
        return 1; // failed to load edm file
      }

      printf("First 6 elements:\n   ");
      for (i=0; i<6; i++) {
        printf("%g, ", voldata[i]);
      }
      printf("\n"); 

      printf("Last 6 elements:\n   ");
      for (i=elements - 6; i<elements; i++) {
        printf("%g, ", voldata[i]);
      }
      printf("\n"); 
    }

    close_edm_read(v);
  }
  return 0;
}

#endif




