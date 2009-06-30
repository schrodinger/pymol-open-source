/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_situsplugin
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
 *      $RCSfile: situsplugin.C,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.15 $       $Date: 2009/04/29 15:45:33 $
 *
 ***************************************************************************/

/* 
 * Situs EM map file reader
 *   map is a cubic lattice
 *   all data is stored in plain ASCII for convenience
 *
 * Format of the file is:
 *   voxel size in Angstroms
 *   coordinates of first voxel
 *   integer X/Y/Z counts
 *   voxels follow in X fastest, Y next fastest, Z slowest
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#if defined(_AIX)
#include <strings.h>
#endif

#define MIN(X,Y) (((X)<(Y))? (X) : (Y))
#define MAX(X,Y) (((X)>(Y))? (X) : (Y))

#define TOLERANCE 1e-4

#ifndef NAN //not a number
  const float NAN = sqrtf(-1.f); //need some kind of portable NAN definition
#endif

#include "molfile_plugin.h"

typedef struct {
  FILE *fd;
  int nsets;
  molfile_volumetric_t *vol;
} situs_t;


static void *open_situs_read(const char *filepath, const char *filetype,
    int *natoms) {
  FILE *fd;
  situs_t *situs;
  float scale;
  int xsize, ysize, zsize;
  float orig[3]; 
  
  fd = fopen(filepath, "r");
  if (!fd) {
    printf("situsplugin) Error opening file.\n");
    return NULL;
  }

  /* get the scale */
  if (fscanf(fd, "%f", &scale) != 1) {;
    printf("situsplugin) Error reading voxel scale.\n");
    return NULL;
  }

  if (fscanf(fd, "%f %f %f", orig, orig+1, orig+2) != 3) {
    printf("situsplugin) Error reading grid origin.\n");
    return NULL;
  }

  /* get the number of grid points */
  if (fscanf(fd, "%d %d %d", &xsize, &ysize, &zsize) != 3) {
    printf("situsplugin) Error reading grid dimensions.\n");
    return NULL;
  }

  /* allocate and initialize the situs structure */
  situs = new situs_t;
  situs->fd = fd;
  situs->vol = NULL;
  *natoms = MOLFILE_NUMATOMS_NONE;
  situs->nsets = 1; /* this file contains only one data set */

  situs->vol = new molfile_volumetric_t[1];
  strcpy(situs->vol[0].dataname, "Situs map");

  /* Set the unit cell origin and basis vectors */
  for (int i=0; i<3; i++) {
    situs->vol[0].origin[i] = orig[i];
    situs->vol[0].xaxis[i] = 0.0;
    situs->vol[0].yaxis[i] = 0.0;
    situs->vol[0].zaxis[i] = 0.0;
  }
  situs->vol[0].xaxis[0] = scale * (xsize-1);
  situs->vol[0].yaxis[1] = scale * (ysize-1);
  situs->vol[0].zaxis[2] = scale * (zsize-1);

  situs->vol[0].xsize = xsize;
  situs->vol[0].ysize = ysize;
  situs->vol[0].zsize = zsize;

  situs->vol[0].has_color = 0; /* Situs maps contain no color information. */

  return situs;
}

static int read_situs_metadata(void *v, int *nsets, 
  molfile_volumetric_t **metadata) {
  situs_t *situs = (situs_t *)v;
  *nsets = situs->nsets; 
  *metadata = situs->vol;  

  return MOLFILE_SUCCESS;
}

static int read_situs_data(void *v, int set, float *datablock,
                         float *colorblock) {
  situs_t *situs = (situs_t *)v;
  FILE *fd = situs->fd;
  int xsize, ysize, zsize, xysize, count, total;

  xsize = situs->vol[0].xsize;
  ysize = situs->vol[0].ysize;
  zsize = situs->vol[0].zsize;
  xysize = xsize * ysize;
  total = xysize * zsize;

  /* Read the values from the file */
  for (count=0; count<total; count++) {
    if (fscanf(fd, "%f", datablock + count) != 1) {
      printf("situsplugin) Failed reading situs map data\n");
      return MOLFILE_ERROR;
    }
  }   

  return MOLFILE_SUCCESS;
}

static void close_situs_read(void *v) {
  situs_t *situs = (situs_t *)v;
  
  fclose(situs->fd);
  if (situs->vol != NULL)
    delete [] situs->vol; 
  delete situs;
}

static void *open_situs_write(const char *filepath, const char *filetype, int natoms) {

  FILE *fd = fopen(filepath, "w");
  if (!fd) {
    fprintf(stderr, "situsplugin) Could not open path '%s' for writing.\n",
      filepath);
  }
  return fd;
}

static void close_situs_write(void *v) {
  if (v) {
    fclose((FILE *)v);
  }
}

/* The Situs format requires the same grid spacing in all directions.
 * The following code allows us to re-sample the map so that we can 
 * write a Situs file when the uniform grid spacing requirement is not met.
 * */

/// return voxel, after safely clamping index to valid range
float situs_voxel_value_safe(int x, int y, int z, const int xsize, const int ysize, const int zsize, const float *data) {
  int xx, yy, zz; 
  xx = (x > 0) ? ((x < xsize) ? x : xsize-1) : 0;
  yy = (y > 0) ? ((y < ysize) ? y : ysize-1) : 0;
  zz = (z > 0) ? ((z < zsize) ? z : zsize-1) : 0;
  int index = zz*xsize*ysize + yy*xsize + xx;
  return data[index];
}

/// return interpolated value from 8 nearest neighbor voxels
float situs_voxel_value_interpolate(float xv, float yv, float zv, const int xsize, const int ysize, const int zsize, const float *data) {
  int x = (int) xv;
  int y = (int) yv;
  int z = (int) zv;
  // fractional offset
  float xf = xv - x;
  float yf = yv - y;
  float zf = zv - z;
  float xlerps[4];
  float ylerps[2];
  float tmp;

  tmp = situs_voxel_value_safe(x, y, z, xsize, ysize, zsize, data);
  xlerps[0] = tmp + xf*(situs_voxel_value_safe(x+1, y, z, xsize, ysize, zsize, data) - tmp);

  tmp = situs_voxel_value_safe(x, y+1, z, xsize, ysize, zsize, data);
  xlerps[1] = tmp + xf*(situs_voxel_value_safe(x+1, y+1, z, xsize, ysize, zsize, data) - tmp);

  tmp = situs_voxel_value_safe(x, y, z+1, xsize, ysize, zsize, data);
  xlerps[2] = tmp + xf*(situs_voxel_value_safe(x+1, y, z+1, xsize, ysize, zsize, data) - tmp);

  tmp = situs_voxel_value_safe(x, y+1, z+1, xsize, ysize, zsize, data);
  xlerps[3] = tmp + xf*(situs_voxel_value_safe(x+1, y+1, z+1, xsize, ysize, zsize, data) - tmp);

  ylerps[0] = xlerps[0] + yf*(xlerps[1] - xlerps[0]);
  ylerps[1] = xlerps[2] + yf*(xlerps[3] - xlerps[2]);

  return ylerps[0] + zf*(ylerps[1] - ylerps[0]);
}

/// return interpolated value of voxel, based on atomic coords.
/// XXX need to account for non-orthog. cells
float situs_voxel_value_interpolate_from_coord(float xpos, float ypos, float zpos, const float *origin, const float *xdelta, const float *ydelta, const float *zdelta, const int xsize, const int ysize, const int zsize, float *data) {
  xpos = (xpos-origin[0])/xdelta[0];
  ypos = (ypos-origin[1])/ydelta[1];
  zpos = (zpos-origin[2])/zdelta[2];
  int gx = (int) xpos; // XXX this is wrong for non-orthog cells.
  int gy = (int) ypos;
  int gz = (int) zpos;
  if (gx < 0 || gx >= xsize) return NAN;
  if (gy < 0 || gy >= ysize) return NAN;
  if (gz < 0 || gz >= zsize) return NAN;
    
  return situs_voxel_value_interpolate(xpos, ypos, zpos, xsize, ysize, zsize, data);

}

static int write_situs_data(void *v, molfile_volumetric_t *metadata, float *datablock, float *colorblock) {

  FILE *fd = (FILE *)v;
  const int xsize = metadata->xsize;
  const int ysize = metadata->ysize;
  const int zsize = metadata->zsize;
  const int xysize = xsize * ysize;

  float xaxis[3], yaxis[3], zaxis[3];
  float xdelta[3], ydelta[3], zdelta[3];
  float origin[3];

  int i, j, k;

  for (i=0; i<3; i++) {
    origin[i] = metadata->origin[i];
    xaxis[i] = metadata->xaxis[i];
    yaxis[i] = metadata->yaxis[i];
    zaxis[i] = metadata->zaxis[i];
    xdelta[i] = xaxis[i]/(xsize - 1);
    ydelta[i] = yaxis[i]/(ysize - 1);
    zdelta[i] = zaxis[i]/(zsize - 1);
  }

  /* Situs format requires an orthogonal cell */
  if (fabs(xaxis[1]) > TOLERANCE || fabs(xaxis[2]) > TOLERANCE ||
      fabs(yaxis[0]) > TOLERANCE || fabs(yaxis[2]) > TOLERANCE ||
      fabs(zaxis[0]) > TOLERANCE || fabs(zaxis[1]) > TOLERANCE) {
    fprintf(stderr, "situsplugin) Could not write situs file: this format requires an orthogonal cell.\n");
    return MOLFILE_ERROR;
  }

  /* Situs format requires the same grid spacing in all dimensions */
  float xres = xdelta[0]*xdelta[0] + xdelta[1]*xdelta[1] + xdelta[2]*xdelta[2];
  float yres = ydelta[0]*ydelta[0] + ydelta[1]*ydelta[1] + ydelta[2]*ydelta[2];
  float zres = zdelta[0]*zdelta[0] + zdelta[1]*zdelta[1] + zdelta[2]*zdelta[2];
  if ( (fabs(xres-yres) > TOLERANCE) || (fabs(xres-zres) > TOLERANCE) ) {

    //fprintf(stderr, "situsplugin) Could not write situs file: this format requires the same grid spacing in all dimensions.\n");
    //return MOLFILE_ERROR;
    
    fprintf(stderr, "situsplugin) Warning: This format requires the same grid spacing in all dimensions. The map will be re-sampled to meet this requirement. The resulting cell may be slightly smaller than the original one.\n");

    // Use the smallest grid spacing for the re-sampled map
    float new_delta = MIN(xdelta[0], ydelta[1]);
    new_delta = MIN(new_delta, zdelta[2]);

    // Calculate the map dimensions for the new grid spacing 
    int new_xsize = int(xaxis[0]/new_delta);
    int new_ysize = int(yaxis[1]/new_delta);
    int new_zsize = int(zaxis[2]/new_delta);
    int new_xysize = new_xsize * new_ysize;
    int new_size = new_xsize * new_ysize * new_zsize;

    // Resample map
    float *new_data = (float *)malloc(3*new_size*sizeof(float));
    for (i=0; i<new_xsize; i++) {
      float xpos = origin[0] + i*new_delta;
      for (j=0; j<new_ysize; j++) {
        float ypos = origin[1] + j*new_delta;
        for (k=0; k<new_zsize; k++) {
          float zpos = origin[2] + k*new_delta;
          new_data[i + j*new_xsize + k*new_xysize] = situs_voxel_value_interpolate_from_coord(xpos, ypos, zpos, origin, xdelta, ydelta, zdelta, xsize, ysize, zsize, datablock);
        }
      }
    }

    /* Write a situs header */
    fprintf(fd, "%g %g %g %g %d %d %d\n\n", new_delta, origin[0], origin[1], origin[2], new_xsize, new_ysize, new_zsize);
  
    /* Write the main data array */
    int numentries = 1;
    for (k=0; k<new_zsize; k++) {
      for (j=0; j<new_ysize; j++) {
        for (i=0; i<new_xsize; i++) {
          fprintf(fd, "%g ", new_data[i + j*new_xsize + k*new_xysize]);
          if (numentries % 10 == 0) fprintf(fd, "\n");
          numentries++;
        }
      }
    }

    free(new_data);

  } else {
  
    /* Write a situs header */
    fprintf(fd, "%g %g %g %g %d %d %d\n\n", xdelta[0], origin[0], origin[1], origin[2], xsize, ysize, zsize);
  
    /* Write the main data array */
    int numentries = 1;
    for (k=0; k<zsize; k++) {
      for (j=0; j<ysize; j++) {
        for (i=0; i<xsize; i++) {
          fprintf(fd, "%g ", datablock[i + j*xsize + k*xysize]);
          if (numentries % 10 == 0) fprintf(fd, "\n");
          numentries++;
        }
      }
    }

  }

  fflush(fd);
  return MOLFILE_SUCCESS;
}

/*
 * Initialization stuff here
 */
static molfile_plugin_t plugin;

VMDPLUGIN_EXTERN int VMDPLUGIN_init(void) {
  memset(&plugin, 0, sizeof(molfile_plugin_t));
  plugin.abiversion = vmdplugin_ABIVERSION;
  plugin.type = MOLFILE_PLUGIN_TYPE;
  plugin.name = "situs";
  plugin.prettyname = "Situs Density Map";
  plugin.author = "John Stone, Leonardo Trabuco";
  plugin.majorv = 1;
  plugin.minorv = 5;
  plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
  plugin.filename_extension = "sit,situs";
  plugin.open_file_read = open_situs_read;
  plugin.read_volumetric_metadata = read_situs_metadata;
  plugin.read_volumetric_data = read_situs_data;
  plugin.close_file_read = close_situs_read;
#if vmdplugin_ABIVERSION > 9
  plugin.open_file_write = open_situs_write;
  plugin.write_volumetric_data = write_situs_data;
  plugin.close_file_write = close_situs_write;
#endif
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_EXTERN int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_EXTERN int VMDPLUGIN_fini(void) { return VMDPLUGIN_SUCCESS; }

