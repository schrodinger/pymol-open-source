/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_edmplugin
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
 *      $RCSfile: edmplugin.C,v $
 *      $Author: ltrabuco $       $Locker:  $             $State: Exp $
 *      $Revision: 1.32 $       $Date: 2009/06/13 18:25:06 $
 *
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "molfile_plugin.h"

#ifndef NAN //not a number
  const float NAN = sqrtf(-1.f); //need some kind of portable NAN definition
#endif

#define TOLERANCE 1e-4

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
    printf("edmplugin) failed to read in title line count\n");
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
    printf("edmplugin) failed to read in box dimensions\n");
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
    printf("edmplugin) failed to read in box lengths and angles\n");
    fclose(edm->fd);
    delete [] edm->vol;
    delete edm;
    return NULL;
  }
  eatline(edm->fd);            // go on to next line

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

  edm->vol[0].xaxis[0] = xaxis[0] * (xsize-1);
  edm->vol[0].xaxis[1] = 0;
  edm->vol[0].xaxis[2] = 0;

  edm->vol[0].yaxis[0] = yaxis[0] * (ysize-1);
  edm->vol[0].yaxis[1] = yaxis[1] * (ysize-1);
  edm->vol[0].yaxis[2] = 0;

  edm->vol[0].zaxis[0] = zaxis[0] * (zsize-1);
  edm->vol[0].zaxis[1] = zaxis[1] * (zsize-1);
  edm->vol[0].zaxis[2] = zaxis[2] * (zsize-1);

  // Check that the EDM file is stored in the "ZYX" format we expect,
  // and return NULL if it is not a supported file type.
  char planeorder[4];
  memset(planeorder, 0, sizeof(planeorder));
  convcnt = fscanf(edm->fd, "%3s", planeorder);
  if (convcnt != 1) {
    printf("edmplugin) failed to read in plane order\n");
    fclose(edm->fd);
    delete [] edm->vol;
    delete edm;
    return NULL;
  }

  if (strcmp(planeorder, "ZYX")) { 
    printf("edmplugin) unsupported plane ordering %s\n", planeorder);
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
        printf("edmplugin) failed reading cell data\n");
        printf("edmplugin) cell %d of %d, slice %d\n", c, nperslice, z);
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
          printf("edmplugin) failed reading cell data\n");
          printf("edmplugin) cell on line %d cell %d, of %d lines\n", l, c, nlines);
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
          printf("edmplugin) failed reading partial line cell data\n");
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
    printf("edmplugin) EOF sentinel != -9999\n");
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

static void *open_edm_write(const char *filepath, const char *filetype, int natoms) {

  FILE *fd = fopen(filepath, "w");
  if (!fd) {
    fprintf(stderr, "edmplugin) Could not open path '%s' for writing.\n",
      filepath);
  }
  return fd;
}

static void close_edm_write(void *v) {
  if (v) {
    fclose((FILE *)v);
  }
}

// XXX - The following interpolation code is duplicated from situsplugin.C
//       (which is in turn duplicated from VMD). The only difference is that
//       we are returning zeroes instead of NANs for out-of-range queries
//       in the code below. These routines should be made available from 
//       a centralized place to all molfile plugins, eliminating code 
//       duplication.

/// return voxel, after safely clamping index to valid range
float edm_voxel_value_safe(int x, int y, int z, const int xsize, const int ysize, const int zsize, const float *data) {
  int xx, yy, zz; 
  xx = (x > 0) ? ((x < xsize) ? x : xsize-1) : 0;
  yy = (y > 0) ? ((y < ysize) ? y : ysize-1) : 0;
  zz = (z > 0) ? ((z < zsize) ? z : zsize-1) : 0;
  int index = zz*xsize*ysize + yy*xsize + xx;
  return data[index];
}

/// return interpolated value from 8 nearest neighbor voxels
float edm_voxel_value_interpolate(float xv, float yv, float zv, const int xsize, const int ysize, const int zsize, const float *data) {
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

  tmp = edm_voxel_value_safe(x, y, z, xsize, ysize, zsize, data);
  xlerps[0] = tmp + xf*(edm_voxel_value_safe(x+1, y, z, xsize, ysize, zsize, data) - tmp);

  tmp = edm_voxel_value_safe(x, y+1, z, xsize, ysize, zsize, data);
  xlerps[1] = tmp + xf*(edm_voxel_value_safe(x+1, y+1, z, xsize, ysize, zsize, data) - tmp);

  tmp = edm_voxel_value_safe(x, y, z+1, xsize, ysize, zsize, data);
  xlerps[2] = tmp + xf*(edm_voxel_value_safe(x+1, y, z+1, xsize, ysize, zsize, data) - tmp);

  tmp = edm_voxel_value_safe(x, y+1, z+1, xsize, ysize, zsize, data);
  xlerps[3] = tmp + xf*(edm_voxel_value_safe(x+1, y+1, z+1, xsize, ysize, zsize, data) - tmp);

  ylerps[0] = xlerps[0] + yf*(xlerps[1] - xlerps[0]);
  ylerps[1] = xlerps[2] + yf*(xlerps[3] - xlerps[2]);

  return ylerps[0] + zf*(ylerps[1] - ylerps[0]);
}

/// return interpolated value of voxel, based on atomic coords.
/// XXX need to account for non-orthog. cells
float edm_voxel_value_interpolate_from_coord(float xpos, float ypos, float zpos, const float *origin, const float *xdelta, const float *ydelta, const float *zdelta, const int xsize, const int ysize, const int zsize, float *data) {
  xpos = (xpos-origin[0])/xdelta[0];
  ypos = (ypos-origin[1])/ydelta[1];
  zpos = (zpos-origin[2])/zdelta[2];
  int gx = (int) xpos; // XXX this is wrong for non-orthog cells.
  int gy = (int) ypos;
  int gz = (int) zpos;

//  if (gx < 0 || gx >= xsize) return NAN;
//  if (gy < 0 || gy >= ysize) return NAN;
//  if (gz < 0 || gz >= zsize) return NAN;

  // Pad with zeroes
  if (gx < 0 || gx >= xsize) return 0;
  if (gy < 0 || gy >= ysize) return 0;
  if (gz < 0 || gz >= zsize) return 0;
    
  return edm_voxel_value_interpolate(xpos, ypos, zpos, xsize, ysize, zsize, data);

}

static int write_edm_data(void *v, molfile_volumetric_t *metadata, float *datablock, float *colorblock) {

  FILE *fd = (FILE *)v;
  const int xsize = metadata->xsize;
  const int ysize = metadata->ysize;
  const int zsize = metadata->zsize;

  float xaxis[3], yaxis[3], zaxis[3];
  float xdelta[3], ydelta[3], zdelta[3];
  float origin[3], porigin[3];

  int i, j, k;

  int na, amin, amax, nb, bmin, bmax, nc, cmin, cmax;
  float a, b, c, alpha, beta, gamma;

  for (i=0; i<3; i++) {
    origin[i] = metadata->origin[i];
    xaxis[i] = metadata->xaxis[i];
    yaxis[i] = metadata->yaxis[i];
    zaxis[i] = metadata->zaxis[i];
    xdelta[i] = xaxis[i]/(xsize - 1);
    ydelta[i] = yaxis[i]/(ysize - 1);
    zdelta[i] = zaxis[i]/(zsize - 1);
  }

  // The origin and the box length must be multiples of the grid spacing,
  // so we pad the box accordingly, keeping the same grid spacing and 
  // resampling the map.

  // For now, we will only write orthonogal cells
  if (fabs(xaxis[1]) > TOLERANCE || fabs(xaxis[2]) > TOLERANCE ||
      fabs(yaxis[0]) > TOLERANCE || fabs(yaxis[2]) > TOLERANCE ||
      fabs(zaxis[0]) > TOLERANCE || fabs(zaxis[1]) > TOLERANCE) {
    fprintf(stderr, "edmplugin) Could not write X-PLOR file: only orthogonal cells are currently supported.\n");
    return MOLFILE_ERROR;
  }

  amin = (int) floorf(origin[0]/xdelta[0]);
  bmin = (int) floorf(origin[1]/ydelta[1]);
  cmin = (int) floorf(origin[2]/zdelta[2]);

  porigin[0] = amin * xdelta[0];
  porigin[1] = bmin * ydelta[1];
  porigin[2] = cmin * zdelta[2];

  amax = (int) ceilf((xaxis[0]+origin[0])/xdelta[0]);
  bmax = (int) ceilf((yaxis[1]+origin[1])/ydelta[1]);
  cmax = (int) ceilf((zaxis[2]+origin[2])/zdelta[2]);

  na = amax - amin + 1;
  nb = bmax - bmin + 1;
  nc = cmax - cmin + 1;

  // The new cell axes may be slightly larger than the original ones
  a = na*xdelta[0];
  b = nb*ydelta[1];
  c = nc*zdelta[2];

  // Assuming the cell is orthognal
  alpha = beta = gamma = 90;

  // Write header
  fprintf(fd,"\n 2 !NTITLE\n"); // number of title lines
  fprintf(fd,"REMARKS FILENAME=\"\"\n");
  fprintf(fd,"REMARKS created by VMD\n");

  // Write the box dimensions and grid spacing deltas
  fprintf(fd,"%d %d %d %d %d %d %d %d %d\n", na, amin, amax, nb, bmin, bmax,
                                             nc, cmin, cmax);
  fprintf(fd,"%g %g %g %g %g %g\n", a, b, c, alpha, beta, gamma);

  // Write plane order
  fprintf(fd,"ZYX\n");

  // Copy voldata to a padded array
  int psize = na*nb*nc;
  int pxysize = na*nb;

  // Resample map
  float *pdata = (float*) malloc(psize*sizeof(float));
  for (i=0; i<na; i++) {
    float xpos = porigin[0] + i*xdelta[0];
    for (j=0; j<nb; j++) {
      float ypos = porigin[1] + j*ydelta[1];
      for (k=0; k<nc; k++) {
        float zpos = porigin[2] + k*zdelta[2];
        pdata[i + j*na + k*pxysize] = edm_voxel_value_interpolate_from_coord(xpos, ypos, zpos, origin, xdelta, ydelta, zdelta, xsize, ysize, zsize, datablock);
      }
    }
  }

  // Write each xy slice separately
  int count = 0;
  for (k=0; k<nc; k++) {
    if (count % 6) fprintf(fd, "\n");
    fprintf(fd, "%8d\n", k);
    count=0;
    for (j=0; j<nb; j++) {
      for (i=0; i<na; i++) {
        fprintf(fd, "%12.5e", pdata[k*pxysize + j*na + i]);
        if (++count % 6 == 0)
          fprintf(fd, "\n");
      }
    }
  }
  if (count % 6) fprintf(fd, "\n");

  // Write footer
  int sentinel = -9999;
  fprintf(fd, "%8d\n", sentinel);

  // Calculate average and standard deviation
  double avg = 0;
  double stddev = 0;
  double sum = 0;
  double sum2 = 0;
  for (i=0; i<psize; i++) {
    sum += pdata[i];
    sum2 += pdata[i]*pdata[i];
  }
  avg = sum/(double)psize;
  stddev = psize/(psize-1) * sqrt(sum2/psize - avg*avg);
  fprintf(fd, "%g %g\n", avg, stddev);

  free(pdata);

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
  plugin.name = "edm";
  plugin.prettyname = "XPLOR Electron Density Map";
  plugin.author = "John Stone, Leonardo Trabuco";
  plugin.majorv = 0;
  plugin.minorv = 8;
  plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
  plugin.filename_extension = "cns,edm,xplor";
  plugin.open_file_read = open_edm_read;
  plugin.read_volumetric_metadata = read_edm_metadata;
  plugin.read_volumetric_data = read_edm_data;
  plugin.close_file_read = close_edm_read;
#if vmdplugin_ABIVERSION > 9
  plugin.open_file_write = open_edm_write;
  plugin.write_volumetric_data = write_edm_data;
  plugin.close_file_write = close_edm_write;
#endif
  return VMDPLUGIN_SUCCESS; 
}

VMDPLUGIN_EXTERN int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_EXTERN int VMDPLUGIN_fini(void) { return VMDPLUGIN_SUCCESS; }



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




