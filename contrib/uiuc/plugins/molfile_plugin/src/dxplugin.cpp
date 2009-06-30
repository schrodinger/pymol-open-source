/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_dxplugin
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
 *      $RCSfile: dxplugin.C,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.40 $       $Date: 2009/04/29 15:45:29 $
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
 * object 3 class array type double rank 0 items { xn*yn*zn } [binary] data follows
 * f1 f2 f3
 * f4 f5 f6 f7 f8 f9
 * .
 * .
 * .
 * object "Dataset name" class field
 
 * Where xn, yn, and zn are the number of data points along each axis;
 * xorg, yorg, and zorg is the origin of the grid, assumed to be in angstroms;
 * xdel, ydel, and zdel are the scaling factors to convert grid units to
 * angstroms.
 *
 * Grid data follows, with a single or multiple values per line (maximum 
 * allowed linelength is hardcoded into the plugin with ~2000 chars), 
 * ordered z fast, y medium, and x slow.
 *
 * Note that the ordering of grid data in VMD's VolumetricData class
 * is transposed, i.e. x changes fastest and z slowest! 
 * The file reading and writing routines perform the transpose.
 *
 * If the optional keyword 'binary' is present, the data is expected to
 * be in binary, native endian, single precision IEEE-754 floating point format.
 *
 * Note: this plugin was written to read the OpenDX files created by the
 * APBS program, and thus supports only files that are writting in this style.
 * the full OpenDX data format is extremely powerful, complex, and flexible.
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
#include "largefiles.h"

#define THISPLUGIN plugin
#include "vmdconio.h"

#define LINESIZE 2040

typedef struct {
  FILE *fd;
  int nsets;
  molfile_volumetric_t *vol;
  int isBinary; 
} dx_t;


// Get a string from a stream, printing any errors that occur
static char *dxgets(char *s, int n, FILE *stream) {
  char *returnVal;

  if (feof(stream)) {
    vmdcon_printf(VMDCON_ERROR, "dxplugin) Unexpected end-of-file.\n");
    return NULL;
  } else if (ferror(stream)) {
    vmdcon_printf(VMDCON_ERROR, "dxplugin) Error reading file.\n");
    return NULL;
  } else {
    returnVal = fgets(s, n, stream);
    if (returnVal == NULL) {
      vmdcon_printf(VMDCON_ERROR, "dxplugin) Error reading line.\n");
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
  int isBinary = 0;
  
  fd = fopen(filepath, "rb");
  if (!fd) {
    vmdcon_printf(VMDCON_ERROR, "dxplugin) Error opening file.\n");
    return NULL;
  }

  /* skip comments */
  do {
    if (dxgets(inbuf, LINESIZE, fd) == NULL) 
      return NULL;
  } while (inbuf[0] == '#');

  /* get the number of grid points */
  if (sscanf(inbuf, "object 1 class gridpositions counts %d %d %d", &xsize, &ysize, &zsize) != 3) {
    vmdcon_printf(VMDCON_ERROR, "dxplugin) Error reading grid dimensions.\n");
    return NULL;
  }

  /* get the cell origin */
  if (dxgets(inbuf, LINESIZE, fd) == NULL) {
    return NULL;
  }
  if (sscanf(inbuf, "origin %e %e %e", orig, orig+1, orig+2) != 3) {
    vmdcon_printf(VMDCON_ERROR, "dxplugin) Error reading grid origin.\n");
    return NULL;
  }

  /* get the cell dimensions */
  if (dxgets(inbuf, LINESIZE, fd) == NULL) {
    return NULL;
  }
  if (sscanf(inbuf, "delta %e %e %e", xdelta, xdelta+1, xdelta+2) != 3) {
    vmdcon_printf(VMDCON_ERROR, "dxplugin) Error reading cell x-dimension.\n");
    return NULL;
  }

  if (dxgets(inbuf, LINESIZE, fd) == NULL) {
    return NULL;
  }
  if (sscanf(inbuf, "delta %e %e %e", ydelta, ydelta+1, ydelta+2) != 3) {
    vmdcon_printf(VMDCON_ERROR, "dxplugin) Error reading cell y-dimension.\n");
    return NULL;
  }

  if (dxgets(inbuf, LINESIZE, fd) == NULL) {
    return NULL;
  }
  if (sscanf(inbuf, "delta %e %e %e", zdelta, zdelta+1, zdelta+2) != 3) {
    vmdcon_printf(VMDCON_ERROR, "dxplugin) Error reading cell z-dimension.\n");
    return NULL;
  }

  /* skip the next line of the header (described at the beginning of
   * the code), which aren't utilized by APBS-produced DX files.  */
  if (dxgets(inbuf, LINESIZE, fd) == NULL) 
    return NULL;
  /* The next line tells us whether to expect ascii or binary format.
   * We scan for the word 'binary' somewhere in the line, and if it's found
   * we assume binary.
   */
  if (dxgets(inbuf, LINESIZE, fd) == NULL)
    return NULL;
  if (strstr(inbuf, "binary")) {
      isBinary = 1;
  }
 
  /* allocate and initialize the dx structure */
  dx = new dx_t;
  dx->fd = fd;
  dx->vol = NULL;
  dx->isBinary = isBinary;
  *natoms = MOLFILE_NUMATOMS_NONE;
  dx->nsets = 1; /* this file contains only one data set */

  dx->vol = new molfile_volumetric_t[1];
  strcpy(dx->vol[0].dataname, "DX map");

  /* Set the unit cell origin and basis vectors */
  for (int i=0; i<3; i++) {
    dx->vol[0].origin[i] = orig[i];
    dx->vol[0].xaxis[i] = xdelta[i] * ((xsize-1 > 0) ? (xsize-1) : 1);
    dx->vol[0].yaxis[i] = ydelta[i] * ((ysize-1 > 0) ? (ysize-1) : 1);
    dx->vol[0].zaxis[i] = zdelta[i] * ((zsize-1 > 0) ? (zsize-1) : 1);
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

static int read_binary_dx_data(dx_t *dx, int set, float *datablock) {
    
  int i, j, k;
  int xsize = dx->vol[0].xsize;
  int ysize = dx->vol[0].ysize;
  int zsize = dx->vol[0].zsize;
  int xysize = xsize * ysize;
  size_t total = xysize * zsize;
  float *tmp = (float *)malloc(total*sizeof(float));
  if (fread(tmp, sizeof(float), total, dx->fd) != total) {
    vmdcon_printf(VMDCON_ERROR, "dxplugin) Failed to read %d binary floats\n", total);
    free(tmp);
    return MOLFILE_ERROR;
  }
  // take the transpose - nasty
  int ind = 0;
  for (i=0; i<xsize; i++) {
      for (j=0; j<ysize; j++) {
          for (k=0; k<zsize; k++) {
              datablock[k * xysize + j*xsize + i] = tmp[ind++];
          }
      }
  }
  free( tmp );
  return MOLFILE_SUCCESS;
}

static int read_dx_data(void *v, int set, float *datablock,
                         float *colorblock) {
  dx_t *dx = (dx_t *)v;
  FILE *fd = dx->fd;
  char inbuf[LINESIZE];
  char *p;
  float grid;
  int x, y, z, xsize, ysize, zsize, xysize, count, total, i, line;

  if (dx->isBinary)
    return read_binary_dx_data(dx, set, datablock);

  xsize = dx->vol[0].xsize;
  ysize = dx->vol[0].ysize;
  zsize = dx->vol[0].zsize;
  xysize = xsize * ysize;
  total = xysize * zsize;

  /* Read the values from the file */
  x = y = z = line = 0;
  for (count = 0; count < total;) {
    ++line;
    p=dxgets(inbuf, LINESIZE, fd);
    if (p == NULL) {
      vmdcon_printf(VMDCON_ERROR, "dxplugin) Error reading grid data.\n");
      vmdcon_printf(VMDCON_ERROR, "dxplugin) line: %d. item: %d/%d. last data: %s\n", 
              line, count, total, inbuf);
      return MOLFILE_ERROR;
    }

    // chop line into whitespace separated tokens and parse them one by one.
    while (*p != '\n' && *p != '\0') {

      // skip over whitespace and try to parse non-blank as number
      while (*p != '\0' && (*p == ' ' || *p == '\t' || *p == '\n')) ++p;
      i = sscanf(p, "%e", &grid);
      if (i < 0) break; // end of line/string. get a new one.
      
      // a 0 return value means non-parsable as number.
      if (i == 0) {
        vmdcon_printf(VMDCON_ERROR, "dxplugin) Error parsing grid data.\n");
        vmdcon_printf(VMDCON_ERROR, "dxplugin) line: %d. item: %d/%d. data %s\n", 
                line, count, total, p);
        return MOLFILE_ERROR;
      }

      // success! add to dataset.
      if (i == 1) {
        ++count;
        datablock[x + y*xsize + z*xysize] = grid;
        z++;
        if (z >= zsize) {
          z = 0; y++;
          if (y >= ysize) {
            y = 0; x++;
          }
        }
      }

      // skip over the parsed text and search for next blank.
      while (*p != '\0' && *p != ' ' && *p != '\t' && *p != '\n') ++p;
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

static void *open_dx_write(const char *filepath, const char *filetype,
        int natoms) {

    FILE *fd = fopen(filepath, "wb");
    if (!fd) {
        vmdcon_printf(VMDCON_ERROR, 
                      "dxplugin) Could not open path '%s' for writing.\n",
                      filepath);
    }
    return fd;
}

static void close_dx_write(void *v) {
    if (v) {
        fclose((FILE *)v);
    }
}

/*
 *
# Data from APBS
# 
# POTENTIAL (kT/e)
# 
object 1 class gridpositions counts 129 129 129
origin -3.075250e+01 -3.848600e+01 -2.908250e+01
delta 4.687500e-01 0.000000e+00 0.000000e+00
delta 0.000000e+00 6.250000e-01 0.000000e+00
delta 0.000000e+00 0.000000e+00 4.687500e-01
object 2 class gridconnections counts 129 129 129
object 3 class array type double rank 0 items 2146689 data follows

*/
static int write_dx_data(void *v, molfile_volumetric_t *metadata,
        float *datablock, float *colorblock) {

    int i, j, k, count;
    FILE *fd = (FILE *)v;
    const int xsize = metadata->xsize;
    const int ysize = metadata->ysize;
    const int zsize = metadata->zsize;
    const int xysize = xsize * ysize;
    const int total = xysize * zsize;

    double xdelta[3], ydelta[3], zdelta[3];
    for (i=0; i<3; i++) {
      xdelta[i] = metadata->xaxis[i]/(xsize - 1);
      ydelta[i] = metadata->yaxis[i]/(ysize - 1);
      zdelta[i] = metadata->zaxis[i]/(zsize - 1);
    }

    fprintf(fd, "# Data from VMD\n");
    fprintf(fd, "# %s\n", metadata->dataname);
    fprintf(fd, "object 1 class gridpositions counts %d %d %d\n",
            xsize, ysize, zsize);
    fprintf(fd, "origin %g %g %g\n", 
            metadata->origin[0], metadata->origin[1], metadata->origin[2]);
    fprintf(fd, "delta %g %g %g\n", 
            xdelta[0], xdelta[1], xdelta[2]);
    fprintf(fd, "delta %g %g %g\n", 
            ydelta[0], ydelta[1], ydelta[2]);
    fprintf(fd, "delta %g %g %g\n", 
            zdelta[0], zdelta[1], zdelta[2]);
    fprintf(fd, "object 2 class gridconnections counts %d %d %d\n",
            xsize, ysize, zsize);

    int useBinary = (getenv("VMDBINARYDX") != NULL);
    fprintf(fd, "object 3 class array type double rank 0 items %d %sdata follows\n", total, useBinary ? "binary " : "");
    count = 0;
    for (i=0; i<xsize; i++) {
        for (j=0; j<ysize; j++) {
            for (k=0; k<zsize; k++) {
                if (useBinary) {
                    fwrite(datablock + k*xysize + j*xsize + i, sizeof(float),
                            1, fd);
                } else {
                    fprintf(fd, "%g ", datablock[k*xysize + j*xsize + i]);
                    if (++count == 3) {
                        fprintf(fd, "\n");
                        count = 0;
                    }
                }
            }
        }
    }
    if (!useBinary && count) 
      fprintf(fd, "\n");

    // Replace any double quotes (") by single quotes (') in the 
    // dataname string to make sure that we don't prematurely
    // terminate the string in the dx file.
    char *squotes = new char[strlen(metadata->dataname)+1];
    strcpy(squotes, metadata->dataname);
    char *s = squotes;
    // while(s=strchr(s, '"')) *s = '\'';
    while(1) { 
      s=strchr(s, '"'); 
      if (s) { 
        *s = '\'';
      } else { 
        break;
      }
    }

    // Print last line
    fprintf(fd, "object \"%s\" class field\n", squotes);
    delete [] squotes;

    fflush(fd);
    return MOLFILE_SUCCESS;
}

/*
 * Initialization stuff here
 */

VMDPLUGIN_EXTERN int VMDPLUGIN_init(void) { 
  memset(&plugin, 0, sizeof(molfile_plugin_t));
  plugin.abiversion = vmdplugin_ABIVERSION;
  plugin.type = MOLFILE_PLUGIN_TYPE;
  plugin.name = "dx";
  plugin.prettyname = "DX";
  plugin.author = "Eamon Caddigan, Justin Gullingsrud, John Stone, Leonardo Trabuco";
  plugin.majorv = 1;
  plugin.minorv = 9;
  plugin.is_reentrant = VMDPLUGIN_THREADUNSAFE;
  plugin.filename_extension = "dx";
  plugin.open_file_read = open_dx_read;
  plugin.read_volumetric_metadata = read_dx_metadata;
  plugin.read_volumetric_data = read_dx_data;
  plugin.close_file_read = close_dx_read;
#if vmdplugin_ABIVERSION > 9
  plugin.open_file_write = open_dx_write;
  plugin.write_volumetric_data = write_dx_data;
  plugin.close_file_write = close_dx_write;
#endif
  return VMDPLUGIN_SUCCESS; 
}

VMDPLUGIN_EXTERN int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_EXTERN int VMDPLUGIN_fini(void) { return VMDPLUGIN_SUCCESS; }

