/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_avsplugin
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
 *      $RCSfile: avsplugin.C,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.17 $       $Date: 2006/02/23 19:36:43 $
 *
 ***************************************************************************/

/*
 * AVS field files
 *
 * XXX - This plugin currently only supports the specific subset of AVS field
 * files that are produced by autodock. 'field' type must be 'uniform',
 * 'data' must be of type 'float', 'ndim' and 'nspace' must be 3, and
 * 'filetype' of all files referenced must be 'ascii'. 
 *
 * XXX - The plugin also expects the values to appear in a certain order,
 * this should definitely be fixed.
 *
 * More info for this format can be found at
 * <http://astronomy.swin.edu.au/~pbourke/geomformats/field/>
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>

#if defined(_AIX)
#include <strings.h>
#endif

#if defined(WIN32) || defined(WIN64)
#define strcasecmp  stricmp
#define strncasecmp strnicmp
#endif

#include "molfile_plugin.h"

#define LINESIZE 256

typedef struct {
  char filename[256]; /* XXX - LAME */
  int filetype, skip, offset, stride;
} datasource_t;

typedef struct {
  int nsets;
  molfile_volumetric_t *vol;
  datasource_t *data;
} avsfield_t;

enum {NONE, ASCII, BINARY, UNFORMATTED};  /* File types */
enum {UNIFORM, IRREGULAR, RECTILINEAR};   /* Field types */
enum {AVSFLOAT};                          /* Data types */

/* Reads lines from the stream into the array pointed to by s.
 * Returns a pointer to s when the first non-comment line is read or NULL on
 * error.
 */
static char *get_string(char *s, int n, FILE *stream) {
  do {
    if (fgets(s, n, stream) == NULL) {
      fprintf(stderr, "Error reading string.\n");
      return NULL;
    }
  } while (s[0] == '#');
  return s;
}

/* Read information from a string and store it into a datasource structure.
 * Returns 0 on success, 1 on error.
 */
static int read_datasource(char *s, datasource_t *data) {
  char *src, *tok, *value;
  src = strdup(s);
  tok = strtok(src, " \t\n");

  /* Load default values -- used if these attributes aren't set */
  data->skip = 0;
  data->offset = 0;
  data->stride = 1;

  /* Load default values -- must be changed */
  data->filename[0] = '\0';
  data->filetype = NONE;

  /* The first word should be "coord" or "variable" */
  if ( (strcasecmp(tok, "coord") != 0) && (strcasecmp(tok, "variable") != 0) ) {
    fprintf(stderr, "Improperly formatted header: expected coord or variable.\n");
    free(src);
    return 1;
  }

  /* Next should be the integer ID of the data source */
  tok = strtok(NULL, " \t\n");
  if (!isdigit(*tok)) {
    fprintf(stderr, "Improperly formatted header: expected ID.\n");
    free(src);
    return 1;
  }

  /* Now read the additional arguments */
  tok = strtok(NULL, " \t\n");
  while(tok) {
    value = strchr(tok, '=');
    if (!value) {
      fprintf(stderr, "Error reading value.\n");
      free(src);
      return 1;
    }
    value++; /* Point to the first character after '=' */

    if (strncasecmp(tok, "file=", value - tok) == 0) {
      /* XXX - This should be changed to something safer */
      strcpy(data->filename, value);
    }
    else if (strncasecmp(tok, "filetype=", value - tok) == 0) {
      /* XXX - For now, only ascii files are recognized. Other possible
       * values are "unformatted" for unformatted Fortran data, and "binary"
       * for raw binary data.
       */
      if (strcasecmp(value, "ascii") == 0) {
        data->filetype = ASCII;
      }
      else {
        fprintf(stderr, "Non-ASCII files are not supported.\n");
        free(src);
        return 1;
      }
    }
    else if (strncasecmp(tok, "skip=", value - tok) == 0) {
      /* XXX - This should probably be more rigorous */
      data->skip = atoi(value);
    }
    else if (strncasecmp(tok, "offset=", value - tok) == 0) {
      /* XXX - This should probably be more rigorous */
      data->offset = atoi(value);
    }
    else if (strncasecmp(tok, "stride=", value - tok) == 0) {
      /* XXX - This should definitely be more rigorous -- we don't want
       * stride set to 0 of the value isn't an integer. */
      data->stride = atoi(value);
    }
    else {
      /* XXX - For now, return with an error if there's an unrecognized
       * argument. This should probably be changed.
       */
      fprintf(stderr, "Unrecognized argument.\n");
      free(src);
      return 1;
    }

    tok = strtok(NULL, " \t\n");
  }

  free(src);

  /* Make sure the filename and filetype have been set */
  if ((data->filename[0] == '\0') || (data->filetype == NONE)) {
    fprintf(stderr, "Filename not set in options.\n");
    return 1;
  }
  
  return 0;
}

static void *open_avsfield_read(const char *filepath, const char *filetype, int *natoms) {
  avsfield_t *avsfield;
  FILE *fd;
  char inbuf[LINESIZE], current_file[256];
  int ndim, nspace, veclen, xsize, ysize, zsize, 
      index, i, coord_count, var_count;
  float value, origin[3], gridlength[3];
  datasource_t *coord = NULL, *variable = NULL;

  fd = fopen(filepath, "rb");
  if (!fd) {
    fprintf(stderr, "Error opening file.\n");
    return NULL;
  }

  /* Check for an AVS file */
  if (fgets(inbuf, LINESIZE, fd) == NULL) {
    fclose(fd);
    fprintf(stderr, "Error reading line.\n");
    return NULL;
  }
  if (strncmp(inbuf, "# AVS", 5) != 0) {
    fclose(fd);
    fprintf(stderr, "Improperly formatted header.\n");
    return NULL;
  }

  /* Check the number of dimensions */
  if (get_string(inbuf, LINESIZE, fd) == NULL) {
    fclose(fd);
    return NULL;
  }
  if (sscanf(inbuf, "ndim=%d", &ndim) != 1) {
    fprintf(stderr, "Error reading ndim.\n");
    fclose(fd);
    return NULL;
  }
  if (ndim != 3) {
    fprintf(stderr, "Error: ndim must be 3.\n");
    fclose(fd);
    return NULL;
  }

  /* Find the size of the grid in grid units */
  if (get_string(inbuf, LINESIZE, fd) == NULL) {
    fclose(fd);
    return NULL;
  }
  if (sscanf(inbuf, "dim1=%d", &xsize) != 1) {
    fprintf(stderr, "Error reading dim1.\n");
    fclose(fd);
    return NULL;
  }
  if (get_string(inbuf, LINESIZE, fd) == NULL) {
    fclose(fd);
    return NULL;
  }
  if (sscanf(inbuf, "dim2=%d", &ysize) != 1) {
    fprintf(stderr, "Error reading dim2.\n");
    fclose(fd);
    return NULL;
  }
  if (get_string(inbuf, LINESIZE, fd) == NULL) {
    fclose(fd);
    return NULL;
  }
  if (sscanf(inbuf, "dim3=%d", &zsize) != 1) {
    fprintf(stderr, "Error reading dim3.\n");
    fclose(fd);
    return NULL;
  }

  /* Check the number of coordinates per point */
  if (get_string(inbuf, LINESIZE, fd) == NULL) {
    fclose(fd);
    return NULL;
  }
  if (sscanf(inbuf, "nspace=%d", &nspace) != 1) {
    fprintf(stderr, "Error reading nspace.\n");
    fclose(fd);
    return NULL;
  }
  if (nspace != 3) {
    fprintf(stderr, "Error: nspace must be 3.\n");
    fclose(fd);
    return NULL;
  }

  /* Find out how many values are stored for each point (the length of the
   * vector for the vector field) */
  if (get_string(inbuf, LINESIZE, fd) == NULL) {
    fclose(fd);
    return NULL;
  }
  if (sscanf(inbuf, "veclen=%d", &veclen) != 1) {
    fprintf(stderr, "Error reading veclen.\n");
    fclose(fd);
    return NULL;
  }

  /* Check that the data type is "float" */
  if (get_string(inbuf, LINESIZE, fd) == NULL) {
    fclose(fd);
    return NULL;
  }
  if (strncmp(inbuf, "data=float", 10) != 0) {
    fprintf(stderr, "Error reading data type.\n");
    fclose(fd);
    return NULL;
  }

  /* Check that the field type is "uniform" */
  if (get_string(inbuf, LINESIZE, fd) == NULL) {
    fclose(fd);
    return NULL;
  }
  if (strncmp(inbuf, "field=uniform", 13) != 0) {
    fprintf(stderr, "Error reading field type.\n");
    fclose(fd);
    return NULL;
  }

  /* Allocate space for the coordinate and variable information 
   * coord is deleted in this fuction, variable is deleted when the plugin
   * is closed.
   */
  coord = new datasource_t[ndim];
  variable = new datasource_t[veclen];

  /* Find the coordinate information */
  for (i = 0; i < ndim; i++) {
    if (get_string(inbuf, LINESIZE, fd) == NULL) {
      delete[] coord;
      fclose(fd);
      return NULL;
    }
    if ( (sscanf(inbuf, "coord %d", &coord_count) != 1) || (coord_count != i+1) ) {
    fprintf(stderr, "Error reading coord count.\n");
      delete[] coord;
      fclose(fd);
      return NULL;
    }
    if (read_datasource(inbuf, &coord[i])) {
      delete[] coord;
      fclose(fd);
      return NULL;
    }
  }

  /* XXX - Ignore the labels (should be used for vol dataname) */
  for (i = 0; i < veclen; i++) {
    if (get_string(inbuf, LINESIZE, fd) == NULL) {
      delete[] coord;
      fclose(fd);
      return NULL;
    }
  }
  
  /* Find the variable information */
  for (i = 0; i < veclen; i++) {
    if (get_string(inbuf, LINESIZE, fd) == NULL) {
      delete[] coord;
      fclose(fd);
      return NULL;
    }
    if ( (sscanf(inbuf, "variable %d", &var_count) != 1) || (var_count != i+1) ) {
      fprintf(stderr, "Error reading variable count.\n");
      delete[] coord;
      fclose(fd);
      return NULL;
    }
    if (read_datasource(inbuf, &variable[i])) {
      delete[] coord;
      fclose(fd);
      return NULL;
    }
  }

  /* Close the AVS file */
  fclose(fd);
  fd = NULL;

  /* Read the coordinate file(s) to find the origin and grid size 
   * XXX - this only works for "uniform" fields
   */
  for (i = 0; i < ndim; i++) {
    if (strcmp(current_file, coord[i].filename) != 0) {
      /* Close the old file if one was open, and open a new one */
      if (fd) {
        fclose(fd);
        fd = NULL;
      }
      strcpy(current_file, coord[i].filename); /* XXX - unsafe */
      fd = fopen(current_file, "rb");
      if (!fd) {
        fprintf(stderr, "Error opening file.\n");
        delete[] coord;
        return NULL;
      }
    }
    else {
      /* Return to the beginning of the file */
      rewind(fd);
    }

    /* Skip the "skip" lines */
    for (index = 0; index < coord[i].skip; index++) {
      if (fgets(inbuf, LINESIZE, fd) == NULL) {
        fprintf(stderr, "Error reading line.\n");
        fclose(fd);
        delete[] coord;
        return NULL;
      }
    }

    /* Skip the "offset" values */
    for (index = 0; index < coord[i].offset; index++) {
      if (fscanf(fd, " %f", &value) != 1) {
        fprintf(stderr, "Error skipping offset.\n");
        fclose(fd);
        delete[] coord;
        return NULL;
      }
    }

    /* Read the origin, skip "stride" values, and read the end */
    if (fscanf(fd, " %f", &value) != 1) {
      fprintf(stderr, "Error reading origin.\n");
      fclose(fd);
      delete[] coord;
      return NULL;
    }
    origin[i] = value;
    for (index = 0; index < coord[i].stride; index++) {
      if (fscanf(fd, " %f", &value) != 1) {
        fprintf(stderr, "Error skipping stride.\n");
        fclose(fd);
        delete[] coord;
        return NULL;
      }
    }
    gridlength[i] = value - origin[i];
  }

  /* Free the coordinates */
  delete[] coord;
  coord = NULL;
  fclose(fd);

  /* Allocate and initialize the avsfield structure */
  avsfield = new avsfield_t;
  avsfield->vol = NULL;
  *natoms = MOLFILE_NUMATOMS_NONE;

  /* AVS field files can have an arbitrary number of sets */
  avsfield->nsets = veclen;
  avsfield->vol = new molfile_volumetric_t[veclen];
  avsfield->data = variable;

  for (i = 0; i < veclen; i++) {
    /* strcpy(avsfield->vol[i].dataname, "AVS Field: "); */
    sprintf(avsfield->vol[i].dataname, "AVS Field: %d", i);

    avsfield->vol[i].origin[0] = origin[0];
    avsfield->vol[i].origin[1] = origin[1];
    avsfield->vol[i].origin[2] = origin[2];

    avsfield->vol[i].xaxis[0] = gridlength[0];
    avsfield->vol[i].xaxis[1] = 0;
    avsfield->vol[i].xaxis[2] = 0;

    avsfield->vol[i].yaxis[0] = 0;
    avsfield->vol[i].yaxis[1] = gridlength[1];
    avsfield->vol[i].yaxis[2] = 0;

    avsfield->vol[i].zaxis[0] = 0;
    avsfield->vol[i].zaxis[1] = 0;
    avsfield->vol[i].zaxis[2] = gridlength[2];

    avsfield->vol[i].xsize = xsize;
    avsfield->vol[i].ysize = ysize;
    avsfield->vol[i].zsize = zsize;

    avsfield->vol[i].has_color = 0;
  }

  return avsfield;
}

static int read_avsfield_metadata(void *v, int *nsets,
  molfile_volumetric_t **metadata) {
  avsfield_t *avsfield = (avsfield_t *)v;
  *nsets = avsfield->nsets;
  *metadata = avsfield->vol;

  return MOLFILE_SUCCESS;
}

static int read_avsfield_data(void *v, int set, float *datablock,
                         float *colorblock) {
  avsfield_t *avsfield = (avsfield_t *)v;
  int skip, offset, stride, count, ndata, index;
  float value, *cellIndex = datablock;
  char inbuf[LINESIZE];
  FILE *fd;
  
  fd = fopen(avsfield->data[set].filename, "rb");
  if (!fd) {
    fprintf(stderr, "Error opening file.\n");
    return MOLFILE_ERROR; 
  }

  skip = avsfield->data[set].skip;
  offset = avsfield->data[set].offset;
  stride = avsfield->data[set].stride;

  count = 0;
  ndata = avsfield->vol[0].xsize * avsfield->vol[0].ysize * avsfield->vol[0].zsize;

  /* Skip the "skip" lines */
  for (index = 0; index < skip; index++) {
    if (fgets(inbuf, LINESIZE, fd) == NULL) {
      fprintf(stderr, "Error skipping lines.\n");
      fclose(fd);
      return MOLFILE_ERROR;
    }
  }

  /* Skip the "offset" values */
  for (index = 0; index < offset; index++) {
    if (fscanf(fd, " %f", &value) != 1) {
      fprintf(stderr, "Error skipping offset.\n");
      fclose(fd);
      return MOLFILE_ERROR;
    }
  }

  while (count < ndata) {
    /* Read a value into the datablock and skip "stride" values */
    if (fscanf(fd, " %f", &value) != 1) {
      fprintf(stderr, "Error reading data.\n");
      fclose(fd);
      return MOLFILE_ERROR;
    }
    *cellIndex = value;
    cellIndex++;
    count++;

    for (index = 0; index < stride-1; index++) {
      if (fscanf(fd, " %f", &value) != 1) {
        fprintf(stderr, "Error skipping stride.\n");
        fclose(fd);
        return MOLFILE_ERROR;
      }
    }
  }

  fclose(fd);
  return MOLFILE_SUCCESS;
}

static void close_avsfield_read(void *v) {
  avsfield_t *avsfield = (avsfield_t *)v;

  if (avsfield->vol != NULL)
    delete [] avsfield->vol;
  if (avsfield->data != NULL)
    delete [] avsfield->data;
  delete avsfield;
}

/*
 * Initialization stuff here
 */
static molfile_plugin_t plugin = {
  vmdplugin_ABIVERSION,   /* ABI version */
  MOLFILE_PLUGIN_TYPE, 	  /* plugin type */
  "fld",                  /* short file format description */
  "AVS Field",            /* pretty file format description */
  "Eamon Caddigan",       /* author(s) */
  0,                      /* major version */
  3,                      /* minor version */
  VMDPLUGIN_THREADSAFE,   /* is reentrant */
  "fld"                   /* filename extension */
};

VMDPLUGIN_EXTERN int VMDPLUGIN_init(void) { return VMDPLUGIN_SUCCESS; }
VMDPLUGIN_EXTERN int VMDPLUGIN_fini(void) { return VMDPLUGIN_SUCCESS; }
VMDPLUGIN_EXTERN int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  plugin.open_file_read = open_avsfield_read;
  plugin.read_volumetric_metadata = read_avsfield_metadata;
  plugin.read_volumetric_data = read_avsfield_data;
  plugin.close_file_read = close_avsfield_read;
  (*cb)(v, (vmdplugin_t *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

