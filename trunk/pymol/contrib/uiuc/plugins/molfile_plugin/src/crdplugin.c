/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_crdplugin
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
 *      $RCSfile: crdplugin.c,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.35 $       $Date: 2008/08/05 20:12:08 $
 *
 ***************************************************************************/

/*
 * TODO: This plugin should probably be merged with the 'rst7' plugin, since
 *       the differences between them are minor, and there's no logical reason
 *       for them to be implemented completely independently as they are now.
 *       The major differences in formatting are in regard to the 6F12.7 (rst7)
 *       versus 10F8.3 (crd) ascii floating point conversion modes. 
 */

#include "largefiles.h"   /* platform dependent 64-bit file I/O defines */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "molfile_plugin.h"

typedef struct {
  FILE *file;
  int has_box;
  int numatoms;
} crddata;
 
static void *open_crd_read(const char *filename, const char *filetype, 
    int *natoms) {
 
  FILE *fd;
  crddata *data;
 
  fd = fopen(filename, "rb");
  if (!fd) return NULL;
  
  /* first line is title, so skip past it */
  while (getc(fd) != '\n');

  /* 
   * CRD's don't store the number of atoms in the timestep, so we assume that
   * the application will determine this for us.  
   */
  data = (crddata *)malloc(sizeof(crddata));
  data->file = fd;
  *natoms = MOLFILE_NUMATOMS_UNKNOWN;
  /* filetype "crd" has no box; filetype "crdbox" does. */
  data->has_box = strcmp(filetype, "crd"); 
  return data;
}

/*
 * CRD files with box info are indistinguishable from regular CRD's.  
 * We regard CRD's with box info as a different file format.
 * CRD's don't tell how many atoms there are in each frame.  We therefore
 * rely on the numatoms field in the molfile_timestep_t parameter.
 */
static int read_crd_timestep(void *mydata, int natoms, molfile_timestep_t *ts) {
  crddata *crd = (crddata *)mydata;
  int i, j;
  float x, y, z;
  float a, b, c;

  /* Read in the atom coordinates */
  for (i=0; i<natoms; i++) {
    j = fscanf(crd->file, "%f %f %f", &x, &y, &z);
    if (j == EOF) {
      return MOLFILE_ERROR;
    } else if (j <= 0) {
      fprintf(stderr, "Problem reading CRD file\n");
      return MOLFILE_ERROR;
    }

    /* only save coords if we're given a valid ts pointer */ 
    /* otherwise assume that VMD wants us to skip it.     */
    if (ts != NULL) {
      ts->coords[3*i  ] = x;
      ts->coords[3*i+1] = y;
      ts->coords[3*i+2] = z;
    }
  }


  /* Read the PBC box info. */
  if (crd->has_box) {
    j = fscanf(crd->file, "%f %f %f", &a, &b, &c);
    if (j == EOF) {
      printf("EOF in box\n");
      return MOLFILE_ERROR;
    } else if (j <= 0) {
      printf("Problem reading box part of CRD file, scanf returned %d\n",j);
      return MOLFILE_ERROR;
    }

    /* only save coords if we're given a valid ts pointer */ 
    /* otherwise assume that VMD wants us to skip it.     */
    if (ts != NULL) {
      ts->A = a;
      ts->B = b;
      ts->C = c;

      /* XXX periodic cell angles are only stored in the PARM file */
      /* we should probably retrieve these from the already-loaded */
      /* molecule when possible.                                   */
      ts->alpha = 90.0;
      ts->beta  = 90.0;
      ts->gamma = 90.0;
    }
  }

  return MOLFILE_SUCCESS;
}
    
static void close_crd_read(void *mydata) {
  crddata *crd = (crddata *)mydata;
  fclose(crd->file);
  free(crd);
}

static void *open_crd_write(const char *path, const char *filetype,
    int natoms) {
  crddata *crd;
  FILE *fd;

  fd = fopen(path, "wb");
  if (!fd) {
    fprintf(stderr, "Could not open file %s for writing\n", path);
    return NULL;
  }
  fprintf(fd, "TITLE : Created by VMD with %d atoms\n", natoms);
  
  crd = (crddata *)malloc(sizeof(crddata));
  crd->file = fd;
  crd->numatoms = natoms;
  crd->has_box = strcmp(filetype, "crd"); 
  return crd;
}    
  
static int write_crd_timestep(void *v, const molfile_timestep_t *ts) {
  crddata *crd = (crddata *)v;
  int i, lfdone;
  const int ndata = crd->numatoms * 3;
  for (i=0; i<ndata; i++) {
    lfdone = 0;
    fprintf(crd->file, "%8.3f", ts->coords[i]);
    if ((i+1) % 10 == 0) {
      fprintf(crd->file, "\n"); 
      lfdone = 1;
    }
  }
  if (!lfdone)
    fprintf(crd->file, "\n"); 
    
  if (crd->has_box) {
    fprintf (crd->file, "%8.3f %8.3f %8.3f\n", ts->A, ts->B, ts->C);
  }

  return MOLFILE_SUCCESS;
}

static void close_crd_write(void *v) {
  crddata *crd = (crddata *)v;
  fclose(crd->file);
  free(crd);
}

/* registration stuff */
    
static molfile_plugin_t plugin;
static molfile_plugin_t crdboxplugin;

VMDPLUGIN_API int VMDPLUGIN_init(void) { 
  memset(&plugin, 0, sizeof(molfile_plugin_t));
  plugin.abiversion = vmdplugin_ABIVERSION;
  plugin.type = MOLFILE_PLUGIN_TYPE;
  plugin.name = "crd";
  plugin.prettyname = "AMBER Coordinates";
  plugin.author = "Justin Gullingsrud, John Stone";
  plugin.majorv = 0;
  plugin.minorv = 7;
  plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
  plugin.filename_extension = "mdcrd,crd";
  plugin.open_file_read = open_crd_read;
  plugin.read_next_timestep = read_crd_timestep;
  plugin.close_file_read = close_crd_read;
  plugin.open_file_write = open_crd_write;
  plugin.write_timestep = write_crd_timestep;
  plugin.close_file_write = close_crd_write;

  memcpy(&crdboxplugin, &plugin, sizeof(molfile_plugin_t));
  crdboxplugin.name = "crdbox";
  crdboxplugin.prettyname = "AMBER Coordinates with Periodic Box";

  return VMDPLUGIN_SUCCESS; 
}

VMDPLUGIN_API int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)(void *)&plugin);
  (*cb)(v, (vmdplugin_t *)(void *)&crdboxplugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_fini(void) { return VMDPLUGIN_SUCCESS; }

