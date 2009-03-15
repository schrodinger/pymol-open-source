/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_tinkerplugin
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
 *      $RCSfile: tinkerplugin.c,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.12 $       $Date: 2007/02/25 20:23:44 $
 *
 ***************************************************************************/

 /*
  * The .arc file is the one coming straight out of tinker.
  * Many frames, looks like a regular tinker, except that in the last columns
  * after x, y and z it has :
  * 
  * atom type (not important for viz)
  * 
  * and then info in a Z-matrix form
  * 
  * atom to which you are bonded, atom to which you are 'angled' and atom to
  * which you are 'torsioned'.
  * 
  */

#include "largefiles.h"   /* platform dependent 64-bit file I/O defines */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "molfile_plugin.h"

typedef struct {
  FILE *file;
  int numatoms;
  char *file_name;
  molfile_atom_t *atomlist;
} tinkerdata;
 
static void *open_tinker_read(const char *filename, const char *filetype, 
                           int *natoms) {
  FILE *fd;
  tinkerdata *data;
  int i;

  fd = fopen(filename, "rb");
  if (!fd) return NULL;
  
  data = (tinkerdata *)malloc(sizeof(tinkerdata));
  data->file = fd;
  data->file_name = strdup(filename);

  /* First line is the number of atoms, followed by other stuff */
  i = fscanf(data->file, "%d", natoms);
  if (i < 1) {
    fprintf(stderr, "\n\nread) ERROR: tinker file '%s' should have the number of atoms in the first line.\n", filename);
    return NULL;
  }
  data->numatoms=*natoms;

  while (getc(fd) != '\n'); /* skip rest of this line */
 
  return data;
}

static int read_tinker_structure(void *mydata, int *optflags, 
                                 molfile_atom_t *atoms) {
  int i, j, atomid;
  char *k;
  float coord;
  molfile_atom_t *atom;
  tinkerdata *data = (tinkerdata *)mydata;

  *optflags = MOLFILE_NOOPTIONS;

  for (i=0; i<data->numatoms; i++) {
    char buffer[1024], fbuffer[1024];
    int atomtypeindex=0;
    k = fgets(fbuffer, 1024, data->file);
    atom = atoms + i;
    j=sscanf(fbuffer, "%d %s %f %f %f %d", &atomid, buffer, &coord, &coord, &coord, &atomtypeindex);
    if (k == NULL) {
      fprintf(stderr, "tinker structure) missing atom(s) in file '%s'\n", data->file_name);
      fprintf(stderr, "tinker structure) expecting '%d' atoms, found only '%d'\n", data->numatoms, i+1);
      return MOLFILE_ERROR;
    } else if (j < 5) {
      fprintf(stderr, "tinker structure) missing type or coordinate(s) in file '%s' for atom '%d'\n", data->file_name, i+1);
      return MOLFILE_ERROR;
    }

    strncpy(atom->name, buffer, sizeof(atom->name));
#if 0
    strncpy(atom->type, atom->name, sizeof(atom->type));
#else
    sprintf(atom->type, "%d", atomtypeindex);
#endif
    atom->resname[0] = '\0';
    atom->resid = 1;
    atom->chain[0] = '\0';
    atom->segid[0] = '\0';
  }

  rewind(data->file);
  return MOLFILE_SUCCESS;
}

static int read_tinker_timestep(void *mydata, int natoms, molfile_timestep_t *ts) {
  int i, j, atomid;
  char atom_name[1024], fbuffer[1024], *k;
  float x, y, z;
  
  tinkerdata *data = (tinkerdata *)mydata;
  
  /* skip over the first line */
  if (NULL == fgets(fbuffer, 1024, data->file))  return MOLFILE_ERROR;

  /* read the coordinates */
  for (i=0; i<natoms; i++) {
    k = fgets(fbuffer, 1024, data->file);

    /* Read in atom type, X, Y, Z, skipping any remaining data fields */
    j = sscanf(fbuffer, "%d %s %f %f %f", &atomid, atom_name, &x, &y, &z);
    if (k == NULL) {
      return MOLFILE_ERROR;
    } else if (j < 5) {
      fprintf(stderr, "tinker timestep) missing type or coordinate(s) in file '%s' for atom '%d'\n",data->file_name,i+1);
      return MOLFILE_ERROR;
    } else if (j >= 5) {
      if (ts != NULL) { 
        /* only save coords if we're given a timestep pointer, */
        /* otherwise assume that VMD wants us to skip past it. */
        ts->coords[3*i  ] = x;
        ts->coords[3*i+1] = y;
        ts->coords[3*i+2] = z;
      }
    } else {
      break;
    }
  }
  
  return MOLFILE_SUCCESS;
}
    
static void close_tinker_read(void *mydata) {
  tinkerdata *data = (tinkerdata *)mydata;
  fclose(data->file);
  free(data->file_name);
  free(data);
}

/* registration stuff */
static molfile_plugin_t plugin;

VMDPLUGIN_API int VMDPLUGIN_init() {
  memset(&plugin, 0, sizeof(molfile_plugin_t));
  plugin.abiversion = vmdplugin_ABIVERSION;
  plugin.type = MOLFILE_PLUGIN_TYPE;
  plugin.name = "tinker";
  plugin.prettyname = "Tinker";
  plugin.author = "John Stone";
  plugin.majorv = 0;
  plugin.minorv = 5;
  plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
  plugin.filename_extension = "arc";
  plugin.open_file_read = open_tinker_read;
  plugin.read_structure = read_tinker_structure;
  plugin.read_next_timestep = read_tinker_timestep;
  plugin.close_file_read = close_tinker_read;
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)(void *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_fini() {
  return VMDPLUGIN_SUCCESS;
}


#ifdef TEST_PLUGIN

int main(int argc, char *argv[]) {
  molfile_timestep_t timestep;
  void *v;
  int natoms;
  int i, nsets, set;

  while (--argc) {
    ++argv;
    v = open_tinker_read(*argv, "tinker", &natoms);
    if (!v) {
      fprintf(stderr, "open_tinker_read failed for file %s\n", *argv);
      return 1;
    }
    fprintf(stderr, "open_tinker_read succeeded for file %s\n", *argv);
    fprintf(stderr, "number of atoms: %d\n", natoms);

    i = 0;
    timestep.coords = (float *)malloc(3*sizeof(float)*natoms);
    while (!read_tinker_timestep(v, natoms, &timestep)) {
      i++;
    }
    fprintf(stderr, "ended read_next_timestep on frame %d\n", i);

    close_tinker_read(v);
  }
  return 0;
}

#endif

