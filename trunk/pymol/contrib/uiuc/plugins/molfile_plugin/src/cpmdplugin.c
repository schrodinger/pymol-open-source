/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_cpmdplugin
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
 *      $RCSfile: cpmdplugin.c,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.11 $       $Date: 2006/02/23 19:36:44 $
 *
 ***************************************************************************/

#include "largefiles.h"   /* platform dependent 64-bit file I/O defines */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "molfile_plugin.h"

typedef struct {
  FILE *file;
  int numatoms;
  const char *file_name;
} cpmddata;
 
static void *open_cpmd_read(const char *filename, const char *filetype, 
                           int *natoms) {
  FILE *fd;
  cpmddata *data;
  char linebuf[255];
  int i, nfi_first, nfi_current, atomcount;
 
  printf("cpmd) trying to open file '%s'\n",filename);
	  
  fd = fopen(filename, "rb");
  if (!fd) return NULL;
  
  data = (cpmddata *)malloc(sizeof(cpmddata));
  data->file = fd;
  data->file_name = filename;

  nfi_first = 0;
  nfi_current = 0;

  /* first column is the current timestep number  */
  fgets(linebuf, 255, fd);
  i = sscanf(linebuf, "%d", &nfi_first);
  if (i < 1) {
    fprintf(stderr, "read) cpmd trajectory file '%s' should have the timestep number "
            "in the first column\n", filename);
    return NULL;
  }
  atomcount   = 0;
  nfi_current = nfi_first;

  /* loop through file until the contents of the first column
     changes, indicating the end of a configuration */
  while ((nfi_first == nfi_current) && !ferror(fd) && !feof(fd)) {
    ++atomcount;
    fgets(linebuf, 255, fd);
    i = sscanf(linebuf, "%d", &nfi_current);
    if (i < 1) {
      fprintf(stderr, "read) cpmd trajectory file '%s' should have the "
              "timestep number in the first column\n", filename);
      return NULL;
    }
  }
  printf("cpmd) found %d atoms in first timestep\n",atomcount);
  *natoms = atomcount;
  data->numatoms=*natoms;

  /* rewind to the beginning for reading the coordinates elsewhere.*/
  rewind(fd);

  return data;
}

static int read_cpmd_structure(void *mydata, int *optflags, 
                              molfile_atom_t *atoms) {
  int i, j;
  char *k;
  float coord;
  molfile_atom_t *atom;
  cpmddata *data = (cpmddata *)mydata;
  
  printf("cpmd) trying to read structure\n");
  *optflags = MOLFILE_NOOPTIONS; /* no optional data */

  for(i=0;i<data->numatoms;i++) {
    char buffer[1024];
    char fbuffer[1024];
    k = fgets(fbuffer, 1024, data->file);
    atom = atoms + i;
    j=sscanf(fbuffer, "%s %f %f %f", buffer, &coord, &coord, &coord);
    if (k == NULL) {
      fprintf(stderr, "cpmd structure) missing atom(s) in file '%s'\n",data->file_name);
      fprintf(stderr, "cpmd structure) expecting '%d' atoms, found only '%d'\n",data->numatoms,i+1);
      return MOLFILE_ERROR;
    } else if (j < 4) {
      fprintf(stderr, "cpmd structure) missing type or coordinate(s) in file '%s' for atom '%d'\n",data->file_name,i+1);
      return MOLFILE_ERROR;
    }

    /* I don't know what to do with these */
    strncpy(atom->name, buffer, sizeof(atom->name));
    strncpy(atom->type, buffer, sizeof(atom->type));
    atom->resname[0] = '\0';
    atom->resid = 1;
    atom->chain[0] = '\0';
    atom->segid[0] = '\0';
    /* skip to the end of line */
  }

  rewind(data->file);
  return MOLFILE_SUCCESS;
}

static int read_cpmd_timestep(void *mydata, int natoms, molfile_timestep_t *ts) {
  int i, j, nfi_first, nfi_current;
  char fbuffer[1024];
  float x, y, z;
  const float bohr=0.529177249;
  char *k;
  
  cpmddata *data = (cpmddata *)mydata;
  nfi_first = nfi_current = -1;
  
  /* read the coordinates */
  for (i=0; i<natoms; i++) {

    k = fgets(fbuffer, 1024, data->file);

    /* read next line if this is a continuation indicator */
    if (strstr(fbuffer, "NEW DATA")) {
      k = fgets(fbuffer, 1024, data->file);
    }
    j = sscanf(fbuffer, "%d %f %f %f", &nfi_current, &x, &y, &z);
    if (nfi_first < 0) nfi_first = nfi_current;
    
    if (k == NULL) {
      return MOLFILE_ERROR;
    } else if (j < 4) {
      fprintf(stderr, "cpmd timestep) missing or illegal data in file"
              " '%s' for atom '%d'\n",data->file_name,i+1);
      return MOLFILE_ERROR;
    } else if (nfi_first != nfi_current) {
      fprintf(stderr, "cpmd timestep) short record in timestep %d of file"
              " '%s' for atom '%d'\n",nfi_first, data->file_name,i+1);
    }
    
    ts->coords[3*i  ] = x*bohr;
    ts->coords[3*i+1] = y*bohr;
    ts->coords[3*i+2] = z*bohr;
  }

  return MOLFILE_SUCCESS;
}
    
static void close_cpmd_read(void *mydata) {
  cpmddata *data = (cpmddata *)mydata;
  
  fclose(data->file);
  free(data);
}


/* registration stuff */
static molfile_plugin_t cpmdplugin = {
  vmdplugin_ABIVERSION,
  MOLFILE_PLUGIN_TYPE,                         /* type */
  "cpmd",                                      /* short name */
  "CPMD",                                      /* pretty name */
  "Axel Kohlmeyer, John E. Stone",             /* author */
  0,                                           /* major version */
  3,                                           /* minor version */
  VMDPLUGIN_THREADSAFE,                        /* is reentrant */
  "cpmd",
  open_cpmd_read,
/*  read_cpmd_structure */
  0,
  0,
  read_cpmd_timestep,
  close_cpmd_read,
/*  open_cpmd_write, */
  0,
  0,
/*  write_cpmd_timestep, */
  0,
/*  close_cpmd_write */
  0
};

VMDPLUGIN_API int VMDPLUGIN_init() {
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)&cpmdplugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_fini() {
  return VMDPLUGIN_SUCCESS;
}

