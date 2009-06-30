/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_mmcif
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
 *      $RCSfile: mmcif.C,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.6 $       $Date: 2009/04/29 15:45:31 $
 *
 ***************************************************************************/

/*
 *  mmCIF molecule file format (a subset of the STAR file format):
 *    http://mmcif.rcsb.org/mmcif-early/background/index.html#1
 *    http://mmcif.rcsb.org/mmcif-early/workshop/mmCIF-tutorials/
 *
 *  mmCIF reserved keywords: data_ loop_ global_ save_ stop_
 *
 *  STAR file syntax and rules (superset of mmCIF):
 *    http://journals.iucr.org/iucr-top/cif/standard/cifstd4.html
 *
 *  STAR syntactic entities: 
 *    Text string: string of characters bounded by blanks, single quotes (') 
 *                 double quotes ("), or by semi-colons (;) as the first 
 *                 character of a line
 *
 *    Data name:   a text string starting with an underline (_) character 
 *
 *    Data item:   a text string not starting with an underline, but preceded 
 *                 by a data name to identify it 
 * 
 *    Data loop:   a list of data names, preceded by loop_ and followed by 
 *                 a repeated list of data items 
 *
 *    Data block:  a collection of data names (looped or not) and data items 
 *                 that are preceded by a data_ code record. A data name must 
 *                 be unique within a data block. A data block is terminated 
 *                 by another data_ statement or the end of file 
 * 
 *    Data file:   a collection of data blocks; the block codes must be 
 *                 unique within a data file  
 *
 *  CIF restrictions to STAR syntax:
 *    http://journals.iucr.org/iucr-top/cif/standard/cifstd5.html
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "molfile_plugin.h"
#include "periodic_table.h"

typedef struct {
  FILE *file;
  int numatoms;
  molfile_atom_t *atomlist;
} mmcifdata;
 
static void *open_mmcif_read(const char *filename, const char *filetype, 
                           int *natoms) {
  // XXX not finished yet
  return NULL;
}

static int read_mmcif_structure(void *mydata, int *optflags, 
                              molfile_atom_t *atoms) {
  // XXX not finished yet 
  return MOLFILE_ERROR;
}

static int read_mmcif_timestep(void *mydata, int natoms, molfile_timestep_t *ts) {
  mmcifdata *data = (mmcifdata *)mydata;

  // XXX not finished yet 
  return MOLFILE_ERROR;
}
    
static void close_mmcif_read(void *mydata) {
  mmcifdata *data = (mmcifdata *)mydata;
  fclose(data->file);
  free(data);
}


static void *open_mmcif_write(const char *filename, const char *filetype, 
                           int natoms) {
  FILE *fd;
  mmcifdata *data;

  fd = fopen(filename, "w");
  if (!fd) { 
    fprintf(stderr, "mmcifplugin) Error: unable to open mmcif file %s for writing\n",
            filename);
    return NULL;
  }
  
  data = (mmcifdata *)malloc(sizeof(mmcifdata));
  data->numatoms = natoms;
  data->file = fd;
  return data;
}

/* registration stuff */
static molfile_plugin_t plugin;

VMDPLUGIN_EXTERN int VMDPLUGIN_init() {
  memset(&plugin, 0, sizeof(molfile_plugin_t));
  plugin.abiversion = vmdplugin_ABIVERSION;
  plugin.type = MOLFILE_PLUGIN_TYPE;
  plugin.name = "mmcif";
  plugin.prettyname = "mmCIF";
  plugin.author = "John Stone";
  plugin.majorv = 0;
  plugin.minorv = 2;
  plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
  plugin.filename_extension = "cif";
  plugin.open_file_read = open_mmcif_read;
  plugin.read_structure = read_mmcif_structure;
  plugin.read_next_timestep = read_mmcif_timestep;
  plugin.close_file_read = close_mmcif_read;
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_EXTERN int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_EXTERN int VMDPLUGIN_fini() {
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
    v = open_mmcif_read(*argv, "mmcif", &natoms);
    if (!v) {
      fprintf(stderr, "open_mmcif_read failed for file %s\n", *argv);
      return 1;
    }
    fprintf(stderr, "open_mmcif_read succeeded for file %s\n", *argv);
    fprintf(stderr, "number of atoms: %d\n", natoms);

    i = 0;
    timestep.coords = (float *)malloc(3*sizeof(float)*natoms);
    while (!read_mmcif_timestep(v, natoms, &timestep)) {
      i++;
    }
    fprintf(stderr, "ended read_next_timestep on frame %d\n", i);

    close_mmcif_read(v);
  }
  return 0;
}

#endif





