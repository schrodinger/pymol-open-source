/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_namdbinplugin
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
 *      $RCSfile: namdbinplugin.c,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.19 $       $Date: 2006/02/23 19:36:45 $
 *
 ***************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>

#include "molfile_plugin.h"

#if INT_MAX == 2147483647
  typedef int namdbin_int32;
#elif SHRT_MAX == 2147483647
  typedef short namdbin_int32;
#elif LONG_MAX == 2147483647
  typedef long namdbin_int32;
#endif

typedef struct {
  FILE *fd;
  int numatoms;
  int wrongendian;
  double *xyz;
} namdbinhandle;

static void *open_namdbin_read(const char *path, const char *filetype, 
    int *natoms) {
  namdbinhandle *namdbin;
  FILE *fd;
  int numatoms;
  namdbin_int32 filen;
  char lenbuf[4];
  char tmpc;

  fd = fopen(path, "rb");
  if (!fd) {
    fprintf(stderr, "Could not open file '%s' for reading.\n", path);
    return NULL;
  }
  namdbin = (namdbinhandle *)malloc(sizeof(namdbinhandle));
  memset(namdbin, 0, sizeof(namdbinhandle));
  fseek(fd,0,SEEK_END);
  numatoms = (ftell(fd)-4)/24;
  if (numatoms < 1) {
    fprintf(stderr, "File '%s' is too short.\n", path);
    fclose(fd);
    free(namdbin);
    return NULL;
  }
  fseek(fd,0,SEEK_SET);
  fread(&filen, sizeof(namdbin_int32), 1, fd);
  if (filen != numatoms) {
    namdbin->wrongendian = 1;
    memcpy(lenbuf, (const char *)&filen, 4);
    tmpc = lenbuf[0]; lenbuf[0] = lenbuf[3]; lenbuf[3] = tmpc;
    tmpc = lenbuf[1]; lenbuf[1] = lenbuf[2]; lenbuf[2] = tmpc;
    memcpy((char *)&filen, lenbuf, 4);
  }
  if (filen != numatoms) {
    fprintf(stderr, "Inconsistent atom count in file '%s'.\n", path);
    fclose(fd);
    free(namdbin);
    return NULL;
  }
  if ( namdbin->wrongendian ) {
    fprintf(stderr, "File '%s' appears to be other-endian.\n", path);
  }
  namdbin->fd = fd;
  namdbin->numatoms = numatoms;
  namdbin->xyz = (double *)malloc(3 * namdbin->numatoms * sizeof(double));
  if (!namdbin->xyz) {
    fprintf(stderr, "Unable to allocate space for %d atoms.\n", namdbin->numatoms);
    fclose(fd);
    free(namdbin);
    return NULL;
  }
  *natoms = namdbin->numatoms;
  return namdbin;
}

static int read_next_timestep(void *v, int natoms, molfile_timestep_t *ts) {
  namdbinhandle *namdbin;
  int i, numatoms;
  char *cdata;
  char tmp0, tmp1, tmp2, tmp3;

  namdbin = (namdbinhandle *)v;
  if (!namdbin->fd) 
    return MOLFILE_ERROR;  /* Done reading frames */

  numatoms = namdbin->numatoms;

  if (fread(namdbin->xyz, sizeof(double), 3 * numatoms, namdbin->fd)
	                         != (size_t)(3 * numatoms)) {
    fprintf(stderr, "Failure reading data from NAMD binary file.\n");
    return MOLFILE_ERROR;
  }

  if (namdbin->wrongendian) {
    fprintf(stderr, "Converting other-endian data from NAMD binary file.\n");
    cdata = (char *) namdbin->xyz;
    for ( i=0; i<3*numatoms; ++i, cdata+=8 ) {
      tmp0 = cdata[0]; tmp1 = cdata[1];
      tmp2 = cdata[2]; tmp3 = cdata[3];
      cdata[0] = cdata[7]; cdata[1] = cdata[6];
      cdata[2] = cdata[5]; cdata[3] = cdata[4];
      cdata[7] = tmp0; cdata[6] = tmp1;
      cdata[5] = tmp2; cdata[4] = tmp3;
    }
  }

  if (ts) {
    for ( i=0; i<numatoms; ++i) {
      ts->coords[3*i] = namdbin->xyz[3*i];
      ts->coords[3*i+1] = namdbin->xyz[3*i+1];
      ts->coords[3*i+2] = namdbin->xyz[3*i+2];
    }
  }
  /*
   * Close the file handle and set to NULL so we know we're done reading 
   */
  fclose(namdbin->fd);
  namdbin->fd = NULL;

  return MOLFILE_SUCCESS;
}
 
static void close_file_read(void *v) {
  namdbinhandle *namdbin = (namdbinhandle *)v;
  if (namdbin->fd)
    fclose(namdbin->fd);
  free(namdbin->xyz);
  free(namdbin);
}

static void *open_namdbin_write(const char *path, const char *filetype, 
    int natoms) {
  namdbinhandle *namdbin;
  FILE *fd;

  fd = fopen(path, "wb");
  if (!fd) {
    fprintf(stderr, "Could not open file %s for writing\n", path);
    return NULL;
  }

  namdbin = (namdbinhandle *)malloc(sizeof(namdbinhandle));
  namdbin->fd = fd;
  namdbin->numatoms = natoms;
  return namdbin;
}

static int write_timestep(void *v, const molfile_timestep_t *ts) {
  
  int i;
  namdbin_int32 myint;
  namdbinhandle *namdbin = (namdbinhandle *)v;
  
  if (!namdbin->fd)
    return MOLFILE_ERROR;
  
  myint = (namdbin_int32)namdbin->numatoms;
  fwrite(&myint, 4, 1, namdbin->fd);

  for (i=0; i<3*namdbin->numatoms; i++) {
    double tmp = ts->coords[i];
    if (fwrite(&tmp, sizeof(double), 1, namdbin->fd) != 1) {
      fprintf(stderr, "Error writing namd binary file\n");
      return MOLFILE_ERROR;
    }
  }
  
  /*
   * Close and NULLify the file handle so we don't write any more frames.
   */
  fclose(namdbin->fd);
  namdbin->fd = NULL;

  return MOLFILE_SUCCESS;
}
       
static void close_file_write(void *v) {
  namdbinhandle *namdbin = (namdbinhandle *)v;
  if (namdbin->fd)
    fclose(namdbin->fd);
  free(namdbin);
}

/*
 * Initialization stuff here
 */

static molfile_plugin_t plugin = {
  vmdplugin_ABIVERSION,         /* ABI verison */
  MOLFILE_PLUGIN_TYPE,          /* type of plugin */
  "namdbin",			/* name of plugin */
  "NAMD Binary Coordinates",    /* name of plugin */
  "James Phillips, Justin Gullingsrud",	 /* author */
  0,				/* major version */
  1,				/* minor version */
  VMDPLUGIN_THREADSAFE,         /* is reentrant */
  "coor",                       /* filename extension */
  open_namdbin_read,
  0,
  0,
  read_next_timestep,
  close_file_read,
  open_namdbin_write,
  0,
  write_timestep,
  close_file_write
};

VMDPLUGIN_API int VMDPLUGIN_init() {
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_fini() {
  return VMDPLUGIN_SUCCESS;
}

  
#ifdef TEST_NAMDBINPLUGIN

int main(int argc, char *argv[]) {
  molfile_header_t header;
  molfile_timestep_t timestep;
  void *v;
  int i;

  while (--argc) {
    ++argv; 
    v = open_namdbin_read(*argv, &header);
    if (!v) {
      fprintf(stderr, "open_namdbin_read failed for file %s\n", *argv);
      return 1;
    }
    timestep.coords = (float *)malloc(3*sizeof(float)*header.numatoms);
    for (i=0; i<header.numsteps; i++) {
      int rc = read_next_timestep(v, &timestep);
      if (rc) {
        fprintf(stderr, "error in read_next_timestep\n");
        return 1;
      }
    }
    close_file_read(v);
  }
  return 0;
}
 
      
#endif  

