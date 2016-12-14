/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_namdbinplugin
#define STATIC_PLUGIN 1

/***************************************************************************
 *cr
 *cr            (C) Copyright 1995-2016 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: namdbinplugin.c,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.23 $       $Date: 2016/11/28 05:01:54 $
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

#define BLOCK 500

typedef struct {
  double xyz[3*BLOCK];
  FILE *fd;
  int numatoms;
  int wrongendian;
} namdbinhandle;

static void *open_namdbin_read(const char *path, const char *filetype, 
    int *natoms) {
  namdbinhandle *namdbin;
  FILE *fd;
  int numatoms;
  namdbin_int32 filen;
  char lenbuf[4];
  char tmpc;

  namdbin = (namdbinhandle *)malloc(sizeof(namdbinhandle));
  if (!namdbin) {
    fprintf(stderr, "Unable to allocate space for read buffer.\n");
    return NULL;
  }
  memset(namdbin, 0, sizeof(namdbinhandle));

  fd = fopen(path, "rb");
  if (!fd) {
    fprintf(stderr, "Could not open file '%s' for reading.\n", path);
    free(namdbin);
    return NULL;
  }
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
  *natoms = namdbin->numatoms;
  return namdbin;
}

static int read_next_timestep(void *v, int natoms, molfile_timestep_t *ts) {
  namdbinhandle *namdbin;
  int i, numatoms;
  char *cdata;

  namdbin = (namdbinhandle *)v;
  if (!namdbin->fd) 
    return MOLFILE_ERROR;  /* Done reading frames */

  numatoms = namdbin->numatoms;

  for (i=0; i<namdbin->numatoms; i+=BLOCK) {
    int j, n;
    n = namdbin->numatoms - i;
    if ( n > BLOCK ) n = BLOCK;

    if (fread(namdbin->xyz, sizeof(double), 3*n, namdbin->fd) != (size_t)(3*n)) {
      fprintf(stderr, "Failure reading data from NAMD binary file.\n");
      return MOLFILE_ERROR;
    }

    if (namdbin->wrongendian) {
      if ( ! i ) fprintf(stderr, "Converting other-endian data from NAMD binary file.\n");
      cdata = (char *) namdbin->xyz;
      for ( j=0; j<3*n; ++j, cdata+=8 ) {
        char tmp0, tmp1, tmp2, tmp3;
        tmp0 = cdata[0]; tmp1 = cdata[1];
        tmp2 = cdata[2]; tmp3 = cdata[3];
        cdata[0] = cdata[7]; cdata[1] = cdata[6];
        cdata[2] = cdata[5]; cdata[3] = cdata[4];
        cdata[7] = tmp0; cdata[6] = tmp1;
        cdata[5] = tmp2; cdata[4] = tmp3;
      }
    }

    if (ts) {
      for (j=0; j<n; ++j) {
        ts->coords[3L*(i+j)  ] = namdbin->xyz[3*j  ];
        ts->coords[3L*(i+j)+1] = namdbin->xyz[3*j+1];
        ts->coords[3L*(i+j)+2] = namdbin->xyz[3*j+2];
      }
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
  free(namdbin);
}

static void *open_namdbin_write(const char *path, const char *filetype, 
    int natoms) {
  namdbinhandle *namdbin;
  FILE *fd;

  namdbin = (namdbinhandle *)malloc(sizeof(namdbinhandle));
  if (!namdbin) {
    fprintf(stderr, "Unable to allocate space for write buffer.\n");
    return NULL;
  }

  fd = fopen(path, "wb");
  if (!fd) {
    fprintf(stderr, "Could not open file %s for writing\n", path);
    free(namdbin);
    return NULL;
  }

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

  for (i=0; i<namdbin->numatoms; i+=BLOCK) {
    double *tmp = namdbin->xyz;
    int j, n;
    n = namdbin->numatoms - i;
    if ( n > BLOCK ) n = BLOCK;
    for (j=0; j<n; ++j) {
      tmp[3*j  ] = ts->coords[3L*(i+j)  ];
      tmp[3*j+1] = ts->coords[3L*(i+j)+1];
      tmp[3*j+2] = ts->coords[3L*(i+j)+2];
    }
    if (fwrite(tmp, sizeof(double), 3*n, namdbin->fd) != 3*n) {
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

static molfile_plugin_t plugin;

VMDPLUGIN_API int VMDPLUGIN_init() {
  memset(&plugin, 0, sizeof(molfile_plugin_t));
  plugin.abiversion = vmdplugin_ABIVERSION;
  plugin.type = MOLFILE_PLUGIN_TYPE;
  plugin.name = "namdbin";
  plugin.prettyname = "NAMD Binary Coordinates";
  plugin.author = "James Phillips, Justin Gullingsrud";
  plugin.majorv = 0;
  plugin.minorv = 2;
  plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
  plugin.filename_extension = "coor";
  plugin.open_file_read = open_namdbin_read;
  plugin.read_next_timestep = read_next_timestep;
  plugin.close_file_read = close_file_read;
  plugin.open_file_write = open_namdbin_write;
  plugin.write_timestep = write_timestep;
  plugin.close_file_write = close_file_write;
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

