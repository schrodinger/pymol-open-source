/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_binposplugin
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
 *      $RCSfile: binposplugin.c,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.10 $       $Date: 2006/05/01 22:34:20 $
 *
 ***************************************************************************/

#include "largefiles.h"   /* platform dependent 64-bit file I/O defines */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>

#include "molfile_plugin.h"

#if INT_MAX == 2147483647
  typedef int binpos_int32;
#elif SHRT_MAX == 2147483647
  typedef short binpos_int32;
#elif LONG_MAX == 2147483647
  typedef long binpos_int32;
#endif

typedef struct {
  FILE *fd;
  int numatoms;
  int wrongendian;
  float *xyz;
} binposhandle;

static void *open_binpos_read(const char *path, const char *filetype, 
    int *natoms) {
  binposhandle *binpos;
  FILE *fd;
  int er=0,point,igarb;
  char lenbuf[4];
  char tmpc;
  char magicchar[5];

  fd = fopen(path, "rb");
  if (!fd) 
   {
    fprintf(stderr, "Could not open file '%s' for reading.\n", path);
    return NULL;
   }
  binpos = (binposhandle *)malloc(sizeof(binposhandle));
  memset(binpos, 0, sizeof(binposhandle));
  fread(magicchar,sizeof(char),4,fd);
  magicchar[4]= '\0' ;
  if(strcmp(magicchar,"fxyz")!=0)
   {
    fprintf(stderr,"not a binpos amber coordinate file\n");
    return NULL;
   }
  fprintf(stderr,"Proceeding to open amber7 binpos coordinate file\n");
  fread(&igarb,sizeof(int),1,fd);
  point=ftell(fd);

/* Check for endianism here*/
  if(igarb>1000000000)
   {
    fprintf(stderr, "File '%s' appears to be other-endian.\n", path);
    binpos->wrongendian = 1;
    memcpy(lenbuf, (const char *)&igarb, 4);
    tmpc = lenbuf[0]; lenbuf[0] = lenbuf[3]; lenbuf[3] = tmpc;
    tmpc = lenbuf[1]; lenbuf[1] = lenbuf[2]; lenbuf[2] = tmpc;
    memcpy((char *)&igarb, lenbuf, 4);

	if((fseek(fd, point, SEEK_SET))!=0)
      {
	  fprintf(stderr,"Endian correction failed. er=%d\n",er);
      return NULL;	
     }
	fseek(fd, point, SEEK_SET);
   }
  binpos->fd = fd;
  binpos->numatoms = igarb;
  binpos->xyz = (float *)malloc(3 * binpos->numatoms * sizeof(float));

  if (!binpos->xyz) {
    fprintf(stderr, "Unable to allocate space for %d atoms.\n", binpos->numatoms);
    fclose(fd);
    free(binpos);
    return NULL;
  }
  *natoms = binpos->numatoms;
  return binpos;
}

static int read_next_timestep(void *v, int natoms, molfile_timestep_t *ts) 
 {
  binposhandle *binpos;
  int i, numatoms,igarb;
  char *cdata;
  char tmpc;

  binpos = (binposhandle *)v;
  if (!binpos->fd) 
    return MOLFILE_ERROR;  /* Done reading frames */

  numatoms = binpos->numatoms;

  if (fread(binpos->xyz, sizeof(float), 3 * numatoms, binpos->fd)
	                         != (size_t)(3 * numatoms)) {
    fprintf(stderr, "Failure reading data from amber7 binary file.\n");
    return MOLFILE_ERROR;
  }

  if (binpos->wrongendian) {

/*For float or single precision endian conversion*/
/*amber7 binpos files are always float not doubles*/
    cdata = (char *) binpos->xyz;
    for ( i=0; i<3*numatoms; ++i, cdata+=4 ) {
    tmpc = cdata[0]; cdata[0] = cdata[3]; cdata[3] = tmpc;
    tmpc = cdata[1]; cdata[1] = cdata[2]; cdata[2] = tmpc;
    }
  }

  if (ts) {
    for ( i=0; i<numatoms; ++i) {
      ts->coords[3*i] = binpos->xyz[3*i];
      ts->coords[3*i+1] = binpos->xyz[3*i+1];
      ts->coords[3*i+2] = binpos->xyz[3*i+2];
    }
  }
  /*
   * Close the file handle and set to NULL so we know we're done reading 
   * 
   */

  if((fread(&igarb,sizeof(int),1,binpos->fd))!=1)
   {
    fclose(binpos->fd);
    binpos->fd = NULL;
   }
  return MOLFILE_SUCCESS;
}
 
static void close_file_read(void *v) {
  binposhandle *binpos = (binposhandle *)v;
  if (binpos->fd)
    fclose(binpos->fd);
  free(binpos->xyz);
  free(binpos);
}

static void *open_binpos_write(const char *path, const char *filetype, 
    int natoms) {
  binposhandle *binpos;
  FILE *fd;

  fd = fopen(path, "wb");
  if (!fd) {
    fprintf(stderr, "Could not open file %s for writing\n", path);
    return NULL;
  }
  fprintf(stderr,"Writing file in current machine endian-ism\n");
  binpos = (binposhandle *)malloc(sizeof(binposhandle));
  binpos->fd = fd;
  binpos->numatoms = natoms;
  fwrite( "fxyz", 4, 1, binpos->fd);
  return binpos;
}

static int write_timestep(void *v, const molfile_timestep_t *ts) {
  
  int i,numatoms;

  binposhandle *binpos = (binposhandle *)v;
  
  if (!binpos->fd)
    return MOLFILE_ERROR;

/*add the number of atoms in between frames*/
  /*myint = (binpos_int32)binpos->numatoms;*/
  numatoms = binpos->numatoms;

  fwrite(&numatoms, 4, 1, binpos->fd);

  for (i=0; i<3*numatoms; i++) 
   {
    float tmp = ts->coords[i];
    if (fwrite(&tmp, sizeof(float), 1, binpos->fd) != 1) {
      fprintf(stderr, "Error writing amber7 binary file\n");
      return MOLFILE_ERROR;
     }
   }

  /*
   * Close and NULLify the file handle so we don't write any more frames.
   */
/*  fclose(binpos->fd);
  binpos->fd = NULL;*/

  return MOLFILE_SUCCESS;
}
       
static void close_file_write(void *v) {
  binposhandle *binpos = (binposhandle *)v;
  if (binpos->fd)
    fclose(binpos->fd);
  free(binpos);
}

/*
 * Initialization stuff here
 */

static molfile_plugin_t plugin;

VMDPLUGIN_API int VMDPLUGIN_init() {
  memset(&plugin, 0, sizeof(molfile_plugin_t));
  plugin.abiversion = vmdplugin_ABIVERSION;
  plugin.type = MOLFILE_PLUGIN_TYPE;
  plugin.name = "binpos";
  plugin.prettyname = "Scripps Binpos";
  plugin.author = "Brian Bennion";
  plugin.majorv = 0;
  plugin.minorv = 4;
  plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
  plugin.filename_extension = "binpos";
  plugin.open_file_read = open_binpos_read;
  plugin.read_next_timestep = read_next_timestep;
  plugin.close_file_read = close_file_read;
  plugin.open_file_write = open_binpos_write;
  plugin.write_timestep = write_timestep;
  plugin.close_file_write = close_file_write;
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)(void *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_fini() {
  return VMDPLUGIN_SUCCESS;
}

