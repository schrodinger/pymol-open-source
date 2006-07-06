/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_rst7plugin
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
 *      $RCSfile: rst7plugin.c,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.16 $       $Date: 2006/02/23 19:36:45 $
 *
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "molfile_plugin.h"

typedef struct {
  FILE *file;
  int has_box;
  int numatoms;
  int count;
  int rstfile;
} rstdata;

static void *open_rst_read(const char *filename, const char *filetype,int *natoms) {
  FILE *fd;
  rstdata *data;
  int numats=0,i,j,point2,kkk=1; 
  char title[82], *field;
  char line[82];
  float x, y, z,a=0.0,b=0.0,c=0.0;
  double  timesteprst;

  /* Amber 7'coord' restart files have a second introduction line with 
   * possibly 2 entries only check for one now...
   * they include three 90.00 ter cards at the end
   * need to fix this, real crd files have atom record but no timestep and no
   * velocity info...arggggg
   */
  fd = fopen(filename, "rb");
  if (!fd) 
    return NULL; /* failure */

  data = (rstdata *)malloc(sizeof(rstdata));
  memset(data, 0, sizeof(rstdata));
  fgets(title, 82, fd);
  printf("Title: %s\n",title);

  fgets(line, 82, fd);
  while (kkk==1) {
    /* try to read first field */
    field = strtok(line, "  ");
    if (field==NULL) {                
      continue; /* no fields at all on this line */
    }
    numats = atoi(field);

    /* try to read second field will be null if not there */
    field = strtok(line, "  ");
    if (field==NULL) {
      kkk=0;
      printf("This file has no velocity info.\n");
    } else {
      timesteprst = strtod(field, NULL);
      printf("This file contains velocity info.\n");
      kkk=0;
    }
  }

  point2=ftell(fd);
  data->file = fd;
  printf("The Restartcrd has %d atoms.\n",numats);
  for (i=0; i<numats; i++) {
    j = fscanf(fd, "%f %f %f", &x, &y, &z);
  }

  j = fscanf(fd, "%f %f %f %f %f %f", &x, &y, &z,&a,&b,&c);
  if (j != EOF) {
    printf("This restartcrd file has box info.\n");
    data->has_box=1;
    if((int)a==90) {
      printf("Box Dimensions are %f  %f  %f  %f  %f  %f\n",x,y,z,a,b,c);
    } else {
      for (i=0; i<numats-2; i++) {
        j = fscanf(fd, "%f %f %f", &x, &y, &z);
      }
      j = fscanf(fd, "%f %f %f %f %f %f", &x, &y, &z,&a,&b,&c);
      if (j != EOF) {
        if((int)a==90) {
          printf("Box Dimensions are %f  %f  %f  %f  %f  %f\n",x,y,z,a,b,c);
        }
      }
    }
  } 

  *natoms=numats;
  data->numatoms=numats;
  data->rstfile=1;
  fseek(fd,point2,SEEK_SET);

  return data;
}

static int read_rst_timestep(void *mydata, int natoms, molfile_timestep_t *ts) {
  rstdata *rst= (rstdata *)mydata;
  int i, j;
  float x, y, z;

  /* check for rst and first read through already taken place */
  if(rst->count==1 && rst->rstfile==1) 
    return MOLFILE_ERROR; 

  for (i=0; i<rst->numatoms; i++)  {
    /* changed to i=1 BB */
    j = fscanf(rst->file, "%f %f %f", &x, &y, &z);
    if (j == EOF) {
      return MOLFILE_ERROR;
    } else if (j <= 0) {
      fprintf(stderr, "Problem reading CRD file\n");
      return MOLFILE_ERROR;
    }
    ts->coords[3*i] = x;
    ts->coords[3*i+1] = y;
    ts->coords[3*i+2] = z;
  }

  /* Read in optional velocity data.  Units are Angstroms per 1/20.455ps. */
  /* XXX Not currently implemented.  This should be added sometime soon.  */
 
  /* Don't Read box info.  No use for it yet.. */
  rst->count++;
  /* printf("rst->count: %d\n",rst->count); */

  return MOLFILE_SUCCESS;
}
    
static void close_rst_read(void *mydata) {
  rstdata *rst= (rstdata *)mydata;
  fclose(rst->file);
  free(rst);
}

static void *open_rst_write(const char *path, const char *filetype, int natoms) {
  /* Not fixed for writing proper rsts yet....BB */
  rstdata *rst;
  FILE *fd;

  fd = fopen(path, "wb");
  if (!fd) {
    fprintf(stderr, "Could not open file %s for writing\n", path);
    return NULL;
  }
  fprintf(fd, "TITLE : Created by VMD with %d atoms\n",natoms);
  
  rst = (rstdata *)malloc(sizeof(rstdata));
  rst->file = fd;
  rst->numatoms = natoms;
  rst->has_box = strcmp(filetype, "rst"); 
  return rst;
}    
  
static int write_rst_timestep(void *v, const molfile_timestep_t *ts) {
  rstdata *rst = (rstdata *)v;
  int i;
  const int ndata = rst->numatoms * 3;
  for (i=0; i<ndata; i++) {
    fprintf(rst->file, "%8.3f", ts->coords[i]);
    if (i % 10 == 0) fprintf(rst->file, "\n"); 
  }
  if (rst->has_box) {
    fprintf (rst->file, "\n0.000 0.000 0.000\n");
  }

  return MOLFILE_SUCCESS;
}

static void close_rst_write(void *v) {
  rstdata *rst = (rstdata *)v;
  fclose(rst->file);
  free(rst);
}

/* registration stuff */
static molfile_plugin_t rst7plugin = {
    vmdplugin_ABIVERSION,   /* ABI version */
    MOLFILE_PLUGIN_TYPE,    /* type of plugin */
    "rst7",                 /* short name of plugin */
    "AMBER7 Restart",       /* pretty name of plugin */
    "Brian Bennion",        /* authors */
    0,                      /* major version */
    2,                      /* minor version */
    VMDPLUGIN_THREADUNSAFE, /* is not reentrant */
    "rst7",                 /* filename extension */
    open_rst_read,
    0,
    0,
    read_rst_timestep,
    close_rst_read,
    open_rst_write,
    0,
    write_rst_timestep,
    close_rst_write
};

VMDPLUGIN_API int VMDPLUGIN_init(){
  return VMDPLUGIN_SUCCESS;
}
VMDPLUGIN_API int VMDPLUGIN_fini(){
  return VMDPLUGIN_SUCCESS;
}
VMDPLUGIN_API int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)&rst7plugin);
  return VMDPLUGIN_SUCCESS;
}
