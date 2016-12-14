/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_rst7plugin
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
 *      $RCSfile: rst7plugin.c,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.20 $       $Date: 2016/11/28 05:01:54 $
 *
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "molfile_plugin.h"
#include "vmdconio.h"

typedef struct {
  FILE *file;
  int has_box;
  int has_vels;
  int numatoms;
  int count;
  int rstfile;
#if vmdplugin_ABIVERSION > 10
  molfile_timestep_metadata_t ts_meta;
#endif
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
#if vmdplugin_ABIVERSION > 10
  data->ts_meta.count = -1;
  data->ts_meta.has_velocities = 0;
#endif

  fgets(title, 82, fd);
  vmdcon_printf(VMDCON_INFO, "rst7plugin) Title: %s\n",title);

  fgets(line, 82, fd);
  while (kkk==1) {
    /* try to read first field */
    field = strtok(line, " \t");
    if (field==NULL) {                
      continue; /* no fields at all on this line */
    }
    numats = atoi(field);

    /* try to read second field will be null if not there */
    field = strtok(NULL, " \t");
    if (field==NULL) {
      kkk=0;
      vmdcon_printf(VMDCON_INFO, "rst7plugin) This file has no velocity info.\n");
      data->has_vels=0;
    } else {
      timesteprst = strtod(field, NULL);
      vmdcon_printf(VMDCON_INFO, "rst7plugin) This file contains velocity info.\n");
      data->has_vels=1;
#if vmdplugin_ABIVERSION > 10
      data->ts_meta.has_velocities = 1;
#endif
      kkk=0;
    }
  }

  point2=ftell(fd);
  data->file = fd;
  vmdcon_printf(VMDCON_INFO, "rst7plugin) The Restartcrd has %d atoms.\n",numats);

  /* skip over coordinate data */
  for (i=0; i<numats; i++) {
    j = fscanf(fd, "%f%f%f", &x, &y, &z);
  }

  /* skip over velocity data, if present */
  if (data->has_vels) {
    for (i=0; i<numats; i++) {
      j = fscanf(fd, "%f%f%f", &x, &y, &z);
    }
  }
  
  j = fscanf(fd, "%f%f%f%f%f%f", &x, &y, &z,&a,&b,&c);
  if (j != EOF) {
    vmdcon_printf(VMDCON_INFO, "rst7plugin) This restartcrd file has box info.\n");
    data->has_box=1;
    vmdcon_printf(VMDCON_INFO, "rst7plugin) Box Dimensions are %f %f %f %f %f %f\n",x,y,z,a,b,c);
  }

  *natoms=numats;
  data->numatoms=numats;
  data->rstfile=1;
  fseek(fd,point2,SEEK_SET);

  return data;
}

#if vmdplugin_ABIVERSION > 10
static int read_timestep_metadata(void *mydata,
                                  molfile_timestep_metadata_t *meta) {
  rstdata *data = (rstdata *)mydata;
  
  meta->count = -1;
  meta->has_velocities = data->ts_meta.has_velocities;
  if (meta->has_velocities) {
    vmdcon_printf(VMDCON_INFO,
                  "rst7plugin) Importing velocities from restart file.\n");
  }
  return MOLFILE_SUCCESS;
}
#endif

static int read_rst_timestep(void *mydata, int natoms, molfile_timestep_t *ts) {
  rstdata *rst= (rstdata *)mydata;
  int i, j;
  float x, y, z, a, b, c;
  
  /* check for rst and first read through already taken place */
  if(rst->count==1 && rst->rstfile==1) 
    return MOLFILE_ERROR; 

  ts->A = ts->B = ts->C = 0.0f;
  ts->alpha = ts->beta = ts->gamma = 90.0f;
  
  for (i=0; i<rst->numatoms; i++)  {
    /* changed to i=1 BB */
    j = fscanf(rst->file, "%f%f%f", &x, &y, &z);
    if (j == EOF) {
      return MOLFILE_ERROR;
    } else if (j <= 0) {
      vmdcon_printf(VMDCON_ERROR, "rst7plugin) Problem reading CRD file\n");
      return MOLFILE_ERROR;
    }
    ts->coords[3*i] = x;
    ts->coords[3*i+1] = y;
    ts->coords[3*i+2] = z;
  }

  if (rst->has_vels) {
    /* Read in optional velocity data.  Units are Angstroms per 1/20.455ps. */
    for (i=0; i<rst->numatoms; i++)  {
      j = fscanf(rst->file, "%f%f%f", &x, &y, &z);
      if (j == EOF) {
        return MOLFILE_ERROR;
      } else if (j <= 0) {
        vmdcon_printf(VMDCON_ERROR, "rst7plugin) Problem reading velocities\n");
        return MOLFILE_ERROR;
      }
#if vmdplugin_ABIVERSION > 10
      if (ts->velocities != NULL) {
        ts->velocities[3*i] = x;
        ts->velocities[3*i+1] = y;
        ts->velocities[3*i+2] = z;
      }
#endif
    }
  }
 
  if (rst->has_box) {
    j = fscanf(rst->file, "%f%f%f%f%f%f", &x, &y, &z, &a, &b, &c);
    if (j == EOF) {
      vmdcon_printf(VMDCON_ERROR, "rst7plugin) Problem reading box data\n");
      return MOLFILE_ERROR;
    }
    ts->A = x;
    ts->B = y;
    ts->C = z;
    ts->alpha = a;
    ts->beta = b;
    ts->gamma = c;
  }
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
  char title[82];
  rstdata *rst;
  FILE *fd;
  int len;

  fd = fopen(path, "wb");
  if (!fd) {
    vmdcon_printf(VMDCON_ERROR, "rst7plugin) Could not open file %s for writing\n", path);
    return NULL;
  }
  /* write out fixed length fortran style string */
  sprintf(title, "TITLE : Created by VMD with %d atoms",natoms);
  len = strlen(title);
  memset(title+len,(int)' ',82-len);
  title[80] = '\n';
  title[81] = '\0';
  fputs(title,fd);

  rst = (rstdata *)malloc(sizeof(rstdata));
  rst->file = fd;
  rst->numatoms = natoms;
  rst->has_box = 1;
  return rst;
}
  
static int write_rst_timestep(void *v, const molfile_timestep_t *ts) {
  rstdata *rst = (rstdata *)v;
  int i;
  const int ndata = rst->numatoms * 3;

#if vmdplugin_ABIVERSION > 10
  if (ts->velocities != NULL) {
    fprintf(rst->file, "%10d %13.7g\n", rst->numatoms, ts->physical_time);
  } else
#endif
    fprintf(rst->file, "%10d\n", rst->numatoms);

  for (i=0; i<ndata; i++) {
    fprintf(rst->file, "%12.7f", ts->coords[i]);
    if ((i+1) % 6 == 0) fprintf(rst->file, "\n"); 
  }
  if (ndata % 6 != 0) fprintf(rst->file,"\n");

#if vmdplugin_ABIVERSION > 10
  if (ts->velocities != NULL) {
    for (i=0; i<ndata; i++) {
      fprintf(rst->file, "%12.7f", ts->velocities[i]);
      if ((i+1) % 6 == 0) fprintf(rst->file, "\n"); 
    }
    if (ndata % 6 != 0) fprintf(rst->file,"\n");
  }
#endif

  fprintf (rst->file, "%12.7f%12.7f%12.7f%12.7f%12.7f%12.7f\n",
           ts->A, ts->B, ts->C, ts->alpha, ts->beta, ts->gamma);

  return MOLFILE_SUCCESS;
}

static void close_rst_write(void *v) {
  rstdata *rst = (rstdata *)v;
  fclose(rst->file);
  free(rst);
}

/* registration stuff */
static molfile_plugin_t plugin;

VMDPLUGIN_API int VMDPLUGIN_init(){
  memset(&plugin, 0, sizeof(molfile_plugin_t));
  plugin.abiversion = vmdplugin_ABIVERSION;
  plugin.type = MOLFILE_PLUGIN_TYPE;
  plugin.name = "rst7";
  plugin.prettyname = "AMBER7 Restart";
  plugin.author = "Brian Bennion, Axel Kohlmeyer";
  plugin.majorv = 0;
  plugin.minorv = 4;
  plugin.is_reentrant = VMDPLUGIN_THREADUNSAFE;
  plugin.filename_extension = "rst7";
  plugin.open_file_read = open_rst_read;
  plugin.read_next_timestep = read_rst_timestep;
#if vmdplugin_ABIVERSION > 10
  plugin.read_timestep_metadata = read_timestep_metadata;
#endif
  plugin.close_file_read = close_rst_read;
  plugin.open_file_write = open_rst_write;
  plugin.write_timestep = write_rst_timestep;
  plugin.close_file_write = close_rst_write;
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_fini(){ return VMDPLUGIN_SUCCESS; }

