/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_dlpolyplugin
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
 *      $RCSfile: dlpolyplugin.c,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.21 $       $Date: 2009/04/29 15:45:29 $
 *
 ***************************************************************************/

/*
 * DLPOLY formatted history file format:
 *   http:

 *   http:

 *
 */

#include "largefiles.h"   /* platform dependent 64-bit file I/O defines */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "molfile_plugin.h"

#ifndef M_PI_2
#define M_PI_2 1.57079632679489661922
#endif

typedef struct {
  FILE *file;
  int dlpolyversion;
  int numatoms;
  int cellwarnflag;
} dlpolydata;
 

static void *open_dlpoly_read(const char *filename, const char *filetype, 
                           int *natoms) {
  FILE *fd;
  dlpolydata *data;
  char fbuffer[4096], buf[4096];
  int dlpolyversion, scancount, nstep, keytrj, atomcount;

  fd = fopen(filename, "rb");
  if (!fd) return NULL;

  if (NULL == fgets(fbuffer, 1024, fd))  
    return NULL;

  if (!strcmp(filetype, "dlpolyhist")) {
    dlpolyversion=2;
  } else if (!strcmp(filetype, "dlpoly3hist")) {
    dlpolyversion=3;
  } else {
    dlpolyversion=2;
  }
 
  /* check to see if the first line is a "timestep" record */ 
  scancount = sscanf(fbuffer, "%s %d %d", buf, &nstep, natoms);
  if (scancount != 3 || strcmp(buf, "timestep") != 0) {
    /* if not a timestep, it might have the normal header on it      */
    /* in which case we'll skip the first line, and parse the second */
    if (NULL == fgets(fbuffer, 1024, fd))  
      return NULL;
    scancount = sscanf(fbuffer, "%d %d %d", &keytrj, &nstep, natoms);
    if (scancount != 3) {
      printf("open_dlpoly_read) unrecognized header record\n");
      return NULL;
    } 

    /* now check the first timestep record for safety */
    if (NULL == fgets(fbuffer, 1024, fd))  
      return NULL;
    scancount = sscanf(fbuffer, "%s %d %d", buf, &nstep, &atomcount);
    if (scancount != 3 || strcmp(buf, "timestep") != 0) {
      printf("open_dlpoly_read) unrecognized timestep record\n");
      return NULL;
    }

    if (atomcount != *natoms) {
      printf("open_dlpoly_read) mismatched atom count\n");
      return NULL;
    }
  }
 
  data = (dlpolydata *) malloc(sizeof(dlpolydata));
  memset(data, 0, sizeof(dlpolydata));
  data->file = fd;
  data->numatoms= *natoms;
  data->cellwarnflag = 0;
  data->dlpolyversion = dlpolyversion;

  rewind(data->file); /* prepare for first read_timestep call */
  return data;
}


static int read_dlpoly_structure(void *mydata, int *optflags,
                              molfile_atom_t *atoms) {
  char fbuffer[4096], buf[4096];
  float rsd, x, y, z;
  int i, nstep, atomcount, keytrj, imcon, scancount, atomid, atomcount2;
  float tstep, mass, charge;
  dlpolydata *data = (dlpolydata *)mydata;
 
  /* we don't read any optional data */
  *optflags = MOLFILE_NOOPTIONS;


  /* if we get nothing, assume we hit end of file */
  if (NULL == fgets(fbuffer, 1024, data->file))  
    return MOLFILE_EOF;

  /* check to see if the first line is a "timestep" record */
  scancount = sscanf(fbuffer, "%s %d %d %d %d %f", buf, 
                     &nstep, &atomcount, &keytrj, &imcon, &tstep);
  if (scancount != 6 || strcmp(buf, "timestep") != 0) {
    /* if not a timestep, it might have the normal header on it      */
    /* in which case we'll skip the first line, and parse the second */
    if (NULL == fgets(fbuffer, 1024, data->file))
      return MOLFILE_EOF;
    scancount = sscanf(fbuffer, "%d %d %d", &keytrj, &nstep, &atomcount);
    if (scancount != 3) {
      printf("dlpoly structure) unrecognized header record\n");
      return MOLFILE_ERROR;
    } 

    /* now check the first timestep record for safety */
    if (NULL == fgets(fbuffer, 1024, data->file))
      return MOLFILE_EOF;
    scancount = sscanf(fbuffer, "%s %d %d %d %d %f", buf, 
                       &nstep, &atomcount2, &keytrj, &imcon, &tstep);
    if (scancount != 6 || strcmp(buf, "timestep") != 0) {
      printf("dlpoly structure) unrecognized timestep record\n");
      return MOLFILE_ERROR;
    }

    /* check atom count */
    if (atomcount != atomcount2) {
      printf("dlpoly structure) header/timestep mismatched atom count\n");
      return MOLFILE_ERROR;
    }
  }

  /* check atom count */
  if (atomcount != data->numatoms) {
    printf("dlpoly structure) mismatched atom count\n");
    return MOLFILE_ERROR;
  }

  /* read periodic cell vectors */
  if (imcon > 0) {
    float xaxis[3];
    float yaxis[3];
    float zaxis[3];

    /* eat the data but don't use it for anything */
    if (3 != fscanf(data->file, "%f %f %f\n", &xaxis[0], &xaxis[1], &xaxis[2])) {
      printf("dlpoly structure) failed reading unit cell basis vectors\n");
      return MOLFILE_ERROR;
    }
    if (3 != fscanf(data->file, "%f %f %f\n", &yaxis[0], &yaxis[1], &yaxis[2])) {
      printf("dlpoly structure) failed reading unit cell basis vectors\n");
      return MOLFILE_ERROR;
    }
    if (3 != fscanf(data->file, "%f %f %f\n", &zaxis[0], &zaxis[1], &zaxis[2])) {
      printf("dlpoly structure) failed reading unit cell basis vectors\n");
      return MOLFILE_ERROR;
    }
  }

  for (i=0; i<data->numatoms; i++) {
    molfile_atom_t *atom = atoms + i;

    /* read the coordinates */
    switch (data->dlpolyversion) {
      case 3:
        if (8 != fscanf(data->file, "%s %d %f %f %f %f %f %f",
               buf, &atomid, &mass, &charge, &rsd, &x, &y, &z)) {
          printf("dlpoly structure v3) failed parsing atom coordinates\n");
          return MOLFILE_ERROR;
        }
        break;

      case 2:
      default: 
        if (7 != fscanf(data->file, "%s %d %f %f %f %f %f",
               buf, &atomid, &mass, &charge, &x, &y, &z)) {
          printf("dlpoly structure v2) failed parsing atom coordinates\n");
          return MOLFILE_ERROR;
        }
        break;
    }

    /* read the velocities */
    if (keytrj > 0) {
      float xv, yv, zv;
      if (3 != fscanf(data->file, "%f %f %f", &xv, &yv, &zv)) {
        printf("dlpoly structure) failed reading atom velocities\n");
        return MOLFILE_ERROR;
      }
    }

    /* read the forces */
    if (keytrj > 1) {
      float xf, yf, zf;
      if (3 != fscanf(data->file, "%f %f %f", &xf, &yf, &zf)) {
        printf("dlpoly structure) failed reading atom forces\n");
        return MOLFILE_ERROR;
      }
    }

    strncpy(atom->name, buf, sizeof(atom->name));
    strncpy(atom->type, atom->name, sizeof(atom->type));
    atom->resname[0] = '\0';
    atom->resid = 1;
    atom->chain[0] = '\0';
    atom->segid[0] = '\0';
  }

  rewind(data->file);
  return MOLFILE_SUCCESS;
}


static int read_dlpoly_timestep(void *mydata, int natoms, molfile_timestep_t *ts) {
  char fbuffer[4096], buf[4096];
  float rsd, x, y, z;
  int i, nstep, atomcount, keytrj, imcon, scancount, atomid, atomcount2;
  float tstep, mass, charge, elapsedps;
  dlpolydata *data = (dlpolydata *)mydata;

  /* if we get nothing, assume we hit end of file */
  if (NULL == fgets(fbuffer, 1024, data->file))  
    return MOLFILE_EOF;

  switch (data->dlpolyversion) {
    case 3:
      scancount = sscanf(fbuffer, "%s %d %d %d %d %f %f", buf, 
                    &nstep, &atomcount, &keytrj, &imcon, &tstep, &elapsedps);
      if (scancount != 7 || strcmp(buf, "timestep") != 0) {
        /* if not a timestep, it might have the normal header on it      */
        /* in which case we'll skip the first line, and parse the second */
        if (NULL == fgets(fbuffer, 1024, data->file))
          return MOLFILE_EOF;
        scancount = sscanf(fbuffer, "%d %d %d", &keytrj, &nstep, &atomcount);
        if (scancount != 3) {
          printf("dlpoly timestep v3) unrecognized header record\n");
          return MOLFILE_ERROR;
        } 

        /* now check the first timestep record for safety */
        if (NULL == fgets(fbuffer, 1024, data->file))
          return MOLFILE_EOF;
        scancount = sscanf(fbuffer, "%s %d %d %d %d %f %f", buf, 
                    &nstep, &atomcount2, &keytrj, &imcon, &tstep, &elapsedps);
        if (scancount != 7 || strcmp(buf, "timestep") != 0) {
          printf("dlpoly timestep v3) unrecognized timestep record\n");
          return MOLFILE_ERROR;
        }

        /* check atom count */
        if (atomcount != atomcount2) {
          printf("dlpoly timestep v3) header/timestep mismatched atom count\n");
          return MOLFILE_ERROR;
        }
      }
      break;

    case 2:
    default: 
      scancount = sscanf(fbuffer, "%s %d %d %d %d %f", buf, 
                         &nstep, &atomcount, &keytrj, &imcon, &tstep);
      if (scancount != 6 || strcmp(buf, "timestep") != 0) {
        /* if not a timestep, it might have the normal header on it      */
        /* in which case we'll skip the first line, and parse the second */
        if (NULL == fgets(fbuffer, 1024, data->file))
          return MOLFILE_EOF;
        scancount = sscanf(fbuffer, "%d %d %d", &keytrj, &nstep, &atomcount);
        if (scancount != 3) {
          printf("dlpoly timestep v2) unrecognized header record\n");
          return MOLFILE_ERROR;
        } 

        /* now check the first timestep record for safety */
        if (NULL == fgets(fbuffer, 1024, data->file))
          return MOLFILE_EOF;
        scancount = sscanf(fbuffer, "%s %d %d %d %d %f", buf, 
                           &nstep, &atomcount2, &keytrj, &imcon, &tstep);
        if (scancount != 6 || strcmp(buf, "timestep") != 0) {
          printf("dlpoly timestep v2) unrecognized timestep record\n");
          return MOLFILE_ERROR;
        }

        /* check atom count */
        if (atomcount != atomcount2) {
          printf("dlpoly timestep v2) header/timestep mismatched atom count\n");
          return MOLFILE_ERROR;
        }
      }
      break;
  }

  /* check atom count */
  if (atomcount != natoms) {
    printf("dlpoly timestep) mismatched atom count\n");
    return MOLFILE_ERROR;
  }

  /* read periodic cell vectors */
  if (imcon > 0) {
    float xaxis[3];
    float yaxis[3];
    float zaxis[3];

    if (3 != fscanf(data->file, "%f %f %f\n", &xaxis[0], &xaxis[1], &xaxis[2])) {
      printf("dlpoly timestep) failed reading unit cell basis vectors\n");
      return MOLFILE_ERROR;
    }

    if (3 != fscanf(data->file, "%f %f %f\n", &yaxis[0], &yaxis[1], &yaxis[2])) {
      printf("dlpoly timestep) failed reading unit cell basis vectors\n");
      return MOLFILE_ERROR;
    }

    if (3 != fscanf(data->file, "%f %f %f\n", &zaxis[0], &zaxis[1], &zaxis[2])) {
      printf("dlpoly timestep) failed reading unit cell basis vectors\n");
      return MOLFILE_ERROR;
    }

    /* check and copy in periodic cell info */
    if (ts != NULL) {
      if (xaxis[1] != 0.0 || xaxis[2] != 0.0 || 
          yaxis[0] != 0.0 || yaxis[2] != 0.0 || 
          zaxis[0] != 0.0 || xaxis[1] != 0.0) {
        if (data->cellwarnflag != 1)
          printf("dlpoly timestep) non-orthogonal DLPOLY periodic cell data unsupported\n");
        data->cellwarnflag = 1;
      } else {
        ts->A = xaxis[0];
        ts->B = yaxis[1];
        ts->C = zaxis[2];
        if (data->cellwarnflag != 2)
          printf("dlpoly timestep) converting DLPOLY periodic cell data\n");
        data->cellwarnflag = 2;
      }
    }
  }

  /* read all per-atom data */
  for (i=0; i<natoms; i++) {
    /* read the coordinates */
    switch (data->dlpolyversion) {
      case 3:
        if (8 != fscanf(data->file, "%s %d %f %f %f %f %f %f",
               buf, &atomid, &mass, &charge, &rsd, &x, &y, &z)) {
          printf("dlpoly timestep v3) failed parsing atom coordinates\n");
          return MOLFILE_ERROR;
        }
        break;

      case 2:
      default: 
        if (7 != fscanf(data->file, "%s %d %f %f %f %f %f",
               buf, &atomid, &mass, &charge, &x, &y, &z)) {
          printf("dlpoly timestep v2) failed parsing atom coordinates\n");
          return MOLFILE_ERROR;
        }
        break;
    }

    /* read the velocities */
    if (keytrj > 0) {
      float xv, yv, zv;
      if (3 != fscanf(data->file, "%f %f %f", &xv, &yv, &zv)) {
        printf("dlpoly timestep) failed parsing atom velocities\n");
        return MOLFILE_ERROR;
      }
    }

    /* read the forces */
    if (keytrj > 1) {
      float xf, yf, zf;
      if (3 != fscanf(data->file, "%f %f %f", &xf, &yf, &zf)) {
        printf("dlpoly timestep) failed parsing atom forces\n");
        return MOLFILE_ERROR;
      }
    }

    /* only save coords if we're given a timestep pointer, */
    /* otherwise assume that VMD wants us to skip past it. */
    if (ts != NULL) { 
      if (atomid > 0 && atomid <= natoms) {
        int addr = 3 * (atomid - 1);
        ts->coords[addr    ] = x;
        ts->coords[addr + 1] = y;
        ts->coords[addr + 2] = z;
      } else {
        fprintf(stderr, "dlpoly timestep) ignoring illegal atom index %d\n", atomid);
      }
    } 
  }

  /* eat remaining spaces and LF character */
  fgets(fbuffer, 1024, data->file);

  return MOLFILE_SUCCESS;
}
    
static void close_dlpoly_read(void *mydata) {
  dlpolydata *data = (dlpolydata *)mydata;
  fclose(data->file);
  free(data);
}


/* registration stuff */
static molfile_plugin_t dlpoly2plugin;
static molfile_plugin_t dlpoly3plugin;

VMDPLUGIN_API int VMDPLUGIN_init() {
  memset(&dlpoly2plugin, 0, sizeof(molfile_plugin_t));
  dlpoly2plugin.abiversion = vmdplugin_ABIVERSION;
  dlpoly2plugin.type = MOLFILE_PLUGIN_TYPE;
  dlpoly2plugin.name = "dlpolyhist";
  dlpoly2plugin.prettyname = "DLPOLY V2 History";
  dlpoly2plugin.author = "John Stone";
  dlpoly2plugin.majorv = 0;
  dlpoly2plugin.minorv = 8;
  dlpoly2plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
  dlpoly2plugin.filename_extension = "dlpolyhist";
  dlpoly2plugin.open_file_read = open_dlpoly_read;
  dlpoly2plugin.read_structure = read_dlpoly_structure;
  dlpoly2plugin.read_next_timestep = read_dlpoly_timestep;
  dlpoly2plugin.close_file_read = close_dlpoly_read;

  memset(&dlpoly3plugin, 0, sizeof(molfile_plugin_t));
  dlpoly3plugin.abiversion = vmdplugin_ABIVERSION;
  dlpoly3plugin.type = MOLFILE_PLUGIN_TYPE;
  dlpoly3plugin.name = "dlpoly3hist";
  dlpoly3plugin.prettyname = "DLPOLY V3 History";
  dlpoly3plugin.author = "John Stone";
  dlpoly3plugin.majorv = 0;
  dlpoly3plugin.minorv = 8;
  dlpoly3plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
  dlpoly3plugin.filename_extension = "dlpolyhist";
  dlpoly3plugin.open_file_read = open_dlpoly_read;
  dlpoly3plugin.read_structure = read_dlpoly_structure;
  dlpoly3plugin.read_next_timestep = read_dlpoly_timestep;
  dlpoly3plugin.close_file_read = close_dlpoly_read;

  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)&dlpoly2plugin);
  (*cb)(v, (vmdplugin_t *)&dlpoly3plugin);
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
    v = open_dlpoly_read(*argv, "dlpoly", &natoms);
    if (!v) {
      fprintf(stderr, "open_dlpoly_read failed for file %s\n", *argv);
      return 1;
    }
    fprintf(stderr, "open_dlpoly_read succeeded for file %s\n", *argv);
    fprintf(stderr, "number of atoms: %d\n", natoms);

    i = 0;
    timestep.coords = (float *)malloc(3*sizeof(float)*natoms);
    while (!read_dlpoly_timestep(v, natoms, &timestep)) {
      i++;
    }
    fprintf(stderr, "ended read_next_timestep on frame %d\n", i);

    close_dlpoly_read(v);
  }
  return 0;
}

#endif

