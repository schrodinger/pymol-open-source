/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_msmsplugin
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
 *      $RCSfile: msmsplugin.C,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.7 $       $Date: 2006/02/23 19:36:45 $
 *
 ***************************************************************************/

/* 
 * Reader for MSMS surface files
 */
#include <stdio.h>
#include <stdlib.h>
//#include <math.h>
#include <string.h>
#include "molfile_plugin.h"

typedef struct {
  FILE *ffd;
  FILE *vfd;
  molfile_graphics_t *graphics;
} msms_t;

static void *open_file_read(const char *filepath, const char *filetype,
    int *natoms) {
  FILE *ffd; // face file
  FILE *vfd; // vertex file
  msms_t *msms;
  char * facefilepath;
  char * vertfilepath;
  char * cp;

  int filenamelen = strlen(filepath);
  facefilepath = (char *) malloc(filenamelen + 10);
  vertfilepath = (char *) malloc(filenamelen + 10);
  strcpy(facefilepath, filepath);
  strcpy(vertfilepath, filepath);

  // require the MSMS output filenames to match what MSMS does
  // If the user selected the .face file or the .vert file, we should
  // be able to cope either way by assigning the right filenames to the
  // right strings and getting the right files opened accordingly.
  cp = strstr(facefilepath, ".face");
  if (cp == NULL) {
    cp = strstr(facefilepath, ".vert");
    if (cp != NULL) {
       strcpy(cp, ".face");
    } else {
      printf("msmsplugin) file names don't match expected MSMS output\n");
      free(facefilepath);
      free(vertfilepath);
      return NULL;
    } 
  }
  cp = strstr(vertfilepath, ".vert");
  if (cp == NULL) {
    cp = strstr(vertfilepath, ".face");
    if (cp != NULL) {
       strcpy(cp, ".vert");
    } else {
      printf("msmsplugin) file names don't match expected MSMS output\n");
      free(facefilepath);
      free(vertfilepath);
      return NULL;
    } 
  }
 
  ffd = fopen(facefilepath, "r");
  vfd = fopen(vertfilepath, "r");
  if (!ffd || !vfd) { 
    printf("msmsplugin) failed to open either the MSMS face or vertex file\n");
    if (ffd) fclose(ffd);
    if (vfd) fclose(vfd);
    free(facefilepath);
    free(vertfilepath);
    return NULL;
  }
  msms = new msms_t;
  msms->ffd = ffd;
  msms->vfd = vfd;
  msms->graphics = NULL;
  *natoms = 0;
  return msms;
}

static int read_rawgraphics(void *v, int *nelem, 
    const molfile_graphics_t **data) {
  msms_t *msms = (msms_t *)v;
  int i, t;
  float tf;
  int facecount=0;
  int vertexcount=0;
  // count number of faces
  while (fscanf(msms->ffd, "%d %d %d %d %d", &t, &t, &t, &t, &t) == 5) 
    facecount++;
  rewind(msms->ffd);
  // count number of vertices
  while (fscanf(msms->vfd, "%f %f %f %f %f %f %d %d %d", 
         &tf, &tf, &tf, &tf, &tf, &tf, &t, &t, &t) == 9)
    vertexcount++;
  rewind(msms->vfd);

  // simple sanity check to insure we have at least one usable triangle
  if (facecount < 1 || vertexcount < 3) 
    return MOLFILE_ERROR;

  // allocate storage for vertex and normal data
  float *vertex = new float[3 * vertexcount];
  float *normal = new float[3 * vertexcount];

  // read in the vertex data
  for (i=0; i<vertexcount; i++) {
    int addr = i * 3;
    int atomid, l0fa, l;
    fscanf(msms->vfd, "%f %f %f %f %f %f %d %d %d",
           &vertex[addr], &vertex[addr+1], &vertex[addr+2], 
           &normal[addr], &normal[addr+1], &normal[addr+2], &l0fa, &atomid, &l);
  }
 
  // allocate the graphics objects, read in the facet data and 
  // copy the vertex coordinates into triangles as necessary 
  msms->graphics = new molfile_graphics_t[2*facecount];

  // read in the facet data
  for (i=0; i<facecount; i++) {
    int v0, v1, v2, surftype, ana;
    // set the graphics object type
    msms->graphics[2*i    ].type = MOLFILE_TRINORM;  
    msms->graphics[2*i + 1].type = MOLFILE_NORMS;

    // read in the next facet
    fscanf(msms->ffd, "%d %d %d %d %d", &v0, &v1, &v2, &surftype, &ana);
    v0--; // convert from 1-based indexing to 0-based indexing
    v1--;
    v2--;

    // copy the triangle vertices
    float *tri = msms->graphics[2*i    ].data;
    float *nrm = msms->graphics[2*i + 1].data;
    memcpy(tri  , vertex+(3*v0), 3*sizeof(float)); 
    memcpy(tri+3, vertex+(3*v1), 3*sizeof(float)); 
    memcpy(tri+6, vertex+(3*v2), 3*sizeof(float)); 
    memcpy(nrm  , normal+(3*v0), 3*sizeof(float)); 
    memcpy(nrm+3, normal+(3*v1), 3*sizeof(float)); 
    memcpy(nrm+6, normal+(3*v2), 3*sizeof(float)); 
  }

  // set the result array pointers
  *nelem = 2*facecount;
  *data = msms->graphics;

  // delete work area storage
  delete [] normal;
  delete [] vertex;

  return MOLFILE_SUCCESS;
}


static void close_file_read(void *v) {
  msms_t *msms = (msms_t *)v;
  fclose(msms->ffd);
  fclose(msms->vfd);
  delete [] msms->graphics;
  delete msms;
}


/*
 * Initialization stuff here
 */
static molfile_plugin_t plugin = {
  vmdplugin_ABIVERSION,               // ABI version
  MOLFILE_PLUGIN_TYPE,                // type of plugin
  "msms",                             // short name of plugin
  "MSMS Surface Mesh",                // pretty name of plugin
  "John Stone",                       // author
  0,                                  // major version
  2,                                  // minor version
  VMDPLUGIN_THREADSAFE,               // is reentrant
  "face,vert"                         // filename extensions 
};

VMDPLUGIN_EXTERN int VMDPLUGIN_init(void) { return VMDPLUGIN_SUCCESS; }
VMDPLUGIN_EXTERN int VMDPLUGIN_fini(void) { return VMDPLUGIN_SUCCESS; }
VMDPLUGIN_EXTERN int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  plugin.open_file_read = open_file_read;
  plugin.read_rawgraphics = read_rawgraphics;
  plugin.close_file_read = close_file_read;
  (*cb)(v, (vmdplugin_t *)&plugin);
  return VMDPLUGIN_SUCCESS;
}


