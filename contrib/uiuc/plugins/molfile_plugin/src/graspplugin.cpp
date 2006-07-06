/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_graspplugin
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
 *      $RCSfile: graspplugin.C,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.17 $       $Date: 2006/02/23 19:36:44 $
 *
 ***************************************************************************/

/* 
 * Reader for GRASP binary surface files
 * The file format is briefly described at these two sites:
 *   http://honiglab.cpmc.columbia.edu/grasp/grasp_contents.html#A.1
 *   http://www.msg.ucsf.edu/local/programs/grasp/html/Appendix%20A.html
 */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "molfile_plugin.h"
#include "endianswap.h"

static int is_little_endian(void) {
  int x=1;
  return *((char *)&x);
}   

typedef struct {
  FILE *fd;
  molfile_graphics_t *graphics;
} grasp_t;

static void *open_file_read(const char *filepath, const char *filetype,
    int *natoms) {
  FILE *fd;
  grasp_t *grasp;
  
  fd = fopen(filepath, "rb");
  if (!fd) 
    return NULL;
  grasp = new grasp_t;
  grasp->fd = fd;
  grasp->graphics = NULL;
  *natoms = 0;
  return grasp;
}

static int read_rawgraphics(void *v, int *nelem, 
    const molfile_graphics_t **data) {
  grasp_t *grasp = (grasp_t *)v;
  FILE *infile = grasp->fd;

  // Reverse engineering is your friend, and combined with FORTRAN code, voila!
  // od -c shows the header starts off:
  // \0  \0  \0   P   f   o   r   m   a   t   =   2
  // and according to ungrasp, this is a 1 for grasp versions 1.0
  // and 1.1, and 2 for grasp version 1.2
  // Also, the header lines are of length 80 characters + 4 header chars
  // + 4 trailer chars
  // The 4 bytes at the beginning/end are standard Fortran array trash
  char trash[4];
#define TRASH fread(trash, 4, 1, infile)
  char line[81];

  // FIRST LINE OF HEADER; contains format type
  TRASH; 
  fread(line, 1, 80, infile); 
  // make sure it says 'format='
  if (strncmp(line, "format=", 7) != 0) {
    fprintf(stderr, "First characters of file don't look like a GRASP file\n");
    return MOLFILE_ERROR;
  }
  TRASH;
   
  // next char should be a 0 or 1
  char gfiletype = line[7];
  if (gfiletype == '1') {
    gfiletype = 1;
  } else if (gfiletype == '2') {
    gfiletype = 2;
  } else {
    fprintf(stderr, "GRASP file is in format %c, but only '1' or '2' is supported\n", gfiletype);
    return MOLFILE_ERROR;
  }

  // SECOND LINE: contains "vertices,accessibles,normals,triangles"
  TRASH; 
  fread(line, 1, 80, infile); 
  TRASH;

  // THIRD LINE: contains (0 or more of)?
  //  "potentials,curvature,distances,gproperty,g2property,vertexcolor
  char line3[81];
  TRASH; 
  fread(line3, 1, 80, infile); 
  line3[80] = 0;
  TRASH;

  // FOURTH LINE stores vertices, triangles, gridsize, lattice spacing
  int nvert, ntriangles, gridsize;
  float lattice;
  TRASH; 
  fread(line, 1, 80, infile); 
  TRASH;
  sscanf(line, "%d%d%d%f", &nvert, &ntriangles, &gridsize, &lattice);

  // FIFTH LINE stores the center (x,y,z) position
  float center[3];
  TRASH; 
  fread(line, 1, 80, infile); 
  TRASH;
  sscanf(line, "%f%f%f", center, center+1, center+2);

  float *vertex = new float[3 * nvert];
  float *access = new float[3 * nvert];
  float *normal = new float[3 * nvert];
  int *triangle = new int[3 * ntriangles];

  if (!vertex || !access || !normal || !triangle) {
    delete [] vertex;
    delete [] access;
    delete [] normal;
    delete [] triangle;
    fprintf(stderr, "Failed vertex/access/normal/triangle allocations.\n");
    return MOLFILE_ERROR;
  }

  // ungrasp says:
  //    if (filetype.eq.1) then integer*2
  //    if (filetype.eq.2) then integer*4

  // And read them in.  Who needs error checking?
  TRASH; 
  fread(vertex, 3 * sizeof(float), nvert, infile); 
  TRASH;
  TRASH; 
  fread(access, 3 * sizeof(float), nvert, infile); 
  TRASH;
  TRASH; 
  fread(normal, 3 * sizeof(float), nvert, infile); 
  TRASH;

  if (is_little_endian()) {
    swap4_aligned(vertex, 3*nvert);
    swap4_aligned(access, 3*nvert);
    swap4_aligned(normal, 3*nvert);
  }

  if (gfiletype == 2) {
    TRASH; 
    fread(triangle, 3 * sizeof(int), ntriangles, infile); 
    TRASH;
    if (is_little_endian())
      swap4_aligned(triangle, 3*ntriangles);
  } else {
#if 1
    int i;
    short *striangle = new short[3 * ntriangles];
    if (!striangle) {
      delete [] vertex;
      delete [] access;
      delete [] normal;
      delete [] triangle;
      fprintf(stderr, "Failed short triangle allocation.\n");
      return MOLFILE_ERROR;
    }

    TRASH;
    fread(striangle, sizeof(short), 3 * ntriangles, infile);
    if (is_little_endian()) swap2_aligned(striangle, 3 * ntriangles);
    for (i=0; i<3*ntriangles; i++) {
      triangle[i] = striangle[i];
    }
    delete [] striangle;  
    TRASH;
#else
    // do it the slow way (converting from short to int)
    int i;
    short tmp[3];
    TRASH;
    for (i=0; i<ntriangles; i++) {
      fread(tmp, sizeof(short), 3, infile);
      if (is_little_endian()) swap2_aligned(tmp, 3);
      triangle[3*i+0] = tmp[0];
      triangle[3*i+1] = tmp[1];
      triangle[3*i+2] = tmp[2];
    }
    TRASH;
#endif
  }

  // And draw things
  grasp->graphics = new molfile_graphics_t[2*ntriangles];
  int vert1, vert2, vert3;

  for (int tri_count = 0; tri_count < ntriangles; tri_count++) {
    vert1 = triangle[3*tri_count+0] - 1;  // from 1-based FORTRAN
    vert2 = triangle[3*tri_count+1] - 1;  // to 0-based C++
    vert3 = triangle[3*tri_count+2] - 1;

    if (vert1 <      0 || vert2 <      0 || vert3 <      0 ||
        vert1 >= nvert || vert2 >= nvert || vert3 >= nvert) {
      printf("graspplugin) Error, out-of-range vertex index, aborting.\n"); 
      delete [] vertex;
      delete [] access;
      delete [] normal;
      delete [] triangle;
      return MOLFILE_ERROR;
    }

    grasp->graphics[2*tri_count  ].type = MOLFILE_TRINORM;
    grasp->graphics[2*tri_count+1].type = MOLFILE_NORMS;

    float *tridata =  grasp->graphics[2*tri_count  ].data;
    float *normdata = grasp->graphics[2*tri_count+1].data;
    memcpy(tridata  , vertex+3*vert1, 3*sizeof(float));
    memcpy(tridata+3, vertex+3*vert2, 3*sizeof(float));
    memcpy(tridata+6, vertex+3*vert3, 3*sizeof(float));
    memcpy(normdata  , normal+3*vert1, 3*sizeof(float));
    memcpy(normdata+3, normal+3*vert2, 3*sizeof(float));
    memcpy(normdata+6, normal+3*vert3, 3*sizeof(float));
  } 

  *nelem = 2*ntriangles;
  *data = grasp->graphics;

  delete [] triangle;
  delete [] normal;
  delete [] access;
  delete [] vertex;

  return MOLFILE_SUCCESS;
}


static void close_file_read(void *v) {
  grasp_t *grasp = (grasp_t *)v;
  fclose(grasp->fd);
  delete [] grasp->graphics;
  delete grasp;
}


/*
 * Initialization stuff here
 */
static molfile_plugin_t plugin = {
  vmdplugin_ABIVERSION,               // ABI version
  MOLFILE_PLUGIN_TYPE,                // type of plugin
  "grasp",                            // short name of plugin
  "GRASP",                            // pretty name of plugin
  "Justin Gullingsrud, John Stone",   // author
  0,                                  // major version
  5,                                  // minor version
  VMDPLUGIN_THREADSAFE,               // is reentrant
  "srf,grasp"                         // filename extension 
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


