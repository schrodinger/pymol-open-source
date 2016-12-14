/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_plyplugin
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
 *      $RCSfile: plyplugin.C,v $
 *      $Author: johns $       $:  $             $State: Exp $
 *      $Revision: 1.8 $       $Date: 2016/11/28 05:01:54 $
 *
 ***************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "molfile_plugin.h"
#include "endianswap.h"

#include "ply.h"
#include "ply_c.h"

/* vertex and face definitions for a polygonal object */

typedef struct Vertex {
  float x,y,z;
  float r,g,b;
  float nx,ny,nz;
  void *other_props;       /* other properties */
} Vertex;

typedef struct Face {
  unsigned char nverts;    /* number of vertex indices in list */
  int *verts;              /* vertex index list */
  void *other_props;       /* other properties */
} Face;

PlyProperty vert_props[] = { /* list of property information for a vertex */
  {"x", Float32, Float32, offsetof(Vertex,x), 0, 0, 0, 0},
  {"y", Float32, Float32, offsetof(Vertex,y), 0, 0, 0, 0},
  {"z", Float32, Float32, offsetof(Vertex,z), 0, 0, 0, 0},
  {"r", Float32, Float32, offsetof(Vertex,r), 0, 0, 0, 0},
  {"g", Float32, Float32, offsetof(Vertex,g), 0, 0, 0, 0},
  {"b", Float32, Float32, offsetof(Vertex,b), 0, 0, 0, 0},
  {"nx", Float32, Float32, offsetof(Vertex,nx), 0, 0, 0, 0},
  {"ny", Float32, Float32, offsetof(Vertex,ny), 0, 0, 0, 0},
  {"nz", Float32, Float32, offsetof(Vertex,nz), 0, 0, 0, 0},
};

PlyProperty face_props[] = { /* list of property information for a face */
  {"vertex_indices", Int32, Int32, offsetof(Face,verts),
   1, Uint8, Uint8, offsetof(Face,nverts)},
  {"vertex_index", Int32, Int32, offsetof(Face,verts),
   1, Uint8, Uint8, offsetof(Face,nverts)},
};



/// plugin data handle
typedef struct {
  FILE *fd;
  molfile_graphics_t *graphics;

  int per_vertex_color;
  int has_normals;
} ply_t;


static void *open_file_read(const char *filepath, const char *filetype,
                            int *natoms) {
  FILE *fd;
  ply_t *ply;
  
  printf("plyplugin) Opening PLY file '%s'\n", filepath);
  fd = fopen(filepath, "rb");
  if (!fd) 
    return NULL;
  ply = new ply_t;
  ply->fd = fd;
  ply->graphics = NULL;
  *natoms = 0;
  return ply;
}


static int read_rawgraphics(void *v, int *nelem, 
                            const molfile_graphics_t **data) {
  ply_t *ply = (ply_t *)v;
  ply->per_vertex_color = 0;
  ply->has_normals = 0;

  int i=0;
  int nverts=0;
  int nfaces=0;
  char *elem_name=NULL;
  Vertex **vlist=NULL;
  Face **flist=NULL;
  PlyOtherProp *vert_other=NULL;
  PlyOtherProp *face_other=NULL;

  printf("plyplugin) Reading PLY file header...\n");
  PlyFile *in_ply = read_ply(ply->fd);

  printf("plyplugin) Processing PLY contents...\n");
  printf("plyplugin) num_elem_types: %d\n", in_ply->num_elem_types);  

  for (i=0; i<in_ply->num_elem_types; i++) {
    int elem_count = 0;

    /* prepare to read the i'th list of elements */
    elem_name = setup_element_read_ply (in_ply, i, &elem_count);

    if (equal_strings ("vertex", elem_name)) {
      int j=0;
      nverts = elem_count;
      printf("plyplugin) reading %d vertex elements...\n", nverts);

      /* create a vertex list to hold all the vertices */
      vlist = (Vertex **) calloc (1, sizeof(Vertex *)*nverts);

      /* set up for getting vertex elements */
      setup_property_ply(in_ply, &vert_props[0]);
      setup_property_ply(in_ply, &vert_props[1]);
      setup_property_ply(in_ply, &vert_props[2]);

      for (j=0; j<in_ply->elems[i]->nprops; j++) {
        PlyProperty *prop;
        prop = in_ply->elems[i]->props[j];
        if (equal_strings ("r", prop->name)) {
          setup_property_ply(in_ply, &vert_props[3]);
          ply->per_vertex_color = 1;
        }
        if (equal_strings ("g", prop->name)) {
          setup_property_ply(in_ply, &vert_props[4]);
          ply->per_vertex_color = 1;
        }
        if (equal_strings ("b", prop->name)) {
          setup_property_ply(in_ply, &vert_props[5]);
          ply->per_vertex_color = 1;
        }
        if (equal_strings ("nx", prop->name)) {
          setup_property_ply(in_ply, &vert_props[6]);
          ply->has_normals = 1;
        }
        if (equal_strings ("ny", prop->name)) {
          setup_property_ply(in_ply, &vert_props[7]);
          ply->has_normals = 1;
        }
        if (equal_strings ("nz", prop->name)) {
          setup_property_ply(in_ply, &vert_props[8]);
          ply->has_normals = 1;
        }
      }

      vert_other = get_other_properties_ply(in_ply,
                                            offsetof(Vertex,other_props));

      /* grab all the vertex elements */
      for (j=0; j<nverts; j++) {
        vlist[j] = (Vertex *) calloc(1, sizeof(Vertex));
        vlist[j]->r = 1;
        vlist[j]->g = 1;
        vlist[j]->b = 1;
        get_element_ply (in_ply, (void *) vlist[j]);
      }

    } else if (equal_strings ("face", elem_name)) {
      int j=0;
      nfaces = elem_count;
      printf("plyplugin) reading %d face elements...\n", nfaces);

      /* create a list to hold all the face elements */
      flist = (Face **) calloc(1, sizeof(Face *)*nfaces);

      /* set up for getting face elements */
      for (j=0; j<in_ply->elems[i]->nprops; j++) {
        PlyProperty *prop;
        prop = in_ply->elems[i]->props[j];
        if (equal_strings ("vertex_indices", prop->name)) {
          setup_property_ply(in_ply, &face_props[0]);
        }
        if (equal_strings ("vertex_index", prop->name)) {
          setup_property_ply(in_ply, &face_props[1]);
        }
      }

      face_other = get_other_properties_ply (in_ply,
                                             offsetof(Face,other_props));

      /* grab all the face elements */
      for (j=0; j<nfaces; j++) {
        flist[j] = (Face *) calloc(1, sizeof(Face));
        get_element_ply(in_ply, (void *) flist[j]);
      }
    } else {
      printf("plyplugin) reading other elements...\n");
      get_other_element_ply(in_ply);
    } 
  }

  printf("plyplugin) freeing PLY structures\n");
  free_ply(in_ply);
  in_ply = NULL;

  printf("plyplugin) generating %d graphics primitives...\n", nfaces); 
  ply->graphics = new molfile_graphics_t[2*nfaces];
  int vert1, vert2, vert3;

  for (i=0; i<nfaces; i++) {
    if (flist[i]->nverts != 3) {
      printf("plyplugin) Found non-triangle facet, aborting.\n");
      return MOLFILE_ERROR;
    }
    vert1 = flist[i]->verts[0];
    vert2 = flist[i]->verts[1];
    vert3 = flist[i]->verts[2];

    if (vert1 <      0  || vert2 <      0  || vert3 <       0 ||
        vert1 >= nverts || vert2 >= nverts || vert3 >= nverts) {
      printf("plyplugin) Error, out-of-range vertex index, aborting.\n"); 
      return MOLFILE_ERROR;
    }

    ply->graphics[i].type = MOLFILE_TRIANGLE;
    float *tridata =  ply->graphics[i].data;
    tridata[0] = vlist[vert1]->x;
    tridata[1] = vlist[vert1]->y;
    tridata[2] = vlist[vert1]->z;
    tridata[3] = vlist[vert2]->x;
    tridata[4] = vlist[vert2]->y;
    tridata[5] = vlist[vert2]->z;
    tridata[6] = vlist[vert3]->x;
    tridata[7] = vlist[vert3]->y;
    tridata[8] = vlist[vert3]->z;
  } 

  *nelem = nfaces;
  *data = ply->graphics;

  printf("plyplugin) freeing ply face list\n");
  for (i=0; i<nfaces; i++) {
    free(flist[i]);
  }
  memset(flist, 0, sizeof(Face *)*nfaces);
  free(flist);
  flist = NULL;

  printf("plyplugin) freeing ply vertex list\n");
  for (i=0; i<nverts; i++) {
    free(vlist[i]);
  }
  memset(vlist, 0, sizeof(float *)*nverts);
  free(vlist);
  vlist=NULL;

  return MOLFILE_SUCCESS;
}


static void close_file_read(void *v) {
  ply_t *ply = (ply_t *)v;
  // close_ply(in_ply);
  fclose(ply->fd);
  
  delete [] ply->graphics;
  delete ply;
}


/*
 * Initialization stuff here
 */
static molfile_plugin_t plugin;

VMDPLUGIN_API int VMDPLUGIN_init(void) {
  memset(&plugin, 0, sizeof(molfile_plugin_t));
  plugin.abiversion = vmdplugin_ABIVERSION;
  plugin.type = MOLFILE_PLUGIN_TYPE;
  plugin.name = "ply";
  plugin.prettyname = "PLY";
  plugin.author = "John Stone";
  plugin.majorv = 0;
  plugin.minorv = 2;
  plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
  plugin.filename_extension = "ply";
  plugin.open_file_read = open_file_read;
  plugin.read_rawgraphics = read_rawgraphics;
  plugin.close_file_read = close_file_read;
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_fini(void) { return VMDPLUGIN_SUCCESS; }



