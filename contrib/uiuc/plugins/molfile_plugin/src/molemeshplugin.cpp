/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_molemeshplugin
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
 *      $RCSfile: molemeshplugin.C,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.3 $       $Date: 2016/11/28 05:01:54 $
 *
 ***************************************************************************/

/*
 *  quadrilaterals are used to represent the geometry. 
 *
 *  ASCII pmesh file from mole 2.0 follows the following format. Files are case 
insensitive
 *  and whitespace is ignored.
 *
 *  number of vertices        ()
 *  vertex v1x v1y v1z        (<v1> is the first vertex)
 *  vertex v2x v2y v2z        (vertices are given in some order)
 *  vertex v3x v3y v3z        (vertices are given as floating-point values)
 *  vertex v4x v4y v4z        (vertices are given as floating-point values)
 *  number of polygons/faces  
 *  face number               (in groups of 4?)
 */

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <string.h>

#if defined(_AIX)
#include <strings.h>
#endif

#if defined(WIN32) || defined(WIN64)
#define strcasecmp stricmp
#endif
 
#include "molfile_plugin.h"

typedef struct graphics_list {
  molfile_graphics_t gItem;
  struct graphics_list *next;
} molfile_graphics_list;

typedef struct {
  FILE *fd;
  molfile_graphics_t *graphics;
} pmesh_t;

static void *open_file_read(const char *filepath, const char *filetype, int *natoms) {
   FILE *fd;
   pmesh_t *pmesh;

   fd = fopen(filepath, "rb");
   if (!fd) {
     fprintf(stderr, "molemeshplugin) Error opening file.\n");
     return NULL;
   }
   pmesh = new pmesh_t;
   pmesh->fd = fd;
   pmesh->graphics = NULL;
   *natoms = 0;
   return pmesh;
} 
static int read_rawgraphics(void *v, int *nelem, const molfile_graphics_t **data) {
    molfile_graphics_list *gListPtr=NULL, *tmpPtr=NULL;
    int i=0, ntriangles=0,j=0;
    int error=0, numVerts=0, numFacets=0, facetType=0, facet=0, tmpFacet[5];
    pmesh_t *pmesh = (pmesh_t *)v;
    FILE *infile = pmesh->fd;
    char line[81];
    float **tmpData;

    // Check your head(er)
    // "number of vertices"
    fgets(line, 80, infile);
    sscanf(line, "%d", &numVerts);
    if (numVerts  < 1) {
      fprintf(stderr, "molespmeshplugin) error: expected \"Positive Number of Vertices\".\n");
      error = 1;
      return MOLFILE_ERROR;
    } else {
     gListPtr = new molfile_graphics_list;
     gListPtr->next = NULL;
     gListPtr->gItem.type = MOLFILE_TRIANGLE;
     ntriangles++;
     tmpPtr = gListPtr;
   }

//Allocate memory for 2d array.  maybe there is a better way to do this.
    tmpData = new float*[numVerts];
    for (int h=0 ; h < numVerts ; h++ ) {
         tmpData[h] = new float[3];
    }

//we know how many vertices there are read them all into tmpData[j][i]
   for ( j=0 ; j < numVerts; j++) {
// "first loop"
// scan each vertex vx, vy, vz  
    i=0;
    fgets(line, 80, infile);
    float t1=0.0f, t2=0.0f, t3=0.0f;
    if ( sscanf(line, "%f %f %f", &t1, &t2, &t3) == 3 ) {
        tmpData[j][i++] = t1;
        tmpData[j][i++] = t2;
        tmpData[j][i++] = t3;
     } else if(ferror(infile)) {
       fprintf(stderr, "molespmeshplugin) error: problem reading file\n");
       error = 1;
       return MOLFILE_ERROR;
     }
   }
//Read in the total number of facets.  For pmesh files the next number is the total number of triangle or quadrilaterals
//If triangles, then the next number after the total is a 4, and the next four numbers after that describe the vertex index for
//each point of the triangle.  If the set begins with a 5 then the next 5 numbers describe the vertex number for
//each point of the quadrilateral.  So far this reader does not support multiple/concatenated mpesh files yet.
   fgets(line, 80, infile);
   sscanf(line, "%d", &numFacets);
//printf("numFacets %d \n",numFacets);
   if (numFacets  < 1) {
      fprintf(stderr, "molespmeshplugin) error: expected \"Positive Number of Facets\".\n");
      error = 1;
      return MOLFILE_ERROR;
    } else {
     gListPtr = new molfile_graphics_list;

     gListPtr->next = NULL;
     gListPtr->gItem.type = MOLFILE_TRIANGLE;
     ntriangles++;
     tmpPtr = gListPtr;
   }
//Read in the facet type 4=triangle 5=quadrilateral. need to do a while loop here, as long as xxxx<numFacets keep reading
   while ( !feof(infile) && ( error == 0 ) ) {
     fgets(line, 80, infile);
     sscanf(line, "%d", &facetType);
     if (facetType == 4) {
//printf("facetype %d \n", facetType);
        int l=0;
        for (int k=0 ; k < facetType-1; k++) {
              fgets(line, 80, infile);
              sscanf(line, "%d", &facet);
              tmpPtr->gItem.data[l++] = tmpData[facet][0];
              tmpPtr->gItem.data[l++] = tmpData[facet][1];
              tmpPtr->gItem.data[l++] = tmpData[facet][2];
         }
     fgets(line, 80, infile); //one more read to keep us in sync
// Create a new list item and initialize it for second triangle.
         tmpPtr->next = new molfile_graphics_list;
         tmpPtr = tmpPtr->next;
         tmpPtr->next = NULL;
         tmpPtr->gItem.type = MOLFILE_TRIANGLE;
         ntriangles++;
     } else if (facetType == 5 ) { 
                for (int k=0 ; k < facetType-1; k++) {
                     fgets(line, 80, infile);
                     sscanf(line, "%d", &facet);
                     tmpFacet[k]=facet;
                }
              
                tmpPtr->gItem.data[0] = tmpData[tmpFacet[0]][0];
                tmpPtr->gItem.data[1] = tmpData[tmpFacet[0]][1];
                tmpPtr->gItem.data[2] = tmpData[tmpFacet[0]][2];

                tmpPtr->gItem.data[3] = tmpData[tmpFacet[1]][0];
                tmpPtr->gItem.data[4] = tmpData[tmpFacet[1]][1];
                tmpPtr->gItem.data[5] = tmpData[tmpFacet[1]][2];

                tmpPtr->gItem.data[6] = tmpData[tmpFacet[2]][0];
                tmpPtr->gItem.data[7] = tmpData[tmpFacet[2]][1];
                tmpPtr->gItem.data[8] = tmpData[tmpFacet[2]][2];

// Create a new list item and initialize it for second triangle.
                tmpPtr->next = new molfile_graphics_list;
                tmpPtr = tmpPtr->next;
                tmpPtr->next = NULL;
                tmpPtr->gItem.type = MOLFILE_TRIANGLE;
                ntriangles++;

                tmpPtr->gItem.data[0] = tmpData[tmpFacet[0]][0];
                tmpPtr->gItem.data[1] = tmpData[tmpFacet[0]][1];
                tmpPtr->gItem.data[2] = tmpData[tmpFacet[0]][2];

                tmpPtr->gItem.data[3] = tmpData[tmpFacet[2]][0];
                tmpPtr->gItem.data[4] = tmpData[tmpFacet[2]][1];
                tmpPtr->gItem.data[5] = tmpData[tmpFacet[2]][2];

                tmpPtr->gItem.data[6] = tmpData[tmpFacet[3]][0];
                tmpPtr->gItem.data[7] = tmpData[tmpFacet[3]][1];
                tmpPtr->gItem.data[8] = tmpData[tmpFacet[3]][2];

// Create a new list item and initialize it for first triangle of the second quadrilateral.
                tmpPtr->next = new molfile_graphics_list;
                tmpPtr = tmpPtr->next;
                tmpPtr->next = NULL;
                tmpPtr->gItem.type = MOLFILE_TRIANGLE;
                ntriangles++;
     } else if ( (facetType != 4 || facetType != 5) && facetType >= 6 ) {
//Find out if this is a concatenated file by reading the next value and testing for a series of three floats.
//If so, free tmpData and reallocate according to this number. then read all the vertices into tmpData and do
//a drop back into the while loop.
//If not, exit with error.
         fgets(line, 80, infile);
         float t1=0.0f, t2=0.0f, t3=0.0f;
         if ( sscanf(line, "%f %f %f", &t1,&t2,&t3) ==3 ) {
//free tmpData
             for (int x=0; x< 3; x++) free(tmpData[x]);
                  free (tmpData);
             numVerts=facetType;
             facetType=0;
//Allocate new size for vertices array
             tmpData = new float*[numVerts];
             for (int h=0 ; h < numVerts ; h++ ) {
                 tmpData[h] = new float[3];
             }

//Read in all new vertices after adding the first test vertex
             tmpData[0][0] = t1;
             tmpData[0][1] = t2;
             tmpData[0][2] = t3;
             for ( j=1 ; j < numVerts; j++) {
                 i=0;
                 fgets(line, 80, infile);
                 float t1=0.0f, t2=0.0f, t3=0.0f;
                 if ( sscanf(line, "%f %f %f", &t1, &t2, &t3) == 3 ) {
                    tmpData[j][i++] = t1;
                    tmpData[j][i++] = t2;
                    tmpData[j][i++] = t3;
                 } else if(ferror(infile)) {
                           fprintf(stderr, "molespmeshplugin) error: problem reading vertices from concatenated file\n");
                           error = 1;
                           break;
                }
             }
         } else if ( feof(infile)  ) {
//end file read gracefully at the last facet
                    break;
         } else { fprintf(stderr, "molespmeshplugin) error: problem reading concatenated file?\n");
                error = 1;
                break;
         }
         fgets(line, 80, infile);
         sscanf(line, "%d", &numFacets);
         if (numFacets  < 1) {
            fprintf(stderr, "molespmeshplugin) error: expected \"Positive Number of Facets\".\n");
           error = 1;
         }
     }
//go back into the while loop with error=0 so that we can reuse the facet code
//for additional meshes
     error = 0;
    }
    // If an error occurred, free the linked list and return MOLFILE_ERROR
    if (error != 0) {
      while (gListPtr != NULL) {
        tmpPtr = gListPtr->next;
        delete gListPtr;
        gListPtr = tmpPtr;
      }
     for (int x=0; x< 3; x++) free(tmpData[x]);
     free (tmpData);
      return MOLFILE_ERROR;
    }

    // Create the array of molfile_graphics_t, and copy the data from the
    // linked list into it, deleting the list as you go.
    pmesh->graphics = new molfile_graphics_t[ntriangles-1];
//    printf("ntriangles %d \n", ntriangles);
    i = 0;
    while (gListPtr != NULL) {
      pmesh->graphics[i] = gListPtr->gItem;
      tmpPtr = gListPtr->next;
      delete gListPtr;
      gListPtr = tmpPtr;
      i++;
    }

    *nelem = ntriangles-1;
    *data = pmesh->graphics;
     for (int x=0; x< 3; x++) free(tmpData[x]);
     free(tmpData);
    return MOLFILE_SUCCESS;
}

static void close_file_read(void *v) {
  pmesh_t *pmesh = (pmesh_t *)v;
  fclose(pmesh->fd);
  if (pmesh->graphics != NULL)
    delete [] pmesh->graphics;
  delete pmesh;
}


/*
 * Initialization stuff here
 */
static molfile_plugin_t plugin;

VMDPLUGIN_API int VMDPLUGIN_init(void) {
  memset(&plugin, 0, sizeof(molfile_plugin_t));
  plugin.abiversion = vmdplugin_ABIVERSION;
  plugin.type = MOLFILE_PLUGIN_TYPE;
  plugin.name = "pmesh";
  plugin.prettyname = "polygon mesh";
  plugin.author = "Brian Bennion";
  plugin.minorv = 0;
  plugin.majorv = 1;
  plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
  plugin.filename_extension = "mesh";
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

