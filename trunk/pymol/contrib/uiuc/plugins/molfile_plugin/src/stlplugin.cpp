/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_stlplugin
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
 *      $RCSfile: stlplugin.C,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.11 $       $Date: 2006/02/23 19:36:45 $
 *
 ***************************************************************************/

/*
 *  STL files are used by most Automatic Fabricators (3D printers). Only 
 *  triangles are used to represent the geometry. 
 *
 *  ASCII STL files follow the following format. Files are case insensitive
 *  and whitespace is ignored.
 *
 *  solid name                ("name" is an arbitrary label for the solid)
 *  facet normal ni nj nk     (<n> is a unit normal of the triangle)
 *  outer loop
 *  vertex v1x v1y v1z        (<v1> is the first vertex)
 *  vertex v2x v2y v2z        (vertices are given in anti-clockwise order)
 *  vertex v3x v3y v3z        (vertices are given as floating-point values)
 *  endloop                   (there is no space in the label "endloop")
 *  endfacet                  (likewise)
 *  ...                       (additional facets are given as above)
 *  endsolid name             (this ends the stl file, same name as above)
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
} stl_t;

static void *open_file_read(const char *filepath, const char *filetype,
    int *natoms) {
  FILE *fd;
  stl_t *stl;
  
  fd = fopen(filepath, "rb");
  if (!fd) {
    fprintf(stderr, "Error opening file.\n");
    return NULL;
  }
  stl = new stl_t;
  stl->fd = fd;
  stl->graphics = NULL;
  *natoms = 0;
  return stl;
}

static int read_rawgraphics(void *v, int *nelem, 
    const molfile_graphics_t **data) {
    molfile_graphics_list *gListPtr=NULL, *tmpPtr=NULL;
    int i=0, ntriangles=0;
    int error=0;
    stl_t *stl = (stl_t *)v;
    FILE *infile = stl->fd;
    char line[81], keyWord[81];

    // Check your head(er)
    // "solid name"
    fgets(line, 80, infile);
    sscanf(line, " %s", keyWord);
    if (strcasecmp(keyWord, "solid") != 0)
    {
      fprintf(stderr, "STL file error: expected \"solid\".\n");
      error = 1;
    }

    // "facet normal ni nj nk"
    fgets(line, 80, infile);
    sscanf(line, " %s", keyWord);
    if (strcasecmp(keyWord, "facet") != 0)
    {
      fprintf(stderr, "STL file error: expected \"facet\".\n");
      error = 1;
    }
    else
    {
      gListPtr = new molfile_graphics_list;
      gListPtr->next = NULL;
      gListPtr->gItem.type = MOLFILE_TRIANGLE;
      ntriangles++;
      tmpPtr = gListPtr;
    }

    while ( !feof(infile) && (error == 0) )
    {
      // "outer loop"
      fgets(line, 80, infile);
      sscanf(line, " %s", keyWord);
      if (strcasecmp(keyWord, "outer") != 0)
      {
        fprintf(stderr, "STL file error: expected \"outer\".\n");
        error = 1;
        break;
      }
      else
      {
        i = 0;
      }
        
      // "vertex vx, vy, vz"
      while (i < 9)
      {
        fgets(line, 80, infile);
        sscanf(line, " %s", keyWord);
        if (strcasecmp(keyWord, "vertex") != 0)
        {
          fprintf(stderr, "STL file error: expected \"vertex\".\n");
          error = 1;
          break; 
        }
        else if ( sscanf(line, " %*s %f %f %f", &(tmpPtr->gItem.data[i++]),
                         &(tmpPtr->gItem.data[i++]), 
                         &(tmpPtr->gItem.data[i++])) != 3 )
        {
          fprintf(stderr, "STL file error: not enough vertices.\n");
          error = 1;
          break;
        }
      }
      if (error != 0) break;

      // "endloop"
      fgets(line, 80, infile);
      sscanf(line, " %s", keyWord);
      if (strcasecmp(keyWord, "endloop") != 0)
      {
        fprintf(stderr, "STL file error: expected \"endloop\".\n");
        error = 1;
        break;
      }
      
      // "endfacet"
      fgets(line, 80, infile);
      sscanf(line, " %s", keyWord);
      if (strcasecmp(keyWord, "endfacet") != 0)
      {
        fprintf(stderr, "STL file error: expected \"endfacet\".\n");
        error = 1;
        break;
      }

      // "endsolid" or "facet normal ni nj nk"
      fgets(line, 80, infile);
      sscanf(line, " %s", keyWord);
      if (strcasecmp(keyWord, "endsolid") == 0)
      {
        break;
      }
      if (strcasecmp(keyWord, "facet") == 0)
      {
        // Create a new list item and initialize it.
        tmpPtr->next = new molfile_graphics_list;
        tmpPtr = tmpPtr->next;
        tmpPtr->next = NULL;
        tmpPtr->gItem.type = MOLFILE_TRIANGLE;
        ntriangles++;
      }
      else
      {
        fprintf(stderr, 
                "STL file error: expected \"facet\" or \"endsolid\".\n");
        error = 1;
        break;
      }

      // file error
      if(ferror(infile))
      {
        fprintf(stderr, "STL file error: problem reading file\n");
        error = 1;
        break;
      }
    }


    // If an error occurred, free the linked list and return MOLFILE_ERROR
    if (error != 0)
    {
      while (gListPtr != NULL)
      {
        tmpPtr = gListPtr->next;
        delete gListPtr;
        gListPtr = tmpPtr;
      }
      return MOLFILE_ERROR;
    }

    // Create the array of molfile_graphics_t, and copy the data from the
    // linked list into it, deleting the list as you go.
    stl->graphics = new molfile_graphics_t[ntriangles];
    i = 0;
    while (gListPtr != NULL)
    {
      stl->graphics[i] = gListPtr->gItem;
      tmpPtr = gListPtr->next;
      delete gListPtr;
      gListPtr = tmpPtr;
      i++;
    }

    *nelem = ntriangles;
    *data = stl->graphics;

    return MOLFILE_SUCCESS;
}

static void close_file_read(void *v) {
  stl_t *stl = (stl_t *)v;
  fclose(stl->fd);
  if (stl->graphics != NULL)
    delete [] stl->graphics;
  delete stl;
}


/*
 * Initialization stuff here
 */
static molfile_plugin_t plugin = {
  vmdplugin_ABIVERSION,   // ABI version
  MOLFILE_PLUGIN_TYPE, 	  // type of plugin
  "stl",                  // name of plugin
  "STL Stereolithography Triangle Mesh", // pretty name of plugin
  "Eamon Caddigan",       // author
  0,                      // major version
  1,                      // minor version
  VMDPLUGIN_THREADSAFE,   // is reentrant
  "stl"                   // filename extension 
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


