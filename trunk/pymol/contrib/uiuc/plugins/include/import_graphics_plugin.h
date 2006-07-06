/***************************************************************************
 *cr
 *cr            (C) Copyright 1995-2002 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: import_graphics_plugin.h,v $
 *      $Author: justin $       $Locker:  $             $State: Exp $
 *      $Revision: 1.4 $       $Date: 2002/07/29 03:45:51 $
 *
 ***************************************************************************/

#ifndef IMPORT_GRAPHICS_PLUGIN_H
#define IMPORT_GRAPHICS_PLUGIN_H

/* 
 * API for C extensions to define a way to import low-level graphics primitives 
 */ 

#include "vmdplugin.h"

/*
 * Define a common plugin type to be used when registering the plugin.
 */
#define IMPORT_GRAPHICS_PLUGIN_TYPE "import graphics"

typedef enum {
	IMPORT_GRAPHICS_LINE_SOLID, IMPORT_GRAPHICS_LINE_DASHED
} import_graphics_linestyle_t;

/* 
 * Application-provided callbacks for specifying graphics primitives.  
 * Items must be maintained in order by the application for the purpose of 
 * coloring; see below.
 */
typedef struct {

  /*
   * Draw a point at the specified location in 3-D space.
   */
  int (* add_point)(void *, const float *x);
  int (* add_triangle)(void *, const float *x1, const float *x2, const float *x3);
  int (* add_trinorm)(void *, const float *x1, const float *x2, const float *x3,
		 const float *n1, const float *n2, const float *n3);
  int (* add_line)(void *, const float *x, const float *y, int line_style, 
		  int width); 
  int (* add_cylinder)(void *, const float *x, const float *y, float radius,
		  int resolution, int filled);
  int (* add_sphere)(void *, const float *x, float rad, int resolution);
  int (* add_text)(void *, const float *x, const char *text, float size);
  /*
   * Color to use for subsequent primitives.  If primitives are added before
   * any call to use_color, the application is free to do whatever it likes.
   */
  int (* use_color)(void *, float r, float g, float b);

  /*
   * Indicate whether the set of primitives is to be lit or not.  Either all
   * or none of the primitives will be lit.
   */ 
  int (* use_materials)(void *, int yes_no);
} import_graphics_cb_t;


/*
 * Main file reader API begins here.  Any function in this struct may be NULL
 * if not implemented by the plugin; the application checks this to determine
 * what functionality is present in the plugin. 
 */ 
typedef struct {
  /*
   * Required header
   */
  vmdplugin_HEAD

  /*
   * Filename extension for this file type.  May be NULL if no filename 
   * extension exists and/or is known.
   */
  const char *filename_extension;

  /* 
   * Try to open the file for reading.  Return an opaque handle, or NULL on
   * failure. filetype should be the name under which this plugin was 
   * registered; this is provided so that plugins can provide the same 
   * function pointer * to handle multiple file types.
   */
  void *(* open_file_read)(const char *filepath, const char *filetype); 
  
  /*
   * Read data and return it to the application in the supplied
   * callbacks.  The first void * is an opaque application handle which 
   * should be passed to all the callbacks in import_cb_t.  The second 
   * void * is the plugin handle returned by open_file_read.  
   */
  int (* read_data)(void *, void *mydata, import_graphics_cb_t *);

  /* 
   * Close the file and release all data.  The handle cannot be reused.
   */
  void (* close_file_read)(void *);

} import_graphics_plugin_t;

#endif

