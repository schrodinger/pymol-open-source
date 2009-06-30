/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_raster3dplugin
#define STATIC_PLUGIN 1


/***************************************************************************
 *cr
 *cr            (C) Copyright 1995-2009 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "molfile_plugin.h"

typedef struct {
  FILE *fd;
  molfile_graphics_t *graphics;
} handle_t;

static void *open_file_read(const char *filepath, const char *filetype,
    int *natoms) {
  FILE *fd;
  handle_t *handle;
  
  fd = fopen(filepath, "rt");
  if (!fd) 
    return NULL;
  handle = new handle_t;
  handle->fd = fd;
  handle->graphics = NULL;
  *natoms = 0;
  return handle;
}

static void next_elem(int &n, int &max, molfile_graphics_t *& graphics) {
  ++n;
  if (n == max) {
    max *= 2;
    graphics = (molfile_graphics_t *)realloc(graphics, max*sizeof(molfile_graphics_t));
  }
}

// gets a line of text from the file, ignoring comments
static int get_line(int &line, char *buf, int len, FILE *f) {
   do {
      line++;
      if (!fgets(buf, len - 1, f)) return 0;
   } while (buf[0] == '#');
   return 1;
}

static int read_rawgraphics(void *v, int *nelem, 
    const molfile_graphics_t **gdata) {

/* TODO
 *
 * 1. doesn't properly render capped cylinders
 * 2. needs to handle object types 15, 19, 9 (affects coordinate system);
 *    needs to first understand VMD's coordinate system.
 * 3. doesn't properly handle object type 8 (material definitions)
 * 4. doesn't properly handle object type 14 (quadrics)
 * 5. doesn't properly support file indirection (using '@')
 * 6. doesn't differentiate between round-ended cylinders (obj 3)
 *    and flat-ended cylinders (obj 5)
 * 7. doesn't handle planes
 * 8. doesn't support text (object types 10-12)
 * 9. doesn't handle depth cueing and transformations in global
 *    properties (object 16)
 */

////  The header contains a lot of information which I ignore
// TITLE (80 chars)
// NTX, NTY (number of tiles -- ignored)
// NPX, NPY (number of points per tile -- ignored)
// SCHEME (antialiasing scheme -- ignored)
// BKGND (background color -- ignored)
// SHADOW (T/F -- ignored)
// IPHONG (phong power -- ignored)
// STRAIT (secondary light source -- ignored)
// AMBIEN (ambient illumination -- ignored)
// SPECLR ignored
// EYEPOS ignored
// SOURCE ignored
// TMAT 4x4 matrix used to view the system
// object mode (1 == triangle, 2 == spheres, or 3 ==  mixed)
// INFMTS  --  FORTRAN input field specifications
//   for input format specifiers 1/6, 2, and 3/5
//   triangle(1) or plane(6) (x1,y1,z1)-(x2,y2,z2)-(x3,y3,z3) (r,g,b)
//   sphere(2)      (x,y,z) r (r,g,b)
//   cylinder(3,5)  (x1,y1,z1) R1 - (x2,y2,z2) R2(is ignored) (r,g,b)
//                      except if one radius is 0, a cone is made
//   I ignore these lines and just read them in as "%f"

   int futureVersion = 0;
   int line = 0;

   int count, i;
   char buffer[200];
   float mat[16];
   FILE *infile;

   handle_t *handle = (handle_t *)v;
   infile = handle->fd;

   int maxelem = 10;
   int n = 0;
   molfile_graphics_t *graphics = (molfile_graphics_t *)malloc(
       maxelem*sizeof(molfile_graphics_t));

   // XXX This should probably be in open_file_read
   if (!get_line(line, buffer, 199, infile)) {
      fprintf(stderr, "raster3dplugin) Error reading file header (line %d)\n",
          line);
      return MOLFILE_ERROR;
   }

   // tell the user info about the file
   for (i = strlen(buffer) - 1; i >= 0 &&
        (buffer[i] == 10 || buffer[i] == 13); i--) buffer[i] = 0;
   printf("raster3dplugin) scene title: '%s'\n", buffer);

   // The next 11 lines of text contain more header information
   // about lighting, anti aliasing, image size, etc. This can be
   // ignored.
   for (count = 0; count < 11; count++) {
      if (!get_line(line, buffer, 199, infile)) {
         fprintf(stderr, 
             "raster3dplugin) error reading file header (line %d)\n", line);
         return MOLFILE_ERROR;
      }
   }

   // Now I have to get the matrix.  This is made nasty since
   // there could be extra text after the first four numbers on a line
   // as in: 1 0 0 0 This is an extra comment
   for (i=0; i<4; i++) {
      get_line(line, buffer, 199, infile);  // read the whole line into a string
      if (sscanf(buffer, "%f %f %f %f",
		 &mat[4*i], &mat[4*i+1], &mat[4*i+2], &mat[4*i+3])<4) {
         fprintf(stderr, "raster3dplugin) invalid format in file (line %d)\n",
             line);
	 return MOLFILE_ERROR;
      }
   }

   get_line(line, buffer, 199, infile);
   if (sscanf(buffer, "%d", &i) < 1) {
     fprintf(stderr, 
         "raster3dplugin) error reading object input mode (line %d)\n", line);
     return MOLFILE_ERROR;
   }

   if (i != 3) {
      fprintf(stderr, 
          "raster3dplugin) the specified file is in an unsupported format\n");
      fprintf(stderr, 
          "(object input mode %d). Aborting.\n", i);
      return MOLFILE_ERROR;
   }

   float data[15];
   float normals[15];
   float tricolors[9];
   float color[9];

      // INFMT/INFMTS input specifiers; these can specifiy Fortran formatted
      // input formats, but are usually "*" meaning free-format input. We
      // give a warning if it's not free-format.
      for (count = 0; count < 3; count++) {
         get_line(line, buffer, 199, infile);
         for (i = strlen(buffer) - 1; i >= 0 && (buffer[i] == 10 || buffer[i] == 13); i--)
            buffer[i] = 0;
         if (strcmp(buffer, "*")) break;
      }
      if (count < 3) {  
        fprintf(stderr, "raster3dplugin) Warning: this file contains input in a nonstandard\n");
        fprintf(stderr, "Fortran format. This is generally not supported, and the read may fail.\n");
      }
      count = 0;

      while (!feof(infile) && !ferror(infile)) {
	 int objtype = -1;

         if (!get_line(line, buffer, 199, infile)) continue;

         if (sscanf(buffer, "%d", &objtype) != 1) {
            fprintf(stderr, "raster3dplugin) bad data in file (line %d)\n",
                line);
            return MOLFILE_ERROR;
         }

	 switch(objtype) {

	    case 1: // triangle
               char buffer2[200];
               int have_normals;
               int have_tricolors;
               long fpos;

               have_normals = 0;
               have_tricolors = 0;

               get_line(line, buffer, 127, infile);
               if (feof(infile)) {
                  //msgErr << "Raster3D input: error reading triangle data (line "
                         //<< line << ")" << sendmsg;
                  return MOLFILE_ERROR;
               }

               if (sscanf(buffer, "%f %f %f %f %f %f %f %f %f %f %f %f",
                          data  , data+1, data+2,
                          data+3, data+4, data+5,
                          data+6, data+7, data+8,
                          data+9, data+10, data+11) < 12) { 
                  //msgErr << "Raster3D input: bad triangle data in file (line "
                         //<< line << "). Will try to continue." << sendmsg;
                  continue;
               }

               while (!feof(infile) && !ferror(infile)) {

                  fpos = ftell(infile);
		  if (!get_line(line, buffer2, 199, infile)) {
                     fseek(infile, fpos, SEEK_SET);
                     break;
                  }

		  if (sscanf(buffer2, "%d", &objtype) != 1) {
                     //msgErr << "Raster3D input: bad data in file (line " << line
                            //<< "). Aborting." << sendmsg;
                     return MOLFILE_ERROR;
                  }

                  switch (objtype) {
                     case 7: // explicit normals

                        if (!get_line(line, buffer2, 199, infile)) {
                           //msgErr << "Raster3D input: read error in file (line " << line
                                  //<< "). Aborting." << sendmsg;
                           return MOLFILE_ERROR;
                        }

                        if (sscanf(buffer2, "%f %f %f %f %f %f %f %f %f",
                                   normals  , normals+1, normals+2,
                                   normals+3, normals+4, normals+5,
                                   normals+6, normals+7, normals+8 ) < 9) { 
                           //msgErr << "Raster3D input: error reading triangle normals (line "
                                  //<< line << "). Will try to continue." << sendmsg;
                           continue;
                        }

                        have_normals = 1;

                        break;

                     case 17: // colors at vertices of a tri-color

                        if (!get_line(line, buffer2, 199, infile)) {
                           //msgErr << "Raster3D input: read error in file (line " << line
                                  //<< "). Aborting." << sendmsg;
                           return MOLFILE_ERROR;
                        }

                        if (sscanf(buffer2, "%f %f %f %f %f %f %f %f %f",
                                   tricolors,  tricolors+ 1, tricolors + 2,
                                   tricolors + 3, tricolors + 4, tricolors + 5,
                                   tricolors + 6, tricolors + 7, tricolors + 8) < 9) {
                           //msgErr << "Raster3D input: error reading vertex colors (line "
                                  //<< line << "). Will try to continue." << sendmsg;
                           continue;
                        }

                        have_tricolors = 1;

                        break;

                     default:

                        fseek(infile, fpos, SEEK_SET);
                        fpos = 0;
                        break;

                  }

                  if (!fpos) break;
               }

               if (ferror(infile)) {
                  //msgErr << "Raster3D input: read error in file (line "
                           //<< line << "). Aborting." << sendmsg;
                  return MOLFILE_ERROR;
               }

                            graphics[n].type = MOLFILE_COLOR;
                            graphics[n].data[0] = sqrt(data[9]);
                            graphics[n].data[1] = sqrt(data[10]);
                            graphics[n].data[2] = sqrt(data[11]);
                            next_elem(n, maxelem, graphics);

               if (have_tricolors) {
                 for (int qq=0; qq<9; qq++) color[qq] = sqrt(tricolors[qq]);
               }

               if (!have_normals && !have_tricolors) {
                 graphics[n].type = MOLFILE_TRIANGLE;
                 memcpy(graphics[n].data, data, 3*sizeof(float));
                 next_elem(n, maxelem, graphics);

	       } else if (have_normals && !have_tricolors) {
           graphics[n].type = MOLFILE_TRINORM;
           memcpy(graphics[n].data, data, 9*sizeof(float));
           next_elem(n, maxelem, graphics);
           graphics[n].type = MOLFILE_NORMS;
           memcpy(graphics[n].data, normals, 9*sizeof(float));
           next_elem(n, maxelem, graphics);
               } else if (have_tricolors && !have_normals) {
#if 0

                  float tmp1[3], tmp2[3], tmp3[3];
                  int j;
                  for (j = 0; j < 3; j++) {
                     tmp1[j] = data[3 + j] - data[j];
                     tmp2[j] = data[6 + j] - data[3 + j];
                  }
                  cross_prod(tmp3, tmp1, tmp2);
                  vec_normalize(tmp3);
                  triclr.putdata(data, data+3, data+6,
                                 tmp3, tmp3, tmp3,
                                 colorIndices[0], colorIndices[1], colorIndices[2], cmdList);
#else
                  // XXX Take the average of the color
                  graphics[n].type = MOLFILE_COLOR;
                  graphics[n].data[0] = (color[0]+color[3]+color[6])/3.0f;
                  graphics[n].data[1] = (color[1]+color[4]+color[7])/3.0f;
                  graphics[n].data[2] = (color[2]+color[5]+color[8])/3.0f;
                  next_elem(n, maxelem, graphics);
                  graphics[n].type = MOLFILE_TRIANGLE;
                  memcpy(graphics[n].data, data, 3*sizeof(float));
                  next_elem(n, maxelem, graphics);

#endif
               } else {

#if 0
                  triclr.putdata(data, data+3, data+6,
                                 normals, normals+3, normals+6,
                                 colorIndices[0], colorIndices[1], colorIndices[2], cmdList);
#else
           graphics[n].type = MOLFILE_TRICOLOR;
           memcpy(graphics[n].data, data, 9*sizeof(float));
           next_elem(n, maxelem, graphics);
           graphics[n].type = MOLFILE_NORMS;
           memcpy(graphics[n].data, normals, 9*sizeof(float));
           next_elem(n, maxelem, graphics);
           graphics[n].type = MOLFILE_COLOR;
           memcpy(graphics[n].data, color, 9*sizeof(float));
           next_elem(n, maxelem, graphics);
#endif
               }

	       break;

	    case 2: // sphere

               if (!get_line(line, buffer, 199, infile)) {
                  //msgErr << "Raster3D input: error reading sphere data (line "
                         //<< line << ")" << sendmsg;
                  return MOLFILE_ERROR;
               }

               if (sscanf(buffer, "%f %f %f %f %f %f %f",
                          data, data + 1, data + 2, data + 3,
                          data + 4, data + 5, data + 6) != 7) {
                  //msgErr << "Raster3D input: bad sphere data in file (line "
                         //<< line << "). Will try to continue." << sendmsg;
                  continue;
               }
                 graphics[n].type = MOLFILE_COLOR;
                 color[0] = sqrt(data[4]);
                 color[1] = sqrt(data[5]);
                 color[2] = sqrt(data[6]);
                 memcpy(graphics[n].data, color, 3*sizeof(float));
                 next_elem(n, maxelem, graphics);
                 graphics[n].type = MOLFILE_SPHERE;
                 graphics[n].size = data[3];
                 graphics[n].style = 12; // XXX hard-coded resolution
                 memcpy(graphics[n].data, data, 3*sizeof(float));
                 next_elem(n, maxelem, graphics);
	       break;

	    case 3: // rounded cylinder
	    case 5: // flat cylinder

               if (!get_line(line, buffer, 199, infile)) {
                  //msgErr << "Raster3D input: error reading cylinder data (line "
                         //<< line << ")" << sendmsg;
                  return MOLFILE_ERROR;
               }

               if (sscanf(buffer, "%f %f %f %f %f %f %f %f %f %f %f",
                          data, data + 1, data + 2, data + 3, data + 4,
                          data + 5, data + 6, data + 7, data + 8, data + 9,
                          data + 10) != 11) {
                  //msgErr << "Raster3D input: bad cylinder data (line "
                         //<< line << "). Will try to continue." << sendmsg;
                  continue;
               }
                 graphics[n].type = MOLFILE_COLOR;
                 color[0] = sqrt(data[8]);
                 color[1] = sqrt(data[9]);
                 color[2] = sqrt(data[10]);
                 memcpy(graphics[n].data, color, 3*sizeof(float));
                 next_elem(n, maxelem, graphics);
                 graphics[n].type = MOLFILE_CYLINDER;
                 graphics[n].size = data[3];
                 graphics[n].style = 12; // XXX hard-coded resolution
                 memcpy(graphics[n].data, data, 3*sizeof(float));
                 memcpy(graphics[n].data+3, data+4, 3*sizeof(float));
                 next_elem(n, maxelem, graphics);

               // cap with spheres for object 3
               if (objtype == 3) {
                 graphics[n].type = MOLFILE_SPHERE;
                 graphics[n].size = data[3];
                 graphics[n].style = 12; // XXX hard-coded resolution
                 memcpy(graphics[n].data, data, 3*sizeof(float));
                 next_elem(n, maxelem, graphics);
                 graphics[n].type = MOLFILE_SPHERE;
                 graphics[n].size = data[3];
                 graphics[n].style = 12; // XXX hard-coded resolution
                 memcpy(graphics[n].data, data+4, 3*sizeof(float));
                 next_elem(n, maxelem, graphics);
               }

	       break;

            case 9:
               break;

            case 6: case 8: case 10: case 11: case 12: case 13: case 15: case 19:

               // Ignore the next line
	       get_line(line, buffer, 199, infile);

               break;

            case 7:
               //msgErr << "Raster3D input: encountered unexpected object 7 (triangle ";
               //msgErr << "vertex normals) (line " << line;
               //msgErr << "). Will try to continue." << sendmsg;
               break;

            case 17:
               //msgErr << "Raster3D input: encountered unexpected object 17 (triangle " << sendmsg;
               //msgErr << "vertex colors) (line " << line << "). Will try to continue." << sendmsg;
               break;

            case 0:  // Raster3D 'force EOF' -- checked later
               break;

            // We encountered an object that is not recognized.
            // Need to warn the user. Future version of R3d perhaps?
	    default:
               if (!futureVersion) {
                  //msgErr << "Raster3D input: encountered unknown object type #"
                         //<< objtype << sendmsg;
                  //msgErr << " (line " << line << "). Future version of Raster3D maybe?"
                         //<< " Will try to continue." << sendmsg;
                  futureVersion = 1;
               }
	       break;

         } // end of switch

         // If this is a Raster3d "force EOF" object, break
         if (objtype == 0) break;

      } // end of while

   if (ferror(infile)) {
      //msgErr << "Raster3D input: read error while reading input file (line "
             //<< line << "). Aborting." << sendmsg;
      return MOLFILE_ERROR;
   }

   // normal exit
   *nelem = n;
   handle->graphics = graphics;
   *gdata = graphics;
   return MOLFILE_SUCCESS;
}

static void close_file_read(void *v) {
  handle_t *handle = (handle_t *)v;
  fclose(handle->fd);
  handle->fd = NULL;
  delete [] handle->graphics;
  handle->graphics = NULL;
  delete handle;
}


/*
 * Initialization stuff here
 */
static molfile_plugin_t plugin;

VMDPLUGIN_EXTERN int VMDPLUGIN_init(void) {
  memset(&plugin, 0, sizeof(molfile_plugin_t));
  plugin.abiversion = vmdplugin_ABIVERSION;
  plugin.type = MOLFILE_PLUGIN_TYPE;
  plugin.name = "raster3d";
  plugin.prettyname = "Raster3d Scene File";
  plugin.author = "Justin Gullingsrud";
  plugin.majorv = 0;
  plugin.minorv = 2;
  plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
  plugin.filename_extension = "r3d";
  plugin.open_file_read = open_file_read;
  plugin.read_rawgraphics = read_rawgraphics;
  plugin.close_file_read = close_file_read;
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_EXTERN int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_EXTERN int VMDPLUGIN_fini(void) { return VMDPLUGIN_SUCCESS; }

