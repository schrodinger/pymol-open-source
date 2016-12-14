/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_offplugin
#define STATIC_PLUGIN 1

/*
 *  Object File Format (.off)
 *
 *  Represent surfaces composed of polygons. It is documented there:
 *  http://people.sc.fsu.edu/~jburkardt/data/off/off.html
 *  http://shape.cs.princeton.edu/benchmark/documentation/off_format.html
 *
 *  Contributed by Francois-Xavier Coudert (fxcoudert@gmail.com)
 *
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

// internal buffer size used in read_rawgraphics()
#define BUFFLEN 1024

static void *open_file_read(const char *filepath, const char *filetype,
                            int *natoms) {
  FILE *f;
  
  f = fopen(filepath, "rb");
  if (!f) {
    fprintf(stderr, "offplugin) Error opening file.\n");
    return NULL;
  }
  *natoms = 0;
  return f;
}


static char *nextNoncommentLine(char *buff, int bufflen, FILE *f) {
  while (1) {
    char *res = fgets(buff, bufflen, f);
    if (!res || (res[0] != '#' && res[0] != '\n' && res[0] != '\r'))
      return res;
  };
}


static void calcNormals (float vert[9], float norm[9]) {
  float x1 = vert[3] - vert[0];
  float y1 = vert[4] - vert[1];
  float z1 = vert[5] - vert[2];
  float x2 = vert[6] - vert[0];
  float y2 = vert[7] - vert[1];
  float z2 = vert[8] - vert[2];
  float nx = y1 * z2 - z1 * y2;
  float ny = z1 * x2 - x1 * z2;
  float nz = x1 * y2 - y1 * x2;
  float n = 1 / sqrtf(nx*nx+ny*ny+nz*nz);
  norm[0] = norm[3] = norm[6] = n * nx;
  norm[1] = norm[4] = norm[7] = n * ny;
  norm[2] = norm[5] = norm[8] = n * nz;
}


static int read_rawgraphics(void *v, int *nelem, 
                            const molfile_graphics_t **data) {
  int i, k, n;
  int nVert, nFaces, nEdges;
  float *vertices = NULL, *vertColors = NULL;
  char *vertHasColor = NULL;
  molfile_graphics_t *graphics = NULL;
  int j=0;

  char buff[BUFFLEN+1];
  FILE *infile = (FILE *)v;

  // First line is the header: "OFF"
  nextNoncommentLine(buff, BUFFLEN, infile);
  if (buff[0] != 'O' || buff[1] != 'F' || buff[2] != 'F') {
    fprintf(stderr, "offplugin) error: expected \"OFF\" header.\n");
    goto error;
  }

  // Second line: numVertices numFaces numEdges
  nextNoncommentLine(buff, BUFFLEN, infile);
  if (sscanf (buff, " %d %d %d", &nVert, &nFaces, &nEdges) < 2 || 
      nVert <= 0 || nFaces <= 0) {
    fprintf(stderr, "offplugin) error: wrong number of elements.\n");
    goto error;
  }

  // Read vertices
  vertices = (float *) calloc (3 * nVert, sizeof(float));
  vertHasColor = (char *) calloc (nVert, sizeof(char));
  vertColors = (float *) calloc (3 * nVert, sizeof(float));
  for (i = 0; i < nVert; i++) {
    nextNoncommentLine(buff, BUFFLEN, infile);
    int n = sscanf (buff, " %g %g %g %g %g %g", 
                    &vertices[3*i], &vertices[3*i+1], &vertices[3*i+2],
                    &vertColors[3*i], &vertColors[3*i+1], &vertColors[3*i+2]);
    if (n != 3 && n != 6) {
      fprintf(stderr, "offplugin) error: not enough data.\n");
      goto error;
    }
    vertHasColor[i] = (n == 6);
  }

  // Read faces
  // We alloc 6 times the memory because:
  //   -- a quadrangle will be transformed into two triangles.
  //   -- each triangle may have color, and then also its norm will be specified
  graphics = (molfile_graphics_t *) calloc(6*nFaces, sizeof(molfile_graphics_t));
  n = 0;
  for (i = 0; i < nFaces; i++) {
    int idx[4];
    float c[3];
    nextNoncommentLine(buff, BUFFLEN, infile);

    if (sscanf (buff, "%d", &k) != 1 || k < 3) {
      fprintf(stderr, "offplugin) error: not enough data.\n");
      goto error;
    }

    if (k > 4) {
      // TODO -- handle polygon decomposition into triangles
      // Follow the algorithm there:
      // http://www.flipcode.com/archives/Efficient_Polygon_Triangulation.shtml
      fprintf(stderr, "offplugin) error: TODO -- handling polygons with more than 4 vertices.\n");
      goto error;
    }

    if (k == 3) {
      j = sscanf (buff, "%d %d %d %d %g %g %g", &k, &idx[0], &idx[1], &idx[2], &c[0], &c[1], &c[2]);
      bool hasColor = ((j == 7) || (vertHasColor[idx[0]] && vertHasColor[idx[1]] && vertHasColor[idx[2]]));

      graphics[n].type = (hasColor ? MOLFILE_TRICOLOR : MOLFILE_TRIANGLE);
      graphics[n].data[0] = vertices[3*idx[0]  ];
      graphics[n].data[1] = vertices[3*idx[0]+1];
      graphics[n].data[2] = vertices[3*idx[0]+2];
      graphics[n].data[3] = vertices[3*idx[1]  ];
      graphics[n].data[4] = vertices[3*idx[1]+1];
      graphics[n].data[5] = vertices[3*idx[1]+2];
      graphics[n].data[6] = vertices[3*idx[2]  ];
      graphics[n].data[7] = vertices[3*idx[2]+1];
      graphics[n].data[8] = vertices[3*idx[2]+2];
      n++;

      if (j == 7) {
        // The facet has a specific color, use it.
        graphics[n].type = MOLFILE_NORMS;
        calcNormals (graphics[n-1].data, graphics[n].data);
        n++;

        graphics[n].type = MOLFILE_COLOR;
        graphics[n].data[0] = graphics[n].data[3] = graphics[n].data[6] = c[0];
        graphics[n].data[1] = graphics[n].data[4] = graphics[n].data[7] = c[1];
        graphics[n].data[2] = graphics[n].data[5] = graphics[n].data[8] = c[2];
        n++;
      } else if (hasColor) {
        // All three vertices have a color attribute
        graphics[n].type = MOLFILE_NORMS;
        calcNormals (graphics[n-1].data, graphics[n].data);
        n++;

        graphics[n].type = MOLFILE_COLOR;
        graphics[n].data[0] = vertColors[3*idx[0]  ];
        graphics[n].data[1] = vertColors[3*idx[0]+1];
        graphics[n].data[2] = vertColors[3*idx[0]+2];
        graphics[n].data[3] = vertColors[3*idx[1]  ];
        graphics[n].data[4] = vertColors[3*idx[1]+1];
        graphics[n].data[5] = vertColors[3*idx[1]+2];
        graphics[n].data[6] = vertColors[3*idx[2]  ];
        graphics[n].data[7] = vertColors[3*idx[2]+1];
        graphics[n].data[8] = vertColors[3*idx[2]+2];
        n++;
      }
    } else if (k == 4) {
      j = sscanf (buff, "%d %d %d %d %d %g %g %g", &k, &idx[0], &idx[1], &idx[2], &idx[3], &c[0], &c[1], &c[2]);
      bool hasColor = ((j == 8) || (vertHasColor[idx[0]] && vertHasColor[idx[1]] && vertHasColor[idx[2]] && vertHasColor[idx[3]]));

      // Split a quadrangle into two triangles
      graphics[n].type = (hasColor ? MOLFILE_TRICOLOR : MOLFILE_TRIANGLE);
      graphics[n].data[0] = vertices[3*idx[0]  ];
      graphics[n].data[1] = vertices[3*idx[0]+1];
      graphics[n].data[2] = vertices[3*idx[0]+2];
      graphics[n].data[3] = vertices[3*idx[1]  ];
      graphics[n].data[4] = vertices[3*idx[1]+1];
      graphics[n].data[5] = vertices[3*idx[1]+2];
      graphics[n].data[6] = vertices[3*idx[2]  ];
      graphics[n].data[7] = vertices[3*idx[2]+1];
      graphics[n].data[8] = vertices[3*idx[2]+2];
      n++;

      if (j == 8) {
        graphics[n].type = MOLFILE_NORMS;
        calcNormals (graphics[n-1].data, graphics[n].data);
        n++;
      
        graphics[n].type = MOLFILE_COLOR;
        graphics[n].data[0] = graphics[n].data[3] = graphics[n].data[6] = c[0];
        graphics[n].data[1] = graphics[n].data[4] = graphics[n].data[7] = c[1];
        graphics[n].data[2] = graphics[n].data[5] = graphics[n].data[8] = c[2];
        n++;
      } else if (hasColor) {
        graphics[n].type = MOLFILE_NORMS;
        calcNormals (graphics[n-1].data, graphics[n].data);
        n++;

        graphics[n].type = MOLFILE_COLOR;
        graphics[n].data[0] = vertColors[3*idx[0]];
        graphics[n].data[1] = vertColors[3*idx[0]+1];
        graphics[n].data[2] = vertColors[3*idx[0]+2];
        graphics[n].data[3] = vertColors[3*idx[1]];
        graphics[n].data[4] = vertColors[3*idx[1]+1];
        graphics[n].data[5] = vertColors[3*idx[1]+2];
        graphics[n].data[6] = vertColors[3*idx[2]];
        graphics[n].data[7] = vertColors[3*idx[2]+1];
        graphics[n].data[8] = vertColors[3*idx[2]+2];
        n++;
      }

      graphics[n].type = (hasColor ? MOLFILE_TRICOLOR : MOLFILE_TRIANGLE);
      graphics[n].data[0] = vertices[3*idx[2]];
      graphics[n].data[1] = vertices[3*idx[2]+1];
      graphics[n].data[2] = vertices[3*idx[2]+2];
      graphics[n].data[3] = vertices[3*idx[3]];
      graphics[n].data[4] = vertices[3*idx[3]+1];
      graphics[n].data[5] = vertices[3*idx[3]+2];
      graphics[n].data[6] = vertices[3*idx[0]];
      graphics[n].data[7] = vertices[3*idx[0]+1];
      graphics[n].data[8] = vertices[3*idx[0]+2];
      n++;

      if (j == 8) {
        graphics[n].type = MOLFILE_NORMS;
        calcNormals (graphics[n-1].data, graphics[n].data);
        n++;

        graphics[n].type = MOLFILE_COLOR;
        graphics[n].data[0] = graphics[n].data[3] = graphics[n].data[6] = c[0];
        graphics[n].data[1] = graphics[n].data[4] = graphics[n].data[7] = c[1];
        graphics[n].data[2] = graphics[n].data[5] = graphics[n].data[8] = c[2];
        n++;
      } else if (hasColor) {
        graphics[n].type = MOLFILE_NORMS;
        calcNormals (graphics[n-1].data, graphics[n].data);
        n++;

        graphics[n].type = MOLFILE_COLOR;
        graphics[n].data[0] = vertColors[3*idx[2]];
        graphics[n].data[1] = vertColors[3*idx[2]+1];
        graphics[n].data[2] = vertColors[3*idx[2]+2];
        graphics[n].data[3] = vertColors[3*idx[3]];
        graphics[n].data[4] = vertColors[3*idx[3]+1];
        graphics[n].data[5] = vertColors[3*idx[3]+2];
        graphics[n].data[6] = vertColors[3*idx[0]];
        graphics[n].data[7] = vertColors[3*idx[0]+1];
        graphics[n].data[8] = vertColors[3*idx[0]+2];
        n++;
      }
    }
  }

  *nelem = n;
  *data = (molfile_graphics_t *) realloc(graphics, n*sizeof(molfile_graphics_t));
  return MOLFILE_SUCCESS;

  // goto jump target for disaster handling: free memory and bail out
  error:
    free (graphics);
    free (vertices);
    return MOLFILE_ERROR;
}


static void close_file_read(void *v) {
  fclose((FILE *)v);
}


/*
 * Initialization stuff here
 */
static molfile_plugin_t plugin;

VMDPLUGIN_API int VMDPLUGIN_init(void) {
  memset(&plugin, 0, sizeof(molfile_plugin_t));
  plugin.abiversion = vmdplugin_ABIVERSION;
  plugin.type = MOLFILE_PLUGIN_TYPE;
  plugin.name = "off";
  plugin.prettyname = "Object File Format (OFF)";
  plugin.author = "Francois-Xavier Coudert";
  plugin.majorv = 0;
  plugin.minorv = 4;
  plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
  plugin.filename_extension = "off";
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

