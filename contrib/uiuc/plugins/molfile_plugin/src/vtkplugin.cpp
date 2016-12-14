/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_vtkplugin
#define STATIC_PLUGIN 1

/***************************************************************************
 *cr
 *cr            (C) Copyright 2007-2011 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: vtkplugin.C,v $
 *      $Author: johns $        $Locker:  $             $State: Exp $
 *      $Revision: 1.5 $      $Date: 2015/10/26 22:14:45 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *   Plugin for reading uniform grids and vector fields written 
 *   in the VTK ASCII format
 ***************************************************************************/

//
// Prototype file reader developed using VTK refs and a couple samples:
//   http://www.vtk.org/Wiki/VTK/Tutorials/3DDataTypes
//   http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
//

// This parser currently expects the data to be formatted like this:
// # vtk DataFile Version X.0
// titletext
// ASCII
// DATASET STRUCTURED_POINTS
// DIMENSIONS XXX YYY ZZZ
// SPACING XX YY ZZ
// ORIGIN XX YY XX
// POINT_DATA n
//
// And either a field like this:
// FIELD fieldname numarrays
// arrayname0 numcomponents numtuples datatype
// val0 val1 ...
//
// Or a list of vectors like this:
// VECTORS dataname datatype
// val0 val1 ...
//

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <string.h>

#if defined(_AIX)
#include <strings.h>
#endif

#if defined(WIN32) || defined(WIN64)
#define strcasecmp  stricmp
#define strncasecmp strnicmp
#endif

#include "molfile_plugin.h"
#include "largefiles.h"

#define THISPLUGIN plugin
#include "vmdconio.h"

#define LINESIZE 2040

typedef struct {
  FILE *fd;
  char title[257]; /// spec says 256 chars max w/ newline termination
  int nsets;
  molfile_volumetric_t *vol;
  int isBinary; 
} vtk_t;


// Get a string from a stream, printing any errors that occur
static char *vtkgets(char *s, int n, FILE *stream) {
  char *returnVal;

  if (feof(stream)) {
    printf("vtkplugin) Unexpected end-of-file.\n");
    return NULL;
  } else if (ferror(stream)) {
    printf("vtkplugin) Error reading file.\n");
    return NULL;
  } else {
    returnVal = fgets(s, n, stream);
    if (returnVal == NULL) {
      printf("vtkplugin) Error reading line.\n");
    }
  }

  return returnVal;
}


static int vtkgetstrcmp(char *s, int n, FILE *stream, const char *cmpstr) {
  char *str = vtkgets(s, n, stream);
  int rc = strncmp(cmpstr, str, strlen(cmpstr));
  if (rc) {
    printf("vtkplugin) found '%s', expected '%s'\n", str, cmpstr);
  }
  return rc;
}
   

static void *open_vtk_read(const char *filepath, const char *filetype,
    int *natoms) {
  FILE *fd;
  vtk_t *vtk;
  char inbuf[LINESIZE];
  int xsize, ysize, zsize;
  float orig[3], xdelta[3], ydelta[3], zdelta[3];
 
  memset(orig, 0, sizeof(orig));
  memset(xdelta, 0, sizeof(xdelta));
  memset(ydelta, 0, sizeof(ydelta));
  memset(zdelta, 0, sizeof(zdelta));
 
  fd = fopen(filepath, "rb");
  if (!fd) {
    printf("vtkplugin) Error opening file.\n");
    return NULL;
  }

  // allocate and initialize the vtk structure
  vtk = new vtk_t;
  memset(vtk, 0, sizeof(vtk_t));
  vtk->fd = fd;
  vtk->vol = NULL;
  vtk->isBinary = 0;
  *natoms = MOLFILE_NUMATOMS_NONE;
  vtk->nsets = 1; /* this file contains only one data set */

  /* skip comments */
  do {
    if (vtkgets(inbuf, LINESIZE, fd) == NULL) 
      return NULL;
  } while (inbuf[0] == '#');

  // read VTK title line 
  printf("vtkplugin) Dataset title: '%s'\n", inbuf);
  strncpy(vtk->title, inbuf, sizeof(vtk->title) - 1);
  vtk->title[256]='\0'; // force-terminate potentially truncated title string

  if (vtkgetstrcmp(inbuf, LINESIZE, fd, "ASCII")) return NULL;
  if (vtkgetstrcmp(inbuf, LINESIZE, fd, "DATASET STRUCTURED_POINTS")) return NULL;

  // get the grid dimensions
  if (vtkgets(inbuf, LINESIZE, fd) == NULL) {
    delete vtk;
    return NULL;
  }
  if (sscanf(inbuf, "DIMENSIONS %d %d %d", &xsize, &ysize, &zsize) != 3) {
    printf("vtkplugin) Error reading grid dimensions!\n");
    delete vtk;
    return NULL;
  }

  // get the grid spacing
  if (vtkgets(inbuf, LINESIZE, fd) == NULL) {
    delete vtk;
    return NULL;
  }
  if (sscanf(inbuf, "SPACING %e %e %e", xdelta, ydelta+1, zdelta+2) != 3) {
    printf("vtkplugin) Error reading cell dimensions!\n");
    delete vtk;
    return NULL;
  }

  // get the grid origin
  if (vtkgets(inbuf, LINESIZE, fd) == NULL) {
    delete vtk;
    return NULL;
  }
  if (sscanf(inbuf, "ORIGIN %e %e %e", orig, orig+1, orig+2) != 3) {
    printf("vtkplugin) Error reading grid origin!\n");
    delete vtk;
    return NULL;
  }

  // get number of grid points
  if (vtkgets(inbuf, LINESIZE, fd) == NULL) {
    delete vtk;
    return NULL;
  }
  int numgridpoints = 0;
  if (sscanf(inbuf, "POINT_DATA %d", &numgridpoints) != 1) {
    printf("vtkplugin) Error reading grid point counts!\n");
    delete vtk;
    return NULL;
  }

  // get field or vector list depending on file format variant we have
  if (vtkgets(inbuf, LINESIZE, fd) == NULL) {
    delete vtk;
    return NULL;
  }

  char tmp[256];
  sscanf(inbuf, "%s", tmp);
  if (!strcmp(tmp, "FIELD")) {
    char fieldname[256];
    int numarrays=0;
    sscanf(inbuf, "FIELD %s %d", fieldname, &numarrays);
    printf("vtkplugin) FIELD: name '%s', %d arrays\n", fieldname, numarrays);

    // eat the array name for field zero
    if (vtkgets(inbuf, LINESIZE, fd) == NULL) {
      delete vtk;
      return NULL;
    }
  } else if (!strcmp(tmp, "VECTORS")) {
    // prepare to eat vectors
    char fieldname[256];
    int numvecs=0;
    sscanf(inbuf, "VECTORS %s %d", fieldname, &numvecs);
    printf("vtkplugin) VECTORS: name '%s', %d arrays\n", fieldname, numvecs);
  } else {
    printf("vtkplugin) Unrecognized file structure, aborting!:\n");
    printf("vtkplugin) line contents: '%s'\n", inbuf);
    delete vtk;
    return NULL;
  } 

  vtk->vol = new molfile_volumetric_t[1];
  memset(vtk->vol, 0, sizeof(molfile_volumetric_t));
  strcpy(vtk->vol[0].dataname, "VTK volumetric map");

  /* Set the unit cell origin and basis vectors */
  for (int i=0; i<3; i++) {
    vtk->vol[0].origin[i] = orig[i];
    vtk->vol[0].xaxis[i] = xdelta[i] * ((xsize-1 > 0) ? (xsize-1) : 1);
    vtk->vol[0].yaxis[i] = ydelta[i] * ((ysize-1 > 0) ? (ysize-1) : 1);
    vtk->vol[0].zaxis[i] = zdelta[i] * ((zsize-1 > 0) ? (zsize-1) : 1);
  }

  vtk->vol[0].xsize = xsize;
  vtk->vol[0].ysize = ysize;
  vtk->vol[0].zsize = zsize;

#if vmdplugin_ABIVERSION > 16
  vtk->vol[0].has_scalar = 1;
  vtk->vol[0].has_gradient = 1;
  vtk->vol[0].has_variance = 0;
#endif
  vtk->vol[0].has_color = 0; // no color data

  return vtk;
}


static int read_vtk_metadata(void *v, int *nsets, 
  molfile_volumetric_t **metadata) {
  vtk_t *vtk = (vtk_t *)v;
  *nsets = vtk->nsets; 
  *metadata = vtk->vol;  

  return MOLFILE_SUCCESS;
}


static int read_vtk_data(void *v, int set, float *datablock,
                         float *colorblock) {
  vtk_t *vtk = (vtk_t *)v;
  FILE *fd = vtk->fd;
  int x, y, z, xsize, ysize, zsize, xysize, total;

  if (vtk->isBinary)
    return MOLFILE_ERROR;

  xsize = vtk->vol[0].xsize;
  ysize = vtk->vol[0].ysize;
  zsize = vtk->vol[0].zsize;
  xysize = xsize * ysize;
  total = xysize * zsize;

  double scalemag = 1.0;
  const char *userscalefactor=getenv("VMDVTKPLUGINSCALEVOXELMAG");
  if (userscalefactor) {
    scalemag = atof(userscalefactor);
    if (scalemag != 0.0) {
      printf("vtkplugin) Applying user scaling factor to voxel scalar/gradient values: %g\n", scalemag); 
    } else {
      printf("vtkplugin) Warning: ignoring user scaling factor due to parse error or zero-value\n"); 
    }
  } else {
    printf("vtkplugin) No user scaling factor set, using scale factor 1.0.\n");
  }

  float maxmag = 0.0f;
  strcpy(vtk->vol[0].dataname, "volgradient");
  for (z=0; z<zsize; z++) {
    for (y=0; y<ysize; y++) {
      for (x=0; x<xsize; x++) {
        double vx, vy, vz;
        fscanf(fd, "%lf %lf %lf", &vx, &vy, &vz);

#if 1
        // XXX hack to allow user override of voxel magnitude
        // during loading...
        vx *= scalemag;
        vy *= scalemag;
        vz *= scalemag;
#endif

        // compute scalar magnitude from vector field
        double mag = sqrt(vx*vx + vy*vy + vz*vz);

        int addr = z*xsize*ysize + y*xsize + x;
        datablock[addr] = mag;

        if (mag > maxmag)
          maxmag = mag;
      }
    }
  }
  printf("vtkplugin) maxmag: %g\n", maxmag);

  return MOLFILE_SUCCESS;
}


#if vmdplugin_ABIVERSION > 16

static int read_vtk_data_ex(void *v, molfile_volumetric_readwrite_t *p) {
  vtk_t *vtk = (vtk_t *)v;
  FILE *fd = vtk->fd;
  int x, y, z, xsize, ysize, zsize, xysize, total;

  if (vtk->isBinary)
    return MOLFILE_ERROR;

  if (!p->scalar || !p->gradient) 
    return MOLFILE_ERROR;

  xsize = vtk->vol[0].xsize;
  ysize = vtk->vol[0].ysize;
  zsize = vtk->vol[0].zsize;
  xysize = xsize * ysize;
  total = xysize * zsize;

  double scalemag = 1.0;
  const char *userscalefactor=getenv("VMDVTKPLUGINSCALEVOXELMAG");
  if (userscalefactor) {
    scalemag = atof(userscalefactor);
    if (scalemag != 0.0) {
      printf("vtkplugin) Applying user scaling factor to voxel scalar/gradient values: %g\n", scalemag); 
    } else {
      printf("vtkplugin) Warning: ignoring user scaling factor due to parse error or zero-value\n"); 
    }
  } else {
    printf("vtkplugin) No user scaling factor set, using scale factor 1.0.\n");
  }

  float maxmag = 0.0f;
  strcpy(vtk->vol[0].dataname, "volgradient");
  for (z=0; z<zsize; z++) {
    for (y=0; y<ysize; y++) {
      for (x=0; x<xsize; x++) {
        double vx, vy, vz;
        fscanf(fd, "%lf %lf %lf", &vx, &vy, &vz);

#if 1
        // XXX hack to allow user override of voxel magnitude
        // during loading...
        vx *= scalemag;
        vy *= scalemag;
        vz *= scalemag;
#endif

        // compute scalar magnitude from vector field
        double mag = sqrt(vx*vx + vy*vy + vz*vz);

        int addr = z*xsize*ysize + y*xsize + x;
        p->scalar[addr] = mag;

        if (mag > maxmag)
          maxmag = mag;

        // assign vector field to gradient map
        // index into vector field of 3-component vectors
        int addr3 = addr *= 3;
        p->gradient[addr3    ] = vx;
        p->gradient[addr3 + 1] = vy;
        p->gradient[addr3 + 2] = vz;
      }
    }
  }
  printf("vtkplugin) maxmag: %g\n", maxmag);

  return MOLFILE_SUCCESS;
}

#endif


static void close_vtk_read(void *v) {
  vtk_t *vtk = (vtk_t *)v;
  
  fclose(vtk->fd);
  if (vtk->vol != NULL)
    delete [] vtk->vol; 
  delete vtk;
}


//
// Initialization stuff here
//
static molfile_plugin_t vtkplugin;

VMDPLUGIN_API int VMDPLUGIN_init(void) {
  memset(&vtkplugin, 0, sizeof(molfile_plugin_t));
  vtkplugin.abiversion = vmdplugin_ABIVERSION;
  vtkplugin.type = MOLFILE_PLUGIN_TYPE;
  vtkplugin.name = "vtk";
  vtkplugin.prettyname = "VTK grid reader";
  vtkplugin.author = "John Stone";
  vtkplugin.majorv = 0;
  vtkplugin.minorv = 2;
  vtkplugin.is_reentrant = VMDPLUGIN_THREADUNSAFE;
  vtkplugin.filename_extension = "vtk";
  vtkplugin.open_file_read = open_vtk_read;
  vtkplugin.read_volumetric_metadata = read_vtk_metadata;
  vtkplugin.read_volumetric_data = read_vtk_data;
#if vmdplugin_ABIVERSION > 16
  vtkplugin.read_volumetric_data_ex = read_vtk_data_ex;
#endif
  vtkplugin.close_file_read = close_vtk_read;
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)&vtkplugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_fini(void) { return VMDPLUGIN_SUCCESS; }


