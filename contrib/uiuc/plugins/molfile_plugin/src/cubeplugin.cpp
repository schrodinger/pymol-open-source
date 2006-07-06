/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_cubeplugin
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
 *      $RCSfile: cubeplugin.C,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.23 $       $Date: 2006/01/05 00:05:52 $
 *
 ***************************************************************************/

//
// Plugin reader for Gaussian "cube" files
//
// Gaussian "cube" file format described here:
//   http://www.gaussian.com/00000430.htm
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "molfile_plugin.h"

#ifndef M_PI_2
#define M_PI_2 1.57079632679489661922
#endif

#include "periodic_table.h"

static const float bohr = 0.529177249;

// A format-independent structure to hold unit cell data
typedef struct {
  float A, B, C, alpha, beta, gamma;
} cube_box;

typedef struct {
  FILE *fd;              // file descriptor
  int nsets;             // number of volume datasets
  int numatoms;          // number of atoms
  bool coord;            // has coordinate data
  long crdpos, datapos;  // seek offsets for coords and data
  char *file_name;       // original filename 
  float *datacache;      // temporary cache of orbital data prior to conversion
  molfile_volumetric_t *vol; // volume data
  float origin[3];       // origin, stored for periodic display hack 
  float rotmat[3][3];    // rotation matrix, stored for periodic display hack
  cube_box box;          // unit cell dimensions
} cube_t;

// Converts box basis vectors to A, B, C, alpha, beta, and gamma.  
// Stores values in cube_box struct, which should be allocated before calling
// this function.
static int cube_readbox(cube_box *box, float *x, float *y, float *z) {
  float A, B, C;

  if (!box) {
    return 1;
  }

  // provide defaults
  box->A = 10.0;
  box->B = 10.0;
  box->C = 10.0;
  box->alpha = 90.0;
  box->beta  = 90.0;
  box->gamma = 90.0;

  // A, B, C are the lengths of the x, y, z vectors, respectively
  A = sqrt( x[0]*x[0] + x[1]*x[1] + x[2]*x[2] );
  B = sqrt( y[0]*y[0] + y[1]*y[1] + y[2]*y[2] );
  C = sqrt( z[0]*z[0] + z[1]*z[1] + z[2]*z[2] );
  if ((A<=0) || (B<=0) || (C<=0)) {
    return 1;
  }
  box->A = A;
  box->B = B;
  box->C = C;

  // gamma, beta, alpha are the angles between the x & y, x & z, y & z
  // vectors, respectively
  box->gamma = acos( (x[0]*y[0]+x[1]*y[1]+x[2]*y[2])/(A*B) ) * 90.0/M_PI_2;
  box->beta = acos( (x[0]*z[0]+x[1]*z[1]+x[2]*z[2])/(A*C) ) * 90.0/M_PI_2;
  box->alpha = acos( (y[0]*z[0]+y[1]*z[1]+y[2]*z[2])/(B*C) ) * 90.0/M_PI_2; 

  return 0;
}

// calculate and store origin and rotation matrix to realign everything later.
static void cube_buildrotmat(cube_t *cube, float *o, float *a, float *b)
{
  // we rotate first around y and z to align a along the x-axis...
  const double len   = sqrt(a[0]*a[0] + a[1]*a[1]);
  const double phi   = atan2((double) a[2], (double) len);
  const double theta = atan2((double) a[1], (double) a[0]);

  const double cph = cos(phi);
  const double cth = cos(theta);
  const double sph = sin(phi);
  const double sth = sin(theta);

  // ...then we rotate around x to put b into the xy-plane.
  const double psi = atan2(-sph*cth*b[0] - sph*sth*b[1] + cph*b[2],-sth*b[0] + cth*b[1]);
  const double cps = cos(psi);
  const double sps = sin(psi);

  const double r[3][3] = { 
    {               cph*cth,                    cph*sth,      sph},
    {-sth*cps - sph*cth*sps,      cth*cps - sph*sth*sps,  cph*sps},
    { sth*sps - sph*cth*cps,     -cth*sps - sph*sth*cps,  cph*cps}
  };

  for (int i=0; i<3; ++i) {
    cube->origin[i] = o[i];
    for (int j=0; j<3; ++j) {
      cube->rotmat[i][j] = r[i][j];
    }
  }
}


static void eatline(FILE * fd) {
  char readbuf[1025];
  fgets(readbuf, 1024, fd);    // go on to next line
}  

// prototype.
static void close_cube_read(void *v);

static void *open_cube_read(const char *filepath, const char *filetype,
    int *natoms) {
  FILE *fd;
  cube_t *cube;
  int xsize, ysize, zsize;
  float a[3], b[3], c[3];
  int i;
 
  fd = fopen(filepath, "rb");
  if (!fd) 
    return NULL;

  cube = new cube_t;
  cube->fd = fd;
  cube->vol = NULL;
  cube->coord = false;
  cube->file_name = strdup(filepath);
  cube->datacache = NULL;

  // initialize origin and rotmat to sensible defaults.
  for (i=0; i<3; ++i) {
      for (int j=0; j<3; ++j) {
          cube->rotmat[i][j] = 0.0;
      }
  }
  for (i=0; i<3; ++i) {
      cube->origin[i] = 0.0;
      cube->rotmat[i][i] = 1.0;
  }

  molfile_volumetric_t voltmpl; // base information for all data sets.

  // read in cube file header information
  char readbuf[256]; 
  fgets(readbuf, 256, cube->fd);    // go on to next line

  // identify this file, and read title string into dataset info
  strcpy(voltmpl.dataname, "Gaussian Cube: ");
  strncat(voltmpl.dataname, readbuf, 240);      // 240 is max space left after
                                                //   "Gaussian Cube: "
  eatline(cube->fd);          // skip second header line 

  // read in number of atoms
  if (fscanf(cube->fd, "%d", &cube->numatoms) != 1) {
    close_cube_read(cube);
    return NULL;
  }

  if (cube->numatoms > 0) {   // density cube file
    cube->nsets = 1;          // this cube file contains only one data set
  } else {
    // cube file with orbitals => multiple densities.
    cube->numatoms = - cube->numatoms;
    cube->nsets = 0;          // we don't know yet how many data sets.
  }
  *natoms = cube->numatoms;

  // read in cube origin
  if ( fscanf(cube->fd, "%f %f %f", 
         &voltmpl.origin[0], 
         &voltmpl.origin[1], 
         &voltmpl.origin[2]) != 3 ) {
    close_cube_read(cube);
    return NULL;
  }

  // read in cube axes and sizes
  if ((fscanf(cube->fd, "%d", &xsize) != 1) ||
      (fscanf(cube->fd, "%f %f %f", &a[0], &a[1], &a[2]) != 3)) {
    close_cube_read(cube);
    return NULL;
  }

  if ((fscanf(cube->fd, "%d", &ysize) != 1) ||
      (fscanf(cube->fd, "%f %f %f", &b[0], &b[1], &b[2]) != 3)) {
    close_cube_read(cube);
    return NULL;
  }

  if ((fscanf(cube->fd, "%d", &zsize) != 1) ||
      (fscanf(cube->fd, "%f %f %f", &c[0], &c[1], &c[2]) != 3)) {
    close_cube_read(cube);
    return NULL;
  }

  // calculate number of samples in each dimension
  voltmpl.xsize = xsize;
  voltmpl.ysize = ysize;
  voltmpl.zsize = zsize;
  voltmpl.has_color = 0;

  eatline(cube->fd);     // skip remaining end of line characters

  // to make the periodic display work, we need to rotate
  // the cellvectors (and the coordinates) in such a way,
  // that the a-vector is collinear with the x-axis and
  // the b-vector is in the xy-plane. 
  cube_buildrotmat(cube, voltmpl.origin, a, b);
  // print warning, if the rotation will be significant:
  if ((fabs((double) a[1]) + fabs((double) a[2]) + fabs((double) b[2]))
      > 0.001) {
    fprintf(stderr, "cubeplugin) WARNING: Coordinates will be rotated to comply \n"
            "cubeplugin) with VMD's conventions for periodic display...\n");
  }
  

  // all dimensional units are always in Bohr
  // so scale axes and origin correctly.
  // NOTE: the angstroms are only allowed in input.
  voltmpl.origin[0] *= bohr; 
  voltmpl.origin[1] *= bohr;
  voltmpl.origin[2] *= bohr;

  // store aligned axes.
  for (i=0; i<3; ++i) {
    voltmpl.xaxis[i] = cube->rotmat[i][0] * a[0] 
      + cube->rotmat[i][1] * a[1] + cube->rotmat[i][2] * a[2];

    voltmpl.yaxis[i] = cube->rotmat[i][0] * b[0] 
      + cube->rotmat[i][1] * b[1] + cube->rotmat[i][2] * b[2];
    
    voltmpl.zaxis[i] = cube->rotmat[i][0] * c[0] 
      + cube->rotmat[i][1] * c[1] + cube->rotmat[i][2] * c[2];
  }

  voltmpl.xaxis[0] *= bohr * xsize;
  voltmpl.xaxis[1] *= bohr * xsize; 
  voltmpl.xaxis[2] *= bohr * xsize; 

  voltmpl.yaxis[0] *= bohr * ysize; 
  voltmpl.yaxis[1] *= bohr * ysize; 
  voltmpl.yaxis[2] *= bohr * ysize; 

  voltmpl.zaxis[0] *= bohr * zsize; 
  voltmpl.zaxis[1] *= bohr * zsize; 
  voltmpl.zaxis[2] *= bohr * zsize; 

#if 1
  /*   As of VMD version 1.8.3, volumetric data points are 
   *   expected to represent the center of a grid box. cube format 
   *   volumetric data represents the value at the edges of the 
   *   grid boxes, so we need to shift the internal origin by half
   *   a grid box diagonal to have the data at the correct position.
   *   This will need to be changed again when the plugin interface 
   *   is updated to explicitly allow point/face-centered data sets.
   */
  voltmpl.origin[0] -= 0.5 * ( voltmpl.xaxis[0] / (double) xsize
			+ voltmpl.yaxis[0] / (double) ysize
			+ voltmpl.zaxis[0] / (double) zsize);
  voltmpl.origin[1] -= 0.5 * ( voltmpl.xaxis[1] / (double) xsize
			+ voltmpl.yaxis[1] / (double) ysize
			+ voltmpl.zaxis[1] / (double) zsize);
  voltmpl.origin[2] -= 0.5 * ( voltmpl.xaxis[2] / (double) xsize
			+ voltmpl.yaxis[2] / (double) ysize
			+ voltmpl.zaxis[2] / (double) zsize);
#endif

#if defined(TEST_PLUGIN)
  printf("cell before rotation:\n");
  printf("a: %12.8f %12.8f %12.8f\n", a[0]*xsize*bohr, a[1]*ysize*bohr, a[2]*zsize*bohr);
  printf("b: %12.8f %12.8f %12.8f\n", b[0]*xsize*bohr, b[1]*ysize*bohr, b[2]*zsize*bohr);
  printf("c: %12.8f %12.8f %12.8f\n", c[0]*xsize*bohr, c[1]*ysize*bohr, c[2]*zsize*bohr);

  printf("cell after rotation:\n");
  printf("x: %12.8f %12.8f %12.8f\n", voltmpl.xaxis[0], voltmpl.xaxis[1], voltmpl.xaxis[2]);
  printf("y: %12.8f %12.8f %12.8f\n", voltmpl.yaxis[0], voltmpl.yaxis[1], voltmpl.yaxis[2]);
  printf("z: %12.8f %12.8f %12.8f\n", voltmpl.zaxis[0], voltmpl.zaxis[1], voltmpl.zaxis[2]);
#endif

  // store the unit cell information for later perusal.
  if (cube_readbox(&(cube->box), voltmpl.xaxis, voltmpl.yaxis, voltmpl.zaxis)) {
	  fprintf(stderr, "cubeplugin) Calculation of unit cell size failed. Trying to continue...\n");
  }

  cube->crdpos = ftell(cube->fd); // and record beginning of coordinates
  // XXX fseek()/ftell() are incompatible with 64-bit LFS I/O implementations, 
  // hope we don't read any files >= 2GB...

  if (cube->nsets >0) { 
    int i;

    // density cube file, copy voltmpl into the cube struct.
    cube->vol = new molfile_volumetric_t[1];
    memcpy(cube->vol, &voltmpl, sizeof(voltmpl));

    // skip over coordinates to find the start of volumetric data
    for (i=0; i < cube->numatoms; i++) {
      eatline(cube->fd);
    }

    cube->datapos = ftell(cube->fd); // and record beginning of data
    // XXX fseek()/ftell() are incompatible with 64-bit LFS I/O, 
    // hope we don't read any files >= 2GB...
  } else {              
    int i;

    // orbital cube file. we now have to read the orbitals section
    // skip over coordinates
    for (i=0; i < cube->numatoms; i++) {
      eatline(cube->fd);
    }
      
    fscanf(cube->fd, "%d", &cube->nsets);
    fprintf(stderr, "\ncubeplugin) found %d orbitals\n", cube->nsets);
    cube->vol = new molfile_volumetric_t[cube->nsets];
    for (i=0; i < cube->nsets; ++i) {
      int orb;
      fscanf(cube->fd, "%d", &orb);
      memcpy(&cube->vol[i], &voltmpl, sizeof(voltmpl));
      sprintf(cube->vol[i].dataname, "Gaussian Cube: Orbital %d", orb);
    }
      
    eatline(cube->fd);        // gobble up rest of line
    cube->datapos = ftell(cube->fd); // and record beginning of data
    // XXX fseek()/ftell() are incompatible with 64-bit LFS I/O, 
    // hope we don't read any files >= 2GB...
  }

  return cube;
}

  
static int read_cube_structure(void *v, int *optflags, molfile_atom_t *atoms) {
  int i, j;
  char *k;
  molfile_atom_t *atom;

  cube_t *cube = (cube_t *)v;

  // go to coordinates
  fseek(cube->fd, cube->crdpos, SEEK_SET);
  // XXX fseek()/ftell() are incompatible with 64-bit LFS I/O implementations, 
  // hope we don't read any files >= 2GB...

  /* set mass and radius from PTE. */
  *optflags = MOLFILE_ATOMICNUMBER | MOLFILE_MASS | MOLFILE_RADIUS; 

  for(i=0;i<cube->numatoms;i++) {
    int idx;
    float chrg, coord;
    char fbuffer[1024];

    atom = atoms + i;

    k = fgets(fbuffer, 1024, cube->fd);
    j=sscanf(fbuffer, "%d %f %f %f %f", &idx, &chrg, &coord, &coord, &coord);
    if (k == NULL) {
      fprintf(stderr, "cube structure) missing atom(s) in "
              "file '%s'\n",cube->file_name);
      fprintf(stderr, "cube structure) expecting '%d' atoms,"
              " found only '%d'\n",cube->numatoms,i+1);
      return MOLFILE_ERROR;
    } else if (j < 5) {
      fprintf(stderr, "cube structure) missing type or coordinate(s)"
              " in file '%s' for atom '%d'\n",cube->file_name,i+1);
      return MOLFILE_ERROR;
    }

    // assign atom symbol to number. flag unknown or dummy atoms with X.
    atom->atomicnumber = idx;
    strncpy(atom->name, get_pte_label(idx), sizeof(atom->name));
    strncpy(atom->type, atom->name, sizeof(atom->type));
    atom->mass = get_pte_mass(idx);
    atom->radius = get_pte_vdw_radius(idx);
    atom->resname[0] = '\0';
    atom->resid = 1;
    atom->chain[0] = '\0';
    atom->segid[0] = '\0';
    /* skip to the end of line */
  }

  return MOLFILE_SUCCESS;
}


static int read_cube_timestep(void *v, int natoms, molfile_timestep_t *ts) {
  int i, j, n;
  char fbuffer[1024];
  float x, y, z;
  char *k;
  
  cube_t *cube = (cube_t *)v;

  // there is only one set of coordinates
  if (cube->coord) return MOLFILE_EOF;
  cube->coord = true;

  // jump to coordinate position
  fseek(cube->fd, cube->crdpos, SEEK_SET);
  // XXX fseek()/ftell() are incompatible with 64-bit LFS I/O implementations, 
  // hope we don't read any files >= 2GB...
 
  /* read the coordinates */
  for (i=0; i<cube->numatoms; i++) {
    k = fgets(fbuffer, 1024, cube->fd);
    j = sscanf(fbuffer, "%*d %*f %f %f %f", &x, &y, &z);
    
    if (k == NULL) {
      return MOLFILE_ERROR;
    } else if (j < 3) {
      fprintf(stderr, "cube timestep) missing type or coordinate(s) in file '%s' for atom '%d'\n",cube->file_name,i+1);
      return MOLFILE_ERROR;
    } else if (j>=3) {
      if (ts != NULL) { 
        // Only save coords if we're given a timestep pointer, 
        // otherwise assume that VMD wants us to skip past it.
        
        // In order to make the periodic display work, we need to
        // rotate the coordinates around the origin by the stored
        // rotation matrix. All coordinates are in Bohr, so they
        // must be converted to Angstrom, too.
        x -= cube->origin[0];
        y -= cube->origin[1];
        z -= cube->origin[2];
        
        for (n=0; n<3; ++n) {
          ts->coords[3*i + n] = bohr*(cube->origin[n] 
                                      + cube->rotmat[n][0] * x
                                      + cube->rotmat[n][1] * y
                                      + cube->rotmat[n][2] * z);
        }
      }
    } else {
      break;
    }
  }
  // set unitcell dimensions from cached data.
  if (ts != NULL) { 
      ts->A = cube->box.A;
      ts->B = cube->box.B;
      ts->C = cube->box.C;
      ts->alpha = cube->box.alpha;
      ts->beta  = cube->box.beta;
      ts->gamma = cube->box.gamma;
  }
  
  return MOLFILE_SUCCESS;
}

static int read_cube_metadata(void *v, int *nsets, 
  molfile_volumetric_t **metadata) {
  cube_t *cube = (cube_t *)v;
  *nsets = cube->nsets; 
  *metadata = cube->vol;  

  return MOLFILE_SUCCESS;
}

static int read_cube_data(void *v, int set, float *datablock, float *colorblock) {
  cube_t *cube = (cube_t *)v;

  fprintf(stderr, "cubeplugin) trying to read cube data set %d\n", set);

  int xsize = cube->vol[set].xsize; 
  int ysize = cube->vol[set].ysize;
  int zsize = cube->vol[set].zsize;
  int xysize = xsize*ysize;
  int nsize = cube->nsets;
  int nzsize = nsize*zsize;
  int nyzsize = nsize*zsize*ysize;
  int x, y, z;

  // go to data
  fseek(cube->fd, cube->datapos, SEEK_SET);
  // XXX fseek()/ftell() are incompatible with 64-bit LFS I/O implementations, 
  // hope we don't read any files >= 2GB...

  // read the data values in 
  if (cube->nsets == 1) { // density cube file
    for (x=0; x<xsize; x++) {
      for (y=0; y<ysize; y++) {
        for (z=0; z<zsize; z++) {
          if (fscanf(cube->fd, "%f", &datablock[z*xysize + y*xsize + x]) != 1) {
              return MOLFILE_ERROR;
          }
        }
      } 
    }
  } else {
    // XXX we may wish to examine this strategy for alternatives that provide
    // the same performance but without the extra copy, but it makes sense 
    // for now.  

    // Since the orbital cube file stores the data orb1(a1,b1,c1), orb2(a1,b1,c1), 
    // ... orbn(a1,b1,c1), orb1(a1,b1,c2), orb2(a1,a1,c2), ..., orbn(ai,bj,ck)
    // we have to cache the whole data set of have any kind of reasonable performance.
    // otherwise we would have to read (and parse!) the whole file over and over again.
    // this way we have to do it only once at the temporary expense of some memory.
    if (cube->datacache == NULL) {
      int points = xsize*ysize*zsize * nsize; // number of grid cells * nsets
      int i;

      // let people know what is going on.
      fprintf(stderr, "\ncubeplugin) creating %d MByte cube orbital cache.\n", 
              (int) (points*sizeof(float)) / 1048576);
      cube->datacache = new float[points];
            
      for (i=0; i < points; ++i) {
        if (fscanf(cube->fd, "%f", &cube->datacache[i]) != 1) {
          return MOLFILE_ERROR;
        }

        // print an ascii progress bar so impatient people do not get scared.
        if ((i % (1048576/sizeof(float))) == 0) {  // one dot per MB.
          fprintf(stderr, "."); 
        }
      }
    }
      
    for (x=0; x<xsize; x++) {
      for (y=0; y<ysize; y++) {
        for (z=0; z<zsize; z++) {
          datablock[z*xysize + y*xsize + x] = cube->datacache[x*nyzsize + y*nzsize + z*nsize + set];
        }
      }
    }
  }

  return MOLFILE_SUCCESS;
}

static void close_cube_read(void *v) {
  cube_t *cube = (cube_t *)v;

  fclose(cube->fd);
  if (cube->vol) {
    delete[] cube->vol; 
  }
  free(cube->file_name);
  if (cube->datacache) { 
    fprintf(stderr, "\ncubeplugin) freeing cube orbital cache.\n");
    delete[] cube->datacache; 
  }  

  delete cube;
}

/*
 * Initialization stuff here
 */
static molfile_plugin_t plugin = {
  vmdplugin_ABIVERSION,   // ABI version
  MOLFILE_PLUGIN_TYPE, 	  // plugin type
  "cube",                 // short file format description
  "Gaussian Cube",        // pretty file format description
  "Axel Kohlmeyer, John E. Stone", // author(s)
  0,                      // major version
  8,                      // minor version
  VMDPLUGIN_THREADSAFE,   // is reentrant
  "cube",                 // filename extension
  open_cube_read,               
  read_cube_structure,
  0,                      // read_bonds
  read_cube_timestep,
  close_cube_read,
  0,                      // open_file_write
  0,                      // write_structure
  0,                      // write_timestep
  0,                      // close_file_write
  read_cube_metadata,
  read_cube_data,
  0                       // read_rawgraphics
};

int VMDPLUGIN_init(void) { return VMDPLUGIN_SUCCESS; }
int VMDPLUGIN_fini(void) { return VMDPLUGIN_SUCCESS; }
int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)&plugin);
  return VMDPLUGIN_SUCCESS;
}



#ifdef TEST_PLUGIN

int main(int argc, char *argv[]) {
  int natoms;
  void *v;
  int i, nsets, set;
  molfile_volumetric_t * meta;

  while (--argc) {
    ++argv;
    v = open_cube_read(*argv, "cube", &natoms);
    if (!v) {
      fprintf(stderr, "open_cube_read failed for file %s\n", *argv);
      return 1;
    }
    fprintf(stderr, "open_cube_read succeeded for file %s\n", *argv);

    // try loading the EDM metadata now
    if (read_cube_metadata(v, &nsets, &meta)) {
      return 1; // failed to load cube file
    }
    fprintf(stderr, "read_cube_metadata succeeded for file %s\n", *argv);

    for (set=0; set<nsets; set++) {
      printf("Loading volume set: %d\n", set);   
      
      int elements = meta[set].xsize * meta[set].ysize * meta[set].zsize;
      printf("   Grid Elements: %d\n", elements);
      printf(" Grid dimensions: X: %d Y: %d Z: %d\n", 
             meta[set].xsize, meta[set].ysize, meta[set].zsize);

      float * voldata = (float *) malloc(sizeof(float) * elements);
      float * coldata = NULL;

      if (meta[set].has_color) {
        coldata = (float *) malloc(sizeof(float) * elements * 3);
      }

      // try loading the data sets now
      if (read_cube_data(v, set, voldata, coldata)) {
        return 1; // failed to load cube file
      }

      printf("First 6 elements:\n   ");
      for (i=0; i<6; i++) {
        printf("%g, ", voldata[i]);
      }
      printf("\n"); 

      printf("Last 6 elements:\n   ");
      for (i=elements - 6; i<elements; i++) {
        printf("%g, ", voldata[i]);
      }
      printf("\n"); 
    }

    close_cube_read(v);
  }
  return 0;
}

#endif
