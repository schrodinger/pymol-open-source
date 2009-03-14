/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_xsfplugin
#define STATIC_PLUGIN 1

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: xsfplugin.C,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.17 $       $Date: 2008/01/09 20:37:47 $
 *
 ***************************************************************************/

//
// Molefile plugin for xsf/axsf format files as created by the 
// Quantum Espresso software package (http://www.quantum-espresso.org) 
// a.k.a. PWScf (http://www.pwscf.org/) and XCrySDen (http://www.xcrysden.org/
//
// a file format description is available at:
// http://www.xcrysden.org/doc/XSF.html
//

#include "largefiles.h"   /* platform dependent 64-bit file I/O defines */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "molfile_plugin.h"

#ifndef M_PI_2
#define M_PI_2 1.57079632679489661922
#endif

#include "periodic_table.h"

static const char *xsf_symtab[] = {
  "(unknown keyword)", "#",
  "BEGIN_INFO", "END_INFO",
  "BEGIN_BLOCK_DATAGRID_2D", "END_BLOCK_DATAGRID_2D",
  "BEGIN_DATAGRID_2D", "END_DATAGRID_2D",
  "BEGIN_BLOCK_DATAGRID_3D", "END_BLOCK_DATAGRID_3D",
  "BEGIN_DATAGRID_3D", "END_DATAGRID_3D",
  "BEGIN_BLOCK_BANDGRID_3D", "END_BLOCK_BANDGRID_3D",
  "ATOMS", "ANIMSTEPS", "BAND",
  "MOLECULE", "POLYMER", "SLAB", "CRYSTAL",
  "PRIMVEC", "CONVVEC", "PRIMCOORD", "CONVCOORD"
};

typedef enum {
  xsf_UNKNOWN = 0,  xsf_COMMENT,
  xsf_BEGINFO,      xsf_ENDINFO, 
  xsf_BEG_2D_BLOCK, xsf_END_2D_BLOCK,
  xsf_BEG_2D_DATA,  xsf_END_2D_DATA,
  xsf_BEG_3D_BLOCK, xsf_END_3D_BLOCK,
  xsf_BEG_3D_DATA,  xsf_END_3D_DATA,
  xsf_BEG_3D_BAND,  xsf_END_3D_BAND,
  xsf_ATOMS, xsf_ANIMSTEPS, xsf_BAND,
  xsf_MOLECULE, xsf_POLYMER, xsf_SLAB, xsf_CRYSTAL,
  xsf_PRIMVEC, xsf_CONVVEC, xsf_PRIMCOORD, xsf_CONVCOORD,
  xsf_NR_KEYWORDS
} xsf_keyword_t;

// list of known alternatives to the keywords above
// last entry has to be an xsf_UNKNOWN.
static const struct {
  const char *name;
  xsf_keyword_t kw;
}  xsf_aliases[] = {
  { "DATAGRID_2D", xsf_BEG_2D_DATA },
  { "DATAGRID_3D", xsf_BEG_3D_DATA },
  { "BEGIN_BLOCK_DATAGRID2D", xsf_BEG_2D_BLOCK },
  { "BEGIN_BLOCK_DATAGRID3D", xsf_BEG_3D_BLOCK },
  { "END_BLOCK_DATAGRID2D", xsf_END_2D_BLOCK },
  { "END_BLOCK_DATAGRID3D", xsf_END_3D_BLOCK },
  { NULL,          xsf_UNKNOWN     }
};

static xsf_keyword_t lookup_keyword(const char* word)
{
  int i, j;
  
  if (word == 0) return xsf_UNKNOWN;
  
  // find start of word.
  j=0;
  for (i=0; i < (int)strlen(word); ++i) {
    j=i;
    if (! isspace(word[i])) break;
  }
  
  for (i=1; i < xsf_NR_KEYWORDS; ++i) {
    if (0 == strncmp(word + j, xsf_symtab[i], strlen(xsf_symtab[i])))
      return (xsf_keyword_t) i;
  }

  // check for known aliases/alternatives
  for (i=0; xsf_aliases[i].kw != xsf_UNKNOWN; ++i) {
    const char *name = xsf_aliases[i].name;

    if (0 == strncmp(word + j, name, strlen(name)))
      return xsf_aliases[i].kw;

  }

  return xsf_UNKNOWN;
}

// A format-independent structure to hold unit cell data
typedef struct {
  float A, B, C, alpha, beta, gamma, cell[3][3];
} xsf_box;

typedef struct {
  FILE *fd;                     // file descriptor
  int nvolsets;                 // number of volumetric datasets
  int numatoms;                 // number of atoms
  int animsteps;                // for comparison.
  int numsteps;                 // number of coordinate sets.
  bool coord;                   // has coordinate data
  char *file_name;              // original filename 
  xsf_keyword_t pbctype;        // type of periodicity (none/polymer/slab/crystal)
  molfile_volumetric_t *vol;    // volume set metadata 
  int numvolmeta;               // number of entries in *vol
  float origin[3];              // origin, stored for periodic display hack 
  float rotmat[3][3];           // rotation matrix, stored for periodic display hack
  float invmat[3][3];           // reciprocal cell matrix (for PBC wrapping).
  xsf_box box;                  // unit cell dimensions (for VMD).
} xsf_t;


// Converts box basis vectors to A, B, C, alpha, beta, and gamma.  
// Stores values in xsf_box struct, which should be allocated before calling
// this function.
static int xsf_readbox(xsf_box *box, float *x, float *y, float *z) {
  float A, B, C;
  int i;

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

  // copy original cell vectors for PBC wrapping.
  for (i=0; i<3; ++i) {
    box->cell[0][i] = x[i];
    box->cell[1][i] = y[i];
    box->cell[2][i] = z[i];
  }
  
  return 0;
}

// calculate and store rotation matrix to realign everything later.
static void xsf_buildrotmat(xsf_t *xsf, float *a, float *b)
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
    for (int j=0; j<3; ++j) {
      xsf->rotmat[i][j] = r[i][j];
    }
  }
}

static void xsf_buildinvmat(xsf_t *xsf, float *a, float *b, float *c)
{
  float det, id;
  
  det = a[0]*b[1]*c[2] + b[0]*c[1]*a[2] + c[0]*a[1]*b[2] 
    - a[0]*c[1]*b[2] - b[0]*a[1]*c[2] - c[0]*b[1]*a[2];
  
  id = 1.0 / det;
  xsf->invmat[0][0] = id * ( b[1]*c[2]-b[2]*c[1] );
  xsf->invmat[1][0] = id * ( a[2]*c[1]-a[1]*c[2] );
  xsf->invmat[2][0] = id * ( a[1]*b[2]-a[2]*b[1] );
  xsf->invmat[0][1] = id * ( b[2]*c[0]-b[0]*c[2] );
  xsf->invmat[1][1] = id * ( a[0]*c[2]-a[2]*c[0] );
  xsf->invmat[2][1] = id * ( a[2]*b[0]-a[0]*b[2] );
  xsf->invmat[0][2] = id * ( b[0]*c[1]-b[1]*c[0] );
  xsf->invmat[1][2] = id * ( a[1]*c[0]-a[0]*c[1] );
  xsf->invmat[2][2] = id * ( a[0]*b[1]-a[1]*b[0] );
}


// read a line and forget the data
static void eatline(FILE *fd) {
  char readbuf[1025];
  fgets(readbuf, 1024, fd);    // go on to next line
}  

static bool xsf_read_cell(FILE *fd, float *a, float *b, float *c)
{
  return (9 == fscanf(fd, "%f%f%f%f%f%f%f%f%f", 
                      &a[0],&a[1],&a[2],
                      &b[0],&b[1],&b[2],
                      &c[0],&c[1],&c[2]));
}

static void close_xsf_read(void *v);

static void *open_xsf_read(const char *filepath, const char *filetype,
                           int *natoms) {
  FILE *fd;
  xsf_t *xsf;
  int i,j;
  
  fd = fopen(filepath, "rb");
  if (!fd) 
    return NULL;

  xsf = new xsf_t;
  xsf->fd = fd;
  xsf->vol = NULL;
  xsf->numvolmeta = 0;
  xsf->coord = false;
  xsf->nvolsets = 0;
  xsf->numatoms = 0;
  xsf->numsteps = 0;
  xsf->file_name = strdup(filepath);
  // this will disable pbc wrapping of coordinates by default
  xsf->pbctype = xsf_MOLECULE; 

  // initialize origin and rotmat to sensible defaults.
  for (i=0; i<3; ++i) {
    for (j=0; j<3; ++j) {
      xsf->rotmat[i][j] = 0.0;
    }
  }
  for (i=0; i<3; ++i) {
    xsf->origin[i] = 0.0;
    xsf->rotmat[i][i] = 1.0;
  }

  // since there can be all kinds of data in the file, 
  // we start by scanning through the whole file and analyse
  // the available info.
  char readbuf[256]; // line buffer
  xsf_keyword_t kw;

  // we loop until can't read anymore.
  do {
    if (NULL == fgets(readbuf, 256, xsf->fd)) break;

    again:
    kw = lookup_keyword(readbuf);
#ifdef TEST_PLUGIN
    fprintf(stderr, "keyword: %d / %s", kw, readbuf);
#endif          
    
    switch (kw) {
      case xsf_ANIMSTEPS:
#ifdef TEST_PLUGIN
      {
        int n;
        if (1 == sscanf(readbuf, "%*s%d", &n)) {
          fprintf(stderr, "ANIMSTEPS: found %d steps\n", n);
        }
      }
#endif          
      break;
        
      case xsf_ATOMS: // no specification for the number of atoms, so we
                      // try to figure them out
        ++ xsf->numsteps;
        if (xsf->numatoms == 0) { // count atoms only, if we don't know how many
          while (fgets(readbuf, 256, xsf->fd)) {
            float x,y,z;
            // the coordinate lines are <index> <x> <y> <z> with optional forces.
            if (3 == sscanf(readbuf, "%*s%f%f%f", &x, &y, &z)) {
              ++ xsf->numatoms;
            } else {
              // we've most likely read the next keyword.
              // reparse buffer.
              goto again;
              break;
            }
          }
#ifdef TEST_PLUGIN
          fprintf(stderr, "ATOMS: found %d atoms\n", xsf->numatoms);
#endif          
        }  else { // skip over the lines
          int n;
          for (n=0; n < xsf->numatoms; ++n) eatline(xsf->fd);
        }
        break;
        
      case xsf_PRIMCOORD: // number of atoms is in the next line
         
        if(fgets(readbuf, 256, xsf->fd)) {
          int mol, mult;
          
          if (xsf->numatoms == 0) {
            if (2 == sscanf(readbuf, "%d%d", &mol, &mult)) {
              xsf->numatoms = mol * mult;
            } else {
              xsf->numatoms = mol;
            }
          }
          // skip over atom coordinates
          int n;
          for (n=0; n < xsf->numatoms; ++n) eatline(xsf->fd);
          ++ xsf->numsteps; 

#ifdef TEST_PLUGIN
          fprintf(stderr, "PRIMCOORD: found %d atoms\n", xsf->numatoms);
#endif          
        }
        break;
        
      case xsf_CONVCOORD: // number of atoms is in the next line
         
        if(fgets(readbuf, 256, xsf->fd)) {
          int mol, mult, num;
          
          num = 0;
          if (2 == sscanf(readbuf, "%d%d", &mol, &mult)) {
            num = mol * mult;
          }
          
          // skip over atom coordinates
          int n;
          for (n=0; n < num; ++n) eatline(xsf->fd);
#ifdef TEST_PLUGIN
          fprintf(stderr, "CONVCOORD: ignoring %d atoms\n", num);
#endif          
        }
        break;

      case xsf_PRIMVEC: // store primitive cell info for rotation of the volumetric data grid vectors
      {
        float a[3], b[3], c[3];
        
        if (xsf_read_cell(xsf->fd, a, b, c)) {
          xsf_buildrotmat(xsf, a, b);
        } else {
          fprintf(stderr, "xsfplugin) WARNING: error reading unit cell. ignoring unit cell info.\n");
        }
      }
      break;
      
      case xsf_CONVVEC: // ignore conventional cells.
      {
        int n;
        for (n=0; n < 3; ++n) eatline(xsf->fd);
      }
      break;

      case xsf_BEG_3D_BLOCK: // analyse number of 3d-data sets
        // ordinarily the parsing of the metadata would be done in read_xsf_metadata()
        // but since we have to move through the whole file to count its contents,
        // it is faster to parse it here and just pass the data later.

        if (xsf->vol == NULL) { // initialize the volume set list with 32 entries
          xsf->numvolmeta = 32;
          xsf->vol = new molfile_volumetric_t[xsf->numvolmeta];
        }
        
        // next line is title, then check for data blocks
        fgets(readbuf, 256, xsf->fd);
        printf("xsfplugin) found grid data block: %s", readbuf);

        do { // loop until we reach the end of the whole block 
             // or run out of data.
          if (NULL == fgets(readbuf, 256, xsf->fd)) break;
          switch (lookup_keyword(readbuf)) {
            case xsf_BEG_3D_DATA: // fallthrough
            {
              int n;
              molfile_volumetric_t *set;
              float a[3], b[3], c[3];
              
              ++ xsf->nvolsets;

              // double the size of the cache for metainfo, if needed
              if (xsf->nvolsets > xsf->numvolmeta) {
                molfile_volumetric_t *ptr = xsf->vol;
                xsf->vol = new molfile_volumetric_t[2 * xsf->numvolmeta];
                memcpy((void *)xsf->vol, (void *)ptr, xsf->numvolmeta*sizeof(molfile_volumetric_t));
                xsf->numvolmeta *= 2;
                delete[] ptr;
              }

              // get a handle to the current volume set meta data
              set = &(xsf->vol[xsf->nvolsets - 1]);
              set->has_color = 0;

              // the begin mark is also the title of the data set.
              // we need the exact name to later find the start of the data set.
              strncpy(set->dataname, readbuf, 255);
              
              // next is the number of grid points, the origin and
              // the spanning vectors of the data block
              fgets(readbuf, 256, xsf->fd);
              sscanf(readbuf, "%d%d%d", &(set->xsize), &(set->ysize), &(set->zsize));
              fgets(readbuf, 256, xsf->fd);
              sscanf(readbuf, "%f%f%f", &(set->origin[0]), &(set->origin[1]), &(set->origin[2]));
              fgets(readbuf, 256, xsf->fd);
              sscanf(readbuf, "%f%f%f", &a[0], &a[1], &a[2]);
              fgets(readbuf, 256, xsf->fd);
              sscanf(readbuf, "%f%f%f", &b[0], &b[1], &b[2]);
              fgets(readbuf, 256, xsf->fd);
              sscanf(readbuf, "%f%f%f", &c[0], &c[1], &c[2]);
              
              // we need to fix up the size of the data points, since xsf file 
              // store the data points at the borders on both sides.
              -- set->xsize; -- set->ysize; -- set->zsize;

              // store the realigned axes.
              for (n=0; n<3; ++n) {
                set->xaxis[n] = xsf->rotmat[n][0] * a[0] 
                  + xsf->rotmat[n][1] * a[1] + xsf->rotmat[n][2] * a[2];
                
                set->yaxis[n] = xsf->rotmat[n][0] * b[0] 
                  + xsf->rotmat[n][1] * b[1] + xsf->rotmat[n][2] * b[2];
    
                set->zaxis[n] = xsf->rotmat[n][0] * c[0] 
                  + xsf->rotmat[n][1] * c[1] + xsf->rotmat[n][2] * c[2];
              }

              do { // loop until we reach the end of the data set
                fgets(readbuf, 256, xsf->fd);
              } while (xsf_END_3D_DATA != lookup_keyword(readbuf));

#if 1
              /*   as of VMD version 1.8.3, volumetric data points are 
	       *   expected to represent the center of a grid box. xsf format 
	       *   volumetric data represents the value at the edges of the 
	       *   grid boxes, so we need to shift the internal origin by half 
	       *   a grid box diagonal to have the data at the correct position 
               *   This will need to be changed again when the plugin interface
	       *   is updated to explicitly allow point/face-centered data sets.
	       */
              set->origin[0] -= 0.5 * ( set->xaxis[0] / (double) set->xsize
                                        + set->yaxis[0] / (double) set->ysize
                                        + set->zaxis[0] / (double) set->zsize);
              set->origin[1] -= 0.5 * ( set->xaxis[1] / (double) set->xsize
                                        + set->yaxis[1] / (double) set->ysize
                                        + set->zaxis[1] / (double) set->zsize);
              set->origin[2] -= 0.5 * ( set->xaxis[2] / (double) set->xsize
                                        + set->yaxis[2] / (double) set->ysize
                                        + set->zaxis[2] / (double) set->zsize);
#endif
            }
            break;

            default:
              break;
          }
        } while (xsf_END_3D_BLOCK != lookup_keyword(readbuf));
        
#ifdef TEST_PLUGIN
        fprintf(stderr, "found %d volumetric data sets\n", xsf->nvolsets);
#endif          
        break;

      case xsf_BEG_2D_BLOCK:
        do { // skip over data
          fgets(readbuf, 256, xsf->fd);
        } while (xsf_END_2D_BLOCK != lookup_keyword(readbuf));
        break;
        
        // periodicity encoding. needed for coordinate wrapping.
      case xsf_MOLECULE: // fallthrough
      case xsf_SLAB:     // fallthrough
      case xsf_POLYMER:  // fallthrough
      case xsf_CRYSTAL:
        xsf->pbctype = kw;
        break;

      case xsf_COMMENT:  // fallthrough
      case xsf_UNKNOWN:  // fallthrough
      default:                  // ignore everything unknown
        break;
        
    }
  } while (! (feof(xsf->fd) || ferror(xsf->fd)));
#ifdef TEST_PLUGIN
  fprintf(stderr, "total of %d coordinate sets\n", xsf->numsteps);
#endif
  rewind(xsf->fd);
  *natoms = xsf->numatoms;
  return xsf;
}

  
static int read_xsf_structure(void *v, int *optflags, molfile_atom_t *atoms) {
  int i;
  xsf_t *xsf = (xsf_t *)v;

  // return immediately if there is no structure in this file.
  if (xsf->numatoms < 1) return MOLFILE_SUCCESS;

  
  // go to beginning of file and find first set of coordinates
  rewind(xsf->fd);

  // we loop until we have found, read and parsed the first 
  // set of atomic coordinates. we only accept ATOMS and PRIMCOORD
  // sections. if we happen to find a PRIMVEC section, too, we use
  // that to set the default for animations.
  char readbuf[1024]; // line buffer
  do {
    if (NULL == fgets(readbuf, 256, xsf->fd)) break;

    switch (lookup_keyword(readbuf)) {

      case xsf_PRIMCOORD: // number of atoms is in the next line. skip.
        eatline(xsf->fd);
        // fallthrough
      case xsf_ATOMS: // no specification for the number of atoms, 
                      // we can use the same parser for both sections.

        /* we set atom mass and VDW radius from the PTE. */
        *optflags = MOLFILE_ATOMICNUMBER | MOLFILE_MASS | MOLFILE_RADIUS;

        for(i=0; i<xsf->numatoms; ++i) {
          int j;
          char *k;
          float coord;
          molfile_atom_t *atom;

          char buffer[1024];
          k = fgets(readbuf, 1024, xsf->fd);
          atom = atoms + i;
          j=sscanf(readbuf, "%s %f %f %f", buffer, &coord, &coord, &coord);
          if (k == NULL) {
            fprintf(stderr, "xsfplugin) structure missing atom(s) in file '%s'\n",xsf->file_name);
            fprintf(stderr, "xsfplugin) expecting '%d' atoms, found only '%d'\n",xsf->numatoms,i+1);
            return MOLFILE_ERROR;
          } else if (j < 4) {
            fprintf(stderr, "xsfplugin) missing type or coordinate(s) in file '%s' for atom '%d'\n",
                    xsf->file_name, i+1);
            return MOLFILE_ERROR;
          }

          /* handle the case if the first item is an ordinal number 
           * from the PTE */
          if (isdigit(buffer[0])) {
            int idx;
            idx = atoi(buffer);
            strncpy(atom->name, get_pte_label(idx), sizeof(atom->name));
            atom->atomicnumber = idx;
            atom->mass = get_pte_mass(idx);
            atom->radius = get_pte_vdw_radius(idx);
          } else {
            int idx;
            strncpy(atom->name, buffer, sizeof(atom->name));
            idx = get_pte_idx(buffer);
            atom->atomicnumber = idx;
            atom->mass = get_pte_mass(idx);
            atom->radius = get_pte_vdw_radius(idx);
          }
          strncpy(atom->type, atom->name, sizeof(atom->type));
          atom->resname[0] = '\0';
          atom->resid = 1;
          atom->chain[0] = '\0';
          atom->segid[0] = '\0';
#ifdef TEST_PLUGIN
          fprintf(stderr,"xsfplugin) atom %4d: %s  / mass= %f\n", i, atom->name, atom->mass);
#endif
        }
        
        // ok. done. rewind once more and get the hell out of here.
        rewind(xsf->fd);
        return MOLFILE_SUCCESS;
        break;
        
        // read primitive cell info.
      case xsf_PRIMVEC: 
      {
        float a[3], b[3], c[3];

        if (xsf_read_cell(xsf->fd, a, b, c)) { // ignore unit cell info,if we cannot parse it.
          xsf_readbox(&(xsf->box), a, b, c);
          xsf_buildrotmat(xsf, a, b);
          // print warning, if the rotation will be significant:
          if ((fabs((double) a[1]) + fabs((double) a[2]) + fabs((double) b[2]))
              > 0.001) {
            fprintf(stderr, "xsfplugin) WARNING: Coordinates will be rotated to comply \n"
                    "xsfplugin) with VMD's conventions for periodic display...\n");
          }
          xsf_buildinvmat(xsf, a, b, c);

#if defined(TEST_PLUGIN)
        printf("cell vectors:\n");
        printf("<a>: %12.8f %12.8f %12.8f\n", a[0], a[1], a[2]);
        printf("<b>: %12.8f %12.8f %12.8f\n", b[0], b[1], b[2]);
        printf("<c>: %12.8f %12.8f %12.8f\n", c[0], c[1], c[2]);
        printf("cell dimensions:\n");
        printf("a= %12.8f   b= %12.8f   c= %12.8f\n", xsf->box.A, xsf->box.B, xsf->box.C);
        printf("alpha= %6.2f  beta= %6.2f  gamma= %6.2f\n", xsf->box.alpha, xsf->box.beta, xsf->box.gamma);
        printf("reciprocal cell vectors:\n");
        printf("i: %12.8f %12.8f %12.8f\n", xsf->invmat[0][0], xsf->invmat[0][1], xsf->invmat[0][2]);
        printf("k: %12.8f %12.8f %12.8f\n", xsf->invmat[1][0], xsf->invmat[1][1], xsf->invmat[1][2]);
        printf("l: %12.8f %12.8f %12.8f\n", xsf->invmat[2][0], xsf->invmat[2][1], xsf->invmat[2][2]);
        printf("cell rotation matrix:\n");
        printf("x: %12.8f %12.8f %12.8f\n", xsf->rotmat[0][0], xsf->rotmat[0][1], xsf->rotmat[0][2]);
        printf("y: %12.8f %12.8f %12.8f\n", xsf->rotmat[1][0], xsf->rotmat[1][1], xsf->rotmat[1][2]);
        printf("z: %12.8f %12.8f %12.8f\n", xsf->rotmat[2][0], xsf->rotmat[2][1], xsf->rotmat[2][2]);
#endif
        }
      }
      break;

      default:                  // ignore everything unknown
        break;
    }
  } while (! (feof(xsf->fd) || ferror(xsf->fd)));

  // if we reach this point, some error must have happened.
  return MOLFILE_ERROR;

}


static int read_xsf_timestep(void *v, int natoms, molfile_timestep_t *ts) {
  int i;
  
  xsf_t *xsf = (xsf_t *)v;

  // we loop until we have found, read and parsed the next.
  // set of atomic coordinates. we only accept ATOMS and PRIMCOORD
  // sections. if we happen to find a PRIMVEC section, too, we use
  // that to reset the cell parameters.

  char readbuf[1024]; // line buffer
  do {
    if (NULL == fgets(readbuf, 256, xsf->fd)) break;

    switch (lookup_keyword(readbuf)) {

      case xsf_PRIMCOORD: // number of atoms is in the next line. skip.
        eatline(xsf->fd);
        // fallthrough
      case xsf_ATOMS: // no specification for the number of atoms, 
                      // we can use the same parser for both sections.

        for(i=0; i<natoms; ++i) {
          int j, n;
          char *k, buffer[1024];
          float x, y, z;
          
          k = fgets(readbuf, 1024, xsf->fd);
          j=sscanf(readbuf, "%s %f %f %f", buffer, &x, &y, &z);

          if (k == NULL) {
            return MOLFILE_ERROR;
          } else if (j < 4) {
            fprintf(stderr, "xsfplugin) missing type or coordinate(s) in file '%s' for atom '%d'\n",
                    xsf->file_name, i+1);
            return MOLFILE_ERROR;
          } else if (j>=3) {
            if (ts != NULL) { 
              // Only save coords if we're given a timestep pointer, 
              // otherwise assume that VMD wants us to skip past it.
              float xf, yf, zf;
              
              // apply periodic boundary conditions.
#ifdef TEST_PLUGIN
                  printf("wrap PBC: before: %12.6f %12.6f %12.6f\n", x, y, z);
#endif                  
              switch(xsf->pbctype) {
                case xsf_CRYSTAL:
                  xf = xsf->invmat[0][0] * x + xsf->invmat[0][1] * y + xsf->invmat[0][2] * z;
                  xf = xf - floor((double)xf);
                  yf = xsf->invmat[1][0] * x + xsf->invmat[1][1] * y + xsf->invmat[1][2] * z;
                  yf = yf - floor((double)yf);
                  zf = xsf->invmat[2][0] * x + xsf->invmat[2][1] * y + xsf->invmat[2][2] * z;
                  zf = zf - floor((double)zf);
                  x = xsf->box.cell[0][0] * xf + xsf->box.cell[0][1] * yf + xsf->box.cell[0][2] * zf;
                  y = xsf->box.cell[1][0] * xf + xsf->box.cell[1][1] * yf + xsf->box.cell[1][2] * zf;
                  z = xsf->box.cell[2][0] * xf + xsf->box.cell[2][1] * yf + xsf->box.cell[2][2] * zf;
                  break;

                case xsf_SLAB:
                  xf = xsf->invmat[0][0] * x + xsf->invmat[0][1] * y + xsf->invmat[0][2] * z;
                  xf = xf - floor((double)xf);
                  yf = xsf->invmat[1][0] * x + xsf->invmat[1][1] * y + xsf->invmat[1][2] * z;
                  yf = yf - floor((double)yf);
                  zf = xsf->invmat[2][0] * x + xsf->invmat[2][1] * y + xsf->invmat[2][2] * z;
                  x = xsf->box.cell[0][0] * xf + xsf->box.cell[0][1] * yf + xsf->box.cell[0][2] * zf;
                  y = xsf->box.cell[1][0] * xf + xsf->box.cell[1][1] * yf + xsf->box.cell[1][2] * zf;
                  z = xsf->box.cell[2][0] * xf + xsf->box.cell[2][1] * yf + xsf->box.cell[2][2] * zf;
                  break;
                  
                case xsf_POLYMER:
                  xf = xsf->invmat[0][0] * x + xsf->invmat[0][1] * y + xsf->invmat[0][2] * z;
                  xf = xf - floor((double)xf);
                  yf = xsf->invmat[1][0] * x + xsf->invmat[1][1] * y + xsf->invmat[1][2] * z;
                  zf = xsf->invmat[2][0] * x + xsf->invmat[2][1] * y + xsf->invmat[2][2] * z;
                  x = xsf->box.cell[0][0] * xf + xsf->box.cell[0][1] * yf + xsf->box.cell[0][2] * zf;
                  y = xsf->box.cell[1][0] * xf + xsf->box.cell[1][1] * yf + xsf->box.cell[1][2] * zf;
                  z = xsf->box.cell[2][0] * xf + xsf->box.cell[2][1] * yf + xsf->box.cell[2][2] * zf;
                  break;
                  
                case xsf_MOLECULE:
                  xf = x;
                  yf = y;
                  zf = z;
                  break;
                  
                default:
                  break;
              }
              
#ifdef TEST_PLUGIN
                  printf("wrap PBC: fract:  %12.6f %12.6f %12.6f\n", xf, yf, zf);
                  printf("wrap PBC: after: %12.6f %12.6f %12.6f\n", x, y, z);
#endif                  

              // In order to make the periodic display work, we need to
              // rotate the coordinates around the origin by the stored
              // rotation matrix. for xsf files the origin is not explicitely
              // so far, but since it is already initialized to (0,0,0) it
              // does no harm, to leave it in here.
              x -= xsf->origin[0];
              y -= xsf->origin[1];
              z -= xsf->origin[2];
        
              for (n=0; n<3; ++n) {
                ts->coords[3*i + n] = (xsf->origin[n] + xsf->rotmat[n][0] * x
                                       + xsf->rotmat[n][1] * y + xsf->rotmat[n][2] * z);
              }
            }
          } else {
            break;
          }
        }
        // (re-)set unitcell dimensions
        if (ts != NULL) { 
          ts->A = xsf->box.A;
          ts->B = xsf->box.B;
          ts->C = xsf->box.C;
          ts->alpha = xsf->box.alpha;
          ts->beta  = xsf->box.beta;
          ts->gamma = xsf->box.gamma;
        }
  
        return MOLFILE_SUCCESS;
        break;
        
        // read and update primitive cell info.
      case xsf_PRIMVEC: 
      {
        float a[3], b[3], c[3];
        
        if (xsf_read_cell(xsf->fd, a, b, c)) {
          xsf_readbox(&(xsf->box), a,b,c);
          xsf_buildrotmat(xsf, a, b);
          // print warning, if the rotation will be significant:
          if ((fabs((double) a[1]) + fabs((double) a[2]) + fabs((double) b[2]))
              > 0.001) {
            fprintf(stderr, "xsfplugin) WARNING: Coordinates will be rotated to comply \n"
                    "xsfplugin) with VMD's conventions for periodic display...\n");
          }
          xsf_buildinvmat(xsf, a, b, c);

#if defined(TEST_PLUGIN)
        printf("new cell vectors:\n");
        printf("<a>: %12.8f %12.8f %12.8f\n", a[0], a[1], a[2]);
        printf("<b>: %12.8f %12.8f %12.8f\n", b[0], b[1], b[2]);
        printf("<c>: %12.8f %12.8f %12.8f\n", c[0], c[1], c[2]);
        printf("new cell dimensions:\n");
        printf("a= %12.8f   b= %12.8f   c= %12.8f\n", xsf->box.A, xsf->box.B, xsf->box.C);
        printf("alpha= %6.2f  beta= %6.2f  gamma= %6.2f\n", xsf->box.alpha, xsf->box.beta, xsf->box.gamma);
        printf("new reciprocal cell vectors:\n");
        printf("i: %12.8f %12.8f %12.8f\n", xsf->invmat[0][0], xsf->invmat[0][1], xsf->invmat[0][2]);
        printf("k: %12.8f %12.8f %12.8f\n", xsf->invmat[1][0], xsf->invmat[1][1], xsf->invmat[1][2]);
        printf("l: %12.8f %12.8f %12.8f\n", xsf->invmat[2][0], xsf->invmat[2][1], xsf->invmat[2][2]);
        printf("new cell rotation matrix:\n");
        printf("x: %12.8f %12.8f %12.8f\n", xsf->rotmat[0][0], xsf->rotmat[0][1], xsf->rotmat[0][2]);
        printf("y: %12.8f %12.8f %12.8f\n", xsf->rotmat[1][0], xsf->rotmat[1][1], xsf->rotmat[1][2]);
        printf("z: %12.8f %12.8f %12.8f\n", xsf->rotmat[2][0], xsf->rotmat[2][1], xsf->rotmat[2][2]);
#endif
        }
      }
      break;

      default:                  // ignore everything unknown
        break;
    }
  } while (! (feof(xsf->fd) || ferror(xsf->fd)));

  // if we reach this point, some error must have happened.
  return MOLFILE_ERROR;
}

static int read_xsf_metadata(void *v, int *nvolsets, 
  molfile_volumetric_t **metadata) {
  xsf_t *xsf = (xsf_t *)v;
  *nvolsets = xsf->nvolsets; 
  *metadata = xsf->vol;  

  return MOLFILE_SUCCESS;
}

static int read_xsf_data(void *v, int set, float *datablock, float *colorblock) {
  xsf_t *xsf = (xsf_t *)v;
  const char *block = xsf->vol[set].dataname;
  
  fprintf(stderr, "xsfplugin) trying to read xsf data set %d: %s\n", set, block);
  
  int xsize = xsf->vol[set].xsize; 
  int ysize = xsf->vol[set].ysize;
  int zsize = xsf->vol[set].zsize;
  int x, y, z;
  int n;
  char readbuf[1024];
  float dummy;
  
  // find data set ...
  rewind(xsf->fd);
  do {
    if (NULL == fgets(readbuf, 1024, xsf->fd)) return MOLFILE_ERROR;
  } while (strncmp(readbuf, block, 1024));
  // ... and skip five more lines to get to the beginning of the data
  eatline(xsf->fd);
  eatline(xsf->fd);
  eatline(xsf->fd);
  eatline(xsf->fd);
  eatline(xsf->fd);
  
  // read in the data values
  // XSF data grids include the points at the border on both sides,
  // so we have to read and skip over those on one of them.
  n = 0;
  for (z=0; z<zsize+1; z++) {
    for (y=0; y<ysize+1; y++) {
      for (x=0; x<xsize+1; x++) {
        if((x>=xsize) || (y>=ysize) || (z>=zsize)) { 
          if (fscanf(xsf->fd, "%f", &dummy) != 1) return MOLFILE_ERROR;
        } else {
          if (fscanf(xsf->fd, "%f", &datablock[n]) != 1) return MOLFILE_ERROR;
          ++ n;
        }
      }
    }
  }
  rewind(xsf->fd);
  return MOLFILE_SUCCESS;
}

static void close_xsf_read(void *v) {
  xsf_t *xsf = (xsf_t *)v;

  fclose(xsf->fd);
  if (xsf->vol) {
    delete[] xsf->vol; 
  }
  free(xsf->file_name);
  delete xsf;
}

/*
 * Initialization stuff here
 */
static molfile_plugin_t plugin;

VMDPLUGIN_EXTERN int VMDPLUGIN_init(void) {
  memset(&plugin, 0, sizeof(molfile_plugin_t));
  plugin.abiversion = vmdplugin_ABIVERSION;
  plugin.type = MOLFILE_PLUGIN_TYPE;
  plugin.name = "xsf";
  plugin.prettyname = "(Animated) XCrySDen Structure File";
  plugin.author = "Axel Kohlmeyer, John Stone";
  plugin.majorv = 0;
  plugin.minorv = 7;
  plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
  plugin.filename_extension = "axsf,xsf";
  plugin.open_file_read = open_xsf_read;
  plugin.read_structure = read_xsf_structure;
  plugin.read_next_timestep =read_xsf_timestep;
  plugin.close_file_read = close_xsf_read;
  plugin.read_volumetric_metadata = read_xsf_metadata;
  plugin.read_volumetric_data = read_xsf_data;
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_EXTERN int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_EXTERN int VMDPLUGIN_fini(void) { return VMDPLUGIN_SUCCESS; }


#ifdef TEST_PLUGIN

int main(int argc, char *argv[]) {
  int natoms, optflags;
  void *v;
  int i, nvolsets, set;
  molfile_volumetric_t * meta;


  printf("got index: %d\n", lookup_keyword("	 ATOMS  "));
  
  while (--argc) {
    ++argv;
    v = open_xsf_read(*argv, "xsf", &natoms);
    if (!v) {
      fprintf(stderr, "open_xsf_read failed for file %s\n", *argv);
      return 1;
    }
    fprintf(stderr, "open_xsf_read succeeded for file %s\n", *argv);
    fprintf(stderr, "input contains %d atoms\n", natoms);

    molfile_atom_t atoms[natoms];
    
    // try reading the structure info
    i = read_xsf_structure(v, &optflags, atoms);
    if (i) {
      fprintf(stderr, "read_xsf_structure failed for file %s\n", *argv);
      return 1;
    }
    fprintf(stderr, "read_xsf_structure succeeded for file %s\n", *argv);

    // try loading the EDM metadata now
    if (read_xsf_metadata(v, &nvolsets, &meta)) {
      return 1; // failed to load xsf file
    }
    fprintf(stderr, "read_xsf_metadata succeeded for file %s\n", *argv);

    for (set=0; set<nvolsets; set++) {
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
      if (read_xsf_data(v, set, voldata, coldata)) {
        return 1; // failed to load xsf file
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

    molfile_timestep_t ts;
    ts.coords = new float[3*natoms];
    
    if (read_xsf_timestep(v, natoms, &ts)) {
      printf("read_xsf_timestep 0: failed\n");
    } else {
      printf("read_xsf_timestep 0: success\n");
    }
    
    if (read_xsf_timestep(v, natoms, &ts)) {
      printf("read_xsf_timestep 1: failed\n");
    } else {
      printf("read_xsf_timestep 1: success\n");
    }
    
    close_xsf_read(v);
  }
  return 0;
}

#endif

