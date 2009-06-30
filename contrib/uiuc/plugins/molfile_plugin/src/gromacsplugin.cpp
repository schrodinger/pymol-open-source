/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_gromacsplugin
#define STATIC_PLUGIN 1

/***************************************************************************
 *cr
 *cr            (C) Copyright 1995-2009 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: gromacsplugin.C,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.46 $       $Date: 2009/04/29 15:45:30 $
 *
 ***************************************************************************/

#include "largefiles.h"   /* platform dependent 64-bit file I/O defines */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "Gromacs.h"
#include "molfile_plugin.h"

#if defined(_AIX)
#include <strings.h>
#endif

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

#if defined(WIN32) || defined(WIN64)
#define strcasecmp stricmp
#endif

typedef struct {
  md_file *mf;
  int natoms;
  int step;
} gmxdata;

static void *open_gro_read(const char *filename, const char *,
    int *natoms) {

    md_file *mf;
    md_header mdh;
    gmxdata *gmx;

    mf = mdio_open(filename, MDFMT_GRO);
    if (!mf) {
        fprintf(stderr, "gromacsplugin) Cannot open file '%s', %s\n",
                filename, mdio_errmsg(mdio_errno()));
        return NULL;
    }

    // read in the header data (careful not to rewind!)
    if (gro_header(mf, mdh.title, MAX_MDIO_TITLE,
    &mdh.timeval, &mdh.natoms, 0) < 0) {
        fprintf(stderr, "gromacsplugin) Cannot read header fromm '%s', %s\n",
                filename, mdio_errmsg(mdio_errno()));
            // XXX should free the file handle...
        return NULL;
    }
    *natoms = mdh.natoms;
    gmx = new gmxdata;
    gmx->mf = mf;
    gmx->natoms = mdh.natoms;
    return gmx;
}

static int read_gro_structure(void *mydata, int *optflags,
    molfile_atom_t *atoms) {

  md_atom ma;
  char buf[MAX_GRO_LINE + 1];
  gmxdata *gmx = (gmxdata *)mydata;
 
  *optflags = MOLFILE_NOOPTIONS; // no optional data

  // read in each atom and add it into the molecule
  for (int i = 0; i < gmx->natoms; i++) {
    molfile_atom_t *atom = atoms+i; 
    if (gro_rec(gmx->mf, &ma) < 0) {
      fprintf(stderr, "gromacsplugin) Error reading atom %d from file, %s\n", 
              i+1, mdio_errmsg(mdio_errno()));
      return MOLFILE_ERROR;
    }
    strcpy(atom->name, ma.atomname);
    strcpy(atom->type, ma.atomname);
    strcpy(atom->resname, ma.resname);
    atom->resid = atoi(ma.resid);
    atom->chain[0] = '\0';
    atom->segid[0] = '\0';
  }
  
  if (mdio_readline(gmx->mf, buf, MAX_GRO_LINE + 1, 0) < 0) {
    fprintf(stderr, "gromacsplugin) Warning, error reading box, %s\n", 
            mdio_errmsg(mdio_errno())); 
  }

  rewind(gmx->mf->f);
  return MOLFILE_SUCCESS;
}

static int read_gro_timestep(void *v, int natoms, molfile_timestep_t *ts) {
  gmxdata *gmx = (gmxdata *)v;
  md_ts mdts;
  memset(&mdts, 0, sizeof(md_ts));
  mdts.natoms = natoms;
 
  if (mdio_timestep(gmx->mf, &mdts) < 0) 
    return MOLFILE_ERROR;
  if (ts) {
    memcpy(ts->coords, mdts.pos, 3 * sizeof(float) * gmx->natoms);
    if (mdts.box) {
      ts->A = mdts.box->A;
      ts->B = mdts.box->B;
      ts->C = mdts.box->C;
      ts->alpha = mdts.box->alpha;
      ts->beta = mdts.box->beta;
      ts->gamma = mdts.box->gamma;
    }
  }
  mdio_tsfree(&mdts);
  return MOLFILE_SUCCESS;
}

static void close_gro_read(void *v) {
  gmxdata *gmx = (gmxdata *)v;
  mdio_close(gmx->mf);
  delete gmx;
}

static void *open_g96_read(const char *filename, const char *,
    int *natoms) {

    md_file *mf;
    md_header mdh;
    char gbuf[MAX_G96_LINE + 1];

    mf = mdio_open(filename, MDFMT_G96);
    if (!mf) {
        fprintf(stderr, "gromacsplugin) Cannot open file '%s', %s\n",
                filename, mdio_errmsg(mdio_errno()));
        return NULL;
    }
    
        // read in the header data
        if (g96_header(mf, mdh.title, MAX_MDIO_TITLE, &mdh.timeval) < 0) {
            fprintf(stderr, "gromacsplugin) Cannot read header from '%s', %s\n",
                    filename, mdio_errmsg(mdio_errno()));
            return NULL;
        }

        // First, look for a timestep block
        if (mdio_readline(mf, gbuf, MAX_G96_LINE + 1) < 0) {
            fprintf(stderr, "gromacsplugin) Cannot read header from '%s', %s\n",
                    filename, mdio_errmsg(mdio_errno()));
            return NULL;
        }
        if (!strcasecmp(gbuf, "TIMESTEP")) {
            // Read in the value line and the END line, and the next
            if (mdio_readline(mf, gbuf, MAX_G96_LINE + 1) < 0 ||
                mdio_readline(mf, gbuf, MAX_G96_LINE + 1) < 0 ||
                mdio_readline(mf, gbuf, MAX_G96_LINE + 1) < 0) {
              fprintf(stderr, "gromacsplugin) Cannot read header from '%s', %s\n",
                      filename, mdio_errmsg(mdio_errno()));
              return NULL;
            }
        }
        if (strcasecmp(gbuf, "POSITION") && strcasecmp(gbuf, "REFPOSITION")) {
          fprintf(stderr, "gromacsplugin) No structure information in file %s\n", filename);
          return NULL;
        }
        *natoms = g96_countatoms(mf); 
       
        gmxdata *gmx = new gmxdata;
        gmx->mf = mf;
        gmx->natoms = *natoms; 
        return gmx;
}

static int read_g96_structure(void *mydata, int *optflags,
    molfile_atom_t *atoms) {
  
    char gbuf[MAX_G96_LINE + 1];
    gmxdata *gmx = (gmxdata *)mydata;
    md_atom ma;
    md_file *mf = gmx->mf;
    
    *optflags = MOLFILE_NOOPTIONS; // no optional data

        for (int i = 0; i < gmx->natoms; i++) {
            molfile_atom_t *atom = atoms+i; 
            if (g96_rec(mf, &ma) < 0) {
                fprintf(stderr, "gromacsplugin) Error reading atom %d from file, %s\n",
                  i+1, mdio_errmsg(mdio_errno()));
                return MOLFILE_ERROR;
            }
            strcpy(atom->name, ma.atomname);
            strcpy(atom->type, ma.atomname);
            strcpy(atom->resname, ma.resname);
            atom->resid = atoi(ma.resid);
            atom->chain[0] = '\0';
            atom->segid[0] = '\0';
        }
        
        if (mdio_readline(mf, gbuf, MAX_G96_LINE + 1) < 0) {
            fprintf(stderr, "gromacsplugin) Warning, error reading END record, %s\n",
                mdio_errmsg(mdio_errno()));
        }

            // ... another problem: there may or may not be a VELOCITY
            // block or a BOX block, so we need to read one line beyond
            // the POSITION block to determine this. If neither VEL. nor
            // BOX are present we've read a line too far and infringed
            // on the next timestep, so we need to keep track of the
            // position now for a possible fseek() later to backtrack.
            long fpos = ftell(mf->f);

            // Now we must read in the velocities and the box, if present
            if (mdio_readline(mf, gbuf, MAX_G96_LINE + 1) >= 0) {

                // Is there a velocity block present ?
                if (!strcasecmp(gbuf, "VELOCITY") || !strcasecmp(gbuf, "VELOCITYRED")) {
                        // Ignore all the coordinates - VMD doesn't use them
                        for (;;) {
                                if (mdio_readline(mf, gbuf, MAX_G96_LINE + 1) < 0)
                                        return MOLFILE_ERROR;
                                if (!strcasecmp(gbuf, "END")) break;
                        }
    
                        // Again, record our position because we may need
                        // to fseek here later if we read too far.
                        fpos = ftell(mf->f);
    
                        // Go ahead and read the next line.
                        if (mdio_readline(mf, gbuf, MAX_G96_LINE + 1) < 0)
                    return MOLFILE_ERROR;
                }
    
                // Is there a box present ?
                if (!strcasecmp(gbuf, "BOX")) {
                        // Ignore the box coordinates at this time.
                        if (mdio_readline(mf, gbuf, MAX_G96_LINE + 1) < 0)
                    return MOLFILE_ERROR;
                        if (mdio_readline(mf, gbuf, MAX_G96_LINE + 1) < 0)
                    return MOLFILE_ERROR;
                        if (strcasecmp(gbuf, "END")) 
                    return MOLFILE_ERROR;
                }
                else {
                        // We have read too far, so fseek back to the
                        // last known safe position so we don't return
                        // with the file pointer set infringing on the
                        // next timestep data.
                        fseek(mf->f, fpos, SEEK_SET);
                }
        }
        else {
            // Go ahead and rewind for good measure
            fseek(mf->f, fpos, SEEK_SET);
        }
        rewind(mf->f);
        return MOLFILE_SUCCESS;
}

static int read_g96_timestep(void *v, int natoms, molfile_timestep_t *ts) {

  gmxdata *gmx = (gmxdata *)v;
  md_ts mdts;
  memset(&mdts, 0, sizeof(md_ts));
  mdts.natoms = natoms;

  if (mdio_timestep(gmx->mf, &mdts) < 0) 
    return MOLFILE_ERROR;
  if (ts) {
    memcpy(ts->coords, mdts.pos, 3 * sizeof(float) * gmx->natoms);
    if (mdts.box) {
      ts->A = mdts.box->A;
      ts->B = mdts.box->B;
      ts->C = mdts.box->C;
      ts->alpha = mdts.box->alpha;
      ts->beta = mdts.box->beta;
      ts->gamma = mdts.box->gamma;
    }
  }
  mdio_tsfree(&mdts);
  return MOLFILE_SUCCESS;
}

static void close_g96_read(void *v) {
  gmxdata *gmx = (gmxdata *)v;
  mdio_close(gmx->mf);
  delete gmx;
}
  

//
// TRR and XTC files
//

static void *open_trr_read(const char *filename, const char *filetype,
    int *natoms) {

    md_file *mf;
    md_header mdh;
    gmxdata *gmx;
    int format;
  
    if (!strcmp(filetype, "trr")) 
      format = MDFMT_TRR;
    else if (!strcmp(filetype, "trj")) 
      format = MDFMT_TRJ;
    else if (!strcmp(filetype, "xtc")) 
      format = MDFMT_XTC;
    else
      return NULL;

    mf = mdio_open(filename, format);
    if (!mf) {
        fprintf(stderr, "gromacsplugin) Cannot open file '%s', %s\n",
                filename, mdio_errmsg(mdio_errno()));
        return NULL;
    }
    if (mdio_header(mf, &mdh) < 0) {
        mdio_close(mf);
        fprintf(stderr, "gromacsplugin) Cannot read header fromm '%s', %s\n",
                filename, mdio_errmsg(mdio_errno()));
        return NULL;
    }
    *natoms = mdh.natoms;
    gmx = new gmxdata;
    gmx->mf = mf;
    gmx->natoms = mdh.natoms;
    return gmx;
}

static int read_trr_timestep(void *v, int natoms, molfile_timestep_t *ts) {
  gmxdata *gmx = (gmxdata *)v;
  md_ts mdts;
  memset(&mdts, 0, sizeof(md_ts));
  mdts.natoms = natoms;

  if (mdio_timestep(gmx->mf, &mdts) < 0) {
    if (mdio_errno() == MDIO_EOF || mdio_errno() == MDIO_IOERROR) {
      // XXX Lame, why does mdio treat IOERROR like EOF?  
      return MOLFILE_ERROR;
    }
    fprintf(stderr, "gromacsplugin) Error reading timestep, %s\n", 
            mdio_errmsg(mdio_errno()));
    return MOLFILE_ERROR;
  }
  if (mdts.natoms != natoms) {
    fprintf(stderr, "gromacsplugin) Timestep in file contains wrong number of atoms\n");
    fprintf(stderr, "gromacsplugin) Found %d, expected %d\n", mdts.natoms, natoms);
    mdio_tsfree(&mdts);
    return MOLFILE_ERROR;
  }
   
  if (ts) {
    memcpy(ts->coords, mdts.pos, 3 * sizeof(float) * gmx->natoms);
    if (mdts.box) {
      ts->A = mdts.box->A;
      ts->B = mdts.box->B;
      ts->C = mdts.box->C;
      ts->alpha = mdts.box->alpha;
      ts->beta = mdts.box->beta;
      ts->gamma = mdts.box->gamma;
    }
  }
  mdio_tsfree(&mdts);
  return MOLFILE_SUCCESS;
}

static void close_trr_read(void *v) {
  gmxdata *gmx = (gmxdata *)v;
  mdio_close(gmx->mf);
  delete gmx;
}

// open file for writing
static void *open_trr_write(const char *filename, const char *filetype,
    int natoms) {

    md_file *mf;
    gmxdata *gmx;
    int format;
  
    if (!strcmp(filetype, "trr")) 
      format = MDFMT_TRR;
    else if (!strcmp(filetype, "xtc")) 
      format = MDFMT_XTC;
    else
      return NULL;

    mf = mdio_open(filename, format, MDIO_WRITE);
    if (!mf) {
        fprintf(stderr, "gromacsplugin) Cannot open file '%s', %s\n",
                filename, mdio_errmsg(mdio_errno()));
        return NULL;
    }
    gmx = new gmxdata;
    gmx->mf = mf;
    gmx->natoms = natoms;
    // set some parameters for the output stream:
    // start at step 0, convert to big-endian, write single precision.
    gmx->step   = 0;
    gmx->mf->rev = host_is_little_endian();
    gmx->mf->prec = sizeof(float);
    return gmx;
}

// write a trr timestep. the file format has a header with each record
static int write_trr_timestep(void *mydata, const molfile_timestep_t *ts) 
{
  const float nm=0.1;
  
  gmxdata *gmx = (gmxdata *)mydata;
  
  // determine and write header from structure info.
  // write trr header. XXX: move this to Gromacs.h ??
  if (gmx->mf->fmt == MDFMT_TRR) {
    int i;

    if ( put_trx_int(gmx->mf, TRX_MAGIC)            // ID
         || put_trx_string(gmx->mf, "GMX_trn_file") // version
         || put_trx_int(gmx->mf, 0)                 // ir_size (ignored)
         || put_trx_int(gmx->mf, 0)                 // e_size (ignored)
         || put_trx_int(gmx->mf, 9*sizeof(float))   // box
         || put_trx_int(gmx->mf, 0)                 // vir_size (ignored)
         || put_trx_int(gmx->mf, 0)                 // pres_size (ignored)
         || put_trx_int(gmx->mf, 0)                 // top_size (ignored)
         || put_trx_int(gmx->mf, 0)                 // sym_size (ignored)
         || put_trx_int(gmx->mf, 3*sizeof(float)*gmx->natoms) // coordinates 
         || put_trx_int(gmx->mf, 0)                 // no velocities
         || put_trx_int(gmx->mf, 0)                 // no forces
         || put_trx_int(gmx->mf, gmx->natoms)       // number of atoms
         || put_trx_int(gmx->mf, gmx->step)         // current step number
         || put_trx_int(gmx->mf, 0)                 // nre (ignored)
         || put_trx_real(gmx->mf, 0.1*gmx->step)    // current time. (dummy value: 0.1)
         || put_trx_real(gmx->mf, 0.0))             // current lambda
      return MOLFILE_ERROR;

    // set up box according to the VMD unitcell conventions.
    // the a-vector is collinear with the x-axis and
    // the b-vector is in the xy-plane. 
    const float sa = sin((double)ts->alpha/180.0*M_PI);
    const float ca = cos((double)ts->alpha/180.0*M_PI);
    const float cb = cos((double)ts->beta/180.0*M_PI);
    const float cg = cos((double)ts->gamma/180.0*M_PI);
    const float sg = sin((double)ts->gamma/180.0*M_PI);
    float box[9];
    box[0] = ts->A;    box[1] = 0.0;      box[2] = 0.0;
    box[3] = ts->B*ca; box[4] = ts->B*sa; box[5] = 0.0;
    box[6] = ts->C*cb; box[7] = ts->C*(ca - cb*cg)/sg;
    box[8] = ts->C*sqrt((double)(1.0 + 2.0*ca*cb*cg 
                                 - ca*ca - cb*cb - cg*cg)/(1.0 - cg*cg));

    for (i=0; i<9; ++i) {
      if (put_trx_real(gmx->mf, box[i]*nm))
        return MOLFILE_ERROR;
    }
#ifdef TEST_TRR_PLUGIN
    fprintf(stderr, "gromacsplugin) box is:\n %f %f %f\n %f %f %f\n %f %f %f\n\n",
            box[0], box[1], box[2], box[3], box[4], box[5], box[6], box[7], box[8]);
#endif

    // write coordinates
    for (i=0; i<(3*gmx->natoms); ++i) {
      if (put_trx_real(gmx->mf, ts->coords[i]*nm))
        return MOLFILE_ERROR;
    }
  } else {
    fprintf(stderr, "gromacsplugin) only .trr is supported for writing\n");
    return MOLFILE_ERROR;
  }

  ++ gmx->step;
  return MOLFILE_SUCCESS;
  }


static void close_trr_write(void *v) {
  gmxdata *gmx = (gmxdata *)v;
  mdio_close(gmx->mf);
  delete gmx;
}

#define GROMACS_PLUGIN_MAJOR_VERSION 0
#define GROMACS_PLUGIN_MINOR_VERSION 9

//
// plugin registration stuff below
//
static molfile_plugin_t gro_plugin = {
  vmdplugin_ABIVERSION,               // ABI version
  MOLFILE_PLUGIN_TYPE,                // type of plugin
  "gro",                              // short name of plugin
  "Gromacs GRO",                      // pretty name of plugin
  "David Norris, Justin Gullingsrud", // authors
  GROMACS_PLUGIN_MAJOR_VERSION,       // major version
  GROMACS_PLUGIN_MINOR_VERSION,       // minor version
  VMDPLUGIN_THREADUNSAFE,             // is not reentrant
  "gro",                              // filename extension
  open_gro_read,
  read_gro_structure,
  0,
  read_gro_timestep,
  close_gro_read,
  0,                                  // open_write
  0,                                  // write_structure
  0,                                  // write_timestep
  0,                                  // close_write
  0,                                  // read_volumetric_metadata
  0,                                  // read_volumetric_data
  0                                   // read_rawgraphics
};

static molfile_plugin_t g96_plugin = {
  vmdplugin_ABIVERSION,               // ABI version
  MOLFILE_PLUGIN_TYPE,                // type of plugin
  "g96",                              // short name of plugin
  "Gromacs g96",                      // pretty name of plugin
  "David Norris, Justin Gullingsrud", // authors
  GROMACS_PLUGIN_MAJOR_VERSION,       // major version
  GROMACS_PLUGIN_MINOR_VERSION,       // minor version
  VMDPLUGIN_THREADUNSAFE,             // is not reentrant
  "g96",                              // filename extension
  open_g96_read,
  read_g96_structure,
  0,
  read_g96_timestep,
  close_g96_read,
  0,                                  // open_write
  0,                                  // write_structure
  0,                                  // write_timestep
  0,                                  // close_write
  0,                                  // read_volumetric_metadata
  0,                                  // read_volumetric_data
  0                                   // read_rawgraphics
};

static molfile_plugin_t trr_plugin = {
  vmdplugin_ABIVERSION,               // ABI version
  MOLFILE_PLUGIN_TYPE,                // type of plugin
  "trr",                              // short name of plugin
  "Gromacs TRR Trajectory",           // pretty name of plugin
  "David Norris, Justin Gullingsrud, Axel Kohlmeyer", // authors
  GROMACS_PLUGIN_MAJOR_VERSION,       // major version
  GROMACS_PLUGIN_MINOR_VERSION,       // minor version
  VMDPLUGIN_THREADUNSAFE,             // is not reentrant
  "trr",                              // filename extension
  open_trr_read,
  0,
  0,
  read_trr_timestep,
  close_trr_read,
  open_trr_write,
  0,                            // write_structure
  write_trr_timestep,
  close_trr_write,
  0,                            // read_volumetric_metadata
  0,                            // read_volumetric_data
  0                             // read_rawgraphics
};

static molfile_plugin_t xtc_plugin = {
  vmdplugin_ABIVERSION,                // ABI version
  MOLFILE_PLUGIN_TYPE,                 // type of plugin
  "xtc",                               // short name of plugin
  "Gromacs XTC Compressed Trajectory", // pretty name of plugin
  "David Norris, Justin Gullingsrud",  // authors
  GROMACS_PLUGIN_MAJOR_VERSION,        // major version
  GROMACS_PLUGIN_MINOR_VERSION,        // minor version
  VMDPLUGIN_THREADUNSAFE,              // is not reentrant
  "xtc",                               // filename extension
  open_trr_read,
  0,
  0,
  read_trr_timestep,
  close_trr_read,
  0,                                  // open_write
  0,                                  // write_structure
  0,                                  // write_timestep
  0,                                  // close_write
  0,                                  // read_volumetric_metadata
  0,                                  // read_volumetric_data
  0                                   // read_rawgraphics
};

static molfile_plugin_t trj_plugin = {
  vmdplugin_ABIVERSION,                // ABI version
  MOLFILE_PLUGIN_TYPE,                 // type of plugin
  "trj",                               // short name of plugin
  "Gromacs TRJ Trajectory", // pretty name of plugin
  "David Norris, Justin Gullingsrud",  // authors
  GROMACS_PLUGIN_MAJOR_VERSION,        // major version
  GROMACS_PLUGIN_MINOR_VERSION,        // minor version
  VMDPLUGIN_THREADUNSAFE,              // is not reentrant
  "trj",                               // filename extension
  open_trr_read,
  0,
  0,
  read_trr_timestep,
  close_trr_read,
  0,                                  // open_write
  0,                                  // write_structure
  0,                                  // write_timestep
  0,                                  // close_write
  0,                                  // read_volumetric_metadata
  0,                                  // read_volumetric_data
  0                                   // read_rawgraphics
};

VMDPLUGIN_EXTERN int VMDPLUGIN_init() { 
  return 0; 
}

VMDPLUGIN_EXTERN int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)&gro_plugin);
  (*cb)(v, (vmdplugin_t *)&g96_plugin);
  (*cb)(v, (vmdplugin_t *)&trr_plugin);
  (*cb)(v, (vmdplugin_t *)&trj_plugin);
  (*cb)(v, (vmdplugin_t *)&xtc_plugin);
  return 0;
}

VMDPLUGIN_EXTERN int VMDPLUGIN_fini() {
  return 0;
}

  
#ifdef TEST_G96_PLUGIN

int main(int argc, char *argv[]) {
  int natoms;
  
  molfile_timestep_t timestep;
  void *v;
  int i;

  if (argc < 2) return 1;
  while (--argc) {
    ++argv;
    v = open_g96_read(*argv, "g96", &natoms);
    if (!v) {
      fprintf(stderr, "open_g96_read failed for file %s\n", *argv);
      return 1;
    }
    timestep.coords = (float *)malloc(3*sizeof(float)*natoms);
    i = 0;
    while(!read_g96_timestep(v, natoms, &timestep)) {
      ++i;
    }
    fprintf(stderr, "ended read_g96_timestep on step %d\n", i);
    free(timestep.coords);
    close_g96_read(v);
  }
  return 0;
}

#endif

#ifdef TEST_TRR_PLUGIN

int main(int argc, char *argv[]) {
  int natoms;
  
  molfile_timestep_t timestep;
  void *v, *w;
  int i;

  if (argc != 3) return 1;
  v = open_trr_read(argv[1], "trr", &natoms);
  if (!v) {
    fprintf(stderr, "open_trr_read failed for file %s\n", argv[1]);
    return 1;
  }
  timestep.coords = (float *)malloc(3*sizeof(float)*natoms);
  w = open_trr_write(argv[2], "trr", natoms);
  if (!w) {
    fprintf(stderr, "open_trr_write failed for file %s\n", argv[2]);
    return 1;
  }
  
  i = 0;
  while(!read_trr_timestep(v, natoms, &timestep)) {
    ++i;
    if (write_trr_timestep(w, &timestep)) {
      fprintf(stderr, "write error\n");
      return 1;
    }
  }
  
  fprintf(stderr, "ended read_trr_timestep on step %d\n", i);
  free(timestep.coords);
  close_trr_read(v);
  close_trr_write(w);
  return 0;
}

#endif

