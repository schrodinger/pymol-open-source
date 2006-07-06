/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_carplugin
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
 *      $RCSfile: carplugin.c,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.12 $       $Date: 2006/02/23 19:36:44 $
 *
 ***************************************************************************/

/*
 * Plugin for Insight II/Discover car (cartesian coordinate file) file format.
 * 
 * TODO: 
 * + 2D PBC info (probably just as a simplified 3D PBC)
 * + HELIX info (not sure how this will be handled)
 *
 * File Header:
 *   !BIOSYM archive 3
 *   HELIX 
 *     (if HELIX information is present)
 *   PBC=ON, PBC=OFF, PBC=2D 
 *     (one of these three choices must be present. PBC cannot be "ON" if
 *     HELIX information is present.)
 * Coordinate Header:
 *   Title (this line may be blank, but must be present)
 *     1-64 title for the system
 *     65-80 energy
 *   ! DATE day month date time year
 *     (day, month, date, time, year are optional)
 *   PBC information if PBC=ON:
 *     1-3 PBC
 *     4-13 a cell vector a in angstroms
 *     14-23 b cell vector b in angstroms
 *     24-33 c cell vector c in angstroms
 *     34-43 alpha cell angle alpha in degrees
 *     44-53 beta cell angle beta in degrees
 *     54-63 gamma cell angle gamma in degrees
 *     64-80 space group name
 *   PBC information if PBC=2D:
 *     1-3 PBC
 *     4-13 k plane vector k in angstroms
 *     14-23 l plane vector l in angstroms
 *     24-33 gamma plane angle gamma in degrees
 *     34-50 plane group name
 * Molecule Data:
 *   If helix info is present:
 *     1-5 HELIX
 *     6-15 sigma in degrees
 *     16-25 d in angstroms
 *     26-35 kappa angle between l axis and helix axis in degrees
 *     36-45 lambda angle between k axis and helix axis in degrees
 *     46-55 Tk fractional position of helix axis along k axis
 *     56-65 Tl fractional position of helix axis along l axis
 *   Atom data:
 *     1-5 atom name
 *     7-20 x Cartesian coordinate of atom in angstroms
 *     22-35 y Cartesian coordinate of atom in angstroms
 *     37-50 z Cartesian coordinate of atom in angstroms
 *     52-55 type of residue containing atom
 *     57-63 residue sequence name relative to beginning of current
 *     molecule, left justified
 *     64-70 potential type of atom left justified
 *     72-73 element symbol
 *     75-80 partial charge on atom
 *   Final line for a given molecule:
 *     1-3 end
 * Final line for the entire molecular system input:
 *     1-3 end
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "molfile_plugin.h"

#define LINESIZE 1024

/* possible values for pbc */
enum {PBC_ON, PBC_OFF, PBC_2D};

typedef struct {
  FILE *file;
  int numatoms, pbc, helix, eof;
  long coord_location;
  molfile_atom_t *atomlist;
} cardata;

/* Parse a line contianing atom data from a car file and store it in the
 * atom structure. Returns 1 on success, 0 on failure.
 */
static int read_car_structure_line(molfile_atom_t *atom, char *line) {
  char name[LINESIZE], type[LINESIZE];
  int resid;
  float charge;

  if (sscanf(line, "%s %*f %*f %*f %*s %d %*s %s %f", name, &resid, type, &charge) 
      != 4)
    return 0;

  /* check the length of the name and type strings before copying. */
  if ( (strlen(name) >= 8) || (strlen(type) >= 8) )
    return 0;
  strcpy(atom->name, name);
  strcpy(atom->type, type);

  /* Check that the resid won't overflow the resname string, then copy it
   * over
   */
  if (resid > 9999999)
    atom->resname[0] = '\0';
  else
    sprintf(atom->resname, "%d", resid);

  atom->resid = resid;

  atom->chain[0] = '\0';
  atom->segid[0] = '\0';

  atom->charge = charge;

  return 1;
}

/* Parse a line contianing atom coordinates from a car file and store them
 * in the array 'coords' in XYZ order. Returns 1 on success, 0 on failure.
 */
static int read_car_coordinates(float *coords, char *line) {
  float x, y, z;

  if (sscanf(line, "%*s %f %f %f %*s %*d %*s %*s %*f", &x, &y, &z) 
      != 3)
    return 0;

  coords[0] = x;
  coords[1] = y;
  coords[2] = z;

  return 1;
}

static void *open_car_read(const char *filename, const char *filetype, 
                           int *natoms) {
  FILE *fd;
  cardata *data;
  char line[LINESIZE];

  fd = fopen(filename, "rb"); if (!fd) return NULL;
  
  data = (cardata *)malloc(sizeof(cardata));
  data->eof = 0;
  data->file = fd;

  /* First line: "!BIOSYM archive n", where n indicates the file format */
  fgets(line, LINESIZE, fd);
  if (strncmp(line, "!BIOSYM archive", 15)) {
    fprintf(stderr, "ERROR) badly formatted/missing header.\n");
    return NULL;
  }

  /* Second line: "HELIX" if helix information is present. Followed by PBC
   * info on the next line: "PBC=ON|OFF|2D". 
   * If HELIX is present, PBC cannot be "ON": return an error if this happens. 
   */
  fgets(line, LINESIZE, fd);
  if (!strncmp(line, "HELIX", 5)) {
    data->helix = 1;
    fgets(line, LINESIZE, fd);
    fprintf(stdout, "WARNING) ignoring helix information.\n");
  }
  else {
    data->helix = 0;
  }

  if (!strncmp(line, "PBC=ON", 6)) {
    data->pbc = PBC_ON;
  }
  else if (!strncmp(line, "PBC=OFF", 7)) {
    data->pbc = PBC_OFF;
  }
  else if (!strncmp(line, "PBC=2D", 6)) {
    data->pbc = PBC_2D;
    fprintf(stdout, "WARNING) ignoring 2D PBC information.\n");
  }
  else {
    fprintf(stderr, "ERROR) badly formatted/missing PBC info.\n");
    return NULL;
  }

  if (data->helix && (data->pbc == PBC_ON)) {
    fprintf(stderr, "ERROR) car file contains helix and 3D PBC information.");
    return NULL;
  }

  /* Next line: title/energy for the system. Skipped. */
  fgets(line, LINESIZE, fd);

  /* Next line: "!DATE [day month date time year]". */
  fgets(line, LINESIZE, fd);
  if (strncmp(line, "!DATE", 5)) {
    fprintf(stderr, "ERROR) badly formatted/missing date.\n");
    return NULL;
  }
    
  /* Store the location of the beginning of the PBC/coordinate data. */
  data->coord_location = ftell(fd);

  /* Skip the PBC and HELIX entries, if present */
  if (data->pbc != PBC_OFF) 
    fgets(line, LINESIZE, fd);
  if (data->helix)
    fgets(line, LINESIZE, fd);

  /* Count the atoms in all molecules*/
  data->numatoms = 0;
  fgets(line, LINESIZE, fd);
  while(strncmp(line, "end", 3)) {
    while(strncmp(line, "end", 3)) {
      data->numatoms++;
      fgets(line, LINESIZE, fd);

      if (feof(fd)) {
        fprintf(stderr, "ERROR) unexpected end-of-file.\n");
        return NULL;
      }
      if (ferror(fd)) {
        fprintf(stderr, "ERROR) error reading car file.\n");
        return NULL;
      }
    }
    fgets(line, LINESIZE, fd);
  }
  *natoms = data->numatoms;

  return data;
}

static int read_car_structure(void *mydata, int *optflags, 
                              molfile_atom_t *atoms) {
  int mol_num;
  char line[LINESIZE];
  molfile_atom_t *atom;
  cardata *data = (cardata *)mydata;

  *optflags = MOLFILE_CHARGE; /* car files contain partial charges */

  /* move to the beginning of the atom data in the file, skipping any PBC or
   * HELIX information that may be present
   */
  fseek(data->file, data->coord_location, SEEK_SET);
  if (data->pbc != PBC_OFF) 
    fgets(line, LINESIZE, data->file);
  if (data->helix)
    fgets(line, LINESIZE, data->file);

  mol_num = 0;
  atom = atoms;
  fgets(line, LINESIZE, data->file);
  /* Loop through all molecules */
  while(strncmp(line, "end", 3)) {
    /* Read the structure for each molecule */
    while(strncmp(line, "end", 3)) {
      if (!read_car_structure_line(atom, line)) {
        fprintf(stderr, "ERROR) badly formatted structure line:\n%s\n", line);
        return MOLFILE_ERROR;
      }

      /* XXX - use the chain name to identify different molecules */
      sprintf(atom->chain, "%d", mol_num);

      atom++;
      
      fgets(line, LINESIZE, data->file);
      if (feof(data->file)) {
        fprintf(stderr, "ERROR) unexpected end-of-file while reading structure.\n");
        return MOLFILE_ERROR;
      }
      if (ferror(data->file)) {
        fprintf(stderr, "ERROR) error reading car file while reading structure.\n");
        return MOLFILE_ERROR;
      }
    }
    fgets(line, LINESIZE, data->file);
    mol_num++;
  }

  return MOLFILE_SUCCESS;
}

static int read_car_timestep(void *mydata, int natoms, molfile_timestep_t *ts) {
  char line[LINESIZE];
  cardata *data = (cardata *)mydata;
  float *coords = NULL;

  /* return if coordinates have been read */
  if (data->eof)
    return MOLFILE_EOF;

  /* move to the beginning of the atom data in the file */
  fseek(data->file, data->coord_location, SEEK_SET);

  /* read the PBC record if present */
  if (data->pbc == PBC_ON) {
    fgets(line, LINESIZE, data->file);

    if (ts) {
      if ( sscanf(line, "PBC %f %f %f %f %f %f %*s", 
                  &ts->A, &ts->B, &ts->C, 
                  &ts->alpha, &ts->beta, &ts->gamma) != 6 ) {
        fprintf(stderr, "ERROR) badly formatted PBC line:\n%s\n", line);
        return MOLFILE_ERROR;
      }
    }
  }
  else if (data->pbc == PBC_2D) {
    /* XXX - Ignore 2D PBC information */
    fgets(line, LINESIZE, data->file);
  }

  /* skip the helix record if present */
  if (data->helix)
    fgets(line, LINESIZE, data->file);

  if (ts)
    coords = ts->coords;

  fgets(line, LINESIZE, data->file);
  /* Loop through all molecules */
  while(strncmp(line, "end", 3)) {
    /* Read the coordinates for each molecule */
    while(strncmp(line, "end", 3)) {
      /* only save coords if we're given a timestep pointer. */
      if (ts) {
        if (!read_car_coordinates(coords, line)) {
          fprintf(stderr, "ERROR) badly formatted coordinate line:\n%s\n", line);
          return MOLFILE_ERROR;
        }
        coords += 3;
      }

      fgets(line, LINESIZE, data->file);
      if (feof(data->file)) {
        fprintf(stderr, "ERROR) unexpected end-of-file while reading coordinates.\n");
        return MOLFILE_ERROR;
      }
      if (ferror(data->file)) {
        fprintf(stderr, "ERROR) file error while reading coordinates.\n");
        return MOLFILE_ERROR;
      }
    }
    fgets(line, LINESIZE, data->file);
  }

  /* set eof since this file contians only one "timestep" */
  data->eof = 1;

  return MOLFILE_SUCCESS;
}
    
static void close_car_read(void *mydata) {
  if (mydata) {
    cardata *data = (cardata *)mydata;
    fclose(data->file);
    free(data);
  }
}


/* registration stuff */
static molfile_plugin_t carplugin = {
  vmdplugin_ABIVERSION,
  MOLFILE_PLUGIN_TYPE,                         /* type */
  "car",                                       /* short name */
  "InsightII car",                             /* pretty name */
  "Eamon Caddigan",                            /* author */
  0,                                           /* major version */
  3,                                           /* minor version */
  VMDPLUGIN_THREADSAFE,                        /* is reentrant */
  "car",
  open_car_read,
  read_car_structure,
  0,
  read_car_timestep,
  close_car_read,
};

VMDPLUGIN_API int VMDPLUGIN_init() {
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)&carplugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_fini() {
  return VMDPLUGIN_SUCCESS;
}

