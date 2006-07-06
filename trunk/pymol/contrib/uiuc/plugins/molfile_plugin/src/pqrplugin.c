/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_pqrplugin
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
 *      $RCSfile: pqrplugin.c,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.12 $       $Date: 2006/02/23 19:36:45 $
 *
 ***************************************************************************/

/*
 * PQR is be a space-delimited variant of the PDB format, with atom radius
 * and charge information replacing the B-factor and occupancy columns, and
 * containing only ATOM records.
 *
 * Most of this code is lifted from pdbplugin and readpdb.h
 */

#include "molfile_plugin.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PQR_RECORD_LENGTH 80

/*  record type defines */
enum {PQR_REMARK, PQR_ATOM, PQR_UNKNOWN, PQR_EOF, PQR_ERROR};

/* 
 * read the next record from the specified pqr file, and put the string
 * found in the given string pointer (the caller must provide adequate (81
 * chars) buffer space); return the type of record found
 */
static int read_pqr_record(FILE *f, char *retStr) {

  char inbuf[PQR_RECORD_LENGTH+2];
  int recType = PQR_UNKNOWN;
 
  /* XXX This PQR record reading code breaks with files that use
   * Mac or DOS style line breaks with ctrl-M characters.  We need
   * to replace the use of fgets() and comparisons against \n with
   * code that properly handles the other cases.
   */
 
  /*  read the next line  */
  if(inbuf != fgets(inbuf, PQR_RECORD_LENGTH+1, f)) {
    retStr[0] = '\0';
    if (feof(f)) {
      recType = PQR_EOF;
    }
    else {
      recType = PQR_ERROR;
    }
  } 
  else {
    /*  remove the newline character, if there is one */
    if(inbuf[strlen(inbuf)-1] == '\n')
      inbuf[strlen(inbuf)-1] = '\0';

    if (!strncmp(inbuf, "ATOM  ", 6)) {
      strcpy(retStr,inbuf);
      recType = PQR_ATOM;
    } 
    else {
      retStr[0] = '\0';
      recType = PQR_UNKNOWN;
    }
  }

  /* read the '\r', if there was one */
  {
    int ch = fgetc(f);
    if (ch != '\r') {
      ungetc(ch, f);
    }
  }
  
  return recType;
}

/*
 * Read the fields in an ATOM line. Return 0 on success, nonzero otherwise.
 */
static int get_pqr_fields(char *record, char *name, char *resname, char *chain,
    char *segname, char *resid,
    float *x, float *y, float *z, float *charge, float *radius) {

  /*
   * PQR files are like pdb files, but not quite.  They're fixed format
   * for the the identifier data up to resid, but after that it's free-format
   * for the floating point data.
   */
  strncpy(name, record+12, 4);    name[4] = 0;
  strncpy(resname, record+17, 4); resname[4] = 0;
  strncpy(chain, record+21, 1);   chain[1] = 0;
  strncpy(resid, record+22, 4);   resid[4] = 0;
  /* XXX what to do about the segid? */
  segname[0] = 0;

  if (sscanf(record+30, "%f%f%f%f%f", x, y, z, charge, radius) != 5) {
    return 1;
  }
  /* success */
  return 0;
}

/*
 * Read the coordinates in an ATOM line. Return 0 on success, nonzero otherwise.
 */
static int get_pqr_coordinates(char *record, float *x, float *y, float *z) {

  if (sscanf(record, "ATOM%*7d  %*4s%*4s%*5s %f%f%f%*f%*f",
             x, y, z) == 3) {
    return 0;
  }

  /* Return an error */
  else {
    return 1;
  }
}

/*
 * Return a default radius for the given atom.
 */
static float get_atom_radius(molfile_atom_t *atom) {
  char *name, *type;
  name = atom->name;
  type = atom->type;

  /* XXX - yeah, this needs to be fixed */
  return 1.0;
}

/*
 * Append a PQR ATOM record to the given file.
 */
static int write_pqr_atom(FILE *fd, int index, const molfile_atom_t *atom, 
                   float x, float y, float z) {
  int rc;

  rc = fprintf(fd, "ATOM  %5d %-4s %s %5d    %8.3f%8.3f%8.3f %.3f %.3f\n",
               index, atom->name, atom->resname, atom->resid, 
               x, y, z, atom->charge, atom->radius);

  return (rc > 0);
}


/*
 * API functions start here
 */

typedef struct {
  FILE *fd;
  int natoms;
  molfile_atom_t *atomlist;
} pqrdata;

static void *open_pqr_read(const char *filepath, const char *filetype, 
    int *natoms) {
  FILE *fd;
  pqrdata *pqr;
  char pqr_record[PQR_RECORD_LENGTH+2];
  int record_type;

  fd = fopen(filepath, "r");
  if (!fd) 
    return NULL;
  pqr = (pqrdata *)malloc(sizeof(pqrdata));
  pqr->fd = fd;

  *natoms=0;
  do {
    record_type = read_pqr_record(pqr->fd, pqr_record);
    if (record_type == PQR_ATOM)
      *natoms += 1;
    else if (record_type == PQR_ERROR) {
      fprintf(stderr, "Error reading PQR file after opening.\n");
      free(pqr);
      return NULL;
    }
    /* else, do nothing */
  } while (record_type != PQR_EOF);

  /* If no atoms were found, this is probably not a PQR file! */
  if (!*natoms) {
    fprintf(stderr, "PQR file '%s' contains no atoms.\n", filepath);
    free(pqr);
    return NULL;
  }

  rewind(pqr->fd);
  pqr->natoms = *natoms;
  return pqr; 
}

/* 
 * Read atom structure, but not coordinate information.
 */
static int read_pqr_structure(void *mydata, int *optflags, 
    molfile_atom_t *atoms) { 
  pqrdata *pqr = (pqrdata *)mydata;
  molfile_atom_t *atom;
  char pqr_record[PQR_RECORD_LENGTH + 1];
  int i, record_type;
  char ridstr[8];
  float newpos;
  long fpos = ftell(pqr->fd);
 
  *optflags = MOLFILE_CHARGE | MOLFILE_RADIUS;

  i = 0;
  do {
    record_type = read_pqr_record(pqr->fd, pqr_record);
    if ( (record_type == PQR_EOF) && (i < pqr->natoms) ) {
      fprintf(stderr, "Unexpected end-of-file.\n");
      return MOLFILE_ERROR;
    } 
    else if (record_type == PQR_ERROR) {
      fprintf(stderr, "Error reading atom coordinates.\n");
      return MOLFILE_ERROR;
    }
    else if (record_type == PQR_ATOM) {
      if (i >= pqr->natoms) {
        fprintf(stderr, "Too many atoms.\n");
        return MOLFILE_ERROR;
      }
      atom = atoms+i;
      get_pqr_fields(pqr_record, atom->name, atom->resname, atom->chain, 
          atom->segid, ridstr,
          &newpos, &newpos, &newpos, &atom->charge, &atom->radius);
      atom->resid = atoi(ridstr);
      strcpy(atom->type, atom->name);
      i++;
    }
    else {
      /* XXX do stuff */
    }
  } while(record_type != PQR_EOF);

  fseek(pqr->fd, fpos, SEEK_SET);

  return MOLFILE_SUCCESS;
}

/* 
 * Read atom coordinates. PQR files contain only a single "timestep".
 */
static int read_pqr_timestep(void *v, int natoms, molfile_timestep_t *ts) {
  pqrdata *pqr = (pqrdata *)v;
  char pqr_record[PQR_RECORD_LENGTH+2];
  int record_type, i;
  float *x, *y, *z;

  if (pqr->natoms == 0) 
    return MOLFILE_ERROR; /* EOF */
  if (ts) {
    x = ts->coords;
    y = x+1;
    z = x+2;
  } else {
    x = y = z = 0;
  }

  i = 0;
  do {
    record_type = read_pqr_record(pqr->fd, pqr_record);
    if ( (record_type == PQR_EOF) && (i < pqr->natoms) ) {
      fprintf(stderr, "Unexpected end-of-file.\n");
      return MOLFILE_ERROR;
    } 
    else if (record_type == PQR_ERROR) {
      fprintf(stderr, "Error reading atom coordinates.\n");
      return MOLFILE_ERROR;
    }
    else if (record_type == PQR_ATOM) {
      /* just get the coordinates, and store them */
      if (ts) {
        get_pqr_coordinates(pqr_record, x, y, z);
        x += 3;
        y += 3;
        z += 3;
        i++;
      } 
    }
    else {
      /* XXX do stuff */
    }
  } while(i < pqr->natoms);

  return MOLFILE_SUCCESS;
}

static void close_pqr_read(void *v) { 
  pqrdata *pqr = (pqrdata *)v;
  if (pqr->fd) {
    fclose(pqr->fd);
    pqr->fd = 0;
  }
  free(pqr);
}

static void *open_pqr_write(const char *path, const char *filetype, 
    int natoms) {
  FILE *fd;
  pqrdata *pqr;
  fd = fopen(path, "w");
  if (!fd) {
    fprintf(stderr, "Unable to open file %s for writing\n", path);
    return NULL;
  }

  pqr = (pqrdata *)malloc(sizeof(pqrdata));
  pqr->fd = fd;
  pqr->natoms = natoms; 
  pqr->atomlist = NULL;
  return pqr;
}
 
static int write_pqr_structure(void *v, int optflags, 
    const molfile_atom_t *atoms) {
  int i;
  pqrdata *pqr = (pqrdata *)v;
  int natoms = pqr->natoms;
  pqr->atomlist = (molfile_atom_t *)malloc(natoms*sizeof(molfile_atom_t));

  memcpy(pqr->atomlist, atoms, natoms*sizeof(molfile_atom_t));

  /* If charge and radius aren't given, we assign defaultvalues. */
  if (!(optflags & MOLFILE_CHARGE)) {
    for (i=0; i<natoms; i++) pqr->atomlist[i].charge = 0.0f;
  }
  if (!(optflags & MOLFILE_RADIUS)) {
    for (i=0; i<natoms; i++) pqr->atomlist[i].radius = 
      get_atom_radius(&(pqr->atomlist[i]));
  }

  return MOLFILE_SUCCESS;
}

static int write_pqr_timestep(void *v, const molfile_timestep_t *ts) {
  pqrdata *pqr = (pqrdata *)v; 
  const molfile_atom_t *atom;
  const float *pos;
  int i;

  if (pqr->natoms == 0)
    return MOLFILE_SUCCESS;

  atom = pqr->atomlist;
  pos = ts->coords;
  for (i=0; i<pqr->natoms; i++) {
    if (!write_pqr_atom(pqr->fd, i+1, atom, pos[0], pos[1], pos[2])) {
      fprintf(stderr, 
        "PQR: Error encoutered writing atom %d; file may be incomplete.\n", 
        i+1);
      return MOLFILE_ERROR;
    }
    atom++;
    pos += 3;
  }

  return MOLFILE_SUCCESS;
}
 
static void close_pqr_write(void *v) {
  pqrdata *pqr = (pqrdata *)v; 
  fclose(pqr->fd);
  free(pqr->atomlist);
  free(pqr);
}

/*
 * Initialization stuff down here
 */

static molfile_plugin_t plugin = {
  vmdplugin_ABIVERSION,         /* ABI version */
  MOLFILE_PLUGIN_TYPE,          /* type */
  "pqr",                        /* short name */
  "PQR",                        /* pretty name */
  "Eamon Caddigan",             /* author */
  0,                            /* major version */
  1,                            /* minor version */
  VMDPLUGIN_THREADSAFE,         /* is_reentrant */
  "pqr",                        /* filename extension */
  open_pqr_read,
  read_pqr_structure,
  0,
  read_pqr_timestep,
  close_pqr_read,
  open_pqr_write,
  write_pqr_structure,
  write_pqr_timestep,
  close_pqr_write
};
 
VMDPLUGIN_API int VMDPLUGIN_init() {
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_fini() {
  return VMDPLUGIN_SUCCESS;
}

