/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_corplugin
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
 *      $RCSfile: corplugin.c,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.29 $       $Date: 2009/04/29 15:45:28 $
 *
 ***************************************************************************/

/*
 * This plugin reads molecular coordinate data stored in 
 * CHARMM Cartesian Coordinate format (ascii text format, not binary).
 *    http:

 *    http:

 */

#include "molfile_plugin.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define COR_RECORD_LENGTH   141 /* 80 */

/* Remove leading and trailing whitespace from the string str of size n */
static void strip_whitespace(char *str, int n) {
  char *beg, *end;

  beg = str;
  end = str + (n-2); /* Point to the last non-null character in the string */

  /* Remove leading whitespace */
  while(beg <= end && *beg == ' ') {
    beg++;
  }

  /* Remove trailing whitespace */
  while(end >= str && *end == ' ') {
    end--;
  }

  if (beg < end) {
    /* Shift the string */
    *(end+1) = '\0';
    memmove(str, beg, (end - beg + 2));
  } else {
    str[0] = '\0';
  }

  return;
}

/* Get a string from a stream, printing any errors that occur */
static char *corgets(char *s, int n, FILE *stream) {
  char *returnVal;

  if (feof(stream)) {
    printf("corplugin) Unexpected end-of-file.\n");
    returnVal = NULL;
  } else if (ferror(stream)) {
    printf("corplugin) Error reading file.\n");
    return NULL;
  } else {
    returnVal = fgets(s, n, stream);
    if (returnVal == NULL) {
      printf("corplugin) Error reading line.\n");
    }
  }

  return returnVal;
}


/* Open the .cor file and skip past the remarks to the first data section.
 * Returns the file pointer, or NULL if error.  Also puts the number of
 * atoms in the molecule into the given integer.  
 */
static FILE *open_cor_file(const char *fname, int *natom, int *iofoext) {
  char inbuf[COR_RECORD_LENGTH+2], header[11];
  FILE *f;

  *natom = MOLFILE_NUMATOMS_NONE;

  if (!fname) {
    printf("corplugin) Error opening file: no filename given.\n");
    return NULL;
  }

  if ((f = fopen(fname, "r")) == NULL) {
    printf("corplugin) Error opening file.\n");
    return NULL;
  }

  /* Read and discard the header */
  do {
    if (fgets(inbuf, COR_RECORD_LENGTH+1, f) == NULL) {
      fclose(f);
      printf("corplugin) Error opening file: cannot read line.\n");
      return NULL;
    }

    if (sscanf(inbuf, "%10c", header) != 1) {
      fclose(f);
      printf("corplugin) Error opening file: improperly formatted line.\n");
      return NULL;
    }

  } while (header[0]=='*');

  /* check for EXT keyword */
  *iofoext = 0 ;
  if (strstr(inbuf,"EXT") != NULL) 
    *iofoext = 1;

  /* check atom count */
  header[10] = '\0';
  *natom = atoi(header);
  if (*natom > 99999) 
    *iofoext = 1;

  if (*iofoext == 1)
    printf("corplugin) Using EXTended CHARMM coordinates file\n");

  return f;
}

/* Read in the next atom info into the given storage areas; this assumes
   that file has already been moved to the beginning of the atom records.
   Returns the serial number of the atom. If there is an error, returns -1.*/
static int get_cor_atom(FILE *f, char *atomName, char *atomType, char
    *resName, char *segName, int *resId, int ioext) {
  char inbuf[COR_RECORD_LENGTH+2], numAtomStr[11], resNStr[11], resIdStr[11];
  char atomNameStr[11], segNameStr[11], resNameStr[11];
  int numAtom;

  if (corgets(inbuf, COR_RECORD_LENGTH+1, f) == NULL) {
    return -1;
  }

  if (strlen(inbuf) < 60) {
    printf("corplugin) Line too short: \n%s\n", inbuf);
    return -1;
  }

  memset(numAtomStr, 0, sizeof(numAtomStr));
  memset(resNStr, 0, sizeof(resNStr));
  memset(resIdStr, 0, sizeof(resIdStr));
  memset(resNameStr, 0, sizeof(resNameStr));
  memset(segNameStr, 0, sizeof(segNameStr));
  memset(atomNameStr, 0, sizeof(atomNameStr));

  /*

    CHARMM has a variety of input formats for the
       read coor card
    command. In order to support all of them the simplest would be to replace 

    if (sscanf(inbuf, "%5c%5c%5c%5c%*10c%*10c%*10c%5c", 

    with

    if (sscanf(inbuf, "%s %s %s %s %*s %*s %*s %s", 

    However this solution has 2 problems:
     1. buffer overruns
     2. it doesn't handle the cases where X,Y,Z values have no spaces in-between

    In this fix we handle only two most important cases, depending on
     the value of IOFOrmat command (either EXTEnded or NOEXtended):
     EXT:    2I10,2X,A8,2X,A8,3F20.10,2X,A8,2X,A8,F20.10
     NOEXT:  2I5,1X,A4,1X,A4,3F10.5,1X,A4,1X,A4,F10.5

    This implementation introduces new flag iofoext in cordata
    structure, which can choose between the 2 formats. Note that
    CHARMM adds EXT keyword in the coordinate and psf files in the
    case of IOFO EXTE command in the input script, or automatically
    when there are 100,000 or more atoms in the system!

  */

  if (ioext == 1 ) {
    if (sscanf(inbuf, "%10c%10c%10c%10c%*20c%*20c%*20c%10c%10c",
               numAtomStr, resNStr, resNameStr, atomNameStr, segNameStr, resIdStr) != 6) {
      printf("corplugin) Improperly formatted line: \n%s\n", inbuf);
      return -1;
    }
    strip_whitespace(resName, sizeof(resName));  /* strip from long original */
    strip_whitespace(atomName, sizeof(atomName));
    strip_whitespace(segName, sizeof(segName));
    memcpy(resName, resNameStr, 7);              /* XXX truncate extra chars */
    memcpy(atomName, atomNameStr, 7);
    memcpy(segName, segNameStr, 7);
    resName[7] = '\0';                           /* NUL terminate strings.. */
    atomName[7] = '\0';
    segName[7] = '\0';
  } else {
    if (sscanf(inbuf, "%5c%5c%5c%5c%*10c%*10c%*10c%5c%5c",
               numAtomStr, resNStr, resName, atomName, segName, resIdStr) != 6) {
      printf("corplugin) Improperly formatted line: \n%s\n", inbuf);
      return -1;
    }
    strip_whitespace(resName, 8);
    strip_whitespace(atomName, 8);
    strip_whitespace(segName, 8);
  }

  numAtom = atoi(numAtomStr);
  *resId = atoi(resIdStr);
  strcpy(atomType, atomName);

  return numAtom;
}


/*
 * API functions
 */

typedef struct {
  FILE *file;
  int numatoms;
  int iofoext;      /* flag for extended CHARMM c31 version support */
} cordata;

static void *open_cor_read(const char *path, const char *filetype, 
    int *natoms) {
  int ioext ;
  FILE *fd;
  cordata *cor;

  if (!(fd = open_cor_file(path, natoms, &ioext))) {
    return NULL;
  } 
  cor = (cordata *) malloc(sizeof(cordata));
  memset(cor, 0, sizeof(cordata));
  cor->numatoms = *natoms;
  cor->file = fd;
  cor->iofoext = ioext ;
  return cor;
}

static int read_cor_structure(void *v, int *optflags, molfile_atom_t *atoms) {
  cordata *data = (cordata *)v;
  int i;
  
  /* we don't read any optional data */
  *optflags = MOLFILE_NOOPTIONS;

  for (i=0; i<data->numatoms; i++) {
    molfile_atom_t *atom = atoms+i; 
    if (get_cor_atom(data->file, atom->name, atom->type, 
                     atom->resname, atom->segid, 
                     &atom->resid, data->iofoext) < 0) {
      printf("corplugin) couldn't read atom %d\n", i);
      return MOLFILE_ERROR;
    }
    atom->chain[0] = atom->segid[0];
    atom->chain[1] = '\0';
  }

  rewind(data->file);
  return MOLFILE_SUCCESS;
}

static int read_cor_timestep(void *v, int natoms, molfile_timestep_t *ts) {
  cordata *cor = (cordata *)v;
  char inbuf[COR_RECORD_LENGTH+2], header[6];
  char xStr[21], yStr[21], zStr[21];
  int i;

  xStr[20] = '\0';
  yStr[20] = '\0';
  zStr[20] = '\0';

  /* Skip the header */
  do {
    /* Return -1 on EOF */
    if (feof(cor->file) || ferror(cor->file) || 
        (fgets(inbuf, COR_RECORD_LENGTH+1, cor->file) == NULL)) {
      return MOLFILE_ERROR;
    }

    if (sscanf(inbuf, " %5c", header) != 1) {
      printf("corplugin) Improperly formatted line.\n");
      return MOLFILE_ERROR;
    }

  } while (header[0]=='*');


  /* read the coordinates */
  for (i = 0; i < natoms; i++) {
    if (corgets(inbuf, COR_RECORD_LENGTH+1, cor->file) == NULL) {
      return MOLFILE_ERROR;
    }
    
    if (cor->iofoext == 1 ) {
      if (sscanf(inbuf, "%*10c%*10c%*10c%*10c%20c%20c%20c%*10c", 
                 xStr, yStr, zStr) != 3) {
        printf("corplugin) Error reading coordinates on line %d.\n%s\n", i, inbuf);
        return MOLFILE_ERROR;
      } else if (ts != NULL) {
        /* We have a timestep -- save the coordinates */
        ts->coords[3*i  ] = atof(xStr);
        ts->coords[3*i+1] = atof(yStr);
        ts->coords[3*i+2] = atof(zStr);
      }
    } else {
      if (sscanf(inbuf, "%*5c%*5c%*5c%*5c%10c%10c%10c%*5c", 
                 xStr, yStr, zStr) != 3) {
        printf("corplugin) Error reading coordinates on line %d.\n%s\n", i, inbuf);
        return MOLFILE_ERROR;
      } else if (ts != NULL) {
        /* We have a timestep -- save the coordinates */
        ts->coords[3*i  ] = atof(xStr);
        ts->coords[3*i+1] = atof(yStr);
        ts->coords[3*i+2] = atof(zStr);
      }
    }
  }

  return MOLFILE_SUCCESS;
}

static void close_cor_read(void *mydata) {
  cordata *data = (cordata *)mydata;
  if (data) {
    if (data->file) fclose(data->file);
    free(data);
  }
}  

/*
 * Initialization stuff down here
 */

static molfile_plugin_t plugin;

VMDPLUGIN_API int VMDPLUGIN_init() {
  memset(&plugin, 0, sizeof(molfile_plugin_t));
  plugin.abiversion = vmdplugin_ABIVERSION;
  plugin.type = MOLFILE_PLUGIN_TYPE;
  plugin.name = "cor";
  plugin.prettyname = "CHARMM Coordinates";
  plugin.author = "Eamon Caddigan, John Stone";
  plugin.majorv = 0;
  plugin.minorv = 8;
  plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
  plugin.filename_extension = "cor";
  plugin.open_file_read = open_cor_read;
  plugin.read_structure = read_cor_structure;
  plugin.read_next_timestep = read_cor_timestep;
  plugin.close_file_read = close_cor_read;
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_fini() {
  return VMDPLUGIN_SUCCESS;
}

