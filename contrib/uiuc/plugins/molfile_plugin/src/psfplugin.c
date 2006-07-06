/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_psfplugin
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
 *      $RCSfile: psfplugin.c,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.34 $       $Date: 2006/03/30 02:53:34 $
 *
 ***************************************************************************/

#include "molfile_plugin.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define PSF_RECORD_LENGTH   160  /* extended to handle Charmm CMAP/CHEQ */ 

typedef struct {
  FILE *fp;
  int numatoms;
  int charmmfmt;  /* whether psf was written in charmm format          */
  int charmmcmap; /* stuff used by charmm for polarizable force fields */
  int charmmcheq; /* stuff used by charmm for polarizable force fields */
  int charmmext;  /* flag used by charmm for IOFOrmat EXTEnded         */
  int nbonds;
  int *from, *to;
} psfdata;



/* Read in the next atom info into the given storage areas; this assumes
   that file has already been moved to the beginning of the atom records.
   Returns the serial number of the atom. If there is an error, returns -1.*/
static int get_psf_atom(FILE *f, char *name, char *atype, char *resname,
   char *segname, int *resid, float *q, float *m, int charmmext) {
  char inbuf[PSF_RECORD_LENGTH+2];
  int i, num;

  if (inbuf != fgets(inbuf, PSF_RECORD_LENGTH+1, f)) {
    return(-1); /* failed to read in an atom */
  }

  if (strlen(inbuf) < 50) {
    fprintf(stderr, "Line too short in psf file: \n%s\n", inbuf);
    return -1;
  }

  num = atoi(inbuf); /* atom index */

  if (charmmext == 1) {
    strncpy(segname, inbuf+11, 7); 
    segname[7] = '\0';
    strncpy(resname, inbuf+29, 7);  
    resname[7] = '\0';
    strncpy(name, inbuf+38, 7);
    name[7] = '\0';
    strncpy(atype, inbuf+46, 7);
    atype[7] = '\0';
    
    /* null terminate any extraneous spaces */
    for (i=7; i >= 0; i--) {
      if (segname[i] == ' ')   
        segname[i] = '\0';
      if (resname[i] == ' ')   
        resname[i] = '\0';
      if (name[i] == ' ')      
        name[i] = '\0';
      if (atype[i] == ' ')     
        atype[i] = '\0';
    }

    *resid = atoi(inbuf+20);
    *q = (float) atof(inbuf+52);
    *m = (float) atof(inbuf+68);
  } else {
    strncpy(segname, inbuf+9, 4); 
    segname[4] = '\0';
    strncpy(resname, inbuf+19, 4);  
    resname[4] = '\0';
    strncpy(name, inbuf+24, 4);  
    name[4] = '\0';
    strncpy(atype, inbuf+29, 4); 
    atype[4] = '\0';

    /* null terminate any extraneous spaces */
    for (i=3; i >= 0; i--) {
      if (segname[i] == ' ')   
        segname[i] = '\0';
      if (resname[i] == ' ')   
        resname[i] = '\0';
      if (name[i] == ' ')      
        name[i] = '\0';
      if (atype[i] == ' ')     
        atype[i] = '\0';
    }

    *resid = atoi(inbuf+13);
    *q = (float) atof(inbuf+34);
    *m = (float) atof(inbuf+50);
  }

#if 0
  /* if this is a Charmm31 PSF file, there may be two extra */
  /* columns containing polarizable force field data.       */
  if (psf->charmmcheq) {
    /* do something to read in these columns here */
  }
#endif

  return num;
}

/* Read in the beginning of the bond information, but don't read in the
   bonds.  Returns the number of bonds in molecule.  If error, returns (-1). */
static int start_psf_bonds(FILE *f) {
  char inbuf[PSF_RECORD_LENGTH+2];
  int nbond = -1;

  /* keep reading the next line until a line with NBOND appears */
  do {
    if (inbuf != fgets(inbuf, PSF_RECORD_LENGTH+1, f)) {
      /* EOF encountered with no NBOND line found ==> error, return (-1) */
      return (-1);
    }
    if (strlen(inbuf) > 0 && strstr(inbuf,"NBOND"))
      nbond = atoi(inbuf);
  } while (nbond == -1);

  return nbond;
}


/* Read in the bond info into the given integer arrays, one for 'from' and
   one for 'to' atoms; remember that .psf files use 1-based indices,
   not 0-based.  Returns 1 if all nbond bonds found; 0 otherwise.  */
static int get_psf_bonds(FILE *f, int nbond, int fromAtom[], int toAtom[], int charmmext) {
  char *bondptr=NULL;
  char inbuf[PSF_RECORD_LENGTH+2];
  int i=0;
  size_t minlinesize;

  while (i < nbond) {
    if ((i % 4) == 0) {
      /* must read next line */
      if (!fgets(inbuf, PSF_RECORD_LENGTH+2, f)) {
        /* early EOF encountered */
        break;
      }
      /* Check that there is enough space in the line we are about to read */
      if (nbond-i >= 4) {
        if(charmmext == 1) minlinesize = 20*4 ; else minlinesize = 16*4;
      } else {
        if(charmmext == 1) minlinesize = 20*(nbond-i); else minlinesize = 16*(nbond-i);
      }
      if (strlen(inbuf) < minlinesize) {
        fprintf(stderr, "Bonds line too short in psf file: \n%s\n", inbuf);
        break;
      }
      bondptr = inbuf;
    }
    if ((fromAtom[i] = atoi(bondptr)) < 1)
      break;
    if(charmmext == 1) bondptr += 10; else bondptr += 8;
    if ((toAtom[i] = atoi(bondptr)) < 1)
      break;
    if(charmmext == 1) bondptr += 10; else bondptr += 8;
    i++;
  }

  return (i == nbond);
}

/*
 * API functions
 */

static void *open_psf_read(const char *path, const char *filetype, 
    int *natoms) {
  FILE *fp;
  char inbuf[PSF_RECORD_LENGTH*8+2];
  psfdata *psf;
  
  /* Open the .psf file and skip past the remarks to the first data section.
   * Returns the file pointer, or NULL if error.  Also puts the number of
   * atoms in the molecule into the given integer.  
   */
  if (!path)
    return NULL;

  if ((fp = fopen(path, "r")) == NULL) {
    fprintf(stderr, "Couldn't open psf file %s\n", path);
    return NULL;
  }

  *natoms = MOLFILE_NUMATOMS_NONE; /* initialize to none */

  psf = (psfdata *) malloc(sizeof(psfdata));
  memset(psf, 0, sizeof(psfdata));
  psf->fp = fp;
  psf->charmmfmt = 0; /* off unless we discover otherwise */
  psf->charmmext = 0; /* off unless we discover otherwise */

  /* read lines until a line with NATOM and without REMARKS appears    */
  do {
    /* be prepared for long lines from CNS remarks */
    if (inbuf != fgets(inbuf, PSF_RECORD_LENGTH*8+1, fp)) {
      /* EOF encountered with no NATOM line found ==> error, return null */
      *natoms = MOLFILE_NUMATOMS_NONE;
      fclose(fp);
      free(psf);
      return NULL;
    }

    if (strlen(inbuf) > 0) {
      if (!strstr(inbuf, "REMARKS")) {
        if (strstr(inbuf, "PSF")) {
          if (strstr(inbuf, "EXT")) {
            psf->charmmfmt = 1; 
            psf->charmmext = 1;      
          }
          if (strstr(inbuf, "CHEQ")) {
            psf->charmmfmt = 1; 
            psf->charmmcheq = 1;      
          }
          if (strstr(inbuf, "CMAP")) {
            psf->charmmfmt = 1; 
            psf->charmmcmap = 1;      
          }
        } else if (strstr(inbuf, "NATOM")) {
          *natoms = atoi(inbuf);
        }
      } 
    }
  } while (*natoms == MOLFILE_NUMATOMS_NONE);

  if (psf->charmmcheq || psf->charmmcmap) {
    printf("psfplugin) Detected a Charmm31 PSF file\n");
  }
  if (psf->charmmext) {
    printf("psfplugin) Detected a Charmm31 PSF EXTEnded file\n");
  }

  psf->numatoms = *natoms;

  return psf;
}

static int read_psf(void *v, int *optflags, molfile_atom_t *atoms) {
  psfdata *psf = (psfdata *)v;
  int i;
  
  /* we read in the optional mass and charge data */
  *optflags = MOLFILE_MASS | MOLFILE_CHARGE;

  for (i=0; i<psf->numatoms; i++) {
    molfile_atom_t *atom = atoms+i; 
    if (get_psf_atom(psf->fp, atom->name, atom->type, 
                     atom->resname, atom->segid, 
                     &atom->resid, &atom->charge, &atom->mass, psf->charmmext) < 0) {
      fprintf(stderr, "couldn't read atom %d\n", i);
      fclose(psf->fp);
      psf->fp = NULL;
      return MOLFILE_ERROR;
    }
    atom->chain[0] = atom->segid[0];
    atom->chain[1] = '\0';
  }

  return MOLFILE_SUCCESS;
}

static int read_bonds(void *v, int *nbonds, int **fromptr, int **toptr, float **bondorder) {
  psfdata *psf = (psfdata *)v;

  /* now read bond data */
  *nbonds = start_psf_bonds(psf->fp);

  if (*nbonds > 0) {
    psf->from = (int *) malloc(*nbonds*sizeof(int));
    psf->to = (int *) malloc(*nbonds*sizeof(int));

    if (!get_psf_bonds(psf->fp, *nbonds, psf->from, psf->to, psf->charmmext)) {
      fclose(psf->fp);
      psf->fp = NULL;
      return MOLFILE_ERROR;
    }
    *fromptr = psf->from;
    *toptr = psf->to;
    *bondorder = NULL; /* PSF files don't provide bond order information */
  } else {
    printf("psfplugin) WARNING: no bonds defined in PSF file.\n");
  }

  return MOLFILE_SUCCESS;
}


static void close_psf_read(void *mydata) {
  psfdata *psf = (psfdata *)mydata;
  if (psf) {
    if (psf->fp != NULL) 
      fclose(psf->fp);
    if (psf->from != NULL) 
      free(psf->from);
    if (psf->to != NULL) 
      free(psf->to);
    free(psf);
  }
}  


static void *open_psf_write(const char *path, const char *filetype,
    int natoms) {
  FILE *fp;
  psfdata *psf;

  fp = fopen(path, "w");
  if (!fp) {
    fprintf(stderr, "Unable to open file %s for writing\n", path);
    return NULL;
  }
  psf = (psfdata *) malloc(sizeof(psfdata));
  psf->fp = fp; 
  psf->numatoms = natoms;
  psf->charmmfmt = 0; /* initialize to off for now */
  psf->nbonds = 0;
  psf->to = NULL;
  psf->from = NULL;

  return psf;
}

static int write_psf_structure(void *v, int optflags,
    const molfile_atom_t *atoms) {
  psfdata *psf = (psfdata *)v;
  const molfile_atom_t *atom;
  int i;

  printf("psfplugin) WARNING: PSF file is incomplete, no angles, dihedrals,\n");
  printf("psfplugin)          or impropers will be written. \n");

  fprintf(psf->fp, "PSF\n\n%8d !NTITLE\n", 1);

  if (psf->charmmfmt) {
    fprintf(psf->fp," REMARKS %s\n","VMD generated structure charmm psf file");

    printf("psfplugin) WARNING: Charmm format PSF file is incomplete, atom type ID\n");
    printf("psfplugin)          codes have been emitted as '0'. \n");
  } else {
    fprintf(psf->fp," REMARKS %s\n","VMD generated structure x-plor psf file");
  }
  fprintf(psf->fp, "\n");

  /* write out total number of atoms */
  fprintf(psf->fp, "%8d !NATOM\n", psf->numatoms);

  /* write out all of the atom records */
  for (i=0; i<psf->numatoms; i++) {
    const char *atomname; 
    atom = &atoms[i];
    atomname = atom->name;

    /* skip any leading space characters given to us by VMD */ 
    while (*atomname == ' ')
      atomname++;

    if (psf->charmmfmt) {
      /* XXX replace hard-coded 0 with proper atom type ID code */
      fprintf(psf->fp, "%8d %-4s %-4d %-4s %-4s %4d %10.6f     %9.4f  %10d\n",
              i+1, atom->segid, atom->resid, atom->resname,
              atomname, /* atom->typeid */ 0, atom->charge, atom->mass, 0);
    } else {
      fprintf(psf->fp, "%8d %-4s %-4d %-4s %-4s %-4s %10.6f     %9.4f  %10d\n",
              i+1, atom->segid, atom->resid, atom->resname,
              atomname, atom->type, atom->charge, atom->mass, 0);
    }
  } 
  fprintf(psf->fp, "\n");

  /* write out bonds if we have bond information */
  if (psf->nbonds > 0 && psf->from != NULL && psf->to != NULL) {
    fprintf(psf->fp, "%8d !NBOND: bonds\n", psf->nbonds);
    for (i=0; i<psf->nbonds; i++) {
      fprintf(psf->fp, "%8d%8d", psf->from[i], psf->to[i]);
      if ((i % 4) == 3) 
        fprintf(psf->fp, "\n");
    }
    if ((i % 4) != 0) 
      fprintf(psf->fp, "\n");
    fprintf(psf->fp, "\n");
  } else {
    fprintf(psf->fp, "%8d !NBOND: bonds\n", 0);
    fprintf(psf->fp, "\n\n");
  }

  fprintf(psf->fp, "%8d !NTHETA: angles\n", 0);
  /* XXX put in code to emit angles here */
  fprintf(psf->fp, "\n\n");

  fprintf(psf->fp, "%8d !NPHI: dihedrals\n", 0);
  /* XXX put in code to emit dihedrals here */
  fprintf(psf->fp, "\n\n");

  fprintf(psf->fp, "%8d !NIMPHI: impropers\n", 0);
  /* XXX put in code to emit impropers here */
  fprintf(psf->fp, "\n\n");

  fprintf(psf->fp, "%8d !NDON: donors\n", 0);
  fprintf(psf->fp, "\n\n");

  fprintf(psf->fp, "%8d !NACC: acceptors\n", 0);
  fprintf(psf->fp, "\n\n");

  fprintf(psf->fp, "%8d !NNB\n\n", 0);
  /* Pad with zeros, one for every atom */
  {
    int i, fullrows;
    fullrows = psf->numatoms/8;
    for (i=0; i<fullrows; ++i)
      fprintf(psf->fp, "%8d%8d%8d%8d%8d%8d%8d%8d\n", 0, 0, 0, 0, 0, 0, 0, 0);
    for (i=psf->numatoms - fullrows*8; i; --i)
      fprintf(psf->fp, "%8d", 0);
  }
  fprintf(psf->fp, "\n\n");

  fprintf(psf->fp, "%8d %7d !NGRP\n%8d%8d%8d\n", 1, 0, 0, 0, 0);
  fprintf(psf->fp, "\n");

  return MOLFILE_SUCCESS;
}

static int write_bonds(void *v, int nbonds, int *fromptr, int *toptr, float *bondorderptr) {
  psfdata *psf = (psfdata *)v;

  /* save bond info until we actually write out the structure file */
  psf->nbonds = nbonds;
  psf->from = (int *) malloc(nbonds * sizeof(int));
  memcpy(psf->from, fromptr, nbonds * sizeof(int));
  psf->to = (int *) malloc(nbonds * sizeof(int));
  memcpy(psf->to, toptr, nbonds * sizeof(int));

  return MOLFILE_SUCCESS;
}

static void close_psf_write(void *v) {
  psfdata *psf = (psfdata *)v;
  fclose(psf->fp);

  /* free bonds if we have them */
  if (psf->from != NULL) 
    free(psf->from);
  if (psf->to != NULL) 
    free(psf->to);
  free(psf);
}


/*
 * Initialization stuff down here
 */

static molfile_plugin_t plugin = {
  vmdplugin_ABIVERSION,              /* ABI version              */
  MOLFILE_PLUGIN_TYPE,               /* type                     */
  "psf",                             /* short name               */
  "CHARMM,NAMD,XPLOR PSF",           /* pretty name              */
  "Justin Gullingsrud, John Stone",  /* author                   */
  0,                                 /* major version            */
  7,                                 /* minor version            */
  VMDPLUGIN_THREADSAFE,              /* is_reentrant             */
  "psf",                             /* filename extension       */
  open_psf_read,                     /* open_file_read           */
  read_psf,                          /* read_structure           */
  read_bonds,                        /* read bond list           */
  0,                                 /* read_next_timestep       */
  close_psf_read,                    /* close_file_read          */
  open_psf_write,                    /* open file for writing    */
  write_psf_structure,               /* write structure          */
  0,                                 /* write timestep           */
  close_psf_write,                   /* close file for writing   */
  0,                                 /* read volumetric metadata */
  0,                                 /* read volumetric data     */
  0,                                 /* read rawgraphics         */
  0,                                 /* read molecule metadata   */
  write_bonds                        /* write bonds              */
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
