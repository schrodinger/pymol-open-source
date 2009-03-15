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
 *      $Revision: 1.56 $       $Date: 2009/02/24 02:25:55 $
 *
 ***************************************************************************/

#include "molfile_plugin.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "fortread.h"

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
  int numangles, *angles;
  int numdihedrals, *dihedrals;
  int numimpropers, *impropers;
  int numcterms, *cterms;
} psfdata;



/* Read in the next atom info into the given storage areas; this assumes
   that file has already been moved to the beginning of the atom records.
   Returns the serial number of the atom. If there is an error, returns -1.*/
static int get_psf_atom(FILE *f, char *name, char *atype, char *resname,
   char *segname, int *resid, float *q, float *m, int charmmext) {
  char inbuf[PSF_RECORD_LENGTH+2];
  int num;

  if (inbuf != fgets(inbuf, PSF_RECORD_LENGTH+1, f)) {
    return(-1); /* failed to read in an atom */
  }

  if (strlen(inbuf) < 50) {
    fprintf(stderr, "Line too short in psf file: \n%s\n", inbuf);
    return -1;
  }

  num = atoi(inbuf); /* atom index */

  if (charmmext == 1) {
    /* CHARMM PSF format is (if CMAP,CHEQ are also included):
     *  '(I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,I4,1X,2G14.6,I8,2G14.6)'
     */
    strnwscpy(segname, inbuf+11, 7);
    strnwscpy(resname, inbuf+29, 7);
    strnwscpy(name, inbuf+38, 7);
    strnwscpy(atype, inbuf+47, 4);
    
    *resid = atoi(inbuf+20);
    *q = (float) atof(inbuf+52);
    *m = (float) atof(inbuf+68);
  } else {
    strnwscpy(segname, inbuf+9, 4);
    strnwscpy(resname, inbuf+19, 4);
    strnwscpy(name, inbuf+24, 4);
    strnwscpy(atype, inbuf+29, 4);

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


/*
 * Read in the beginning of the bond/angle/dihed/etc information,
 * but don't read in the data itself.  Returns the number of the record type
 * for the molecule.  If error, returns (-1). 
 */
static int psf_start_block(FILE *file, const char *blockname) {
  char inbuf[PSF_RECORD_LENGTH+2];
  int nrec = -1;

  /* keep reading the next line until a line with blockname appears */
  do {
    if(inbuf != fgets(inbuf, PSF_RECORD_LENGTH+1, file)) {
      /* EOF encountered with no blockname line found ==> error, return (-1) */
      return (-1);
    }
    if(strlen(inbuf) > 0 && strstr(inbuf, blockname))
      nrec = atoi(inbuf);
  } while (nrec == -1);

  return nrec;
}


/* Read in the bond info into the given integer arrays, one for 'from' and
   one for 'to' atoms; remember that .psf files use 1-based indices,
   not 0-based.  Returns 1 if all nbond bonds found; 0 otherwise.  */
static int psf_get_bonds(FILE *f, int nbond, int fromAtom[], int toAtom[], int charmmext) {
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


static int read_bonds(void *v, int *nbonds, int **fromptr, int **toptr, 
                      float **bondorder, int **bondtype, 
                      int *nbondtypes, char ***bondtypename) {
  psfdata *psf = (psfdata *)v;

  *nbonds = psf_start_block(psf->fp, "NBOND"); /* get bond count */

  if (*nbonds > 0) {
    psf->from = (int *) malloc(*nbonds*sizeof(int));
    psf->to = (int *) malloc(*nbonds*sizeof(int));

    if (!psf_get_bonds(psf->fp, *nbonds, psf->from, psf->to, psf->charmmext)) {
      fclose(psf->fp);
      psf->fp = NULL;
      return MOLFILE_ERROR;
    }
    *fromptr = psf->from;
    *toptr = psf->to;
    *bondorder = NULL; /* PSF files don't provide bond order or type information */
    *bondtype = NULL;
    *nbondtypes = 0;
    *bondtypename = NULL;
  } else {
    *fromptr = NULL;
    *toptr = NULL;
    *bondorder = NULL; /* PSF files don't provide bond order or type information */
    *bondtype = NULL;
    *nbondtypes = 0;
    *bondtypename = NULL;
    printf("psfplugin) WARNING: no bonds defined in PSF file.\n");
  }

  return MOLFILE_SUCCESS;
}


static int psf_get_angles(FILE *f, int n, int *angles) {
  char inbuf[PSF_RECORD_LENGTH+2];
  char *bondptr = NULL;
  int i=0;
  while (i<n) {
    if((i % 3) == 0) {
      /* must read next line */
      if(!fgets(inbuf,PSF_RECORD_LENGTH+2,f)) {
        /* early EOF encountered */
        break;
      }
      bondptr = inbuf;
    }
    if((angles[3*i] = atoi(bondptr)) < 1)
      break;
    bondptr += 8;
    if((angles[3*i+1] = atoi(bondptr)) < 1)
      break;
    bondptr += 8;
    if((angles[3*i+2] = atoi(bondptr)) < 1)
      break;
    bondptr += 8;
    i++;
  }

  return (i != n);
}


static int psf_get_dihedrals_impropers(FILE *f, int n, int *dihedrals) {
  char inbuf[PSF_RECORD_LENGTH+2];
  char *bondptr = NULL;
  int i=0;
  while (i<n) {
    if((i % 2) == 0) {
      /* must read next line */
      if(!fgets(inbuf,PSF_RECORD_LENGTH+2,f)) {
        /* early EOF encountered */
        break;
      }
      bondptr = inbuf;
    }
    if((dihedrals[4*i] = atoi(bondptr)) < 1)
      break;
    bondptr += 8;
    if((dihedrals[4*i+1] = atoi(bondptr)) < 1)
      break;
    bondptr += 8;
    if((dihedrals[4*i+2] = atoi(bondptr)) < 1)
      break;
    bondptr += 8;
    if((dihedrals[4*i+3] = atoi(bondptr)) < 1)
      break;
    bondptr += 8;
    i++;
  }

  return (i != n);
}

#if vmdplugin_ABIVERSION > 14
static int read_angles(void *v, int *numangles, int **angles, 
                       int **angletypes, int *numangletypes, 
                       char ***angletypenames, int *numdihedrals,
                       int **dihedrals, int **dihedraltypes, 
                       int *numdihedraltypes, char ***dihedraltypenames,
                       int *numimpropers, int **impropers, 
                       int **impropertypes, int *numimpropertypes, 
                       char ***impropertypenames, int *numcterms, 
                       int **cterms, int *ctermcols, int *ctermrows) {
  psfdata *psf = (psfdata *)v;

  /* initialize data to zero */
  *numangles         = 0;
  *angles            = NULL;
  *angletypes        = NULL;
  *numangletypes     = 0;
  *angletypenames    = NULL;
  *numdihedrals      = 0;
  *dihedrals         = NULL;
  *dihedraltypes     = NULL;
  *numdihedraltypes  = 0;
  *dihedraltypenames = NULL;
  *numimpropers      = 0;
  *impropers         = NULL;
  *impropertypes     = NULL;
  *numimpropertypes  = 0;
  *impropertypenames = NULL;
  *numcterms         = 0;
  *cterms            = NULL;
  *ctermrows         = 0;
  *ctermcols         = 0;

  psf->numangles    = psf_start_block(psf->fp, "NTHETA"); /* get angle count */
  if (psf->numangles > 0) {
    psf->angles = (int *) malloc(3*psf->numangles*sizeof(int));
    psf_get_angles(psf->fp, psf->numangles, psf->angles);
  } else {
    printf("psfplugin) WARNING: no angles defined in PSF file.\n");
  }
 
  psf->numdihedrals = psf_start_block(psf->fp, "NPHI");   /* get dihed count */
  if (psf->numdihedrals > 0) {
    psf->dihedrals = (int *) malloc(4*psf->numdihedrals*sizeof(int));
    psf_get_dihedrals_impropers(psf->fp, psf->numdihedrals, psf->dihedrals);
  } else {
    printf("psfplugin) WARNING: no dihedrals defined in PSF file.\n");
  }
 
  psf->numimpropers = psf_start_block(psf->fp, "NIMPHI"); /* get imprp count */
  if (psf->numimpropers > 0) {
    psf->impropers = (int *) malloc(4*psf->numimpropers*sizeof(int));
    psf_get_dihedrals_impropers(psf->fp, psf->numimpropers, psf->impropers);
  } else {
    printf("psfplugin) WARNING: no impropers defined in PSF file.\n");
  }

  psf->numcterms = psf_start_block(psf->fp, "NCRTERM"); /* get cmap count */
  if (psf->numcterms > 0) {
    psf->cterms = (int *) malloc(8*psf->numcterms*sizeof(int));

    /* same format as dihedrals, but double the number of terms */
    psf_get_dihedrals_impropers(psf->fp, psf->numcterms * 2, psf->cterms);
  } else {
    printf("psfplugin) no cross-terms defined in PSF file.\n");
  }

  *numangles = psf->numangles;
  *angles = psf->angles;

  *numdihedrals = psf->numdihedrals;
  *dihedrals = psf->dihedrals;

  *numimpropers = psf->numimpropers;
  *impropers = psf->impropers;

  *numcterms = psf->numcterms;
  *cterms = psf->cterms;

  *ctermcols = 0;
  *ctermrows = 0;

  return MOLFILE_SUCCESS;
}
#else
static int read_angles(void *v,
               int *numangles,    int **angles,    double **angleforces,
               int *numdihedrals, int **dihedrals, double **dihedralforces,
               int *numimpropers, int **impropers, double **improperforces,
               int *numcterms,    int **cterms,
               int *ctermcols,    int *ctermrows,  double **ctermforces) {
  psfdata *psf = (psfdata *)v;

  psf->numangles    = psf_start_block(psf->fp, "NTHETA"); /* get angle count */
  if (psf->numangles > 0) {
    psf->angles = (int *) malloc(3*psf->numangles*sizeof(int));
    psf_get_angles(psf->fp, psf->numangles, psf->angles);
  } else {
    printf("psfplugin) WARNING: no angles defined in PSF file.\n");
  }
 
  psf->numdihedrals = psf_start_block(psf->fp, "NPHI");   /* get dihed count */
  if (psf->numdihedrals > 0) {
    psf->dihedrals = (int *) malloc(4*psf->numdihedrals*sizeof(int));
    psf_get_dihedrals_impropers(psf->fp, psf->numdihedrals, psf->dihedrals);
  } else {
    printf("psfplugin) WARNING: no dihedrals defined in PSF file.\n");
  }
 
  psf->numimpropers = psf_start_block(psf->fp, "NIMPHI"); /* get imprp count */
  if (psf->numimpropers > 0) {
    psf->impropers = (int *) malloc(4*psf->numimpropers*sizeof(int));
    psf_get_dihedrals_impropers(psf->fp, psf->numimpropers, psf->impropers);
  } else {
    printf("psfplugin) WARNING: no impropers defined in PSF file.\n");
  }

  psf->numcterms = psf_start_block(psf->fp, "NCRTERM"); /* get cmap count */
  if (psf->numcterms > 0) {
    psf->cterms = (int *) malloc(8*psf->numcterms*sizeof(int));

    /* same format as dihedrals, but double the number of terms */
    psf_get_dihedrals_impropers(psf->fp, psf->numcterms * 2, psf->cterms);
  } else {
    printf("psfplugin) no cross-terms defined in PSF file.\n");
  }

  *numangles = psf->numangles;
  *angles = psf->angles;
  *angleforces = NULL;

  *numdihedrals = psf->numdihedrals;
  *dihedrals = psf->dihedrals;
  *dihedralforces = NULL;

  *numimpropers = psf->numimpropers;
  *impropers = psf->impropers;
  *improperforces = NULL;

  *numcterms = psf->numcterms;
  *cterms = psf->cterms;

  *ctermcols = 0;
  *ctermrows = 0;
  *ctermforces = NULL;

  return MOLFILE_SUCCESS;
}
#endif

static void close_psf_read(void *mydata) {
  psfdata *psf = (psfdata *)mydata;
  if (psf) {
    if (psf->fp != NULL) 
      fclose(psf->fp);

    /* free bond data */
    if (psf->from != NULL) 
      free(psf->from);
    if (psf->to != NULL) 
      free(psf->to);

    /* free angle data */
    if (psf->angles != NULL)
      free(psf->angles);
    if (psf->dihedrals != NULL)
      free(psf->dihedrals);
    if (psf->impropers != NULL)
      free(psf->impropers);

    /* free cross-term data */
    if (psf->cterms != NULL)
      free(psf->cterms);

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
  memset(psf, 0, sizeof(psfdata));
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
  int i, fullrows;

  if (psf->numcterms == 0) {
    fprintf(psf->fp, "PSF\n\n%8d !NTITLE\n", 1);
  } else {
    fprintf(psf->fp, "PSF CMAP\n\n%8d !NTITLE\n", 1);
  }

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

  if (psf->numangles == 0 && psf->numdihedrals == 0 && psf->numimpropers == 0 && psf->numcterms == 0) {
    printf("psfplugin) WARNING: PSF file is incomplete, no angles, dihedrals,\n");
    printf("psfplugin)          impropers, or cross-terms will be written. \n");

    fprintf(psf->fp, "%8d !NTHETA: angles\n\n\n", 0);
    fprintf(psf->fp, "%8d !NPHI: dihedrals\n\n\n", 0);
    fprintf(psf->fp, "%8d !NIMPHI: impropers\n\n\n", 0);
  } else {
    int i, numinline;

    printf("psfplugin) Writing angles/dihedrals/impropers...\n");

    fprintf(psf->fp, "%8d !NTHETA: angles\n", psf->numangles);
    for (numinline=0,i=0; i<psf->numangles; i++) {
      if ( numinline == 3 ) { fprintf(psf->fp, "\n");  numinline = 0; }
      fprintf(psf->fp, " %7d %7d %7d", 
              psf->angles[i*3], psf->angles[i*3+1], psf->angles[i*3+2]);
      numinline++;
    }
    fprintf(psf->fp, "\n\n");

    fprintf(psf->fp, "%8d !NPHI: dihedrals\n", psf->numdihedrals);
    for (numinline=0,i=0; i<psf->numdihedrals; i++) {
      if ( numinline == 2 ) { fprintf(psf->fp, "\n");  numinline = 0; }
      fprintf(psf->fp, " %7d %7d %7d %7d", 
              psf->dihedrals[i*4], psf->dihedrals[i*4+1], 
              psf->dihedrals[i*4+2], psf->dihedrals[i*4+3]);
      numinline++;
    }
    fprintf(psf->fp, "\n\n");

    fprintf(psf->fp, "%8d !NIMPHI: impropers\n", psf->numimpropers);
    for (numinline=0,i=0; i<psf->numimpropers; i++) {
      if ( numinline == 2 ) { fprintf(psf->fp, "\n");  numinline = 0; }
      fprintf(psf->fp, " %7d %7d %7d %7d",
              psf->impropers[i*4  ], psf->impropers[i*4+1], 
              psf->impropers[i*4+2], psf->impropers[i*4+3]);
      numinline++;
    }
    fprintf(psf->fp, "\n\n");
  }


  /*
   * write out empty donor/acceptor records since we don't
   * presently make use of this information 
   */
  fprintf(psf->fp, "%8d !NDON: donors\n\n\n", 0);
  fprintf(psf->fp, "%8d !NACC: acceptors\n\n\n", 0);
  fprintf(psf->fp, "%8d !NNB\n\n", 0);

  /* Pad with zeros, one for every atom */
  fullrows = psf->numatoms/8;
  for (i=0; i<fullrows; ++i)
    fprintf(psf->fp, "%8d%8d%8d%8d%8d%8d%8d%8d\n", 0, 0, 0, 0, 0, 0, 0, 0);
  for (i=psf->numatoms - fullrows*8; i; --i)
    fprintf(psf->fp, "%8d", 0);
  fprintf(psf->fp, "\n\n");
  fprintf(psf->fp, "%8d %7d !NGRP\n%8d%8d%8d\n\n", 1, 0, 0, 0, 0);


  /* write out cross-terms */
  if (psf->numcterms > 0) {
    fprintf(psf->fp, "%8d !NCRTERM: cross-terms\n", psf->numcterms);
    for (i=0; i<psf->numcterms; i++) {
      fprintf(psf->fp, " %7d %7d %7d %7d %7d %7d %7d %7d\n",
              psf->cterms[i*8  ], psf->cterms[i*8+1], 
              psf->cterms[i*8+2], psf->cterms[i*8+3],
              psf->cterms[i*8+4], psf->cterms[i*8+5],
              psf->cterms[i*8+6], psf->cterms[i*8+7]);
    }
    fprintf(psf->fp, "\n\n");
  }

  return MOLFILE_SUCCESS;
}

static int write_bonds(void *v, int nbonds, int *fromptr, int *toptr, 
                       float *bondorderptr, int *bondtype, 
                       int nbondtypes, char **bondtypename) {
  psfdata *psf = (psfdata *)v;

  /* save info until we actually write out the structure file */
  psf->nbonds = nbonds;
  psf->from = (int *) malloc(nbonds * sizeof(int));
  memcpy(psf->from, fromptr, nbonds * sizeof(int));
  psf->to = (int *) malloc(nbonds * sizeof(int));
  memcpy(psf->to, toptr, nbonds * sizeof(int));

  return MOLFILE_SUCCESS;
}

#if vmdplugin_ABIVERSION > 14
static int write_angles(void * v, int numangles, const int *angles,
                        const int *angletypes, int numangletypes,
                        const char **angletypenames, int numdihedrals, 
                        const int *dihedrals, const int *dihedraltype,
                        int numdihedraltypes, const char **dihedraltypenames,
                        int numimpropers, const int *impropers, 
                        const int *impropertypes, int numimpropertypes, 
                        const char **impropertypenames, int numcterms, 
                        const int *cterms, int ctermcols, int ctermrows) {
  psfdata *psf = (psfdata *)v;

  /* save info until we actually write out the structure file */
  psf->numangles = numangles;
  psf->numdihedrals = numdihedrals;
  psf->numimpropers = numimpropers;
  psf->numcterms = numcterms;

  psf->angles = (int *) malloc(3*psf->numangles*sizeof(int));
  memcpy(psf->angles, angles, 3*psf->numangles*sizeof(int));

  psf->dihedrals = (int *) malloc(4*psf->numdihedrals*sizeof(int));
  memcpy(psf->dihedrals, dihedrals, 4*psf->numdihedrals*sizeof(int));

  psf->impropers = (int *) malloc(4*psf->numimpropers*sizeof(int));
  memcpy(psf->impropers, impropers, 4*psf->numimpropers*sizeof(int));

  psf->cterms = (int *) malloc(8*psf->numcterms*sizeof(int));
  memcpy(psf->cterms, cterms, 8*psf->numcterms*sizeof(int));

  return MOLFILE_SUCCESS;
}
#else
static int write_angles(void * v,
        int numangles,    const int *angles,    const double *angleforces,
        int numdihedrals, const int *dihedrals, const double *dihedralforces,
        int numimpropers, const int *impropers, const double *improperforces,
        int numcterms,   const int *cterms,
        int ctermcols, int ctermrows, const double *ctermforces) {
  psfdata *psf = (psfdata *)v;

  /* save info until we actually write out the structure file */
  psf->numangles = numangles;
  psf->numdihedrals = numdihedrals;
  psf->numimpropers = numimpropers;
  psf->numcterms = numcterms;

  psf->angles = (int *) malloc(3*psf->numangles*sizeof(int));
  memcpy(psf->angles, angles, 3*psf->numangles*sizeof(int));

  psf->dihedrals = (int *) malloc(4*psf->numdihedrals*sizeof(int));
  memcpy(psf->dihedrals, dihedrals, 4*psf->numdihedrals*sizeof(int));

  psf->impropers = (int *) malloc(4*psf->numimpropers*sizeof(int));
  memcpy(psf->impropers, impropers, 4*psf->numimpropers*sizeof(int));

  psf->cterms = (int *) malloc(8*psf->numcterms*sizeof(int));
  memcpy(psf->cterms, cterms, 8*psf->numcterms*sizeof(int));

  return MOLFILE_SUCCESS;
}
#endif

static void close_psf_write(void *v) {
  psfdata *psf = (psfdata *)v;
  fclose(psf->fp);

  /* free bonds if we have them */
  if (psf->from != NULL) 
    free(psf->from);
  if (psf->to != NULL) 
    free(psf->to);

  /* free angles if we have them */
  if (psf->angles)
    free(psf->angles);
  if (psf->dihedrals)
    free(psf->dihedrals);
  if (psf->impropers)
    free(psf->impropers);

  /* free cross-terms if we have them */
  if (psf->cterms)
    free(psf->cterms);

  free(psf);
}


/*
 * Initialization stuff down here
 */

static molfile_plugin_t plugin;

VMDPLUGIN_API int VMDPLUGIN_init() {
  memset(&plugin, 0, sizeof(molfile_plugin_t));
  plugin.abiversion = vmdplugin_ABIVERSION;
  plugin.type = MOLFILE_PLUGIN_TYPE;
  plugin.name = "psf";
  plugin.prettyname = "CHARMM,NAMD,XPLOR PSF";
  plugin.author = "Justin Gullingsrud, John Stone";
  plugin.majorv = 1;
  plugin.minorv = 3;
  plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
  plugin.filename_extension = "psf";
  plugin.open_file_read = open_psf_read;
  plugin.read_structure = read_psf;
  plugin.read_bonds = read_bonds;
#if vmdplugin_ABIVERSION > 9
  plugin.read_angles = read_angles;
#endif
  plugin.close_file_read = close_psf_read;
  plugin.open_file_write = open_psf_write;
  plugin.write_structure = write_psf_structure;
  plugin.close_file_write = close_psf_write;
  plugin.write_bonds = write_bonds;
#if vmdplugin_ABIVERSION > 9
  plugin.write_angles = write_angles;
#endif
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)(void *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_fini() {
  return VMDPLUGIN_SUCCESS;
}
