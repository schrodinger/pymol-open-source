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
 *      $RCSfile: basissetplugin.h,v $
 *      $Author: saam $       $Locker:  $             $State: Exp $
 *      $Revision: 1.1 $       $Date: 2009/02/20 22:37:14 $
 *
 ***************************************************************************/
/*******************************************************************
 * 
 *  headerfile for the gamessplugin
 *
 *  
 ******************************************************************/

#ifndef BASISSETPLUGIN_H
#define BASISSETPLUGIN_H

#include <stdio.h>
#include "molfile_plugin.h"


/* maximum number of Gaussian basis functions;
 * 1000 seems to be a proper upper limit for now;
 * state-of-the-art simulation could do more, hence
 * maybe increase to 5000 later */
#define MAXBASISFUNCTIONS 1000

/* define macros for true/false to make code 
 * look somewhat nicer; the macro DONE signals
 * that we're done with reading an should return
 * with what we have */
#define FALSE 0
#define TRUE  1



typedef struct {
  float exponent;
  float contraction_coeff;
} prim_t;

typedef struct {
  int numprims;
  int symmetry;     /* S, P, D, F, ...
                      * just for convenience when retrieving info */
  int wave_offset;   /* index into wave_function array */
  prim_t *prim;      /* array of primitives */
} shell_t;

/* Basis set definition for one atom */
typedef struct {
  char name[11];  /* atom name or type */
  int atomicnum;  /* atomic number (nuclear charge) */
  int numshells;
  shell_t *shell;
} basis_atom_t;


/* structure for storing temporary values read in 
 * from the gamess output file */
typedef struct 
{
  char type [11]; /* atom name or type */

  int atomicnum;  /* atomic number (nuclear charge) */

  float x,y,z;    /* coordinates of atom */
} qm_atom_t;



/* main gamess plugin data structure */
typedef struct 
{
  FILE *file;
  int numatoms;

  char basis_string[BUFSIZ]; /* basis name as "nice" string */

  char *file_name;

  /* this array of floats stores the contraction coefficients
   * and exponents for the basis functions:
   * { exp(1), c-coeff(1), exp(2), c-coeff(2), .... }
   * This holds also for double-zeta basis functions with
   * exp(i) = exp(j) and c-coeff(i) != c-coeff(j). */
  float *basis;

  basis_atom_t *basis_set;

  int num_basis_funcs;

  int num_basis_atoms;

  /* atomic number per atom in basis set */
  int *atomicnum_per_basisatom;

  /* number of shells per atom in basis set */
  int *num_shells_per_atom;

  /* the total number of atomic shells */
  int num_shells;

  /* number of primitives in shell i */
  int *num_prim_per_shell;

  /* symmetry type of each shell */
  int *shell_symmetry; 

} gamessdata;


#endif
