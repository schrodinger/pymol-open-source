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
 *      $RCSfile: gamessplugin.h,v $
 *      $Author: saam $       $Locker:  $             $State: Exp $
 *      $Revision: 1.73 $       $Date: 2009/06/27 01:09:57 $
 *
 ***************************************************************************/
/*******************************************************************
 * 
 *  headerfile for the gamessplugin
 *
 *  
 ******************************************************************/

#ifndef GAMESSPLUGIN_H
#define GAMESSPLUGIN_H

#include <stdio.h>
#include "molfile_plugin.h"

/* define macros for true/false to make code 
 * look somewhat nicer; the macro DONE signals
 * that we're done with reading an should return
 * with what we have */
#define FALSE 0
#define TRUE  1

#define NONE  0

/* macros describing the RUNTYP */
#define ENERGY   1
#define OPTIMIZE 2
#define SADPOINT 3
#define HESSIAN  4
#define SURFACE  5
#define GRADIENT 6


/* macros defining the SCFTYP */
#define RHF   1
#define UHF   2
#define ROHF  3
#define GVB   4
#define MCSCF 5

/* macros defining CITYP */
#define UNKNOWN -1
#define CIS   1
#define ALDET 2
#define ORMAS 3
#define GUGA  4
#define FSOCI 5
#define GENCI 6

/* Basis set definition for a primitive */
typedef struct {
  float exponent;
  float contraction_coeff;
} prim_t;

/* Basis set definition for a shell */
typedef struct {
  int numprims;      /* number of primitives in this shell */
  int symmetry;      /* S, P, D, F, ...
                      * just for convenience when retrieving info */
  int wave_offset;   /* index into wave_function array */
  prim_t *prim;      /* array of primitives */
} shell_t;

/* Basis set definition for an atom */
typedef struct {
  char name[11];  /* atom name or type */
  int atomicnum;  /* atomic number (nuclear charge) */
  int numshells;  /* number of shells for this atom */
  shell_t *shell; /* array of shells */
} basis_atom_t;


/* structure for storing temporary values read in 
 * from the gamess output file */
typedef struct 
{
  char type [11]; /* atom name or type */

  int atomicnum;  /* atomic number (nuclear charge) */

  float x,y,z;    /* coordinates of atom */
} qm_atom_t;


typedef struct {
  int   type;           /**< CANONICAL, LOCALIZED, OTHER */
  int   spin;           /**< 0 for alpha, 1 for beta */
  int   exci;           /**< 0 for ground state, 1,2,3,... for excited states */
  char info[MOLFILE_BUFSIZ]; /**< string for additional type info */

  int   num_orbitals;   /**< number of orbitals that was really 
                         *   present in the output for this step */
  int   num_coeffs;     /**< number of coefficients per orbital */
  int   have_energies;  /**< number of orbital energies */
  int   have_occup;     /**< number of occupancies */
  float *wave_coeffs;   /**< expansion coefficients for wavefunction in the
                         *   form {orbital1(c1),orbital1(c2),.....,orbitalM(cN)} */
  float *orb_energies;  /**< list of orbital energies for wavefunction */
  float *occupancies;   /**< orbital occupancies */
} qm_wavefunction_t;


typedef struct {
  qm_wavefunction_t *wave;
  int     numwave;          /* number of wavefunctions for this ts */
  float  *gradient;         /* energy gradient for each atom */
  int     num_scfiter;      /* number of SCF iterations */

  double *scfenergies;      /* scfenergies per trajectory point */
  double *mulliken_charges; /* per-atom Mulliken charges */
  double *lowdin_charges;   /* per-atom Lowdin charges */

  double *esp_charges;      /* per-atom esp charges */
  double *npa_charges;      /* per-atom npa charges */
} qm_timestep_t;


/* main gamess plugin data structure */
typedef struct 
{
  FILE *file;       /* the file we are reading */

  int numatoms;     /* number of atoms in structure */
  int runtype;      /* type of calculation 
                     * (ENERGY, OPTIMIZE, GRADIENT, ...) */
  int scftype;      /* UHF, RHF, ROHF, ... */
  int dfttype;      /* NONE, B3LYP, ...,   */
  int citype;       /* NONE, GUGA, ...     */

  int mplevel;      /* Moller-Plesset perturbation level */

  char gbasis[10];  /* GBASIS of GAMESS run */

  char basis_string[BUFSIZ]; /* basis name as "nice" string */

  char runtitle[BUFSIZ];  /* title of gamess run */

  char geometry[BUFSIZ];  /* either UNIQUE, CART or ZMP/ZMTMPC */
  char guess[BUFSIZ];     /* type of guess method used */

  char version_string[BUFSIZ]; /* GAMESS version used for run */
  int  version; /* here we track the GAMESS versions, since the
                 * file format has changed with 
                 * version 27 JUN 2005 (R2);
                 * version = 1  : pre-27 JUN 2005 (R2)
                 * version = 2  : 27 JUN 2005 (R2)
                 * version = 0  : this we might set if we
                 *                detect an unsupported 
                 *                version and then bomb out */

  int have_pcgamess; /* this flag is set to 1 if the output
                      * file is recognized as a PC Gamess output
                      * file; we might need to introduce a few
                      * switches in the code depending on if
                      * the log file is plain Gamess or PC Gamess
                      */
  
  int  nproc;          /* Number processors used */
  char memory[256];    /* Amount of memory used, e.g. 1Gb */

  int totalcharge;     /* Total charge of the system */
  int multiplicity;    /* Multiplicity of the system */
  int num_electrons;   /* Number of electrons */



  int max_opt_steps;   /* Max. number of geom. opt. steps */
  float opt_tol;       /* gradient convergence tolerance,
                        * in Hartree/Bohr. */

  /* arrays with atom charges */
  double *mulliken_charges; 

  double *esp_charges;
  int   have_mulliken; 
  int   have_esp; 


  /******************************************************
   * normal modes
   *****************************************************/

  int have_normal_modes; /* TRUE/FALSE flag indicating if we
                          * could properly read normal modes,
                          * wavenumbers and intensities. */

  int  nimag;          /* Number of imaginary frequencies */
  int *nimag_modes;    /* List of imaginary modes */

  float *wavenumbers;  /* rotational and translational DoF 
                        * are included, but can be removed due
                        * to their zero frequencies */
  float *intensities;  /* Intensities of spectral lines */

  float *normal_modes; /* the normal modes themselves */


  /******************************************************
   * internal coordinate stuff
   *****************************************************/

  int have_internals;  /* TRUE/FALSE flag indicating if we
                        * could properly read the internal
                        * coordinates + internal hessian */

  int have_cart_hessian; /* TRUE/FALSE flag indicating if the
                          * cartesian Hessian matrix could
                          * be read from the output file */

  int nintcoords;    /* Number of internal coordinates */
  int nbonds;        /* Number of bonds */
  int nangles;       /* Number of angles */
  int ndiheds;       /* Number of dihedrals */
  int nimprops;      /* Number of impropers */

  int *bonds;        /* bond list (atom tuples) */
  int *angles;       /* angle list (atom triples) */
  int *dihedrals;    /* dihedral list (atom quadrupels) */
  int *impropers;    /* improper list (atom quadrupels) */

  double *internal_coordinates; /* value of internal coordinates */ 
  
  /* the order of force constants has to match the internal
   * coordinates in *bonds, *angles, *dihedrals */

  double *bond_force_const;     /* force constant for bonds */
  double *angle_force_const;    /* force constant for angles */
  double *dihedral_force_const; /* force constant for dihedrals */
  double *improper_force_const; /* force constant for impropers */

  /*******************************************************
   * Hessian matrices
   *******************************************************/

  double *carthessian;  /* Hessian matrix in cartesian coordinates,
                         * dimension (3*numatoms)*(3*numatoms),
                         * single array of floats 
                         * (row(1),row(2),...,row(numatoms))
                         */

  double *inthessian;  /* Hessian matrix in internal coordinates,
                        * dimension nintcoords*nintcoords,
                        * single array of floats 
                        * (row(1),row(2),...,row(nintcoords))
                        */


  /*********************************************************
   * Basis set data
   *********************************************************/

  /* this array of floats stores the contraction coefficients
   * and exponents for the basis functions:
   * { exp(1), c-coeff(1), exp(2), c-coeff(2), .... }
   * This holds also for double-zeta basis functions with
   * exp(i) = exp(j) and c-coeff(i) != c-coeff(j). */
  float *basis;

  /* hierarchical basis set structures for each atom */
  basis_atom_t *basis_set;

  /* number of uncontracted basis functions in basis array */
  int num_basis_funcs;

  /* number of atoms listed in basis set */
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

  /* number of occupied spin A and B orbitals */
  int num_occupied_A;
  int num_occupied_B;


  /* Max. size of the wave_function array per orbital.
   * I.e. this is also the number of contracted
   * cartesian gaussian basis functions or the size
   * of the secular equation.
   * While the actual # of MOs present can be different
   * for each frame, this is the maximum number of 
   * possible occupied and virtual orbitals. */
  int wavef_size;

  /* Array of length 3*num_wave_f containing the exponents 
   * describing the cartesian components of the angular momentum. 
   * E.g. S={0 0 0}, Px={1 0 0}, Dxy={1 1 0}, or Fyyz={0 2 1}. */
  int *angular_momentum;


  /* the structure qm_atom_t was defined to read in data from
   * the GAMESS output file and store it temporarily;
   * it is then copied into the VMD specific arrays at the
   * appropriate point in time;
   * this was partially implemented since the output file does
   * not, e.g., contain the number of atoms per se. One rather
   * has to count them by hand - at that point one could as 
   * well already read in the initial coordinates, atom types ...
   * which is not really supported by the way the VMD provided
   * function are arranged....this implementation could of
   * course be changed later..... */
  qm_atom_t *initatoms;

  /* per timestep data like wavefunctions and scf iterations */
  qm_timestep_t *qm_timestep;

  /* this flag tells if scf cycle and the geometry search converged */
  int opt_status;


  /* number of trajectory frames; single point corresponds to 1 */
  int num_frames;
  int num_frames_sent;
  int num_frames_read;

  /* flag to indicate wether we are done with reading frames */
  int trajectory_done;

  /* file pointers to the beginning of each trajectory frame */
  long *filepos_array;

  /* file pointer to the beginning of final section printed after
   * the last trajectory frame */
  long end_of_traj;

} gamessdata;


#endif
