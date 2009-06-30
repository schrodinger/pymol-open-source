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
 *      $RCSfile: gaussianplugin.h,v $
 *      $Author: akohlmey $       $Locker:  $             $State: Exp $
 *      $Revision: 1.9 $       $Date: 2009/06/21 23:39:31 $
 *
 ***************************************************************************/
/*******************************************************************
 * 
 *  headerfile for the gaussianplugin
 *  
 ******************************************************************/

#ifndef GAUSSIANPLUGIN_H
#define GAUSSIANPLUGIN_H

#include <stdio.h>
#include "molfile_plugin.h"

/* define macros for true/false to make code 
* look somewhat nicer; the macro DONE signals
* that we're done with reading an should return
* with what we have */
#define FALSE 0
#define TRUE  1

/** macros describing the RUNTYP */
#define RUNTYP_UNKNOWN    0     /**< not set.  */
#define RUNTYP_ENERGY     1     /**< total energy, single point run  */
#define RUNTYP_OPTIMIZE   2     /**< geometry optimization  */
#define RUNTYP_SADPOINT   3     /**< transition state search  */
#define RUNTYP_HESSIAN    4     /**< frequency calculation  */
#define RUNTYP_SURFACE    5     /**< potential energy scan  */
#define RUNTYP_DYNAMICS   6     /**< molecular dynamics or monte carlo  */
#define RUNTYP_PROPERTIES 7     /**< wavefunction analysis  */

/** macros defining the SCFTYP */
#define SCFTYP_UNKNOWN 0        /**< not set. */
#define SCFTYP_RHF   1          /**< closed shell or restricted wfn.  */
#define SCFTYP_UHF   2          /**< open shell or unrestricted wfn.  */
#define SCFTYP_ROHF  3          /**< restricted open shell wfn.  */
#define SCFTYP_GVB   4          /**< generalized valence bond wfn. */
#define SCFTYP_MCSCF 5          /**< multi-configuration SCF.  */
#define SCFTYP_FF    6          /**< force field calculation.  */

/** XXX: check if this can go away. */
#define MAX_NUM_WAVE 10

typedef struct {
  float exponent;
  float contraction_coeff;
} prim_t;


typedef struct {
  int symmetry;                 /* symmetry of this shell (S, P, D, ...) */
  int numprims;                 /* number of primitive basis functions */
  prim_t *prim;                 /* array of primitives */
} shell_t;


/** Basis set definition for one atom */
typedef struct {
  char name[20];                /* atom name or type (pseudopotentialfile?) */
  int atomicnum;                /* atomic number (nuclear charge XXX) */
  int numshells;
  shell_t *shell;
} basis_atom_t;


/** structure for storing temporary values read in 
 * from the gaussian output file */
typedef struct 
{
  char type [8];                /* atom type H,N,O ..... */
  int atomicnum;                /* index in PTE. */
  float x,y,z;                  /* coordinates of atom i */
} qm_atom_t;


typedef struct {
  int   idtag;              /**< unique tag to identify this wavefunction over the trajectory */
  int   type;               /**< CANONICAL, LOCALIZED, OTHER */
  int   spin;               /**< 0 for alpha, 1 for beta */
  int   cartesian;          /**< 1 if cartesian, 0 if pure  */
  int   excitation;         /**< 0 for ground state, 1,2,3,... for excited states */
  char info[MOLFILE_BUFSIZ]; /**< string for additional type info */

  int   num_orbitals;       /**< number of orbitals that were really 
                             *   present in the output for this step */
  int   *orb_indices;       /**< list of orbital indices for wavefunction */
  float *occupancies;       /**< list of orbital occupancies for wavefunction */
  float *orb_energies;      /**< list of orbital energies for wavefunction */
  float *wave_coeffs;       /**< expansion coefficients for wavefunction in the
                             *   form {orbital1(c1),orbital1(c2),.....,orbitalM(cN)} */
} qm_wavefunction_t;


typedef struct {
  qm_wavefunction_t *wave;  /**<  */
  int numwave;              /**< number of independent wavefunctions  */
  int num_orbitals;         /**< number of orbitals per wavefunction  */
  int num_scfiter;          /**< number of SCF iterations for this step */
  double *scfenergies;      /* scfenergies per trajectory point 
                             * XXX: how about post-HF calculations?
                             *      we have the HF, MP2, CASSCF, CCS,
                             *      CCSD, CCSD(T) energies... */

  double *mulliken_charges; /* per-atom Mulliken charges */
  double *lowdin_charges;   /* per-atom Lowdin charges */
  double *esp_charges;      /* per-atom esp charges */
  double *npa_charges;      /* per-atom npa charges */

  float *gradient;          /* energy gradient for each atom */
} qm_timestep_t;


/** main gaussian plugin data structure */
typedef struct 
{
  FILE *file;
  int numatoms;
  int runtyp;        /* run type of this job as int for internal use */
  char gbasis[20];   /* canonicalized basis set string
                        for loading and external basis set. */

  char basis_string[MOLFILE_BUFSIZ]; /* basis name as "nice" string */

  char runtitle[MOLFILE_BUFSIZ];  /* title/info section of run */

  char geometry[MOLFILE_BUFSIZ];  /* either UNIQUE, CART or ZMP/ZMTMPC */
  char guess[MOLFILE_BUFSIZ];    /* type of guess method used */

  char version_string[MOLFILE_BUFSIZ]; /* full version string of binary used for run */
  int  version;  /* gaussian version code;
                  * canonicalized and set up for integer comparisons.
                  * here we keep track the exact Gaussian version,
                  * since the log file format keeps changing all
                  * the time. 
                  * Format is Year/Revision/Patchlevel: YYYYRRPP
                  * with:
                  * YYYY: g94 -> 1994, g98 -> 1998, g03 -> 2003
                  *   RR: A -> 1, B -> 2, C ->3, ...
                  *   PP: 01, 02, 03, ...
                  *
                  *   Example: G03RevB.04 -> 20030204
                  *
                  *   version = 0  => unknown/unreadable version. 
                  */

  char *file_name;

  /******************************************************
   * new API functions
   *****************************************************/

  int  scftyp;              /* UHF, RHF, ROHF, as in for 
                             * internal use*/
  char scftyp_string[MOLFILE_BUFSIZ]; /* scftyp as string. XXX: remove */
  int totalcharge;              /* Total charge of the system */
  int multiplicity;             /* Multiplicity of the system */
  int num_electrons;            /* Number of electrons */
  int nimag;                    /* Number of imaginary frequencies */
  int *nimag_modes;             /* List of imaginary modes */

  float *wavenumbers;          /* rotational and translational DoF are
                                 * included, but can be removed due to
                                 * their zero frequencies */
  float *intensities;          /* Intensities of spectral lines */

  float *normal_modes;         /* the normal modes themselves */

  int nproc;                    /* Number processors used */
  int memory;                   /* Amount of memory used in MBytes */

  int have_wavefunction; /** TRUE/FALSE flag indicating if we should
                          *  try searching for wavefunction data.
                          *  Gaussian needs IOP(6/7=3) to do this.
                          *  if it is not set, we can save time searching.
                          *  should help a lot with large log files... */

  int have_basis;        /** TRUE/FALSE flag initially indicating if we 
                          *  should try searching for basis set data.
                          *  Gaussian needs GFINPUT to print the data in a 
                          *  form that we can parse. If it is not set, we 
                          *  can save time searching. Should help a lot 
                          *  with large log files. If there is no basis in
                          *  the log, we try reading the basis set from
                          *  a local database. After initial parse it indicates
                          *  whether basis set data is available.
                          */

  int have_cart_basis;   /** flag indicating if we have a cartesian
                          *  basis set. this is additive.
                          *  0 = none (can only visualize up to P)
                          *  1 = have cartesian (6 instead of 5) d-functions)
                          *  2 = have cartesian (10 instead of 7) f-functions)
                          *  4 = have cartesian (XX instead of X) g-functions)
                          *  basis transferred to VMD is always cartesian,
                          *  so wavefunctions expanded in pure atomic orbitals
                          *  will have to be converted.
                          */

  /* arrays with atom charges */
  double *mulliken_charges; 
  int   have_mulliken; 

  int have_normal_modes; /** TRUE/FALSE flag indicating if we
			  * could properly read normal modes,
			  * wavenumbers and intensities. */

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
  
  /*******************************************************
   * end internal coordinate stuff
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
   * END OF NEW API data members
   *********************************************************/

  /* this array of floats stores the contraction coefficients
   * and exponents for the basis functions:
   * { exp(1), c-coeff(1), exp(2), c-coeff(2), .... }
   * This holds also for double-zeta basis functions with
   * exp(i) = exp(j) and c-coeff(i) != c-coeff(j). */
  float *basis;

  basis_atom_t *basis_set;

  int num_basis_funcs;          /** total number of primitive basis functions */
  int num_basis_atoms;          /** total number of atoms with basis functions */

  /** atomic number per atom in basis set */
  int *atomicnum_per_basisatom;

  /** the total number of atomic shells */
  int num_shells;
  /** number of shells per atom i */
  int *num_shells_per_atom;
  /** number of primitives in shell i */
  int *num_prim_per_shell;
  /** symmetry type of each shell */
  int *shell_symmetry; 

  /** number of occupied spin alpha and beta orbitals */
  int occ_orbitals_A;
  int occ_orbitals_B;

  /** Max. rank of the wave_function matrix.
   * While the actual number of MOs present can be different
   * for each frame, this is the maximum number of 
   * possible occupied and virtual orbitals. */
  int wavef_size;
  int num_orbitals;      /* number of orbitals in file. can be less if pure WF. */

  /** Array of length 3*num_wave_f containing the exponents 
   *  describing the cartesian components of the angular momentum. 
   *  E.g. S={0 0 0}, Px={1 0 0}, Dxy={1 1 0}, or Fyyz={0 2 1}. */
  int *angular_momentum;

  /** this flag tells if the geometry search converged */
  int converged;
  int opt_status;

  /* the structure qm_atom_t was defined to read in data from
   * the Gaussian output file and store it temporarily;
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

  /** initial cell info (in cpmd notation)":
   * a, b/a, c/a, cos(alpha), cos(beta), cos(gamma) */
  float initcell[6];

  /** per timestep data like wavefunctions and scf iterations */
  qm_timestep_t *qm_timestep;

  /** number of trajectory points; single point corresponds to 1 */
  int num_frames;
  int num_frames_sent;
  int num_frames_read;

  int end_of_trajectory;

} gaussiandata;

#endif
