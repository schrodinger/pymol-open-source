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
 *      $RCSfile: qmplugin.h,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.18 $       $Date: 2011/06/21 05:44:25 $
 *
 ***************************************************************************/
/*******************************************************************
 * 
 *  Data structures and utility functions for plugins
 *  reading logfiles from QM packages.
 *  
 ******************************************************************/

#ifndef QMPLUGIN_H
#define QMPLUGIN_H

#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include "molfile_plugin.h"


#define FALSE 0
#define TRUE  1

#define NONE  0

/* macros for shell types */
#define UNK_SHELL -666
#define SPD_SHELL   -11
#define SP_SHELL    -10
#define SPD_D_SHELL -5
#define SPD_P_SHELL -4
#define SPD_S_SHELL -3
#define SP_S_SHELL -2
#define SP_P_SHELL -1
#define S_SHELL 0
#define P_SHELL 1
#define D_SHELL 2
#define F_SHELL 3
#define G_SHELL 4
#define H_SHELL 5

#define SPIN_ALPHA  0
#define SPIN_BETA   1

/* XXX the following macros should better be in molfileplugin.h
 * since the same macros must be defined in VMD too. */


/* macros defining type of CI method (CITYP in GAMESS) */
#define CI_UNKNOWN -1
#define CI_NONE     0
#define CI_CIS      1
#define CI_ALDET    2
#define CI_ORMAS    3
#define CI_GUGA     4
#define CI_FSOCI    5
#define CI_GENCI    6

/* Basis set definition for a primitive */
typedef struct {
  float exponent;
  float contraction_coeff;
} prim_t;

/* Basis set definition for a shell */
typedef struct {
  int numprims;      /* number of primitives in this shell */
  int type;          /* S, P, D, F, ...
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


/* Atoms */
typedef struct 
{
  char type [11]; /* atom name or type */

  int atomicnum;  /* atomic number (nuclear charge) */

  float x,y,z;    /* coordinates of atom */
} qm_atom_t;


/* Wave function */
typedef struct {
  int   type;           /**< CANONICAL, LOCALIZED, OTHER */
  int   spin;           /**< 0 for alpha, 1 for beta */
  int   exci;           /**< 0 for ground state, 1,2,3,... for excited states */
  int   mult;           /**< spin multiplicity of the electronic state */
  char info[MOLFILE_BUFSIZ]; /**< string for additional type info */

  int   num_orbitals;   /**< number of orbitals that was really 
                         *   present in the output for this step */
  int   num_coeffs;     /**< number of coefficients per orbital  */
  int   has_orben;      /**< flag for orbital energies    */
  int   has_occup;      /**< flag for orbital occupancies */
  double energy;        /**< total energy for this state */
  float *wave_coeffs;   /**< expansion coefficients for wavefunction in the
                         *   form {orbital1(c1),orbital1(c2),.....,orbitalM(cN)} */
  float *orb_energies;    /**< array of orbital energies*/
  float *orb_occupancies; /**< array of orbital occupancies */
} qm_wavefunction_t;


/* Timestep specific data.
 * (Note that atoms are treated separately) */
typedef struct {
  qm_wavefunction_t *wave;
  int     numwave;          /* number of wavefunctions for this ts */
  float  *gradient;         /* energy gradient for each atom */
  int     num_scfiter;      /* number of SCF iterations */

  double *scfenergies;      /* scfenergies per trajectory point */

  /* arrays with atom charges */
  double *mulliken_charges; /* per-atom Mulliken charges */
  double *lowdin_charges;   /* per-atom Lowdin charges */
  double *esp_charges;      /* per-atom ESP charges */
  int   have_mulliken; 
  int   have_lowdin; 
  int   have_esp; 
} qm_timestep_t;


/* Main QM plugin data structure */
typedef struct 
{
  /* File format specific data.
   * This pointer must be cast to the according type by the plugin.
   * Typically, this will be a struct containing various data that
   * are needed by different functions during the file parsing process
   * and cannot be sent through the molfile_plugin QM interface since 
   * the underlying data types are specific to the file format being
   * read. */
  void *format_specific_data; 

  FILE *file;       /* the file we are reading */


  /******************************************************
   * calculation metadata (input data)
   *****************************************************/

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


  int  nproc;          /* Number processors used */
  char memory[256];    /* Amount of memory used, e.g. 1Gb */

  int totalcharge;     /* Total charge of the system */
  int multiplicity;    /* Multiplicity of the system */
  int num_electrons;   /* Number of electrons */

  char pointgroup[BUFSIZ]; /* Symmetry point group */
  int naxis;
  int order;               /* Order of highest axis */

  int mcscf_num_core;  /* Number of MCSCF core orbitals 
                        * (determines # valid orb energies
                        * for MCSCF natural orbitals) */

  int max_opt_steps;   /* Max. number of geom. opt. steps */
  float opt_tol;       /* gradient convergence tolerance,
                        * in Hartree/Bohr. */


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

  /* type of each shell */
  int *shell_types; 

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


  /******************************************************
   * normal modes
   *****************************************************/

  int have_normal_modes; /* TRUE/FALSE flag indicating if we
                          * could properly read normal modes,
                          * wavenumbers and intensities. */

  int  nimag;          /* Number of imaginary frequencies */
  int *imag_modes;     /* List of imaginary modes */

  float *wavenumbers;  /* rotational and translational DoF 
                        * are included, but can be removed due
                        * to their zero frequencies */
  float *intensities;  /* Intensities of spectral lines */

  float *normal_modes; /* the normal modes themselves */


  /******************************************************
   * internal coordinate stuff
   *****************************************************/

  int have_internals;  /* TRUE/FALSE flag indicating if we
                        * have internal coordinates */

  int have_cart_hessian; /* TRUE/FALSE flag indicating if we
                          * have a cartesian Hessian matrix */

  int have_int_hessian;  /* TRUE/FALSE flag indicating if we
                          * have a Hessian matrix in internal
                          * coordinates */

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
  
  /* XXX GAMESS */
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

  /*******************************************************
   * Trajectory related data
   *******************************************************/

  /* per timestep data like wavefunctions and scf iterations */
  qm_timestep_t *qm_timestep;

  /* array of atoms for the current timestep */
  qm_atom_t *atoms;

  /* flag to tell if SCF cycle and the geometry search converged.
   * XXX should distinguish between SCF and geometry convergence
   * in separate flags */
  int status;

  /* number of trajectory frames: */
  int num_frames;       /* total # frames */
  int num_frames_read;  /* # frames read in so far */
  int num_frames_sent;  /* # frames read sent to VMD so far */

  /* flag to indicate wether we are done with reading frames */
  int trajectory_done;

  /* file positions of the beginning of each trajectory frame */
  long *filepos_array;

  /* file position indicator for the beginning of final section 
   * printed after the last trajectory frame */
  long end_of_traj;

} qmdata_t;


/* #######################################################
 *
 * Function declarations
 *
 * ####################################################### */

/* Expand a set of symmetry unique atoms by creating images of
 * the atoms so that we have the full coordinate set. */
static int symmetry_expand(qm_atom_t **atoms, int numunique, int natoms, 
                            char *pg, int naxis);

/* Skip n lines at a time */
static void eatline(FILE * fd, int n);

/* Advance to the next non-white line */
static void eatwhitelines(FILE *fd);

/* Trim leading whitespaces from string */
static char* trimleft(char *);

/* Trim trailing whitespaces from string */
static char* trimright(char *);

/* Return 1 if the string consists of whitespace only */
static int iswhiteline(char *s);

/* Convert a string to upper case */
static char *strtoupper(char *s);


/* Place file pointer AFTER the line containing one of
 * the keystrings. */
static int pass_keyline(FILE *file, const char *keystring,
			const char *keystring2);

/* Place file pointer AT THE BEGINNING of the line containing 
 * a keystring. The keystrings are specified as a list of
 * const char* function arguments. The last argument must be 
 * NULL in order to terminate the list. */
static int goto_keyline(FILE *file, ...);

/* Check wether keystring1 occurs before keystring2 and
 * jumps back to beginning of search */
static int have_keyline(FILE *file, const char *keystring1,
                        const char *keystring2);


/* Print the current line but don't advance the file pointer. */
static void thisline(FILE *file);

/* Print next nonempty, nonwhite line but do not advance
 * the file pointer. */
static void whereami(FILE *file);



/* #######################################################
 *
 * Function definitions
 *
 * ####################################################### */


/*********************************************************
 *
 * Allocates and initiates qmdata_t structure.
 *
 *********************************************************/
static qmdata_t* init_qmdata(qmdata_t *data) {
  /* allocate memory for main data structure */
  data = (qmdata_t *)calloc(1,sizeof(qmdata_t));
  if (data == NULL) return NULL;

  data->runtype = NONE;
  data->scftype = NONE;
  data->dfttype = NONE;
  data->citype  = NONE;
  data->status = MOLFILE_QMSTATUS_UNKNOWN;
  data->trajectory_done   = FALSE;
  data->have_internals    = FALSE;
  data->have_int_hessian  = FALSE;
  data->have_cart_hessian = FALSE;
  data->have_normal_modes = FALSE;

  /* initialize some of the character arrays */
  memset(data->basis_string,0,sizeof(data->basis_string));
  memset(data->version_string,0,sizeof(data->version_string));
  memset(data->memory,0,sizeof(data->memory));

  return data;
}


/*********************************************************
 *
 * functions to manipulate the wavefunction array 
 * in qm_timestep_t.
 *
 *********************************************************/

/* Increase wavefunction array in ts by one. */
static qm_wavefunction_t* add_wavefunction(qm_timestep_t *ts) {
  if (ts->numwave) {
    /* Add a new wavefunction */
    ts->wave = (qm_wavefunction_t *)realloc(ts->wave,
                        (ts->numwave+1)*sizeof(qm_wavefunction_t));
    memset(&ts->wave[ts->numwave], 0, sizeof(qm_wavefunction_t));
    ts->numwave++;
  } else {
    /* We have no wavefunction for this timestep so create one */
    ts->wave = (qm_wavefunction_t *)calloc(1, sizeof(qm_wavefunction_t));
    ts->numwave = 1;
  }

  return &ts->wave[ts->numwave-1];
}


/* Replace the n-th wavefunction in ts with the last
 * one and decrease the array length by one. */
static void replace_wavefunction(qm_timestep_t *ts, int n) {
  if (ts->numwave>=2 && n>=0 && n<ts->numwave-1) {
    qm_wavefunction_t *w1, *w2;
    w2 = &ts->wave[n];
    w1 = &ts->wave[ts->numwave-1];
    free(w2->wave_coeffs);
    free(w2->orb_energies);
    free(w2->orb_occupancies);
    memcpy(w2, w1, sizeof(qm_wavefunction_t));
    ts->wave = (qm_wavefunction_t *) realloc(ts->wave,
                   (ts->numwave-1)*sizeof(qm_wavefunction_t));
    ts->numwave--;
  }
}


/* Delete the last wavefunction in ts */
static void del_wavefunction(qm_timestep_t *ts) {
  if (ts->numwave) {
    qm_wavefunction_t *w;
    w = &ts->wave[ts->numwave-1];
    free(w->wave_coeffs);
    free(w->orb_energies);
    free(w->orb_occupancies);
    ts->numwave--;
    ts->wave = (qm_wavefunction_t *)realloc(ts->wave,
                        ts->numwave*sizeof(qm_wavefunction_t));    
  }
}

/* Translate angular momentum string representation into
 * a triplet angular momentum exponents for X, Y, Z.
 * Angular momentum strings are a more human readable way
 * to specify the order of coefficients for wavefunctions.
 * So if the order of coefficients for D-shells implied
 * in the QM file is xx, yy, zz, xy, xz, yz then you can
 * use subsequent calls to this function in order to 
 * translate the strings into the more obscure but machine
 * friendly angular momentum exponents required by the 
 * molfile_plugin interface.
 *
 * Example translations:
 * S   --> {0 0 0}
 * X   --> {1 0 0}
 * Y   --> {0 1 0}
 * Z   --> {0 0 1}
 * XX  --> {2 0 0}
 * XY  --> {1 1 0}
 * YYZ --> {0 2 1}
 */
static void angular_momentum_expon(int  *ang_mom_expon,
                                   char *ang_mom_str) {
  int i;
  int xexp=0, yexp=0, zexp=0;

  for (i=0; i<strlen(ang_mom_str); i++) {
    switch (toupper(ang_mom_str[i])) {
    case 'X':
      xexp++;
      break;
    case 'Y':
      yexp++;
      break;
    case 'Z':
      zexp++;
      break;
    }
  }
  ang_mom_expon[0] = xexp;
  ang_mom_expon[1] = yexp;
  ang_mom_expon[2] = zexp;
}

/******************************************************
 *
 * matrix and vector functions
 *
 ******************************************************/

/* Degree-to-Radians and Radians-to-Degrees Conversion macros */
#define VMD_PI      3.14159265358979323846
#define DEGTORAD(a)     (a*VMD_PI/180.0)
#define RADTODEG(a)     (a*180.0/VMD_PI)

/* clears the matrix (resets it to identity) */
static void identity(float mat[16]) {
  memset(mat, 0, 16*sizeof(float));
  mat[0]=1.0f;
  mat[5]=1.0f;
  mat[10]=1.0f;
  mat[15]=1.0f;
}


/* Print a matrix for debugging purpose */
static void print_matrix4(const float mat[16]) {
  int i, j;
  for (i=0; i<4; i++) {
    for (j=0; j<4; j++) {
      printf("%f ", mat[4*j+i]);
    }
    printf("\n");
  }
  printf("\n");
}


/* premultiply the matrix by the given matrix */
static void multmatrix(const float *m1, float m2[16]) {
  int i, j;
  float tmp[4];

  for (j=0; j<4; j++) {
    tmp[0] = m2[j];
    tmp[1] = m2[4+j];
    tmp[2] = m2[8+j]; 
    tmp[3] = m2[12+j];
    for (i=0; i<4; i++) {
      m2[4*i+j] = m1[4*i]*tmp[0] + m1[4*i+1]*tmp[1] +
        m1[4*i+2]*tmp[2] + m1[4*i+3]*tmp[3]; 
    }
  } 
}


/* performs a rotation around an axis (char == 'x', 'y', or 'z')
 * angle is in degrees */
static void rot(float a, char axis, float mat[16]) {
  double angle;
  float m[16];
  identity(m);			/* create identity matrix */

  angle = (double)DEGTORAD(a);

  if (axis == 'x') {
    m[0]  = 1.f;
    m[5]  = (float)cos(angle);
    m[10] = m[5];
    m[6]  = (float)sin(angle);
    m[9]  = -m[6];
  } else if (axis == 'y') {
    m[0]  = (float)cos(angle);
    m[5]  = 1.f;
    m[10] = m[0];
    m[2]  = (float) -sin(angle);
    m[8]  = -m[2];
  } else if (axis == 'z') {
    m[0]  = (float)cos(angle);
    m[5]  = m[0];
    m[10] = 1.f;
    m[1]  = (float)sin(angle);
    m[4]  = -m[1];
  }

  multmatrix(m, mat);
}


/* scale a matrix */
static void scale(float s, float m[16]) {
  float t[16];
  identity(t);		/* create identity matrix */
  t[0]  = s;
  t[5]  = s;
  t[10] = s;
  multmatrix(t, m);
}


/* reflect through mirror plane */
static void mirror(char axis, float m[16]) {
  scale(-1.f, m);
  rot(180, axis, m);
}


/* multiplies a 3D point (first arg) by the Matrix, returns in second arg */
static void multpoint3d(const float *mat, const float opoint[3], float npoint[3]) {
  float tmp[3];
  float itmp3 = 1.0f / (opoint[0]*mat[3] + opoint[1]*mat[7] +
                        opoint[2]*mat[11] + mat[15]);
  tmp[0] = itmp3*opoint[0];
  tmp[1] = itmp3*opoint[1];
  tmp[2] = itmp3*opoint[2];
  npoint[0]=tmp[0]*mat[0] + tmp[1]*mat[4] + tmp[2]*mat[ 8] + itmp3*mat[12];
  npoint[1]=tmp[0]*mat[1] + tmp[1]*mat[5] + tmp[2]*mat[ 9] + itmp3*mat[13];
  npoint[2]=tmp[0]*mat[2] + tmp[1]*mat[6] + tmp[2]*mat[10] + itmp3*mat[14];
}


/* subtract 3rd vector from 2nd and put into 1st
 * in other words, a = b - c  */
static void vec_sub(float *a, const float *b, const float *c) {
  a[0]=b[0]-c[0];
  a[1]=b[1]-c[1];
  a[2]=b[2]-c[2];
}


/* length of vector */
static float norm(const float *vect) {
#if defined(_MSC_VER)
  return sqrt(vect[0]*vect[0] + vect[1]*vect[1] + vect[2]*vect[2]);
#else
  return sqrtf(vect[0]*vect[0] + vect[1]*vect[1] + vect[2]*vect[2]);
#endif
}


/******************************************************
 *
 * Symmetry
 *
 ******************************************************/

/* Expand a set of symmetry unique atoms by creating images of
 * the atoms so that we have the full coordinate set.
 * The atom array, the number of unique atoms, the total 
 * number of atoms, the pointgroup string and the order of
 * the highest axis mut be provided. The atoms array will be resized
 * and extended with the image atoms accordingly. */
static int symmetry_expand(qm_atom_t **atoms, int numunique, int natoms, 
                            char *pg, int naxis) {
  int i, j;
  float *unique, *image, *expanded;
  int *indexmap;
  int numexp;
  float m[16];

  printf("gamessplugin) Expanding %d coordinates for pointgroup %s, naxis=%d\n",
         numunique, pg, naxis);

  unique = calloc(3*numunique, sizeof(float));
  image  = calloc(3*numunique, sizeof(float));
  indexmap = calloc(numunique, sizeof(int));
  for (i=0; i<numunique; i++) {
    unique[3*i  ] = (*atoms)[i].x;
    unique[3*i+1] = (*atoms)[i].y;
    unique[3*i+2] = (*atoms)[i].z;
    indexmap[i] = i;
/*     printf("unique[%d]={%f %f %f}\n", i, unique[3*i], */
/*            unique[3*i+1],unique[3*i+2]); */
  }

  /* Define the generating symmetry operations. */
  /* XXX this is just a start, 
   *     many molecules with higher symmetries will need 
   *     two generating operations, e.g. a rotation and a
   *     reflection. */
  identity(m);
  if (!strcmp(pg, "CI")) {
    scale(-1.f, m);
  }
  else if (!strcmp(pg, "CS")) {
    mirror('y', m);
  }
  else if (!strcmp(pg, "CN")) {
    rot(360.f/naxis, 'z', m);
  }
  else if (!strcmp(pg, "CNV")) {
    rot(360.f/naxis, 'z', m);
  }
  else if (!strcmp(pg, "CNH")) {
    rot(360.f/naxis, 'z', m);
  }
  else if (!strcmp(pg, "DNH")) {
    rot(360.f/naxis, 'z', m);
  }

  for (i=0; i<numunique; i++) {
    multpoint3d(m, &unique[3*i], &image[3*i]);
/*     printf("image[%d]={%f %f %f}\n", i, image[3*i  ], */
/*            image[3*i+1], image[3*i+2]); */
  }

  expanded = calloc(3*numunique, sizeof(float));
  memcpy(expanded, unique, 3*numunique* sizeof(float));
  numexp=numunique;
  for (i=0; i<numunique; i++) {
    int found = 0;
    for (j=0; j<numunique; j++) {
      float d[3];
      vec_sub(d, &image[3*i], &unique[3*j]);
      /* printf("%d,%d norm(d)=%f\n", i, j, norm(d)); */
      if (norm(d)<0.001) {
        found = 1;
        break;
      }
    }

    if (!found) {
      expanded = realloc((float*)expanded, 3*(numexp+1)*sizeof(float));
      indexmap = realloc((int*)indexmap, (numexp+1)*sizeof(int));
      expanded[3*numexp  ] = image[3*i  ];
      expanded[3*numexp+1] = image[3*i+1];
      expanded[3*numexp+2] = image[3*i+2];
      indexmap[numexp] = i;
      numexp++;
    }
  }

/*   for (i=0; i<numexp; i++) { */
/*     printf("expanded[%d]={%f %f %f}\n", i, expanded[3*i], */
/*            expanded[3*i+1], expanded[3*i+2]); */
/*   } */

  free(unique);
  free(image);

  if (natoms && numexp!=natoms) {
    printf("gamessplugin) Couldn't expand symmetry unique atoms.\n");
    free(expanded);
    free(indexmap);
    return FALSE;
  }

  /* XXX handling the natoms==0 case is futile unless we return 
   *     the new number of atoms, so that the caller knows how
   *     many atoms there are. */
  if (!natoms)
    *atoms = calloc(numexp, sizeof(qm_atom_t));
  else 
    *atoms = realloc((qm_atom_t*)(*atoms), numexp*sizeof(qm_atom_t));

  for (i=numunique; i<numexp; i++) {
    (*atoms)[i].x = expanded[3*i  ];
    (*atoms)[i].y = expanded[3*i+1];
    (*atoms)[i].z = expanded[3*i+2];
    (*atoms)[i].atomicnum = (*atoms)[indexmap[i]].atomicnum;
    strncpy((*atoms)[i].type, (*atoms)[indexmap[i]].type, 10);
  }

  free(expanded);
  free(indexmap);
  return TRUE;
}


/******************************************************
 *
 * string functions and other parsing helpers
 *
 ******************************************************/

/* skip n lines at a time */
static void eatline(FILE * fd, int n)
{
  int i;
  for (i=0; i<n; i++) {
    char readbuf[1025];
    fgets(readbuf, 1024, fd);
  }
}


/* Skip all following consecutive white lines. */
static void eatwhitelines(FILE *file) {
  char buffer[BUFSIZ];
  long filepos;
  filepos = ftell(file);
  while (fgets(buffer, sizeof(buffer), file)) {
    if (strlen(trimright(buffer))) {
      fseek(file, filepos, SEEK_SET);
      break;
    }
    filepos = ftell(file);
  }
}


/* Returns a pointer to the first non-whitespace
 * character in a string.
 * The input string s must be null-terminated. 
 */
static char* trimleft(char* the_string)
{
  char *new_string = the_string;
  while ( (*new_string=='\n' || *new_string==' ' || *new_string=='\t') && 
	  (*new_string != '\0'))
  {
    new_string++;
  }

  return new_string;
}


/* Places a NULL character after the last non-whitespace
 * character in a string and returns the pointer to the
 * beginning of the string. 
 * The input string s must be null-terminated.
 */
static char* trimright(char* s)
{
  int i;
  for (i=strlen(s)-1; i>=0; i--) {
    if (!isspace(s[i])) break;
  }
  s[i+1] = '\0';

  return s;
}


/* Return 1 if the string contains only whitespace. */
static int iswhiteline(char *s) {
  return (!strlen(trimleft(s)));
}


/* Convert a string to upper case */
static char *strtoupper(char *s) {
  int i;
  int sz = strlen(s);

  if (s != NULL) {
    for(i=0; i<sz; i++)
      s[i] = toupper(s[i]);
  }

  return s;
}



/* Places file pointer AFTER the line containing one of the
 * keystrings. If keystring2 is NULL then it will be ignored. */
static int pass_keyline(FILE *file, const char *keystring,
        const char *keystring2) {
  char buffer[BUFSIZ];
  char *line;
  int found = 0;
  long filepos;
  filepos = ftell(file);

  do {
    if (!fgets(buffer, sizeof(buffer), file)) {
      fseek(file, filepos, SEEK_SET);
      return 0;
    }
    line = trimleft(buffer);
    if (strstr(line, keystring)) {
      found = 1;
      break;
    }
    else if (keystring2 && strstr(line, keystring2)) {
      found = 2;
      break;
    }
  } while (1);
    
  if (!found) {
    fseek(file, filepos, SEEK_SET);
  }

  return found;
}


/* Advances the file pointer until the first appearance
 * of a line containing one of the given keystrings.
 * The keystrings are specified as a list of const char*
 * function arguments. The last argument must be NULL
 * signifying the end of the keystring argument list.
 * Returns the 1-based index of the keystring from the
 * argument list for which a match was found.
 * The file pointer will be placed AT THE BEGINNING of
 * the line containing the matched keystring.
 * If no keystring was found before EOF then the file
 * is rewound to the position where the search started
 * and 0 is returned.
 */
static int goto_keyline(FILE *file, ...) {
  char buffer[BUFSIZ];
  const char *keystring;
  int found=0, loop=1;
  long filepos, curline;
  va_list argptr;
  filepos = ftell(file);

  /* loop over lines */
  while (loop) {
    int narg = 0;
    curline = ftell(file);
    if (!fgets(buffer, sizeof(buffer), file)) {
      fseek(file, filepos, SEEK_SET);

      return 0;
    }

    /* loop over the list of search strings */
    va_start(argptr, file);
    while ((keystring = va_arg(argptr, const char*)) != NULL) {
      /* search for keystring in line buffer */
      if (strstr(buffer, keystring)) {
        found = narg+1;
        /* rewind to beginning of current line */
        fseek(file, curline, SEEK_SET);
        loop = 0;
        break;
      }
      narg++;
    }
    va_end (argptr);
  };
    
  if (!found) {
    /* no match, rewind to beginning of search */
    fseek(file, filepos, SEEK_SET);
  }

  return found;
}

/* Check wether keystring1 occurs before keystring2 and
 * jumps back to beginning of search */
static int have_keyline(FILE *file, const char *keystring1,
        const char *keystring2) {
  int found;
  long filepos;
  filepos = ftell(file);

  found = pass_keyline(file, keystring1, keystring2);

  fseek(file, filepos, SEEK_SET);
  return found;
}


/* Some MOLDEN files use Fortran-style notation where 'D' is used instead of 
 * 'E' as exponential character. 
 * If you pass this function a C string that contains floating point numbers
 * in Fortran-style scientific notation, it will find the "D" characters
 * that correspond to just the floating point numbers and fix them by
 * changing them to "E".  Once converted, one should be able to
 * call scanf() on the string to parse it normally.  The code modifies
 * the string in-place, which should work well for the molden plugin. */

static int fpexpftoc(char *ftocstr) {
  int convcnt = 0;
  int len = strlen(ftocstr);

  // Compute string length, minus three chars, since any floating point 
  // number in scientific notation is going to have at least a "-" or "+",
  // and one or more integer digits following the "E" or "D" character for
  // exponent notation.
  int lenm2 = len - 2;

  // Replace "D" exponential characters with "E", so we can parse with 
  // ANSI stdio routines.
  // Our loop starts at the second character in the string since the shortest
  // possible FP number in scientific notation would be "1D+01" or something
  // similar.
  int i;
  for (i=1; i<lenm2; i++) {
    // Check for a 'D' character.  If we find one, then we look to see that the
    // preceding character is numeric, that the following character is +/-,
    // and that the character following +/- is also numeric.  If all of the
    // necessary conditions are true then we replace the 'D' with an 'E'.
    // The strict checking should allow us to safely convert entire strings
    // of characters where there are also integers, and other kinds of strings
    // on the same line.
    if (ftocstr[i] == 'D' &&
        ((ftocstr[i+1] == '-') || (ftocstr[i+1] == '+')) &&
        ((ftocstr[i-1] >= '0') && (ftocstr[i-1] <= '9')) &&
        ((ftocstr[i+2] >= '0') && (ftocstr[i+2] <= '9'))) {
      ftocstr[i] = 'E';
      convcnt++;
    }
  }

  return convcnt;
}




/*****************************************************
 *
 *   Debugging functions
 *
 *****************************************************/

/* Print the current line in the format
 * "HERE) <contents of line>".
 * Doesn't advance the file pointer.
 */
static void thisline(FILE *file) {
  char buffer[BUFSIZ];
  long filepos;
  filepos = ftell(file);
  if (!fgets(buffer, sizeof(buffer), file)) {
    if (feof(file)) printf("HERE) EOF\n");
    else printf("HERE) ????\n");
    return;
  }
  printf("HERE) %s\n", buffer);
  fseek(file, filepos, SEEK_SET);
}


/* Print all lines up to and including the next line that contains
 * non-whitespace characters. The file pointer will be restored to
 * its previous position. 
 * Output format: "HERE) <contents of line>".
 * If the end of file has been reached print "HERE) EOF".
 */
static void whereami(FILE *file) {
  char buffer[BUFSIZ];
  char *line;
  long filepos;
  filepos = ftell(file);
  do {
    if (!fgets(buffer, sizeof(buffer), file)) {
      if (feof(file)) printf("HERE) EOF\n");
      else printf("HERE) ????\n");
      return;
    }
    line = trimleft(buffer);
    printf("HERE) %s\n", buffer);
  } while (!strlen(line));

  fseek(file, filepos, SEEK_SET);
}


#endif
