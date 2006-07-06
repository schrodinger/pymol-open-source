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
 *      $RCSfile: gamessplugin.h,v $
 *      $Author: markus $       $Locker:  $             $State: Exp $
 *      $Revision: 1.32 $       $Date: 2006/03/06 16:16:50 $
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
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <time.h>
#include <math.h>
#include "molfile_plugin.h"


/* in order to be able to reserve the proper
 * amount of temporary arrays I have to define an
 * upper limit for the number of atoms in the QM
 * system; 1000 atoms should be sufficient for all
 * but abnoxiously large systems; hopefully there will
 * eventually be a more elegant way to circumvent 
 * this */
#define MAXQMATOMS 1000


/* maximum number of Gaussian basis functions;
 * 1000 seems to be a proper upper limit for now;
 * state-of-the-art simulation could do more, hence
 * maybe increase to 5000 later */
#define MAXBASISFUNCTIONS 1000


/* numerical representation of pi, from math.h */
#define MY_PI 3.1415926535897932384626433832795029L


/* convert Bohr to Angstrom */
#define BOHR_TO_ANGS 0.529177


/* convert Hartree into kcal/mol */
#define HARTREE_TO_KCAL 627.503


/* maximum number of points for the orbital
 * grid that the grid optimization code is 
 * allowed to generate */
#define MAX_GRIDPOINTS 20000


/* define macros for true/false to make code 
 * look somewhat nicer; the macro DONE signals
 * that we're done with reading an should return
 * with what we have */
#define FALSE 0
#define TRUE  1


/* macros describing the RUNTYP */
#define ENERGY   1
#define OPTIMIZE 2
#define SADPOINT 3
#define HESSIAN  4


/* macros defining the SCFTYP */
#define RHF   1
#define UHF   2
#define ROHF  3
#define GVB   4
#define MCSCF 5


/* this routine is the main gamess log file
 * parser responsible for static, i.e. 
 * non-trajectory information */
static int parse_gamess_log_static(void *, int *);


/* this routine checks if the current run is an
 * actual GAMESS run; returns true/false */
static int have_gamess(void *);


/* this routine extracts the GBASIS; returns
 * true/false */
static int get_gbasis(void *);


/* this routine reads the contrl group and
 * checks the RUNTYP terminating the plugin
 * if it encounters an unsupported one */
static int check_contrl(void *);


/* this routine prints the current date and time,
 * always useful to have that available, in my 
 * opinion at least. */
static void get_time(char *);


/* helper routine to chop spaces/newlines off
 * a C character string 
 *
 * TODO: This function is horrible and should
 *       be replaced by a cleaner solutions */
char* chop_string_all(char *);


/* helper routine to chop newlines off
 * a C character string 
 * 
 * TODO: This function is horrible and should
 *       be replaced by a cleaner solutions */
char* chop_string_nl(char *);


/* this routine renormalizes the orbital
 * coefficients read in from the input file */
float renorm_coefficient(float, float, char);


/* routine to determine the run title of the
 * GAMESS run */
static int get_runtitle(void *);


/* the function get_initial_info provides the atom number,
 * coordinates, and atom types and stores them
 * temporarily. */ 
static int get_initial_info (void *);


/* the function get_basis we also parse the basis function section to
 * determine the number of basis functions, contraction
 * coefficients. For Pople/Huzinga style basis sets
 * this numbers are in principle fixed, and could hence
 * be provided by the the plugin itself; however, the user might
 * define his own basis/contraction coeffients and hence reading
 * them from the input file seem to be somewhat more general. */
static int get_basis (void *);


/* this function reads the number of processors requested */
static int get_proc_mem(void *);


/* read in the guess options */
static int get_guess(void *);


/* this function reads the Gaussian Basis Set information
 * for an individual atom in the output file */
static int atomic_basis(int, void *, float *, char *, int *, int*);


/* this function parses the input file for the final
 * wavefunction and stores it in the appropriate arrays; */
static int get_wavefunction(void *);


/* this function parses the input file and reads the
 * number of orbitals; in the case of UHF this might
 * be two numbers, namely the number of A and B orbitals;
 * it also read the number of GAUSSIAN basis functions */
static int get_num_orbitals(void *);


/* this function is the main driver for computing and
 * handling orbital grid data
 */
static int orbital_grid_driver(void *);


/* this short test routine checks if the 
 * wavefunction/orbital stuff is supported/possible for 
 * the current GBASIS */
static int have_supported_gbasis(void *);


/* this function parses all the stored wavefunctions
 * for the HOMO which is the orbital with the smallest
 * negative energy; NOTE: selecting the HOMO in such
 * a manner might fail for some reasons, but simply is
 * currently the easiest way to implement this.
 * In the future we probably need some sort of GUI such
 * that the user can select whatever orbital they want
 * to look at not just the HOMO !! */
static int find_homo(void *);


/* this subroutine determines the cartesian origin
 * and the dimensions of the system under consideration */
static int get_system_dimensions(void *);


/* given the system dimensions as determined by the
 * get_system_dimension call it evaluates the value
 * of the orbital at the points determined by the 
 * computational grid */
static int calculate_orbital(void *);


/* this subroutine scans the output file for
 * the trajectory information */
static int get_trajectory(void *, molfile_timestep_t *, int);


/* For runtyp=HESSIAN, this subroutine scans the file for 
 * the hessian matrix in internal coordinates 
 * as well as the internal coordinate information */
static int get_int_coords(void *);


/* For runtyp=HESSIAN, this subroutine scans the file for 
 * the cartesian hessian matrix */ 
static int get_cart_hessian(void *);


/* For runtyp=HESSIAN, this subroutine reads the frequencies
 * and intensities of the normal modes */
static int get_normal_modes(void *);


/* this function calculates the value of the wavefunction
 *  * corresponding to a particular orbital at grid point
 *   * grid_x, grid_y, grid_z */
float orbital_at_grid_xyz(void*, float*, float, float, 
    float, float);


/* this function animates a given normal mode by means of
 * generating mod_num_frames frames away from the equilibrium
 * structure in a direction given by the hessiane */
static int animate_normal_mode(void*, unsigned int);


/* this function generates animated frames for normal
 * mode mode_to_animate */
static int initialize_animated_mode(void*);


/* structure for storing temporary values read in 
 * from the gamess output file */
typedef struct 
{
  char type [8]; /* atom type H,N,O ..... */

  float charge; /* array containing charge of atom i */

  float x,y,z; /* array containing the coordinate of
		* atom i*/
} gamess_temp;


/* structure for storing an animated normal mode */
typedef struct
{
  double *mode_frames; /* coordinate array containing
			    frames of an animated normal mode */

  unsigned int mode_num_frames; /* number of frames when animating 
				   modes */

  unsigned int current_mode_frame; /* tracker of current frame 
				      of an animated mode */

  double mode_scaling; /* scaling factor used for animated modes */
} mode_data;


/* main gamess plugin data structure */
typedef struct 
{
  FILE *file;
  int numatoms;
  int runtyp;   /* RUNTYP of GAMESS as int for internal use */
  char runtyp_string[BUFSIZ];  /* RUNTYP as string */  
  char gbasis[10];   /* GBASIS of GAMESS run */

  char basis_string[BUFSIZ]; /* basis name as "nice" string */
  int ngauss;        /* number of gaussian function for 
			pople style and STO basis */
  int npfunc;        /* number of p,d,f and diffuse funtions used */ 
  int ndfunc;        
  int nffunc;
  int diffs;
  int diffsp;

  char runtitle[BUFSIZ];  /* title of gamess run */

  char geometry[BUFSIZ];  /* either UNIQUE, CART or ZMP/ZMTMPC */
  char guess[BUFSIZ];    /* type of guess method used */

  char version_string[BUFSIZ]; /* GAMESS version used for run */
  int  version;  /* here we track the GAMESS versions, since the
		  * file format has changed with 
		  * version 27 JUN 2005 (R2);
		  * version = 1  : pre-27 JUN 2005 (R2)
		  * version = 2  : 27 JUN 2005 (R2)
		  * version = 0  : this we might set if we
		  *                detect an unsupported 
		  *                version and then bomb out */


  char *file_name;

  /******************************************************
   * new API functions
   *****************************************************/

  int  scftyp;              /* UHF, RHF, ROHF, as in for 
			       internal use*/
  char scftyp_string[BUFSIZ]; /* scftyp as string */
  int *atomic_number;       /* atomic numbers of molecule elements */
  /* char usertitle[80]; */ /* 80 chars for an user comment */
  /* char fromcheckfile[80]; */ /* mother checkpoint file */
  int totalcharge;          /* Total charge of the system */
  int multiplicity;         /* Multiplicity of the system */
  int num_electrons;        /* Number of electrons */
  int  nimag;           /* Number of imaginary frequencies */
  int *nimag_modes;           /* List of imaginary modes */

  int num_scfenergies;     /* number of SCF energies */
  double *scfenergies;     /* Converged SCF energies */

  double *wavenumbers; /* rotational and translational DoF 
                               are included, but can be removed due
			       to their zero frequencies */
  double *intensities; /* Intensities of spectral lines */

  double *normal_modes; /* the normal modes themselves */

  int  nproc;           /* Number processors used */
  char memory[256];     /* Amount of memory used, e.g. 1Gb */
  /* char checkfile[256]; */ /* The checkpoint file */

  /* GAUSSIAN specific */
  /* char route[1000]; */  /* line: "#..." */
  /* char geometry[256]; *//* options of the Geom keyword */
  /* char guess[256];  */  /* options of the Guess keyword */

  /* arrays with atom charges */
  double *mulliken_charges; 
  /* float *mullikengroup; */
  double *esp_charges;
  /* float *npacharges; */
  int   have_mulliken; 
  int   have_esp; 
  /* int   have_npa; */

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

  double *bond_force_const;  /* force constant for bonds */
  double *angle_force_const; /* force constant for angles */
  double *dihedral_force_const;  /* force constant for dihedrals */
  double *improper_force_const;  /* force constant for impropers */

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

  /* float** get_cartesian_hessian(units);
  float** get_internal_hessian(units);
  char* get_elements();
  char* get_names(); */

  /* NBO stuff */
  /* alpha and beta 
  int have_nbo; 
  int** lonepairs;
  int** singlebonds;
  int** doublebonds;
  int** triplebonds; */ 


  /*********************************************************
   * END OF NEW API data members
   *********************************************************/

  float *system_dimensions; /* stores the minmax xyz dimensions
                            of the system */
  float *system_center; /* stores the geometric center of the
			 * system */
  float *orbital_grid; /*this is a large array storing the 
                        * value of the orbitals at each grid
                        * point; the values are arranged as
                        * {xmin,ymin,zmin,xmin,ymin,zmin+i,
                        *  ....xmin,ymin+i,zmin,.....}
                        *  were i is distance between grid
                        *  points */
  int num_gridpoints;  /* number of gridpoints in the orbital
                        *  grid */


  /* this variable flags the presence (1) or
   * absence (!=1) of volumetric data */
  int have_volumetric;
 

  /* this struct holds the volumetric metadata
   * for the grid used to display the orbitals */
  molfile_volumetric_t *vol;


  /* this array of floats stores the contraction coefficients
   * and exponents for the basis functions:
   * { exp(1), c-coeff(1), exp(2), c-coeff(2), .... }
   * This holds also for double-zeta basis functions with
   * exp(i) = exp(j) and c-coeff(i) != c-coeff(j). 
   * The array basis_counter holds the number of basis functions
   * per atom i, in order for them to be assignable to the 
   * proper real space positions later on. 
   * the integer num_basis_funcs holds the number of elementary
   * basis functions present in the basis array */
  float *basis;
  int *basis_counter;
  int num_basis_funcs;


  /* the array atom_shells stores the number of shells per
   * atom i */
  int *atomic_shells;


  /* the array shell_primitives contains the number of 
   * primitives in shell i */
  int *shell_primitives;


  /* the array orbital_symmetry stores the symmetry type 
   * (S,L,D,..) of each (exponent, c-coeff) couple */
  char *orbital_symmetry; 


  /* the array wave_function contains the expansion coefficients
   * for the wavefunction in the form:
   * orbital1(c1,c2,c3,...),orbital2(c1,c2,c3,....) ...... */
  float *wave_function;


  /* the array orbital_energy contains the energies of 
   * all orbitals */
  float *orbital_energy;


  /* these two variable store the number of A and B orbitals */
  int num_orbitals_A;
  int num_orbitals_B;


  /* here we store the orbital index of the HOMO */
  int homo_index;

  /* this variable stores the number of GAUSSIAN basis functions */
  int num_gauss_basis_funcs;


  /* this flags signals if we were successful in reading in 
   * wavefunction information ( got_wavefunction = 1) or
   * not ( got_wavefunction = 0). This will help to avoid
   * doing wavefunction type analysis on arrays that have
   * not been properly initialized with useful data */
  int got_wavefunction;
 

  /* this variable stores the number of orbital read from the
   * output file */
  int orbital_counter;


  /* the structure gamess_temp was defined to read in data from
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
  gamess_temp *temporary;


  /* pointer to a structure that keeps track-of an animated
   * normal mode */
  mode_data *animated_mode;


  /* flag to indicate wether to read a single point or 
   * a trajectory */
  int have_trajectory;

  /* number of trajectory points; single point corresponds
   * to 1 */
  int num_traj_points;


} gamessdata;

/* this is currently a hack and provides a
 * way to access the data read by the
 * plugin from the tcl interface */
/*gamessdata* tcl_pointer; */


/* this will skip one line at a time */
static void eatline(FILE * fd)
{
  char readbuf[1025];
  fgets(readbuf, 1024, fd);
}

float calculate_wavefunction(void*, float*, float, float, 
    float, float);

#endif
