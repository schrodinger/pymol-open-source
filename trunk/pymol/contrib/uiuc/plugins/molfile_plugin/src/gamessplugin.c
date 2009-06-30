/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_gamessplugin
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
 *      $RCSfile: gamessplugin.c,v $
 *      $Author: saam $       $Locker:  $             $State: Exp $
 *      $Revision: 1.144 $       $Date: 2009/06/27 00:47:00 $
 *
 ***************************************************************************/

/* *******************************************************
 *
 *          G A M E S S     P L U G I N 
 *
 * This plugin allows VMD to read GAMESS log files
 * currently only single point geometries and trajectories
 * for optimizations, saddle point runs are supported 
 *
 * ********************************************************/

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <time.h>
#include <math.h>

#include "gamessplugin.h"
#include "unit_conversion.h"
 
#define ANGSTROM 0
#define BOHR     1
#define SPIN_ALPHA  0
#define SPIN_BETA   1

/* #define DEBUGGING 1 */

/*
 * Error reporting macro for use in DEBUG mode
 */

#define GAMESS_DEBUG
#ifdef GAMESS_DEBUG
#define PRINTERR fprintf(stderr, "\n In file %s, line %d: \n %s \n \n", \
                            __FILE__, __LINE__, strerror(errno))
#else
#define PRINTERR (void)(0)
#endif

/*
 * Error reporting macro for the multiple fgets calls in
 * the code
 */
#define GET_LINE(x,y) if (!fgets(x, sizeof(x), y)) return FALSE

#define UNK_SHELL -666
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

#define FOUND   1
#define STOPPED 2


/* ######################################################## */
/* declaration/documentation of internal (static) functions */
/* ######################################################## */

/* this routine is the main gamess log file
 * parser responsible for static, i.e. 
 * non-trajectory information */
static int parse_static_data(gamessdata *, int *);

static void print_input_data(gamessdata *);

/* this routine checks if the current run is an
 * actual GAMESS run; returns true/false */
static int have_gamess(gamessdata *);


/* this function reads the number of processors requested */
static int get_proc_mem(gamessdata *);


/* Parse the $BASIS options*/
static int get_basis_options(gamessdata *);


/* Determine the run title line */
static int get_runtitle(gamessdata *);

/* Read the input atom definitions and geometry */
static int get_input_structure(gamessdata *data);

/* Read basis set and orbital statistics such as
 * # of shells, # of A/B orbitals, # of electrons, 
 * multiplicity and total charge */
static int get_basis_stats(gamessdata *);

/* Read the contrl group and check for
 * supported RUNTYPes. Terminate the plugin
 * if an unsupported one is encountered. */
static int get_contrl(gamessdata *);

/* Read input parameters regarding calculation of 
 * certain molecular properties such as electrostatic
 * moments and the MEP. */
static int get_properties_input(gamessdata *);

/* Read symmetry point group and highest axis */
static int get_symmetry(gamessdata *);

/* read in the $GUESS options */
static int get_guess_options(gamessdata *);

/* the function get_initial_info provides the atom number,
 * coordinates, and atom types and stores them
 * temporarily. */ 
static int get_final_properties (gamessdata *);

static int get_coordinates(FILE *file, qm_atom_t **atoms, int unit,
                           int *numatoms);


/* the function get_basis we also parse the basis function section to
 * determine the number of basis functions, contraction
 * coefficients. For Pople/Huzinga style basis sets
 * this numbers are in principle fixed, and could hence
 * be provided by the the plugin itself; however, the user might
 * define his own basis/contraction coeffients and hence reading
 * them from the input file seem to be somewhat more general. */
static int get_basis (gamessdata *);


/* read all primitives for the current shell */
static int read_shell_primitives(gamessdata *, prim_t **prim,
                                 char *shellsymm, int icoeff);

/* convert shell symmetry type from char to int */
static int shellsymm_int(char symm);

/* Populate the flat arrays containing the basis set data */
static int fill_basis_arrays(gamessdata *);

static int read_first_frame(gamessdata *);

/* this subroutine scans the output file for
 * the trajectory information */
static int get_traj_frame(gamessdata *, qm_atom_t *, int);


/* returns 1 if the optimization has converged */
static int analyze_traj(gamessdata *);

/* read the number of scf iterations and the scf energies
 * for the current timestep. */
static int get_scfdata(gamessdata *, qm_timestep_t *);


/* this function parses the input file for the final
 * wavefunction and stores it in the appropriate arrays; */
static int get_wavefunction(gamessdata *, qm_timestep_t *, qm_wavefunction_t *);
static int read_coeff_block(FILE *file, int wavef_size,
                              float *wave_coeff, int *angular_momentum);

static int read_localized_orbitals(gamessdata *data);

static int get_population(gamessdata *, qm_timestep_t *);

static int get_gradient(gamessdata *, qm_timestep_t *ts);

/* Read ESP charges. */
static int get_esp_charges(gamessdata *data);

/* For runtyp=HESSIAN, this subroutine scans the file for 
 * the hessian matrix in internal coordinates 
 * as well as the internal coordinate information */
static int get_int_coords(gamessdata *);


/* For runtyp=HESSIAN, this subroutine scans the file for 
 * the cartesian hessian matrix */ 
static int get_cart_hessian(gamessdata *);


/* For runtyp=HESSIAN, this subroutine reads the frequencies
 * and intensities of the normal modes */
static int get_normal_modes(gamessdata *);


/* helper routine to chop spaces/newlines off
 * a C character string 
 *
 * TODO: This function is horrible and should
 *       be replaced by a cleaner solution */
static char* chop_string_all(char *);

static char* trimleft(char *);
static char* trimright(char *);
static void eatwhitelines(FILE *fd);

static int goto_keyline(FILE *file, const char *keystring,
                        const char *stopstring);
static int pass_keyline(FILE *file, const char *keystring,
                          const char *stopstring);
static int goto_keystring2(FILE *file, const char *keystring,
        const char *stopstring1, const char *stopstring2);
static void whereami(FILE *file);

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

/* helper routine to chop newlines off
 * a C character string 
 * 
 * TODO: This function is horrible and should
 *       be replaced by a cleaner solution */
static char* chop_string_nl(char *);


/* skip n line at a time */
static void eatline(FILE * fd, int n)
{
  int i;
  for (i=0; i<n; i++) {
    char readbuf[1025];
    fgets(readbuf, 1024, fd);
  }
}


/* Increase wavefunction array in ts by one. */
static qm_wavefunction_t* add_wavefunction(qm_timestep_t *ts) {
  if (ts->numwave) {
    /* Add a new wavefunction */
    ts->wave = (qm_wavefunction_t *)realloc(ts->wave,
                        (ts->numwave+1)*sizeof(qm_wavefunction_t));
    ts->numwave++;
  } else {
    /* We have no wavefunction for this timestep so create one */
    ts->wave = (qm_wavefunction_t *)calloc(1, sizeof(qm_wavefunction_t));
    ts->numwave = 1;
  }
  memset(&ts->wave[ts->numwave-1], 0, sizeof(qm_wavefunction_t));

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
    free(w2->occupancies);
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
    free(w->occupancies);
    ts->numwave--;
    ts->wave = (qm_wavefunction_t *)realloc(ts->wave,
                        ts->numwave*sizeof(qm_wavefunction_t));    
  }
}


/* ######################################################## */
/* Functions that are needed by the molfile_plugin          */
/* interface to provide VMD with the parsed data            */
/* ######################################################## */


/***************************************************************
 *
 * Called by VMD to open the GAMESS logfile and get the number
 * of atoms.
 * We are also reading all the static (i.e. non-trajectory)
 * data here since we have to parse a bit to get the atom count
 * anyway. These data will then be provided to VMD by
 * read_gamess_metadata() and read_gamess_rundata().
 *
 * *************************************************************/
static void *open_gamess_read(const char *filename, 
                  const char *filetype, int *natoms) {

  FILE *fd;
  gamessdata *data;
  
  /* open the input file */
  fd = fopen(filename, "rb");
 
  if (!fd) {
    PRINTERR;
    return NULL;
  }

  /* allocate memory for main data structure */
  data = (gamessdata *)calloc(1,sizeof(gamessdata));

  /* make sure memory was allocated properly */
  if (data == NULL) {
    PRINTERR;
    return NULL;
  }

  data->runtype = NONE;
  data->scftype = NONE;
  data->dfttype = NONE;
  data->citype  = NONE;
  data->num_shells = 0;
  data->num_basis_funcs = 0;
  data->num_basis_atoms = 0;
  data->num_frames = 0;
  data->num_frames_sent = 0;
  data->num_frames_read = 0;
  data->trajectory_done = FALSE;
  data->opt_status = MOLFILE_QM_STATUS_UNKNOWN;
  data->have_internals = FALSE;
  data->have_cart_hessian = FALSE;
  data->have_normal_modes = FALSE;
  data->nimag = 0;
  data->num_electrons = 0;

  data->nintcoords = 0;
  data->nbonds = 0;
  data->nangles = 0;
  data->ndiheds = 0;
  data->nimprops = 0;

  data->version = 0;
  data->have_pcgamess = 0;

  /* initialize some of the character arrays */
  memset(data->basis_string,0,sizeof(data->basis_string));
  memset(data->version_string,0,sizeof(data->version_string));
  memset(data->memory,0,sizeof(data->memory));

  /* store file pointer and filename in gamess struct */
  data->file = fd;


  /* check if the file is GAMESS format; if yes
   * parse it, if no exit */
  if (have_gamess(data)==TRUE) {
    /* if we're dealing with an unsupported GAMESS
     * version, we better quit */
    if (data->version==0) {
      printf("gamessplugin) GAMESS version %s not supported. \n",
             data->version_string);
      printf("gamessplugin) .... bombing out! Sorry :( \n");
      return NULL;
    }

    /* get the non-trajectory information from the log file */    
    if (parse_static_data(data, natoms) == FALSE) 
      return NULL;
  }
  else {
    return NULL;
  }

  return data;
}


/************************************************************
 * 
 * Provide VMD with the structure of the molecule, i.e the
 * atoms coordinates names, etc.
 *
 *************************************************************/
static int read_gamess_structure(void *mydata, int *optflags, 
                      molfile_atom_t *atoms) 
{
  gamessdata *data = (gamessdata *)mydata;
  qm_atom_t *cur_atom;
  molfile_atom_t *atom;
  int i = 0;
 
  *optflags = MOLFILE_ATOMICNUMBER; /* no optional data */
  if (data->have_mulliken) 
    *optflags |= MOLFILE_CHARGE;

  /* all the information I need has already been read in
   * via the initial scan and I simply need to copy 
   * everything from the temporary arrays into the 
   * proper VMD arrays.
   * Since there are no atom names in the GAMESS output
   * I use the atom type here --- maybe there is a better
   * way to do this !!?? */

  /* get initial pointer for atom array */
  cur_atom = data->initatoms;

  for(i=0; i<data->numatoms; i++) {
    atom = atoms+i;
    strncpy(atom->name, cur_atom->type, sizeof(atom->name)); 
    strncpy(atom->type, cur_atom->type, sizeof(atom->type));
    strncpy(atom->resname,"", sizeof(atom->resname)); 
    atom->resid = 1;
    atom->chain[0] = '\0';
    atom->segid[0] = '\0';
    atom->atomicnumber = cur_atom->atomicnum;
#ifdef DEBUGGING
    printf("gamessplugin) atomicnum[%d] = %d\n", i, atom->atomicnumber);
#endif

    if (data->have_mulliken)
      atom->charge = data->qm_timestep->mulliken_charges[i];

    cur_atom++;
  }
 
  return MOLFILE_SUCCESS; 
}


/*****************************************************
 *
 * provide VMD with the sizes of the QM related
 * data structure arrays that need to be made
 * available
 *
 *****************************************************/
static int read_gamess_metadata(void *mydata, 
    molfile_qm_metadata_t *metadata) {

  gamessdata *data = (gamessdata *)mydata;

  if (data->runtype == HESSIAN) {
    metadata->ncart = (3*data->numatoms);
    metadata->nimag = data->nimag;             
   
    if (data->have_internals) {
      metadata->nintcoords = data->nintcoords; 
    } else {
      metadata->nintcoords = 0;
    }
  }
  else {
    metadata->ncart = 0;
    metadata->nimag = 0;
    metadata->nintcoords = 0;
  }

  /* orbital data */
  metadata->num_basis_funcs = data->num_basis_funcs;
  metadata->num_basis_atoms = data->num_basis_atoms;
  metadata->num_shells      = data->num_shells;
  metadata->wavef_size      = data->wavef_size;  

#if vmdplugin_ABIVERSION > 11
  /* system and run info */
  metadata->have_sysinfo = 1;

  /* charges */
  metadata->have_esp = data->have_esp;

  /* hessian info */
  metadata->have_carthessian = data->have_cart_hessian;
  metadata->have_internals   = data->have_internals;

  /* normal mode info */
  metadata->have_normalmodes = data->have_normal_modes;
#endif

  return MOLFILE_SUCCESS;
}


/******************************************************
 * 
 * Provide VMD with the static (i.e. non-trajectory)
 * data. That means we are filling the molfile_plugin
 * data structures.
 *
 ******************************************************/
static int read_gamess_rundata(void *mydata, 
                               molfile_qm_t *qm_data) {

  gamessdata *data = (gamessdata *)mydata;
  int i, j;
  int ncart;
  molfile_qm_hessian_t *hessian_data = &qm_data->hess;
  molfile_qm_basis_t   *basis_data   = &qm_data->basis;
  molfile_qm_sysinfo_t *sys_data     = &qm_data->run;

  /* fill in molfile_qm_hessian_t */
  if (data->runtype == HESSIAN) {
    ncart = 3*data->numatoms;

    /* Hessian matrix in cartesian coordinates */
    if (data->have_cart_hessian) {
      for (i=0; i<ncart; i++) {
        for (j=0; j<=i; j++) {
          hessian_data->carthessian[ncart*i+j] = data->carthessian[ncart*i+j];
          hessian_data->carthessian[ncart*j+i] = data->carthessian[ncart*i+j];
        }
      }
    }

    /* Hessian matrix in internal coordinates */
    if (data->have_internals) {
      for (i=0; i<(data->nintcoords)*(data->nintcoords); i++) {
        hessian_data->inthessian[i] = data->inthessian[i];
      }
    }

    /* wavenumbers, intensities, normal modes */
    if (data->have_normal_modes) {
      for (i=0; i<ncart*ncart; i++) {
        hessian_data->normalmodes[i] = data->normal_modes[i];
      }
      for (i=0; i<ncart; i++) {
        hessian_data->wavenumbers[i] = data->wavenumbers[i];
        hessian_data->intensities[i] = data->intensities[i];
      }
    }

    /* imaginary modes */
    for (i=0; i<data->nimag; i++) {
      hessian_data->imag_modes[i] = data->nimag_modes[i];
    }
  }

  /* fill in molfile_qm_sysinfo_t */
  sys_data->runtype = data->runtype;
  sys_data->scftype = data->scftype;
  sys_data->nproc = data->nproc;
  sys_data->num_electrons = data->num_electrons;
  sys_data->totalcharge = data->totalcharge;
  sys_data->multiplicity = data->multiplicity;
  sys_data->num_occupied_A = data->num_occupied_A;
  sys_data->num_occupied_B = data->num_occupied_B;

#if vmdplugin_ABIVERSION > 11
  sys_data->status = data->opt_status;

  if (data->have_esp) {
    for (i=0; i<data->numatoms; i++) {
      sys_data->esp_charges[i] = data->esp_charges[i];
    }
  }
#endif

  strncpy(sys_data->basis_string, data->basis_string,
          sizeof(sys_data->basis_string));

  sys_data->memory = 0; /* XXX fixme */

  strncpy(sys_data->runtitle, data->runtitle, sizeof(sys_data->runtitle));
  strncpy(sys_data->geometry, data->geometry, sizeof(sys_data->geometry));
  strncpy(sys_data->version_string, data->version_string,
          sizeof(sys_data->version_string));

#if vmdplugin_ABIVERSION > 11
  /* fill in molfile_qm_basis_t */
  if (data->num_basis_funcs) {
    for (i=0; i<data->num_basis_atoms; i++) {
      basis_data->num_shells_per_atom[i] = data->num_shells_per_atom[i];
      basis_data->atomic_number[i] = data->atomicnum_per_basisatom[i];
    }
    
    for (i=0; i<data->num_shells; i++) {
      basis_data->num_prim_per_shell[i] = data->num_prim_per_shell[i];
      basis_data->shell_symmetry[i] = data->shell_symmetry[i];
    }
    
    for (i=0; i<2*data->num_basis_funcs; i++) {
      basis_data->basis[i] = data->basis[i];
    }

    for (i=0; i<3*data->wavef_size; i++) {
      basis_data->angular_momentum[i] = data->angular_momentum[i];
    }
  }
#endif
 
  return MOLFILE_SUCCESS;
}




#if vmdplugin_ABIVERSION > 11

/***********************************************************
 * Provide non-QM metadata for next timestep. 
 * Required by the plugin interface.
 ***********************************************************/
static int read_timestep_metadata(void *mydata,
                                  molfile_timestep_metadata_t *meta) {
  meta->count = -1;
  meta->has_velocities = 0;

  return MOLFILE_SUCCESS;
}

/***********************************************************
 * Provide QM metadata for next timestep. 
 * This actually triggers reading the entire next timestep
 * since we have to parse the whole timestep anyway in order
 * to get the metadata. So we store the read data locally
 * and hand them to VMD when requested by read_timestep().
 *
 ***********************************************************/
static int read_qm_timestep_metadata(void *mydata,
                                    molfile_qm_timestep_metadata_t *meta) {
  int have = 0;

  gamessdata *data = (gamessdata *)mydata;
  printf("gamessplugin) read_qm_timestep_metadata(): %i/%i/%i\n",
         data->num_frames, 
         data->num_frames_read,
         data->num_frames_sent);

  meta->count = -1; /* Don't know the number of frames yet */
  meta->has_gradient = 0;

  if (data->num_frames_read > data->num_frames_sent) {
    have = 1;
  }
  else if (data->num_frames_read < data->num_frames) {
    printf("gamessplugin) Probing timestep %i\n", data->num_frames_read);

    have = get_traj_frame(data, data->initatoms, data->numatoms);
  }

  if (have) {
    int i;
    qm_timestep_t *cur_ts;

    /* get a pointer to the current qm timestep */
    cur_ts = data->qm_timestep+data->num_frames_sent;
    printf("gamessplugin) Approved timestep %i\n", data->num_frames_sent);
    
    meta->num_scfiter  = cur_ts->num_scfiter;

    for (i=0; (i<MOLFILE_MAXWAVEPERTS && i<cur_ts->numwave); i++) {
      char typestr[MOLFILE_BUFSIZ];
      strcpy(typestr, cur_ts->wave[i].info);
      
      printf("gamessplugin) num_orbitals_per_wavef[%d/%d]=%3d type=%s\n",
             i+1, cur_ts->numwave, cur_ts->wave[i].num_orbitals, typestr);
      meta->num_orbitals_per_wavef[i] = cur_ts->wave[i].num_orbitals;
    }
    meta->num_wavef  = cur_ts->numwave;
    meta->wavef_size = data->wavef_size;

  } else {
    meta->num_scfiter  = 0;
    meta->num_orbitals_per_wavef[0] = 0;
    meta->num_wavef = 0;
    meta->wavef_size = 0;

    data->trajectory_done = TRUE;
  }

  return MOLFILE_SUCCESS;
}


/***********************************************************
 *
 * This function provides the data of the next timestep.
 * Here we actually don't read the data from file, that had
 * to be done already upon calling read_timestep_metadata().
 * Instead we copy the stuff from the local data structure
 * into the one's provided by VMD.
 *
 ***********************************************************/
static int read_timestep(void *mydata, int natoms, 
       molfile_timestep_t *ts, molfile_qm_metadata_t *qm_metadata,
			 molfile_qm_timestep_t *qm_ts) 
{
  gamessdata *data = (gamessdata *)mydata;
  qm_atom_t *cur_atom;
  int i = 0;
  qm_timestep_t *cur_ts;

  if (data->trajectory_done == TRUE) return MOLFILE_ERROR;

  printf("gamessplugin) Sending timestep %i\n", data->num_frames_sent);

  /* initialize pointer for temporary arrays */
  cur_atom = data->initatoms; 
  
  /* copy the coordinates */
  for(i=0; i<natoms; i++) {
    ts->coords[3*i  ] = cur_atom->x;
    ts->coords[3*i+1] = cur_atom->y;
    ts->coords[3*i+2] = cur_atom->z; 
    cur_atom++;
  }    
  
  /* get a convenient pointer to the current qm timestep */
  cur_ts = data->qm_timestep+data->num_frames_sent;

  /* store the SCF energies */
  for (i=0; i<cur_ts->num_scfiter; i++) {
    qm_ts->scfenergies[i] = cur_ts->scfenergies[i];
  }

  /* store the wave function and orbital energies */
  if (cur_ts->wave) {
    for (i=0; i<cur_ts->numwave; i++) {
      qm_wavefunction_t *wave = &cur_ts->wave[i];
      qm_ts->wave[i].type        = wave->type;
      qm_ts->wave[i].spin        = wave->spin;
      qm_ts->wave[i].excitation  = wave->exci;
      strncpy(qm_ts->wave[i].info, wave->info, MOLFILE_BUFSIZ);

      if (wave->wave_coeffs && wave->orb_energies) {
        memcpy(qm_ts->wave[i].wave_coeffs, wave->wave_coeffs,
               wave->num_orbitals*data->wavef_size*sizeof(float));
        memcpy(qm_ts->wave[i].orbital_energies, wave->orb_energies,
               wave->num_orbitals*sizeof(float));
      }
    }
  }

  if (data->runtype == ENERGY || data->runtype == HESSIAN) {
    /* We have only a single point */
    data->trajectory_done = TRUE;
  }

  data->num_frames_sent++;

  return MOLFILE_SUCCESS;
}
#endif



/**********************************************************
 *
 * clean up when done and free all the memory do avoid
 * memory leaks
 *
 **********************************************************/
static void close_gamess_read(void *mydata) {

  gamessdata *data = (gamessdata *)mydata;
  int i, j;
  fclose(data->file);

  free(data->initatoms);
  free(data->basis);
  free(data->shell_symmetry);
  free(data->atomicnum_per_basisatom);
  free(data->num_shells_per_atom);
  free(data->num_prim_per_shell);
  free(data->mulliken_charges);
  free(data->esp_charges);
  free(data->bonds);
  free(data->angles);
  free(data->dihedrals);
  free(data->impropers);
  free(data->internal_coordinates);
  free(data->bond_force_const);
  free(data->angle_force_const);
  free(data->dihedral_force_const);
  free(data->improper_force_const);
  free(data->inthessian);
  free(data->carthessian);
  free(data->wavenumbers);
  free(data->intensities);
  free(data->normal_modes);
  free(data->angular_momentum);
  free(data->filepos_array);

  if (data->basis_set) {
    for(i=0; i<data->num_basis_atoms; i++) {
      for (j=0; j<data->basis_set[i].numshells; j++) {
        free(data->basis_set[i].shell[j].prim);
      }
      free(data->basis_set[i].shell);
    } 
    free(data->basis_set);
  }

  for (i=0; i<data->num_frames; i++) {
    free(data->qm_timestep[i].scfenergies);
    free(data->qm_timestep[i].gradient);
    free(data->qm_timestep[i].mulliken_charges);
    for (j=0; j<data->qm_timestep[i].numwave; j++) {
      free(data->qm_timestep[i].wave[j].wave_coeffs);
      free(data->qm_timestep[i].wave[j].orb_energies);
      free(data->qm_timestep[i].wave[j].occupancies);
    }
    free(data->qm_timestep[i].wave);
  }
  free(data->qm_timestep);
  
  free(data);
}

/* ####################################################### */
/*             End of API functions                        */
/* The following functions actually do the file parsing.   */
/* ####################################################### */



/********************************************************
 *
 * Main gamess log file parser responsible for static,  
 * i.e. non-trajectory information.
 *
 ********************************************************/
static int parse_static_data(gamessdata *data, int *natoms) 
{
  /* Read # of procs and amount of requested memory */
  if (!get_proc_mem(data))        return FALSE;

  /* Read the $BASIS options */
  if (!get_basis_options(data))   return FALSE;

  /* Read the run title */
  if (!get_runtitle(data))        return FALSE;

  /* Read the input atom definitions and geometry */
  if (!get_input_structure(data)) return FALSE;

  /* Read the basis set */
  if (!get_basis(data))           return FALSE; 

  /* Read the number of orbitals, electrons, 
   * charge, multiplicity, ... */
  if (!get_basis_stats(data))     return FALSE;

  /* Read the $CONTRL group; here we determine the
   * RUNTYP and bomb out if we don't support
   * it; also read in some control flow related
   * info as well as the COORD type */
  if (!get_contrl(data))          return FALSE;

  /* Read input parameters regarding calculation of 
   * certain molecular properties such as electrostatic
   * moments and the MEP. */
  if (!get_properties_input(data)) return FALSE;

  /* Read symmetry point group and highest axis */
  if (!get_symmetry(data))         return FALSE;

  /* Read the $GUESS options */
  if (!get_guess_options(data))    return FALSE;

  /* Find the end of the trajectory and count the
   * frames on the way.
   * If no regular end is found we won't find any
   * properties to read and return. */
  if (!analyze_traj(data)) {
    printf("gamessplugin) WARNING: Truncated or abnormally terminated file!\n\n");
  }



  /* provide VMD with the proper number of atoms */
  *natoms = data->numatoms;

  /* Read the first frame*/
  read_first_frame(data);


  /* Read the properties at the end of a calculation */
  get_final_properties(data);

  printf("gamessplugin) num_frames_read = %i\n", data->num_frames_read);
  printf("gamessplugin) num_frames_sent = %i\n", data->num_frames_sent);

#ifdef DEBUGGING 
  /* Test print the parsed data in same format as logfile */
  print_input_data(data);
#endif

  return TRUE;
}


#define TORF(x) (x ? "T" : "F")

static void print_input_data(gamessdata *data) {
  int i, j, k;
  int primcount=0;
  int shellcount=0;

  printf("\nDATA READ FROM FILE:\n\n");
  printf(" %10s WORDS OF MEMORY AVAILABLE\n", data->memory);
  printf("\n");
  printf("     BASIS OPTIONS\n");
  printf("     -------------\n");
  printf("%s\n", data->basis_string);
  printf("\n\n\n");
  printf("     RUN TITLE\n");
  printf("     ---------\n");
  printf(" %s\n", data->runtitle);
  printf("\n");
  printf(" THE POINT GROUP OF THE MOLECULE IS %s\n", "XXX");
  printf(" THE ORDER OF THE PRINCIPAL AXIS IS %5i\n", 0);
  printf("\n");
  printf(" YOUR FULLY SUBSTITUTED Z-MATRIX IS\n");
  printf("\n");
  printf(" THE MOMENTS OF INERTIA ARE (AMU-ANGSTROM**2)\n");
  printf(" IXX=%10.3f   IYY=%10.3f   IZZ=%10.3f\n", 0.0, 0.0, 0.0);
  printf("\n");
  printf(" ATOM      ATOMIC                      COORDINATES (BOHR)\n");
  printf("           CHARGE         X                   Y                   Z\n");
  for (i=0; i<data->numatoms; i++) {
    printf(" %-8s %6d", data->initatoms[i].type, data->initatoms[i].atomicnum);
    
    printf("%17.10f",   ANGS_TO_BOHR*data->initatoms[i].x);
    printf("%20.10f",   ANGS_TO_BOHR*data->initatoms[i].y);
    printf("%20.10f\n", ANGS_TO_BOHR*data->initatoms[i].z);
  }
  printf("\n");
  printf("     ATOMIC BASIS SET\n");
  printf("     ----------------\n");
  printf(" THE CONTRACTED PRIMITIVE FUNCTIONS HAVE BEEN UNNORMALIZED\n");
  printf(" THE CONTRACTED BASIS FUNCTIONS ARE NOW NORMALIZED TO UNITY\n");
  printf("\n");
  printf("  SHELL TYPE  PRIMITIVE        EXPONENT          CONTRACTION COEFFICIENT(S)\n");
  printf("\n");

#if 0
  for (i=0; i<data->numatoms; i++) {
    printf("%-8s\n\n", data->initatoms[i].type);
    printf("\n");
    printf("nshells=%d\n", data->num_shells_per_atom[i]);

    for (j=0; j<data->num_shells_per_atom[i]; j++) {
      printf("nprim=%d\n", data->num_prim_per_shell[shellcount]);

      for (k=0; k<data->num_prim_per_shell[shellcount]; k++) {
        printf("%6d   %d %7d %22f%22f\n", j, data->shell_symmetry[shellcount],
               primcount+1, data->basis[2*primcount], data->basis[2*primcount+1]);
        primcount++;
      }

      printf("\n");
      shellcount++;
    }
  }
#endif
  printf("gamessplugin) =================================================================\n");
  for (i=0; i<data->num_basis_atoms; i++) {
    printf("%-8s (%10s)\n\n", data->initatoms[i].type, data->basis_set[i].name);
    printf("\n");

    for (j=0; j<data->basis_set[i].numshells; j++) {

      for (k=0; k<data->basis_set[i].shell[j].numprims; k++) {
        printf("%6d   %d %7d %22f%22f\n", j,
               data->basis_set[i].shell[j].symmetry,
               primcount+1,
               data->basis_set[i].shell[j].prim[k].exponent,
               data->basis_set[i].shell[j].prim[k].contraction_coeff);
        primcount++;
      }

      printf("\n");
      shellcount++;
    }
  }
  printf("\n");
  printf(" TOTAL NUMBER OF BASIS SET SHELLS             =%5d\n", data->num_shells);
  printf(" NUMBER OF CARTESIAN GAUSSIAN BASIS FUNCTIONS =%5d\n", data->wavef_size);
  printf(" NUMBER OF ELECTRONS                          =%5d\n", data->num_electrons);
  printf(" CHARGE OF MOLECULE                           =%5d\n", data->totalcharge);
  printf(" SPIN MULTIPLICITY                            =%5d\n", data->multiplicity);
  printf(" NUMBER OF OCCUPIED ORBITALS (ALPHA)          =%5d\n", data->num_occupied_A);
  printf(" NUMBER OF OCCUPIED ORBITALS (BETA )          =%5d\n", data->num_occupied_B);
  printf(" TOTAL NUMBER OF ATOMS                        =%5i\n", data->numatoms);
  printf("\n");
}



/**********************************************************
 *
 * this subroutine checks if the provided files is
 * actually a GAMESS file;
 *
 **********************************************************/
static int have_gamess(gamessdata *data) 
{
  char word[3][BUFSIZ];
  char buffer[BUFSIZ];
  int day, year;
  char month[BUFSIZ], rev[BUFSIZ];
  int i = 0;
 
  buffer[0] = '\0';
  for (i=0; i<3; i++) word[i][0] = '\0';


  /* check if the file is GAMESS format 
   * for now I just read line by line until 
   * the word GAMESS appears                */
  do {
    GET_LINE(buffer, data->file);
    sscanf(buffer,"%s %s %s",&word[0][0],&word[1][0],&word[2][0]);

  } while( strcmp(&word[1][0],"GAMESS") || 
           strcmp(&word[2][0],"VERSION") );

  /* extract the version number if possible; otherwise
   * return empty string */
  if (strstr(buffer,"=") != NULL) {
    strncpy(data->version_string,strstr(buffer,"=")+2,16);
    data->version_string[16] = '\0';
  }
  
  /* determine if we're dealing with pre-"27 JUN 2005"
   * version */
  sscanf(data->version_string,"%d %s %d %s",&day, month, &year, rev);
  
  if ( ( year >= 2006 ) ||
       ( year == 2005 && !strcmp(month,"JUN") ) ||
       ( year == 2005 && !strcmp(month,"NOV") ) ||
       ( year == 2005 && !strcmp(month,"DEC") ) )
  {
    data->version = 2;
  }
  else { 
    data->version = 1;
  }

  /* scan the next 20 lines and see if we're dealing with
   * PC Gamess output */
  for (i=0; i<20; i++) {
    GET_LINE(buffer, data->file);
    sscanf(buffer,"%s %s %s",&word[0][0],&word[1][0],&word[2][0]);

    if (!strcmp(&word[1][0],"PC") && !strcmp(&word[2][0],"GAMESS"))
    {
      data->have_pcgamess = 1;
      break;
    }
  }

  /* short messsage to stdout */
  if ( data->have_pcgamess == 1 ) {
    printf("gamessplugin) Detected PC GAMESS output :)\n");
  }
  else {
    printf("gamessplugin) Detected GAMESS format :)\n");
  }
  printf("gamessplugin) GAMESS version = %s \n", 
         data->version_string);

  return TRUE;
}


/**********************************************************
 *
 * this subroutine reads the number of procs and the amount
 * of memory requested
 *
 **********************************************************/
static int get_proc_mem(gamessdata *data) {

  char word[4][BUFSIZ];
  char buffer[BUFSIZ];
  char *temp;
  int nproc;
  int i;

  buffer[0] = '\0';
  for (i=0; i<3; i++) word[i][0] = '\0';

  rewind(data->file);

  
  /* scan for the number of processors; here we need
   * distinguish between vanilla Gamess and PC Gamess */
  if (data->have_pcgamess == 1) {
    /* for now we fake ncpu = 1 until we know exactly
     * how the output format looks like */
    nproc = 1;
  }
  else {
    do {
      GET_LINE(buffer, data->file);
      sscanf(buffer,"%s %d %s",&word[0][0],&nproc,&word[1][0]);

      if (!strcmp(&word[0][0],"Initiating") &&
          !strcmp(&word[1][0],"compute")) {
        break;
      }

      /* Some versions */
      if (!strcmp(&word[0][0],"PARALLEL") &&
          !strcmp(&word[0][1],"RUNNING")) {
        sscanf(buffer,"%*s %*s %*s %*s %d %*s",&nproc);
        break;
      }

    } while (strcmp(&word[0][0],"ECHO") || 
             strcmp(&word[1][0],"THE") );
  }
  
  /* store the number of processors */
  data->nproc = nproc;

  
  /* scan for the amount of memory requested */
  do {
    GET_LINE(buffer, data->file);
    sscanf(buffer,"%s %s",&word[0][0],&word[1][0]);

  } while( strcmp(&word[0][0],"$SYSTEM") || 
           strcmp(&word[1][0],"OPTIONS") );

  eatline(data->file, 1);


  /* next line contains the amount of memory requested,
   * vanilla Gamess and PC Gamess need separate treatment */
  if (data->have_pcgamess == 1) {
    GET_LINE(buffer, data->file);

    /* store it */
    if ((temp = strstr(buffer,"MEMORY=")+8)==NULL) return FALSE;
    strncpy(data->memory,chop_string_all(temp),sizeof(data->memory));
  }
  else {
    GET_LINE(buffer, data->file);
    sscanf(buffer,"%s %s %s",&word[0][0],&word[1][0],&word[2][0]);

    /* store it */
    strncpy(data->memory,&word[2][0],sizeof(data->memory));
  }

  printf("gamessplugin) GAMESS used %d compute processes \n", nproc);
  printf("gamessplugin) GAMESS used %s words of memory \n", data->memory);

  return TRUE;
}


/**********************************************************
 *
 * Extract the $BASIS options
 *
 **********************************************************/
static int get_basis_options(gamessdata *data) {

  char word[3][BUFSIZ];
  char buffer[BUFSIZ];
  char diffuse[BUFSIZ];
  char polarization[BUFSIZ];
  int i = 0;

  buffer[0] = '\0';
  diffuse[0] = '\0';
  polarization[0] = '\0';
  for (i=0; i<3; i++) word[i][0] = '\0';

  /* to be safe let's rewind the file */
  rewind(data->file);

  /* start scanning */
  if (!goto_keyline(data->file, "BASIS OPTIONS", "RUN TITLE")) {
    /* No Basis options section found
     * (basis was entered explicitly) */
    return TRUE;
  }

  eatline(data->file, 2);


  /* the first string in the current line contains the
   * GBASIS used; copy it over into the gbasis variable
   * of gamessdata */
  GET_LINE(buffer, data->file);
  sscanf(buffer,"%s %s %s", &word[0][0], &word[1][0], &word[2][0]);
 
  strncpy(data->gbasis, (&word[0][0])+7, sizeof(data->gbasis));


  /* in case we're using a pople style basis set, i.e. 
   * GBASIS=N311,N31,N21 or STO we also scan for the number 
   * of gaussians, as well as p,d,f and diffuse functions
   * and use this info to assemble a "basis set string" */
  if ( !strncmp(data->gbasis,"N311",sizeof(data->gbasis)) ||
       !strncmp(data->gbasis,"N31",sizeof(data->gbasis)) ||
       !strncmp(data->gbasis,"N21",sizeof(data->gbasis)) ||
       !strncmp(data->gbasis,"STO",sizeof(data->gbasis)) ) 
  {
    int ngauss, npfunc, ndfunc, nffunc;
    int diffs=FALSE, diffsp=FALSE;
    char torf;

    /* word[2] read previously should still contain the
     * number of gaussians used */
    ngauss = atoi(&word[2][0]);


    /* the next line gives us the d,f and diffuse sp
     * functions */
    GET_LINE(buffer, data->file);
    sscanf(buffer,"%*s %d %*s %d %*s %c", &ndfunc, &nffunc, &torf);

    /* convert GAMESS' .TRUE./.FALSE. for DIFFSP into 1/0 */
    if (torf=='T') diffsp = TRUE;


    /* the next line gives us the p and diffuse s
     * functions */
    GET_LINE(buffer, data->file);
    sscanf(buffer,"%*s %d %*s %c", &npfunc, &torf);

    /* convert GAMESS' .TRUE./.FALSE. for DIFFSP into 1/0 */
    if (torf=='T') diffs = TRUE;


    /* now we need some logic to assemble this info into
     * some nice looking string a la "6-31G*" */

    /* create the diffuse function string */
    if (diffs && diffsp) {
    	strncpy(diffuse,"++",sizeof(diffuse));
    }
    else if (diffsp) {
    	strncpy(diffuse,"+",sizeof(diffuse));
    }
    else {
    	strncpy(diffuse,"",sizeof(diffuse));
    }


    /* create the polarization function string */
    if (npfunc>0 && ndfunc>0 && nffunc>0) {
      sprintf(polarization, "(%dp,%dd,%df)", npfunc, ndfunc, nffunc);
    } else if (npfunc>0 && ndfunc>0) {
      sprintf(polarization, "(%dp,%dd)", npfunc, ndfunc);
    }
    else if (npfunc>0) {
      sprintf(polarization, "(%dp)", npfunc);
    }
    else if (ndfunc>0) {
      sprintf(polarization, "(%dd)", ndfunc);
    } 
    else {
      strncpy(polarization, "", sizeof(polarization));
    } 

    /* assemble the bits */ 
    if (!strcmp(data->gbasis, "STO")) {
      sprintf(data->basis_string, "STO-%dG%s%s",
              ngauss, diffuse, polarization);
    }
    else {
      sprintf(data->basis_string, "%d-%s%sG%s",
              ngauss, (data->gbasis+1), diffuse, 
              polarization);
    }      
  }

  /* cc-pVnZ and cc-pCVnZ */
  else if (!strncmp(data->gbasis, "CC",  2)) {
    strcpy(data->basis_string, "cc-p");
    if (strlen(data->gbasis)==4 && data->gbasis[3]=='C') {
      strcat(data->basis_string, "C");
    }
    strcat(data->basis_string, "V");
    strncat(data->basis_string, &data->gbasis[2], 1);
    strcat(data->basis_string, "Z");
  }

  /* aug-cc-pVnZ and aug-cc-pCVnZ */
  else if (!strncmp(data->gbasis, "ACC", 3)) {
    strcpy(data->basis_string, "aug-cc-p");
    if (strlen(data->gbasis)==5 && data->gbasis[4]=='C') {
      strcat(data->basis_string, "C");
    }
    strcat(data->basis_string, "V");
    strncat(data->basis_string, &data->gbasis[3], 1);
    strcat(data->basis_string, "Z");
  }

  /* for non Pople style basis sets we just use the GBASIS
   * for the basis string;
   * TODO: make the basis_string more comprehensive for non
   *       pople-style basis sets */
  else {
    strncpy(data->basis_string,data->gbasis,
            sizeof(data->basis_string));
  }

  return TRUE;
}



/**********************************************************
 *
 * Extract the run title line
 *
 **********************************************************/
static int get_runtitle(gamessdata *data) {

  char word[2][BUFSIZ];
  char buffer[BUFSIZ];

  /* initialize arrays */
  word[0][0] = '\0';
  word[1][0] = '\0';
  buffer[0] = '\0';

  /* look for RUN TITLE section */
  do {
    GET_LINE(buffer, data->file);
    sscanf(buffer,"%s %s",&word[0][0],&word[1][0]);

  } while (strcmp(&word[0][0],"RUN") || 
           strcmp(&word[1][0],"TITLE"));

  eatline(data->file, 1);

  GET_LINE(buffer, data->file);
  strncpy(data->runtitle,chop_string_nl(buffer),sizeof(data->runtitle));

  return TRUE;
} 


/* Read the input atom definitions and geometry */
static int get_input_structure(gamessdata *data) {
  char buffer[BUFSIZ];
  char word[4][BUFSIZ];
  int i, numatoms;
  int bohr;

  buffer[0] = '\0';
  for (i=0; i<4; i++) word[i][0] = '\0';

  do {
    GET_LINE(buffer, data->file);
    sscanf(buffer,"%s %s %s %s",&word[0][0],&word[1][0],
           &word[2][0],&word[3][0]);

  } while(strcmp(&word[0][0],"ATOM") || 
          strcmp(&word[1][0],"ATOMIC"));

  
  /* test if coordinate units are Bohr */
  bohr = (!strcmp(&word[3][0],"(BOHR)"));

  eatline(data->file, 1);

  /* we don't know the number of atoms yet */
  numatoms=-1;  

  /* Read the coordinate block */
  if (get_coordinates(data->file, &data->initatoms, bohr, &numatoms))
    data->num_frames_read = 0;
  else return FALSE;


  /* store number of atoms in data structure */
  data->numatoms = numatoms;

  return TRUE; 
}


/**********************************************************
 *
 * Read data from the $CONTRL card
 *
 **********************************************************/
static int get_contrl(gamessdata *data) {

  char word[3][BUFSIZ];
  char buffer[BUFSIZ];
  char *temp;

  word[0][0] = '\0';
  word[1][0] = '\0';
  word[2][0] = '\0';
  buffer[0] = '\0';


  /* start scanning; currently we support
   * RUNTYP = ENERGY, OPTIMIZE, SADPOINT, HESSIAN, SURFACE */
  do {
    GET_LINE(buffer, data->file);
    sscanf(buffer,"%s %s",&word[0][0],&word[1][0]);

  } while( strcmp(&word[0][0],"\044CONTRL") || 
           strcmp(&word[1][0],"OPTIONS") );

  eatline(data->file, 1);

  /* current line contains SCFTYP, RUNTYP, EXETYP info; scan it */
  GET_LINE(buffer, data->file);
  sscanf(buffer,"%s %s",&word[0][0],&word[1][0]);

  /* check for supported RUNTYPs */
  if      (!strcmp(&word[1][0],"RUNTYP=ENERGY")) {
    data->runtype = ENERGY;
  }
  else if (!strcmp(&word[1][0],"RUNTYP=OPTIMIZE")) {
    data->runtype = OPTIMIZE;
  }
  else if (!strcmp(&word[1][0],"RUNTYP=SADPOINT")) {
    data->runtype = SADPOINT;
  }
  else if (!strcmp(&word[1][0],"RUNTYP=HESSIAN")) {
    data->runtype = HESSIAN;
  }
  else if (!strcmp(&word[1][0],"RUNTYP=SURFACE")) {
    data->runtype = SURFACE;
  }
  else if (!strcmp(&word[1][0],"RUNTYP=GRADIENT")) {
    data->runtype = GRADIENT;
  }
  else {
    printf("gamessplugin) The %s is currently not supported \n",
           &word[1][0]);
  /*   return FALSE; */
  }
  printf("gamessplugin) File generated via %s \n",&word[1][0]);


  /* determine SCFTYP */
  if (!strcmp(&word[0][0],"SCFTYP=RHF")) {
    data->scftype = RHF;
  }
  else if (!strcmp(&word[0][0],"SCFTYP=UHF")) {
    data->scftype = UHF;
  }
  else if (!strcmp(&word[0][0],"SCFTYP=ROHF")) {
    data->scftype = ROHF;
  }
  else if (!strcmp(&word[0][0],"SCFTYP=GVB")) {
    data->scftype = GVB;
  }
  else if (!strcmp(&word[0][0],"SCFTYP=MCSCF")) {
    data->scftype = MCSCF;
  }
  else if (!strcmp(&word[0][0],"SCFTYP=NONE")) {
    data->scftype = NONE;
  }
  else {
    /* if we don't find a supported SCFTYP we bomb out; this
     * might be a little drastic */
    printf("gamessplugin) %s is currently not supported \n",
           &word[0][0]);
    return FALSE;
  }
  printf("gamessplugin) Type of wavefunction used %s \n",
         &word[0][0]);

  /* scan for MPLEVL, CITYP, CCTYP, VBTYP info; */
  GET_LINE(buffer, data->file);
  sscanf(buffer,"%s %s %*s %s",&word[0][0],&word[1][0],&word[2][0]);

  if (!strcmp(&word[0][0],"MPLEVL=")) {
    /* Moller-Plesset perturbation level */
    printf("gamessplugin) MP perturbation level %s \n",&word[1][0]);
    data->mplevel = atoi(&word[1][0]);

    /* determine CITYP */
    if      (!strcmp(&word[2][0],"=NONE"))  data->citype = NONE;
    else if (!strcmp(&word[2][0],"=CIS"))   data->citype = CIS;
    else if (!strcmp(&word[2][0],"=ALDET")) data->citype = ALDET;
    else if (!strcmp(&word[2][0],"=ORMAS")) data->citype = ORMAS;
    else if (!strcmp(&word[2][0],"=GUGA"))  data->citype = GUGA;
    else if (!strcmp(&word[2][0],"=FSOCI")) data->citype = FSOCI;
    else if (!strcmp(&word[2][0],"=GENCI")) data->citype = GENCI;
    else                                    data->citype = UNKNOWN;
    printf("gamessplugin) CI method %s \n",&word[2][1]);

    GET_LINE(buffer, data->file);
    sscanf(buffer,"%s %s",&word[0][0],&word[1][0]);
  }

  /* scan for DFTTYP, TDDFT info; */
  if (!strncmp(&word[0][0],"DFTTYP=", 7)) {
    printf("gamessplugin) Density functional used is %s \n",&word[0][7]);
    GET_LINE(buffer, data->file);
  }


  /* find the coordinate type in next line */
  while ( (temp=strstr(buffer,"COORD =")) == NULL ) {
    GET_LINE(buffer, data->file);;
  }
  strncpy(data->geometry,chop_string_all(temp+7),
          sizeof(data->geometry)); 
  printf("gamessplugin) Coordinate type used is %s \n", data->geometry);

  return TRUE;
}

/* Read input parameters regarding calculation of 
 * certain molecular properties such as electrostatic
 * moments and the MEP. */
static int get_properties_input(gamessdata *data) {
  /* TODO!! */
  return TRUE;
}

/* Read symmetry point group and highest axis */
static int get_symmetry(gamessdata *data) {
  /* TODO!! */
  /* This could be lumped together with 
   * get_guess_options() where we deal with orbital
   * symmetry and which comes right after the pointgroup
   * stuff */
  return TRUE;
}


static int read_first_frame(gamessdata *data) {
  /* the angular momentum is populated in get_wavefunction 
   * which is called by get_traj_frame(). We have obtained
   * the array size wavef_size already from the basis set
   * statistics */
  data->angular_momentum = (int*)calloc(3*data->wavef_size, sizeof(int));

  /* Try reading the first frame. 
   * If there is only one frame then also read the
   * final wavefunction. */
  if (!get_traj_frame(data, data->initatoms, data->numatoms)) {
    return FALSE;
  }

  return TRUE;
}

/******************************************************
 *
 * Reads the info printed after the geometry search
 * has finished or whatever analysis was done in a 
 * single point run, e.g. ESP charges, Hessian, etc.
 * Rewinds to the beginning of the search when done,
 * because we read this part at in the initial phase
 * and might have to look for additional timesteps
 * later.
 *
 ******************************************************/
static int get_final_properties(gamessdata *data) {
  long filepos;
  filepos = ftell(data->file);

  /* Go to end of trajectory */
  fseek(data->file, data->end_of_traj, SEEK_SET);


  if (get_esp_charges(data)) {
    printf("gamessplugin) ESP charges found!\n");
  }

  if (data->runtype == GRADIENT) {
    /* get_singlepoint_gradient(data); */
  }


  if (data->runtype == HESSIAN) {
    /* try reading the hessian matrix in internal and
     * cartesian coordinates as well as the internal
     * coordinates together with their associated
     * force constants */
    if (!get_int_coords(data)) {
      printf("gamessplugin) \n");
      printf("gamessplugin) Could not determine the internal \n");
      printf("gamessplugin) coordinate and Hessian info!! \n");
      printf("gamessplugin) \n");
    }
    
    if (!get_cart_hessian(data)) {
      printf("gamessplugin) \n");
      printf("gamessplugin) Could not determine the cartesian \n");
      printf("gamessplugin) Hessian matrix!! \n");
      printf("gamessplugin) \n");
    }

    /* read the wavenumbers, intensities of the normal modes 
     * as well as the modes themselves */
    if (!get_normal_modes(data)) {
      printf("gamessplugin) \n");
      printf("gamessplugin) Could not scan the normal modes,\n");
      printf("gamessplugin) \n");
    }
  }

  /* Read localized orbitals if there are any */
  read_localized_orbitals(data);


  fseek(data->file, filepos, SEEK_SET);
  return TRUE; 
}


static int read_localized_orbitals(gamessdata *data) {
  int i;
  qm_timestep_t *ts;
  qm_wavefunction_t *wavef;

  /* Move past the listing of the canonical MOs */
  goto_keyline(data->file, "ENERGY COMPONENTS", NULL);
  eatline(data->file, 1);

  ts = data->qm_timestep + data->num_frames-1;

  for (i=0; i<2; i++) {
    wavef = add_wavefunction(ts);

    if (get_wavefunction(data, ts, wavef) == FALSE ||
        (wavef->type!=MOLFILE_WAVE_BOYS &&
         wavef->type!=MOLFILE_WAVE_PIPEK &&
         wavef->type!=MOLFILE_WAVE_RUEDEN)) {
      del_wavefunction(ts);
      return FALSE;
    }
    else {
      char typestr[64];
      if (wavef->spin==SPIN_ALPHA) {
        strcpy(typestr, "alpha");
      }
      else if (wavef->spin==SPIN_BETA) {
        strcpy(typestr, "beta");
      }
      printf("gamessplugin) Localized orbitals (%s) found for timestep %d\n",
             typestr, data->num_frames-1);
    }
  }

  return TRUE;
}


static int get_coordinates(FILE *file, qm_atom_t **atoms, int unit,
                           int *numatoms) {
  int i = 0;
  int growarray = 0;

  if (*numatoms<0) {
    *atoms = (qm_atom_t*)calloc(1, sizeof(qm_atom_t));
    growarray = 1;
  }

  /* Read in the coordinates until an empty line is reached.
   * We expect 5 entries per line */
  while (1) {
    char buffer[BUFSIZ];
    char atname[BUFSIZ];
    float atomicnum;
    float x,y,z;
    int n;
    qm_atom_t *atm;

    GET_LINE(buffer, file);
    n = sscanf(buffer,"%s %f %f %f %f",atname,&atomicnum,&x,&y,&z);

    if (n!=5) break;

    if (growarray && i>0) {
      *atoms = (qm_atom_t*)realloc(*atoms, (i+1)*sizeof(qm_atom_t));
    }
    atm = (*atoms)+i;

    strncpy(atm->type, atname, sizeof(atm->type));
    atm->atomicnum = floor(atomicnum+0.5); /* nuclear charge */
    
    /* if coordinates are in Bohr convert them to Angstrom */
    if (unit==BOHR) {
      x *= BOHR_TO_ANGS;
      y *= BOHR_TO_ANGS;
      z *= BOHR_TO_ANGS;
    }
    
    atm->x = x;
    atm->y = y;
    atm->z = z; 
    i++;
  }

  /* If file is broken off in the middle of the coordinate block 
   * we cannot use this frame. */
  if (*numatoms>=0 && *numatoms!=i) return FALSE;

  (*numatoms) = i;
  return TRUE;
}


/********************************************************
 *
 * Read basis set and orbital statistics such as
 * # of shells, # of A/B orbitals, # of electrons, 
 * multiplicity and total charge
 *
 ********************************************************/
static int get_basis_stats(gamessdata *data) {

  char buffer[BUFSIZ];
  char word[7][BUFSIZ];
  int i;

  buffer[0] = '\0';
  for (i=0; i<7; i++) word[i][0] = '\0';

  /* look for the orbital/charge/... info section */
  pass_keyline(data->file, "TOTAL NUMBER OF BASIS", NULL);


  /* go ahead reading the info */
  GET_LINE(buffer, data->file);
  sscanf(buffer,"%*s %*s %*s %*s %*s %*s %*s %d",
         &(data->wavef_size));

  /* read the number of electrons */
  GET_LINE(buffer, data->file);
  sscanf(buffer,"%*s %*s %*s %*s %d",
         &(data->num_electrons));

  /* read the charge of the molecule */
  GET_LINE(buffer, data->file);
  sscanf(buffer,"%*s %*s %*s %*s %d",
         &(data->totalcharge));

  /* read the multiplicity of the molecule */
  GET_LINE(buffer, data->file);
  sscanf(buffer,"%*s %*s %*s %d",
         &(data->multiplicity));

  /* read number of A orbitals */
  GET_LINE(buffer, data->file);

  /* note the different number of items per line for A/B orbitals
   * due to "(ALPHA)" and "(BETA )" !! */
  sscanf(buffer,"%*s %*s %*s %*s %*s %*s %d",
         &(data->num_occupied_A)); 
    
  /* read number of B orbitals */
  GET_LINE(buffer, data->file);
  sscanf(buffer,"%*s %*s %*s %*s %*s %*s %*s %d",
         &(data->num_occupied_B)); 


  printf("gamessplugin) Number of Electrons: %d \n",
      data->num_electrons);

  printf("gamessplugin) Charge of Molecule : %d \n",
      data->totalcharge);

  printf("gamessplugin) Multiplicity of Wavefunction: %d \n",
      data->multiplicity);

  printf("gamessplugin) Number of occupied A / B orbitals: %d / %d \n",\
      data->num_occupied_A, data->num_occupied_B);

  printf("gamessplugin) Number of gaussian basis functions: %d \n",\
      data->wavef_size);

 
  return TRUE;
}



/*******************************************************
 *
 * Reads in the $GUESS options.
 *
 * ******************************************************/
static int get_guess_options(gamessdata *data)
{
  char word[2][BUFSIZ];
  char buffer[BUFSIZ];
  int i = 0;
  long filepos;
  filepos = ftell(data->file);

  /* initialize buffers */
  buffer[0] = '\0';
  for (i=0; i<2; i++) word[i][0] = '\0';

  /* parse for GUESS field */
  if (!pass_keyline(data->file, "GUESS OPTIONS", NULL)) {
    printf("gamessplugin) No GUESS OPTIONS found!\n");
    fseek(data->file, filepos, SEEK_SET);
    return FALSE;
  }

  /* next line contains all we need */
  eatline(data->file, 1);
  GET_LINE(buffer, data->file);
  sscanf(buffer,"%s %s",&word[0][0],&word[1][0]);

  /* store it */
  strncpy(data->guess,(&word[1][0])+1,sizeof(data->guess));

  printf("gamessplugin) Run was performed with GUESS = %s \n",
	  data->guess);

 /* Since this block occurs in the middle of first frame
  * we need to rewind. */
  fseek(data->file, filepos, SEEK_SET);

  return TRUE;
}



/*******************************************************
 *
 * this function reads in the basis set data 
 *
 * ******************************************************/
int get_basis(gamessdata *data) {

  char buffer[BUFSIZ];
  char word[4][BUFSIZ];
  int i = 0; 
  int success = 0;
  int numread, numshells;
  shell_t *shell;
  long filepos;

  if (!strcmp(data->gbasis, "MNDO") ||
      !strcmp(data->gbasis, "AM1")  ||
      !strcmp(data->gbasis, "PM3")) {
    /* Semiempirical methods are based on STOs.
     * The only parameter we need for orbital rendering
     * are the exponents zeta for S, P, D,... shells for
     * each atom. Since GAMESS doesn't print these values
     * we skip reading the basis set but and hardcode the
     * parameters in tables in VMD. */
    return MOLFILE_ERROR;
  }

  /* initialize buffers */
  buffer[0] = '\0';
  for (i=0; i<3; i++) word[i][0] = '\0';
  
  /* Search for "ATOMIC BASIS SET" line */
  /* XXX what is a good stopstring? */
  if (!pass_keyline(data->file, "ATOMIC BASIS SET", NULL))
    printf("gamessplugin) No basis set found!\n");


  /* skip the next 6 lines */
  for(i=0; i<5; i++) { eatline(data->file, 1); }

  /* Allocate space for the basis for all atoms */
  /* When the molecule is symmetric the actual number atoms with
   * a basis set could be smaller */
  data->basis_set = (basis_atom_t*)calloc(data->numatoms, sizeof(basis_atom_t));


  i = 0; /* basis atom counter */

  do {
    prim_t *prim = NULL;
    char shellsymm;
    int numprim = 0;
    int icoeff = 0;
    filepos = ftell(data->file);
    GET_LINE(buffer, data->file);
      
    /* Count the number of relevant words in the line. */
    numread = sscanf(buffer,"%s %s %s %s",&word[0][0], &word[1][0],
           &word[2][0], &word[3][0]);

    switch (numread) {
      case 1:
        /* Next atom */
        strcpy(data->basis_set[i].name, &word[0][0]);

        /* skip initial blank line */
        eatline(data->file, 1);

        /* read the basis set for the current atom */
        shell = (shell_t*)calloc(1, sizeof(shell_t)); 
        numshells = 0;

        do {
          filepos = ftell(data->file);
          numprim = read_shell_primitives(data, &prim, &shellsymm, icoeff);

          if (numprim>0) {
            /* make sure we have eiter S, L, P, D, F or G shells */
            if ( (shellsymm!='S' && shellsymm!='L' && shellsymm!='P' && 
                  shellsymm!='D' && shellsymm!='F' && shellsymm!='G') ) {
              printf("gamessplugin) WARNING ... %c shells are not supported \n", shellsymm);
            }
            
            /* create new shell */
            if (numshells) {
              shell = (shell_t*)realloc(shell, (numshells+1)*sizeof(shell_t));
            }
            shell[numshells].numprims = numprim;
            shell[numshells].symmetry = shellsymm_int(shellsymm);
            shell[numshells].prim = prim;
            data->num_basis_funcs += numprim;

            /* We split L-shells into one S and one P-shell.
             * I.e. for L-shells we have to go back read the shell again
             * this time using the second contraction coefficients.
             * We use L and M instead of S and P for the shell symmetry
             * in order to be able to distinguish SP-type shells. */
            if (shellsymm=='L' && !icoeff) {
              fseek(data->file, filepos, SEEK_SET);
              icoeff++;
            } else if (shellsymm=='L' && icoeff) {
              shell[numshells].symmetry = SP_P_SHELL;
              icoeff = 0;
            }

            numshells++;
          }
        } while (numprim);

        /* store shells in atom */
        data->basis_set[i].numshells = numshells;
        data->basis_set[i].shell = shell;

        /* store the total number of basis functions */
        data->num_shells += numshells;
        i++;

        /* go back one line so that we can read the name of the
         * next atom */
        fseek(data->file, filepos, SEEK_SET);

        break;

      case 4:
        /* this is the very end of the basis set */
        if (!strcmp(&word[0][0],"TOTAL")  &&
            !strcmp(&word[1][0],"NUMBER") && 
            !strcmp(&word[2][0],"OF")     &&
            !strcmp(&word[3][0],"BASIS")) {
          success = 1;
          /* go back one line so that get_basis_stats()
             can use this line as a keystring. */
          fseek(data->file, filepos, SEEK_SET);
        }
        break;
    }

  } while (!success);


  printf("gamessplugin) Parsed %d uncontracted basis functions for %d atoms.\n",
         data->num_basis_funcs, i);

  data->num_basis_atoms = i;


  /* allocate and populate flat arrays needed for molfileplugin */
  return fill_basis_arrays(data);
}


/**************************************************
 *
 * Convert shell symmetry type from char to int.
 *
 ************************************************ */
static int shellsymm_int(char symm) {
  int shell_symmetry;

  switch (symm) {
    case 'L':
      shell_symmetry = SP_S_SHELL;
      break;
    case 'M':
      shell_symmetry = SP_P_SHELL;
      break;
    case 'S':
      shell_symmetry = S_SHELL;
      break;
    case 'P':
      shell_symmetry = P_SHELL;
      break;
    case 'D':
      shell_symmetry = D_SHELL;
      break;
    case 'F':
      shell_symmetry = F_SHELL;
      break;
    case 'G':
      shell_symmetry = G_SHELL;
      break;
    default:
      shell_symmetry = UNK_SHELL;
      break;
  }

  return shell_symmetry;
}



/******************************************************
 *
 * Populate the flat arrays containing the basis
 * set data.
 *
 ******************************************************/
static int fill_basis_arrays(gamessdata *data) {
  int i, j, k;
  int shellcount = 0;
  int primcount = 0;
  int found = 0;

  float *basis;
  int *num_shells_per_atom;
  int *num_prim_per_shell;
  int *shell_symmetry;
  int *atomicnum_per_basisatom;

  /* Count the total number of primitives which
   * determines the size of the basis array. */
  for(i=0; i<data->num_basis_atoms; i++) {
    for (j=0; j<data->basis_set[i].numshells; j++) {
      primcount += data->basis_set[i].shell[j].numprims;
    }
  }

  /* reserve space for pointer to array containing basis
   * info, i.e. contraction coeficients and expansion 
   * coefficients; need 2 entries per basis function, i.e.
   * exponent and contraction coefficient; also,
   * allocate space for the array holding the orbital symmetry
   * information per primitive Gaussian.
   * Finally, initialize the arrays holding the number of 
   * shells per atom and the number of primitives per shell*/
  basis = (float *)calloc(2*primcount,sizeof(float));

  /* make sure memory was allocated properly */
  if (basis == NULL) {
    PRINTERR;
    return MOLFILE_ERROR;
  }

  shell_symmetry = (int *)calloc(data->num_shells, sizeof(int));
  
  /* make sure memory was allocated properly */
  if (shell_symmetry == NULL) {
    PRINTERR; 
    return MOLFILE_ERROR;
  }

  num_shells_per_atom = (int *)calloc(data->num_basis_atoms, sizeof(int));

  /* make sure memory was allocated properly */
  if (num_shells_per_atom == NULL) {
    PRINTERR; 
    return MOLFILE_ERROR;
  }

  num_prim_per_shell = (int *)calloc(data->num_shells, sizeof(int));

  /* make sure memory was allocated properly */
  if (num_prim_per_shell == NULL) {
    PRINTERR;
    return MOLFILE_ERROR;
  }

  atomicnum_per_basisatom = (int *)calloc(data->num_basis_atoms, sizeof(int));

  /* make sure memory was allocated properly */
  if (atomicnum_per_basisatom == NULL) {
    PRINTERR;
    return MOLFILE_ERROR;
  }


  /* store pointers in struct gamessdata */
  data->basis = basis;
  data->shell_symmetry = shell_symmetry;
  data->num_shells_per_atom = num_shells_per_atom;
  data->num_prim_per_shell = num_prim_per_shell;
  data->atomicnum_per_basisatom = atomicnum_per_basisatom;

  primcount = 0;
  for (i=0; i<data->num_basis_atoms; i++) {
    /* assign atomic number */
    for(j=0; j<data->numatoms; j++) {
      if (!strcmp(data->initatoms[j].type, data->basis_set[i].name)) {
        found = 1;
        break;
      }
    }
    if (!found) {
      printf("gamessplugin) ERROR: Couldn't find atomic number for basis set atom %s",
             data->initatoms[j].type);
      return MOLFILE_ERROR;
    }
    data->basis_set[i].atomicnum = data->initatoms[j].atomicnum;
    atomicnum_per_basisatom[i] = data->initatoms[j].atomicnum;

    num_shells_per_atom[i] = data->basis_set[i].numshells;

    for (j=0; j<data->basis_set[i].numshells; j++) {
      shell_symmetry[shellcount] = data->basis_set[i].shell[j].symmetry;
      num_prim_per_shell[shellcount] = data->basis_set[i].shell[j].numprims;

      for (k=0; k<data->basis_set[i].shell[j].numprims; k++) {
        basis[2*primcount  ] = data->basis_set[i].shell[j].prim[k].exponent;
        basis[2*primcount+1] = data->basis_set[i].shell[j].prim[k].contraction_coeff;
        primcount++;
      }
      shellcount++;
    }
  } 

  return TRUE;
}


/******************************************************
 *
 * read all primitives for the current shell
 *
 ******************************************************/
static int read_shell_primitives(gamessdata *data, prim_t **prim, char *shellsymm,
                                 int icoeff) {
  char buffer[BUFSIZ];
  float exponent = 0.0; 
  float contract[2] = {0.0, 0.0};
  int shell, success;
  int primcounter = 0;
  (*prim) = (prim_t*)calloc(1, sizeof(prim_t));

  do {
    GET_LINE(buffer, data->file);
      success = sscanf(buffer,"%i %c %*s %f %f %f", &shell,
                       shellsymm,
                       &exponent, &contract[0], &contract[1]); 

    /* store in basis array and increase the counter */ 
    switch (success) {
      case 4:
        if (primcounter) {
          *prim = (prim_t*)realloc(*prim, (primcounter+1)*sizeof(prim_t));
        }

        /* store exponent */
        (*prim)[primcounter].exponent = exponent;
          
        /* store coefficient */
        (*prim)[primcounter].contraction_coeff = contract[0];
        
        primcounter++;
        break;

      case 5:
        if (primcounter) {
          *prim = (prim_t*)realloc(*prim, (primcounter+1)*sizeof(prim_t));
        }

        /* store exponent */
        (*prim)[primcounter].exponent = exponent;
          
        /* store coefficient */
        (*prim)[primcounter].contraction_coeff = contract[icoeff];
        
        primcounter++;
        break;

      case -1:
        /* otherwise it's an empty line which represents the end of the shell */
        break;

      case 1:
        /* the user had given the next atom a numeric name */
        success = -1;
        break;
    }

  } while(success>0);

  if (!primcounter) free(*prim);

  return primcounter;
}



/******************************************************
 *
 * this function extracts the trajectory information
 * from the output file
 *
 * *****************************************************/
static int get_traj_frame(gamessdata *data, qm_atom_t *atoms,
                          int natoms) {
  qm_timestep_t *cur_ts;
  qm_wavefunction_t *wavef_alpha;
  char buffer[BUFSIZ];
  char word[BUFSIZ];
  buffer[0] = '\0';
  word[0]   = '\0';

  printf("gamessplugin) Timestep %i:\n", data->num_frames_read);
  printf("gamessplugin) ============\n");

  fseek(data->file, data->filepos_array[data->num_frames_read], SEEK_SET);


  if (data->runtype==OPTIMIZE || data->runtype==SADPOINT) {

    pass_keyline(data->file, "COORDINATES OF ALL ATOMS", NULL);
    eatline(data->file, 2);

    if (!get_coordinates(data->file, &data->initatoms, ANGSTROM, &natoms)) {
      printf("gamessplugin) Couldn't find coordinates for timestep %i\n", data->num_frames_read);
    }

  }
  else if (data->runtype==SURFACE) {
    if (!pass_keyline(data->file, "---- SURFACE MAPPING GEOMETRY", NULL)) {
      return FALSE;
    }
  }

 
  /* get a convenient pointer to the current qm timestep */
  cur_ts = data->qm_timestep + data->num_frames_read;

  /* read the SCF energies */
  if (get_scfdata(data, cur_ts) == FALSE) {
    printf("gamessplugin) Couldn't find SCF iterations for timestep %i\n", data->num_frames_read);
  }

  /* Allocate memory for new wavefunction */
  wavef_alpha = add_wavefunction(cur_ts);

  /* Try to read wavefunction and orbital energies */
  if (get_wavefunction(data, cur_ts, wavef_alpha) == FALSE) {
    /* Free the last wavefunction again. */
    del_wavefunction(cur_ts);
    printf("gamessplugin) No canonical wavefunction present for timestep %i\n", data->num_frames_read);

  } else {
    char typestr[64];
    strcpy(typestr, wavef_alpha->info);
    if (data->scftype==UHF) {
      strcat(typestr, ", alpha");
    }
    printf("gamessplugin) Wavefunction (%s) found for timestep %i\n", typestr, data->num_frames_read);

    if (data->scftype==UHF || data->scftype==GVB || data->scftype==MCSCF) {
      /* Try to read second wavefunction
       * (spin beta or GI orbitals or MCSCF optimized orbs) */
      qm_wavefunction_t *wavef_beta;
      wavef_beta = add_wavefunction(cur_ts);

      if (get_wavefunction(data, cur_ts, wavef_beta) == FALSE) {
        /* Free the last wavefunction again. */
        del_wavefunction(cur_ts);

        printf("gamessplugin) No beta wavefunction present for timestep %i\n", data->num_frames_read);
      } else {
        strcpy(typestr, wavef_beta->info);
        if (data->scftype==UHF) {
          strcat(typestr, ", beta");
        }
        printf("gamessplugin) Wavefunction (%s)  found for timestep %i\n", typestr, data->num_frames_read);
      }
    }
  }

  if (get_population(data, cur_ts)) {
    printf("gamessplugin) Mulliken charges found\n");
  }

  if (data->citype!=NONE) {
    if (pass_keyline(data->file, "CI DENSITY MATRIX AND NATURAL ORBITALS",
                       "GRADIENT (HARTREE/BOHR)")) {
      int i, numstates=0;
      qm_wavefunction_t *wave_ci;
      goto_keyline(data->file, "NUMBER OF STATES", NULL);
      GET_LINE(buffer, data->file);
      trimleft(buffer);
      sscanf(buffer, " NUMBER OF STATES = %d", &numstates);
      printf("gamessplugin) Number of CI states = %d\n", numstates);

      for (i=0; i<numstates; i++) {
        float cienergy = 0.f;
        goto_keyline(data->file, "CI EIGENSTATE", NULL);
        GET_LINE(buffer, data->file);
        sscanf(buffer,"%*s %*s %*d %*s %*s %*s %f", &cienergy);
        printf("gamessplugin) CI energy[%d] = %f\n", i, cienergy);

        wave_ci = add_wavefunction(cur_ts);

        if (get_wavefunction(data, cur_ts, wave_ci) == FALSE) {
          del_wavefunction(cur_ts);      
        }
        else {
          int j, canon =-1;
          wave_ci->exci = i;
          printf("gamessplugin) Found CI natural orbitals for timestep %d\n", data->num_frames_read);

          for (j=0; j<cur_ts->numwave; j++) {
            if (cur_ts->wave[j].type==MOLFILE_WAVE_CANON &&
                cur_ts->wave[j].spin==wave_ci->spin) {
              canon = j;
              break;
            }
          }
          if (canon>=0) {
            /* Replace existing canonical wavefunction for this step */
            replace_wavefunction(cur_ts, canon);
          }
        }
      }
    }
  }


  if (get_gradient(data, cur_ts) == FALSE) {
    printf("gamessplugin) No energy gradient present for timestep %i\n", data->num_frames_read);
  } else {
    printf("gamessplugin) Energy gradient found for timestep %i\n", data->num_frames_read);
  }

  /* If this is the last frame of the trajectory and the file
   * wasn't truncated and the program didn't terminate
   * abnormally then read the final wavefunction. */
  if ((data->runtype == OPTIMIZE || data->runtype == SADPOINT) &&
      (data->num_frames_read+1 == data->num_frames &&
       (data->opt_status == MOLFILE_QM_STATUS_UNKNOWN || 
        data->opt_status == MOLFILE_QM_OPT_CONVERGED ||
        data->opt_status == MOLFILE_QM_OPT_NOT_CONV))) {
    qm_wavefunction_t *wave_final;

    /* We need to jump over the end of the trajectory because 
     * this is also the keystring for get_wavefunction() to
     * bail out. */
    if (data->opt_status == MOLFILE_QM_OPT_CONVERGED || 
        data->opt_status == MOLFILE_QM_OPT_NOT_CONV) {
      fseek(data->file, data->end_of_traj, SEEK_SET);
    }

    /* Try to read final wavefunction and orbital energies
     * Any wavefunction previously stored for the final
     * timestep will be overwritten. */
    wave_final = add_wavefunction(cur_ts);

    if (get_wavefunction(data, cur_ts, wave_final) == FALSE) {
      printf("gamessplugin) No final wavefunction present for timestep %d\n", data->num_frames_read);
      /* Free the the pointer an the contents of the last wfn */
      del_wavefunction(cur_ts);

    } else {
      char typestr[MOLFILE_BUFSIZ];

      /* if there exists a canonical wavefunction of the same spin
       * we'll replace it */
      if (cur_ts->numwave>1 && wave_final->type==MOLFILE_WAVE_CANON) {
        int i, found =-1;
        for (i=0; i<cur_ts->numwave; i++) {
          if (cur_ts->wave[i].type==wave_final->type &&
              cur_ts->wave[i].spin==wave_final->spin &&
              cur_ts->wave[i].exci==wave_final->exci &&
              !strncmp(cur_ts->wave[i].info, wave_final->info, MOLFILE_BUFSIZ)) {
            found = i;
            break;
          }
        }
        if (found>=0) {
          /* If the new wavefunction has more orbitals we 
           * replace the old one for this step. */
          if (wave_final->num_orbitals > 
              cur_ts->wave[found].num_orbitals) {
            /* Replace existing wavefunction for this step */
            replace_wavefunction(cur_ts, found);
          } else {
            /* Delete last wavefunction again */
            del_wavefunction(cur_ts);
          }
          wave_final = &cur_ts->wave[cur_ts->numwave-1];
        }
      }


      strcpy(typestr, wave_final->info);
      if (data->scftype==UHF) {
        strcat(typestr, ", alpha");
      }

      printf("gamessplugin) Final wavefunction (%s) found for timestep %d\n",
             typestr, data->num_frames_read);
    }

  }

  data->num_frames_read++;

  return TRUE;
}


/* Analyze the trajectory.
 * Read the parameters comtrolling geometry search and
 * find the end of the trajectory, couinting the frames
 * on the way. Store the filepointer for the beginning of
 * each frame in *filepos_array. */
static int analyze_traj(gamessdata *data) {
  char buffer[BUFSIZ], nserch[BUFSIZ];
  char *line;
  long filepos;
  filepos = ftell(data->file);

  data->filepos_array = (long* )calloc(1, sizeof(long ));


  if (data->runtype==OPTIMIZE || data->runtype==SADPOINT) {
    pass_keyline(data->file,
                   "PARAMETERS CONTROLLING GEOMETRY SEARCH", NULL);
    eatline(data->file, 2);

    GET_LINE(buffer, data->file);
    sscanf(buffer, "NSTEP  = %i", &data->max_opt_steps);
    eatline(data->file, 3);
    GET_LINE(buffer, data->file);
    sscanf(buffer, "OPTTOL = %f", &data->opt_tol);

    /* The $STATP options are followed by the coordinates 
     * but we can skip them here because we rewind after
     * get_guess_options() and try to read them in
     * get_traj_frame(). */
  }
  else if (data->runtype==SURFACE) {
    if (pass_keyline(data->file,
                        "POTENTIAL SURFACE MAP INPUT", NULL)) {
      
      int coord1[2];
      int mplevel1=-1, nstep1;
      float origin1, disp1;
      char runtype1[BUFSIZ];
      char scftype1[BUFSIZ];
      char dfttyp1[BUFSIZ];
      char *tmp;
        
      eatline(data->file, 1);

      GET_LINE(buffer, data->file);
      /* int n=sscanf(buffer, "JOB 1 IS RUNTYP=%s SCFTYP=%s CITYP=%*s",
         runtype1, scftype1); */
      tmp = strstr(buffer, "RUNTYP=") + 7;
      sscanf(tmp, "%s", runtype1);
      tmp = strstr(buffer, "SCFTYP=") + 7;
      sscanf(tmp, "%s", scftype1);
      printf("gamessplugin) JOB 1 IS RUNTYP=%s SCFTYP=%s\n", runtype1, scftype1);

      GET_LINE(buffer, data->file);
      sscanf(buffer, "MPLEVL= %i", &mplevel1);
      tmp = strstr(buffer, "DFTTYP=") + 7;
      sscanf(tmp, "%s", dfttyp1);
      printf("gamessplugin) MPLEVL= %i DFTTYP=%s\n", mplevel1, dfttyp1);

      GET_LINE(buffer, data->file);
      sscanf(buffer, "COORD 1 LYING ALONG ATOM PAIR %i %i",
             coord1, coord1+1);
      GET_LINE(buffer, data->file);
      tmp = strstr(buffer, "ORIGIN=") + 7;
      sscanf(tmp, "%f", &origin1);
      tmp = strstr(buffer, "DISPLACEMENT=") + 13;
      sscanf(tmp, "%f", &disp1);
      tmp = strstr(buffer, "AND") + 3;
      sscanf(tmp, "%i STEPS.", &nstep1);
      printf("gamessplugin) origin=%f, displacement=%f nstep=%i\n", origin1, disp1, nstep1);
    }
  }
  else {
    /* We have just one frame */
    data->num_frames = 1;
    pass_keyline(data->file, "1 ELECTRON INTEGRALS",
                 "ENERGY COMPONENTS");
    data->filepos_array[0] = ftell(data->file);

    /* Check wether SCF has converged */
    if (pass_keyline(data->file, "SCF IS UNCONVERGED, TOO MANY ITERATIONS",
                     "ENERGY COMPONENTS")==FOUND) {
      printf("gamessplugin) SCF IS UNCONVERGED, TOO MANY ITERATIONS\n");
      data->opt_status = MOLFILE_QM_SCF_NOT_CONV;
    } else {
      data->opt_status = MOLFILE_QM_OPT_CONVERGED;
      fseek(data->file, data->filepos_array[0], SEEK_SET);
    }

    pass_keyline(data->file, "ENERGY COMPONENTS", NULL);
    data->end_of_traj = ftell(data->file);

    /* Allocate memory for the frame */
    data->qm_timestep = (qm_timestep_t *)calloc(1, sizeof(qm_timestep_t));
    memset(data->qm_timestep, 0, sizeof(qm_timestep_t));
    
    return TRUE;
  }

  printf("gamessplugin) Analyzing trajectory...\n");
  data->opt_status = MOLFILE_QM_STATUS_UNKNOWN;

  while (1) {
    if (!fgets(buffer, sizeof(buffer), data->file)) break;
    line = trimleft(buffer);

    /* at this point we have to distinguish between
     * pre="27 JUN 2005 (R2)" and "27 JUN 2005 (R2)"
     * versions since the output format for geometry
     * optimizations has changed */
    if (data->version==1) {
      strcpy(nserch, "1NSERCH=");
    }
    else if (data->version==2) {
      strcpy(nserch, "BEGINNING GEOMETRY SEARCH POINT NSERCH=");
    }

    if (strstr(line, nserch) ||
        strstr(line, "---- SURFACE MAPPING GEOMETRY")) {
      printf("gamessplugin) %s", line);

      if (data->num_frames > 0) {
        data->filepos_array = (long*)realloc(data->filepos_array,
                                (data->num_frames+1)*sizeof(long));
      }
      data->filepos_array[data->num_frames] = ftell(data->file);

      /* Make sure that we have at least a complete coordinate
         block in order to consider this a new frame. */
      if (pass_keyline(data->file, "COORDINATES OF",
               "BEGINNING GEOMETRY SEARCH POINT NSERCH=") == FOUND) {
        if (pass_keyline(data->file, "INTERNUCLEAR DISTANCES",
                           "1 ELECTRON INTEGRALS") ||
            pass_keyline(data->file, "1 ELECTRON INTEGRALS",
               "BEGINNING GEOMETRY SEARCH POINT NSERCH=")) {
          data->num_frames++;
        }
      }
    }
    else if (strstr(line, "***** EQUILIBRIUM GEOMETRY LOCATED") ||
             strstr(line, "... DONE WITH POTENTIAL SURFACE SCAN")) {
      printf("gamessplugin) ==== End of trajectory. ====\n");
      data->opt_status = MOLFILE_QM_OPT_CONVERGED;
      break;
    }
    else if (strstr(line, "***** FAILURE TO LOCATE STATIONARY POINT,")) {
      printf("gamessplugin) %s\n", line);
      if (strstr(strchr(line, ','), "SCF HAS NOT CONVERGED")) {
        data->opt_status = MOLFILE_QM_SCF_NOT_CONV;
        break;
      }
      else if (strstr(strchr(line, ','), "TOO MANY STEPS TAKEN")) {
        data->opt_status = MOLFILE_QM_OPT_NOT_CONV;
        break;
      }
    }
  }

  data->end_of_traj = ftell(data->file);
  fseek(data->file, filepos, SEEK_SET);

  if (data->opt_status == MOLFILE_QM_STATUS_UNKNOWN) {
    /* We didn't find any of the regular key strings,
     * the run was most likely broken off and we have an
     * incomplete file. */
    data->opt_status = MOLFILE_QM_FILE_TRUNCATED;
  }


  /* Allocate memory for all frames */
  data->qm_timestep = (qm_timestep_t *)calloc(data->num_frames,
                                              sizeof(qm_timestep_t));
  memset(data->qm_timestep, 0, data->num_frames*sizeof(qm_timestep_t));


  if (data->opt_status == MOLFILE_QM_SCF_NOT_CONV ||
      data->opt_status == MOLFILE_QM_FILE_TRUNCATED) {
    return FALSE;  
  }

  return TRUE;
}


/***************************************************************
 *
 * Read the number of scf iterations and the scf energies
 * for the current timestep. 
 * Assumes that the file pointer is somewhere before this:
 * ITER EX DEM     TOTAL ENERGY        E CHANGE  DENSITY CHANGE    DIIS ERROR
 * 1  0  0      -39.7266993475   -39.7266993475   0.000000118   0.000000000
 * 2  1  0      -39.7266991566     0.0000001909   0.000000032   0.000000000
 * ...
 * then it reads the block up to the next blank line.
 * The second argument is a pointer to the qm timestep you want to
 * store the data in. Memory for the scfenergies will be
 * allocated.
 *
 ***************************************************************/
static int get_scfdata(gamessdata *data, qm_timestep_t *ts) {
  char buffer[BUFSIZ];
  char word[3][BUFSIZ];
  long filepos;
  int i;
  int numread, numiter=0, dum;
  char *line;

  buffer[0] = '\0';
  for (i=0; i<3; ++i) word[i][0] = '\0';

  /* look for SCF iteration energies */
  do {
    GET_LINE(buffer, data->file);
    sscanf(buffer,"%s %s %s", &word[0][0], &word[1][0],
	   &word[2][0]);
  } while (!(!strcmp(&word[0][0],"ITER") &&
             (!strcmp(&word[1][0],"EX") || !strcmp(&word[1][0],"TOTAL"))));

  
  /* store current file position since we first have to count iterations */
  filepos = ftell(data->file);

  /* read until the next blank line count the iterations */
  do {
    GET_LINE(buffer, data->file);
    line = trimleft(buffer);
    numread = sscanf(line,"%i %i %*i %*f", &dum, &dum);
    if (numread==2) numiter++;
  } while (strlen(line));

  printf("gamessplugin) %i SCF iterations\n", numiter);

  /* go back and read energies */
  fseek(data->file, filepos, SEEK_SET);
  

  /* allocate memory for scfenergy array */
  ts->scfenergies = (double *)calloc(numiter,sizeof(double));
  
  i=0;
  do {
    GET_LINE(buffer, data->file);
    line = trimleft(buffer);
    numread = sscanf(line,"%i %i %*i %*f", &dum, &dum);
    if (numread==2) {
      sscanf(buffer,"%*i %*i %*i %lf", ts->scfenergies+i);
      i++;
    }
  } while (strlen(line));

#if 0
  for (i=0; i<numiter; i++) {
    printf("scfenergies[%i] = %f\n", i, ts->scfenergies[i]);
  }
#endif

  ts->num_scfiter = numiter;
  
  return TRUE;
}


/*********************************************************
 *
 * this function reads the actual wavefunction, which is
 * punched at the end of the log file
 *
 **********************************************************/
static int get_wavefunction(gamessdata *data, qm_timestep_t *ts, qm_wavefunction_t *wf)
{
  float *orb_energies, orben[5];
  float *wave_coeff;
  char buffer[BUFSIZ];
  char word[6][BUFSIZ];
  int num_orbitals = 0;
  int i = 0, num_values = 0;
  long filepos;
  char *line;
  int have_orben = 0;
  int n[5];

  buffer[0] = '\0';
  for (i=0; i<6; i++) word[i][0] = '\0';

  if (wf == NULL) {
    PRINTERR;	    
    return FALSE;
  }

  wf->type = MOLFILE_WAVE_UNKNOWN;
  wf->spin = SPIN_ALPHA;
  wf->exci = 0;
  strncpy(wf->info, "unknown", MOLFILE_BUFSIZ);

  /*
   * Scan for something like this:

          ------------------
          MOLECULAR ORBITALS     <<--- sometimes EIGENVECTORS
          ------------------

                      1          2          3          4          5
                  -11.0297    -0.9121    -0.5205    -0.5205    -0.5205  <<-- orbital energies
                     A          A          A          A          A   
    1  C  1  S    0.991925   0.221431   0.000006  -0.000001   0.000002
    2  C  1  S    0.038356  -0.627585  -0.000021   0.000003  -0.000006
    3  C  1  X    0.000000  -0.000004   0.338169  -0.030481  -0.460283
     ...

                     6          7          8          9
                    0.7192     0.7192     0.7193     0.7611
                     A          A          A          A   
    1  C  1  S    0.000028   0.000012   0.000092   0.252320
    2  C  1  S   -0.000183  -0.000077  -0.000594  -1.632834
    3  C  1  X   -0.890147   0.062618   0.654017  -0.000154
      ...


     ----------------------------------------------------------------
     PROPERTY VALUES FOR THE RHF   SELF-CONSISTENT FIELD WAVEFUNCTION
     ----------------------------------------------------------------
  * 
  */

  /* Remember position in order to go back if no wave function was found */
  filepos = ftell(data->file);

  do {
    GET_LINE(buffer, data->file);

    line = trimleft(buffer);
    if      (strstr(line, "----- ALPHA SET -----")) {
      wf->type = MOLFILE_WAVE_CANON;
      strncpy(wf->info, "canonical", MOLFILE_BUFSIZ);
      pass_keyline(data->file, "EIGENVECTORS", NULL);
    }
    else if (strstr(line, "----- BETA SET -----")) {
      wf->type = MOLFILE_WAVE_CANON;
      wf->spin = SPIN_BETA;
      strncpy(wf->info, "canonical", MOLFILE_BUFSIZ);
      pass_keyline(data->file, "EIGENVECTORS", NULL);
    }
    else if (strstr(line, "***** BETA ORBITAL LOCALIZATION *****")) {
      wf->spin = SPIN_BETA;
    }
    else if (strstr(line, "EIGENVECTORS")==line) {
      wf->type = MOLFILE_WAVE_CANON;
      strncpy(wf->info, "canonical", MOLFILE_BUFSIZ);
    }
    else if (strstr(line, "MOLECULAR ORBITALS") &&
             !strstr(line, "LOCALIZED")) {
      wf->type = MOLFILE_WAVE_CANON;
      strncpy(wf->info, "canonical", MOLFILE_BUFSIZ);
    }
    else if (strstr(line, "LOCALIZED ORBITALS ARE")) {
      if (strstr(line, "BOYS"))   {
        wf->type = MOLFILE_WAVE_BOYS;
        strncpy(wf->info, "Boys localized", MOLFILE_BUFSIZ);
      }
      else if (strstr(line, "RUEDEN")) {
        wf->type = MOLFILE_WAVE_RUEDEN;
        strncpy(wf->info, "Ruedenberg localized", MOLFILE_BUFSIZ);
      }
      else if (strstr(line, "PIPEK-MEZEY")) {
        wf->type = MOLFILE_WAVE_PIPEK;
        strncpy(wf->info, "Pipek-Mezey localized", MOLFILE_BUFSIZ);
      }
    }
    else if (strstr(line, "GI ORBITALS")) {
      wf->type = MOLFILE_WAVE_GEMINAL;
      strncpy(wf->info, "GVB geminal pairs", MOLFILE_BUFSIZ);
    }
    else if (strstr(line, "MCSCF NATURAL ORBITALS")) {
      wf->type = MOLFILE_WAVE_MCSCFNAT;
      strncpy(wf->info, "MCSCF natural orbitals", MOLFILE_BUFSIZ);
    }
    else if (strstr(line, "MCSCF OPTIMIZED ORBITALS")) {
      wf->type = MOLFILE_WAVE_MCSCFOPT;
      strncpy(wf->info, "MCSCF optimized orbitals", MOLFILE_BUFSIZ);
    }
    else if (strstr(line, "NATURAL ORBITALS IN ATOMIC")==line) {
      wf->type = MOLFILE_WAVE_CINATUR;
      strncpy(wf->info, "CI natural orbitals", MOLFILE_BUFSIZ);
    }
  } while(wf->type==MOLFILE_WAVE_UNKNOWN &&
          !strstr(line, "ENERGY COMPONENTS") &&
          !strstr(line, "***** EQUILIBRIUM GEOMETRY LOCATED") &&
          !strstr(line, "**** THE GEOMETRY SEARCH IS NOT CONVERGED!"));


  /* If we reach the last line of the rhf section without finding 
   * one of the keywords marking the beginning of a wavefunction
   * table then we return.*/
  if (wf->type==MOLFILE_WAVE_UNKNOWN) {
    fseek(data->file, filepos, SEEK_SET);
    return FALSE;
  }

  /* Reserve space for arrays storing wavefunction and orbital
   * energies. we reserve the max. space (num_orbitals==wavef_size)
   * for now and realloc later if, we have less orbitals. */
  wave_coeff = (float *)calloc(data->wavef_size*data->wavef_size,
                                  sizeof(float)); 

  if (wave_coeff == NULL) {
    PRINTERR;	    
    return FALSE;
  }
  
  orb_energies = (float *)calloc(data->wavef_size, sizeof(float));

  if (orb_energies == NULL) {
    free(wave_coeff);
    PRINTERR; 
    return FALSE;
  }

  /* store the pointers */
  wf->wave_coeffs  = wave_coeff;
  wf->orb_energies = orb_energies;

  /* skip the next line which here is typically "-------" */
  eatline(data->file, 1);

  while (1) {
    int over=0;
    float coeff[5];

    if (wf->type == MOLFILE_WAVE_GEMINAL) {
      /* Skip over "PAIR x" header line */
      pass_keyline(data->file, "PAIR ", NULL);
    }

    eatwhitelines(data->file);

    filepos = ftell(data->file);
    GET_LINE(buffer, data->file);
    num_values = sscanf(buffer, "%d %d %d %d %d",
                          &n[0], &n[1], &n[2], &n[3], &n[4]);

    /* If we didn't find this line then coefficient table
     * is complete. */
    if (!num_values) {
      fseek(data->file, filepos, SEEK_SET);
      break;
    }

    eatwhitelines(data->file);

    /* Read first line of orbital energies */
    filepos = ftell(data->file);
    GET_LINE(buffer, data->file);
    have_orben = sscanf(buffer,"%f %f %f %f %f", &orben[0],
                        &orben[1], &orben[2], &orben[3], &orben[4]);
    /* Make sure this is not the first line containing coeffs */
    i = sscanf(buffer, " 1 %*s 1 %*s %f %f %f %f %f",
               &coeff[0], &coeff[1], &coeff[2], &coeff[3], &coeff[4]);
    if (i==num_values) have_orben = 0;

    if (have_orben) {
      /* store the orbital energies in the appropriate arrays 
       * read them until we encounter an empty string */
      for(i=0; i<num_values; i++) {
        *(orb_energies+i) = orben[i];
      }
      
      /* increase orbital energy pointer */
      orb_energies = orb_energies+5;
    }      
    else {
      /* No orbital energies present, go back one line */
      fseek(data->file, filepos, SEEK_SET);
    }

    num_orbitals += num_values;

    /* Find first line containing coefficients */
    filepos = ftell(data->file);
    while (fgets(buffer, sizeof(buffer), data->file)) {
      trimleft(buffer);
      if (strstr(line, "ENERGY COMPONENTS") ||
          strstr(line, "---") ||
          strstr(line, "...")) {
        over = 1; break;
      }
      i = sscanf(buffer, " 1 %*s 1 %*s %f %f %f %f %f",
             &coeff[0], &coeff[1], &coeff[2], &coeff[3], &coeff[4]);
      if (i==num_values) break;
      filepos = ftell(data->file);
    }
    fseek(data->file, filepos, SEEK_SET);

    if (over) break;


    /* Read the wave function coefficient block for up to 5
     * orbitals per line. */
    read_coeff_block(data->file, data->wavef_size,
                     wave_coeff, data->angular_momentum);


    /* move wavefunction pointer to start of next five orbitals */
    if (wf->type == MOLFILE_WAVE_GEMINAL) {
      wave_coeff = wave_coeff + 2*data->wavef_size;
    } else {
      wave_coeff = wave_coeff + 5*data->wavef_size;
    }
  }


  if (!num_orbitals) return FALSE;

  /* resize the array to the actual number of read orbitals */
  if (data->wavef_size!=num_orbitals) {
    wf->orb_energies = (float *)realloc(wf->orb_energies,
          num_orbitals*sizeof(float));

    wf->wave_coeffs  = (float *)realloc(wf->wave_coeffs, data->wavef_size*
					num_orbitals*sizeof(float)); 
  }


  /* store the number of orbitals read in */
  wf->num_orbitals  = num_orbitals;
  wf->have_energies = TRUE;
  wf->have_occup    = FALSE;

  printf("gamessplugin) Number of orbitals scanned: %d \n",\
             wf->num_orbitals);

  return TRUE;
}


/* Read the wave function coefficient block for up to 5
 * orbitals per line:
 *  1  C  1  S    0.989835   0.155361   0.000000  -0.214258   0.000000
 *  2  C  1  S    0.046228  -0.548915   0.000000   0.645267   0.000000
 *  3  C  1  X    0.000000   0.000000   1.030974   0.000000   0.000000
 */
static int read_coeff_block(FILE *file, int wavef_size,
                            float *wave_coeff, int *angular_momentum) {
  int i, j;
  int truncated = 0;
  char buffer[BUFSIZ];
  /* read in the wavefunction */
  /* XXX This will fail for truncated files */
  for (i=0; i<wavef_size; i++) {
    int xexp=0, yexp=0, zexp=0;
    char symm[BUFSIZ];
    float coeff[5];
    int num_values = 0;
    
    GET_LINE(buffer, file);
    
    /* read in the wavefunction coefficients for 5
     * orbitals at a time line by line */
    num_values = sscanf(buffer,"%*5i%*4s%*2i%4s %f %f %f %f %f", 
                        symm, &coeff[0], &coeff[1], &coeff[2],
                        &coeff[3], &coeff[4]);
    
    if (num_values==0) {
      /* The file must have been truncated! */
      truncated = 1;
      break;
    }
    
    for (j=0; j<strlen(symm); j++) {
      switch (symm[j]) {
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
    angular_momentum[3*i  ] = xexp;
    angular_momentum[3*i+1] = yexp;
    angular_momentum[3*i+2] = zexp;
    
    /* each orbital has data->wavef_size entries, 
     * hence we have to use this number as offset when storing 
     * them in groups of five */
    for (j=0 ; j<num_values-1; j++) {
      wave_coeff[j*wavef_size+i] = coeff[j];
    }
  }
  
  if (truncated) return 0;
  
  return 1;
}

/* Read the population analysis section.
 * Currently we parse only the Mulliken charges
 * but we might want to add support for populations
 * and for Lowdin analysis. */
static int get_population(gamessdata *data, qm_timestep_t *ts) {
  int i;
  char buffer[BUFSIZ];
  data->have_mulliken = FALSE;

  if (pass_keyline(data->file,
                     "TOTAL MULLIKEN AND LOWDIN ATOMIC POPULATIONS",
                     "NSERCH=") != FOUND) {
    return FALSE;
  }

  /* Read Mulliken charges if present */
  ts->mulliken_charges = 
    (double *)calloc(data->numatoms, sizeof(double));

  if (!ts->mulliken_charges) {
    PRINTERR; 
    return FALSE;
  }
  
  eatline(data->file, 1);
  
  for (i=0; i<data->numatoms; i++) {
    int n;
    float mullpop, mullcharge, lowpop, lowcharge;
    GET_LINE(buffer, data->file);
    n = sscanf(buffer,"%*i %*s %f %f %f %f",
                   &mullpop, &mullcharge, &lowpop, &lowcharge);
    if (n!=4) return FALSE;
    ts->mulliken_charges[i] = mullcharge;
  }

  if (i!=data->numatoms) return FALSE;

  data->have_mulliken = TRUE;
  return TRUE;
}


/* Read ESP charges.
 * XXX Right now we don't distinguish between different type of
 * ESP-style charges (CHELPG, CONNOLLY, GEODESIC). 
 * This could be solved by reading in the PTSEL keyword in
 * the $PDC group. */
static int get_esp_charges(gamessdata *data) {
  int i;
  char buffer[BUFSIZ];
  data->have_esp = FALSE;

  if (pass_keyline(data->file,
           "ATOM                CHARGE    E.S.D.",
           "...... END OF PROPERTY EVALUATION ") != FOUND) {
    return FALSE;
  }


  /* Read ESP charges if present */
  data->esp_charges = 
    (double *)calloc(data->numatoms, sizeof(double));

  if (data->esp_charges == NULL) {
    PRINTERR; 
    return FALSE;
  }

  eatline(data->file, 1);

  for (i=0; i<data->numatoms; i++) {
    int n;
    double charge;
    GET_LINE(buffer, data->file);
    n = sscanf(buffer,"%*s %lf ", &charge);
    if (n!=2) return FALSE; /* XXX n!=1 ?? */
    data->esp_charges[i] = charge;
  }

  if (i!=data->numatoms) return FALSE;

  data->have_esp = TRUE;
  return TRUE;
}


static int get_gradient(gamessdata *data, qm_timestep_t *ts) {
  char buffer[BUFSIZ];
  float dx, dy, dz;
  long filepos;
  int i=0;
  int numread;

  buffer[0] = '\0';

  /* remember position in order to go back if no forces were found */
  filepos = ftell(data->file);

  /* look for GRADIENT section */
  if (goto_keystring2(data->file, "GRADIENT (HARTREE",
                "***** EQUILIBRIUM GEOMETRY LOCATED", 
                " BEGINNING GEOMETRY SEARCH") != FOUND) {
    fseek(data->file, filepos, SEEK_SET);
    return FALSE;
  }

  eatline(data->file, 3);

  ts->gradient = (float *)calloc(3*data->numatoms, sizeof(float));

  if (ts->gradient == NULL) {
    PRINTERR;	    
    fseek(data->file, filepos, SEEK_SET);
    return FALSE;
  }

  /* read the gradient table */
  do {
    GET_LINE(buffer, data->file);
    numread = sscanf(buffer, "%i %*s %*f %f %f %f", &i, &dx, &dy, &dz);
    if (numread==4) {
      ts->gradient[3*(i-1)  ] = dx;
      ts->gradient[3*(i-1)+1] = dy;
      ts->gradient[3*(i-1)+2] = dz;
    }

  } while(numread==4);

  if (i!=data->numatoms) {
    printf("gamessplugin) Number of gradients != number of atoms!\n");
    fseek(data->file, filepos, SEEK_SET);
    return FALSE;
  }

  return TRUE;
}



/***********************************************************
 *
 * this function reads in the wavenumbers and intensities of
 * the normal modes
 *
 * *********************************************************/
static int get_normal_modes(gamessdata *data) {

  char word[4][BUFSIZ];
  char buffer[BUFSIZ];
  int i = 0, k = 0, j = 0;
  int remaining_columns;
  double entry[6]; 
  char separator;
  char *item;
  int counter;


  /* initialize arrays */
  buffer[0] = '\0';
  memset(entry, 0, sizeof(entry));
  for (i=0; i<4; i++) word[i][0] = '\0';

    
  /* reserve memory for dynamically allocated data
   * arrays */
  item = (char *)calloc(BUFSIZ,sizeof(char));

  if (item==NULL){
    PRINTERR;
    return FALSE;
  }


  data->wavenumbers = 
    (float *)calloc(data->numatoms*3,sizeof(float));

  if (data->wavenumbers==NULL) {
    PRINTERR;
    return FALSE;
  }


  data->intensities = 
    (float *)calloc(data->numatoms*3,sizeof(float));

  if (data->intensities==NULL) {
    PRINTERR;
    return FALSE;
  }


  data->normal_modes = 
    (float *)calloc((data->numatoms*3)*(data->numatoms*3),
		     sizeof(float));

  if (data->normal_modes==NULL) {
    PRINTERR;
    return FALSE;
  }


  /* look for FREQUENCYs and IR INTENSITIES */
  for (i=0; i<data->numatoms*3/5; i++) {
    do {
      GET_LINE(buffer, data->file);
      sscanf(buffer,"%s ",&word[0][0]);

    } while(strcmp(&word[0][0],"FREQUENCY:"));


    /* scan the frequencies; 
     * this requires some care since there might
     * be imaginary modes present which leads to
     * an additional char per imaginary mode per
     * line */

    /* check for imaginary modes */
    if (strchr(buffer,'I') != NULL) {
      separator = ' ';
      counter = 0;

      /* read all line elements into individual strings */

      /* skip first entry "FREQUENCY" */
      item = strtok(buffer,&separator);
      
      /* start going through the string */
      while ((item = strtok(NULL,&separator)) != NULL) {
        /* check if item is 'I'; if yes, mark previous mode
         * as imaginary; otherwise save the mode */
        if (*item=='I') data->nimag++;

        else {
          /* save only the first 5 modes - there NEVER should
           * be more in any case, but just to make sure
           * we don't overrun the array */
          if (counter<5) {
            *(data->wavenumbers+(i*5)+counter) = atof(item);
            counter++;
          }
        }
      }
    }

    /* no imaginary mode, reading is straightforward */
    else {
      sscanf(buffer,"%s %lf %lf %lf %lf %lf",&word[0][0],
             &entry[0],&entry[1],&entry[2],&entry[3],&entry[4]); 
    
      for (k=0; k<5; k++) {
        *(data->wavenumbers+(i*5)+k) = entry[k]; 
      }
    }

    eatline(data->file, 1);

    /* next line contains the IR INTENSITIES */
    GET_LINE(buffer, data->file);

    /* scan the IR INTENSITIES */
    sscanf(buffer,"%s %s %lf %lf %lf %lf %lf",&word[0][0],&word[1][0],
           &entry[0],&entry[1],&entry[2],&entry[3],&entry[4]);
 
    for (k=0; k<5; k++) {
      *(data->intensities+(i*5)+k) = entry[k]; 
    }

    eatline(data->file, 1);

    /* read the following five modes */
    for (k=0; k<data->numatoms; k++) {
      /* x */
      GET_LINE(buffer, data->file);
      sscanf(buffer,"%s %s %s %lf %lf %lf %lf %lf",&word[0][0], 
             &word[1][0], &word[2][0], &entry[0], &entry[1], &entry[2],
             &entry[3], &entry[4]);

      for (j=0; j<5; j++) {
        *(data->normal_modes+(3*k)+((i*5+j)*3*data->numatoms)) = 
          entry[j];
      }

      /* y */
      GET_LINE(buffer, data->file);
      sscanf(buffer,"%s %lf %lf %lf %lf %lf",&word[0][0], &entry[0],
             &entry[1],&entry[2], &entry[3],&entry[4]);

      for (j=0; j<5; j++) {
        *(data->normal_modes+(3*k+1)+((i*5+j)*3*data->numatoms)) = 
          entry[j];
      }

      /* z */
      GET_LINE(buffer, data->file);
      sscanf(buffer,"%s %lf %lf %lf %lf %lf",&word[0][0], &entry[0],
             &entry[1], &entry[2], &entry[3],&entry[4]);

      for (j=0; j<5; j++) {
        *(data->normal_modes+(3*k+2)+((i*5+j)*3*data->numatoms)) = 
          entry[j];
      }
    }
  }


  /* read the remaining columns */
  if ( (remaining_columns=(data->numatoms*3)%5) != 0 ) { 

    /* move to next set of modes */
    do {
      GET_LINE(buffer, data->file);
      sscanf(buffer,"%s ",&word[0][0]);

    } while(strcmp(&word[0][0],"FREQUENCY:"));


    /* scan the frequencies; 
     * this requires some care since there might
     * be imaginary modes present which leads to
     * an additional char per imaginary mode per
     * line */

    /* check for imaginary modes */
    if ( strchr(buffer,'I') != NULL ) {
      /* initialize */
      separator = ' ';
      counter = 0;

      /* read all line elements into individual strings */

      /* skip first entry "FREQUENCY" */
      item = strtok(buffer,&separator);

      
      /* start going through the string */
      while ( (item = strtok(NULL,&separator)) != NULL ) 
      {
        /* check if item is 'I'; if yes, mark previous mode
         * as imaginary; otherwise save the mode */
        if ( *item == 'I' ) data->nimag++;

        else {
          /* save only the first remaining columns modes - 
           * there NEVER should
           * be more in any case, but just to make sure
           * we don't overrun the array */
          if (counter<remaining_columns) {
            *(data->wavenumbers+(i*5)+counter) = atof(item);
            counter++;
          }
        }
      } 
    }

    /* no imaginary mode, reading is straightforward */
    else {
      sscanf(buffer,"%s %lf %lf %lf %lf %lf",&word[0][0],
             &entry[0],&entry[1],&entry[2],&entry[3],&entry[4]); 
 
      for (k=0; k<remaining_columns; k++) {
        *(data->wavenumbers+(i*5)+k) = entry[k]; 
      }
    }

    eatline(data->file, 1);
    
    /* next line contains the IR INTENSITIES */
    GET_LINE(buffer, data->file);

    /* scan the IR INTENSITIES */
    sscanf(buffer,"%s %s %lf %lf %lf %lf %lf",&word[0][0],&word[1][0],
           &entry[0],&entry[1],&entry[2],&entry[3],&entry[4]);
 
    for (k=0; k<remaining_columns; k++) {
      *(data->intensities+(i*5)+k) = entry[k];
    }

    eatline(data->file, 1);

    /* read the following five modes */
    for (k=0; k<data->numatoms; k++) {
      /* x */
      GET_LINE(buffer, data->file);
      sscanf(buffer,"%s %s %s %lf %lf %lf %lf %lf",&word[0][0], 
             &word[1][0], &word[2][0], &entry[0], &entry[1], &entry[2],
             &entry[3], &entry[4]);

      for (j=0; j<remaining_columns; j++) {
        *(data->normal_modes+(3*k)+((i*5+j)*3*data->numatoms)) = 
          entry[j];
      }
      
      /* y */
      GET_LINE(buffer, data->file);
      sscanf(buffer,"%s %lf %lf %lf %lf %lf",&word[0][0], &entry[0],
             &entry[1],&entry[2], &entry[3],&entry[4]);
      
      for (j=0; j<remaining_columns; j++) {
        *(data->normal_modes+(3*k+1)+((i*5+j)*3*data->numatoms)) = 
          entry[j];
      }

      /* z */
      GET_LINE(buffer, data->file);
      sscanf(buffer,"%s %lf %lf %lf %lf %lf",&word[0][0], &entry[0],
             &entry[1], &entry[2], &entry[3],&entry[4]);
      
      for (j=0; j<remaining_columns; j++) {
        *(data->normal_modes+(3*k+2)+((i*5+j)*3*data->numatoms)) = 
          entry[j];
      }
    }
  }

  data->have_normal_modes = TRUE;

  /* release memory that is not needed any more */
  free(item);

  /* print brief message */
  printf("gamessplugin) Successfully scanned normal modes \n");

  return TRUE;
}



/***********************************************************
 *
 * this function reads in the cartesian hessian matrix 
 *
 * *********************************************************/
static int get_cart_hessian(gamessdata *data)
{
  char word[4][BUFSIZ];
  char buffer[BUFSIZ];
  char dummy; 
  int i,j,k;
  float entry[6]; 

  buffer[0] = '\0';
  memset(entry, 0, sizeof(entry));
  for (i=0; i<4; i++) word[i][0] = '\0';


  /* at this point we need to rewind the file, since
   * in case that there is no internal Hessian stuff the
   * previous call to get_int_coords scanned the file
   * until EOF */
  rewind(data->file);


  /* look for CARTESIAN FORCE CONSTANT MATRIX */
  do {
    GET_LINE(buffer, data->file);
    sscanf(buffer,"%s %s %s %s",&word[0][0],&word[1][0],
           &word[2][0],&word[3][0]);

  } while(strcmp(&word[0][0],"CARTESIAN") || 
          strcmp(&word[1][0],"FORCE") ||
          strcmp(&word[2][0],"CONSTANT") ||
          strcmp(&word[3][0],"MATRIX"));


  /* skip next 5 lines */
  for (i=0; i<5; i++) eatline(data->file, 1);


  /* reserve memory for array; 
   * NOTE: this is a lower triangular matrix, but for now
   * we save it in an square matrix of dim(3Nx3N) to 
   * facilitate element access */
  data->carthessian = 
    (double *)calloc((data->numatoms*3)*(data->numatoms*3),
		     sizeof(double));

  
  /* make sure memory was allocated properly */
  if (data->carthessian == NULL) {
    PRINTERR;
    return FALSE;
  }


  /* start scanning; the cartesian hessian matrix is a lower
   * triangular matrix, organized in rows of 6 */

  /* read blocks with complete rows of 6 */
  for (i=0; i<(int)(data->numatoms/2); i++) {
    for (j=0; j<(data->numatoms*3)-(i*6); j++) {
      GET_LINE(buffer, data->file);
 
      if (j%3==0) {
        sscanf(buffer,"%s %s %c %f %f %f %f %f %f",
               &word[0][0],&word[1][0],&dummy,&entry[0],&entry[1],
               &entry[2],&entry[3],&entry[4],&entry[5]);
      }
      else {
        sscanf(buffer,"%1s %f %f %f %f %f %f",
               &dummy,&entry[0],&entry[1],&entry[2],&entry[3],&entry[4],
               &entry[5]);
      }


      /* save entries (lower triangular matrix) in a 
       * square matrix */
      for (k=0; k<=(j<5 ? j : 5); k++) {
        *(data->carthessian+((j+(i*6))*3*data->numatoms)+
          (k+(i*6))) = entry[k];
      }
    }

    /* skip the three line separating the matrix entries */
    eatline(data->file, 4);
  }

  printf("gamessplugin) Scanned Hessian in CARTESIAN coordinates\n");

  data->have_cart_hessian = TRUE;

  return TRUE;
}
  
  
  
/***********************************************************
 *
 * this function reads in the hessian in internal coordinates
 * as well as the internal coordinates 
 *
 * *********************************************************/
static int get_int_coords(gamessdata *data) {

  char word[5][BUFSIZ];
  char buffer[BUFSIZ];
  long position;
  int first, second, third, fourth;
  double value;
  double hess[5];
  int i = 0, j = 0, k = 0, l = 0;
  int n, dummy, remaining_blocks;

  buffer[0] = '\0';
  memset(hess, 0, sizeof(hess));
  for (i=0; i<5; i++) word[i][0] = '\0';

  /* look for list of INTERNAL COORDINATES */
  do {
    GET_LINE(buffer, data->file);
    sscanf(buffer,"%s %s ",&word[0][0],&word[1][0]);

  } while(strcmp(&word[0][0],"INTERNAL") || 
          strcmp(&word[1][0],"COORDINATES"));

  
  /* skip next 5 lines */
  for (i=0; i<5; i++) eatline(data->file, 1);

  /* remember current position so we can jump back */
  position = ftell(data->file);

  /* scan the next line */
  GET_LINE(buffer, data->file);
  n = sscanf(buffer,"%s %s", &word[0][0], &word[1][0]); 

  /* read line by line */
  while (n!=-1) {
    /* start counting the number of internal coordinates */
    data->nintcoords++;

    /* count the number of bonds, angles, dihedrals */
    if (!strcmp(&word[1][0],"STRETCH")) {
      data->nbonds++;
    }
    else if (!strcmp(&word[1][0],"BEND")) {
      data->nangles++;
    }
    else if (!strcmp(&word[1][0],"TORSION")) {
      data->ndiheds++;
    }
    else if (!strcmp(&word[1][0],"PLA.BEND")) {
      data->nimprops++;
    }

    /* scan next line */
    GET_LINE(buffer, data->file);
    n = sscanf(buffer,"%s %s", &word[0][0], &word[1][0]); 
  }

  /* now that we know the number of bonds, angles, etc.
   * we can read and store the internal coordinates */
  fseek(data->file,position,SEEK_SET);


  /* reserve memory for the arrays storing the internal
   * coordinates and their values */
  data->bonds = (int *)calloc(2*data->nbonds,sizeof(int));
  data->angles = (int *)calloc(3*data->nangles,sizeof(int));
  data->dihedrals = (int *)calloc(4*data->ndiheds,sizeof(int));
  data->impropers = (int *)calloc(4*data->nimprops,sizeof(int));
  data->internal_coordinates = (double *)calloc(data->nintcoords,
	sizeof(double));


  /* check if we have sufficient memory available */
  if ( (data->bonds == NULL) || 
       (data->angles == NULL) ||
       (data->dihedrals == NULL) || 
       (data->internal_coordinates == NULL)) 
  {
    PRINTERR; 
    return FALSE;
  }


  /* now start going through the internal coordinates
   * and save them in the appropriate arrays; here
   * I drop all safety check since we went through
   * this part of the file already once and should
   * be good */
 
  /* scan the STRETCHES */
  for (i=0; i<data->nbonds; i++) {
    GET_LINE(buffer, data->file);
    sscanf(buffer,"%s %s %d %d %lf", &word[0][0], &word[1][0], 
           &first, &second, &value);

    *(data->bonds+2*i) = first;
    *(data->bonds+2*i+1) = second;
    *(data->internal_coordinates+i) = value;
  }

  /* scan the BENDS */
  for (j=0; j<data->nangles; j++) {
    GET_LINE(buffer, data->file);
    sscanf(buffer,"%s %s %d %d %d %lf", &word[0][0], &word[1][0], 
           &first, &second, &third, &value);

    *(data->angles+3*j) = first;
    *(data->angles+3*j+1) = second;
    *(data->angles+3*j+2) = third;
    *(data->internal_coordinates+i+j) = value;
  }

  /* scan the TORSIONS */
  for (k=0; k<data->ndiheds; k++) {
    GET_LINE(buffer, data->file);
    sscanf(buffer,"%s %s %d %d %d %d %lf", &word[0][0], &word[1][0],
           &first, &second, &third, &fourth, &value);

    *(data->dihedrals+4*k) = first;
    *(data->dihedrals+4*k+1) = second;
    *(data->dihedrals+4*k+2) = third;
    *(data->dihedrals+4*k+3) = fourth;
    *(data->internal_coordinates+i+j+k) = value;
  }

  /* scan the IMPROPERS */
  for (l=0; l<data->nimprops; l++) {
    GET_LINE(buffer, data->file);
    sscanf(buffer,"%s %s %d %d %d %d %lf", &word[0][0], &word[1][0],
           &first, &second, &third, &fourth, &value);

    *(data->impropers+4*l) = first;
    *(data->impropers+4*l+1) = second;
    *(data->impropers+4*l+2) = third;
    *(data->impropers+4*l+3) = fourth;
    *(data->internal_coordinates+i+j+k+l) = value;
  }

  printf("gamessplugin) Scanned %d INTERNAL coordinates \n",
         data->nintcoords);
  printf("gamessplugin)    %d BONDS \n",data->nbonds);
  printf("gamessplugin)    %d ANGLES \n",data->nangles);
  printf("gamessplugin)    %d DIHEDRALS \n",data->ndiheds);
  printf("gamessplugin)    %d IMPROPERS \n",data->nimprops);


  /* next read in the hessian in internal coordinates;
   * we would expect the matrix immediately after the
   * internal coordinates in the output files;
   * we check this first */
  eatline(data->file, 1);


  GET_LINE(buffer, data->file);
  sscanf(buffer,"%s %s %s %s %s", &word[0][0], &word[1][0], 
      &word[2][0], &word[3][0], &word[4][0]);

  if ( strcmp(&word[0][0],"HESSIAN") || 
       strcmp(&word[1][0],"MATRIX") ||
       strcmp(&word[3][0],"INTERNAL") || 
       strcmp(&word[4][0],"COORDINATES")) 
  {
    /* apparently we are out of luck - no Hessian in internal
     * coordinates -- GOOD BYE :) */
    return FALSE;
  }
 

  /* skip to the hessian arrays */
  while (sscanf(buffer,"%d %lf %lf %lf %lf %lf", 
                &dummy, &hess[0], &hess[1], &hess[2], 
                &hess[3], &hess[4]) != 6) {
    GET_LINE(buffer, data->file);
  }

  
  /* reserve memory for inthessian array */
  data->inthessian = 
    (double *)calloc((data->nintcoords)*(data->nintcoords),
		     sizeof(double));


  /* make sure memory was allocated properly */
  if (data->inthessian == NULL) {
    PRINTERR;
    return FALSE;
  }


  /* start scanning; GAMESS organized the output of the
   * internal HESSIAN in rows of 5 */

  /* read blocks with complete rows of 5 */
  for (i=0; i<(int)(data->nintcoords/5); i++) {
    for (j=0; j<data->nintcoords; j++) {
      sscanf(buffer,"%d %lf %lf %lf %lf %lf", &dummy, &hess[0], 
             &hess[1], &hess[2], &hess[3], &hess[4]);

      /* save entries */
      for (k=0; k<5; k++) { 
        *(data->inthessian+(j*data->nintcoords)+(i*5)+k) = hess[k];
      }

      /* next line */
      GET_LINE(buffer, data->file);
    }

    /* skip the two lines separating the matrix entries 
     * and scan next line */
    eatline(data->file, 2);

    GET_LINE(buffer, data->file);
  }

  
  /* read the remaining block with less then 5 rows
   * if present */
  remaining_blocks = data->nintcoords%5;
  
  if (remaining_blocks!=0) {
    for (j=0; j<data->nintcoords; j++) {
      sscanf(buffer,"%d %lf %lf %lf %lf %lf", &dummy, &hess[0], 
             &hess[1], &hess[2], &hess[3], &hess[4]);

      for (k=0; k<remaining_blocks; k++) { 
        *(data->inthessian+(j*data->nintcoords)+(i*5)+k) = hess[k];
      }

      GET_LINE(buffer, data->file);
    }
  }

  printf("gamessplugin) Scanned Hessian in INTERNAL coordinates\n");

  /* finally, dump the diagonal elements of the hessian into the
   * force constant arrays, after converting the units 
   * appropriately;
   * BONDS are in HARTREE/BOHR**2
   * ANGLES,DIHEDRALS,IMPROPERS are in HARTREE/RADIAN**2 */
  
  /* allocate dynamic arrays */
  data->bond_force_const = 
    (double *)calloc(data->nbonds,sizeof(double));
 
  if (data->bond_force_const==NULL) {
    PRINTERR;
    return FALSE;
  }


  data->angle_force_const =
    (double *)calloc(data->nangles,sizeof(double));

  if (data->angle_force_const==NULL) {
    PRINTERR;
    return FALSE;
  }


  data->dihedral_force_const =
    (double *)calloc(data->ndiheds,sizeof(double));

  if (data->dihedral_force_const==NULL) {
    PRINTERR;
    return FALSE;
  }


  data->improper_force_const =
    (double *)calloc(data->nimprops,sizeof(double));

  if (data->improper_force_const==NULL) {
    PRINTERR;
    return FALSE;
  }

  /* scan the bonds */
  for (i=0; i<data->nbonds; i++) {
    *(data->bond_force_const + i) = 
      *(data->inthessian+(i*data->nintcoords)+i) 
      * HARTREE_TO_KCAL / BOHR_TO_ANGS / BOHR_TO_ANGS;

    printf("%3d (BOND) %2d - %2d : %f (CHARMM) %f \n",i, 
           *(data->bonds+2*i), *(data->bonds+2*i+1),
           *(data->bond_force_const +i),
           *(data->bond_force_const +i)*0.5);
  }
  
  /* scan the angles */
  for (j=i; j<i+(data->nangles); j++) {
    *(data->angle_force_const + (j-i)) = 
      *(data->inthessian+(j*data->nintcoords)+j) 
      * HARTREE_TO_KCAL;
    
    printf("%3d (ANGLE) %2d - %2d - %2d : %f (CHARMM) %f \n",j,
           *(data->angles+3*(j-i)), *(data->angles+3*(j-i)+1), 
           *(data->angles+3*(j-i)+2), 
           *(data->angle_force_const + (j-i)),
           *(data->angle_force_const + (j-i))*0.5);
  }

  /* scan the dihedrals */
  for (k=j; k<j+(data->ndiheds); k++) {
    *(data->dihedral_force_const + (k-j)) = 
      *(data->inthessian+(k*data->nintcoords)+k)
      * HARTREE_TO_KCAL;
    
    printf("%3d (DIHEDRAL) %2d - %2d - %2d - %2d : %f \n",k,
           *(data->dihedrals+4*(k-j)), *(data->dihedrals+4*(k-j)+1),
           *(data->dihedrals+4*(k-j)+2), *(data->dihedrals+4*(k-j)+3),
           *(data->dihedral_force_const + (k-j))); 
  }

  /* scan the impropers */
  for (l=k; l<k+(data->nimprops); l++) {
    *(data->improper_force_const + (l-k)) = 
      *(data->inthessian+(l*data->nintcoords)+l)
      * HARTREE_TO_KCAL;
    
    printf("%3d (IMPROPERS) %2d - %2d - %2d - %2d : %f \n",l,
           *(data->impropers+4*(l-k)), *(data->impropers+4*(l-k)+1),
           *(data->impropers+4*(l-k)+2), *(data->impropers+4*(l-k)+3),
           *(data->improper_force_const + (l-k)));
  }


  data->have_internals = TRUE;
  return TRUE;
}


#if 0

/************************************************************
 *
 * this function animates a given normal mode by means of
 * generating mod_num_frames frames away from the equilibrium
 * structure in a direction given by the hessiane 
 *
 ************************************************************/
static int animate_normal_mode(gamessdata *data, int mode)
{
  mode_data *animated_mode = data->animated_mode;
  float *normal_modes = data->normal_modes;
  float scale = animated_mode->mode_scaling;
  int i = 0, k = 0; 
  int l = 0, m = 0;
  int natoms = data->numatoms;
  int num_frames = animated_mode->mode_num_frames;

  /* first sweep to max of interval */
  for ( k = 0; k < num_frames+1; ++k)
  {
    for ( i = 0; i < natoms; ++i)
    {
      *(animated_mode->mode_frames+(k*natoms*3)+(3*i)) = 
	  (data->initatoms+i)->x * (1+( k*scale * 
	  (*(normal_modes+(mode*natoms*3)+(3*i)))));

      *(animated_mode->mode_frames+(k*natoms*3)+(3*i+1)) = 
	  (data->initatoms+i)->y * (1+( k*scale * 
	  (*(normal_modes+(mode*natoms*3)+(3*i+1)))));

      *(animated_mode->mode_frames+(k*natoms*3)+(3*i+2)) = 
	  (data->initatoms+i)->z * (1+( k*scale * 
	  (*(normal_modes+(mode*natoms*3)+(3*i+2)))));
    }
  }


  /* second sweep all the way back to min of interval */
  for ( l = 0; l < 2*num_frames+1; ++l)
  {
    for ( i = 0; i < natoms; ++i)
    {
      *(animated_mode->mode_frames+((l+k)*natoms*3)+(3*i)) = 
	  (data->initatoms+i)->x * (1+((int)(num_frames-l)*scale * 
	  (*(normal_modes+(mode*natoms*3)+(3*i)))));

      *(animated_mode->mode_frames+((l+k)*natoms*3)+(3*i+1)) = 
	  (data->initatoms+i)->y * (1+((int)(num_frames-l)*scale * 
	  (*(normal_modes+(mode*natoms*3)+(3*i+1)))));

      *(animated_mode->mode_frames+((l+k)*natoms*3)+(3*i+2)) = 
	  (data->initatoms+i)->z * (1+((int)(num_frames-l)*scale * 
	  (*(normal_modes+(mode*natoms*3)+(3*i+2)))));
    }
  }


  /* third sweep back to the native starting structure */
  for ( m = 0; m < num_frames+1; ++m)
  {
    for ( i = 0; i < natoms; ++i)
    {
      *(animated_mode->mode_frames+((l+k+m)*natoms*3)+(3*i)) = 
	  (data->initatoms+i)->x * (1+((int)(m-num_frames)*scale * 
	  (*(normal_modes+(mode*natoms*3)+(3*i)))));

      *(animated_mode->mode_frames+((l+k+m)*natoms*3)+(3*i+1)) = 
	  (data->initatoms+i)->y * (1+((int)(m-num_frames)*scale * 
	  (*(normal_modes+(mode*natoms*3)+(3*i+1)))));

      *(animated_mode->mode_frames+((l+k+m)*natoms*3)+(3*i+2)) = 
	  (data->initatoms+i)->z * (1+((int)(m-num_frames)*scale * 
	  (*(normal_modes+(mode*natoms*3)+(3*i+2)))));
    }
  }

  printf("gamessplugin) Successfully animated mode %d \n", mode);

  return TRUE;
}
#endif


/********************************************************
 *
 * this function removes trailing spaces/newlines off
 * a character array
 *
 * FIXME: This function is dangerous; if the string is 
 * screwed up and doesn't contain any char that the while
 * loop is looking for we get a nice segfault :(
 ********************************************************/
static char* chop_string_all(char* the_string)
{
  int i = 0;

  while (the_string[i]!='\n' && the_string[i]!=' ' && 
         the_string[i]!='\0') {
    i++;
  }

  the_string[i] = '\0';

  return the_string;
}


/********************************************************
 *
 * this function returns a pointer to the first non-whitespace
 * character in a string.
 * The c-string must be null-terminated.
 *
 ********************************************************/
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

static char* trimright(char* s)
{
  int i;
  for (i=strlen(s)-1; i>=0; i--) {
    if (!isspace(s[i])) break;
  }
  s[i+1] = '\0';

  return s;
}

/* Advances the file pointer until the first appearance
 * of a line beginning with the given keystring. Leading
 * whitespace in the lines are ignored, the keystring 
 * should not begin with whitespace otherwise no match
 * will be found.
 * If stopstring is encountered before the keystring 
 * the file is rewound to the position where the search
 * started. If stopstring is NULL then the search stops
 * at EOF. */
static int goto_keyline(FILE *file, const char *keystring,
        const char *stopstring) {
  char buffer[BUFSIZ];
  char *line;
  int found = 0;
  long filepos, curline;
  filepos = ftell(file);

  do {
    curline = ftell(file);
    if (!fgets(buffer, sizeof(buffer), file)) {
      fseek(file, filepos, SEEK_SET);
      return 0;
    }
    line = trimleft(buffer);
    if (strstr(line, keystring)) {
      found = 1;
      fseek(file, curline, SEEK_SET);
      break;
    }
  } while (!stopstring || !strstr(line, stopstring));
    
  if (!found) {
    fseek(file, filepos, SEEK_SET);
    return STOPPED;
  }

  return FOUND;
}

/* places file pointer AFTER the line containing the string */
static int pass_keyline(FILE *file, const char *keystring,
        const char *stopstring) {
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
  } while (!stopstring || !strstr(line, stopstring));
    
  if (!found) {
    fseek(file, filepos, SEEK_SET);
    return STOPPED;
  }

  return FOUND;
}

static int goto_keystring2(FILE *file, const char *keystring,
        const char *stopstring1, const char *stopstring2) {
  char buffer[BUFSIZ];
  char *line;
  int found = 0;
  long filepos;
  filepos = ftell(file);

  do {
    if (!fgets(buffer, sizeof(buffer), file)) break;
    line = trimleft(buffer);
    if (strstr(line, keystring)) {
      found = 1;
      break;
    }
  } while (!stopstring1 || !strstr(line, stopstring1) || 
           !stopstring2 || !strstr(line, stopstring2));
    
  if (!found) {
    fseek(file, filepos, SEEK_SET);
    return 0;
  }

  return 1;
}

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
  } while (!strlen(line));

  printf("HERE) %s\n", buffer);
  fseek(file, filepos, SEEK_SET);
}

/********************************************************
 *
 * this function removes trailing newlines off
 * a character array
 *
 ********************************************************/
static char* chop_string_nl(char* the_string)
{
  int i = 0;

  while (the_string[i]!='\n' && the_string[i]!='\0') i++;

  the_string[i] = '\0';

  return the_string;
}



/*************************************************************
 *
 * plugin registration 
 *
 **************************************************************/
static molfile_plugin_t plugin;

VMDPLUGIN_API int VMDPLUGIN_init(void) {
  memset(&plugin, 0, sizeof(molfile_plugin_t));
  plugin.abiversion = vmdplugin_ABIVERSION;
  plugin.type = MOLFILE_PLUGIN_TYPE;
  plugin.name = "gamess";
  plugin.prettyname = "GAMESS";
  plugin.author = "Markus Dittrich, Jan Saam";
  plugin.majorv = 0;
  plugin.minorv = 11;
  plugin.is_reentrant = VMDPLUGIN_THREADUNSAFE;
  plugin.filename_extension = "log";
  plugin.open_file_read = open_gamess_read;
  plugin.read_structure = read_gamess_structure;
  plugin.close_file_read = close_gamess_read;

  plugin.read_qm_metadata = read_gamess_metadata;
  plugin.read_qm_rundata  = read_gamess_rundata;

#if vmdplugin_ABIVERSION > 11
  plugin.read_timestep_metadata    = read_timestep_metadata;
  plugin.read_qm_timestep_metadata = read_qm_timestep_metadata;
  plugin.read_timestep = read_timestep;
#endif

  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_fini(void) {
  return VMDPLUGIN_SUCCESS;
}
