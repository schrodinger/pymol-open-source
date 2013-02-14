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
 *      $Author: johanstr $       $Locker:  $             $State: Exp $
 *      $Revision: 1.207 $       $Date: 2012/01/24 00:59:06 $
 *
 ***************************************************************************/

/* *******************************************************
 *
 *          G A M E S S     P L U G I N 
 *
 * Read GAMESS-US log files.
 * This file reader is fairly robust and should be able to 
 * parse GAMESS logfiles from a wide range of versions
 * (2000-2009 tested). 
 * From the variety of data that can be found in GAMESS
 * logfiles we read coordinates, basis set, wavefunctions, 
 * gradients, hessian, charges, frequencies and a number 
 * of data describing the type of calculation.
 *
 * ********************************************************/


/**********************************************************

 
FUNCTION CALL CHAIN
===================

Below is an overview about the hierarchy of function calls.
Not all functions or possible branches are listed but the
general picture is covered.

The top level consists of functions defined by the 
molfile_plugin interface: First VMD calls open_gamess_read()
then it requests information about atoms and topology using 
read_gamess_structure(). Note that no coordinates are provided
at this step. Timestep independent data such as info about the
calculation method or the basis set are provided through 
read_gamess_rundata(). Since the allocation of these
arrays populated for this purpose is done by VMD rather than
the plugin, VMD needs to know the sizes beforehand. This is 
achieved by calling read_qm_metadata() before read_qm_rundata().

Next, in order to obtain the info for all timesteps VMD will call
read_timestep_metadata(), read_qm_timestep_metadata(), and
read_timestep() repeatedly until read_timestep() returns
MOLFILE_ERROR indicating the end of the trajectory.
Here too, the metadata are transferred before the main chunk
of data.

Finally, VMD calls close_gamess_read() which frees the temporary
memory used by the plugin and closes the file.


VMD
 |__ open_gamess_read()
 |    |__ parse_static_data()
 |        |__ get_proc_mem()
 |        |__ get_basis_options()
 |        |__ get_runtitle()
 |        |__ get_contrl()
 |        |__ get_input_structure()
 |        |    |__ get_coordinates()
 |        |
 |        |__ get_basis()
 |        |    |__ read_shell_primitives()
 |        |    |__ fill_basis_arrays()
 |        |
 |        |__ get_basis_stats()
 |        |__ get_properties_input()
 |        |__ get_int_coords()
 |        |__ get_symmetry()
 |        |__ get_guess_options()
 |        |__ get_mcscf()
 |        |__ analyze_traj()
 |        |__ read_first_frame()
 |        |    |__get_traj_frame()
 |        |        |__get_coordinates()
 |        |        |__get_scfdata()
 |        |        |__check_add_wavefunctions()
 |        |        |   |__add_wavefunction()
 |        |        |   |__get_wavefunction()
 |        |        |   |   |__read_coeff_block()
 |        |        |   |       |__angular_momentum_expon()
 |        |        |   |
 |        |        |   |__del_wavefunction()
 |        |        |
 |        |        |__get_population()
 |        |        |__get_gradient()
 |        |
 |        |__ get_final_properties()
 |             |__ get_population()
 |             |__ get_esp_charges()
 |             |__ get_final_gradient()
 |             |__ get_int_hessian()
 |             |__ get_cart_hessian()
 |             |__ get_normal_modes()
 |             |__ read_localized_orbitals()
 |        
 |
 |__ read_gamess_structure()
 |
 |__ read_gamess_metadata()
 |
 |__ read_gamess_rundata()
 |
 |
 |  DO:                         <----.
 |__ read_timestep_metadata()        |
 |__ read_qm_timestep_metadata()     |
 |__ read_timestep()                 |
 |                                   |
 |  WHILE(FRAMES)               -----' 
 |
 |
 |__ close_gamess_read()

 

PARSING STRATEGY
================

Because we potentially have to read quite a bit into the
logfile in order to obtain the number of atoms required by 
open_gamess_read(), we just parse the whole file and store 
all timestep independent data. This process is managed by 
parse_static_data(). Some of the static data are at the end
of the file. Function analyze_traj() find the end of the
trajectory and records the file pointer for the beginning
of each frame. Thus, reading the frames when they are requested 
later by read_timestep() is much faster without having to  
keep large amountss of data in memory.


**********************************************************/

#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <time.h>

#include "qmplugin.h"
#include "unit_conversion.h"
#include "periodic_table.h"

#define ANGSTROM 0
#define BOHR     1


/* #define DEBUGGING 1 */

/*
 * Error reporting macro for use in DEBUG mode
 */

#ifdef DEBUGGING
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


#define NOTFOUND 0
#define FOUND    1
#define STOPPED  2


/* Data specific to parsing GAMESS files */
typedef struct 
{
  int  version; /* here we track the GAMESS versions, since the
                 * file format has changed with 
                 * version 27 JUN 2005 (R2);
                 * version = 1  : pre-27 JUN 2005 (R2)
                 * version = 2  : 27 JUN 2005 (R2) or later
                 * version = 0  : version unrecognized;
                 *                file will not be read */

  int have_pcgamess; /* this flag is set to 1 if the output
                      * file is recognized as a PC Gamess output
                      * file; we might need to introduce a few
                      * switches in the code depending on if
                      * the log file is plain Gamess or PC Gamess
                      */

  int have_fmo;      /* set to 1 if this is an FMO calculation */

} gmsdata;


/* ######################################################## */
/* declaration/documentation of internal (static) functions */
/* ######################################################## */

/* this routine is the main gamess log file
 * parser responsible for static, i.e. 
 * non-trajectory information */
static int parse_static_data(qmdata_t *, int *);

/* for debugging */
static void print_input_data(qmdata_t *);

/* this routine checks if the current run is an
 * actual GAMESS run; returns true/false */
static int have_gamess(qmdata_t *, gmsdata *);

/* this function reads the number of processors requested */
static int get_proc_mem(qmdata_t *, gmsdata *);

/* Parse the $BASIS options*/
static int get_basis_options(qmdata_t *);

/* Determine the run title line */
static int get_runtitle(qmdata_t *);

/* Read the input atom definitions and geometry */
static int get_input_structure(qmdata_t *data, gmsdata *gms);

/* Read basis set and orbital statistics such as
 * # of shells, # of A/B orbitals, # of electrons, 
 * multiplicity and total charge */
static int get_basis_stats(qmdata_t *);

/* Read the contrl group for firefly calc and check for
 * supported RUNTYPes. Terminate the plugin
 * if an unsupported one is encountered. */
static int get_contrl_firefly(qmdata_t *);

/* Read the contrl group and check for
 * supported RUNTYPes. Terminate the plugin
 * if an unsupported one is encountered. */
static int get_contrl(qmdata_t *);

/* Read input parameters regarding calculation of 
 * certain molecular properties such as electrostatic
 * moments and the MEP. */
static int get_properties_input(qmdata_t *);

/* Read symmetry point group and highest axis */
static int get_symmetry(qmdata_t *);

/* Read in the $GUESS options */
static int get_guess_options(qmdata_t *);

/* Read MCSCF data */
static int get_mcscf(qmdata_t *data);

/* the function get_initial_info provides the atom number,
 * coordinates, and atom types and stores them
 * temporarily. */ 
static int get_final_properties (qmdata_t *);

/* Read atom coordinate block */
static int get_coordinates(FILE *file, qm_atom_t **atoms, int unit,
                           int *numatoms);

/* Read coordinates from $FMOXYZ section in the INPUT CARD
 * listing at the beginning of the file. */
static int get_fmoxyz(FILE *file, qm_atom_t **atoms, int unit,
                      int *numatoms);

/* Read the basis set data */
static int get_basis (qmdata_t *);

/* Read all primitives for the current shell */
static int read_shell_primitives(qmdata_t *, prim_t **prim,
                                 char *shelltype, int icoeff, int pcgamess);

/* convert shell type from char to int */
static int shelltype_int(char type);

/* Populate the flat arrays containing the basis set data */
static int fill_basis_arrays(qmdata_t *);

/* Read the first trajectory frame. */
static int read_first_frame(qmdata_t *);

/* Read next trajectory frame. */
static int get_traj_frame(qmdata_t *, qm_atom_t *, int);

/* returns 1 if the optimization has converged */
static int analyze_traj(qmdata_t *, gmsdata *);

/* read the number of scf iterations and the scf energies
 * for the current timestep. */
static int get_scfdata(qmdata_t *, qm_timestep_t *);

/* Reads a set of wavefunctions for the current timestep.
 * (typically alpha and beta spin wavefunctions) */
static int check_add_wavefunctions(qmdata_t *data,
                                   qm_timestep_t *ts);

/* Parse the wavefunction. */
static int get_wavefunction(qmdata_t *, qm_timestep_t *, qm_wavefunction_t *);

/* Read the wavefunction coefficients. */
static int read_coeff_block(FILE *file, int wavef_size,
                              float *wave_coeff, int *angular_momentum);

/* Read localized orbitals (Boys/Ruedenberg/Pipek) */
static int read_localized_orbitals(qmdata_t *data);

/* Read population analysis (Mulliken and Lowdin charges) */
static int get_population(qmdata_t *, qm_timestep_t *);

/* Read the energy gradient for each atom. */
static int get_gradient(qmdata_t *, qm_timestep_t *ts);

/* Read energy gradient from final traj step. */
static int get_final_gradient(qmdata_t *, qm_timestep_t *ts);

/* Read ESP charges. */
static int get_esp_charges(qmdata_t *data);

/* For runtyp=HESSIAN, this subroutine scans the file for 
 * the hessian matrix in internal coordinates 
 * as well as the internal coordinate information */
static int get_int_coords(qmdata_t *);

/* get Hessian matrix in internal coordinates */
static int get_int_hessian(qmdata_t *);

/* For runtyp=HESSIAN, this subroutine scans the file for 
 * the cartesian hessian matrix */ 
static int get_cart_hessian(qmdata_t *);

/* For runtyp=HESSIAN, this subroutine reads the frequencies
 * and intensities of the normal modes */
static int get_normal_modes(qmdata_t *);



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
  qmdata_t *data;
  gmsdata *gms;

  /* open the input file */
  fd = fopen(filename, "rb");
 
  if (!fd) {
    PRINTERR;
    return NULL;
  }

  /* allocate and initialize main data structure */
  data = init_qmdata(data);

  /* make sure memory was allocated properly */
  if (data == NULL) {
    PRINTERR;
    return NULL;
  }

  /* allocate GAMESS specific data */
  gms = (gmsdata *)calloc(1,sizeof(gmsdata));
  data->format_specific_data = gms;

  gms->version = 0;
  gms->have_pcgamess = 0;
  gms->have_fmo = 0;


  /* store file pointer and filename in gamess struct */
  data->file = fd;


  /* check if the file is GAMESS format; if yes
   * parse it, if no exit */
  if (have_gamess(data, gms)==TRUE) {
    if (gms->have_pcgamess) {
      printf("gamessplugin) Warning: PC GAMESS/FIREFLY is not yet fully supported!\n");
    

    }
    /* if we're dealing with an unsupported GAMESS
     * version, we better quit */
    if (gms->version==0) {
      printf("gamessplugin) GAMESS version %s not supported. \n",
             data->version_string);
      return NULL;
    }

    /* get the non-trajectory information from the log file */    
    if (parse_static_data(data, natoms) == FALSE) 
      return NULL;
  }
  else {
    printf("gamessplugin) This seems to not be a GAMESS logfile.\n");
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
  qmdata_t *data = (qmdata_t *)mydata;
  qm_atom_t *cur_atom;
  molfile_atom_t *atom;
  int i = 0;
 
  /* optional atomic number provided */
  *optflags = MOLFILE_ATOMICNUMBER;
  /*  if (data->have_mulliken) 
   *optflags |= MOLFILE_CHARGE;*/

  /* all the information I need has already been read in
   * via the initial scan and I simply need to copy 
   * everything from the temporary arrays into the 
   * proper VMD arrays.
   * Since there are no atom names in the GAMESS output
   * I use the atom type here --- maybe there is a better
   * way to do this !!?? */

  /* get initial pointer for atom array */
  cur_atom = data->atoms;

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

    /* if (data->have_mulliken)
      atom->charge = data->qm_timestep->mulliken_charges[i];
    */
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

  qmdata_t *data = (qmdata_t *)mydata;

  if (data->runtype == MOLFILE_RUNTYPE_HESSIAN) {
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

  /* hessian info */
  metadata->have_carthessian = data->have_cart_hessian;
  metadata->have_inthessian  = data->have_int_hessian;

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

  qmdata_t *data = (qmdata_t *)mydata;
  int i, j;
  int ncart;
  molfile_qm_hessian_t *hessian_data = &qm_data->hess;
  molfile_qm_basis_t   *basis_data   = &qm_data->basis;
  molfile_qm_sysinfo_t *sys_data     = &qm_data->run;

  /* fill in molfile_qm_hessian_t */
  if (data->runtype == MOLFILE_RUNTYPE_HESSIAN) {
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
    if (data->have_int_hessian) {
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
      /*printf("imag_modes[%d]=%d\n", i, data->imag_modes[i]);*/
      hessian_data->imag_modes[i] = data->imag_modes[i];
    }
  }

  /* fill in molfile_qm_sysinfo_t */
  sys_data->runtype = data->runtype;
  sys_data->scftype = data->scftype;
  sys_data->nproc   = data->nproc;
  sys_data->num_electrons  = data->num_electrons;
  sys_data->totalcharge    = data->totalcharge;
  sys_data->num_occupied_A = data->num_occupied_A;
  sys_data->num_occupied_B = data->num_occupied_B;
  sys_data->status         = data->status;


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
      basis_data->shell_types[i] = data->shell_types[i];
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

  qmdata_t *data = (qmdata_t *)mydata;

  meta->count = -1; /* Don't know the number of frames yet */

  if (data->num_frames_read > data->num_frames_sent) {
    have = 1;
  }
  else if (data->num_frames_read < data->num_frames) {
    /*printf("gamessplugin) Probing timestep %d\n", data->num_frames_read);*/

    have = get_traj_frame(data, data->atoms, data->numatoms);
  }

  if (have) {
    int i;
    qm_timestep_t *cur_ts;

    /* get a pointer to the current qm timestep */
    cur_ts = data->qm_timestep+data->num_frames_sent;
    
    for (i=0; (i<MOLFILE_MAXWAVEPERTS && i<cur_ts->numwave); i++) {
      meta->num_orbitals_per_wavef[i] = cur_ts->wave[i].num_orbitals;
      meta->has_occup_per_wavef[i]    = cur_ts->wave[i].has_occup;
      meta->has_orben_per_wavef[i]    = cur_ts->wave[i].has_orben;
    }
    meta->wavef_size      = data->wavef_size;
    meta->num_wavef       = cur_ts->numwave;
    meta->num_scfiter     = cur_ts->num_scfiter;
    meta->num_charge_sets = cur_ts->have_mulliken +
      cur_ts->have_lowdin + cur_ts->have_esp;
    if (cur_ts->gradient) meta->has_gradient = TRUE;

  } else {
    meta->has_gradient = FALSE;
    meta->num_scfiter  = 0;
    meta->num_orbitals_per_wavef[0] = 0;
    meta->has_occup_per_wavef[0] = FALSE;
    meta->num_wavef = 0;
    meta->wavef_size = 0;
    meta->num_charge_sets = 0;
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
  qmdata_t *data = (qmdata_t *)mydata;
  qm_timestep_t *cur_ts;
  int offset;
  int i = 0;
  int num_charge_sets = 0;

  if (data->trajectory_done == TRUE) return MOLFILE_ERROR;


  /* copy the coordinates */
  for (i=0; i<natoms; i++) {
    ts->coords[3*i  ] = data->atoms[i].x;
    ts->coords[3*i+1] = data->atoms[i].y;
    ts->coords[3*i+2] = data->atoms[i].z; 
  }    
  
  /* get a convenient pointer to the current qm timestep */
  cur_ts = data->qm_timestep+data->num_frames_sent;

  /* store the SCF energies */
  for (i=0; i<cur_ts->num_scfiter; i++) {
    qm_ts->scfenergies[i] = cur_ts->scfenergies[i];
  }

  /* store gradients */
  if (cur_ts->gradient) {
    for (i=0; i<3*natoms; i++) {
      qm_ts->gradient[i] = cur_ts->gradient[i];
    }
  }

  /* store charge sets*/
  if (cur_ts->have_mulliken) {
    offset = num_charge_sets*data->numatoms;
    for (i=0; i<data->numatoms; i++) {
      qm_ts->charges[offset+i] = cur_ts->mulliken_charges[i];
    }
    qm_ts->charge_types[num_charge_sets] = MOLFILE_QMCHARGE_MULLIKEN;
    num_charge_sets++;
  }

  if (cur_ts->have_lowdin) {
    offset = num_charge_sets*data->numatoms;
    for (i=0; i<data->numatoms; i++) {
      qm_ts->charges[offset+i] = cur_ts->lowdin_charges[i];
    }
    qm_ts->charge_types[num_charge_sets] = MOLFILE_QMCHARGE_LOWDIN;
    num_charge_sets++;
  }
  if (cur_ts->have_esp) {
    offset = num_charge_sets*data->numatoms;
    for (i=0; i<data->numatoms; i++) {
      qm_ts->charges[offset+i] = cur_ts->esp_charges[i];
    }
    qm_ts->charge_types[num_charge_sets] = MOLFILE_QMCHARGE_ESP;
    num_charge_sets++;
  }


  /* store the wave function and orbital energies */
  if (cur_ts->wave) {
    for (i=0; i<cur_ts->numwave; i++) {
      qm_wavefunction_t *wave = &cur_ts->wave[i];
      qm_ts->wave[i].type         = wave->type;
      qm_ts->wave[i].spin         = wave->spin;
      qm_ts->wave[i].excitation   = wave->exci;
      qm_ts->wave[i].multiplicity = wave->mult;
      qm_ts->wave[i].energy       = wave->energy;
      strncpy(qm_ts->wave[i].info, wave->info, MOLFILE_BUFSIZ);

      if (wave->wave_coeffs) {
        memcpy(qm_ts->wave[i].wave_coeffs, wave->wave_coeffs,
               wave->num_orbitals*data->wavef_size*sizeof(float));
      }
      if (wave->orb_energies) {
        memcpy(qm_ts->wave[i].orbital_energies, wave->orb_energies,
               wave->num_orbitals*sizeof(float));
      }
      if (wave->has_occup) {
        memcpy(qm_ts->wave[i].occupancies, wave->orb_occupancies,
               wave->num_orbitals*sizeof(float));
      }
    }
  }

  if (data->runtype == MOLFILE_RUNTYPE_ENERGY || 
      data->runtype == MOLFILE_RUNTYPE_HESSIAN) {
    /* We have only a single point */
    data->trajectory_done = TRUE;
  }

  data->num_frames_sent++;

  return MOLFILE_SUCCESS;
}
#endif



/**********************************************************
 *
 * close file and free memory
 *
 **********************************************************/
static void close_gamess_read(void *mydata) {

  qmdata_t *data = (qmdata_t *)mydata;
  int i, j;
  fclose(data->file);

  free(data->atoms);
  free(data->basis);
  free(data->shell_types);
  free(data->atomicnum_per_basisatom);
  free(data->num_shells_per_atom);
  free(data->num_prim_per_shell);
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
  free(data->imag_modes);
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
    free(data->qm_timestep[i].lowdin_charges);
    free(data->qm_timestep[i].esp_charges);
    for (j=0; j<data->qm_timestep[i].numwave; j++) {
      free(data->qm_timestep[i].wave[j].wave_coeffs);
      free(data->qm_timestep[i].wave[j].orb_energies);
      free(data->qm_timestep[i].wave[j].orb_occupancies);
    }
    free(data->qm_timestep[i].wave);
  }
  free(data->qm_timestep);
  free(data->format_specific_data);
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
static int parse_static_data(qmdata_t *data, int *natoms) 
{
  /* Cast GAMESS specific data */
  gmsdata *gms = (gmsdata *)data->format_specific_data;


  /* Read # of procs and amount of requested memory */
  get_proc_mem(data, gms);

  /* Read the $BASIS options */
  if (!get_basis_options(data))   return FALSE;

  /* Read the run title */
  if (!get_runtitle(data))        return FALSE;

  /* Read the $CONTRL group;
   * It actually appears after the basis stats in the
   * output file but we jump ahead and read it here
   * because we need the units (ANGS/BOHR) before
   * reading the input structure. */
  if (gms->have_pcgamess){
    if (!get_contrl_firefly(data))          return FALSE;
  }
  else {
    if (!get_contrl(data))          return FALSE;
  }

  /* Read the input atom definitions and geometry */
  if (!get_input_structure(data, gms)) return FALSE;

  /* Read the basis set */
  if (!get_basis(data))           return FALSE; 

  /* Read the number of orbitals, electrons, 
   * charge, multiplicity, ... */
  if (!get_basis_stats(data))     return FALSE;

  /* Read input parameters regarding calculation of 
   * certain molecular properties such as electrostatic
   * moments and the MEP. */
  if (!get_properties_input(data)) return FALSE;

  /* Read internal coordinates */
  get_int_coords(data);

  /* Read symmetry point group and highest axis */
  if (!get_symmetry(data))         return FALSE;

  /* Read the $GUESS options */
  get_guess_options(data);

  if (data->scftype==MOLFILE_SCFTYPE_MCSCF) {
    if (!get_mcscf(data))          return FALSE;
  }

  /* Find the end of the trajectory and count the
   * frames on the way.
   * If no regular end is found we won't find any
   * properties to read and return. */
  if (!analyze_traj(data, gms)) {
    printf("gamessplugin) WARNING: Truncated or abnormally terminated file!\n\n");
  }


  /* provide VMD with the proper number of atoms */
  *natoms = data->numatoms;

  /* Read the first frame*/
  read_first_frame(data);

  /* Read the properties at the end of a calculation */
  get_final_properties(data);

#ifdef DEBUGGING 
  printf("gamessplugin) num_frames_read = %d\n", data->num_frames_read);
  printf("gamessplugin) num_frames_sent = %d\n", data->num_frames_sent);

  /* Test print the parsed data in same format as logfile */
  print_input_data(data);
#endif

  return TRUE;
}


#define TORF(x) (x ? "T" : "F")

static void print_input_data(qmdata_t *data) {
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
    printf(" %-8s %6d", data->atoms[i].type, data->atoms[i].atomicnum);
    
    printf("%17.10f",   ANGS_TO_BOHR*data->atoms[i].x);
    printf("%20.10f",   ANGS_TO_BOHR*data->atoms[i].y);
    printf("%20.10f\n", ANGS_TO_BOHR*data->atoms[i].z);
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
    printf("%-8s\n\n", data->atoms[i].type);
    printf("\n");
    printf("nshells=%d\n", data->num_shells_per_atom[i]);

    for (j=0; j<data->num_shells_per_atom[i]; j++) {
      printf("nprim=%d\n", data->num_prim_per_shell[shellcount]);

      for (k=0; k<data->num_prim_per_shell[shellcount]; k++) {
        printf("%6d   %d %7d %22f%22f\n", j, data->shell_types[shellcount],
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
    printf("%-8s (%10s)\n\n", data->atoms[i].type, data->basis_set[i].name);
    printf("\n");

    for (j=0; j<data->basis_set[i].numshells; j++) {

      for (k=0; k<data->basis_set[i].shell[j].numprims; k++) {
        printf("%6d   %d %7d %22f%22f\n", j,
               data->basis_set[i].shell[j].type,
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
static int have_gamess(qmdata_t *data, gmsdata *gms) 
{
  char word[3][BUFSIZ];
  char buffer[BUFSIZ];
  char versionstr[BUFSIZ];
  int day, year;
  char month[BUFSIZ], rev[BUFSIZ];
  int i = 0;
  int program;

  buffer[0] = '\0';
  for (i=0; i<3; i++) word[i][0] = '\0';


  /* check if the file is GAMESS format */
  program = goto_keyline(data->file,
                          "PC GAMESS version",
                          "GAMESS VERSION =", 
                          "Firefly version",NULL);
  if (program==1) {
    gms->have_pcgamess = 1;
    gms->version = 1;
    strcpy(data->version_string, "PC GAMESS ");
  } else if (program==2) {
    gms->have_pcgamess = 0;
    strcpy(data->version_string, "GAMESS ");
  } else if (program==3) {
    gms->have_pcgamess = 1;
    gms->version = 1;
    strcpy(data->version_string, "Firefly ");
  } else {
    printf("gamessplugin) This is no GAMESS/PCGAMESS logfile!\n");
    return FALSE;
  }

  GET_LINE(buffer, data->file);

  if (gms->have_pcgamess) {
    if (strstr(buffer,"version") != NULL) {
      strncpy(versionstr, strstr(buffer,"version")+8, 16);
      *strchr(versionstr, ' ') = '\0';
    }
  } else {
    /* extract the version number if possible; otherwise
     * return empty string */
    if (strstr(buffer,"=") != NULL) {
      strncpy(versionstr, strstr(buffer,"=")+2, 16);
      versionstr[16] = '\0';
    }
    
    /* determine if we're dealing with pre-"27 JUN 2005"
     * version */
    sscanf(versionstr, "%d %s %d %s", &day, month, &year, rev);
    
    if ( ( year >= 2006 ) ||
         ( year == 2005 && !strcmp(month,"JUN") ) ||
         ( year == 2005 && !strcmp(month,"NOV") ) ||
         ( year == 2005 && !strcmp(month,"DEC") ) )
      {
        gms->version = 2;
      } else { 
        gms->version = 1;
      }
  }

  strcat(data->version_string, versionstr);

  printf("gamessplugin) Version = %s\n", 
         data->version_string);

  return TRUE;
}


/**********************************************************
 *
 * this subroutine reads the number of procs and the amount
 * of memory requested
 *
 **********************************************************/
static int get_proc_mem(qmdata_t *data, gmsdata *gms) {

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
  if (gms->have_pcgamess == 1) {
    /* XXX for now we fake ncpu = 1 until we know exactly
     *     how the output format looks like */
    nproc = 1;
    do {
      GET_LINE(buffer, data->file);
      sscanf(buffer,"%s %d %s",&word[0][0],&nproc,&word[1][0]);
        if (!strcmp(&word[0][0],"PARALLEL") &&
               !strcmp(&word[0][1],"RUNNING")) {
            sscanf(buffer,"%*s %*s %*s %*s %*s %d %*s %*s",&nproc);
            break;
        }
      } while (strcmp(&word[0][0],"ECHO") || 
               strcmp(&word[1][0],"THE") );

  }
  else {
    do {
      GET_LINE(buffer, data->file);
      sscanf(buffer,"%s %d %s",&word[0][0],&nproc,&word[1][0]);

      if (!strcmp(&word[0][0],"Initiating") &&
          !strcmp(&word[1][0],"compute")) {
        break;
      }

      else if (!strcmp(&word[0][0],"Initiating") &&
               !strcmp(&word[1][0],"processes")) {
        break;
      }

      /* Some versions */
      else if (!strcmp(&word[0][0],"PARALLEL") &&
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
  if (gms->have_pcgamess == 1) {
    GET_LINE(buffer, data->file);

    /* store it */
    if ((temp = strstr(buffer,"MEMORY=")+8)==NULL) return FALSE;
    strncpy(data->memory,trimright(temp),sizeof(data->memory));
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
static int get_basis_options(qmdata_t *data) {

  char buffer[BUFSIZ];
  char diffuse[BUFSIZ];
  char polarization[BUFSIZ];
  int ngauss;

  buffer[0] = '\0';
  diffuse[0] = '\0';
  polarization[0] = '\0';

  /* The $BASIS section is somewhere at the beginning of
   * the file. Rewind to be sure not to miss it. */
  rewind(data->file);

  /* start scanning */
  if (pass_keyline(data->file, "BASIS OPTIONS",
                    "RUN TITLE") != FOUND) {
    /* No Basis options section found
     * (basis was entered explicitly) */
    return TRUE;
  }

  eatline(data->file, 1);


  /* the first string in the current line contains the
   * GBASIS used; copy it over into the gbasis variable
   * of qmdata_t */
  GET_LINE(buffer, data->file);
  sscanf(buffer," GBASIS=%s IGAUSS= %d", data->gbasis, &ngauss);
 

  /* in case we're using a pople style basis set, i.e. 
   * GBASIS=N311,N31,N21 or STO we also scan for the number 
   * of gaussians, as well as p,d,f and diffuse functions
   * and use this info to assemble a "basis set string" */
  if ( !strncmp(data->gbasis,"N311",sizeof(data->gbasis)) ||
       !strncmp(data->gbasis,"N31",sizeof(data->gbasis)) ||
       !strncmp(data->gbasis,"N21",sizeof(data->gbasis)) ||
       !strncmp(data->gbasis,"STO",sizeof(data->gbasis)) ) 
  {
    int npfunc, ndfunc, nffunc;
    int diffs=FALSE, diffsp=FALSE;
    char torf;

    /* the next line gives us the d,f and diffuse sp
     * functions */
    GET_LINE(buffer, data->file);
    if (sscanf(buffer," NDFUNC= %d NFFUNC= %d DIFFSP= %c",
               &ndfunc, &nffunc, &torf) != 3) {
      sscanf(buffer," NDFUNC= %d DIFFSP= %c", &ndfunc, &torf);
    }

    /* convert GAMESS' .TRUE./.FALSE. for DIFFSP into 1/0 */
    if (torf=='T') diffsp = TRUE;


    /* the next line gives us the p and diffuse s
     * functions */
    GET_LINE(buffer, data->file);
    sscanf(buffer," NPFUNC= %d DIFFS= %c", &npfunc, &torf);

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
static int get_runtitle(qmdata_t *data) {

  char buffer[BUFSIZ];

  if (pass_keyline(data->file, "RUN TITLE",
                   "THE POINT GROUP") != FOUND) {
    /* This is most likely a broken file, but who knows. 
     * Since we don't really care about the title string
     * we go on here without error. */
    data->runtitle[0] = '\0';
    return TRUE;
  }

  GET_LINE(buffer, data->file);
  strncpy(data->runtitle,trimright(buffer),sizeof(data->runtitle));

  return TRUE;
} 


/**********************************************************
 *
 * Read the input atom definitions and geometry
 *
 **********************************************************/
static int get_input_structure(qmdata_t *data, gmsdata *gms) {
  char buffer[BUFSIZ];
  char units[BUFSIZ];
  int numatoms = -1;
  int bohr;
  long filepos;
  filepos = ftell(data->file);

  /* See if we find the "ATOM      ATOMIC ..." line before 
   * any of the three stopstrings mrking the beginning of
   * possible following sections. */
  if (goto_keyline(data->file,
         "ATOM      ATOMIC                      COORDINATES (",
         "INTERNUCLEAR DISTANCES",
         "ATOMIC BASIS SET",
         "$CONTRL OPTIONS", NULL) == FOUND) {

    GET_LINE(buffer, data->file);
    sscanf(buffer, " ATOM      ATOMIC  %*s  %s", units);
    eatline(data->file, 1);

  } else {
    /* This is probably an FMO calc.; if so, set flag. */
    fseek(data->file, filepos, SEEK_SET);
    if (pass_keyline(data->file,
                     "The Fragment Molecular Orbital (FMO) method.", 
                     NULL)) {
      gms->have_fmo = 1;
      printf("gamessplugin) Fragment Molecular Orbital (FMO) method.\n");
    }

    /* We didn't find the normal input section.
     * Let's see i we can find coordinates for the first
     * frame of a trajectory. */
    fseek(data->file, filepos, SEEK_SET);
    if (pass_keyline(data->file,
                      "BEGINNING GEOMETRY SEARCH POINT NSERCH=   0",
                      NULL) &&
        goto_keyline(data->file, "COORDINATES OF ALL ATOMS", NULL)) {
      GET_LINE(buffer, data->file);
      sscanf(buffer, " COORDINATES OF ALL ATOMS ARE %s", units);
      eatline(data->file, 2);

    } else {
      /* As last resort look for FMO coordinates in the
       * INPUT CARD section: */

      /* But first we have to get the units from the $CONTRL section */
      rewind(data->file);
      if (!pass_keyline(data->file, "$CONTRL OPTIONS", NULL)) {
        printf("gamessplugin) Missing $CONTRL OPTIONS section!\n");
        return FALSE;
      }
      goto_keyline(data->file, "UNITS =", NULL);
      GET_LINE(buffer, data->file);
      sscanf(strstr(buffer, "UNITS ="), "%s", units);
      bohr = !strcmp(units, "BOHR");
      
      /* Find beginning of $FMOXYZ input card */
      rewind(data->file);
      if (!pass_keyline(data->file, "INPUT CARD> $fmoxyz", 
                        "INPUT CARD> $FMOXYZ")) {
        printf("gamessplugin) No atom coordinates found!\n");
        return FALSE;
      }
            
      /* Read the $FMOXYZ coordinates */     
      if (!get_fmoxyz(data->file, &data->atoms, bohr, &numatoms)) {
        printf("gamessplugin) Could not read coordinates from $FMOXYZ input!\n");
        return FALSE;
      } else {
        printf("gamessplugin) Fragment Molecular Orbital (FMO) method.\n");
        gms->have_fmo = 1;
        data->numatoms = numatoms;
        return TRUE;
      }
    }
  }

  /* If we reached this point we have found either a regular
   * input coordinate block or the coordinate block of the
   * first trajectory frame. */

  /* test if coordinate units are Bohr */
  bohr = !strcmp(units, "(BOHR)");

  /* Read the coordinate block */
  if (get_coordinates(data->file, &data->atoms, bohr, &numatoms))
    data->num_frames_read = 0;
  else {
    printf("gamessplugin) Bad atom coordinate block!\n");
    return FALSE;
  }

  fseek(data->file, filepos, SEEK_SET);

  /* store number of atoms in data structure */
  data->numatoms = numatoms;

  return TRUE; 
}


/**********************************************************
 *
 * Read an atom coordinate block.
 *
 * Example:
 *  F         9.0    3.04259     -0.07605       0.00000
 *  N         7.0    0.03017      0.38347       0.00000
 * ...
 *
 **********************************************************/
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
    float x,y,z, dum;
    int n;
    qm_atom_t *atm;

    GET_LINE(buffer, file);

    /* For FMO there is an additional atom index in the
     * second column. Try both variants: */
    n = sscanf(buffer,"%s %f %f %f %f %f",atname,&dum,&atomicnum,&x,&y,&z);
    if (n!=6) {
      n = sscanf(buffer,"%s %f %f %f %f",atname,&atomicnum,&x,&y,&z);
    }
    if (n!=5 && n!=6) break;

    if (growarray && i>0) {
      *atoms = (qm_atom_t*)realloc(*atoms, (i+1)*sizeof(qm_atom_t));
    }
    atm = (*atoms)+i;

    strncpy(atm->type, atname, sizeof(atm->type));
    atm->atomicnum = floor(atomicnum+0.5); /* nuclear charge */
    /*printf("coor: %s %d %f %f %f\n", atm->type, atm->atomicnum, x, y, z);*/
   
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
  if (*numatoms>=0 && *numatoms!=i) {
    (*numatoms) = i;
    return FALSE;
  }

  (*numatoms) = i;
  return TRUE;
}


/* Read coordinates from $FMOXYZ section in the INPUT CARD
 * listing at the beginning of the file. This is a method
 * of last resort used only for FMO calculations where no
 * coordinates are printed. */
static int get_fmoxyz(FILE *file, qm_atom_t **atoms, int unit,
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
    char atname[BUFSIZ], element[BUFSIZ];
    float x,y,z;
    int n;
    qm_atom_t *atm;

    GET_LINE(buffer, file);

    /* skip " INPUT CARD>" at the beginning of the line */
    n = sscanf(buffer+12,"%s %s %f %f %f",atname,element,&x,&y,&z);

    if (n!=5) break;

    if (growarray && i>0) {
      *atoms = (qm_atom_t*)realloc(*atoms, (i+1)*sizeof(qm_atom_t));
    }
    atm = (*atoms)+i;

    strncpy(atm->type, atname, sizeof(atm->type));
    if (isalpha(element[0]))
      atm->atomicnum = get_pte_idx_from_string(element);
    else if (isdigit(element[0])) 
      atm->atomicnum = floor(element[0]+0.5); /* nuclear charge */
    else break;

    /* If coordinates are in Bohr convert them to Angstrom */
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


/**********************************************************
 *
 * Read data from the $CONTRL card from FIREFLY
 *
 **********************************************************/
static int get_contrl_firefly(qmdata_t *data) {

  char word[3][BUFSIZ];
  char buffer[BUFSIZ];
  char *temp;
  long filepos;
  filepos = ftell(data->file);

  word[0][0] = '\0';
  word[1][0] = '\0';
  word[2][0] = '\0';
  buffer[0] = '\0';

  


  /* start scanning; currently we support
   * RUNTYP = ENERGY, OPTIMIZE, SADPOINT, HESSIAN, SURFACE */
  if (!pass_keyline(data->file, "$CONTRL OPTIONS", NULL)) {
    fseek(data->file, filepos, SEEK_SET);
    return FALSE;
  }

  eatline(data->file, 1);

  /* current line contains SCFTYP, RUNTYP, EXETYP info; scan it */
  GET_LINE(buffer, data->file);
  sscanf(buffer,"%s %s",&word[0][0],&word[1][0]);

  /* check for supported RUNTYPs */
  if      (!strcmp(&word[1][0],"RUNTYP=ENERGY")) {
    data->runtype = MOLFILE_RUNTYPE_ENERGY;
  }
  else if (!strcmp(&word[1][0],"RUNTYP=OPTIMIZE")) {
    data->runtype = MOLFILE_RUNTYPE_OPTIMIZE;
  }
  else if (!strcmp(&word[1][0],"RUNTYP=SADPOINT")) {
    data->runtype = MOLFILE_RUNTYPE_SADPOINT;
  }
  else if (!strcmp(&word[1][0],"RUNTYP=HESSIAN")) {
    data->runtype = MOLFILE_RUNTYPE_HESSIAN;
  }
  else if (!strcmp(&word[1][0],"RUNTYP=SURFACE")) {
    data->runtype = MOLFILE_RUNTYPE_SURFACE;
  }
  else if (!strcmp(&word[1][0],"RUNTYP=GRADIENT")) {
    data->runtype = MOLFILE_RUNTYPE_GRADIENT;
  }
  else if (!strcmp(&word[1][0],"RUNTYP=MEX")) {
    data->runtype = MOLFILE_RUNTYPE_MEX;
  }
  else {
#ifdef DEBUGGING
    printf("gamessplugin) The %s is currently not supported \n",
           &word[1][0]);
#endif
    data->runtype = MOLFILE_RUNTYPE_UNKNOWN;
  }
  printf("gamessplugin) File generated via %s \n",&word[1][0]);


  /* determine SCFTYP */
  if (!strcmp(&word[0][0],"SCFTYP=RHF")) {
    data->scftype = MOLFILE_SCFTYPE_RHF;
  }
  else if (!strcmp(&word[0][0],"SCFTYP=UHF")) {
    data->scftype = MOLFILE_SCFTYPE_UHF;
  }
  else if (!strcmp(&word[0][0],"SCFTYP=ROHF")) {
    data->scftype = MOLFILE_SCFTYPE_ROHF;
  }
  else if (!strcmp(&word[0][0],"SCFTYP=GVB")) {
    data->scftype = MOLFILE_SCFTYPE_GVB;
  }
  else if (!strcmp(&word[0][0],"SCFTYP=MCSCF")) {
    data->scftype = MOLFILE_SCFTYPE_MCSCF;
  }
  else if (!strcmp(&word[0][0],"SCFTYP=NONE")) {
    data->scftype = MOLFILE_SCFTYPE_NONE;
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

  /* scan for MPLEVL, LOCAL and UNITS; */
  GET_LINE(buffer, data->file);
  sscanf(buffer,"%s %s %*s %s",&word[0][0],&word[1][0],&word[2][0]);
  data->mplevel = atoi(&word[1][0]);

  /* scan for MULT, ICHARG and MAXIT; */
  GET_LINE(buffer, data->file);
  


  /* find the coordinate type in next line */
  while ( (temp=strstr(buffer,"COORD =")) == NULL ) {
    GET_LINE(buffer, data->file);;
  }
  strncpy(data->geometry, trimright(temp+7), sizeof(data->geometry)); 
  printf("gamessplugin) Coordinate type used is %s \n", data->geometry);

  while ( (temp=strstr(buffer,"CITYP =")) == NULL ) {
    GET_LINE(buffer, data->file);;
  }
  strncpy(buffer, trimright(temp+7), 8); 

  /* determine CITYP */
  if      (!strcmp(buffer,"NONE"))  data->citype = CI_NONE;
  else if (!strcmp(buffer,"CIS"))   data->citype = CI_CIS;
  else if (!strcmp(buffer,"ALDET")) data->citype = CI_ALDET;
  else if (!strcmp(buffer,"ORMAS")) data->citype = CI_ORMAS;
  else if (!strcmp(buffer,"GUGA"))  data->citype = CI_GUGA;
  else if (!strcmp(buffer,"FSOCI")) data->citype = CI_FSOCI;
  else if (!strcmp(buffer,"GENCI")) data->citype = CI_GENCI;
  else                                    data->citype = CI_UNKNOWN;
  printf("gamessplugin) CI method %s \n",buffer);

  GET_LINE(buffer, data->file);
  sscanf(buffer,"%s %*s",&word[0][0]);

  /* scan for DFTTYP, TDDFT info; */
  if (!strncmp(&word[0][0],"DFTTYP=", 7)) {
    printf("gamessplugin) Density functional used is %s \n",&word[0][7]);
    GET_LINE(buffer, data->file);
  }



  fseek(data->file, filepos, SEEK_SET);
  return TRUE;
}


/**********************************************************
 *
 * Read data from the $CONTRL card from GAMESS
 *
 **********************************************************/
static int get_contrl(qmdata_t *data) {

  char word[3][BUFSIZ];
  char buffer[BUFSIZ];
  char *temp;
  long filepos;
  filepos = ftell(data->file);

  word[0][0] = '\0';
  word[1][0] = '\0';
  word[2][0] = '\0';
  buffer[0] = '\0';


  /* start scanning; currently we support
   * RUNTYP = ENERGY, OPTIMIZE, SADPOINT, HESSIAN, SURFACE */
  if (!pass_keyline(data->file, "$CONTRL OPTIONS", NULL)) {
    fseek(data->file, filepos, SEEK_SET);
    return FALSE;
  }

  eatline(data->file, 1);

  /* current line contains SCFTYP, RUNTYP, EXETYP info; scan it */
  GET_LINE(buffer, data->file);
  sscanf(buffer,"%s %s",&word[0][0],&word[1][0]);

  /* check for supported RUNTYPs */
  if      (!strcmp(&word[1][0],"RUNTYP=ENERGY")) {
    data->runtype = MOLFILE_RUNTYPE_ENERGY;
  }
  else if (!strcmp(&word[1][0],"RUNTYP=OPTIMIZE")) {
    data->runtype = MOLFILE_RUNTYPE_OPTIMIZE;
  }
  else if (!strcmp(&word[1][0],"RUNTYP=SADPOINT")) {
    data->runtype = MOLFILE_RUNTYPE_SADPOINT;
  }
  else if (!strcmp(&word[1][0],"RUNTYP=HESSIAN")) {
    data->runtype = MOLFILE_RUNTYPE_HESSIAN;
  }
  else if (!strcmp(&word[1][0],"RUNTYP=SURFACE")) {
    data->runtype = MOLFILE_RUNTYPE_SURFACE;
  }
  else if (!strcmp(&word[1][0],"RUNTYP=GRADIENT")) {
    data->runtype = MOLFILE_RUNTYPE_GRADIENT;
  }
  else if (!strcmp(&word[1][0],"RUNTYP=MEX")) {
    data->runtype = MOLFILE_RUNTYPE_MEX;
  }
  else {
#ifdef DEBUGGING
    printf("gamessplugin) The %s is currently not supported \n",
           &word[1][0]);
#endif
    data->runtype = MOLFILE_RUNTYPE_UNKNOWN;
  }
  printf("gamessplugin) File generated via %s \n",&word[1][0]);


  /* determine SCFTYP */
  if (!strcmp(&word[0][0],"SCFTYP=RHF")) {
    data->scftype = MOLFILE_SCFTYPE_RHF;
  }
  else if (!strcmp(&word[0][0],"SCFTYP=UHF")) {
    data->scftype = MOLFILE_SCFTYPE_UHF;
  }
  else if (!strcmp(&word[0][0],"SCFTYP=ROHF")) {
    data->scftype = MOLFILE_SCFTYPE_ROHF;
  }
  else if (!strcmp(&word[0][0],"SCFTYP=GVB")) {
    data->scftype = MOLFILE_SCFTYPE_GVB;
  }
  else if (!strcmp(&word[0][0],"SCFTYP=MCSCF")) {
    data->scftype = MOLFILE_SCFTYPE_MCSCF;
  }
  else if (!strcmp(&word[0][0],"SCFTYP=NONE")) {
    data->scftype = MOLFILE_SCFTYPE_NONE;
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
    if      (!strcmp(&word[2][0],"=NONE"))  data->citype = CI_NONE;
    else if (!strcmp(&word[2][0],"=CIS"))   data->citype = CI_CIS;
    else if (!strcmp(&word[2][0],"=ALDET")) data->citype = CI_ALDET;
    else if (!strcmp(&word[2][0],"=ORMAS")) data->citype = CI_ORMAS;
    else if (!strcmp(&word[2][0],"=GUGA"))  data->citype = CI_GUGA;
    else if (!strcmp(&word[2][0],"=FSOCI")) data->citype = CI_FSOCI;
    else if (!strcmp(&word[2][0],"=GENCI")) data->citype = CI_GENCI;
    else                                    data->citype = CI_UNKNOWN;
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
  strncpy(data->geometry, trimright(temp+7), sizeof(data->geometry)); 
  printf("gamessplugin) Coordinate type used is %s \n", data->geometry);

  fseek(data->file, filepos, SEEK_SET);
  return TRUE;
}


/* Read input parameters regarding calculation of 
 * certain molecular properties such as electrostatic
 * moments and the MEP. */
static int get_properties_input(qmdata_t *data) {
  /* TODO!! */
  return TRUE;
}


/* Read symmetry point group and highest axis.
 * Currently these values are not used yet, but if this
 * section was not found the file is corrupt. */
static int get_symmetry(qmdata_t *data) {
  char buffer[BUFSIZ];
  char *sep, tmp[BUFSIZ];

  

  

  long filepos = ftell(data->file);

  if (goto_keyline(data->file, "THE POINT GROUP IS",
		    "1 ELECTRON INTEGRALS", NULL) != FOUND) {
    printf("gamessplugin) No symmetry info found!\n");
    return FALSE;
  }

  GET_LINE(buffer, data->file);
  sscanf(buffer," THE POINT GROUP IS %s", data->pointgroup);

  sep = strchr(data->pointgroup, ',');
  if (sep) *sep = '\0';
  trimright(data->pointgroup);

  sep = strstr(buffer, "NAXIS=") + 6;
  strncpy(tmp, sep, 2); tmp[2] = '\0';
  data->naxis = atoi(tmp);

  sep = strstr(buffer, "ORDER=") + 6;
  sscanf(sep, "%d", &data->order);

  printf("gamessplugin) Point group = %s, naxis = %d, order = %d\n",
         data->pointgroup, data->naxis, data->order);

  fseek(data->file, filepos, SEEK_SET);

  return TRUE;
}


/* Read MCSCF input data */
static int get_mcscf(qmdata_t *data) {
  gmsdata *gms = (gmsdata *)data->format_specific_data;
  char buffer[BUFSIZ];
  char *temp;
  long filepos;
  int tmp;

  filepos = ftell(data->file);

  if (gms->have_pcgamess){
      if (pass_keyline(data->file,"XMCQDPT INPUT PARAMETERS",
                        "DONE SETTING UP THE RUN") != FOUND) {
         

         if(pass_keyline(data->file, "MCSCF CALCULATION",
                       "ITER     TOTAL ENERGY") != FOUND)
                return FALSE;

          if (goto_keyline(data->file, "-CORE-    -INTERNAL-  -EXTERNAL-",
                           "ITER     TOTAL ENERGY", NULL) != FOUND)
            return FALSE;

          while ( (temp=strstr(buffer,"NFZC=")) == NULL ) {
            GET_LINE(buffer, data->file);
          }
          strncpy(buffer, trimright(temp+6), 5); 
          sscanf(buffer, "%d", &data->mcscf_num_core);

          while ( (temp=strstr(buffer,"NMCC=")) == NULL ) {
            GET_LINE(buffer, data->file);
          }
          strncpy(buffer, trimright(temp+6), 5); 
          sscanf(buffer, "%d", &tmp);
          data-> mcscf_num_core += tmp;
          printf("gamessplugin) Number of MCSCF core orbitals = %d\n",
             data->mcscf_num_core);
      }
      else{
          

          while ( (temp=strstr(buffer,"# OF FROZEN CORE ORBITALS")) == NULL ) {
            GET_LINE(buffer, data->file);
          }
          sscanf(buffer, "%*s %*s %*s %*s %*s %*s %d",&data->mcscf_num_core);

          GET_LINE(buffer,data->file);
          sscanf(buffer, "%*s %*s %*s %*s %*s %*s %d",&tmp);
          data->mcscf_num_core += tmp;
          printf("gamessplugin) Number of MCSCF core orbitals = %d\n",
             data->mcscf_num_core);
          printf("gamessplugin) XMCQDPT2 not supported.\n");
          

          data->scftype = MOLFILE_SCFTYPE_NONE;

      } 
  } 
  else {
      if (pass_keyline(data->file, "MCSCF CALCULATION",
                       "ITER     TOTAL ENERGY") != FOUND)
        return FALSE;

      if (goto_keyline(data->file, "NUMBER OF CORE ORBITALS",
                       "ITER     TOTAL ENERGY", NULL) != FOUND)
        return FALSE;

      GET_LINE(buffer, data->file);
      sscanf(buffer," NUMBER OF CORE ORBITALS          = %d",
             &data->mcscf_num_core);

      printf("gamessplugin) Number of MCSCF core orbitals = %d\n",
         data->mcscf_num_core);
  }

  fseek(data->file, filepos, SEEK_SET);
  return TRUE;
}


/* Read the first trajectory frame. */
static int read_first_frame(qmdata_t *data) {
  /* The angular momentum is populated in get_wavefunction 
   * which is called by get_traj_frame(). We have obtained
   * the array size wavef_size already from the basis set
   * statistics */
  data->angular_momentum = (int*)calloc(3*data->wavef_size, sizeof(int));

  /* Try reading the first frame. 
   * If there is only one frame then also read the
   * final wavefunction. */
  if (!get_traj_frame(data, data->atoms, data->numatoms)) {
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
static int get_final_properties(qmdata_t *data) {
  qm_timestep_t *ts;
  long filepos;
  filepos = ftell(data->file);
  ts = data->qm_timestep + data->num_frames-1;

  /* Go to end of trajectory */
  fseek(data->file, data->end_of_traj, SEEK_SET);

  printf("gamessplugin) Reading final properties section (timestep %d):\n",
         data->num_frames-1);
  printf("gamessplugin) ===============================================\n");

  /* Read population analysis (Mulliken and Lowdin charges),
   * but only if wasn't read already (while parsing the first
   * timestep). */
  if (!ts->have_mulliken && get_population(data, ts)) {
    printf("gamessplugin) Mulliken charges found\n");
  }

  if (get_esp_charges(data)) {
    printf("gamessplugin) ESP charges found\n");
  }

  if (data->runtype == MOLFILE_RUNTYPE_GRADIENT ||
      data->runtype == MOLFILE_RUNTYPE_HESSIAN) {
    if (get_final_gradient(data, ts)) {
      printf("gamessplugin) Final gradient found\n");
    }
  }


  if (data->runtype == MOLFILE_RUNTYPE_HESSIAN) {
    /* try reading the hessian matrix in internal and
     * cartesian coordinates as well as the internal
     * coordinates together with their associated
     * force constants */
    
    if (!get_int_hessian(data)) {
      printf("gamessplugin) No internal Hessian matrix found.\n");
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
      printf("gamessplugin) No normal modes found.\n");
      printf("gamessplugin) \n");
    }
  }

  /* Read localized orbitals if there are any */
  read_localized_orbitals(data);


  fseek(data->file, filepos, SEEK_SET);
  return TRUE; 
}


/* Read localized orbitals (Boys/Ruedenberg/Pipek) */
static int read_localized_orbitals(qmdata_t *data) {
  int i;
  qm_timestep_t *ts;
  qm_wavefunction_t *wavef;

  /* Move past the listing of the canonical MOs */
  pass_keyline(data->file, "ENERGY COMPONENTS", NULL);

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
      wavef->mult = data->multiplicity;
      wavef->energy = ts->scfenergies[ts->num_scfiter-1];

      printf("gamessplugin) Localized orbitals (%s) found for timestep %d\n",
             typestr, data->num_frames-1);
    }
  }

  return TRUE;
}



/********************************************************
 *
 * Read basis set and orbital statistics such as
 * # of shells, # of A/B orbitals, # of electrons, 
 * multiplicity and total charge
 *
 ********************************************************/
static int get_basis_stats(qmdata_t *data) {

  gmsdata *gms = (gmsdata *)data->format_specific_data;

  char buffer[BUFSIZ];
  char word[7][BUFSIZ];
  int i;

  buffer[0] = '\0';
  for (i=0; i<7; i++) word[i][0] = '\0';

  /* look for the orbital/charge/... info section */
  if(gms->have_pcgamess){
    if (!pass_keyline(data->file, "TOTAL NUMBER OF SHELLS", NULL)){
        printf("ERROR!\n");
        return FALSE;
     }
      /* # cartesian gaussian function = wavefunction size */
      GET_LINE(buffer, data->file);
      sscanf(buffer,"%*s %*s %*s %*s %*s %*s %d",
             &(data->wavef_size));

      /* number of electrons */
      GET_LINE(buffer, data->file);
      sscanf(buffer,"%*s %*s %*s %*s %d",
             &(data->num_electrons));

      /* charge of the molecule */
      GET_LINE(buffer, data->file);
      sscanf(buffer,"%*s %*s %*s %*s %d",
             &(data->totalcharge));

      /* Multiplicity of the molecule.
       * Multiplicity is actually defined per wavefunction
       * but in some cases where there's only one wavefunction
       * in the output or they all have the same multiplicity
       * it will not be printed with the wavefunction.
       * Thus we use this one as default value. */
      GET_LINE(buffer, data->file);
      sscanf(buffer,"%*s %*s %*s %d",
             &(data->multiplicity));

      /* number of A orbitals */
      /* Note the different number of items per line for A/B orbitals
       * due to "(ALPHA)" and "(BETA )" !! */
      GET_LINE(buffer, data->file);
      sscanf(buffer,"%*s %*s %*s %*s %*s %*s %d",
             &(data->num_occupied_A)); 
        
      /* number of B orbitals */
      GET_LINE(buffer, data->file);
      sscanf(buffer,"%*s %*s %*s %*s %*s %*s %*s %d",
             &(data->num_occupied_B)); 

  }
  else {
    if (!pass_keyline(data->file, "TOTAL NUMBER OF BASIS", NULL))
        return FALSE;

      /* # cartesian gaussian function = wavefunction size */
      GET_LINE(buffer, data->file);
      sscanf(buffer,"%*s %*s %*s %*s %*s %*s %*s %d",
             &(data->wavef_size));

      /* number of electrons */
      GET_LINE(buffer, data->file);
      sscanf(buffer,"%*s %*s %*s %*s %d",
             &(data->num_electrons));

      /* charge of the molecule */
      GET_LINE(buffer, data->file);
      sscanf(buffer,"%*s %*s %*s %*s %d",
             &(data->totalcharge));

      /* Multiplicity of the molecule.
       * Multiplicity is actually defined per wavefunction
       * but in some cases where there's only one wavefunction
       * in the output or they all have the same multiplicity
       * it will not be printed with the wavefunction.
       * Thus we use this one as default value. */
      GET_LINE(buffer, data->file);
      sscanf(buffer,"%*s %*s %*s %d",
             &(data->multiplicity));

      /* number of A orbitals */
      /* Note the different number of items per line for A/B orbitals
       * due to "(ALPHA)" and "(BETA )" !! */
      GET_LINE(buffer, data->file);
      sscanf(buffer,"%*s %*s %*s %*s %*s %*s %d",
             &(data->num_occupied_A)); 
        
      /* number of B orbitals */
      GET_LINE(buffer, data->file);
      sscanf(buffer,"%*s %*s %*s %*s %*s %*s %*s %d",
             &(data->num_occupied_B)); 
  }


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
static int get_guess_options(qmdata_t *data)
{
  char word[BUFSIZ];
  char buffer[BUFSIZ];
  long filepos;
  filepos = ftell(data->file);

  /* initialize buffers */
  buffer[0] = '\0';
  word[0]   = '\0';

  /* parse for GUESS field */
  if (pass_keyline(data->file, "GUESS OPTIONS",
                   "2 ELECTRON INTEGRALS") != FOUND) {
    printf("gamessplugin) No GUESS OPTIONS found.\n");
    fseek(data->file, filepos, SEEK_SET);

    /* This section id not mandatory, there are a few
       calculation types the don't print it, so we
       always return TRUE*/
    return TRUE;
  }

  /* next line contains all we need */
  eatline(data->file, 1);
  GET_LINE(buffer, data->file);
  sscanf(buffer," GUESS %s NORB",&word[0]);

  /* the first character is '=', we skip it */
  strncpy(data->guess,&word[1], sizeof(data->guess));

  printf("gamessplugin) Run was performed with GUESS = %s \n",
	  data->guess);

 /* Since this block occurs in the middle of first frame
  * we need to rewind. */
  fseek(data->file, filepos, SEEK_SET);

  return TRUE;
}



/*******************************************************
 *
 * Read the basis set data into hierarchical structures.
 * L-shells are expanded into an S-shell and a P-shell.
 *
 * ******************************************************/
int get_basis(qmdata_t *data) {

  gmsdata *gms = (gmsdata *)data->format_specific_data;

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
     * we skip reading the basis set and hardcode the
     * parameters in tables in VMD. */
    return TRUE;
  }

  /* Search for "ATOMIC BASIS SET" line */
  if (pass_keyline(data->file, "ATOMIC BASIS SET", 
                    "$CONTRL OPTIONS") !=FOUND ) {
    printf("gamessplugin) No basis set found!\n");

    return FALSE;
  }

  /* initialize buffers */
  buffer[0] = '\0';
  for (i=0; i<3; i++) word[i][0] = '\0';
  

  /* skip the next 5 lines */
  eatline(data->file, 5);

  /* Allocate space for the basis for all atoms */
  /* When the molecule is symmetric the actual number atoms with
   * a basis set could be smaller */
  data->basis_set = (basis_atom_t*)calloc(data->numatoms, sizeof(basis_atom_t));


  i = 0; /* basis atom counter */

  do {
    prim_t *prim = NULL;
    char shelltype;
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
          numprim = read_shell_primitives(data, &prim, &shelltype, icoeff, gms->have_pcgamess);

          if (numprim>0) {
            /* make sure we have eiter S, L, P, D, F or G shells */
            if ( (shelltype!='S' && shelltype!='L' && shelltype!='P' && 
                  shelltype!='D' && shelltype!='F' && shelltype!='G') ) {
              printf("gamessplugin) WARNING ... %c shells are not supported \n", shelltype);
            }
            
            /* create new shell */
            if (numshells) {
              shell = (shell_t*)realloc(shell, (numshells+1)*sizeof(shell_t));
            }
            shell[numshells].numprims = numprim;
            /* assign a numeric shell type */
            shell[numshells].type = shelltype_int(shelltype);
            shell[numshells].prim = prim;
            data->num_basis_funcs += numprim;

            /* We split L-shells into one S and one P-shell.
             * I.e. for L-shells we have to go back read the shell again
             * this time using the second contraction coefficients. */
            if (shelltype=='L' && !icoeff) {
              fseek(data->file, filepos, SEEK_SET);
              icoeff++;
            } else if (shelltype=='L' && icoeff) {
              shell[numshells].type = SP_P_SHELL;
              icoeff = 0;  /* reset the counter */
            }

            numshells++;
          }
        } while (numprim);

        /* store shells in atom */
        data->basis_set[i].numshells = numshells;
        data->basis_set[i].shell = shell;

        /* Update total number of basis functions */
        data->num_shells += numshells;
        i++;

        /* go back one line so that we can read the name of the
         * next atom */
        fseek(data->file, filepos, SEEK_SET);

        break;

      case 4:
        /* this is the very end of the basis set */
        if(gms->have_pcgamess){
            if (!strcmp(&word[0][0],"TOTAL")  &&
                !strcmp(&word[1][0],"NUMBER") && 
                !strcmp(&word[2][0],"OF")     &&
                !strcmp(&word[3][0],"SHELLS")) {
              success = 1;
              /* go back one line so that get_basis_stats()
                 can use this line as a keystring. */
              fseek(data->file, filepos, SEEK_SET);
            }
        }
        else {
            if (!strcmp(&word[0][0],"TOTAL")  &&
                !strcmp(&word[1][0],"NUMBER") && 
                !strcmp(&word[2][0],"OF")     &&
                !strcmp(&word[3][0],"BASIS")) {
              success = 1;
              /* go back one line so that get_basis_stats()
                 can use this line as a keystring. */
              fseek(data->file, filepos, SEEK_SET);
            }
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
 * Convert shell type from char to int.
 *
 ************************************************ */
static int shelltype_int(char type) {
  int shelltype;

  switch (type) {
    case 'L':
      /* SP_P shells are assigned in get_basis() */
      shelltype = SP_S_SHELL;
      break;
    case 'S':
      shelltype = S_SHELL;
      break;
    case 'P':
      shelltype = P_SHELL;
      break;
    case 'D':
      shelltype = D_SHELL;
      break;
    case 'F':
      shelltype = F_SHELL;
      break;
    case 'G':
      shelltype = G_SHELL;
      break;
    default:
      shelltype = UNK_SHELL;
      break;
  }

  return shelltype;
}



/******************************************************
 *
 * Populate the flat arrays containing the basis
 * set data.
 *
 ******************************************************/
static int fill_basis_arrays(qmdata_t *data) {
  gmsdata *gms = (gmsdata *)data->format_specific_data;
  int i, j, k;
  int shellcount = 0;
  int primcount = 0;

  float *basis;
  int *num_shells_per_atom;
  int *num_prim_per_shell;
  int *shell_types;
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
    return FALSE;
  }

  shell_types = (int *)calloc(data->num_shells, sizeof(int));
  
  /* make sure memory was allocated properly */
  if (shell_types == NULL) {
    PRINTERR; 
    return FALSE;
  }

  num_shells_per_atom = (int *)calloc(data->num_basis_atoms, sizeof(int));

  /* make sure memory was allocated properly */
  if (num_shells_per_atom == NULL) {
    PRINTERR; 
    return FALSE;
  }

  num_prim_per_shell = (int *)calloc(data->num_shells, sizeof(int));

  /* make sure memory was allocated properly */
  if (num_prim_per_shell == NULL) {
    PRINTERR;
    return FALSE;
  }

  atomicnum_per_basisatom = (int *)calloc(data->num_basis_atoms, sizeof(int));

  /* make sure memory was allocated properly */
  if (atomicnum_per_basisatom == NULL) {
    PRINTERR;
    return FALSE;
  }


  /* store pointers in struct qmdata_t */
  data->basis = basis;
  data->shell_types = shell_types;
  data->num_shells_per_atom = num_shells_per_atom;
  data->num_prim_per_shell = num_prim_per_shell;
  data->atomicnum_per_basisatom = atomicnum_per_basisatom;

  /* Go through all basis set atoms and try to assign the
   * atomic numbers. The basis set atoms are specified by 
   * name strings (the same as in the coordinate section,
   * except for FMO calcs.) and we try to match the names
   * from the two lists. The basis set atom list is symmetry
   * unique while the coordinate atom list is complete.*/
  primcount = 0;
  for (i=0; i<data->num_basis_atoms; i++) {
    int found = 0;

    /* For this basis atom find a matching atom from the
     * coordinate atom list. */
    for(j=0; j<data->numatoms; j++) {
      char basisname[BUFSIZ];
      strcpy(basisname, data->basis_set[i].name);

      /* for FMO calculations we have to strip the "-n" tail
       * of the basis atom name. */
      if (gms->have_fmo) {
        *strchr(basisname, '-') = '\0';
      }

      if (!strcmp(data->atoms[j].type, basisname)) {
        found = 1;
        break;
      }
    }
    if (!found) {
      printf("gamessplugin) WARNING: Couldn't find atomic number for basis set atom %s\n",
             data->basis_set[i].name);
      data->basis_set[i].atomicnum = 0;
      atomicnum_per_basisatom[i] = 0;
    } else {
      /* assign atomic number */
      data->basis_set[i].atomicnum = data->atoms[j].atomicnum;
      atomicnum_per_basisatom[i]   = data->atoms[j].atomicnum;
    }
    num_shells_per_atom[i] = data->basis_set[i].numshells;

    for (j=0; j<data->basis_set[i].numshells; j++) {
      shell_types[shellcount] = data->basis_set[i].shell[j].type;
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
static int read_shell_primitives(qmdata_t *data, prim_t **prim, char *shelltype,
                                 int icoeff, int pcgamess) {
  char buffer[BUFSIZ];
  float exponent = 0.0; 
  float contract[2] = {0.0, 0.0};
  int shell, success;
  int primcounter = 0;
  (*prim) = (prim_t*)calloc(1, sizeof(prim_t));

  do {
    GET_LINE(buffer, data->file);
      if (pcgamess)
        success = sscanf(buffer,"%d %c %*s %f %f %*s %*s %f", &shell,
                       shelltype,
                       &exponent, &contract[0], &contract[1]); 

      else
        success = sscanf(buffer,"%d %c %*s %f %f %f", &shell,
                       shelltype,
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
static int get_traj_frame(qmdata_t *data, qm_atom_t *atoms,
                          int natoms) {
  gmsdata *gms = (gmsdata *)data->format_specific_data;
  qm_timestep_t *cur_ts;
  char buffer[BUFSIZ];
  char word[BUFSIZ];
  int units;
  buffer[0] = '\0';
  word[0]   = '\0';

  printf("gamessplugin) Timestep %d:\n", data->num_frames_read);
  printf("gamessplugin) ============\n");

  fseek(data->file, data->filepos_array[data->num_frames_read], SEEK_SET);

  /* Read the coordinate block */
  if (data->runtype==MOLFILE_RUNTYPE_OPTIMIZE ||
      data->runtype==MOLFILE_RUNTYPE_SADPOINT) {
    goto_keyline(data->file, "COORDINATES OF ALL ATOMS", NULL);
    /* get the units */
    GET_LINE(buffer, data->file);
    sscanf(buffer, " COORDINATES OF ALL ATOMS ARE %s", word);
    units = !strcmp(word, "(BOHR)");
    eatline(data->file, 2);

    if (!get_coordinates(data->file, &data->atoms, units, &natoms)) {
      printf("gamessplugin) Couldn't find coordinates for timestep %d\n", data->num_frames_read);
    }
  }
  else if (data->runtype==MOLFILE_RUNTYPE_SURFACE) {
    if (pass_keyline(data->file, "HAS ENERGY VALUE",
                     "...... END OF ONE-ELECTRON INTEGRALS ......")
        == FOUND) {
      /* Read the coordinate block following 
       * ---- SURFACE MAPPING GEOMETRY ---- */
      int i, n;
      for (i=0; i<natoms; i++) {
        char atname[BUFSIZ];
        float x,y,z;
        GET_LINE(buffer, data->file);
        n = sscanf(buffer,"%s %f %f %f", atname, &x,&y,&z);
        if (n!=4 || strcmp(atname, data->atoms[i].type)) break;
        data->atoms[i].x = x;
        data->atoms[i].y = y;
        data->atoms[i].z = z;
      }
      if (i!=natoms) {
        printf("gamessplugin) Couldn't read surface mapping geometry for timestep %d\n", data->num_frames_read);
      }
    }
    else {
      /* Read the coordinate block following 
       * ATOM      ATOMIC                      COORDINATES (BOHR) */
      goto_keyline(data->file, "ATOM      ATOMIC", NULL);
      /* get the units */
      GET_LINE(buffer, data->file);
      sscanf(buffer, " ATOM      ATOMIC                      COORDINATES %s", word);
      units = !strcmp(word, "(BOHR)");
      eatline(data->file, 1);
      
      if (!get_coordinates(data->file, &data->atoms, units, &natoms)) {
        printf("gamessplugin) Couldn't find coordinates for timestep %d\n", data->num_frames_read);
      }
    }
  }
  /* XXX could merge this with OPTIMIZE/SADPOINT */
  else if (data->runtype==MOLFILE_RUNTYPE_MEX) {
    int numuniqueatoms = natoms;
    goto_keyline(data->file, "COORDINATES OF SYMMETRY UNIQUE ATOMS", NULL);
    /* get the units */
    GET_LINE(buffer, data->file);
    sscanf(buffer, " COORDINATES OF SYMMETRY UNIQUE ATOMS ARE %s", word);
    units = !strcmp(word, "(BOHR)");
    eatline(data->file, 2);
    if (!get_coordinates(data->file, &data->atoms, units, &numuniqueatoms)) {
      printf("gamessplugin) Expanding symmetry unique coordinates for timestep %d\n", data->num_frames_read);

      /* Create images of symmetry unique atoms so that we have
       * the full coordinate set. */
      symmetry_expand(&data->atoms, numuniqueatoms, natoms,
                      data->pointgroup, data->naxis);
    }
  }

  /* For FMO calculations we read the coordinates only
   * because the wavefunctions are printed per fragment
   * and VMD requires that there's a wavefunction present
   * for each atom.
   * A possible workaround would be to pad the wavefunctions
   * accordingly and add a wavefunction for each fragment. */
  if (gms->have_fmo) {
    data->num_frames_read++;
    return TRUE;
  }

  /* get a convenient pointer to the current qm timestep */
  cur_ts = data->qm_timestep + data->num_frames_read;

  /* read the SCF energies */
  if (get_scfdata(data, cur_ts) == FALSE) {
    printf("gamessplugin) Couldn't find SCF iterations for timestep %d\n",
           data->num_frames_read);
  }

  /* Try reading canonical alpha/beta wavefunction */
  check_add_wavefunctions(data, cur_ts);


  /* Read population analysis (Mulliken and Lowdin charges)
   * only if wasn't read already while parsing the final
   * property section. Otherwise we would potentially 
   * overwrite the data with empty fields. */
  if (!cur_ts->have_mulliken &&
      get_population(data, cur_ts)) {
    printf("gamessplugin) Mulliken/Loewdin charges found\n");
  }

  if (data->citype==CI_GUGA) {
    if (pass_keyline(data->file, "CI DENSITY MATRIX AND NATURAL ORBITALS",
                       "GRADIENT (HARTREE/BOHR)")) {
      int i, numstates=0, state;
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
        sscanf(buffer," CI EIGENSTATE %d %*s %*s %*s %f", &state, &cienergy);
        printf("gamessplugin) CI energy[%d] = %f\n", state-1, cienergy);

        wave_ci = add_wavefunction(cur_ts);

        if (get_wavefunction(data, cur_ts, wave_ci) == FALSE) {
          del_wavefunction(cur_ts);
          break;
        }
        else {
          /* canon =-1;*/
          wave_ci->exci = state-1;
          wave_ci->energy = cienergy;
          wave_ci->mult = data->multiplicity;

          printf("gamessplugin) Found %d CI natural orbitals for excited state %d, mult=%d\n",
                 wave_ci->num_orbitals, state-1, wave_ci->mult);
        }
      }
    }
  }
  else if (data->citype==CI_CIS) {
    if (pass_keyline(data->file,
                     "USING DAVIDSON ALGORITHM TO FIND CIS EIGENVALUES",
                     NULL)) {
      int i, numstates=0, state;
      qm_wavefunction_t *wave_ci;
      float *state_energies, *state_spinquant;
      goto_keyline(data->file, "NUMBER OF STATES REQUESTED", NULL);
      GET_LINE(buffer, data->file);
      trimleft(buffer);
      sscanf(buffer, " NUMBER OF STATES REQUESTED = %d", &numstates);
      printf("gamessplugin) Number of CIS states = %d\n", numstates);

      /* For CIS only the wavefunction for the excited state 
       * (specified by IROOT in $CIS group) will be printed.
       * Here we read in the energies for all states, and store 
       * the energy for the selected state in its wavefunction */
      state_energies  = calloc(numstates+1, sizeof(float));
      state_spinquant = calloc(numstates+1, sizeof(float));
      goto_keyline(data->file, "RHF REFERENCE ENERGY  =", NULL);
      GET_LINE(buffer, data->file);
      trimleft(buffer);
      sscanf(buffer, " RHF REFERENCE ENERGY  = %f", &state_energies[0]);
      state_spinquant[0] = 1.f;

      for (i=1; i<=numstates; i++) {
        goto_keyline(data->file, "EXCITED STATE", NULL);
        GET_LINE(buffer, data->file);
        trimleft(buffer);
        sscanf(buffer, " EXCITED STATE %*d  ENERGY= %f  S = %f",
               &state_energies[i], &state_spinquant[i]);
      }

      goto_keyline(data->file,
                   "CIS NATURAL ORBITAL OCCUPATION NUMBERS FOR EXCITED STATE",
                   NULL);
      GET_LINE(buffer, data->file);
      trimleft(buffer);
      sscanf(buffer,
             " CIS NATURAL ORBITAL OCCUPATION NUMBERS FOR EXCITED STATE %d",
             &state);
      
      wave_ci = add_wavefunction(cur_ts);

      if (get_wavefunction(data, cur_ts, wave_ci) == FALSE) {
        del_wavefunction(cur_ts);
      }
      else {
        wave_ci->exci = state;
        wave_ci->energy = state_energies[state];
        wave_ci->mult = 2*(int)state_spinquant[state]+1;
        printf("gamessplugin) Found %d CIS natural orbitals for excited state %d\n",
               wave_ci->num_orbitals, state);
      }

      free(state_energies);
      free(state_spinquant);
    }
  }


  /* Read the energy gradients (=forces on atoms) */
  if (get_gradient(data, cur_ts)) {
    printf("gamessplugin) Energy gradient found.\n");
  }


  /* If this is the last frame of the trajectory and the file
   * wasn't truncated and the program didn't terminate
   * abnormally then read the final wavefunction. */
  if ((data->runtype == MOLFILE_RUNTYPE_OPTIMIZE ||
       data->runtype == MOLFILE_RUNTYPE_SADPOINT) &&
      (data->num_frames_read+1 == data->num_frames &&
       (data->status == MOLFILE_QMSTATUS_UNKNOWN || 
        data->status == MOLFILE_QMSTATUS_OPT_CONV ||
        data->status == MOLFILE_QMSTATUS_OPT_NOT_CONV))) {

    /* We need to jump over the end of the trajectory because 
     * this is also the keystring for get_wavefunction() to
     * bail out. */
    if (data->status == MOLFILE_QMSTATUS_OPT_CONV || 
        data->status == MOLFILE_QMSTATUS_OPT_NOT_CONV) {
      fseek(data->file, data->end_of_traj, SEEK_SET);
    }

    /* Try to read final wavefunction and orbital energies
     * A preexisting canonical wavefunction for this timestep
     * with the same characteristics (spin, exci, info) will
     * be overwritten by the final wavefuntion if it has more
     * orbitals. */
    check_add_wavefunctions(data, cur_ts);
  }


  /* For MCSCF optimized orbitals no occupancies are given
   * but since their occupancies are identical to the ones
   * from natural orbitals we can use those. The natural 
   * orbitals are always listed right before the optimized
   * ones so we simply copy the data over. */
  if (cur_ts->numwave>=2 &&
      cur_ts->wave[cur_ts->numwave-1].type==MOLFILE_WAVE_MCSCFOPT &&
      cur_ts->wave[cur_ts->numwave-2].type==MOLFILE_WAVE_MCSCFNAT) {
    int i;
    qm_wavefunction_t *waveopt = &cur_ts->wave[cur_ts->numwave-1];
    qm_wavefunction_t *wavenat = &cur_ts->wave[cur_ts->numwave-2];
    waveopt->orb_occupancies = (float *)calloc(waveopt->num_orbitals,
                                               sizeof(float));
    /* Only the core and active natural orbitals are listed. 
     * We copy the occupancies for these orbitals and pad the
     * rest with zeros. */
    for (i=0; i<wavenat->num_orbitals; i++) {
      waveopt->orb_occupancies[i] = wavenat->orb_occupancies[i];
    }
    for (i=wavenat->num_orbitals; i<waveopt->num_orbitals; i++) {
      waveopt->orb_occupancies[i] = 0.f;
    }
    waveopt->has_occup = TRUE;
  }


  data->num_frames_read++;

  return TRUE;
}


/* Analyze the trajectory.
 * Read the parameters controlling geometry search and
 * find the end of the trajectory, couinting the frames
 * on the way. Store the filepointer for the beginning of
 * each frame in *filepos_array. */
static int analyze_traj(qmdata_t *data, gmsdata *gms) {
  char buffer[BUFSIZ], nserch[BUFSIZ];
  char *line;
  long filepos;
  filepos = ftell(data->file);

  data->filepos_array = (long* )calloc(1, sizeof(long ));

  if (data->runtype==MOLFILE_RUNTYPE_OPTIMIZE ||
      data->runtype==MOLFILE_RUNTYPE_SADPOINT) {
    pass_keyline(data->file,
                   "PARAMETERS CONTROLLING GEOMETRY SEARCH", NULL);
    eatline(data->file, 2);

    GET_LINE(buffer, data->file);
    sscanf(buffer, "NSTEP  = %d", &data->max_opt_steps);
    eatline(data->file, 3);
    GET_LINE(buffer, data->file);
    sscanf(buffer, "OPTTOL = %f", &data->opt_tol);

    /* The $STATP options are followed by the coordinates 
     * but we can skip them here because we rewind after
     * get_guess_options() and try to read them in
     * get_traj_frame(). */
  }
  else if (data->runtype==MOLFILE_RUNTYPE_SURFACE) {
    if (pass_keyline(data->file,
                     "POTENTIAL SURFACE MAP INPUT", NULL)) {
      
      int coord1[2];
      int mplevel1=-1, mplevel2=-1, nstep1;
      float origin1, disp1;
      char runtype1[BUFSIZ], runtype2[BUFSIZ];
      char scftype1[BUFSIZ], scftype2[BUFSIZ];
      char dfttype1[BUFSIZ], dfttype2[BUFSIZ];
      char citype1[BUFSIZ],  citype2[BUFSIZ];
      char cctype1[BUFSIZ],  cctype2[BUFSIZ];
      char *tmp;
      int n;
        
      eatline(data->file, 1);

      GET_LINE(buffer, data->file);
      n=sscanf(buffer, " JOB 1 IS RUNTYP=%s SCFTYP=%s CITYP=%s",
               runtype1, scftype1, citype1);
      if (n==3) {
        GET_LINE(buffer, data->file);
        sscanf(buffer, " MPLEVL= %d CCTYP=%s, DFTTYP=%s\n",
               &mplevel1, dfttype1, cctype1);
        GET_LINE(buffer, data->file);
      }
      n=sscanf(buffer, " JOB 2 IS RUNTYP=%s SCFTYP=%s CITYP=%s",
               runtype2, scftype2, citype2);
      if (n==3) {
        GET_LINE(buffer, data->file);
        sscanf(buffer, " MPLEVL= %d CCTYP=%s, DFTTYP=%s\n",
               &mplevel2, dfttype2, cctype2);
        GET_LINE(buffer, data->file);
      }

      sscanf(buffer, " COORD 1 LYING ALONG ATOM PAIR %d %d",
             coord1, coord1+1);
      GET_LINE(buffer, data->file);
      tmp = strstr(buffer, "ORIGIN=") + 7;
      sscanf(tmp, "%f", &origin1);
      tmp = strstr(buffer, "DISPLACEMENT=") + 13;
      sscanf(tmp, "%f", &disp1);
      tmp = strstr(buffer, "AND") + 3;
      sscanf(tmp, "%d STEPS.", &nstep1);
      printf("gamessplugin) origin=%f, displacement=%f nstep=%d\n", origin1, disp1, nstep1);
    }
  }
  else if (data->runtype==MOLFILE_RUNTYPE_MEX) {
    char scftype1[BUFSIZ];
    char scftype2[BUFSIZ];
    rewind(data->file);
    if (!pass_keyline(data->file, "$MEX OPTIONS", NULL)) {
      printf("gamessplugin) No $MEX OPTIONS found!\n");
      return FALSE;
    }
    eatline(data->file, 2);
    GET_LINE(buffer, data->file);
    sscanf(strstr(buffer, "SCF1    =")+7, "%s", scftype1);
    sscanf(strstr(buffer, "SCF2   =")+7, "%s", scftype2);
    printf("gamessplugin) MEX SCF1=%s SCF2=%s\n", scftype1, scftype2);

  }
  else {
    /* We have just one frame */
    data->num_frames = 1;
    pass_keyline(data->file, "1 ELECTRON INTEGRALS",
                 "ENERGY COMPONENTS");
    data->filepos_array[0] = ftell(data->file);

    /* Check wether SCF has converged */
    if (pass_keyline(data->file,
                     "SCF IS UNCONVERGED, TOO MANY ITERATIONS",
                     "ENERGY COMPONENTS")==FOUND) {
      printf("gamessplugin) SCF IS UNCONVERGED, TOO MANY ITERATIONS\n");
      data->status = MOLFILE_QMSTATUS_SCF_NOT_CONV;
    } else {
      data->status = MOLFILE_QMSTATUS_OPT_CONV;
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
  data->status = MOLFILE_QMSTATUS_UNKNOWN;

  while (1) {
    if (!fgets(buffer, sizeof(buffer), data->file)) break;
    line = trimleft(buffer);

    /* at this point we have to distinguish between
     * pre="27 JUN 2005 (R2)" and "27 JUN 2005 (R2)"
     * versions since the output format for geometry
     * optimizations has changed */
    if (gms->version==1) {
      strcpy(nserch, "1NSERCH=");
    }
    else if (gms->version==2) {
      strcpy(nserch, "BEGINNING GEOMETRY SEARCH POINT NSERCH=");
    }

    if (strstr(line, nserch) ||
        strstr(line, "---- SURFACE MAPPING GEOMETRY") ||
        strstr(line, "MINIMUM ENERGY CROSSING POINT SEARCH") ||
        (data->runtype==MOLFILE_RUNTYPE_MEX && strstr(line, "NSERCH=")==line)) {
      printf("gamessplugin) %s", line);

      if (data->num_frames > 0) {
        data->filepos_array = (long*)realloc(data->filepos_array,
                                (data->num_frames+1)*sizeof(long));
      }
      data->filepos_array[data->num_frames] = ftell(data->file);
      if (data->runtype==MOLFILE_RUNTYPE_SURFACE) {
        int ret = goto_keyline(data->file,
                               "ATOM      ATOMIC", "HAS ENERGY VALUE",
                               "---- SURFACE MAPPING GEOMETRY ----", NULL);
        if (ret>0 && ret<3 &&
            (have_keyline(data->file, "...... END OF ONE-ELECTRON INTEGRALS ......",
                          "---- SURFACE MAPPING GEOMETRY ----") ||
             have_keyline(data->file, "... DONE WITH POTENTIAL SURFACE SCAN",
                          "---- SURFACE MAPPING GEOMETRY ----"))) {
          data->num_frames++;          
        }
      }
      else if (pass_keyline(data->file, "COORDINATES OF",
                            "BEGINNING GEOMETRY SEARCH POINT NSERCH=")==FOUND)
      {
        /* Make sure that we have at least a complete coordinate
           block in order to consider this a new frame. */
        if (have_keyline(data->file, "INTERNUCLEAR DISTANCES",
                         "1 ELECTRON INTEGRALS") ||
            have_keyline(data->file, "1 ELECTRON INTEGRALS",
                         "BEGINNING GEOMETRY SEARCH POINT NSERCH=")) {
          data->num_frames++;
        }
      }
    }
    else if (strstr(line, "***** EQUILIBRIUM GEOMETRY LOCATED") ||
             strstr(line, "... DONE WITH POTENTIAL SURFACE SCAN")) {
      printf("gamessplugin) ==== End of trajectory (%d frames) ====\n",
             data->num_frames);
      data->status = MOLFILE_QMSTATUS_OPT_CONV;
      break;
    }
    else if (strstr(line, "***** FAILURE TO LOCATE STATIONARY POINT,")) {
      printf("gamessplugin) %s\n", line);
      if (strstr(strchr(line, ','), "SCF HAS NOT CONVERGED")) {
        data->status = MOLFILE_QMSTATUS_SCF_NOT_CONV;
        break;
      }
      else if (strstr(strchr(line, ','), "TOO MANY STEPS TAKEN")) {
        data->status = MOLFILE_QMSTATUS_OPT_NOT_CONV;
        break;
      }
    }
  }
  
  data->end_of_traj = ftell(data->file);
  fseek(data->file, filepos, SEEK_SET);

  if (data->status == MOLFILE_QMSTATUS_UNKNOWN) {
    /* We didn't find any of the regular key strings,
     * the run was most likely broken off and we have an
     * incomplete file. */
    data->status = MOLFILE_QMSTATUS_FILE_TRUNCATED;
  }


  /* Allocate memory for all frames */
  data->qm_timestep = (qm_timestep_t *)calloc(data->num_frames,
                                              sizeof(qm_timestep_t));
  memset(data->qm_timestep, 0, data->num_frames*sizeof(qm_timestep_t));


  if (data->status == MOLFILE_QMSTATUS_SCF_NOT_CONV ||
      data->status == MOLFILE_QMSTATUS_FILE_TRUNCATED) {
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
static int get_scfdata(qmdata_t *data, qm_timestep_t *ts) {
  char buffer[BUFSIZ];
  char word[3][BUFSIZ];
  long filepos;
  int i, epos = -1;
  int numread, numiter=0, dum, dum2;
  char *line;
  float dumf;
  filepos = ftell(data->file);

  for (i=0; i<3; i++) word[i][0] = '\0';

  if (!goto_keyline(data->file, "ITER EX", "ITER     TOTAL",
                     "ITER    TOTAL", NULL)) {
    fseek(data->file, filepos, SEEK_SET);
    ts->num_scfiter = 0;
    return FALSE;
  }

  /* determine in which column the energy is stored */
  GET_LINE(buffer, data->file);
  numread = sscanf(buffer, "%*s %s %s %s",
                   &word[0][0], &word[1][0], &word[2][0]);
  for (i=0; i<numread; i++) {
    if (!strcmp(&word[i][0], "TOTAL")) epos = i+1;
  }
   
  if (epos<0) {
    fseek(data->file, filepos, SEEK_SET);
    ts->num_scfiter = 0;
    return FALSE;
  }

  /* store current file position since we first have to count
   * the iterations */
  filepos = ftell(data->file);

  /* read until the next blank line and count the iterations */
  do {
    GET_LINE(buffer, data->file);
    line = trimleft(buffer);
    numread = sscanf(line,"%d %d %*d %*f", &dum, &dum2);
    if (numread==2) numiter++;
  } while (strlen(line)>2);

  printf("gamessplugin) %d SCF iterations\n", numiter);

  /* go back and read energies */
  fseek(data->file, filepos, SEEK_SET);
  

  /* allocate memory for scfenergy array */
  ts->scfenergies = (double *)calloc(numiter,sizeof(double));
  
  i=0;
  do {
    GET_LINE(buffer, data->file);
    line = trimleft(buffer);
    numread = sscanf(line,"%d %f %*i %*f", &dum, &dumf);
    if (numread==2) {
      switch (epos) {
      case 1:
        sscanf(buffer,"%*d %lf", ts->scfenergies+i);
        break;
      case 2:
        sscanf(buffer,"%*d %*d %lf", ts->scfenergies+i);
        break;
      case 3:
        sscanf(buffer,"%*d %*d %*d %lf", ts->scfenergies+i);
        break;
      }
      i++;
    }
  } while (strlen(line)>2);

#if 0
  for (i=0; i<numiter; i++) {
    printf("scfenergies[%d] = %f\n", i, ts->scfenergies[i]);
  }
#endif

  ts->num_scfiter = numiter;
  
  return TRUE;
}


/*********************************************************
 *
 * Reads a set of wavefunctions for the current timestep.
 * These are typically the alpha and beta spin wavefunctions
 * or the MCSCF natural and optimized orbitals or the GVB
 * canonical orbitals and geminal pairs.
 *
 **********************************************************/
static int check_add_wavefunctions(qmdata_t *data,
                                   qm_timestep_t *ts) {
  qm_wavefunction_t *wavef;
  int i, n=1;

  if (data->scftype==MOLFILE_SCFTYPE_UHF || 
      data->scftype==MOLFILE_SCFTYPE_GVB ||
      data->scftype==MOLFILE_SCFTYPE_MCSCF) {
    /* Try to read second wavefunction
     * (spin beta or GI orbitals or MCSCF optimized orbs) */
    n = 2;
  }

  for (i=0; i<n; i++) {
    /* Allocate memory for new wavefunction */
    wavef = add_wavefunction(ts);

    /* Try to read wavefunction and orbital energies */
    if (get_wavefunction(data, ts, wavef) == FALSE) {
      /* Free the last wavefunction again. */
      del_wavefunction(ts);
#ifdef DEBUGGING
      printf("gamessplugin) No canonical wavefunction present for timestep %d\n", data->num_frames_read);
#endif
      break;

    } else {
      char action[32];
      char spinstr[32];
      strcpy(spinstr, "");
      if (data->scftype==MOLFILE_SCFTYPE_UHF) {
        if (wavef->spin==SPIN_BETA) {
          strcat(spinstr, "spin  beta, ");
        } else {
          strcat(spinstr, "spin alpha, ");
        }
      }
      
      /* The last SCF energy is the energy of this electronic state */
      if (ts->scfenergies) {
        wavef->energy = ts->scfenergies[ts->num_scfiter-1];
      } else {
        wavef->energy = 0.f;
      }
      
      /* Multiplicity */
      wavef->mult = data->multiplicity;
      

      /* String telling wether wavefunction was added, updated
       * or ignored. */
      strcpy(action, "added");

      /* If there exists a canonical wavefunction of the same spin
       * we'll replace it */
      if (ts->numwave>1 && wavef->type==MOLFILE_WAVE_CANON) {
        int i, found =-1;
        for (i=0; i<ts->numwave-1; i++) {
          if (ts->wave[i].type==wavef->type &&
              ts->wave[i].spin==wavef->spin &&
              ts->wave[i].exci==wavef->exci &&
              !strncmp(ts->wave[i].info, wavef->info, MOLFILE_BUFSIZ)) {
            found = i;
            break;
          }
        }
        if (found>=0) {
          /* If the new wavefunction has more orbitals we 
           * replace the old one for this step. */
          if (wavef->num_orbitals > 
              ts->wave[found].num_orbitals) {
            /* Replace existing wavefunction for this step */
            replace_wavefunction(ts, found);
            sprintf(action, "%d updated", found);
          } else {
            /* Delete last wavefunction again */
            del_wavefunction(ts);
            sprintf(action, "matching %d ignored", found);
          }
          wavef = &ts->wave[ts->numwave-1];
        }
      }

      printf("gamessplugin) Wavefunction %s (%s):\n", action, wavef->info);
      printf("gamessplugin)   %d orbitals, %sexcitation %d, multiplicity %d\n",
             wavef->num_orbitals, spinstr, wavef->exci, wavef->mult);
    }
  }

  return i;
}


/*********************************************************
 *
 * Finds the next wavefunction, determines its type by
 * analyzing the keystring and reads in the wavefunction
 * coefficients.
 *
 **********************************************************/
static int get_wavefunction(qmdata_t *data, qm_timestep_t *ts,
                            qm_wavefunction_t *wf)
{
  float *orb_enocc;
  float *wave_coeff;
  char buffer[BUFSIZ];
  char word[6][BUFSIZ];
  int num_orbitals = 0;
  int i = 0, num_values = 0;
  long filepos;
  char *line;
  int have_orbenocc = 0;
  int n[5];

  buffer[0] = '\0';
  for (i=0; i<6; i++) word[i][0] = '\0';

  if (wf == NULL) {
    PRINTERR;	    
    return FALSE;
  }

  wf->has_occup = FALSE;
  wf->has_orben = FALSE;
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
                  -11.0297    -0.9121    -0.5205    -0.5205    -0.5205  <<-- orbital energies (or occupancies)                     A          A          A          A          A   
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

    line = trimleft(trimright(buffer));

    if      (!strcmp(line, "----- ALPHA SET -----")) {
      wf->type = MOLFILE_WAVE_CANON;
      strncpy(wf->info, "canonical", MOLFILE_BUFSIZ);
      pass_keyline(data->file, "EIGENVECTORS", NULL);
    }
    else if (!strcmp(line, "----- BETA SET -----")) {
      wf->type = MOLFILE_WAVE_CANON;
      wf->spin = SPIN_BETA;
      strncpy(wf->info, "canonical", MOLFILE_BUFSIZ);
      pass_keyline(data->file, "EIGENVECTORS", NULL);
    }
    else if (!strcmp(line, "****** BETA ORBITAL LOCALIZATION *****")) {
      wf->spin = SPIN_BETA;
    }
    else if (!strcmp(line, "EIGENVECTORS")) {
      wf->type = MOLFILE_WAVE_CANON;
      strncpy(wf->info, "canonical", MOLFILE_BUFSIZ);
    }
    else if (!strcmp(line, "MOLECULAR ORBITALS")) {
      wf->type = MOLFILE_WAVE_CANON;
      strncpy(wf->info, "canonical", MOLFILE_BUFSIZ);
    }
    else if (!strcmp(line, "THE BOYS LOCALIZED ORBITALS ARE")) {
      wf->type = MOLFILE_WAVE_BOYS;
      strncpy(wf->info, "Boys localized", MOLFILE_BUFSIZ);
    }
    else if (!strcmp(line, "THE PIPEK-MEZEY POPULATION LOCALIZED ORBITALS ARE")) {
      wf->type = MOLFILE_WAVE_PIPEK;
      strncpy(wf->info, "Pipek-Mezey localized", MOLFILE_BUFSIZ);
    }
    else if (!strcmp(line, "EDMISTON-RUEDENBERG ENERGY LOCALIZED ORBITALS")) {
      wf->type = MOLFILE_WAVE_RUEDEN;
      strncpy(wf->info, "Ruedenberg localized", MOLFILE_BUFSIZ);
    }
    else if (!strcmp(line, "GI ORBITALS")) {
      wf->type = MOLFILE_WAVE_GEMINAL;
      strncpy(wf->info, "GVB geminal pairs", MOLFILE_BUFSIZ);
    }
    else if (!strcmp(line, "MCSCF NATURAL ORBITALS")) {
      wf->type = MOLFILE_WAVE_MCSCFNAT;
      strncpy(wf->info, "MCSCF natural orbitals", MOLFILE_BUFSIZ);
    }
    else if (!strcmp(line, "MCSCF OPTIMIZED ORBITALS")) {
      wf->type = MOLFILE_WAVE_MCSCFOPT;
      strncpy(wf->info, "MCSCF optimized orbitals", MOLFILE_BUFSIZ);
    }
    else if (!strcmp(line, "NATURAL ORBITALS IN ATOMIC ORBITAL BASIS")) {
      wf->type = MOLFILE_WAVE_CINATUR;
      strncpy(wf->info, "CI natural orbitals", MOLFILE_BUFSIZ);
    }
    else if (!strcmp(line, "CIS NATURAL ORBITALS")) {
      wf->type = MOLFILE_WAVE_CINATUR;
      strncpy(wf->info, "CIS natural orbitals", MOLFILE_BUFSIZ);
    }
    

    else if (!strcmp(line, "-MCHF- NATURAL ORBITALS")) {
      wf->type = MOLFILE_WAVE_MCSCFNAT;
      strncpy(wf->info, "MCSCF natural orbitals", MOLFILE_BUFSIZ);
    }
    else if (!strcmp(line, "-MCHF- OPTIMIZED ORBITALS")) {
      wf->type = MOLFILE_WAVE_MCSCFOPT;
      strncpy(wf->info, "MCSCF optimized orbitals", MOLFILE_BUFSIZ);
    }
    else if (!strcmp(line, "ZERO-ORDER QDPT NATURAL ORBITALS")){
      

    }

  } while(wf->type==MOLFILE_WAVE_UNKNOWN &&
          strcmp(line, "ENERGY COMPONENTS") &&
          strcmp(line, "***** EQUILIBRIUM GEOMETRY LOCATED *****") &&
          strcmp(line, "**** THE GEOMETRY SEARCH IS NOT CONVERGED! ****"));

  /* If we reach the last line of the rhf section without finding 
   * one of the keywords marking the beginning of a wavefunction
   * table then we return.*/
  if (wf->type==MOLFILE_WAVE_UNKNOWN) {
#ifdef DEBUGGING
    printf("gamessplugin) get_wavefunction(): No wavefunction found!\n");
#endif
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

  /* orbital energies/occupancies */  
  orb_enocc = (float *)calloc(data->wavef_size, sizeof(float));

  if (orb_enocc == NULL) {
    free(orb_enocc);
    PRINTERR; 
    return FALSE;
  }


  /* store the coeficient pointer */
  wf->wave_coeffs  = wave_coeff;

  /* depending on the wavefunction type the line after the
     orbital index stores the orbital occupancies or energies */
  if (wf->type == MOLFILE_WAVE_CINATUR  ||
      wf->type == MOLFILE_WAVE_MCSCFNAT) {
    wf->orb_occupancies = orb_enocc;
    wf->has_occup = TRUE;
  } else {
    wf->orb_energies    = orb_enocc;
    wf->has_orben = TRUE;
  }

  /* skip the next line which here is typically "-------" */
  eatline(data->file, 1);


  while (1) {
    int nr, over=0;
    float coeff[5], enocc[5];

    if (wf->type == MOLFILE_WAVE_GEMINAL) {
      /* Skip over "PAIR x" header line */
      pass_keyline(data->file, "PAIR ", NULL);
    }

    eatwhitelines(data->file);
    filepos = ftell(data->file);

    /* Parse the orbital indexes */
    GET_LINE(buffer, data->file);
    num_values = sscanf(buffer, "%d %d %d %d %d",
                          &n[0], &n[1], &n[2], &n[3], &n[4]);

    /* If there are no orbital indexes then this must be the
     * end of the wavefunction coefficient table. */
    if (!num_values) {
      fseek(data->file, filepos, SEEK_SET);
      break;
    }

    eatwhitelines(data->file);

    /* Read first line of orbital energies/occupancies */
    filepos = ftell(data->file);
    GET_LINE(buffer, data->file);
    have_orbenocc = sscanf(buffer,"%f %f %f %f %f", &enocc[0],
                        &enocc[1], &enocc[2], &enocc[3], &enocc[4]);

    /* Make sure this is not the first line containing coeffs */
    nr = sscanf(buffer, " 1 %*s 1 %*s %f %f %f %f %f",
               &coeff[0], &coeff[1], &coeff[2], &coeff[3], &coeff[4]);
    if (nr==num_values) have_orbenocc = 0;

    if (have_orbenocc) {
      /* store the orbital energies in the appropriate arrays 
       * read them until we encounter an empty string */
      for(i=0; i<num_values; i++) {
        orb_enocc[i] = enocc[i];
      }
      

      /* If we are in the first block we have to distinguish 
         between energies and occupancies */
      if (wf->type  == MOLFILE_WAVE_MCSCFNAT &&
          orb_enocc == wf->orb_occupancies   &&
          enocc[0] <= 0.f) {
        wf->orb_occupancies = NULL;
        wf->has_occup = FALSE;
        wf->orb_energies    = orb_enocc;
        wf->has_orben = TRUE;
      }

      /* increase orbital energy pointer */
      orb_enocc = orb_enocc+5;
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

      nr = sscanf(buffer, " 1 %*s 1 %*s %f %f %f %f %f",
             &coeff[0], &coeff[1], &coeff[2], &coeff[3], &coeff[4]);
      if (nr==num_values) break;
      filepos = ftell(data->file);
    }
    fseek(data->file, filepos, SEEK_SET);

    if (over) break;


    /* Read the wave function coefficient block for up to 5
     * orbitals per line. */
    if (!read_coeff_block(data->file, data->wavef_size,
                          wave_coeff, data->angular_momentum)) {
      printf("gamessplugin) Wavefunction coefficient block truncated or ill formatted!\n");
      data->status = MOLFILE_QMSTATUS_FILE_TRUNCATED;
      return FALSE;
    }


    /* move wavefunction pointer to start of next five orbitals */
    if (wf->type == MOLFILE_WAVE_GEMINAL) {
      wave_coeff = wave_coeff + 2*data->wavef_size;
    } else {
      wave_coeff = wave_coeff + 5*data->wavef_size;
    }
  }


  if (!num_orbitals) {
    printf("gamessplugin) No orbitals in wavefunction!\n");
    return FALSE;
  }

  /* resize the array to the actual number of read orbitals */
  if (data->wavef_size!=num_orbitals) {

    if (wf->has_occup) {
      wf->orb_occupancies = (float *)realloc(wf->orb_occupancies,
               num_orbitals*sizeof(float));
    }
    if (wf->has_orben) {
      wf->orb_energies = (float *)realloc(wf->orb_energies,
               num_orbitals*sizeof(float));
    }

    wf->wave_coeffs  = (float *)realloc(wf->wave_coeffs, data->wavef_size*
					num_orbitals*sizeof(float)); 
  }

  /* In case MCSCF natural orbitals are present, then GAMESS
     prints the orbital energy for the core orbitals and the
     occupancy for the other orbitals. We zero out the non-core
     energies to prevent confusion. The orbital occupancies 
     are read separately elsewhere. */
  if (wf->type == MOLFILE_WAVE_MCSCFNAT &&
      wf->has_orben == TRUE ) {
    wf->orb_occupancies = (float *)calloc(num_orbitals, sizeof(float));

    for (i=0; i<data->mcscf_num_core; i++) {
      wf->orb_occupancies[i] = 2.f;
    }

    for (i=data->mcscf_num_core; i<num_orbitals; i++) {
      wf->orb_occupancies[i] = wf->orb_energies[i];
      wf->orb_energies[i] = 0.f;
    }

    wf->has_occup = TRUE;
  }


  /* store the number of orbitals read in */
  wf->num_orbitals  = num_orbitals;

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

  /* Read a line with coefficients for up to 5 orbitals
   * for each cartesian basis function. */
  for (i=0; i<wavef_size; i++) {
    char type[BUFSIZ];
    float coeff[5];
    int num_values = 0;
    
    GET_LINE(buffer, file);
    
    /* read in the wavefunction coefficients for 5
     * orbitals at a time line by line */
    num_values = sscanf(buffer,"%*5i%*4s%*2i%4s %f %f %f %f %f", 
                        type, &coeff[0], &coeff[1], &coeff[2],
                        &coeff[3], &coeff[4]);
    
    if (num_values==0) {
      /* The file must have been truncated! */
      truncated = 1;
      break;
    }

    angular_momentum_expon(&angular_momentum[3*i], type);
   
    /* Each orbital has data->wavef_size entries, 
     * hence we have to use this number as offset when storing 
     * them in groups of five. */
    for (j=0 ; j<num_values-1; j++) {
      wave_coeff[j*wavef_size+i] = coeff[j];
    }
  }
  
  if (truncated) return 0;
  
  return 1;
}

/* Read the population analysis section.
 * Currently we parse only the Mulliken and Lowdin charges
 * but we might want to add support for population analysis. */
static int get_population(qmdata_t *data, qm_timestep_t *ts) {
  int i;
  char buffer[BUFSIZ];
  long filepos;
  ts->have_mulliken = FALSE;
  ts->have_lowdin   = FALSE;
  filepos = ftell(data->file);

  if (pass_keyline(data->file,
                     "TOTAL MULLIKEN AND LOWDIN ATOMIC POPULATIONS",
                     "NSERCH=") != FOUND) {
    fseek(data->file, filepos, SEEK_SET);
    return FALSE;
  }

  /* Read Mulliken charges if present */
  ts->mulliken_charges = 
    (double *)calloc(data->numatoms, sizeof(double));

  if (!ts->mulliken_charges) {
    PRINTERR; 
    return FALSE;
  }

  ts->lowdin_charges = 
    (double *)calloc(data->numatoms, sizeof(double));

  if (!ts->lowdin_charges) {
    free(ts->mulliken_charges);
    ts->mulliken_charges = NULL;
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
    if (n!=4) {
      free(ts->mulliken_charges);
      free(ts->lowdin_charges);
      ts->mulliken_charges = NULL;
      ts->lowdin_charges   = NULL;
      return FALSE;
    }
    ts->mulliken_charges[i] = mullcharge;
    ts->lowdin_charges[i]   = lowcharge;
  }

  if (i!=data->numatoms) {
    free(ts->mulliken_charges);
    free(ts->lowdin_charges);
    ts->mulliken_charges = NULL;
    ts->lowdin_charges   = NULL;
    return FALSE;
  }

  ts->have_mulliken = TRUE;
  ts->have_lowdin   = TRUE;
  return TRUE;
}


/* Read ESP charges.
 * XXX Right now we don't distinguish between different type of
 * ESP-style charges (CHELPG, CONNOLLY, GEODESIC). 
 * This could be solved by reading in the PTSEL keyword in
 * the $PDC group. */
static int get_esp_charges(qmdata_t *data) {
  int i;
  char buffer[BUFSIZ];
  long filepos;
  /* Store charges in last timestep */
  qm_timestep_t *ts = &data->qm_timestep[data->num_frames-1];

  ts->have_esp = FALSE;
  filepos = ftell(data->file);

  if (pass_keyline(data->file,
           "ATOM                CHARGE    E.S.D.",
           "...... END OF PROPERTY EVALUATION ") != FOUND) {
    fseek(data->file, filepos, SEEK_SET);
    return FALSE;
  }

  /* Read ESP charges if present */
  ts->esp_charges = 
    (double *)calloc(data->numatoms, sizeof(double));

  if (ts->esp_charges == NULL) {
    PRINTERR; 
    return FALSE;
  }

  eatline(data->file, 1);

  for (i=0; i<data->numatoms; i++) {
    int n;
    double charge;
    GET_LINE(buffer, data->file);
    n = sscanf(buffer,"%*s %lf ", &charge);
    if (n!=1) return FALSE;
    ts->esp_charges[i] = charge;
  }

  if (i!=data->numatoms) {
    
    return FALSE;
  }

  ts->have_esp = TRUE;
  return TRUE;
}


/* Read the energy gradient (=force) for each atom */
static int get_gradient(qmdata_t *data, qm_timestep_t *ts) {
  int numgrad=0;
  int numread;
  char buffer[BUFSIZ];
  long filepos;

  buffer[0] = '\0';

  /* remember position in order to go back if no forces were found */
  filepos = ftell(data->file);

  /* look for GRADIENT section */
  if (goto_keyline(data->file, "GRADIENT (HARTREE",
                "***** EQUILIBRIUM GEOMETRY LOCATED", 
                " BEGINNING GEOMETRY SEARCH", NULL) != FOUND) {
    fseek(data->file, filepos, SEEK_SET);
    return FALSE;
  }

  eatline(data->file, 4);

  ts->gradient = (float *)calloc(3*data->numatoms, sizeof(float));

  if (ts->gradient == NULL) {
    PRINTERR;	    
    fseek(data->file, filepos, SEEK_SET);
    return FALSE;
  }

  /* read the gradient table */
  do {
    int i;
    float dx, dy, dz;
    GET_LINE(buffer, data->file);
    numread = sscanf(buffer, "%d %*s %*f %f %f %f", &i, &dx, &dy, &dz);
    if (numread==4) {
      ts->gradient[3*(i-1)  ] = dx;
      ts->gradient[3*(i-1)+1] = dy;
      ts->gradient[3*(i-1)+2] = dz;
      numgrad++;
    }
  } while(numread==4);

  if (numgrad!=data->numatoms) {
    printf("gamessplugin) Found %d gradients for %d atoms!\n",
           numgrad, data->numatoms);
    fseek(data->file, filepos, SEEK_SET);
    return FALSE;
  }

  return TRUE;
}


/* Read energy gradients from final step
 * (different format in final step) */
static int get_final_gradient(qmdata_t *data, qm_timestep_t *ts) {
  int numgrad=0;
  int numread;
  char buffer[BUFSIZ];
  long filepos;

  /* remember position in order to go back at the end */
  filepos = ftell(data->file);

  if (pass_keyline(data->file,
                   "ATOM                 E'X", NULL) != FOUND) {
    fseek(data->file, filepos, SEEK_SET);
    return FALSE;
  }

  ts->gradient = (float *)calloc(3*data->numatoms, sizeof(float));

  if (ts->gradient == NULL) {
    PRINTERR;	    
    fseek(data->file, filepos, SEEK_SET);
    return FALSE;
  }

  /* read the gradient table */
  do {
    int i;
    float dx, dy, dz;
    GET_LINE(buffer, data->file);
    numread = sscanf(buffer, "%d %*s %f %f %f", &i, &dx, &dy, &dz);
    if (numread==4) {
      ts->gradient[3*(i-1)  ] = dx;
      ts->gradient[3*(i-1)+1] = dy;
      ts->gradient[3*(i-1)+2] = dz;
      numgrad++;
    }
  } while(numread==4);

  /* go back to where search started */
  fseek(data->file, filepos, SEEK_SET);

  if (numgrad!=data->numatoms) {
    printf("gamessplugin) Number of gradients != number of atoms!\n");
    return FALSE;
  }

  return TRUE;
}


/***********************************************************
 *
 * Read in wavenumbers and intensities of the normal modes
 *
 **********************************************************/
static int get_normal_modes(qmdata_t *data) {
  char buffer[BUFSIZ];
  int i = 0, k = 0, j = 0;
  double entry[6]; 
  char *token;

  if (!pass_keyline(data->file, "NORMAL COORDINATE ANALYSIS", NULL)) {
    return FALSE;
  }

  /* initialize array */
  memset(entry, 0, sizeof(entry));

    
  /* allocate memory for arrays */
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

  data->imag_modes = 
    (int *)calloc(data->numatoms*3,sizeof(int));
  if (data->imag_modes==NULL) {
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


  /* Example:
   *                       1           2           3           4           5
   *    FREQUENCY:      4805.54 I   1793.00 I   1317.43 I      2.13        1.99
   * REDUCED MASS:      1.09875     1.00796     1.18876     1.04129     1.03921
   * IR INTENSITY:    105.04109     2.80788     3.07217     0.01689     0.01587
   *
   * 1   CARBON       X  0.00000000  0.00000000 -0.11767262 -0.05407389  0.00000000
   *                  Y  0.00000000 -0.00345200  0.00000000  0.00000000 -0.05241437
   *                  Z -0.08676597  0.00000000  0.00000000  0.00000000  0.00000000
   */
  for (i=0; i<ceil(data->numatoms*3/5.f); i++) {
    int numread = 0;

    if (!goto_keyline(data->file, "FREQUENCY:", NULL)) {
      break;
    }
    GET_LINE(buffer, data->file);

    /* Scan the frequencies; 
     * If there are imaginary modes present then the
     * frequency is followed by the 'I' which represents
     * an additional char token in the line. */

    /* Skip first token "FREQUENCY:" */
    token = strtok(buffer, " \t\r\n");
    
    /* Walk through the remaining tokens */
    while ((token = strtok(NULL, " \t\r\n")) != NULL) {
      /* Check if token is 'I'.
       * If yes, mark previous mode as imaginary. */
      if (*token=='I') {
        data->imag_modes[data->nimag] = numread-1;
        data->nimag++;
      } else {
        /* save only the first 5 modes - there NEVER should
         * be more in any case, but just to make sure
         * we don't overrun the array */
        if (numread<5) {
          data->wavenumbers[i*5+numread] = atof(token);
          numread++;
        }
      }
    }

    eatline(data->file, 1);

    /* Read the IR INTENSITIES */
    GET_LINE(buffer, data->file);
    numread = sscanf(buffer,"%*s %*s %lf %lf %lf %lf %lf", &entry[0],
                     &entry[1], &entry[2], &entry[3], &entry[4]);
 
    for (k=0; k<numread; k++) {
      data->intensities[i*5+k] = entry[k]; 
    }

    eatline(data->file, 1);

    /* Read the normal mode vectors */
    for (k=0; k<data->numatoms; k++) {
      /* x */
      GET_LINE(buffer, data->file);
      numread = sscanf(buffer,"%*s %*s %*s %lf %lf %lf %lf %lf",
             &entry[0], &entry[1], &entry[2], &entry[3], &entry[4]);

      for (j=0; j<numread; j++) {
        data->normal_modes[3*k + (i*5+j)*3*data->numatoms] = 
          entry[j];
      }

      /* y */
      GET_LINE(buffer, data->file);
      numread = sscanf(buffer,"%*s %lf %lf %lf %lf %lf", &entry[0],
                       &entry[1],&entry[2], &entry[3],&entry[4]);

      for (j=0; j<numread; j++) {
        data->normal_modes[(3*k+1) + (i*5+j)*3*data->numatoms] =
          entry[j];
      }

      /* z */
      GET_LINE(buffer, data->file);
      numread = sscanf(buffer,"%*s %lf %lf %lf %lf %lf", &entry[0],
                       &entry[1], &entry[2], &entry[3],&entry[4]);

      for (j=0; j<numread; j++) {
        data->normal_modes[(3*k+2) + (i*5+j)*3*data->numatoms] = 
          entry[j];
      }
    }
  }


  /* Chop unused part of imag_modes array */
  data->imag_modes = 
    (int *)realloc(data->imag_modes, data->nimag*sizeof(int));

/*   free(token); */

  data->have_normal_modes = TRUE;
  printf("gamessplugin) Successfully scanned normal modes (%d imag.)\n", data->nimag);

  return TRUE;
}



/***********************************************************
 *
 * Read the cartesian hessian matrix 
 * XXX Does not read blocks with less than 6 entries correctly!
 *
 * *********************************************************/
static int get_cart_hessian(qmdata_t *data)
{
  char buffer[BUFSIZ];
  int i,j,k;
  float entry[6]; 

  buffer[0] = '\0';
  memset(entry, 0, sizeof(entry));

  /* at this point we need to rewind the file, since
   * in case that there is no internal Hessian stuff the
   * previous call to get_int_coords scanned the file
   * until EOF */
  rewind(data->file);

  if (pass_keyline(data->file,
                   "CARTESIAN FORCE CONSTANT MATRIX",
                   NULL) != FOUND) {
    return FALSE;
  }

  /* skip next 5 lines */
  eatline(data->file, 5);


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
  for (i=0; i<(int)ceil(data->numatoms/2.f); i++) {
    for (j=0; j<(data->numatoms*3)-(i*6); j++) {
      GET_LINE(buffer, data->file);
 
      if (j%3==0) {
        sscanf(buffer,"%*s %*s %*c %f %f %f %f %f %f",
               &entry[0],&entry[1],
               &entry[2],&entry[3],&entry[4],&entry[5]);
      }
      else {
        sscanf(buffer,"%*1s %f %f %f %f %f %f",
               &entry[0],&entry[1],&entry[2],&entry[3],&entry[4],
               &entry[5]);
      }


      /* save entries (lower triangular matrix) in a 
       * square matrix */
      for (k=0; k<=(j<5 ? j : 5); k++) {
        data->carthessian[(j+i*6)*3*data->numatoms + (k+i*6)] =
          entry[k];
      }
    }

    /* skip the four lines separating the data blocks */
    eatline(data->file, 4);
  }

  printf("gamessplugin) Scanned Hessian in CARTESIAN coordinates\n");

  data->have_cart_hessian = TRUE;

  return TRUE;
}
  
  
  
/***********************************************************
 *
 * Read the internal coordinates and rewind to the file
 * position where we started the search.
 *
 **********************************************************/
static int get_int_coords(qmdata_t *data) {

  char word[BUFSIZ];
  char buffer[BUFSIZ];
  long filepos, beginning;
  int first, second, third, fourth;
  double value;
  int n, i = 0, j = 0, k = 0, l = 0;

  /* remember current filepos so we can jump back */
  beginning = ftell(data->file);

  if (pass_keyline(data->file, "INTERNAL COORDINATES",
                   "1 ELECTRON INTEGRALS") != FOUND) {
    printf("gamessplugin) No internal coordinates found.\n");
    fseek(data->file, beginning, SEEK_SET);
    return FALSE;
  }

  /* skip next 5 lines */
  eatline(data->file, 5);

  /* remember current filepos so we can jump back */
  filepos = ftell(data->file);

  /* scan the next line */
  GET_LINE(buffer, data->file);
  n = sscanf(buffer,"%*s %s", word); 

  /* read line by line */
  while (n!=-1) {
    /* start counting the number of internal coordinates */
    data->nintcoords++;

    /* count the number of bonds, angles, dihedrals */
    if (!strcmp(word,"STRETCH")) {
      data->nbonds++;
    }
    else if (!strcmp(word,"BEND")) {
      data->nangles++;
    }
    else if (!strcmp(word,"TORSION")) {
      data->ndiheds++;
    }
    else if (!strcmp(word,"PLA.BEND")) {
      data->nimprops++;
    }

    /* scan next line */
    GET_LINE(buffer, data->file);
    n = sscanf(buffer,"%*s %s", word); 
  }

  /* now that we know the number of bonds, angles, etc.
   * we can read and store the internal coordinates */
  fseek(data->file, filepos, SEEK_SET);


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
    sscanf(buffer,"%*s %*s %d %d %lf", &first, &second, &value);

    *(data->bonds+2*i)   = first;
    *(data->bonds+2*i+1) = second;
    *(data->internal_coordinates+i) = value;
  }

  /* scan the BENDS */
  for (j=0; j<data->nangles; j++) {
    GET_LINE(buffer, data->file);
    sscanf(buffer,"%*s %*s %d %d %d %lf",
           &first, &second, &third, &value);

    *(data->angles+3*j)   = first;
    *(data->angles+3*j+1) = second;
    *(data->angles+3*j+2) = third;
    *(data->internal_coordinates+i+j) = value;
  }

  /* scan the TORSIONS */
  for (k=0; k<data->ndiheds; k++) {
    GET_LINE(buffer, data->file);
    sscanf(buffer,"%*s %*s %d %d %d %d %lf",
           &first, &second, &third, &fourth, &value);

    *(data->dihedrals+4*k)   = first;
    *(data->dihedrals+4*k+1) = second;
    *(data->dihedrals+4*k+2) = third;
    *(data->dihedrals+4*k+3) = fourth;
    *(data->internal_coordinates+i+j+k) = value;
  }

  /* scan the IMPROPERS */
  for (l=0; l<data->nimprops; l++) {
    GET_LINE(buffer, data->file);
    sscanf(buffer,"%*s %*s %d %d %d %d %lf",
           &first, &second, &third, &fourth, &value);

    *(data->impropers+4*l)   = first;
    *(data->impropers+4*l+1) = second;
    *(data->impropers+4*l+2) = third;
    *(data->impropers+4*l+3) = fourth;
    *(data->internal_coordinates+i+j+k+l) = value;
  }

  /* Since the internal coordinate section can appear
   * before or after the symmetry section we have to
   * jump back to the beginning of the search. */
  fseek(data->file, beginning, SEEK_SET);

  printf("gamessplugin) Scanned %d INTERNAL coordinates \n",
         data->nintcoords);
  printf("gamessplugin)    %d BONDS \n",data->nbonds);
  printf("gamessplugin)    %d ANGLES \n",data->nangles);
  printf("gamessplugin)    %d DIHEDRALS \n",data->ndiheds);
  printf("gamessplugin)    %d IMPROPERS \n",data->nimprops);

  data->have_internals = TRUE;
  return TRUE;
}



/***********************************************************
 *
 * Read the the Hessian in internal coordinates
 *
 **********************************************************/
static int get_int_hessian(qmdata_t *data) {
  char buffer[BUFSIZ];
  double hess[5];
  int i = 0, j = 0, k = 0, l = 0;

  memset(hess, 0, sizeof(hess));

  if (pass_keyline(data->file,
                   "HESSIAN MATRIX IN INTERNAL COORDINATES",
                   "ENERGY GRADIENT") != FOUND) {
    return FALSE;
  }
  if (pass_keyline(data->file,
                     "UNITS ARE HARTREE/",
                     "ENERGY GRADIENT") != FOUND) {
    return FALSE;
  }

  eatline(data->file, 3);
  
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
  for (i=0; i<(int)ceil(data->nintcoords/5.f); i++) {
    for (j=0; j<data->nintcoords; j++) {
      int numread = 0;

      GET_LINE(buffer, data->file);
      numread = sscanf(buffer,"%*d %lf %lf %lf %lf %lf", &hess[0], 
             &hess[1], &hess[2], &hess[3], &hess[4]);

      /* save entries */
      for (k=0; k<numread; k++) { 
        data->inthessian[j*data->nintcoords + i*5+k] = hess[k];
      }
    }

    /* skip the two lines separating the matrix entries 
     * and scan next line */
    eatline(data->file, 2);

    GET_LINE(buffer, data->file);
  }

#if 0
  /* read the remaining block with less then 5 rows
   * if present */
  remaining_blocks = data->nintcoords%5;
  
  if (remaining_blocks!=0) {
    for (j=0; j<data->nintcoords; j++) {
      GET_LINE(buffer, data->file);
      sscanf(buffer,"%*d %lf %lf %lf %lf %lf", &hess[0], 
             &hess[1], &hess[2], &hess[3], &hess[4]);

      for (k=0; k<remaining_blocks; k++) { 
        *(data->inthessian+(j*data->nintcoords)+(i*5)+k) = hess[k];
      }
    }
  }
#endif

  printf("gamessplugin) Scanned Hessian in INTERNAL coordinates\n");

  /* finally, dump the diagonal elements of the hessian into the
   * force constant arrays, after converting the units 
   * appropriately;
   * BONDS are in HARTREE/BOHR**2
   * ANGLES,DIHEDRALS,IMPROPERS are in HARTREE/RADIAN**2 */
  
  /* allocate dynamic arrays */
  data->bond_force_const = 
    (double *)calloc(data->nbonds, sizeof(double));

  if (data->bond_force_const==NULL) {
    PRINTERR;
    return FALSE;
  }


  data->angle_force_const =
    (double *)calloc(data->nangles, sizeof(double));

  if (data->angle_force_const==NULL) {
    PRINTERR;
    return FALSE;
  }


  data->dihedral_force_const =
    (double *)calloc(data->ndiheds, sizeof(double));

  if (data->dihedral_force_const==NULL) {
    PRINTERR;
    return FALSE;
  }


  data->improper_force_const =
    (double *)calloc(data->nimprops, sizeof(double));

  if (data->improper_force_const==NULL) {
    PRINTERR;
    return FALSE;
  }

  /* scan the bonds */
  for (i=0; i<data->nbonds; i++) {
    data->bond_force_const[i] = 
      data->inthessian[(i*data->nintcoords)+i] * 
      HARTREE_TO_KCAL / BOHR_TO_ANGS / BOHR_TO_ANGS;

    printf("%3d (BOND) %2d - %2d : %f\n", i, 
           data->bonds[2*i], data->bonds[2*i+1],
           data->bond_force_const[i]);
  }
  
  /* scan the angles */
  for (j=i; j<i+(data->nangles); j++) {
    data->angle_force_const[j-i] = 
      data->inthessian[j*data->nintcoords + j] * HARTREE_TO_KCAL;
    
    printf("%3d (ANGLE) %2d - %2d - %2d : %f\n", j,
           data->angles[3*(j-i)], data->angles[3*(j-i)+1], 
           data->angles[3*(j-i)+2], 
           data->angle_force_const[j-i]);
  }

  /* scan the dihedrals */
  for (k=j; k<j+(data->ndiheds); k++) {
    data->dihedral_force_const[k-j] = 
      data->inthessian[k*data->nintcoords + k] * HARTREE_TO_KCAL;
    
    printf("%3d (DIHEDRAL) %2d - %2d - %2d - %2d : %f \n", k,
           data->dihedrals[4*(k-j)  ], data->dihedrals[4*(k-j)+1],
           data->dihedrals[4*(k-j)+2], data->dihedrals[4*(k-j)+3],
           data->dihedral_force_const[k-j]);
  }

  /* scan the impropers */
  for (l=k; l<k+(data->nimprops); l++) {
    data->improper_force_const[l-k] = 
      data->inthessian[l*data->nintcoords + l] * HARTREE_TO_KCAL;
    
    printf("%3d (IMPROPERS) %2d - %2d - %2d - %2d : %f \n", l,
           data->impropers[4*(l-k)  ], data->impropers[4*(l-k)+1],
           data->impropers[4*(l-k)+2], data->impropers[4*(l-k)+3],
           data->improper_force_const[l-k]);
  }

  data->have_int_hessian = TRUE;
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
static int animate_normal_mode(qmdata_t *data, int mode)
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
	  (data->atoms+i)->x * (1+( k*scale * 
	  (*(normal_modes+(mode*natoms*3)+(3*i)))));

      *(animated_mode->mode_frames+(k*natoms*3)+(3*i+1)) = 
	  (data->atoms+i)->y * (1+( k*scale * 
	  (*(normal_modes+(mode*natoms*3)+(3*i+1)))));

      *(animated_mode->mode_frames+(k*natoms*3)+(3*i+2)) = 
	  (data->atoms+i)->z * (1+( k*scale * 
	  (*(normal_modes+(mode*natoms*3)+(3*i+2)))));
    }
  }


  /* second sweep all the way back to min of interval */
  for ( l = 0; l < 2*num_frames+1; ++l)
  {
    for ( i = 0; i < natoms; ++i)
    {
      *(animated_mode->mode_frames+((l+k)*natoms*3)+(3*i)) = 
	  (data->atoms+i)->x * (1+((int)(num_frames-l)*scale * 
	  (*(normal_modes+(mode*natoms*3)+(3*i)))));

      *(animated_mode->mode_frames+((l+k)*natoms*3)+(3*i+1)) = 
	  (data->atoms+i)->y * (1+((int)(num_frames-l)*scale * 
	  (*(normal_modes+(mode*natoms*3)+(3*i+1)))));

      *(animated_mode->mode_frames+((l+k)*natoms*3)+(3*i+2)) = 
	  (data->atoms+i)->z * (1+((int)(num_frames-l)*scale * 
	  (*(normal_modes+(mode*natoms*3)+(3*i+2)))));
    }
  }


  /* third sweep back to the native starting structure */
  for ( m = 0; m < num_frames+1; ++m)
  {
    for ( i = 0; i < natoms; ++i)
    {
      *(animated_mode->mode_frames+((l+k+m)*natoms*3)+(3*i)) = 
	  (data->atoms+i)->x * (1+((int)(m-num_frames)*scale * 
	  (*(normal_modes+(mode*natoms*3)+(3*i)))));

      *(animated_mode->mode_frames+((l+k+m)*natoms*3)+(3*i+1)) = 
	  (data->atoms+i)->y * (1+((int)(m-num_frames)*scale * 
	  (*(normal_modes+(mode*natoms*3)+(3*i+1)))));

      *(animated_mode->mode_frames+((l+k+m)*natoms*3)+(3*i+2)) = 
	  (data->atoms+i)->z * (1+((int)(m-num_frames)*scale * 
	  (*(normal_modes+(mode*natoms*3)+(3*i+2)))));
    }
  }

  printf("gamessplugin) Successfully animated mode %d \n", mode);

  return TRUE;
}
#endif



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
  plugin.author = "Jan Saam, Markus Dittrich, Johan Strumpfer";
  plugin.majorv = 1;
  plugin.minorv = 0;
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
