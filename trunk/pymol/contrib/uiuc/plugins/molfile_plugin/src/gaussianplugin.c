/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_gaussianplugin
#define STATIC_PLUGIN 1

/***************************************************************************
 *cr
 *cr            (C) Copyright 1995-2009 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/


/* *******************************************************
 *
 *          Gaussian Plugin
 *
 * This plugin allows VMD to read Gaussian log files.
 * It is modeled after the corresponding GAMESS plugin.
 *
 * ********************************************************/
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <time.h>
#include <math.h>

#include "gaussianplugin.h"
#include "periodic_table.h"
#include "unit_conversion.h"

#define THISPLUGIN plugin
#include "vmdconio.h"
 
/*
 * Error reporting macro for use in DEBUG mode
 */
#define GAUSSIAN_DEBUG 0
#define GAUSSIAN_BASIS_DEBUG 1

#if GAUSSIAN_DEBUG
#define PRINTERR vmdcon_printf(VMDCON_ERROR,                            \
                               "\n In file %s, line %d: \n %s \n \n",   \
                               __FILE__, __LINE__, strerror(errno))
#else
#define PRINTERR (void)(0)
#endif

/*
 * Error reporting macro for the multiple fgets calls in
 * the code
 */
#if GAUSSIAN_DEBUG
#define GET_LINE(x,y) if (!fgets(x, sizeof(x), y)) return FALSE ;   \
    else fprintf(stderr,"%s:%d %s",__FILE__, __LINE__, x)
#else
#define GET_LINE(x,y) if (!fgets(x, sizeof(x), y)) return FALSE
#endif

/* make sure pointers are NULLed after free(3)ing them. */
#define SAFE_FREE(ptr) free(ptr); ptr=NULL
/* calloc with test of success */
#define SAFE_CALLOC(ptr,type,count)                                 \
  ptr = (type *)calloc(count,sizeof(type));                         \
  if (ptr == NULL) {                                                \
    PRINTERR;                                                       \
    return MOLFILE_ERROR;                                           \
  }                                                                 \
  

#define UNK_SHELL -666
#define SPD_D_SHELL -5
#define SPD_P_SHELL -4
#define SPD_S_SHELL -3
#define SP_S_SHELL  -2
#define SP_P_SHELL  -1
#define S_SHELL 0
#define P_SHELL 1
#define D_SHELL 2
#define F_SHELL 3
#define G_SHELL 4
#define H_SHELL 5
#define I_SHELL 6

#define STATUS_UNKNOWN       -1
#define STATUS_CONVERGED      0
#define STATUS_SCF_NOT_CONV   1
#define STATUS_TOO_MANY_STEPS 2
#define STATUS_BROKEN_OFF     3

static const char *runtypes[] = { 
  "(unknown)", "ENERGY", "OPTIMIZE", "SADPOINT", "HESSIAN", 
  "SURFACE", "DYNAMICS", "PROPERTIES", NULL};

static const char *scftypes[] = { 
  "(unknown)", "RHF", "UHF", "ROHF", "GVB", "MCSCF", "FF", NULL };


/* ######################################################## */
/* declaration/documentation of internal (static) functions */
/* ######################################################## */

/** Top level gaussian log file parser. Responsible 
 *  for static, i.e. non-trajectory information. */
static int parse_static_data(gaussiandata *, int *);

/** Check if the current run is an actual Gaussian run; 
 *  returns true/false */
static int have_gaussian(gaussiandata *);

/** Get number of processors requested and amount of memory used. */
static int get_proc_mem(gaussiandata *);

/** Get basis set statistics */
static int get_basis_options(gaussiandata *);

/** Determine the run's title text. */
static int get_runtitle(gaussiandata *);

/** Read the input atom definitions.  */
static int get_input_structure(gaussiandata *);
/** Read coordinates */
static int get_coordinates(FILE *file, qm_atom_t *atoms, int numatoms);
/** Read internal orientation coordinates */
static int get_int_coordinates(FILE *file, qm_atom_t *atoms, int numatoms);

/** Read and parse the Route section of the input. */
static int get_contrl(gaussiandata *);

/** the function get_initial_info provides the atom number,
 *  coordinates, and atom types and stores them temporarily. */ 
static int get_final_info(gaussiandata *);

/* the function get_basis we also parse the basis function section to
 * determine the number of basis functions, contraction
 * coefficients. For Pople/Huzinga style basis sets
 * this numbers are in principle fixed, and could hence
 * be provided by the the plugin itself; however, the user might
 * define his own basis/contraction coeffients and hence reading
 * them from the input file seem to be somewhat more general. */
static int get_basis(gaussiandata *);

/* this function replaces the basis set data from the log file with
 * with the equivalent data read from an internal basis set data base.
 * for simplicity we use the same format as gaussian. the basis set
 * data is expected to be in $VMDDIR/basis/<basis-set-name>.gbs. */
static int get_internal_basis(gaussiandata *);

/* convert shell symmetry type from char to int */
static int shellsymm_int(char *symm);

/* Populate the flat arrays containing the basis set data */
static int fill_basis_arrays(gaussiandata *);

static int read_first_frame(gaussiandata *);

/* this subroutine scans the output file for
 * the trajectory information */
static int get_traj_frame(gaussiandata *);

/* returns 1 if the optimization has converged */
static int find_traj_end(gaussiandata *);

/* this function parses the input file for the final
 * wavefunction and stores it in the appropriate arrays; */
static int get_wavefunction(gaussiandata *, qm_timestep_t *);

/** read in mulliken charges */
static int get_population(gaussiandata *, qm_timestep_t *);

/* turn fortran double precision 'D' exponents into c parsable 'E's */
static void fix_fortran_exp(char *string) {
  while (*string) {
    if ( (*string == 'D') || (*string == 'd')) *string='e';
    ++string;
  }
}

/* ######################################################## */
/* Functions that are needed by the molfile_plugin          */
/* interface to provide VMD with the parsed data            */
/* ######################################################## */


/***************************************************************
 *
 * Called by VMD to open the Gaussian logfile and get the number
 * of atoms.
 * We are also reading all the static (i.e. non-trajectory)
 * data here since we have to parse a bit to get the atom count
 * anyway. These data will then be provided to VMD by
 * read_gaussian_metadata() and read_gaussian_rundata().
 *
 * *************************************************************/
static void *open_gaussian_read(const char *filename, 
                  const char *filetype, int *natoms) {

  FILE *fd;
  gaussiandata *data;
  
  /* open the input file */
  fd = fopen(filename, "rb");
  if (fd == NULL) {
    PRINTERR;
    return NULL;
  }

  /* set up main data structure */
  data = (gaussiandata *) calloc(1,sizeof(gaussiandata));
  if (data == NULL) return NULL;
  
  data->runtyp = RUNTYP_UNKNOWN;
  data->scftyp = SCFTYP_UNKNOWN;
  data->file = fd;
  data->file_name = strdup(filename);

  /* check if the file is Gaussian format; 
   * if yes parse it, if not exit */
  if (have_gaussian(data)==TRUE) {
    /* if we're dealing with an unsupported Gaussian
     * version, we better quit. so far we tested g98(g94?),g03 */
    if ((data->version < 19940000) || (data->version > 20040000)) {
      vmdcon_printf(VMDCON_ERROR,
                    "gaussianplugin) Gaussian version %s is not "
                    "(yet) supported. Bailing out.\n",
                    data->version_string);
      free(data);
      return NULL;
    }

    /* get the non-trajectory information from the log file */    
    if (parse_static_data(data, natoms) == FALSE) {
      free(data);
      return NULL;
    }
  }
  else {
    free(data);
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
static int read_gaussian_structure(void *mydata, int *optflags, 
                      molfile_atom_t *atoms) 
{
  gaussiandata *data = (gaussiandata *)mydata;
  qm_atom_t *cur_atom;
  molfile_atom_t *atom;
  int i = 0;
 
  /* optional data from PTE */
  *optflags = MOLFILE_ATOMICNUMBER | MOLFILE_MASS | MOLFILE_RADIUS;

  if (data->have_mulliken) 
    *optflags |= MOLFILE_CHARGE;

  /* all the information I need has already been read in
   * via the initial scan and I simply need to copy 
   * everything from the temporary arrays into the 
   * proper VMD arrays. */

  /* get initial pointer for atom array */
  cur_atom = data->initatoms;

  for(i=0; i<data->numatoms; i++) {
    atom = atoms+i;
    strncpy(atom->name, cur_atom->type, sizeof(atom->name)); 
    strncpy(atom->type, get_pte_label(cur_atom->atomicnum), sizeof(atom->type));
    strncpy(atom->resname,"QM", sizeof(atom->resname)); 
    atom->resid = 1;
    atom->chain[0] = '\0';
    atom->segid[0] = '\0';
    atom->atomicnumber = cur_atom->atomicnum;
    atom->radius = get_pte_vdw_radius(cur_atom->atomicnum);
    /* XXX; check on isotopes. should be possible to read. */
    atom->mass   = get_pte_mass(cur_atom->atomicnum);  
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
static int read_gaussian_metadata(void *mydata, 
    molfile_qm_metadata_t *gaussian_metadata) {

  gaussiandata *data = (gaussiandata *)mydata;

  gaussian_metadata->ncart = 0;
  gaussian_metadata->nimag = 0;
  gaussian_metadata->nintcoords = 0;

  /* orbital data */
  gaussian_metadata->num_basis_funcs = data->num_basis_funcs;
  gaussian_metadata->num_shells      = data->num_shells;
  gaussian_metadata->wavef_size      = data->wavef_size;  

  /* trajectory information */
  /* XXX needed here or better in read_qm_timestep_metadata()? */
  gaussian_metadata->num_traj_points = data->num_frames;

#if vmdplugin_ABIVERSION > 11
  /* charges */
  gaussian_metadata->have_esp = 0;
  gaussian_metadata->have_carthessian = 0;
  gaussian_metadata->have_internals   = 0;
  gaussian_metadata->have_normalmodes = FALSE;
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
static int read_gaussian_rundata(void *mydata, 
                               molfile_qm_t *qm_data) {

  gaussiandata *data = (gaussiandata *)mydata;
  int i;

  molfile_qm_basis_t   *basis_data   = &qm_data->basis;
  molfile_qm_sysinfo_t *sys_data     = &qm_data->run;

  /* fill in molfile_qm_sysinfo_t */
  sys_data->nproc = data->nproc;
  sys_data->memory = data->memory; 
  sys_data->runtyp = data->runtyp;
  sys_data->scftyp = data->scftyp;
  sys_data->totalcharge = data->totalcharge;
  sys_data->multiplicity = data->multiplicity;
/*   sys_data->num_wave_f = data->num_wave_f; */
  sys_data->num_electrons = data->num_electrons;
  sys_data->num_orbitals_A = data->num_orbitals_A;
  sys_data->num_orbitals_B = data->num_orbitals_B;

  strncpy(sys_data->basis_string, data->basis_string,
          sizeof(sys_data->basis_string));
  
  strncpy(sys_data->runtitle, data->runtitle, sizeof(sys_data->runtitle));
  strncpy(sys_data->geometry, data->geometry, sizeof(sys_data->geometry));
  strncpy(sys_data->version_string, data->version_string,
          sizeof(sys_data->version_string));

#if vmdplugin_ABIVERSION > 11
  /* fill in molfile_qm_basis_t */
  if (data->num_basis_funcs) {
    for (i=0; i<data->numatoms; i++) {
      basis_data->num_shells_per_atom[i] = data->num_shells_per_atom[i];
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
  gaussiandata *data = (gaussiandata *)mydata;
#if GAUSSIAN_DEBUG
  vmdcon_printf(VMDCON_INFO, 
                "gaussianplugin) read_qm_timestep_metadata(): %d/%d/%d\n",
                data->num_frames, 
                data->num_frames_read,
                data->num_frames_sent);
#endif

  meta->count = -1; /* Don't know the number of frames yet */
  meta->has_gradient = 0;

  if (data->num_frames_read > data->num_frames_sent) {
    have = 1;
  } else if (data->num_frames_read < data->num_frames) {
#if GAUSSIAN_DEBUG
    vmdcon_printf(VMDCON_INFO,
                  "gaussianplugin) Probing timestep %d\n",
                  data->num_frames_read);
#endif
    have = get_traj_frame(data);
  }

  if (have) {
    /* get a pointer to the current qm timestep */
    qm_timestep_t *cur_qm_ts = data->qm_timestep+data->num_frames_sent;
#if GAUSSIAN_DEBUG
    vmdcon_printf(VMDCON_INFO,
                  "gaussianplugin) Approved timestep %d\n", 
                  data->num_frames_sent);
#endif
    meta->num_scfiter  = cur_qm_ts->num_scfiter;
    meta->num_orbitals_per_wavef[0] = cur_qm_ts->orbital_counter;
    meta->wavef_size = data->wavef_size;

  } else {
    meta->num_scfiter  = 0;
    meta->num_orbitals_per_wavef[0] = 0;
    meta->wavef_size = 0;

    data->end_of_trajectory = TRUE;
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
  gaussiandata *data = (gaussiandata *)mydata;
  qm_atom_t *cur_atom;
  int i = 0;
  qm_timestep_t *cur_qm_ts;

  if (data->end_of_trajectory == TRUE) return MOLFILE_ERROR;

#if GAUSSIAN_DEBUG
  vmdcon_printf(VMDCON_INFO,
                "gaussianplugin) Sending timestep %d\n", 
                data->num_frames_sent);
#endif

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
  cur_qm_ts = data->qm_timestep+data->num_frames_sent;

  /* store the SCF energies */
  for (i=0; i<cur_qm_ts->num_scfiter; i++) {
    qm_ts->scfenergies[i] = cur_qm_ts->scfenergies[i];
  }

  /* store the wave function and orbital energies */
  if (cur_qm_ts->wave_function && cur_qm_ts->orbital_counter) {
    memcpy(qm_ts->wave_function, cur_qm_ts->wave_function,
	    cur_qm_ts->orbital_counter*data->wavef_size*sizeof(float));
    
    memcpy(qm_ts->orbital_energies, cur_qm_ts->orbital_energies,
	    cur_qm_ts->orbital_counter*sizeof(float));
  }

  if (data->runtyp == RUNTYP_ENERGY || data->runtyp == RUNTYP_HESSIAN) {
    /* We have only a single point */
    data->end_of_trajectory = TRUE;
  }

  data->num_frames_sent++;

  return MOLFILE_SUCCESS;
}
#endif



/** Clean up when done and free all memory 
 *  to avoid memory leaks.
 *
 **********************************************************/
static void close_gaussian_read(void *mydata) {

  gaussiandata *data = (gaussiandata *)mydata;
  int i, j;
  fclose(data->file);

  free(data->file_name);
  free(data->initatoms);
  free(data->basis);
  free(data->shell_symmetry);
  free(data->num_shells_per_atom);
  free(data->num_prim_per_shell);
  free(data->mulliken_charges);
  free(data->internal_coordinates);
  free(data->wavenumbers);
  free(data->intensities);
  free(data->normal_modes);
  free(data->angular_momentum);

  if (data->basis_set) {
    for(i=0; i<data->numatoms; i++) {
      for (j=0; j<data->basis_set[i].numshells; j++) {
        free(data->basis_set[i].shell[j].prim);
      }
      free(data->basis_set[i].shell);
    } 
    free(data->basis_set);
  }

  for (i=0; i<data->num_frames_read; i++) {
    free(data->qm_timestep[i].scfenergies);
    free(data->qm_timestep[i].orbital_energies);
    free(data->qm_timestep[i].wave_function);
    free(data->qm_timestep[i].gradient);
    free(data->qm_timestep[i].mulliken_charges);
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
 * Main gaussian log file parser responsible for static,  
 * i.e. non-trajectory information.
 *
 ********************************************************/
static int parse_static_data(gaussiandata *data, int *natoms) 
{
  /* Read # of procs and amount of requested memory */
  if (!get_proc_mem(data))        return FALSE;

  /* Read the basis options */
  if (!get_basis_options(data))   return FALSE;

  /* Read the route section and try to determine
   * the job type. exit if unsupported. */
  if (!get_contrl(data))          return FALSE;

  /* Read the run title */
  if (!get_runtitle(data))        return FALSE;

  /* Read the input atom definitions and geometry */
  if (!get_input_structure(data)) return FALSE;
  /* provide VMD with the proper number of atoms */
  *natoms = data->numatoms;
  /* read first set of coordinates and basis set data */
  read_first_frame(data);
  
  /* */
  get_final_info(data);

  vmdcon_printf(VMDCON_INFO, "gaussianplugin) found %d QM data frames.\n", data->num_frames);
#if GAUSSIAN_DEBUG
  vmdcon_printf(VMDCON_INFO, "gaussianplugin) num_frames_read = %d\n", data->num_frames_read);
  vmdcon_printf(VMDCON_INFO, "gaussianplugin) num_frames_sent = %d\n", data->num_frames_sent);
#endif
  return TRUE;
}

/**********************************************************
 *
 * this subroutine checks if the provided files is
 * actually a Gaussian file;
 *
 **********************************************************/
static int have_gaussian(gaussiandata *data) 
{
  char word[4][MOLFILE_BUFSIZ];
  char buffer[BUFSIZ];
  char *ptr;
  int i = 0;
 
  buffer[0] = '\0';
  for (i=0; i<3; i++) word[i][0] = '\0';

  /* check if the file is Gaussian format 
   * Gaussian output typically begins with:
   * 'Entering Gaussian System' */
  i=0; /* check only the first 100 lines */
  do {
    GET_LINE(buffer, data->file);
    sscanf(buffer,"%s%s%s",word[0],word[1],word[2]);
    ++i;
  } while( (strcmp(word[0],"Entering") || 
            strcmp(word[1],"Gaussian") || 
            strcmp(word[2],"System,")) && (i<100) );
  if (i>=100) return FALSE;
  vmdcon_printf(VMDCON_INFO, "gaussianplugin) Analyzing Gaussian log file: %s\n",data->file_name);
  
  /* now read on until we find the block of text with encoded version
   * number and compile date. */
  i=0; /* check only the next 100 lines */
  do {
    GET_LINE(buffer, data->file);
    buffer[20] = '\0';          /* length of the version block varies. */
  } while ( (i<100) && strcmp(buffer,
      " *******************"));
  if (i>=100) return FALSE;

  /* one more line... */
  GET_LINE(buffer, data->file);
  sscanf(buffer,"%s%s%s%s",word[0],word[1],
         word[2],word[3]);

  data->version = atoi(word[1]) * 10000;
  if (data->version > 700000) {
      data->version += 19000000;
  } else {
      data->version += 20000000;
  }
  strcpy(data->version_string,word[2]);

  ptr=strrchr(word[2],'-');
  if (ptr != NULL) { 
      /* extract revision and patchlevel from G##Rev%.## word.*/
      ptr += 7;
      i =  (*ptr) - 'A' + 1;
      data->version += i*100;
      ++ptr; ++ptr;
      data->version += atoi(ptr);
  }
  
  vmdcon_printf(VMDCON_INFO, 
                "gaussianplugin) Gaussian version = %s  (Version code: %d)\n",
                data->version_string, data->version);
  vmdcon_printf(VMDCON_INFO,
                "gaussianplugin) Compiled on      = %s \n", word[3]);


  return TRUE;
}


/**********************************************************
 *
 * this subroutine reads the number of procs and the amount
 * of memory requested
 *
 **********************************************************/
static int get_proc_mem(gaussiandata *data) {

  char word[5][MOLFILE_BUFSIZ];
  char buffer[BUFSIZ];
  int nproc,maxmem;
  int i, link;

  buffer[0] = '\0';
  for (i=0; i<3; i++) word[i][0] = '\0';

  rewind(data->file);

  /* set some defaults */
  nproc = 1;
  maxmem = -1;
  
  do {
      GET_LINE(buffer, data->file);
      sscanf(buffer,"%s%s%s%*s%s%*s%*s%*s%*s%*s%s",
             word[0],word[1],word[2],word[3],word[4]);

      /* max dynamical memory. */ 
      if ((strcmp(word[0],"Leave") == 0) &&
          (strcmp(word[1],"Link")  == 0)) {
        link = atoi(word[2]);
        /* gaussian uses real*8 words internally. convert to MByte. */
        if (link > 1) maxmem=atoi(word[4])/128/1024;
      }

      /* number of SMP cpus */
      if ( (strcmp(word[0],"Will") == 0) &&
           (strcmp(word[1],"use")  == 0) &&
           (strcmp(word[2],"up")  == 0) ) {
        nproc = atoi(word[3]);
      }

      /* detect if have read too far */ 
      if (((strcmp(word[0],"Standard") == 0) ||
           (strcmp(word[0],"Z-Matrix") == 0) ||
           (strcmp(word[0],"Input") == 0) ) &&
          (strcmp(word[1],"orientation:")  == 0)) {
        /* can not detect memory used */
        maxmem=0;
      }
  } while (maxmem < 0);

  /* store findings */
  data->nproc = nproc;
  data->memory = maxmem;
  if (maxmem) 
    vmdcon_printf(VMDCON_INFO, 
                  "gaussianplugin) Gaussian used %2d SMP process(es), "
                  "% 6d Mbytes of memory \n", nproc, maxmem);

  return TRUE;
}


/**********************************************************
 *
 * Extract basis set options
 *
 **********************************************************/
static int get_basis_options(gaussiandata *data) {

  char word[5][MOLFILE_BUFSIZ];
  char buffer[BUFSIZ], *ptr;
  int i = 0, nfunc, nprim, nume_a, nume_b;

  buffer[0] = '\0';
  for (i=0; i<3; i++) word[i][0] = '\0';

  /* to be safe let's rewind the file */
  rewind(data->file);

  /* scanning for basis set string */
  nume_a=-1;
  do {
    GET_LINE(buffer, data->file);
    i=sscanf(buffer,"%s%s%s",word[0],word[1],word[2]);
    if (i==3) {

      if ( (strcmp(word[0],"Standard") == 0) &&
           (strcmp(word[1],"basis:") == 0) ) {

        ptr = &buffer[0] + strlen(buffer) - 1;
        while(*ptr==' ') --ptr;
        *ptr='\0';
        strncpy(data->gbasis, word[2], 10);
        strncpy(data->basis_string, buffer+17, MOLFILE_BUFSIZ);

        /* make sure gbasis is uppercase. */
        ptr=data->gbasis;
        for (;*ptr;++ptr) *ptr=toupper(*ptr);

      } else if ( (strcmp(word[0],"General") == 0) &&
           (strcmp(word[1],"basis") == 0) ) {

        /* General basis read from cards */

        ptr = &buffer[0] + strlen(buffer) - 1;
        while(*ptr==' ') --ptr;
        *ptr='\0';

        strncpy(data->gbasis, "GEN", 4);
        strncpy(data->basis_string, buffer, MOLFILE_BUFSIZ);
      
      } else {

        i=sscanf(buffer,"%d%s%s%d%s",&nfunc,word[0],word[1],
                 &nprim,word[2]);
        if (i==5) {
          if ((strncmp(word[0],"basis",5) == 0)     &&
              (strncmp(word[1],"functions",9) == 0) &&
              (strncmp(word[2],"primitive",9) == 0) ) {
            GET_LINE(buffer, data->file);
            sscanf(buffer,"%d%s%s%d%s",&nume_a,word[0],word[1],
                   &nume_b,word[2]);
          }
        }
      }
    }
  } while (nume_a < 0);

  vmdcon_printf(VMDCON_INFO, 
                "gaussianplugin) with %d basis functions, %d primitive gaussians\n",
                nfunc, nprim);
  vmdcon_printf(VMDCON_INFO, 
                "gaussianplugin) with %d alpha electrons, %d beta electron\n",
                nume_a, nume_b);

  data->wavef_size     = nfunc;
  data->num_orbitals_A = nume_a;
  data->num_orbitals_B = nume_b;
  data->num_electrons  = nume_a + nume_b;
  data->multiplicity   = (nume_a+nume_b)/2-(nume_a>nume_b?nume_b:nume_a)+1;
  
  return TRUE;
}


/**********************************************************
 *
 * Extract the run title line
 *
 **********************************************************/
static int get_runtitle(gaussiandata *data) {

  char buffer[BUFSIZ];
  char *temp;
  size_t offs;
  
  /* look for RUN TITLE section */
  do {
    GET_LINE(buffer, data->file);
  } while (strncmp(buffer," -", 2) );
  
  GET_LINE(buffer, data->file);
  do {
    char s; 

    offs=strlen(buffer);
    temp=&buffer[offs];

    while (*temp == '\n' || *temp == '\r' || *temp == '\0') --temp;
    s = *temp;
    fgets(temp,BUFSIZ-offs,data->file);
    *temp = s;
  } while ( strncmp(temp+1,"--",2) );

  offs=strlen(buffer);
  temp=&buffer[offs-1];
  while (*temp == '-' || *temp == '\n' || *temp == '\r') --temp;
  ++temp;  *temp='\0';
  strncpy(data->runtitle,buffer,sizeof(data->runtitle));
  
  return TRUE;
} 


/* Read the input atom definitions and geometry */
static int get_input_structure(gaussiandata *data) {
  char buffer[BUFSIZ];
  char word[4][MOLFILE_BUFSIZ];
  int i, numatoms;

  buffer[0] = '\0';
  for (i=0; i<4; i++) word[i][0] = '\0';
  
  GET_LINE(buffer, data->file);
  sscanf(buffer,"%s%s",word[0],word[1]);
  
  if ( ( (strcmp(word[0],"Symbolic") == 0) &&
         (strcmp(word[1],"Z-matrix:") == 0) ) ||
       ( (strcmp(word[0],"Redundant") == 0) &&
         (strcmp(word[1],"internal") == 0)  ) || 
       ( (strcmp(word[0],"Z-Matrix") == 0) &&
         (strcmp(word[1],"taken") == 0) ) ) {

    /* skip over line with checkpoint file name */
    if ( ( (strcmp(word[0],"Redundant") == 0) &&
           (strcmp(word[1],"internal") == 0)  ) || 
         ( (strcmp(word[0],"Z-Matrix") == 0) &&
           (strcmp(word[1],"taken") == 0) ) ) {
      GET_LINE(buffer, data->file);
    }
      
    /* charge and multiplicity */
    GET_LINE(buffer, data->file);
    sscanf(buffer,"%*s%*s%d%*s%*s%d", &(data->totalcharge),
           &(data->multiplicity));

    /* parse coordinates from input */
    numatoms=0;
    data->initatoms=NULL;
    i=1;
    do {
      char *ptr;
      qm_atom_t *atm;
      
      i=1;
      GET_LINE(buffer, data->file);
      ptr = strtok(buffer," ,\t\n");
      
      if (ptr == NULL) break;
      if (strcmp(ptr, "Variables:") == 0) break;
      if (strcmp(ptr, "Recover") == 0) break;
      if (strcmp(ptr, "The") == 0) break;
      if (strcmp(ptr, "Leave") == 0) break;
      
      /* for now we only read the atom label */
      data->initatoms=realloc(data->initatoms,(numatoms+1)*sizeof(qm_atom_t));
      atm = data->initatoms + numatoms;
      strncpy(atm->type, ptr, sizeof(atm->type));
      atm->atomicnum=get_pte_idx(ptr);
      ++numatoms;
    } while ( i >= 0 );
    /* TODO */
  } else {
    vmdcon_printf(VMDCON_ERROR,
                  "gaussianplugin) ERROR, cannot parse input coordinates.\n");
    return FALSE;
  }
  
  /* store number of atoms in data structure */
  data->numatoms = numatoms;
  vmdcon_printf(VMDCON_INFO, 
                "gaussianplugin) Atoms: %d   Charge: %d   Multiplicity: %d\n", 
                numatoms, data->totalcharge, data->multiplicity);
  return TRUE; 
}


/**********************************************************
 *
 * Read data from the Route section.
 *
 * XXX: there currently do not seem to be any provisions
 * XXX: for multi-step jobs, i.e. jobs that either have
 * XXX: multiple run types in one calculations or perform 
 * XXX: several calculation steps subsequently...
 *
 **********************************************************/
static int get_contrl(gaussiandata *data) {

  char buffer[BUFSIZ];
  const char *vmdbasis;
  char *temp;
  size_t offs;

  buffer[0] = '\0';
  rewind(data->file);
  
  /* try to find route section. */
  do {
    GET_LINE(buffer, data->file);
  } while( strncmp(buffer," #", 2) );

  /* append to string in buffer until we reach the end of 
   * the route section. this way we'll have the whole info
   * in one long string.
   * Gaussian writes this out with format (x,A70), so 
   * joining lines need some special tricks.  */
  do {
    char s; 

    offs=strlen(buffer);
    temp=&buffer[offs];

    while (*temp == '\n' || *temp == '\r' || *temp == '\0') --temp;
    s = *temp;
    fgets(temp,BUFSIZ-offs,data->file);
    *temp = s;
  } while ( strncmp(temp+1,"--",2) );

  offs=strlen(buffer);
  temp=&buffer[offs-1];
  while (*temp == '-' || *temp == '\n' || *temp == '\r') --temp;

  /* some more magic is required.
   * make sure we are zero terminated and have a trailing blank.
   * the latter allows to search for keywords with strstr without
   * getting false positives on keyword strings that are contained
   * in other keywords */
  ++temp;  *temp=' '; ++temp; *temp='\0';

  /* convert to upper case */
  temp=&buffer[0];
  while (*temp++) *temp = toupper(*temp);

  /* now we are ready to look for useful information */

  /* will we have wavefunction output? */
  if (strstr(buffer," IOP(6/7=3) ")) {
    data->have_wavefunction=TRUE;
  } else {
    data->have_wavefunction=FALSE;
  }

  /* is there a basis set in "input format" */
  if (strstr(buffer," GFINPUT ")) {
    data->have_basis=TRUE;
  } else {
    data->have_basis=FALSE;
  }

  /* cartesian d-basis functions */
  if (strstr(buffer," 6D ")) {
    data->have_cart_basis |= 1;
  }
  /* cartesian f-basis functions */
  if (strstr(buffer," 10F ")) {
    data->have_cart_basis |= 2;
  }
  /* pure d-basis functions */
  if (strstr(buffer," 5D ")) {
    data->have_cart_basis &= ~1;
  }
  /* pure f-basis functions */
  if (strstr(buffer," 7F ")) {
    data->have_cart_basis &= ~2;
  }

  /* find scf calculation type. XXX: might be safer to detect from output. */
  if ((strstr(buffer," ROHF/")) ||
      (strstr(buffer," ROHF ")) ||
      (strstr(buffer," ROMP"))) {
    data->scftyp = SCFTYP_ROHF;
  } else if (data->multiplicity != 1) {
    data->scftyp = SCFTYP_UHF;
  } else {
    data->scftyp = SCFTYP_RHF;
  }
        
  /* for semi-empirical, we set the basis to valence-STO-3G. */
  if ((strstr(buffer," AM1/")) ||
      (strstr(buffer," AM1 ")) ||
      (strstr(buffer," PM3/")) ||
      (strstr(buffer," PM3 ")) ||
      (strstr(buffer," MNDO/")) ||
      (strstr(buffer," MNDO "))) {

    vmdbasis = getenv("VMDDEFBASISSET");
    if (vmdbasis == NULL) 
      vmdbasis = "VSTO-3G";
    
    if(strlen(data->gbasis) == 0)
      strncpy(data->gbasis, vmdbasis, sizeof(data->gbasis));

    if(strlen(data->basis_string) == 0) {
      strncpy(data->basis_string, "Internal ", sizeof(data->basis_string));
      strncat(data->basis_string, vmdbasis, sizeof(data->basis_string) - 10);

      if(data->have_cart_basis & 1)
        strcat(data->basis_string," 6D");
      else
        strcat(data->basis_string," 5D");
      if(data->have_cart_basis & 2)
        strcat(data->basis_string," 10F");
      else
        strcat(data->basis_string," 7F");
    }
  }

  /* find run type. XXX: might be safer to detect from output. */
  data->runtyp = RUNTYP_ENERGY;  /* default */
  if ((strstr(buffer," FOPT ")) || 
      (strstr(buffer," FOPT=")) ||
      (strstr(buffer," FOPT(")) ||
      (strstr(buffer," OPT=")) ||
      (strstr(buffer," OPT(")) ||
      (strstr(buffer," OPT "))) {
    data->runtyp = RUNTYP_OPTIMIZE;
  } 
  if (strstr(buffer," FREQ ")) {
    data->runtyp = RUNTYP_HESSIAN;
  }
  if (strstr(buffer," SCAN ")) {
    data->runtyp = RUNTYP_SURFACE;
  }
        
  /* make sure we have some basis set string
   * and allow to set it from within VMD. */
  vmdbasis = getenv("VMDDEFBASISSET");
  if(strlen(data->gbasis) == 0) {
    if (vmdbasis == NULL) {
      strncpy(data->gbasis, "(unknown)", sizeof(data->gbasis));
      strncpy(data->basis_string, "(unknown)", sizeof(data->basis_string));
    } else {
      strncpy(data->gbasis, vmdbasis, sizeof(data->gbasis));
      strncpy(data->basis_string, "Internal ", sizeof(data->basis_string));
      strncat(data->basis_string, vmdbasis, sizeof(data->basis_string) - 10);

      if(data->have_cart_basis & 1)
        strcat(data->basis_string," 6D");
      else
        strcat(data->basis_string," 5D");
      if(data->have_cart_basis & 2)
        strcat(data->basis_string," 10F");
      else
        strcat(data->basis_string," 7F");
    }
  }

  /* TODO: add more. e.g. Opt=(ModRedundant), IRC */

  /* print (some of) our findings */
  vmdcon_printf(VMDCON_INFO, "gaussianplugin) Run-type: %s, SCF-type: %s\n",
                runtypes[data->runtyp], scftypes[data->scftyp]);
  vmdcon_printf(VMDCON_INFO, 
                "gaussianplugin) using %s basis set.\n", data->basis_string);
  
  return TRUE;
}

static int read_first_frame(gaussiandata *data) {
  data->qm_timestep = NULL;

  /* the angular momentum is populated in get_wavefunction 
   * which is called by get_traj_frame(). We have obtained
   * the array size wavef_size already from the basis set
   * statistics */
  vmdcon_printf(VMDCON_INFO, 
                "gaussianplugin) allocating data for %d wavefunctions\n",
                data->wavef_size);
  SAFE_CALLOC(data->angular_momentum,int,3*data->wavef_size);

  if (!get_traj_frame(data)) {
    return FALSE;
  }

  data->num_frames = 1;
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
static int get_final_info(gaussiandata *data) {
  long filepos;
  filepos = ftell(data->file);

  if (data->runtyp == RUNTYP_OPTIMIZE || 
      data->runtyp == RUNTYP_SADPOINT ||
      data->runtyp == RUNTYP_SURFACE) {
    /* Try to advance to the end of the geometry
     * optimization. If no regular end is found we
     * won't find any propertiies to read and return. */
    if (!find_traj_end(data)) return FALSE;
  }

  if (data->runtyp == RUNTYP_HESSIAN || data->runtyp == RUNTYP_SURFACE) {
    vmdcon_printf(VMDCON_WARN, "gaussianplugin) this run type is not fully supported\n");
  }
    
  fseek(data->file, filepos, SEEK_SET);
  return TRUE; 
}


static int get_coordinates(FILE *file, qm_atom_t *atoms, int numatoms) {
  char buffer[BUFSIZ];
  int atomicnum;
  float x,y,z;
  int i,n;

  /* we look for:
                              Input orientation:                          
 ---------------------------------------------------------------------
 Center     Atomic     Atomic              Coordinates (Angstroms)
 Number     Number      Type              X           Y           Z
 ---------------------------------------------------------------------
   */
  do {
    GET_LINE(buffer, file);
  } while ((strstr(buffer,"Input orientation:") == NULL) &&
           (strstr(buffer,"Z-Matrix orientation:") == NULL));
  GET_LINE(buffer, file);
  GET_LINE(buffer, file);
  GET_LINE(buffer, file);
  GET_LINE(buffer, file);
  
  for (i=0; i < numatoms; ++i) {
    GET_LINE(buffer, file);
    n = sscanf(buffer,"%*d%d%*d%f%f%f",&atomicnum,&x,&y,&z);
    if (n!=4) return FALSE;
    atoms[i].x=x;
    atoms[i].y=y;
    atoms[i].z=z;
  }
  return TRUE;
}

/** same as get_coordinates, but we look for coordinates 
    in internal orientation */
static int get_int_coordinates(FILE *file, qm_atom_t *atoms, int numatoms) {
  char buffer[BUFSIZ];
  int atomicnum;
  float x,y,z;
  int i,n;

  /* we look for:
                              Standard orientation:                          
 ---------------------------------------------------------------------
 Center     Atomic     Atomic              Coordinates (Angstroms)
 Number     Number      Type              X           Y           Z
 ---------------------------------------------------------------------
   */
  do {
    GET_LINE(buffer, file);
  } while ((strstr(buffer,"Standard orientation:") == NULL) &&
           (strstr(buffer,"Z-Matrix orientation:") == NULL));
  GET_LINE(buffer, file);
  GET_LINE(buffer, file);
  GET_LINE(buffer, file);
  GET_LINE(buffer, file);
  
  for (i=0; i < numatoms; ++i) {
    GET_LINE(buffer, file);
    n = sscanf(buffer,"%*d%d%*d%f%f%f",&atomicnum,&x,&y,&z);
    if (n!=4) return FALSE;
    atoms[i].x=x;
    atoms[i].y=y;
    atoms[i].z=z;
  }
  return TRUE;
}



/*******************************************************
 *
 * this function reads in the basis set data 
 *
 * ******************************************************/
/* typical data looks like this: 
 ****
   <atomnum> 0
   <shellsymm> <#prims> <scalefactor> 0.0
   <#prims> lines of <coeff(i)> <exp(i>
   <shellsymm> <#prims> <scalefactor> 0.0
   <#prims> lines of <coeff(i)> <exp(i>
    ...
 ****    

  1 0
 S   6 1.00       0.000000000000
      0.1941330000D+05  0.1851598923D-02
      0.2909420000D+04  0.1420619174D-01
      0.6613640000D+03  0.6999945928D-01
      0.1857590000D+03  0.2400788603D+00
      0.5919430000D+02  0.4847617180D+00
      0.2003100000D+02  0.3351998050D+00
 SP   6 1.00       0.000000000000
      0.3394780000D+03 -0.2782170105D-02  0.4564616191D-02
      0.8101010000D+02 -0.3604990135D-01  0.3369357188D-01
      0.2587800000D+02 -0.1166310044D+00  0.1397548834D+00
      0.9452210000D+01  0.9683280364D-01  0.3393617168D+00
      0.3665660000D+01  0.6144180231D+00  0.4509206237D+00
      0.1467460000D+01  0.4037980152D+00  0.2385858009D+00
 */

int get_basis(gaussiandata *data) {

  char buffer[BUFSIZ];
  char word[3][MOLFILE_BUFSIZ];
  int i; 

  /* no point in searching through the log file,
   * if we cannot have this information */
  if (!data->have_basis) return FALSE;
  
  /* search for characteristic line. but only in the next 1000 lines */
  i=0;
  if(data->version < 20030000) {            /* g98 */
    do {
      GET_LINE(buffer, data->file);
      sscanf(buffer,"%s%s",word[0],word[1]);
      ++i;
    } while( (strcmp(word[0],"Basis") || 
              strcmp(word[1],"set"))  && (i<1000) );
  } else {                                  /* g03 */
    do {
      GET_LINE(buffer, data->file);
      sscanf(buffer,"%s%s",word[0],word[1]);
      ++i;
    } while( (strcmp(word[0],"AO") || 
              strcmp(word[1],"basis"))  && (i<1000) );
  }
  if (i>=1000) {
    /* flag that we have no basis set data */
    data->num_shells_per_atom=NULL;
    data->have_basis=FALSE;
    return MOLFILE_ERROR;
  }
  
  /* Allocate space for the basis for all atoms */
  /* When the molecule is symmetric the actual number atoms with
   * a basis set could be smaller */
  SAFE_CALLOC(data->basis_set,basis_atom_t,data->numatoms);

  for (i=0; i< data->numatoms; ++i) {
    int numshells, numprim;
    int numread, ishell;
    float scalef;
    shell_t *shell;
    
    /* this line should be '<atomindex> 0'. gaussian allows much more
     * flexible input (it is a zero terminated list of indices or labels), 
     * but upon writing to the log files there is exactly one basis set 
     * specification per atom. */
    GET_LINE(buffer, data->file);
    if (atoi(buffer) != (i+1) ) {
      vmdcon_printf(VMDCON_WARN, 
                    "gaussianplugin) basis set atom counter mismatch: %d vs %s\n",
                    atoi(buffer), i+1);
      SAFE_FREE(data->basis_set);
      data->num_shells_per_atom=NULL;
      return MOLFILE_ERROR;
    }
    strncpy(data->basis_set[i].name,data->gbasis,10);
    
    numshells=0;
    shell=NULL;
    data->num_shells_per_atom=NULL;
    GET_LINE(buffer, data->file);
    while( strncmp(buffer," ****",5) ) {
      int n;

      numread=sscanf(buffer,"%s%d%f",word[0],&numprim,&scalef);
      if (numread == 3) {
#if GAUSSIAN_DEBUG && GAUSSIAN_BASIS_DEBUG
        vmdcon_printf(VMDCON_INFO, "gaussianplugin) atom: %d, element: %s, shell: %d "
                      "%s-type shell, %d primitives, scalefactor %f\n", i,
                      get_pte_label(data->initatoms[i].atomicnum), numshells+1, 
                      word[0], numprim, scalef);
#endif
        ;
      } else {
        vmdcon_printf(VMDCON_WARN, 
                      "gaussianplugin) basis set parse error: %s",buffer);
        free(data->basis_set);
        data->have_basis=FALSE;
        data->basis_set=NULL;
        data->num_shells_per_atom=NULL;
        return FALSE;
      }
      ishell=numshells;
      ++numshells;
      shell=realloc(shell,numshells*sizeof(shell_t));
      shell[ishell].numprims=numprim;
      shell[ishell].symmetry=shellsymm_int(word[0]);
      shell[ishell].prim = (prim_t *)calloc(numprim,sizeof(prim_t));
      data->num_basis_funcs += numprim;
      if (shell[ishell].symmetry == SP_S_SHELL) {
        ++numshells;
        shell=realloc(shell,numshells*sizeof(shell_t));
        shell[ishell+1].numprims=numprim;
        shell[ishell+1].symmetry=SP_P_SHELL;
        shell[ishell+1].prim = (prim_t *)calloc(numprim,sizeof(prim_t));
        data->num_basis_funcs += numprim;
      }

      for (n=0; n<numprim; ++n) {
        GET_LINE(buffer, data->file);
        sscanf(buffer,"%s%s%s", word[0],word[1],word[2]);
        fix_fortran_exp(word[0]);
        shell[ishell].prim[n].exponent=atof(word[0])*scalef*scalef;
        fix_fortran_exp(word[1]);
        shell[ishell].prim[n].contraction_coeff=atof(word[1]);
        if (shell[ishell].symmetry == SP_S_SHELL) {
          shell[ishell+1].prim[n].exponent=shell[ishell].prim[n].exponent;
          fix_fortran_exp(word[2]);
          shell[ishell+1].prim[n].contraction_coeff=atof(word[2]);
        }
      }
      GET_LINE(buffer, data->file);
    }

    /* store shells in atom */
    data->basis_set[i].numshells = numshells;
    data->basis_set[i].shell = shell;
    
    /* store the total number of basis functions */
    data->num_shells += numshells;
  }

  vmdcon_printf(VMDCON_INFO, 
                "gaussianplugin) Parsed %d uncontracted basis functions. \n",
                data->num_basis_funcs);

  /* allocate and populate flat arrays needed for molfileplugin */
  return fill_basis_arrays(data);
}


/*******************************************************
 *
 * this function reads in the basis set data from 
 * <basis>.gbs or $VMDDIR/basis/<basis>.gbs
 *
 * ******************************************************/
int get_internal_basis(gaussiandata *data) {

  char *vmddir = NULL;
  FILE *fp;
  char buffer[BUFSIZ];
  char word[3][MOLFILE_BUFSIZ];
  char filepath[256];
  int  i,n; 

  /* no point in adding a basis set if we already have this information */
  if (data->have_basis) return TRUE;

  /* try to open basis set database file. */
  sprintf(filepath,"%s.gbs",data->gbasis);
  fp=fopen(filepath,"rb");
  if (fp == NULL) {
      vmddir=getenv("VMDDIR");
      if (vmddir == NULL) {
          vmddir="/usr/local/lib/vmd";
      }
      sprintf(filepath,"%s/basis/%s.gbs",vmddir,data->gbasis);
      fp=fopen(filepath,"rb");
  }
  
  if (fp == NULL) {
    vmdcon_printf(VMDCON_ERROR, "gaussianplugin) failed to read basis set "
                  "from data base file %s\n", filepath);
    data->num_shells_per_atom=NULL;
    data->have_basis=FALSE;
    return FALSE;
  } else {
    vmdcon_printf(VMDCON_INFO, "gaussianplugin) reading basis set "
                  "from data base file %s\n", filepath);
  }
  
  /* Allocate space for the basis for all atoms */
  /* When the molecule is symmetric the actual number atoms with
   * a basis set could be smaller */
  SAFE_CALLOC(data->basis_set,basis_atom_t,data->numatoms);

  for (i=0; i < data->numatoms; ++i) {
    int numshells, numprim;
    int numread, ishell;
    float scalef;
    shell_t *shell;
    
    /* search for the characteristic first line starting with '****'. */
    rewind(fp);
    do {
      fgets(buffer, sizeof(buffer), fp);
      sscanf(buffer,"%s%s",word[0],word[1]);
    } while(strcmp(word[0],"****"));
  
    /* search for an entry for the current atom in the format '<name> 0'. */
    do {
      fgets(buffer, sizeof(buffer), fp);
      if (feof(fp)) {
        free(data->basis_set);
        data->basis_set=NULL;
        data->num_shells_per_atom=NULL;
        data->have_basis=FALSE;
        vmdcon_printf(VMDCON_ERROR, "gaussianplugin) EOF in data base "
                      "file %s while looking for element %s.\n", filepath, 
                      get_pte_label(data->initatoms[i].atomicnum), buffer);
        fclose(fp);
        return FALSE;
      }
      n=sscanf(buffer,"%s%s",word[0],word[1]);
    } while ( (n != 2) || strcmp(word[1],"0") ||
              (strcmp(word[0],get_pte_label(data->initatoms[i].atomicnum))) );
    
    strncpy(data->basis_set[i].name,data->gbasis,sizeof(data->gbasis));
    
    numshells=0;
    shell=NULL;

    /* read basis set until end of element */
    do {

      fgets(buffer, sizeof(buffer), fp);
      if (strstr(buffer,"****")) break;
      if (ferror(fp)) {
        vmdcon_printf(VMDCON_ERROR, "gaussianplugin) read error in data "
                      "base file %s while reading basis of element %s.\n", 
                      filepath, get_pte_label(data->initatoms[i].atomicnum));
        free(data->basis_set);
        data->basis_set=NULL;
        return FALSE;
      }

      numread=sscanf(buffer,"%s%d%f",word[0],&numprim,&scalef);
      if (numread == 3) {
#if GAUSSIAN_DEBUG && GAUSSIAN_BASIS_DEBUG
        vmdcon_printf(VMDCON_INFO, "gaussianplugin) atom: %d, element: %s, shell: %d "
                      "%s-type shell, %d primitives, scalefactor %f\n", i,
                      get_pte_label(data->initatoms[i].atomicnum), numshells+1, 
                      word[0], numprim, scalef);
#endif
        ;
      } else {
        vmdcon_printf(VMDCON_WARN, 
                      "gaussianplugin) basis set parse error: %s",buffer);
        free(data->basis_set);
        data->basis_set=NULL;
        return FALSE;
      }

      ishell=numshells;
      ++numshells;
      shell=realloc(shell,numshells*sizeof(shell_t));
      shell[ishell].numprims=numprim;
      shell[ishell].symmetry=shellsymm_int(word[0]);
      shell[ishell].prim = (prim_t *)calloc(numprim,sizeof(prim_t));
      data->num_basis_funcs += numprim;
      if (shell[ishell].symmetry == SP_S_SHELL) {
        ++numshells;
        shell=realloc(shell,numshells*sizeof(shell_t));
        shell[ishell+1].numprims=numprim;
        shell[ishell+1].symmetry=SP_P_SHELL;
        shell[ishell+1].prim = (prim_t *)calloc(numprim,sizeof(prim_t));
        data->num_basis_funcs += numprim;
      }

      for (n=0; n<numprim; ++n) {
        fgets(buffer, sizeof(buffer), fp);
        if (ferror(fp)) {
          vmdcon_printf(VMDCON_ERROR, "gaussianplugin) read error in data "
                        "base file %s while reading basis of element %s.\n", 
                        filepath, get_pte_label(data->initatoms[i].atomicnum));
          free(data->basis_set);
          data->basis_set=NULL;
          return FALSE;
        }
        sscanf(buffer,"%s%s%s", word[0],word[1],word[2]);
        fix_fortran_exp(word[0]);
        shell[ishell].prim[n].exponent=atof(word[0])*scalef*scalef;
        fix_fortran_exp(word[1]);
        shell[ishell].prim[n].contraction_coeff=atof(word[1]);
        if (shell[ishell].symmetry == SP_S_SHELL) {
          shell[ishell+1].prim[n].exponent=shell[ishell].prim[n].exponent;
          fix_fortran_exp(word[2]);
          shell[ishell+1].prim[n].contraction_coeff=atof(word[2]);
        }
      }
    } while(1);

  
    /* store shells in atom */
    data->basis_set[i].numshells = numshells;
    data->basis_set[i].shell = shell;
    
    /* store the total number of basis functions */
    data->num_shells += numshells;
  }

  vmdcon_printf(VMDCON_INFO, 
                "gaussianplugin) Parsed %d uncontracted basis functions. \n",
                data->num_basis_funcs);

  /* allocate and populate flat arrays needed for molfileplugin */
  data->have_basis = TRUE;
  return fill_basis_arrays(data);
}


/**************************************************
 *
 * Convert shell symmetry type from char to int.
 *
 ************************************************ */
static int shellsymm_int(char *symm) {
  int shell_symmetry;

  switch (toupper(symm[0])) {
    case 'S':
      if (symm[1] == '\0') {
        shell_symmetry = S_SHELL;
      } else if (toupper(symm[1]) == 'P') {
        if (symm[2] == '\0') {
          shell_symmetry = SP_S_SHELL;
        } else if (toupper(symm[1]) == 'D') {
          shell_symmetry = SPD_S_SHELL;
        } else {
          shell_symmetry = UNK_SHELL;
        } 
      } else {
        shell_symmetry = UNK_SHELL;
      }
      break;
    case 'L':
      shell_symmetry = SP_S_SHELL;
      break;
    case 'M': 
      shell_symmetry = SP_P_SHELL;
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


/** Populate the flat arrays containing the basis set data. */
static int fill_basis_arrays(gaussiandata *data) {
  int i, j, k;
  int shellcount = 0;
  int primcount = 0 ;
  float *basis;
  int *num_shells_per_atom;
  int *num_prim_per_shell;
  int *shell_symmetry;

  /* reserve space for pointer to array containing basis
   * info, i.e. contraction coeficients and expansion 
   * coefficients; need 2 entries per primitive gaussian, i.e.
   * exponent and contraction coefficient; also,
   * allocate space for the array holding the orbital symmetry
   * information per primitive Gaussian.
   * Finally, initialize the arrays holding the number of 
   * shells per atom and the number of primitives per shell*/
  SAFE_CALLOC(basis,float,2*data->num_basis_funcs);
  SAFE_CALLOC(shell_symmetry,int,data->num_shells);
  SAFE_CALLOC(num_shells_per_atom,int,data->numatoms);
  SAFE_CALLOC(num_prim_per_shell,int,data->num_shells);

  /* place pointers into struct gaussiandata */
  data->basis = basis;
  data->shell_symmetry = shell_symmetry;
  data->num_shells_per_atom = num_shells_per_atom;
  data->num_prim_per_shell = num_prim_per_shell;

  for(i=0; i<data->numatoms; i++) {
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



/** this function extracts the trajectory information
 *  from the output file */
static int get_traj_frame(gaussiandata *data) {
  qm_timestep_t *cur_qm_ts;
  long fpos;
  char buffer[BUFSIZ];
  char word[MOLFILE_BUFSIZ];  
  buffer[0] = '\0';
  word[0] = '\0';

#if GAUSSIAN_DEBUG
  vmdcon_printf(VMDCON_INFO, 
                "gaussianplugin) Timestep %d: =======================\n", 
                data->num_frames_read);
#endif

  fpos=ftell(data->file);       /* XXX: */
  if (!get_coordinates(data->file, data->initatoms, data->numatoms)) {
    vmdcon_printf(VMDCON_WARN, 
                  "gaussianplugin) Couldn't find input orientation coordinates"
                  " for timestep %d\n", data->num_frames_read);
    vmdcon_printf(VMDCON_WARN, 
                  "gaussianplugin) Trying internal coordinates instead.\n");

    fseek(data->file, fpos, SEEK_SET); /* XXX: */
    if(!get_int_coordinates(data->file, data->initatoms, data->numatoms)) {
      vmdcon_printf(VMDCON_ERROR, 
                    "gaussianplugin) Couldn't find any coordinates.\n", 
                    data->num_frames_read);
      return FALSE;
    } 
  }
  
  /* allocate more memory for the timestep array */
  data->qm_timestep = 
    (qm_timestep_t *)realloc(data->qm_timestep, 
                             (data->num_frames_read+1)*sizeof(qm_timestep_t));

  /* get a convenient pointer to the current qm timestep */
  cur_qm_ts = data->qm_timestep+data->num_frames_read;
  memset(cur_qm_ts, 0, sizeof(qm_timestep_t));

  /* Read the basis set. if available */
  if (data->num_frames_read == 0) {
    get_basis(data);
    get_internal_basis(data);
  }
  
  /* Try to read wavefunction and orbital energies */
  if (get_wavefunction(data, cur_qm_ts) == FALSE) {
    vmdcon_printf(VMDCON_WARN, "gaussianplugin) No wavefunction present for timestep %d\n", data->num_frames_read);
  } else {
    vmdcon_printf(VMDCON_INFO, "gaussianplugin) Wavefunction found for timestep %d\n", data->num_frames_read);
  }

#if 0
  if (get_population(data, cur_qm_ts)) {
    vmdcon_printf(VMDCON_INFO, "gaussianplugin) Mulliken charges found\n");
  }
#endif

  data->num_frames_read++;

  return TRUE;
}


/* Look for the "EQUILIBRIUM GEOMETRY LOCATED" line thereby
 * advancing the file pointer so that the final info block
 * can be parsed.
 * If we don't find this line before the next geometry
 * the file pointer will be set back to where the search
 * started. */
static int find_traj_end(gaussiandata *data) {
  char buffer[BUFSIZ];
  long filepos;
  filepos = ftell(data->file);

  while (1) {
    if (!fgets(buffer, sizeof(buffer), data->file)) break;

    if (strstr(buffer, "Berny optimization.")) {
          data->num_frames++;
    } else if (strstr(buffer, "Optimization completed.")) {
#if GAUSSIAN_DEBUG
      vmdcon_printf(VMDCON_INFO, "gaussianplugin) ==== End of trajectory. ====\n");
#endif
      data->opt_status = STATUS_CONVERGED;
      return TRUE;
    } else if (strstr(buffer, "Optimization stopped.")) {
#if GAUSSIAN_DEBUG
      vmdcon_printf(VMDCON_INFO, "gaussianplugin) ==== End of trajectory. ====\n");
#endif
      data->opt_status = STATUS_TOO_MANY_STEPS;
      return TRUE;
    } else if (strstr(buffer, "Convergence failure -- run terminated.")) {
      data->opt_status = STATUS_SCF_NOT_CONV;
      return FALSE;
    } 
  }

  /* We didn't find any of the regular key strings,
   * the run was most likely broken off and we have an
   * incomplete file. */
  data->opt_status = STATUS_BROKEN_OFF;

  fseek(data->file, filepos, SEEK_SET);
  return FALSE;  
}

/*********************************************************
 *
 * this function reads the actual wavefunction, which is
 * punched at the end of the log file
 *
 **********************************************************/
static int get_wavefunction(gaussiandata *data, qm_timestep_t *ts)
{
  float *orbital_energies;
  float *wave_function;
  char buffer[BUFSIZ];
  char word[6][MOLFILE_BUFSIZ];
  int orbital_counter = 0;
  int i = 0, j = 0, num_values = 0, num_orbs =0;
  long filepos;
  int uhfmult = 1;

  /* no point in searching for wavefunction info, if there cannot be any. */
  if (!data->have_wavefunction) return FALSE;
  
  buffer[0] = '\0';
  for (i=0; i<6; i++) word[i][0] = '\0';

  /*
   * the following output requires  IOP(6/7=3) in the route section.
   *
   * Scan for something like this (g98 closed shell):
     Molecular Orbital Coefficients
                           1         2         3         4         5
                           O         O         O         O         O
     EIGENVALUES --    -7.61925  -7.61868  -0.88486  -0.63877  -0.55750
   1 1   B  1S          0.71623  -0.69260  -0.13195  -0.12199   0.00000
   2        2S          0.01884  -0.01744   0.20755   0.19054   0.00000
   3        2PX         0.00000   0.00000   0.00000   0.00000   0.24930

   * or this (g03 closed shell):

     Molecular Orbital Coefficients
                           1         2         3         4         5
                       (SGG)--O  (SGU)--O  (SGG)--O  (PIU)--O  (PIU)--O
     EIGENVALUES --    -1.30558  -1.16543  -0.63335  -0.56287  -0.51096
   1 1   O  1S          0.68764   0.70168   0.16476   0.00000   0.00000
   2        1PX         0.00000   0.00000   0.00000   0.70711   0.00000
   3        1PY         0.00000   0.00000   0.00000   0.00000   0.70711


   * or this (g03 open shell):

   Alpha Molecular Orbital Coefficients
                           1         2         3         4         5
                       (SGG)--O  (SGU)--O  (SGG)--O  (SGU)--O  (PIU)--O
     EIGENVALUES --   -20.82303 -20.82274  -1.55752  -1.30463  -0.78343
   1 1   O  1S          0.70309   0.70312  -0.15839  -0.17089   0.00000
   2        2S          0.01545   0.01540   0.38998   0.42214   0.00000


   */

  /* remember position in order to go back if no wave function was found */
  filepos = ftell(data->file);

  do {
    GET_LINE(buffer, data->file);
    /* check if we searched too far. XXX need more tests here. */
    if (strstr(buffer,"Input orientation") ) {
      fseek(data->file, filepos, SEEK_SET);
      return FALSE;
    }
  } while(!strstr(buffer,"Molecular Orbital Coefficients"));
  sscanf(buffer,"%s",word[0]);
  if (strcmp(word[0],"Alpha") == 0) uhfmult=2;

  /* Reserve space for arrays storing wavefunction and orbital
   * energies. For the wavefunction we reserve the number
   * of alpha orbitals squared, for unrestricted twice as much.
   * Accordingly, for the energies I use the number of A orbitals 
   * In gaussian the number of alpha and beta orbitals is always 
   * the same. */
  SAFE_CALLOC(wave_function,float,uhfmult*data->wavef_size*data->wavef_size);
  SAFE_CALLOC(orbital_energies,float,uhfmult*data->wavef_size);

  /* store the pointers in gaussiandata */
  ts->wave_function    = wave_function;
  ts->orbital_energies = orbital_energies;

  while (orbital_counter < uhfmult*data->wavef_size) {
    /* read up to line of orbital energies */
    GET_LINE(buffer, data->file);           /* orbital index */
    GET_LINE(buffer, data->file);           /* XXX: occ/virt (+ orbital symm) */
    GET_LINE(buffer, data->file);           /* eigenvalues */
    num_orbs = sscanf(buffer,"%*s%*s%s%s%s%s%s",word[0],
                      word[1],word[2],word[3],word[4]);
    
    /* store the orbital energies in the appropriate arrays 
     * read them until we encounter an empty string */
    for(i=0; i<num_orbs; i++) 
      *(orbital_energies+i) = atof(word[i]);
      
    /* increase orbital energy pointer */
    orbital_energies = orbital_energies+num_orbs;
    orbital_counter += num_orbs;

    /* read in the wavefunction */
    for (i=0; i<data->wavef_size; i++) {
      int xexp=0, yexp=0, zexp=0;

      /* read in the wavefunction coefficients for up 
       * to 5 orbitals at a time line by line */
      GET_LINE(buffer, data->file);
      num_values = sscanf(buffer+12,"%4s%s%s%s%s%s", 
                          word[0], word[1], word[2],
                          word[3], word[4], word[5]);

      /* handle magenetic quantum number. in cartesian basis the
       * labels are: S, PX, PY, PZ, DXX, DXY, DXZ, DYY, DYZ, DZZ, ...*/
      for (j=1; j<strlen(word[0]); j++) {
        switch (word[0][j]) {
          case 'X':
            xexp++;
            break;
          case 'Y':
            yexp++;
            break;
          case 'Z':
            zexp++;
            break;
            /* if we have pure d/f-functions the nomenclature changes to 
             * 'D 0', 'D-1', 'D+1', 'D-2', 'D+2' */
          case '+': /* fallthrough */
          case '-': /* fallthrough */
          case '0': /* fallthrough */
          case '1': /* fallthrough */
          case '2': /* fallthrough */
            vmdcon_printf(VMDCON_ERROR, "gaussianplugin) pure basis function "
                          "detected: '%s'. bailing out...\n", word[0]);
            free(ts->wave_function);
            ts->wave_function=NULL;
            free(ts->orbital_energies);
            ts->orbital_energies=NULL;
            return FALSE;
            break;
          default:
            /* do nothing */
            break;
        }
      }
      data->angular_momentum[3*i  ] = xexp;
      data->angular_momentum[3*i+1] = yexp;
      data->angular_momentum[3*i+2] = zexp;
#if GAUSSIAN_DEBUG && GAUSSIAN_BASIS_DEBUG
      fprintf(stderr,"%s:%d orbital %d/%d  function %d/%s: %d %d %d\n", 
              __FILE__, __LINE__, num_orbs, orbital_counter, 
              i, word[0], xexp, yexp, zexp);
#endif

      /* each orbital has data->wavef_size entries, 
       * hence we have to use this number as offset when storing 
       * them in groups of five */
      for (j=0 ; j<num_values-1; j++) {
        wave_function[j*data->wavef_size+i] = atof(&word[j+1][0]);
      }
    }
    /* move wavefunction pointer to start of next five orbitals */
    wave_function = wave_function + num_orbs*data->wavef_size;
    if ((uhfmult==2) && (orbital_counter==data->wavef_size)) {
      GET_LINE(buffer, data->file);
    }
  }

  /* store the number of orbitals read in */
  ts->orbital_counter = orbital_counter;

  vmdcon_printf(VMDCON_INFO, 
                "gaussianplugin) Number of orbitals scanned: %d \n",
                ts->orbital_counter);

  return TRUE;
}

/* Read the population analysis section.
 * Currently we parse only the Mulliken charges
 * but we might want to add support for populations
 * and for Lowdin analysis. */
static int get_population(gaussiandata *data, qm_timestep_t *ts) {
  int i;
  char buffer[BUFSIZ];
  data->have_mulliken = FALSE;


  /* Read Mulliken charges if present */
  ts->mulliken_charges = 
    (double *)calloc(data->numatoms, sizeof(double));

  if (!ts->mulliken_charges) {
    PRINTERR; 
    return FALSE;
  }
  
  for (i=0; i<data->numatoms; i++) {
    int n;
    float mullpop, mullcharge, lowpop, lowcharge;
    GET_LINE(buffer, data->file);
    n = sscanf(buffer,"%*i%*s%f%f%f%f",
                   &mullpop, &mullcharge, &lowpop, &lowcharge);
    if (n!=4) return FALSE;
    ts->mulliken_charges[i] = mullcharge;
  }

  if (i!=data->numatoms) return FALSE;

  data->have_mulliken = TRUE;
  return TRUE;
}


/*************************************************************
 *
 * plugin registration 
 *
 **************************************************************/

VMDPLUGIN_API int VMDPLUGIN_init(void) {
  memset(&plugin, 0, sizeof(molfile_plugin_t));
  plugin.abiversion = vmdplugin_ABIVERSION;
  plugin.type = MOLFILE_PLUGIN_TYPE;
  plugin.name = "gaussian";
  plugin.prettyname = "Gaussian Logfile (g94,g98,g03)";
  plugin.author = "Axel Kohlmeyer, Markus Dittrich, Jan Saam";
  plugin.majorv = 0;
  plugin.minorv = 1;
  plugin.is_reentrant = VMDPLUGIN_THREADUNSAFE;
  plugin.filename_extension = "log";
  plugin.open_file_read = open_gaussian_read;
  plugin.read_structure = read_gaussian_structure;
  plugin.close_file_read = close_gaussian_read;

  plugin.read_qm_metadata = read_gaussian_metadata;
  plugin.read_qm_rundata  = read_gaussian_rundata;

#if vmdplugin_ABIVERSION > 11
  plugin.read_timestep_metadata    = read_timestep_metadata;
  plugin.read_qm_timestep_metadata = read_qm_timestep_metadata;
  plugin.read_timestep = read_timestep;
#endif

  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)(void *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_fini(void) {
  return VMDPLUGIN_SUCCESS;
}


#ifdef TEST_PLUGIN

int main(int argc, char *argv[]) {
  int numatoms, i, optflags;
  molfile_atom_t *atoms;
  molfile_timestep_t timestep;
  molfile_qm_metadata_t qm_metadata;
  molfile_qm_timestep_t qm_ts;
    
  void *v;

  while (--argc) {
    ++argv;
    v = open_gaussian_read(*argv, "log", &numatoms);
    if (!v) {
      fprintf(stderr, "open_gaussian_read failed for file %s\n", *argv);
      return 1;
    }
    fprintf(stderr, "open_gaussian_read succeeded for file %s\n", *argv);
    fprintf(stderr, "number of atoms: %d\n", numatoms);
    atoms = (molfile_atom_t *)malloc(sizeof(molfile_atom_t)*numatoms);
    read_gaussian_structure(v,&optflags, atoms);

    i = 0;
    timestep.coords = (float *)malloc(3*sizeof(float)*numatoms);
    while (!read_timestep(v, numatoms, &timestep, &qm_metadata, &qm_ts)) {
      i++;
    }
    fprintf(stderr, "ended read_timestep on frame %d\n", i);
    close_gaussian_read(v);
  }
  return 0;
}
#endif
