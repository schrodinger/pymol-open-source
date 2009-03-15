/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_cpmdlogplugin
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
 *          CPMD output file Plugin
 *
 * This plugin allows VMD to read CPMD output files.
 * So far only supports PROPERTIES and specially
 * modified MOLECULAR DYNAMICS BO runs.
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
#define CPMDLOG_DEBUG 0
#define CPMDLOG_BASIS_DEBUG 1
#if CPMDLOG_DEBUG
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
#if CPMDLOG_DEBUG
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

/** Top level CPMD log file parser. Responsible 
 *  for static, i.e. non-trajectory information. */
static int parse_static_data(gaussiandata *, int *);

/** Check if the current run is an actual CPMD run; 
 *  returns true/false */
static int have_cpmd(gaussiandata *);

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

/* count the number of readable QM timesteps 
 * and collect other information about the
 * total trajectory. */
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
 * Called by VMD to open the CPMD logfile and get the number
 * of atoms.
 * We are also reading all the static (i.e. non-trajectory)
 * data here since we have to parse a bit to get the atom count
 * anyway. These data will then be provided to VMD by
 * read_cpmdlog_metadata() and read_cpmdlog_rundata().
 *
 * *************************************************************/
static void *open_cpmdlog_read(const char *filename, 
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

  /* check if the file is CPMD format; 
   * if yes parse it, if not exit */
  if (have_cpmd(data)==TRUE) {
    /* if we're dealing with an unsupported CPMD
     * version, we better quit. so far we can test 3.9.x-3.13.x */
    if ((data->version < 30900) || (data->version > 40000)) {
      vmdcon_printf(VMDCON_ERROR,
                    "cpmdlogplugin) CPMD version %s is not "
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
static int read_cpmdlog_structure(void *mydata, int *optflags, 
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
static int read_cpmdlog_metadata(void *mydata, 
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
static int read_cpmdlog_rundata(void *mydata, 
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
/*   sys_data->wavef_size = data->wavef_size; */
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

#if CPMDLOG_DEBUG
  vmdcon_printf(VMDCON_INFO, 
                "cpmdlogplugin) read_qm_timestep_metadata(): %d/%d/%d\n",
                data->num_frames, 
                data->num_frames_read,
                data->num_frames_sent);
#endif

  meta->count = -1; /* Don't know the number of frames yet */
  meta->has_gradient = 0;

  if (data->num_frames_read > data->num_frames_sent) {
    have = 1;
  } else if (data->num_frames_read < data->num_frames) {
#if CPMDLOG_DEBUG
    vmdcon_printf(VMDCON_INFO,
                  "cpmdlogplugin) Probing timestep %d\n", 
                  data->num_frames_read);
#endif
    have = get_traj_frame(data);
  }

  if (have) {
    /* get a pointer to the current qm timestep */
    qm_timestep_t *cur_qm_ts = data->qm_timestep+data->num_frames_sent;
#if CPMDLOG_DEBUG
    vmdcon_printf(VMDCON_INFO,
                  "cpmdlogplugin) Approved timestep %d\n", 
                  data->num_frames_sent);
#endif
    meta->num_scfiter  = 0;
    meta->num_orbitals_per_wavef[0] = cur_qm_ts->orbital_counter;
    meta->wavef_size = data->wavef_size;

  } else {
    meta->num_scfiter = 0;
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

#if CPMDLOG_DEBUG
  vmdcon_printf(VMDCON_INFO,
                "cpmdlogplugin) Sending timestep %d\n", 
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
static void close_cpmdlog_read(void *mydata) {

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

/*! count number of QM dataset frames. */
static int find_traj_end(gaussiandata *data) {
  char buffer[BUFSIZ];
  long filepos;
  filepos = ftell(data->file);

  while (1) {
    if (!fgets(buffer, sizeof(buffer), data->file)) break;

    if (strstr(buffer, "PROJECTION COORDINATES")) {
      data->num_frames++;
    } 
  }
  data->opt_status = STATUS_UNKNOWN;

  fseek(data->file, filepos, SEEK_SET);
  return FALSE;  
}


static int get_final_info(gaussiandata *data) {
  long filepos;
  filepos = ftell(data->file);

  if (data->runtyp == RUNTYP_OPTIMIZE || 
      data->runtyp == RUNTYP_DYNAMICS) {
    /* Try to advance to the end of the geometry
     * optimization or MD. If no regular end is found we
     * won't find any propertiies to read and return. */
    if (!find_traj_end(data)) return FALSE;
  }

#if 0
  if (get_esp_charges(data)) {
    vmdcon_printf(VMDCON_INFO, "gaussianplugin) ESP charges found!\n");
  }
#endif

  fseek(data->file, filepos, SEEK_SET);
  return TRUE; 
}



/********************************************************
 *
 * Main gaussian log file parser responsible for static,  
 * i.e. non-trajectory information.
 *
 ********************************************************/
static int parse_static_data(gaussiandata *data, int *natoms) 
{
  char buffer[BUFSIZ];
  char word[4][MOLFILE_BUFSIZ];
  char *vmdbasis;
  int  i,n;
  int  numatoms, numstates, numelectrons, totalcharge;

  buffer[0] = '\0';
  
  /* set some defaults */
  data->scftyp = SCFTYP_RHF;
  data->runtyp = RUNTYP_UNKNOWN;
  numelectrons = 0;
  data->totalcharge = 0;
  data->multiplicity = 1;
  data->have_basis=FALSE;       
  /* CPMD never outputs basis set info, so we use 
   * VSTO-6G unless overridden by environment */
  vmdbasis = getenv("VMDDEFBASISSET");
  if (vmdbasis == NULL) 
    vmdbasis = "VSTO-6G";

  strncpy(data->gbasis, vmdbasis, sizeof(data->gbasis));
  strncpy(data->basis_string, "Internal ", sizeof(data->basis_string));
  strncat(data->basis_string, vmdbasis, sizeof(data->basis_string) - 10);


  /* try to find job type parameters within the next 100 lines.*/
  for (i=0; i<100; ++i) {
    GET_LINE(buffer, data->file);
    sscanf(buffer,"%s%s%s",word[0],word[1],word[2]);
    if ( (strcmp(word[0],"CALCULATE" ) == 0 &&
          strcmp(word[1],"SOME"      ) == 0 && 
          strcmp(word[2],"PROPERTIES") == 0 ) ) {
      data->runtyp = RUNTYP_PROPERTIES;
      break;
    } else if ( (strcmp(word[0],"GEOMETRY"    ) == 0 &&
                 strcmp(word[1],"OPTIMIZATION") == 0 ) ) {
      data->runtyp = RUNTYP_OPTIMIZE;
      break;
    } else if ( (strcmp(word[1],"MOLECULAR"    ) == 0 &&
                 strcmp(word[2],"DYNAMICS") == 0 ) ) {
      data->runtyp = RUNTYP_DYNAMICS;
      break;
      
    }
  }

  /* XXX: add support for other types later? */
  if ((data->runtyp != RUNTYP_PROPERTIES) &&
      (data->runtyp != RUNTYP_DYNAMICS) ) {
    vmdcon_printf(VMDCON_ERROR, "cpmdlogplugin) Run-type: %s is currently "
                  "not supported by this plugin.\n", runtypes[data->runtyp]);
    return FALSE;
  }
  
  if ((data->runtyp != RUNTYP_PROPERTIES) && (data->version < 31303)) {
    vmdcon_printf(VMDCON_ERROR, "cpmdlogplugin) Run-type: %s is currently "
                  "not supported for outputs of this CPMD version.\n", 
                  runtypes[data->runtyp]);
    return FALSE;
  }
  
  /* scavange for more setup information */
  do {

    GET_LINE(buffer, data->file);
    n = sscanf(buffer,"%s%s%s%s",word[0],word[1],word[2],word[3]);

    /* empty line */
    if (n < 0) continue;

    /* atom types and initial coordinates */
    if ( (strstr(word[0],"************") != 0 &&
          strcmp(word[1],"ATOMS"       ) == 0 && 
          strstr(word[2],"************") != 0 ) ) {
      /* NR TYPE X(bohr) Y(bohr) Z(bohr) MBL */
      GET_LINE(buffer, data->file);

      numatoms=0;
      data->initatoms=NULL;
      while (1) {
        qm_atom_t *atm;
      
        GET_LINE(buffer, data->file);
        /* end of ATOMS block */
        if (strstr(buffer, "*************************") != NULL)
          break;
        
        data->initatoms=realloc(data->initatoms,(numatoms+1)*sizeof(qm_atom_t));
        atm = data->initatoms + numatoms;
        
        n=sscanf(buffer,"%*d%s%g%g%g", atm->type, &atm->x, &atm->y, &atm->z);
        if (n != 4) {
          free(data->initatoms);
          data->initatoms=NULL;
          vmdcon_printf(VMDCON_ERROR, "cpmdlogplugin) Failed to parse initial" 
                        " coordinates. Stopping.\n");
          return FALSE;
        }
        atm->atomicnum=get_pte_idx(atm->type);
        atm->x *= BOHR_TO_ANGS;
        atm->y *= BOHR_TO_ANGS;
        atm->z *= BOHR_TO_ANGS;
        ++numatoms;
      };
      data->numatoms = numatoms;
      *natoms = numatoms;
      
    } else if ( (strcmp(word[0],"NUMBER" ) == 0 &&
                 strcmp(word[2],"STATES:") == 0 ) ) {
      numstates=atoi(word[3]);

    } else if ( (strcmp(word[0],"NUMBER" ) == 0 &&
                 strcmp(word[2],"ELECTRONS:") == 0 ) ) {
      numelectrons=(int) (atof(word[3]) + 0.5); /* XXX */
      data->num_electrons = numelectrons;
      
    } else if (strcmp(word[0],"CHARGE:" ) == 0 ) {
      totalcharge =(int) (atof(word[1]) + 0.5); /* XXX */
      data->totalcharge = totalcharge;
      
    } else if (strcmp(word[0],"OCCUPATION" ) == 0 ) {
      ; /* XXX. */
      
    } else if ( (strcmp(word[0],"CELL" ) == 0 &&
                 strcmp(word[1],"DIMENSION:") == 0 ) ) {
      n=sscanf(buffer,"%*s%*s%g%g%g%g%g%g",data->initcell,data->initcell+1, 
               data->initcell+2,data->initcell+3,data->initcell+4,data->initcell+5);
      
    /*  
        PROJECT WAVEFUNCTION ON ATOMIC ORBITALS EVERY           10 STEPS
    */
    } else if ( (strcmp(word[0],"PROJECT" ) == 0 &&
                 strcmp(word[1],"WAVEFUNCTION") == 0 ) ) {
        data->have_wavefunction=1;
    /*
      DIPOLE MOMENT CALCULATION 
      STORE DIPOLE MOMENTS EVERY                          10 STEPS
      WANNIER FUNCTION DYNAMICS 
    */
    } else if ( (strcmp(word[0],"WANNIER" ) == 0 &&
                 strcmp(word[1],"FUNCTION") == 0 ) ) {
        if (data->have_wavefunction) {
            data->have_wavefunction=2;
        }
    }
  } while( strcmp(word[0],"INITIALIZATION") || 
           strcmp(word[1],"TIME:") );

  if (numstates >= numelectrons) 
    data->scftyp = SCFTYP_UHF;

  switch (data->scftyp) {
    case SCFTYP_RHF:
      data->num_orbitals_A = numstates;
      data->num_orbitals_B = 0;
      break;
    case SCFTYP_UHF:            /* XXX: this is most likely wrong. check! */
      data->num_orbitals_A = numstates/2;
      data->num_orbitals_B = numstates/2;
      break;
    default:
      break;
  }
  
  vmdcon_printf(VMDCON_INFO, 
                "cpmdlogplugin) Atoms: %d   Charge: %d   Multiplicity: %d\n", 
                data->numatoms, data->totalcharge, data->multiplicity);

  vmdcon_printf(VMDCON_INFO, "cpmdlogplugin) Run-type: %s, SCF-type: %s\n",
                runtypes[data->runtyp], scftypes[data->scftyp]);
  vmdcon_printf(VMDCON_INFO, 
                "cpmdlogplugin) using %s basis set.\n", data->basis_string);

  read_first_frame(data);

  get_final_info(data);
  
  vmdcon_printf(VMDCON_INFO, "cpmdlogplugin) found %d QM data frames.\n", data->num_frames);
#if CPMDLOG_DEBUG
  vmdcon_printf(VMDCON_INFO, "cpmdlogplugin) num_frames_read = %d\n", 
                data->num_frames_read);
  vmdcon_printf(VMDCON_INFO, "cpmdlogplugin) num_frames_sent = %d\n", 
                data->num_frames_sent);
#endif
  return TRUE;
}

/**********************************************************
 *
 * this subroutine checks if the provided files is
 * actually a CPMD file and gathers its version code.
 *
 **********************************************************/
static int have_cpmd(gaussiandata *data) 
{
  char word[4][MOLFILE_BUFSIZ];
  char buffer[BUFSIZ];
  char *ptr;
  int i = 0;
 
  buffer[0] = '\0';
  for (i=0; i<3; i++) word[i][0] = '\0';

  /* check if the file is CPMD format 
   * CPMD output typically begins with:
   *  'PROGRAM CPMD STARTED'
   */
  i=0; /* check only the first 100 lines */
  do {
    GET_LINE(buffer, data->file);
    sscanf(buffer,"%s%s%s",word[0],word[1],word[2]);
    ++i;
  } while( (strcmp(word[0],"PROGRAM") || 
            strcmp(word[1],"CPMD") || 
            strcmp(word[2],"STARTED")) && (i<100) );
  if (i>=100) return FALSE;
  vmdcon_printf(VMDCON_INFO, "cpmdlogplugin) Analyzing CPMD log file: %s\n",data->file_name);
  
  /* now read on until we find the block of text with encoded version
   * number and compile date. */
  i=0; /* check only the next 100 lines */
  do {
    GET_LINE(buffer, data->file);
    sscanf(buffer,"%s%s",word[0],word[1]);
  } while ( (i<100) &&  (strcmp(word[0],"VERSION")) );
  if (i>=100) return FALSE;

  strcpy(data->version_string,word[1]);

  /* now split version number strings */
  ptr=strtok(word[1],"._");
  data->version = 10000*atoi(ptr);
  ptr=strtok(NULL,"._");
  data->version += 100*atoi(ptr);
  ptr=strtok(NULL,"._");
  data->version += atoi(ptr);
  
  vmdcon_printf(VMDCON_INFO, 
                "cpmdlogplugin) CPMD version = %s  (Version code: %d)\n",
                data->version_string, data->version);
  return TRUE;
}

/*! read coordinates for old-style population analysis */
static int read_first_frame(gaussiandata *data) {
  
  data->qm_timestep = NULL;

  /* Read the basis set. */
  get_internal_basis(data);

  /* the angular momentum is populated in get_wavefunction 
   * which is called by get_traj_frame(). We have obtained
   * the array size wavef_size already from the basis set
   * statistics */
  vmdcon_printf(VMDCON_INFO, 
                "cpmdlogplugin) Allocating data for %d wavefunctions\n",
                data->wavef_size);
  SAFE_CALLOC(data->angular_momentum,int,3*data->wavef_size);

  if (data->version < 31303) {

    vmdcon_printf(VMDCON_WARN, 
                  "cpmdlogplugin) ##################################"
                  "####################################\n");
    vmdcon_printf(VMDCON_WARN, 
                  "cpmdlogplugin) This version of CPMD does not print "
                  "the actual coordinates on properties runs.\n");
    vmdcon_printf(VMDCON_WARN, 
                  "cpmdlogplugin) Using initial coordinates from the "
                  "input file instead.\n");
    vmdcon_printf(VMDCON_WARN,
                  "cpmdlogplugin) These coodinates are most likely "
                  "inconsistent with the projection.\n Try the following:\n\n"
                  "set cur [atomselect top all]\n"
                  "set new [atomselect [mol new GEOMETRY.xyz] all]\n"
                  "$cur set {x y z} [$new get {x y z}]\n"
                  "$cur delete; mol delete [$new molid]; $new delete\n\n");
    vmdcon_printf(VMDCON_WARN, 
                  "cpmdlogplugin) ##################################"
                  "####################################\n");
    /* Try to read wavefunction and orbital energies */
    SAFE_CALLOC(data->qm_timestep,qm_timestep_t,1);
    if (get_wavefunction(data, data->qm_timestep) == FALSE) {
      vmdcon_printf(VMDCON_INFO, "cpmdlogplugin) No wavefunction present for timestep %d\n", data->num_frames_read);
      free(data->qm_timestep);
      data->qm_timestep=NULL;
    } else {
      vmdcon_printf(VMDCON_INFO, "cpmdlogplugin) Wavefunction found.\n");
      data->num_frames_read=1;
      data->num_frames=1;
    }
  } else {
    /* don't copy coordinates yet. we have the newer version
     * that prints them right before each projection. */
    data->num_frames = 0;
  }
  
  return TRUE;
}


/*******************************************************
 *
 * this function reads in the basis set data from 
 * <basis>.gbs or $VMDDIR/basis/<basis>.gbs
 *
 * XXX: this is the same function as in gaussianplugin
 * ******************************************************/
int get_internal_basis(gaussiandata *data) {

  char *vmddir=NULL;
  FILE *fp;
  char buffer[BUFSIZ];
  char word[3][MOLFILE_BUFSIZ];
  char filepath[256];
  int  i,n, wavef_size; 

  /* no point in adding a basis set if we already have this information */
  if (data->have_basis) return TRUE;

  /* try to open basis set database file. a file in the current
   * directory takes priority over what is shipped with VMD. */
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
    vmdcon_printf(VMDCON_ERROR, "cpmdlogplugin) failed to read basis set "
                  "from data base file %s\n", filepath);
    data->num_shells_per_atom=NULL;
    data->have_basis=FALSE;
    return FALSE;
  } else {
    vmdcon_printf(VMDCON_INFO, "cpmdlogplugin) reading basis set "
                  "from data base file %s\n", filepath);
  }
  data->wavef_size=0;
  
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
        vmdcon_printf(VMDCON_ERROR, "cpmdlogplugin) EOF in data base "
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
        vmdcon_printf(VMDCON_ERROR, "cpmdlogplugin) read error in data "
                      "base file %s while reading basis of element %s.\n", 
                      filepath, get_pte_label(data->initatoms[i].atomicnum));
        free(data->basis_set);
        data->basis_set=NULL;
        return FALSE;
      }

      numread=sscanf(buffer,"%s%d%f",word[0],&numprim,&scalef);
      if (numread == 3) {
#if CPMDLOG_DEBUG && CPMDLOG_BASIS_DEBUG
        vmdcon_printf(VMDCON_INFO, "cpmdlogplugin) atom: %d, element: %s, shell: %d "
                      "%s-type shell, %d primitives, scalefactor %f\n", i,
                      get_pte_label(data->initatoms[i].atomicnum), numshells+1, 
                      word[0], numprim, scalef);
#endif
        ;
      } else {
        vmdcon_printf(VMDCON_ERROR, 
                      "cpmdlogplugin) basis set parse error: %s",buffer);
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

      switch(shell[ishell].symmetry) {
        case S_SHELL:
          data->wavef_size += 1;
          break;
        case P_SHELL:
          data->wavef_size += 3;
          break;
        case SP_S_SHELL:
          data->wavef_size += 4;
          break;
        case D_SHELL:
          data->wavef_size += 5;  /* XXX: handle pure vs. cartesian */
          break;
        case SPD_S_SHELL:
          data->wavef_size += 9;  /* XXX: handle pure vs. cartesian */
          break;
        case F_SHELL:
          data->wavef_size += 7;  /* XXX: handle pure vs. cartesian */
          break;
        default:
          break;
      }
      
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
          vmdcon_printf(VMDCON_ERROR, "cpmdlogplugin) read error in data "
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
                "cpmdlogplugin) Parsed %d uncontracted basis functions. \n", 
                data->num_basis_funcs);

  /* allocate and populate flat arrays needed for molfileplugin */
  data->have_basis = TRUE;
  return fill_basis_arrays(data);
}


/**************************************************
 * Convert shell symmetry type from char to int.
 * XXX: same function exists in gaussianplugin.
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
   * information per primitive gaussian.
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
  char buffer[BUFSIZ];
  char word[4][MOLFILE_BUFSIZ];
  int n,i;
  

  buffer[0] = '\0';
  word[0][0] = '\0';
  word[1][0] = '\0';
  word[2][0] = '\0';
  word[3][0] = '\0';

#if CPMDLOG_DEBUG
  vmdcon_printf(VMDCON_INFO, 
                "cpmdlogplugin) Timestep %d: ====\n", 
                data->num_frames_read);
#endif

  /* allocate more memory for the timestep array */
  data->qm_timestep = 
    (qm_timestep_t *)realloc(data->qm_timestep, 
                             (data->num_frames_read+1)*sizeof(qm_timestep_t));

  /* get a convenient pointer to the current qm timestep */
  cur_qm_ts = data->qm_timestep+data->num_frames_read;
  memset(cur_qm_ts, 0, sizeof(qm_timestep_t));

  /* search for data */
  while (1) {

    GET_LINE(buffer,data->file);
    n = sscanf(buffer,"%s%s%s%s",word[0],word[1],word[2],word[3]);
    
    /* empty or uninteresting line */
    if (n < 3) continue;

    /* coordinates relevant for projection.
     * this needs a CPMD version > 3.13.2
 ********************* PROJECTION COORDINATES ********************
   NR   TYPE        X(bohr)        Y(bohr)        Z(bohr)    
    1      B       8.536664      10.525423      10.202766
    2      B      10.195353       8.831898      12.528872
    3      H      10.467571       8.964546       9.816927
    4      H       8.327684      10.371750      12.680045
    5      H      12.222474       9.664531      13.239192
    6      H       9.579094       6.722748      12.777072
    7      H       6.516216       9.875605       9.562778
    8      H       9.404536      12.649641       9.952504
 ****************************************************************
 */
    if ( (strstr(word[0],"************") != 0 &&
          strcmp(word[1],"PROJECTION"  ) == 0 && 
          strcmp(word[2],"COORDINATES" ) == 0 && 
          strstr(word[3],"************") != 0 ) ) {

      /* NR TYPE X(bohr) Y(bohr) Z(bohr) */
      GET_LINE(buffer, data->file);

      if (data->initatoms==NULL) {
        vmdcon_printf(VMDCON_ERROR,"why is initatoms NULL?\n");
        return FALSE;
      }

      for (i=0; i < data->numatoms; ++i) {
        qm_atom_t *atm;
      
        GET_LINE(buffer, data->file);
        atm = data->initatoms + i;
        
        n=sscanf(buffer,"%*d%*s%g%g%g", &atm->x, &atm->y, &atm->z);
        
        if (n != 3) {
          vmdcon_printf(VMDCON_ERROR, "cpmdlogplugin) Failed to parse "
                        "projection coordinates. Stopping.\n");
          return FALSE;
        }
        atm->atomicnum=get_pte_idx(atm->type);
        atm->x *= BOHR_TO_ANGS;
        atm->y *= BOHR_TO_ANGS;
        atm->z *= BOHR_TO_ANGS;
      }
      buffer[0] = '\0';
#if CPMDLOG_DEBUG
      vmdcon_printf(VMDCON_INFO, "cpmdlogplugin) Coordinates found for timestep %d\n", data->num_frames_read);
#endif
    } else if ( (strcmp(word[0],"WAVEFUNCTIONS") == 0) &&
         (strcmp(word[2],"ATOMIC"       ) == 0) &&
         (strcmp(word[3],"ORBITAL"      ) == 0) ) {

      /* Try to read wavefunction and orbital energies */
      if (get_wavefunction(data, cur_qm_ts) == FALSE) {
          vmdcon_printf(VMDCON_WARN, "cpmdlogplugin) No wavefunction present for timestep %d\n", data->num_frames_read);
          /* XXX: add flag to ignore wfn */
      }
      if (data->have_wavefunction == 2) {
#if 0          
          /* Try to read localized wavefunctions  */
          /* XXX: used offset or something. */
          if (get_wavefunction(data, cur_qm_ts+XXX) == FALSE) {
              vmdcon_printf(VMDCON_WARN, "cpmdlogplugin) No localized wavefunction present for timestep %d\n", data->num_frames_read);
          }
#endif
      }
      buffer[0] = '\0';
      break; /* XXX */
    } else if ( (strcmp(word[0],"POPULATION") == 0) &&
                (strcmp(word[1],"ANALYSIS"  ) == 0) &&
                (strcmp(word[3],"PROJECTED" ) == 0) ) {
#if 0
      if (get_population(data, cur_qm_ts)) {
        vmdcon_printf(VMDCON_INFO, "cpmdlogplugin) Mulliken charges found\n");
      }
#endif
      buffer[0] = '\0';
      
    } else if ( strcmp(word[0],"*"     ) == 0 &&
                strcmp(word[1],"TIMING") == 0 && 
                strcmp(word[2],"*"     ) == 0 ) {
      data->end_of_trajectory=FALSE;
      break;
    } else if (feof(data->file)) {
      data->end_of_trajectory=TRUE;
      break;
    }
  }

  /* next timestep or end of file */
  data->num_frames_read++;

  return TRUE;
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
#define ORBSPERBLOCK 8
  char word[ORBSPERBLOCK+1][MOLFILE_BUFSIZ];
  int orbital_counter = 0;
  int i = 0, j = 0, num_values = 0;
  int num_orbs = data->num_orbitals_A + data->num_orbitals_B;

  buffer[0] = '\0';
  for (i=0; i<ORBSPERBLOCK+1; i++) word[i][0] = '\0';
  /*
   * Scan for something like this:
      ORBITAL      1       2       3       4       5       6       7       8
  COMPLETNESS    0.972   0.951   0.972   0.972   0.951   0.972   0.951   0.972
  OCCUPATION     2.000   2.000   2.000   2.000   2.000   2.000   2.000   2.000
  1   B  S      -0.064  -0.163   0.064   0.365   0.163  -0.365   0.163  -0.365
         Px      0.000  -0.170   0.000   0.000  -0.170   0.000  -0.170   0.000
         Pz      0.086   0.227  -0.086   0.181  -0.227  -0.181  -0.227  -0.181
         Py      0.043   0.000   0.043   0.411   0.000   0.411   0.000   0.411
  2   B  S       0.365  -0.163  -0.365  -0.064   0.163   0.064   0.163   0.064
         Px      0.000  -0.170   0.000   0.000  -0.170   0.000  -0.170   0.000
         Pz     -0.181  -0.227   0.181  -0.086   0.227   0.086   0.227   0.086
         Py     -0.411   0.000  -0.411  -0.043   0.000  -0.043   0.000  -0.043
  3   H  S      -0.049   0.200   0.049  -0.049   0.540   0.049   0.540   0.049
  4   H  S      -0.049  -0.540   0.049  -0.049  -0.200   0.049  -0.200   0.049
  5   H  S       0.520   0.078   0.041   0.056  -0.078   0.031  -0.078   0.031
    ...
   */


  /* Reserve space for arrays storing wavefunction and orbital
   * energies. For the wavefunction we reserve the number
   * of alpha orbitals squared, for unrestricted twice as much.
   * Accordingly, for the energies I use the number of A orbitals 
   * In gaussian the number of alpha and beta orbitals is always 
   * the same. */
  SAFE_CALLOC(wave_function,float,data->wavef_size*num_orbs);
  SAFE_CALLOC(orbital_energies,float,num_orbs);
  ts->wave_function    = wave_function;
  ts->orbital_energies = orbital_energies;

  while (orbital_counter < num_orbs) {
    /* read up to line of orbital indices */
    do {
      GET_LINE(buffer, data->file);
      sscanf(buffer,"%s",word[0]);
    } while (strcmp(word[0],"ORBITAL"));

    GET_LINE(buffer, data->file);           /* completeness */
    num_orbs = sscanf(buffer,"%*s%s%s%s%s%s%s%s%s",word[0],word[1],
                      word[2],word[3],word[4],word[5],word[6],word[7]);

    /* we don't have orbital energies here, but we
     * can use it for the completenes parameter. */
    for(i=0; i<num_orbs; i++) 
      *(orbital_energies+i) = atof(word[i]);

    GET_LINE(buffer, data->file);           /* occupation */
             
    /* step orbital energy pointer */
    orbital_energies = orbital_energies+num_orbs;
    orbital_counter += num_orbs;

    /* read in the wavefunction */
    for (i=0; i<data->wavef_size; i++) {
      int xexp=0, yexp=0, zexp=0;

      /* read in the wavefunction coefficients for up 
       * to 8 orbitals at a time line by line */
      GET_LINE(buffer, data->file);
      num_values = sscanf(buffer+7,"%4s%s%s%s%s%s%s%s%s", word[0], 
                          word[1], word[2],word[3], word[4], 
                          word[5], word[6], word[7], word[8]);

      /* handle magenetic quantum number. in cartesian basis the
       * labels are: S, PX, PY, PZ, DXX, DXY, DXZ, DYY, DYZ, DZZ, ...*/
      for (j=1; j<strlen(word[0]); j++) {
        switch (toupper(word[0][j])) {
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
            vmdcon_printf(VMDCON_ERROR, "cpmdlogplugin) pure basis function "
                          "detected: '%s'. bailing out...\n", word[0]);
            free(ts->wave_function);
            wave_function=NULL;
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
#if CPMDLOG_DEBUG && CPMDLOG_BASIS_DEBUG
      fprintf(stderr,"%s:%d orbital %d/%d  function %d/%s: %d %d %d\n", 
              __FILE__, __LINE__, num_orbs, orbital_counter, 
              i, word[0], xexp, yexp, zexp);
#endif

      /* each orbital has data->wavef_size entries, 
       * hence we have to use this number as offset when storing 
       * them in groups of five */
      for (j=0 ; j<num_values-1; j++) {
        wave_function[j*data->wavef_size+i] = atof(word[j+1]);
      }
    }
    /* move wavefunction pointer to start of next five orbitals */
    wave_function = wave_function + num_orbs*data->wavef_size;
    /* XXX: FIXME test with open shell run */
    if ((data->scftyp == SCFTYP_UHF) && (orbital_counter==data->wavef_size)) {
      GET_LINE(buffer, data->file);
    }
  }

  /* store the number of orbitals read in */
  ts->orbital_counter = orbital_counter;

#if CPMDLOG_DEBUG
  vmdcon_printf(VMDCON_INFO, 
                "cpmdlogplugin) Number of orbitals scanned: %d \n",
                ts->orbital_counter);
#endif
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
  plugin.name = "cpmdlog";
  plugin.prettyname = "CPMD 3.x output ";
  plugin.author = "Axel Kohlmeyer";
  plugin.majorv = 0;
  plugin.minorv = 0;
  plugin.is_reentrant = VMDPLUGIN_THREADUNSAFE;
  plugin.filename_extension = "out";
  plugin.open_file_read = open_cpmdlog_read;
  plugin.read_structure = read_cpmdlog_structure;
  plugin.close_file_read = close_cpmdlog_read;

  plugin.read_qm_metadata = read_cpmdlog_metadata;
  plugin.read_qm_rundata  = read_cpmdlog_rundata;

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
    v = open_cpmdlog_read(*argv, "log", &numatoms);
    if (!v) {
      fprintf(stderr, "open_cpmdlog_read failed for file %s\n", *argv);
      return 1;
    }
    fprintf(stderr, "open_cpmdlog_read succeeded for file %s\n", *argv);
    fprintf(stderr, "number of atoms: %d\n", numatoms);
    atoms = (molfile_atom_t *)malloc(sizeof(molfile_atom_t)*numatoms);
    read_cpmdlog_structure(v,&optflags, atoms);

    i = 0;
    timestep.coords = (float *)malloc(3*sizeof(float)*numatoms);
    while (!read_timestep(v, numatoms, &timestep, &qm_metadata, &qm_ts)) {
      i++;
    }
    fprintf(stderr, "ended read_timestep on frame %d\n", i);
    close_cpmdlog_read(v);
  }
  return 0;
}

#endif
