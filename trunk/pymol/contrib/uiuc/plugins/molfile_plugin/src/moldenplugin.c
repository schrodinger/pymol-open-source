/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_moldenplugin
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
 *      $RCSfile: moldenplugin.c,v $
 *      $Author: saam $       $Locker:  $             $State: Exp $
 *      $Revision: 1.32 $       $Date: 2011/06/20 20:54:58 $
 *
 ***************************************************************************/

/* This is a plugin that will read input from a MOLDEN
** generated output file 
** some more details will go here soon 
** NOTE: The current version of the plugin relies
** on the fact that the [Atom] field comes before
** the [Geometries] field */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "molfile_plugin.h"
#include "unit_conversion.h"
#include "periodic_table.h"
#include "qmplugin.h"


#define ALLOCATE(array, type, size) \
  array = (type *)calloc(size, sizeof(type)); \
  if (array == NULL) { \
    fprintf(stderr, "moldenplugin) Memory allocation for %s failed!\n", #array); \
    return FALSE; \
  }

#define GET_LINE(x,y) if (!fgets(x, sizeof(x), y)) return FALSE


/* I could use a flag already present in qmdata_t to 
 * indicate a trajectory, but I'm using moldendata_t 
 * to demonstrate the use of 
 * void *format_specific_data;
 * in qmdata_t as a means to store data specific to 
 * the plugin. 
 */
typedef struct {
  long filepos_atoms;   /* [ATOMS] section */
  long filepos_geomxyz; /* [GEOMETRIES] XYZ section */
  long filepos_gto;     /* [GTO] section */
  long filepos_mo;      /* [MO] section */
  char units[16];
  int coordsonly;
} moldendata_t;



/* Read the basis set data */
static int get_basis (qmdata_t *);
static int shelltype_int(char *type);
static int fill_basis_arrays(qmdata_t *data);

static int read_geom_block(qmdata_t *data);
static int read_molecular_orbitals(qmdata_t *data);
static int read_wave_coeffs(FILE *file, qm_wavefunction_t *wave);
static int count_orbitals(qmdata_t *data);


/*********************************************************
 *
 * Open file and return number of atoms.
 *
 * In order to determine the # atom we have to read into
 * the [ATOMS]/[GEOMETRIES] sections where the atoms are
 * defined.
 * There can be either [ATOMS] or [GEOMETRIES] or both.
 * While [GEOMETRIES] has the # atoms as its first line,
 * we actually have to count lines in [ATOMS].
 *
 * We assume no particular order for the sections and scan
 * the entire file for the according keywords. In order to
 * save time when the section contents are actually read
 * we store the file pointers to the beginning of each
 * section in moldendata_t.
 *
 *********************************************************/
static void *open_molden_read(const char *filename,
                              const char *filetype,
                              int *natoms) {
  FILE *fd;
  qmdata_t *data;
  moldendata_t *moldendata;
  char buffer[1024];
  char keystring[20];

  fd = fopen(filename, "rb");
  if (!fd) return NULL;
  
  /* allocate memory for main QM data structure */
  data = init_qmdata(data);
  if (!data) return NULL;

  data->file = fd;

  /* allocate GAMESS specific data */
  moldendata = (moldendata_t *)calloc(1, sizeof(moldendata_t));
  if (!moldendata) return NULL;

  data->format_specific_data = moldendata;


  /* Read first line */
  if (!fgets(buffer,1024,data->file)) return NULL;

  /* Check if the file is MOLDEN format */
  if (!strcmp(strtoupper(trimleft(trimright(buffer))), "[MOLDEN FORMAT]")) {
    printf("moldenplugin) Detected MOLDEN file format!\n");
  } else {
    printf("moldenplugin) The file is not in MOLDEN format!\n");
    return NULL;
  }

  eatwhitelines(data->file);


  /* Identify the different sections. */
  while (fgets(buffer,1024,data->file)) {

    /* Get key string (ignoring empty lines) */
    if (!sscanf(buffer, "%s", keystring)) continue;

    /* Quick test avoiding the uppercase transformation */
    if (keystring[0]!='[') continue;

    /* Make keystring upper case */
    strtoupper(keystring);

    if (!strcmp(keystring, "[5D]") || !strcmp(keystring, "[5D7F]") ||
	!strcmp(keystring, "[7F]") || !strcmp(keystring, "[5D10F]") ||
	!strcmp(keystring, "[9G]")) {
      printf("moldenplugin) Spherical harmonic basis found %s. \n", keystring);
      printf("moldenplugin)   Currently VMD handles only basis sets/wave functions\n");
      printf("moldenplugin)   with cartesian Gaussian functions.\n");
      printf("moldenplugin)   Loading coordinates only.\n");
      moldendata->coordsonly = 1;
    }

    if (!strcmp(keystring, "[ATOMS]")) {
      char *s;
      long prevline=ftell(fd);
      printf("moldenplugin) Found [ATOMS] section ...\n");
      moldendata->filepos_atoms = ftell(data->file);

      if (!sscanf(buffer, "%*s %s", moldendata->units)) {
        printf("moldenplugin) Missing units in [ATOMS] section!\n");
        return NULL;
      }

      /* start counting the atoms; 
       * read until I hit the first line that starts with a "["
       * bracket */      
      (*natoms) = 0; 
      s = fgets(buffer, 1024, fd);

      /* Here we assume that the [ATOMS] section goes
       * on until the next empty line or another section
       * starts, i.e. there is a "[" or we encounter EOF */
      while (trimleft(buffer)[0]!='[' && s!=NULL && !iswhiteline(buffer)) {
        (*natoms)++;
        prevline = ftell(fd);
        s = fgets(buffer, 1024, fd);
      }
      data->numatoms = *natoms;
      data->num_frames = 1;

      /* Go back to the previous line, it might contain
       * the next keyword */
      fseek(data->file, prevline, SEEK_SET);
    }

    else if (!strcmp(keystring, "[GEOMETRIES]")) {
      if (!strcmp(trimright(buffer), "[GEOMETRIES] XYZ")) {
        printf("moldenplugin) Found [GEOMETRIES] XYZ section ...\n");

        moldendata->filepos_geomxyz = ftell(data->file);

        /* The first line of the XYZ type [GEOMETRIES] input
         * contains the number of atoms. */
        if (fscanf(data->file, "%d", natoms) != 1) {
          printf("moldenplugin) No # atoms found in [GEOMTERIES] section!\n");
          return NULL;
        }
        data->numatoms = *natoms;

        /* Jump back to the beginning of the section */
        fseek(data->file, moldendata->filepos_geomxyz, SEEK_SET);

        /* Count # frames */
        data->num_frames = 0;
        do {
          int natm = 0;
          fscanf(data->file, "%d", &natm);
          if (natm!=data->numatoms) break;
          eatline(data->file, 1);

          data->filepos_array = (long*)realloc(data->filepos_array,
                                               (data->num_frames+1)*sizeof(long));
          data->filepos_array[data->num_frames] = ftell(data->file);

          /* Skip title line + numatoms lines */
          eatline(data->file, 1+data->numatoms);
          if (feof(data->file)) break;

          data->num_frames++;
        } while (1);

        printf("moldenplugin) Found %d frames\n", data->num_frames);
      } else if (!strcmp(trimright(buffer), "[GEOMETRIES] ZMAT")) {
        printf("moldenplugin) [GEOMETRIES] ZMAT not supported!\n");
      }
    }

    else if (!strcmp(keystring,"[GTO]")) {
      printf("moldenplugin) Found [GTO] section ...\n");
      moldendata->filepos_gto = ftell(data->file);
    }

    else if (!strcmp(keystring,"[MO]")) {
      printf("moldenplugin) Found [MO] section ...\n");
      moldendata->filepos_mo = ftell(data->file);
    }

  };
  
  return data;
}



/*********************************************************
 *
 * Read geometry from file
 *
 * The [ATOMS] section provides atom name, atomic number
 * and coordinates while the [GEOMETRIES] XYZ section
 * provides atom type and coordinates. Trajectories can
 * only be specified using [GEOMETRIES].
 * In case we only have [GEOMETRIES] the atomic number
 * will have to be deduced from the atom name.
 *
 *********************************************************/
static int read_molden_structure(void *mydata, int *optflags, 
                                 molfile_atom_t *atoms) 
{
  int i;
  char buffer[1024];
  char atname[1024];
  int num, atomicnum;
  molfile_atom_t *atom;
  qmdata_t *data = (qmdata_t *)mydata;
  moldendata_t *moldendata = (moldendata_t *)data->format_specific_data;

  ALLOCATE(data->atoms, qm_atom_t, data->numatoms);

  /* atomic number is provided by plugin.
   * (This is required for QM plugins!) */
  *optflags = MOLFILE_ATOMICNUMBER;

  /* [ATOMS] section */
  if (moldendata->filepos_atoms) { 
    float unitfac = 1.f;

    /* If the units are given in AU we have to convert them.       */
    /* Note: Also recognize parenthesized units emitted by Molcas. */
    if (!strcmp(moldendata->units, "AU") ||
        !strcmp(moldendata->units, "(AU)")) {
      unitfac = BOHR_TO_ANGS;
    }
    
    /* Jump to beginning of [ATOMS] section. */
    fseek(data->file, moldendata->filepos_atoms, SEEK_SET);

    /* Read the atom types, names, atomic numbers
     * as well as x,y,z coordinates */
    for (i=0; i<data->numatoms; i++) {
      float x,y,z;
      atom = atoms+i;
      
      if (!fgets(buffer,1024,data->file)) return MOLFILE_ERROR;
      
      sscanf(buffer,"%s %d %d %f %f %f", atname, &num,
             &atomicnum, &x, &y, &z);

      /* populate data structure for VMD */
      strncpy(atom->name, atname,     sizeof(atom->name)); 
      strncpy(atom->type, atom->name, sizeof(atom->type));
      atom->atomicnumber = atomicnum;
      atom->resname[0] = '\0';
      atom->resid = 1;
      atom->chain[0] = '\0';
      atom->segid[0] = '\0';

      /* keep local copy */
      strncpy(data->atoms[i].type, atname, sizeof(data->atoms[i].type));
      data->atoms[i].atomicnum = atomicnum;
      data->atoms[i].x = x*unitfac;
      data->atoms[i].y = y*unitfac;
      data->atoms[i].z = z*unitfac;
    }
    data->num_frames_read = 1;

    return MOLFILE_SUCCESS;
  }

  /* [GEOMETRIES] XYZ section */
  if (moldendata->filepos_geomxyz) {

    /* Jump to beginning of [GEOMETRIES] section. */
    fseek(data->file, moldendata->filepos_geomxyz, SEEK_SET);
    eatline(data->file, 2);

    /* Read block from file */
    for (i=0; i<data->numatoms; i++) {
      atom = atoms+i;
      if (!fgets(buffer,1024,data->file)) return MOLFILE_ERROR;
      sscanf(buffer,"%s %*f %*f %*f", atname);

      strncpy(atom->type, atname, sizeof(atom->type));
      strncpy(atom->name, atname, sizeof(atom->name)); 
      atom->atomicnumber = get_pte_idx_from_string(atname);
      atom->resname[0] = '\0';
      atom->resid = 1;
      atom->chain[0] = '\0';
      atom->segid[0] = '\0';
      data->atoms[i].atomicnum = atom->atomicnumber;
    }
    data->num_frames_read = 0;

    return MOLFILE_SUCCESS;
  }

  printf("Sorry, could not obtain structure information \n");
  printf("from either the [ATOMS] or [GEOMETRIES] section! \n");
  printf("Please check your MOLDEN output file! \n"); 
  return MOLFILE_ERROR; 
}


/***********************************************************
 *
 * Read atoms for one frame from [GEOMETRIES} section.
 *
 ***********************************************************/
static int read_geom_block(qmdata_t *data) {
  int i;
  char buffer[1024];
  float x,y,z;

  /* Skip title line */
  eatline(data->file, 1);

  for (i=0; i<data->numatoms; i++) {
    if (!fgets(buffer,1024,data->file)) return 0;
    sscanf(buffer,"%*s %f %f %f", &x, &y, &z);
    data->atoms[i].x = x;
    data->atoms[i].y = y;
    data->atoms[i].z = z;
  }

  return 1;
}


/***********************************************************
 *
 * Provide non-QM metadata for next timestep. 
 * Required by the plugin interface.
 *
 ***********************************************************/
static int read_timestep_metadata(void *mydata,
                                  molfile_timestep_metadata_t *meta) {
  meta->count = -1;
  meta->has_velocities = 0;

  return MOLFILE_SUCCESS;
}


/***********************************************************
 *
 * We are not reading the coefficients themselves,
 * because that could require a large amount of memory.
 *
 ***********************************************************/
static int read_qm_timestep_metadata(void *mydata,
                                    molfile_qm_timestep_metadata_t *meta) {
  qmdata_t *data = (qmdata_t *)mydata;
  moldendata_t *moldendata = (moldendata_t *)data->format_specific_data;

  if (data->num_frames_sent >= data->num_frames) {
    /* All frames were sent. */
    return MOLFILE_ERROR;
  }

  /* Can't send metadata if only coordinates were read */
  if (moldendata->coordsonly) return MOLFILE_ERROR;

  /* Count the number of cartesian basis functions in 
     the basis set */
  if (data->num_frames_sent == data->num_frames-1) {
    int i;
    qm_timestep_t *cur_ts;

    if (!count_orbitals(data)) return MOLFILE_ERROR;

    /* get a pointer to the current qm timestep */
    cur_ts = data->qm_timestep;
    
    for (i=0; (i<MOLFILE_MAXWAVEPERTS && i<cur_ts->numwave); i++) {
      meta->num_orbitals_per_wavef[i] = cur_ts->wave[i].num_orbitals;
      meta->has_occup_per_wavef[i]    = cur_ts->wave[i].has_occup;
      meta->has_orben_per_wavef[i]    = cur_ts->wave[i].has_orben;
    }
    meta->wavef_size   = data->wavef_size;
    meta->num_wavef    = cur_ts->numwave;
    meta->num_scfiter  = cur_ts->num_scfiter;
    meta->has_gradient = FALSE;
    meta->num_charge_sets = 0;
  }

  return MOLFILE_SUCCESS;
}



/***********************************************************
 *
 * Provides VMD with the data of the next timestep.
 *
 ***********************************************************/
static int read_timestep(void *mydata, int natoms, 
       molfile_timestep_t *ts, molfile_qm_metadata_t *qm_metadata,
                         molfile_qm_timestep_t *qm_ts) {
  int i;
  qmdata_t *data = (qmdata_t *)mydata;
  qm_timestep_t *cur_ts;

  if (data->num_frames_sent >= data->num_frames) {
    /* All frames were sent. */
    return MOLFILE_ERROR;
  }

  if (data->num_frames_sent == data->num_frames_read) {
    /* Read next coordinate block from file */
    fseek(data->file, data->filepos_array[data->num_frames_read], SEEK_SET);
    read_geom_block(data);

    /*printf("moldenplugin) Read frame %d\n", data->num_frames_read); */
    data->num_frames_read++;
  }


  /* Copy the coordinates */
  for (i=0; i<natoms; i++) {
    ts->coords[3*i  ] = data->atoms[i].x;
    ts->coords[3*i+1] = data->atoms[i].y;
    ts->coords[3*i+2] = data->atoms[i].z; 
  }

  /*printf("moldenplugin) Sent frame %d\n", data->num_frames_sent); */
  data->num_frames_sent++;

  /* In MOLDEN the MOs are listed only for the last frame */
  if (data->num_frames_sent == data->num_frames) {

    /* get a convenient pointer to the current qm timestep */
    cur_ts = data->qm_timestep;

    read_molecular_orbitals(data);

    /* store the wave function and orbital energies */
    if (cur_ts != NULL && cur_ts->wave != NULL) {
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

  }

  return MOLFILE_SUCCESS;
}
  

/*****************************************************
 *
 * Provide VMD with the sizes of the QM related
 * data structure arrays that need to be made
 * available.
 * Since we cannot determine the basis set meta data
 * without parsing the whole basis set section, we
 * read all basis set data here. The data is stored
 * in the qmdata_t structure for later use in
 * read_molden_rundata().
 *
 *****************************************************/
static int read_molden_metadata(void *mydata, 
    molfile_qm_metadata_t *metadata) {

  qmdata_t *data;
  moldendata_t *moldendata;
  data = (qmdata_t *)mydata;
  moldendata = (moldendata_t *)data->format_specific_data;


  metadata->ncart = 0;
  metadata->nimag = 0;
  metadata->nintcoords = 0;

  metadata->have_sysinfo = 0;
  metadata->have_carthessian = 0;
  metadata->have_inthessian = 0;
  metadata->have_normalmodes = 0;

  metadata->num_basis_funcs = 0;
  metadata->num_basis_atoms = 0;
  metadata->num_shells = 0;
  metadata->wavef_size = 0;

  if (!moldendata->coordsonly) {
    /* Read the basis set */
    if (!get_basis(data)) return MOLFILE_ERROR; 

    /* orbital + basis set data */
    metadata->num_basis_funcs = data->num_basis_funcs;
    metadata->num_basis_atoms = data->num_basis_atoms;
    metadata->num_shells      = data->num_shells;
    metadata->wavef_size      = data->wavef_size;  
  }

  return MOLFILE_SUCCESS;
}


/******************************************************
 * 
 * Provide VMD with the static (i.e. non-trajectory)
 * data. That means we are filling the molfile_plugin
 * data structures.
 *
 ******************************************************/
static int read_molden_rundata(void *mydata, 
                               molfile_qm_t *qm_data) {
  qmdata_t *data = (qmdata_t *)mydata;
  int i;
  molfile_qm_hessian_t *hessian_data;
  molfile_qm_basis_t   *basis_data;
  molfile_qm_sysinfo_t *sys_data;

  if (!qm_data) return MOLFILE_ERROR;


  hessian_data = &qm_data->hess;
  basis_data   = &qm_data->basis;
  sys_data     = &qm_data->run;

  sys_data->num_electrons = data->num_electrons;
  sys_data->totalcharge = data->totalcharge;

  /* Populate basis set data */
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

    /* If we have MOs in the file we must provide the 
     * angular momentum exponents */
    if (data->angular_momentum) {
      for (i=0; i<3*data->wavef_size; i++) {
        basis_data->angular_momentum[i] = data->angular_momentum[i];
      }
    }
  }

  /* fill in molfile_qm_sysinfo_t */
  /*sys_data->runtype = data->runtype;
  sys_data->scftype = data->scftype;
  sys_data->nproc   = data->nproc;
  sys_data->num_electrons  = data->num_electrons;
  sys_data->totalcharge    = data->totalcharge;
  sys_data->num_occupied_A = data->num_occupied_A;
  sys_data->num_occupied_B = data->num_occupied_B;
  sys_data->status         = data->opt_status;
  */
  return MOLFILE_SUCCESS;
}


/**********************************************************
 *
 * close file and free memory
 *
 **********************************************************/
static void close_molden_read(void *mydata) {
  int i, j;
  qmdata_t *data = (qmdata_t *)mydata;

  fclose(data->file);

  free(data->atoms);
  free(data->basis);
  free(data->shell_types);
  free(data->atomicnum_per_basisatom);
  free(data->num_shells_per_atom);
  free(data->num_prim_per_shell);
  free(data->angular_momentum);

  if (data->basis_set) {
    for(i=0; i<data->num_basis_atoms; i++) {
      for (j=0; j<data->basis_set[i].numshells; j++) {
        free(data->basis_set[i].shell[j].prim);
      }
      free(data->basis_set[i].shell);
    } 
    free(data->basis_set);
  }

  free(data->format_specific_data);
  free(data->filepos_array);

  if (data->qm_timestep != NULL) {
    for (j=0; j<data->qm_timestep[0].numwave; j++) {
      free(data->qm_timestep[0].wave[j].wave_coeffs);
      free(data->qm_timestep[0].wave[j].orb_energies);
      free(data->qm_timestep[0].wave[j].orb_occupancies);
    }
    free(data->qm_timestep[0].wave);
    free(data->qm_timestep);
  } else {
    printf("close_molden_read(): NULL qm_timestep!\n");
  }

  free(data);
}



/* ####################################################### */
/*             End of API functions                        */
/* The following functions actually do the file parsing.   */
/* ####################################################### */


/******************************************************
 * 
 * Format specification of the basis-set consisting of 
 * contracted Gaussian Type Orbitals.
 *
 * [GTO]
 * atom_sequence_number1 0
 * shell_label number_of_primitives 1.00
 * exponent_prim_1 contraction_coeff_1 (contraction_coeff_sp_1)
 * ...
 * <empty line>
 * atom_sequence_number2 0
 * shell_label number_of_primitives 1.00
 * exponent_prim_1 contraction_coeff_1 (contraction_coeff_1)
 * ...
 * <empty line>
 *
 * recognized shell_labels: s, p, d, f, sp, g
 *
 * For 'sp' shells two contraction coefficients must be given,
 * for both s and p functions. 
 * The 0 on the shell_number line and the 1.00 on the
 * shell_label line are no longer functional and can be
 * ignored.
 *
 * All workings with the [GTO] keyword are in Atomic Units.
 *
 ******************************************************/


/*******************************************************
 *
 * Read the basis set data
 *
 * Format example:
 * [GTO]
 *   1 0
 *  s    2 1.00
 *   0.2738503300E+02  0.4301284983E+00
 *   0.4874522100E+01  0.6789135305E+00
 *  sp   2 1.00
 *   0.1136748200E+01  0.4947176920E-01  0.5115407076E+00
 *   0.2883094000E+00  0.9637824081E+00  0.6128198961E+00
 *
 *   2 0
 *  s    2 1.00
 *   0.1309756400E+01  0.4301284983E+00
 *   0.2331360000E+00  0.6789135305E+00
 *  ...
 *
 * qmdata_t provides hierarchical data structures for 
 * the basis set which are convenient for parsing. 
 * The molfile_plugin interface, however, requires flat
 * arrays, so after reading is done we have to populate
 * the according arrays using fill_basis_arrays().
 * 
 *******************************************************/
static int get_basis(qmdata_t *data) {
  char buffer[1024];
  char shelltype[1024];
  int atomid, numprims;
  int i, j=0;
  int numshells;
  moldendata_t *moldendata = (moldendata_t *)data->format_specific_data;

  /* XXX already initialized in open_molden_read() */
  data->num_shells = 0;
  data->num_basis_funcs = 0;
  data->num_basis_atoms = 0;

  /* initialize basis set the character array */
  memset(data->basis_string, 0, sizeof(data->basis_string));


  /* Place file pointer on line after the [GTO] keyword. */
  fseek(data->file, moldendata->filepos_gto, SEEK_SET);
 
  /* Allocate memory for the basis of all atoms. */
  ALLOCATE(data->basis_set, basis_atom_t, data->numatoms);

  /* Loop over all atoms. */
  for (i=0; i<data->numatoms; i++) {

    if (!fgets(buffer,1024,data->file)) return FALSE;
    sscanf(buffer,"%d %*d", &atomid);

    numshells = 0;
    data->basis_set[i].shell = NULL;

    /* Read an unknown number of shells */
    while (1) {
      shell_t *shell, *shell2=NULL;

      if (!fgets(buffer,1024,data->file)) return FALSE;
      
      /* Empty line signifies beginning of next atom */
      if (!strlen(trimleft(buffer))) break;
      
      /* Get shell type (s, p, d, f, g, sp) and # of primitives */
      sscanf(buffer,"%s %d %*f", shelltype, &numprims);


      /* Add new shell(s). */
      if (!strcmp(shelltype, "sp")) {
        /* Two new shells for SP */
        data->basis_set[i].shell =
          (shell_t*)realloc(data->basis_set[i].shell,
                            (numshells+2)*sizeof(shell_t));
      } else {
        /* One new shell for non-SP */
        data->basis_set[i].shell =
          (shell_t*)realloc(data->basis_set[i].shell,
                            (numshells+1)*sizeof(shell_t));
      }

      shell  = &(data->basis_set[i].shell[numshells]);
      memset(shell, 0, sizeof(shell_t));
      shell->numprims = numprims;
      shell->type     = shelltype_int(shelltype);
      shell->prim     = (prim_t*)calloc(numprims, sizeof(prim_t));

      /* If this is an SP-shell we have to add as separate 
       * S-shell and P-shell. */
      if (!strcmp(shelltype, "sp")) {
        shell->type      = SP_S_SHELL;
        shell2 = &(data->basis_set[i].shell[numshells+1]);
        shell2->numprims = numprims;
        shell2->type     = SP_P_SHELL;
        shell2->prim     = (prim_t*)calloc(numprims, sizeof(prim_t));
      }

      /* Loop over the primitives */
      for (j=0; j<numprims; j++) {
        int nr;
        double expon=0.f, coeff1, coeff2=0.f;
	char s_expon[128], s_coeff1[128], s_coeff2[128];
	if (!fgets(buffer,1024,data->file)) return FALSE;

	/* MOLDEN writes the basis set coefficients using Fortran style notation 
	 * where the exponential character is 'D' instead of 'E'. Other packages 
	 * adhere to C-style notation. Unfortunately sscanf() won't recognize 
	 * Fortran-style numbers. Therefore we have to read the line as string 
	 * first, convert the numbers by replacing the 'D' and then extract the 
	 * floats using sscanf(). */
	fpexpftoc(buffer);
	nr = sscanf(buffer,"%lf %lf %lf", &expon, &coeff1, &coeff2);
        if (nr<2) {
          printf("moldenplugin) Bad format in [GTO] section\n");
          return FALSE;
        }
        shell->prim[j].exponent = expon;
        shell->prim[j].contraction_coeff = coeff1;

        /* P-shell component of SP-shell */
        if (!strcmp(shelltype, "sp")) {
          if (nr!=3) {
            printf("moldenplugin) Bad SP-shell format in [GTO] section\n");
            return FALSE;
          }
          shell2->prim[j].exponent = expon;
          shell2->prim[j].contraction_coeff = coeff2;
        }
      }

      /* Update # uncontracted basis functions */
      data->num_basis_funcs += numprims;

      numshells++;

      /* Account for SP-shells */
      if (!strcmp(shelltype, "sp")) {
        numshells++;
        data->num_basis_funcs += numprims;
      }
    }

    /* Store # shells for current atom */
    data->basis_set[i].numshells = numshells;

    /* Update total number of basis functions */
    data->num_shells += numshells;
  }

  /* As far as I know in Molden format the basis has
   * to be specified for each individual atom, even
   * if the types are the same. */
  data->num_basis_atoms = data->numatoms;

  /* allocate and populate flat arrays needed for molfileplugin */
  fill_basis_arrays(data);

  /* Count the number of cartesian basis functions in 
   * the basis set */
  data->wavef_size = 0;
  for (i=0; i<data->num_shells; i++) {
    switch (data->shell_types[i]) {
    case S_SHELL:
    case SP_S_SHELL:
      data->wavef_size += 1;
      break;
    case P_SHELL:
    case SP_P_SHELL:
      data->wavef_size += 3;
      break;
    case D_SHELL:
      data->wavef_size += 6;
      break;
    case F_SHELL:
      data->wavef_size += 10;
      break;
    case G_SHELL:
      data->wavef_size += 15;
      break;
    }
  }


  /* If we have MOs in the file we must provide the 
   * angular momentum exponents.
   * The order of P, D, F en G functions is as follows:

   *  5D: D 0, D+1, D-1, D+2, D-2
   *  6D: xx, yy, zz, xy, xz, yz

   *  7F: F 0, F+1, F-1, F+2, F-2, F+3, F-3
   * 10F: xxx, yyy, zzz, xyy, xxy, xxz, xzz, yzz, yyz, xyz

   *  9G: G 0, G+1, G-1, G+2, G-2, G+3, G-3, G+4, G-4
   * 15G: xxxx yyyy zzzz xxxy xxxz yyyx yyyz zzzx zzzy,
   *      xxyy xxzz yyzz xxyz yyxz zzxy
   */
  ALLOCATE(data->angular_momentum, int, 3*data->wavef_size);

  for (i=0; i<data->num_shells; i++) {
    switch (data->shell_types[i]) {
    case P_SHELL:
    case SP_P_SHELL:
      angular_momentum_expon(&data->angular_momentum[j  ], "x");
      angular_momentum_expon(&data->angular_momentum[j+3], "y");
      angular_momentum_expon(&data->angular_momentum[j+6], "z");
      j += 9;
      break;
    case D_SHELL:
      angular_momentum_expon(&data->angular_momentum[j   ], "xx");
      angular_momentum_expon(&data->angular_momentum[j+3 ], "yy");
      angular_momentum_expon(&data->angular_momentum[j+6 ], "zz");
      angular_momentum_expon(&data->angular_momentum[j+9 ], "xy");
      angular_momentum_expon(&data->angular_momentum[j+12], "xz");
      angular_momentum_expon(&data->angular_momentum[j+15], "yz");
      j += 18;
      break;
    case F_SHELL:
      angular_momentum_expon(&data->angular_momentum[j   ], "xxx");
      angular_momentum_expon(&data->angular_momentum[j+3 ], "yyy");
      angular_momentum_expon(&data->angular_momentum[j+6 ], "zzz");
      angular_momentum_expon(&data->angular_momentum[j+9 ], "xyy");
      angular_momentum_expon(&data->angular_momentum[j+12], "xxy");
      angular_momentum_expon(&data->angular_momentum[j+15], "xxz");
      angular_momentum_expon(&data->angular_momentum[j+18], "xzz");
      angular_momentum_expon(&data->angular_momentum[j+21], "yzz");
      angular_momentum_expon(&data->angular_momentum[j+24], "yyz");
      angular_momentum_expon(&data->angular_momentum[j+27], "xyz");
      j += 30;
      break;
    case G_SHELL:
      angular_momentum_expon(&data->angular_momentum[j   ], "xxxx");
      angular_momentum_expon(&data->angular_momentum[j+3 ], "yyyy");
      angular_momentum_expon(&data->angular_momentum[j+6 ], "zzzz");
      angular_momentum_expon(&data->angular_momentum[j+9 ], "xxxy");
      angular_momentum_expon(&data->angular_momentum[j+12], "xxxz");
      angular_momentum_expon(&data->angular_momentum[j+15], "yyyx");
      angular_momentum_expon(&data->angular_momentum[j+18], "yyyz");
      angular_momentum_expon(&data->angular_momentum[j+21], "zzzx");
      angular_momentum_expon(&data->angular_momentum[j+24], "zzzy");
      angular_momentum_expon(&data->angular_momentum[j+27], "xxyy");
      angular_momentum_expon(&data->angular_momentum[j+30], "xxzz");
      angular_momentum_expon(&data->angular_momentum[j+33], "yyzz");
      angular_momentum_expon(&data->angular_momentum[j+36], "xxyz");
      angular_momentum_expon(&data->angular_momentum[j+39], "yyxz");
      angular_momentum_expon(&data->angular_momentum[j+42], "zzxy");
      j += 45;
      break;
    }
  }

  return TRUE;
}


/******************************************************
 *
 * Populate the flat arrays containing the basis
 * set data.
 *
 ******************************************************/
static int fill_basis_arrays(qmdata_t *data) {
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
  data->num_basis_funcs = primcount;

  /* reserve space for pointer to array containing basis
   * info, i.e. contraction coeficients and expansion 
   * coefficients; need 2 entries per basis function, i.e.
   * exponent and contraction coefficient; also,
   * allocate space for the array holding the orbital symmetry
   * information per primitive Gaussian.
   * Finally, initialize the arrays holding the number of 
   * shells per atom and the number of primitives per shell*/
  ALLOCATE(basis,                   float, 2*primcount);
  ALLOCATE(shell_types,             int,   data->num_shells);
  ALLOCATE(num_shells_per_atom,     int,   data->num_basis_atoms);
  ALLOCATE(num_prim_per_shell,      int,   data->num_shells);
  ALLOCATE(atomicnum_per_basisatom, int,   data->num_basis_atoms);



  /* store pointers in struct qmdata_t */
  data->basis = basis;
  data->shell_types = shell_types;
  data->num_shells_per_atom = num_shells_per_atom;
  data->num_prim_per_shell  = num_prim_per_shell;
  data->atomicnum_per_basisatom = atomicnum_per_basisatom;

  primcount = 0;
  for (i=0; i<data->num_basis_atoms; i++) {
    /* assign atomic number */
    data->basis_set[i].atomicnum = data->atoms[i].atomicnum;
    atomicnum_per_basisatom[i]   = data->atoms[i].atomicnum;

    num_shells_per_atom[i] = data->basis_set[i].numshells;

    for (j=0; j<data->basis_set[i].numshells; j++) {
      shell_types[shellcount]        = data->basis_set[i].shell[j].type;
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



/**************************************************
 *
 * Convert shell type from char to int.
 * Note that SP_P shells are assigned in get_basis()
 *
 ************************************************ */
static int shelltype_int(char *type) {
  int shelltype;
  if      (!strcmp(type, "sp")) shelltype = SP_SHELL;
  else if (!strcmp(type, "s"))  shelltype = S_SHELL;
  else if (!strcmp(type, "p"))  shelltype = P_SHELL;
  else if (!strcmp(type, "d"))  shelltype = D_SHELL;
  else if (!strcmp(type, "f"))  shelltype = F_SHELL;
  else if (!strcmp(type, "g"))  shelltype = G_SHELL;
  else shelltype = UNK_SHELL;
  
  return shelltype;
}


/***********************************************************
 *
 * Parse through the [MO] section and count orbitals and
 * the number of wavefunction coefficients per orbital.
 *
 * Format example of [MO] section:
 *  [MO]
 *  Sym=   1a         <-- this line is optional
 * Ene=   -15.5764
 *  Spin= Alpha
 * Occup=    2.00000
 *    1        -0.00000435
 *    2         0.00005919
 *    ...
 *
 ***********************************************************/
static int count_orbitals(qmdata_t *data) {
  int nr;
  int num_wave_coeff=0;
  float orbenergy, occu;
  char spin[1024];
  qm_wavefunction_t *wave;
  moldendata_t *moldendata = (moldendata_t *)data->format_specific_data;


  /* Place file pointer after [MO] keyword in line containing "Spin". */
  fseek(data->file, moldendata->filepos_mo, SEEK_SET);
  if (!goto_keyline(data->file, "Spin=", NULL)) {
    printf("moldenplugin) Couldn't find keyword 'Spin' in [MO] section!\n");
    return FALSE;
  }
 
  nr = fscanf(data->file, " Spin= %s\n", spin);
  eatline(data->file, 1);

  /* The first wavefunction should have spin alpha */
  strtoupper(spin);
  if (strcmp(spin, "ALPHA")) return FALSE;

  /* Count wavefunction coefficients */
  while (1) {
    int nr, atomid;
    char buffer[1024];
    if (!fgets(buffer,1024,data->file)) return FALSE;
    nr = sscanf(buffer,"%d %*f", &atomid);
    if (nr==0) break;
    num_wave_coeff++;
  }


  if (data->wavef_size && 
      data->wavef_size != num_wave_coeff) {
    printf("moldenplugin) No match between # wavefunction coefficients\n");
    printf("moldenplugin) and # cart. basis functions in basis set!\n");
    return FALSE;
  }

  /* Allocate memory for the qm_timestep frame */
  data->qm_timestep = (qm_timestep_t *)calloc(1, sizeof(qm_timestep_t));

  /* Add wavefunction for spin alpha */
  wave = add_wavefunction(data->qm_timestep);

  wave->spin = SPIN_ALPHA;
  wave->type = MOLFILE_WAVE_UNKNOWN;
  wave->exci = 0;
  wave->mult = 1;
  wave->num_coeffs = num_wave_coeff;

  /* Place file pointer on line after the [MO] keyword. */
  fseek(data->file, moldendata->filepos_mo, SEEK_SET);

  /* Count orbitals */
  while (1) {
    nr  = fscanf(data->file, " Ene= %f\n", &orbenergy);
    nr += fscanf(data->file, " Spin= %s\n", spin);
    nr += fscanf(data->file, " Occup= %f\n", &occu);

    eatline(data->file, num_wave_coeff);
    if (nr!=3 || toupper(spin[0])!='A') break;
    wave->num_orbitals++;
  }


  /* Add wavefunction for spin beta */
  if (!strcmp(strtoupper(spin), "BETA")) {
    wave = add_wavefunction(data->qm_timestep);
    wave->spin = SPIN_BETA;
    wave->type = MOLFILE_WAVE_UNKNOWN;
    wave->exci = 0;
    wave->mult = 1;
    wave->num_coeffs = num_wave_coeff;
    wave->num_orbitals = 1;

    while (1) {
      nr  = fscanf(data->file, " Ene= %f\n", &orbenergy);
      nr += fscanf(data->file, " Spin= %s\n", spin);
      nr += fscanf(data->file, " Occup= %f\n", &occu);

      eatline(data->file, num_wave_coeff);
      if (nr!=3 || toupper(spin[0])!='B' ||
          wave->num_orbitals>=num_wave_coeff) break;
      wave->num_orbitals++;
    }
  }

  return TRUE;
}


static int read_molecular_orbitals(qmdata_t *data) {
  moldendata_t *moldendata = (moldendata_t *)data->format_specific_data;
  qm_wavefunction_t *wave;

  if (!data->qm_timestep || moldendata->coordsonly) return FALSE;

  /* Place file pointer on line after the [MO] keyword. */
  fseek(data->file, moldendata->filepos_mo, SEEK_SET);

  wave = &data->qm_timestep->wave[0];
  ALLOCATE(wave->wave_coeffs, float, wave->num_coeffs*wave->num_orbitals);
  /* printf("num_coeffs   = %d\n", wave->num_coeffs);
  printf("num_orbitals = %d\n", wave->num_orbitals);
  printf("num_wave     = %d\n", data->qm_timestep->numwave);
  */

  /* Read wavefunction coefficients for spin alpha */
  if (!read_wave_coeffs(data->file, wave)) return FALSE;

  if (data->qm_timestep->numwave==1) return TRUE;

  /* Read wavefunction coefficients for spin beta */
  wave = &data->qm_timestep->wave[1];
  ALLOCATE(wave->wave_coeffs, float, wave->num_coeffs*wave->num_orbitals);

  if (!read_wave_coeffs(data->file, wave)) return FALSE;

  return TRUE;
}

static int read_wave_coeffs(FILE *file, qm_wavefunction_t *wave) {
  int i, j, nr;
  char buffer[1024];
  float *wave_coeffs = wave->wave_coeffs;

  for (i=0; i<wave->num_orbitals; i++) {
    eatline(file, 3);
    for (j=0; j<wave->num_coeffs; j++) {
      int atomid;
      if (!fgets(buffer,1024,file)) return FALSE;
      nr = sscanf(buffer,"%d %f", &atomid, &wave_coeffs[i*wave->num_coeffs+j]);
      /*printf("%d,%d: %d %f\n", i, j, atomid, wave_coeffs[i*wave->num_coeffs+j]);*/
      if (nr==0) {
        printf("moldenplugin) Error reading wavefunction coefficients!\n");
        return FALSE;
      }
    }
  }

  return TRUE;
}

/*************************************************************
 *
 * plugin registration 
 *
 **************************************************************/
static molfile_plugin_t plugin;

VMDPLUGIN_API int VMDPLUGIN_init() {
  memset(&plugin, 0, sizeof(molfile_plugin_t));
  plugin.abiversion = vmdplugin_ABIVERSION;
  plugin.type = MOLFILE_PLUGIN_TYPE;
  plugin.name = "molden";
  plugin.prettyname = "Molden";
  plugin.author = "Markus Dittrich, Jan Saam";
  plugin.majorv = 0;
  plugin.minorv = 5;
  plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
  plugin.filename_extension = "molden";
  plugin.open_file_read = open_molden_read;
  plugin.read_structure = read_molden_structure;

  plugin.read_timestep_metadata    = read_timestep_metadata;
  plugin.read_timestep             = read_timestep;
  plugin.read_qm_timestep_metadata = read_qm_timestep_metadata;

  plugin.read_qm_metadata = read_molden_metadata;
  plugin.read_qm_rundata  = read_molden_rundata;

  plugin.close_file_read = close_molden_read;
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_fini() {
  return VMDPLUGIN_SUCCESS;
}

