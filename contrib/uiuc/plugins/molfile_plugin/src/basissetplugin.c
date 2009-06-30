/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_basissetplugin
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
 *      $RCSfile: basissetplugin.c,v $
 *      $Author: saam $       $Locker:  $             $State: Exp $
 *      $Revision: 1.4 $       $Date: 2009/06/27 00:45:27 $
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
#define SPIN_ALPHA 0
#define SPIN_BETA  1

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

#define NUM_ELEMENTS 10



static const char *elements[] = { 
  "(unknown)", "HYDROGEN", "HELIUM", "LITHIUM", "BERYLLIUM", "BORON",
  "CARBON", "NITROGEN", "OXYGEN", "FLUORINE", "NEON"};



/* ######################################################## */
/* declaration/documentation of internal (static) functions */
/* ######################################################## */

static void print_input_data(gamessdata *);


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

static char* trimleft(char *);
static int goto_keystring(FILE *file, const char *keystring,
                                const char *stopstring);


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
static void *open_basis_read(const char *filename, 
                  const char *filetype, int *natoms) {

  FILE *fd;
  gamessdata *data;

  printf("open_basis_read();\n");
  
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

  data->num_shells = 0;
  data->num_basis_funcs = 0;
  data->num_basis_atoms = 0;

  /* initialize some of the character arrays */
  memset(data->basis_string,0,sizeof(data->basis_string));

  /* store file pointer and filename in gamess struct */
  data->file = fd;

  /* Read the basis set */
  if (!get_basis(data)) return NULL; 


  /* provide VMD with the proper number of atoms */
  *natoms = 0;

  /* Test print the parsed data in same format as logfile */
  print_input_data(data);

  return data;
}




/*****************************************************
 *
 * provide VMD with the sizes of the QM related
 * data structure arrays that need to be made
 * available
 *
 *****************************************************/
static int read_basis_metadata(void *mydata, 
    molfile_qm_metadata_t *metadata) {

  gamessdata *data = (gamessdata *)mydata;

  metadata->ncart = 0;
  metadata->nimag = 0;
  metadata->nintcoords = 0;

  metadata->have_sysinfo = 0;
  metadata->have_esp = 0;
  metadata->have_npa = 0;
  metadata->have_carthessian = 0;
  metadata->have_internals = 0;
  metadata->have_normalmodes = 0;

  /* orbital + basis set data */
  metadata->num_basis_funcs = data->num_basis_funcs;
  metadata->num_basis_atoms = data->num_basis_atoms;
  metadata->num_shells      = data->num_shells;
  metadata->wavef_size      = 0;  

  return MOLFILE_SUCCESS;
}


/******************************************************
 * 
 * Provide VMD with the static (i.e. non-trajectory)
 * data. That means we are filling the molfile_plugin
 * data structures.
 *
 ******************************************************/
static int read_basis_rundata(void *mydata, 
                               molfile_qm_t *qm_data) {

  gamessdata *data = (gamessdata *)mydata;
  int i;
  molfile_qm_basis_t   *basis_data   = &qm_data->basis;
/*   molfile_qm_sysinfo_t *sys_data     = &qm_data->run; */

/*   strncpy(sys_data->basis_string, data->basis_string, */
/*           sizeof(sys_data->basis_string)); */


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
  }
#endif
 
  return MOLFILE_SUCCESS;
}



/**********************************************************
 *
 * clean up when done and free all the memory do avoid
 * memory leaks
 *
 **********************************************************/
static void close_basis_read(void *mydata) {

  gamessdata *data = (gamessdata *)mydata;
  int i, j;
  fclose(data->file);

  free(data->initatoms);
  free(data->basis);
  free(data->shell_symmetry);
  free(data->atomicnum_per_basisatom);
  free(data->num_shells_per_atom);
  free(data->num_prim_per_shell);
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

  free(data);
}

/* ####################################################### */
/*             End of API functions                        */
/* The following functions actually do the file parsing.   */
/* ####################################################### */


#define TORF(x) (x ? "T" : "F")

static void print_input_data(gamessdata *data) {
  int i, j, k;
  int primcount=0;
  int shellcount=0;

/*   printf("\n"); */
/*   printf("     BASIS OPTIONS\n"); */
/*   printf("     -------------\n"); */
/*   printf("%s\n", data->basis_string); */
/*   printf("\n\n\n"); */
  printf("\n");
  printf("     ATOMIC BASIS SET\n");
  printf("     ----------------\n");
  printf(" THE CONTRACTED PRIMITIVE FUNCTIONS HAVE BEEN UNNORMALIZED\n");
  printf(" THE CONTRACTED BASIS FUNCTIONS ARE NOW NORMALIZED TO UNITY\n");
  printf("\n");
  printf("  SHELL TYPE  PRIMITIVE        EXPONENT          CONTRACTION COEFFICIENT(S)\n");
  printf("\n");

  printf("gamessplugin) =================================================================\n");
  for (i=0; i<data->num_basis_atoms; i++) {
    printf("%-8d (%10s)\n\n", data->basis_set[i].atomicnum, data->basis_set[i].name);
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
  printf(" TOTAL NUMBER OF ATOMS                        =%5i\n", data->numatoms);
  printf("\n");
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

  printf("get_basis()\n");

  /* initialize buffers */
  buffer[0] = '\0';
  for (i=0; i<3; i++) word[i][0] = '\0';
  
  if (!goto_keystring(data->file, "$DATA", NULL))
    printf("gamessplugin) No basis set found!\n");


  /* Allocate space for the basis for all atoms */
  /* When the molecule is symmetric the actual number atoms with
   * a basis set could be smaller */
  data->basis_set = (basis_atom_t*)calloc(1, sizeof(basis_atom_t));


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
    printf("%s",buffer);
    if (!strcmp(&word[0][0], "$END")) break;

    switch (numread) {
      case 1:
        /* Next atom */
        if (i>0) {
          data->basis_set = (basis_atom_t*)realloc(data->basis_set, (i+1)*sizeof(basis_atom_t));
        }

        strcpy(data->basis_set[i].name, &word[0][0]);


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
    int j;
    /* assign atomic number from element name */
    data->basis_set[i].atomicnum = 0;
    for (j=0; j<NUM_ELEMENTS; j++) {
      if (!strcmp(elements[j], data->basis_set[i].name))
        data->basis_set[i].atomicnum = j;
    }
    printf("%d %s\n", data->basis_set[i].atomicnum, data->basis_set[i].name);
    atomicnum_per_basisatom[i] = data->basis_set[i].atomicnum;

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
  int i, success;
  int primcounter = 0, nprim = 0;;

  GET_LINE(buffer, data->file);
  success = sscanf(buffer,"%c %d", shellsymm, &nprim);

  (*prim) = (prim_t*)calloc(nprim, sizeof(prim_t));

  for (i=0; i<nprim; i++) {
    GET_LINE(buffer, data->file);
    success = sscanf(buffer,"%*d %f %f %f",
                       &exponent, &contract[0], &contract[1]); 
printf("%s",buffer);

    /* store in basis array and increase the counter */ 
    switch (success) {
      case 2:
        /* store exponent */
        (*prim)[i].exponent = exponent;
          
        /* store coefficient */
        (*prim)[i].contraction_coeff = contract[0];

        primcounter++;
        break;

      case 3:
        /* store exponent */
        (*prim)[i].exponent = exponent;
          
        /* store coefficient */
        (*prim)[i].contraction_coeff = contract[icoeff];
        
        primcounter++;
        break;

      case -1:
        /* otherwise it's an empty line which represents the end of the shell */
        break;

      case 1:
        /* the user had given the next atom a numeric name */
        break;
    }

  }

  if (!primcounter) free(*prim);

  return primcounter;
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


/* Advances the file pointer until the first appearance
 * of a line beginning with the given keystring. Leading
 * whitespace in the lines are ignored, the keystring 
 * should not begin with whitespace otherwise no match
 * will be found.
 * If stopstring is encountered before the keystring 
 * the file is rewound to the position where the search
 * started. If stopstring is NULL then the search stops
 * at EOF. */
static int goto_keystring(FILE *file, const char *keystring,
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
  plugin.name = "basisset";
  plugin.prettyname = "Basis Set";
  plugin.author = "Jan Saam";
  plugin.majorv = 0;
  plugin.minorv = 1;
  plugin.is_reentrant = VMDPLUGIN_THREADUNSAFE;
  plugin.filename_extension = "basis";
  plugin.open_file_read  = open_basis_read;
  plugin.close_file_read = close_basis_read;
  plugin.read_structure = NULL;

  plugin.read_qm_metadata = read_basis_metadata;
  plugin.read_qm_rundata  = read_basis_rundata;

#if vmdplugin_ABIVERSION > 11
  plugin.read_timestep_metadata    = NULL;
  plugin.read_qm_timestep_metadata = NULL;
  plugin.read_timestep = NULL;
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
