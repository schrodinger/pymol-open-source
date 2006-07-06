/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_gamessplugin
#define STATIC_PLUGIN 1

/***************************************************************************
 *cr
 *cr            (C) Copyright 1995-2006 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
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

#include "gamessplugin.h"

/* 
 * pre-processor macro to activate grid code 
 */
#ifdef GRID_ACTIVE
#define GRID 1
#else
#define GRID 0
#endif

/*
 * Error reporting macro for use in DEBUG mode
 */
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
#define ERR_FALSE(x) if ( x == NULL ) return FALSE;



/***************************************************************
 *
 * subroutine doing the initial reading of the GAMESS
 * input file
 *
 * *************************************************************/
static void *open_gamess_read(const char *filename, 
                  const char *filetype, int *natoms) {

  FILE *fd;
  gamessdata *data;
  char mytime[BUFSIZ];

  /* initialize array */
  mytime[0] = '\0';

  
  /* open the input file */
  fd = fopen(filename, "rb");
 
  /* make sure we were able to load the file properly*/
  if (!fd) 
  {
    PRINTERR;
    return NULL;
  }

  /* obtain the current time */
  get_time(mytime);

  /* print out a short warning message -- we're very alpha */
  printf("\n");
  printf("     *************************************************\n");
  printf("     ***       GAMESSPLUGIN for VMD                ***\n");
  printf("     *** %s      ***\n", mytime);
  printf("     *** email bugs to vmd@ks.uiuc.edu             ***\n");
  printf("     *** v0.3.0     04/03/2006 (c) Markus Dittrich ***\n");
  printf("     *************************************************\n");
  printf("\n");

  /* allocate memory for main data structure */
  data = (gamessdata *)calloc(1,sizeof(gamessdata));

  /* make sure memory was allocated properly */
  if (data == NULL) 
  {
    PRINTERR;
    return NULL;
  }


  /**********************************
   * initialize some variables
   *********************************/

  /* volumetric */
  data->have_volumetric = 0; 

  /* runtyp */
  data->runtyp = 0;

  /* have_trajectory flag to 1 */
  data->have_trajectory = 1; 

  /* reset trajectory counter */
  data->num_traj_points = 0;

  /* presence of wavefunction */
  data->got_wavefunction = FALSE;

  /* orbital index of HOMO */
  data->homo_index = 0;

  /* flag indicating the presence of
   * internal coordinates */
  data->have_internals = FALSE;

  /* flag indicating the presence of
   * the cartesian Hessian */
  data->have_cart_hessian = FALSE;

  /* initialize the number of imaginary 
   * modes for HESSIAN type runs */
  data->nimag = 0;

  /* initialize the GAMESS version */
  data->version = 0;

  /* initialize the basis set info */
  data->ngauss = 0;
  data->npfunc = 0;
  data->ndfunc = 0;
  data->nffunc = 0;
  data->diffs  = 0;
  data->diffsp = 0;

  /* initialize the number of scfenergies */
  data->num_scfenergies = 0;

  /* initialize some of the character arrays */
  strncpy(data->runtyp_string,"\0",sizeof(data->runtyp_string));
  strncpy(data->basis_string,"\0",sizeof(data->basis_string));
  strncpy(data->version_string,"\0",sizeof(data->version_string));
  strncpy(data->scftyp_string,"\0",sizeof(data->scftyp_string));
  strncpy(data->memory,"\0",sizeof(data->memory));

  
  /***********************************
   * done intializing
   **********************************/


  /* store file pointer and filename in gamess
   * struct */
  data->file = fd;
  data->file_name = strdup(filename);


  /* check if the file is GAMESS format; if yes
   * parse it, if no exit */
  if ( have_gamess(data) == TRUE ) 
  {
    /* if we're dealing with an unsupported GAMESS
     * version, we better quit */
    if ( data->version == 0 )
    {
      printf("gamessplugin> GAMESS version %s not supported. \n",
          data->version_string);
      printf("gamessplugin> .... bombing out! Sorry :( \n");
      return NULL;
    }

    /* get the "static" information from the log file */    
    if ( parse_gamess_log_static(data,natoms) == FALSE ) 
      return NULL;
  }
  else 
  {
    return NULL;
  }


  /* done with the gamess plugin */
  rewind(fd);
  printf("\n");


  /* copy pointer also in bogus tcl data pointer */
  /* tcl_pointer = data; */

  return data;
}


/************************************************************
 * 
 * subroutine reading in the structure of the molecule
 * in the GAMESS output file
 *
 *************************************************************/
static int read_gamess_structure(void *mydata, int *optflags, 
                      molfile_atom_t *atoms) 
{
  gamessdata *data = (gamessdata *)mydata;
  gamess_temp *temp_atom;
  molfile_atom_t *atom;
  unsigned int i = 0;
 
  *optflags = MOLFILE_NOOPTIONS; /* no optional data */

  /* all the information I need has already been read in
   * via the initial scan and I simply need to copy 
   * everything from the temporary arrays into the 
   * proper VMD arrays.
   * Since there are no atom names in the GAMESS output
   * I use the atom type here --- maybe there is a better
   * way to do this !!?? */

  /* get initial pointer for temp arrays */
  temp_atom = data->temporary;

  for(i=0;i<data->numatoms;i++)
  {
    atom=atoms+i;
    strncpy(atom->name,temp_atom->type,sizeof(atom->name)); 
    strncpy(atom->type,temp_atom->type,sizeof(atom->type));
    strncpy(atom->resname,"",sizeof(atom->resname)); 
    atom->resid = 1;
    atom->charge = temp_atom->charge;
    atom->chain[0] = '\0';
    atom->segid[0] = '\0';
   

    temp_atom++;
  }
 
  return MOLFILE_SUCCESS; 
}


/***********************************************************
 *
 * this function reads in the information for the next
 * timestep
 *
 ***********************************************************/
static int read_next_timestep(void *mydata, int natoms, 
    molfile_timestep_t *ts) 
{
  gamessdata *data = (gamessdata *)mydata;
  gamess_temp *temp_atom;
  unsigned int i = 0;

#ifdef ANIMATE_MODE
  mode_data *animated_mode = data->animated_mode;
#endif
  
  /* Now, if we are not dealing with a trajectory (a single 
   * * point conformation is considered one element of a trajectory
   */
  if (data->have_trajectory == 0) return MOLFILE_ERROR;
 

  /* in the case of RUNTYP=ENERGY read data from temp
   * arrays and that's it;
   * in the case of RUNTYP=OPTIMIZE and SADPOINT just start with
   * the first 1NSERCH entry */
  if ( data->runtyp == ENERGY ) 
  {
     /* initialize pointer for temporary arrays */
     temp_atom = data->temporary; 
 

     /* now copy the initial coordinates from the temporary
      * arrays into the VMD ones */
     for(i=0;i<natoms;i++)  {

       ts->coords[3*i  ] = temp_atom->x;
       ts->coords[3*i+1] = temp_atom->y;
       ts->coords[3*i+2] = temp_atom->z; 

       temp_atom++;
     }    

     /* Since we are only dealing with a single point which has
      * just been read we have to make sure to let VMD know 
      * and also register the number of trajectory points as 1 */
     data->have_trajectory = 0;

     return MOLFILE_SUCCESS;
  }

#ifndef ANIMATE_MODE
  else if ( data->runtyp == HESSIAN ) 
  {
     /* initialize pointer for temporary arrays */
     temp_atom = data->temporary; 
 

     /* now copy the initial coordinates from the temporary
      * arrays into the VMD ones */
     for(i=0;i<natoms;i++)  {

       ts->coords[3*i  ] = temp_atom->x;
       ts->coords[3*i+1] = temp_atom->y;
       ts->coords[3*i+2] = temp_atom->z; 

       temp_atom++;
     }    

     /* Since we are only dealing with a single point which has
      * just been read we have to make sure to let VMD know 
      * and also register the number of trajectory points as 1 */
     data->have_trajectory = 0;

     return MOLFILE_SUCCESS;
  }
#endif

#ifdef ANIMATE_MODE
  else if ( data->runtyp == HESSIAN ) 
  {
    /* read until we run out of frames for the current mode */
    if ( animated_mode->current_mode_frame < 
	      4*animated_mode->mode_num_frames+3 )
    {
      /* now copy the initial coordinates from the temporary
       	* arrays into the VMD ones */
      for(i=0; i< natoms; i++)  
      {
	ts->coords[3*i  ] = *(animated_mode->mode_frames +
	    (animated_mode->current_mode_frame*3*natoms)+3*i);

	ts->coords[3*i+1] = *(animated_mode->mode_frames +
	    (animated_mode->current_mode_frame*3*natoms)+3*i+1);

	ts->coords[3*i+2] = *(animated_mode->mode_frames +
	    (animated_mode->current_mode_frame*3*natoms)+3*i+2);
      }    

      /* increase the trajectory and mode counter */ 
      (animated_mode->current_mode_frame)++;
    }
    else 
    {
      return MOLFILE_ERROR;
    }

    return MOLFILE_SUCCESS;
  }
#endif

  else if ( data->runtyp == OPTIMIZE || data->runtyp == SADPOINT ) 
  {
    /* call the routine scanning the output file for the
     * trajectory */
    if ( get_trajectory(data,ts,natoms) == FALSE) {

      /* before we return let's check if we have read
       * any trajectory points whatsoever; if not
       * this could mean that the file might be
       * truncated, and in this case we dump the initial
       * coordinate info, such that the user has something
       * to look at; */
      if ( data->num_traj_points == 0 ) {

	/* initialize pointer for temporary arrays */
	temp_atom = data->temporary; 
 
	
	/* now copy the initial coordinates from the temporary
	* arrays into the VMD ones */
	for(i=0;i<natoms;i++)  {

	  ts->coords[3*i  ] = temp_atom->x;
	  ts->coords[3*i+1] = temp_atom->y;
	  ts->coords[3*i+2] = temp_atom->z; 

	  temp_atom++; }


	/* make sure we don't continue reading the
	 * trajectory */
        data->have_trajectory = 0;


	/* that's all we can do */
	return MOLFILE_SUCCESS;
      }


      return MOLFILE_ERROR;
    }


    /* increase the trajectory counter */
    (data->num_traj_points)++;
   

    /* success, it seems :) */
    return MOLFILE_SUCCESS;
  }


  /* all other cases should have been rejected allready
   * in open_gamess_read;
   * but just to make sure we return MOLFILE_ERROR 
   * by default */
  return MOLFILE_ERROR;
}


/*********************************************************
 *
 * this subroutine makes the volumetric orbital metadata
 * available to VMD;
 *
 ********************************************************/
static int read_orbital_metadata( void *mydata, int *nsets, 
                          molfile_volumetric_t **metadata) 
{
    gamessdata *data = (gamessdata *)mydata;

    /* first of all we have to test the presence of
     * volumetric data; if not set number of volumetric
     * data sets to zero and return */

    if ( data->have_volumetric != 1 ) 
    {
      *nsets = 0;
      return MOLFILE_SUCCESS;
    }

    /* we only have one set of volumetric data */
    *nsets = 1; 

    /* all the necessary data structures have been filled
     * in the function get_system_dimensions and we only
     * need to copy the pointer. */
    *metadata = data->vol;

    return MOLFILE_SUCCESS;
}


/**********************************************************
 *
 * this subroutine makes the volumetric data available 
 * to VMD
 *
 *********************************************************/
static int read_orbital_data( void *mydata, int set, 
    float *datablock, float *colorblock) 
{
  gamessdata *data = (gamessdata *)mydata;
  gamess_temp *temp_data;
  molfile_volumetric_t *vol_metadata;
  unsigned int i = 0;
  int numxyz;

  /* get pointer of temporary data storage */
  temp_data = data->temporary; 

  /* get pointer to volumetric metadata */
  vol_metadata = data->vol;

  /* move array containing volumetric over to datablock 
   * the array contains xsize*ysize*zsize elements */
  numxyz = ( vol_metadata->xsize * vol_metadata->ysize * 
             vol_metadata->zsize );

  for ( i = 0 ; i < numxyz ; ++i )
  {
    *(datablock + i) = *(data->orbital_grid+i);
  }

  return MOLFILE_SUCCESS;
}


/**********************************************************
 *
 * clean up when done and free all the memory do avoid
 * memory leaks
 *
 **********************************************************/
static void close_gamess_read(void *mydata) {

  gamessdata *data = (gamessdata *)mydata;
  fclose(data->file);

  /* free memory */
  free(data->temporary);
  free(data->basis);
  free(data->system_dimensions);
  free(data->system_center); 
  free(data->orbital_grid);
  free(data->basis_counter);
  free(data->wave_function);
  free(data->orbital_symmetry);
  free(data->atomic_shells);
  free(data->shell_primitives);
  free(data->orbital_energy);
  free(data->atomic_number);
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
  free(data->vol);
  free(data->scfenergies);
  free(data->wavenumbers);
  free(data->intensities);
  free(data->normal_modes);
  
  if ( data->runtyp == HESSIAN )
  {
#ifdef ANIMATE_MODE
    free(data->animated_mode->mode_frames);
    free(data->animated_mode);
#endif
  }

  free(data);
}


/*************************************************************
 *
 * registration stuff 
 *
 **************************************************************/
static molfile_plugin_t gamessplugin = {
  vmdplugin_ABIVERSION,
  MOLFILE_PLUGIN_TYPE,                         /* type */
  "gamess",                                    /* short name */
  "GAMESS",                                    /* pretty name */
  "Markus Dittrich",                           /* author */
  0,                                           /* major version */
  3,                                           /* minor version */
  VMDPLUGIN_THREADSAFE,                        /* is reentrant */
  "log",
  open_gamess_read,
  read_gamess_structure,
  0,
  read_next_timestep,
  close_gamess_read,
  0,
  0,
  0,
  0,
  read_orbital_metadata,
  read_orbital_data
};


VMDPLUGIN_API int VMDPLUGIN_init(void) {
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)&gamessplugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_fini(void) {
  return VMDPLUGIN_SUCCESS;
}




/********************************************************
 *
 * THE NEXT SECTION CONTAINs THE GAMESS OUTPUT FILE
 * SPECIFIC SUBROUTINES
 *
 ********************************************************/





/********************************************************
 *
 * this routine is the main gamess log file
 * parser responsible for static, i.e. 
 * non-trajectory information 
 *
 ********************************************************/
int parse_gamess_log_static(void *mydata, int *natoms) 
{
  gamessdata *data = (gamessdata *)mydata;

  /* determine the number of procs and the amount
   * of memory requested */
  if ( get_proc_mem(data) == FALSE ) return FALSE;


  /* determine the GBASIS used for the run */
  if ( get_gbasis(data) == FALSE ) return FALSE;


  /* read the run title */
  if ( get_runtitle(data) == FALSE ) return FALSE;


  /* read the contrl group; here we determine the
   * RUNTYP and bomb out if we don't support
   * it; also read in some controlflow related
   * info as well as the COORD type */
  if ( check_contrl(data) == FALSE ) return FALSE;


  /* now call routine get_initial_info to obtain 
   * the number of atoms, the atom types, and the
   * coordinates, which will be stored in gamess_temp */
  if ( get_initial_info(data) == FALSE ) return FALSE;


  /* provide VMD with the proper number of atoms */
  *natoms = data->numatoms;


  /* read in the number of orbitals, charge, 
   * multiplicity, ... */
  if ( get_num_orbitals(data) == FALSE ) return FALSE;


  /* check if the wavefunction and orbital stuff
   * is suported for current GBASIS; if not stop
   * here and return data pointer */
  if ( have_supported_gbasis(data) == FALSE) return TRUE;


  /* read in the guess options */
  if ( get_guess(data) == FALSE ) return FALSE;

  
  /* read in the final wavefunction; do this only for 
   * RUNTYP=ENERGY for now. For the other runs I have to
   * put in some more infrastructure to make sure to
   * correlate each of the possible many geometries
   * during e.g. an optimization run with the correstpoding
   * wavefunction 
   * for now we also have to restrict reading the wavefunctions
   * for singlets only, since otherwise the punched wavefunction
   * output differs;
   * At least the call to generate orbital grid should
   * not be done automagically but rather after a TCL call
   * only */
  if ( data->runtyp == ENERGY 
	 && data->num_orbitals_A == data->num_orbitals_B 
	 && GRID)
  {
    /* read in the basis set information */
    if ( get_basis(data) == FALSE ) return FALSE; 


    /* read the wavefunction */
    if ( get_wavefunction(data) == FALSE) return FALSE;  


    /* figure out which one of the orbitals is the HOMO */
    if ( find_homo(data) == FALSE) return FALSE; 

 
    /* generate the grid containing the volumetric data
     * for an orbital */
    if( orbital_grid_driver(data) == FALSE ) return FALSE;  

  }

  return TRUE;
}



/******************************************************
 *
 * this function extracts the trajectory information
 * from the output file
 *
 * *****************************************************/
int get_trajectory(void *mydata, molfile_timestep_t *ts,
                   int natoms) 
{
  gamessdata *data = (gamessdata *)mydata;
  char buffer[BUFSIZ];
  char word[4][BUFSIZ];
  char type[8];
  unsigned int i = 0;
  float charge,x,y,z;
  double scfenergy = 0.0;

  /* zero out the word arrays */
  buffer[0] = '\0';
  type[0] = '\0';
  for ( i = 0; i < 4; ++i) word[i][0] = '\0';


  /* search for the first coordinate frame during optimization */
  ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )


  /* at this point we have to distinguish between
   * pre="27 JUN 2005 (R2)" and "27 JUN 2005 (R2)"
   * versions since the output format for geometry
   * optimizations has changed */
  if ( data->version == 1)
  {
     /* scan for the next optimized geometry punched in the 
      * logfile */
     do 
     {
       sscanf(buffer,"%s %s",&word[0][0],&word[1][0]); ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )

     } while (strcmp(&word[0][0],"1NSERCH="));
  }
  else if ( data->version == 2 )
  {
     /* scan for the next optimized geometry punched in the 
      * logfile */
     do 
     {
       sscanf(buffer,"%s %s %s",&word[0][0],&word[1][0],
	    &word[2][0]); 
       ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )

     } while (strcmp(&word[0][0],"BEGINNING") || 
	      strcmp(&word[1][0],"GEOMETRY")  ||
	      strcmp(&word[2][0],"SEARCH"));
  }


  /* skip the next three lines */
  eatline(data->file);
  eatline(data->file);
  eatline(data->file);
  

  /* read the actual coordinates */
  for(i=0;i<natoms;i++) 
  {
    ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )
    sscanf(buffer,"%s %f %f %f %f",type,&charge,&x,&y,&z);  

    /* store the coordinates */
    ts->coords[3*i  ] = x; 
    ts->coords[3*i+1] = y;
    ts->coords[3*i+2] = z; 
  }     


  /* now we look for the SCF energy of this trajectory 
   * point */
   do 
   {
    ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )
    sscanf(buffer,"%s %s %s %lf",&word[0][0],&word[1][0],
	&word[2][0],&scfenergy); 
 
    } while (strcmp(&word[2][0],"ENERGY=") && 
	     strcmp(&word[1][0],"ENERGY=")); 
  
  /* allocate more memory for scfenergy array */
  data->scfenergies = (double *)realloc(data->scfenergies,
	(data->num_traj_points+1)*sizeof(double));


  /* append current scfenergy to array and increase counter */
  *(data->scfenergies+data->num_traj_points) = scfenergy;
  data->num_scfenergies++;


  /* done with this trajectory point */ 
  return TRUE;
}



/******************************************************
 *
 * this function performs an initial file read
 * to check its validity and determine the number
 * of atoms in the system 
 *
 * *****************************************************/
int get_initial_info(void *mydata)
{
  gamessdata *data = (gamessdata *)mydata;
  gamess_temp *atoms;
  char buffer[BUFSIZ];
  char word[4][BUFSIZ];
  char atname[BUFSIZ];
  float charge;
  float x,y,z;
  unsigned int numatoms;
  unsigned int bohr;
  unsigned int i,n;
  unsigned int have_normal_modes = TRUE;
  double scfenergy;
  char *status;

  /* initialize buffers */
  buffer[0] = '\0';
  atname[0] = '\0';
  for ( i = 0; i < 4; ++i ) word[i][0] = '\0';
 

  /* alocate temporary data structure */
  atoms = (gamess_temp *)calloc(MAXQMATOMS,sizeof(gamess_temp));

  /* make sure memory was allocated properly */
  if ( atoms == NULL) 
  {
    PRINTERR; 
    return FALSE;
  }


  /* save pointer */
  data->temporary = atoms;


  /* counts the number of atoms in input file */
  numatoms=0;  

  /* look for the initial coordinate section 
   * NOTE: I also have to check if the coordinates are
   * in Angstrom or Bohr and convert them to A in the
   * latter case */
  rewind(data->file);

  do 
  {
    ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )
    sscanf(buffer,"%s %s %s %s",&word[0][0],&word[1][0],
	&word[2][0],&word[3][0]);

  } while(strcmp(&word[0][0],"ATOM") || 
          strcmp(&word[1][0],"ATOMIC"));

  
  /* test if coordinate units are Bohr */
  bohr = (!strcmp(&word[3][0],"(BOHR)"));


  /* skip next line */
  eatline(data->file);


  /* now read in the coordinates until an empty line is being
   * encoutered */
  ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )
  n = sscanf(buffer,"%s %f %f %f %f",atname,&charge,&x,&y,&z);

  /* we expect 5 entries per line */
  while( n == 5 ) 
  {
    strncpy(atoms->type,atname,sizeof(atoms->type));
    atoms->charge = charge;


    /* if coordinates are in Bohr convert them to Angstrom here */
    if (bohr) 
    {
      x = BOHR_TO_ANGS * x;
      y = BOHR_TO_ANGS * y;
      z = BOHR_TO_ANGS * z;
    }

    atoms->x = x;
    atoms->y = y;
    atoms->z = z; 


    atoms++;
    numatoms++;

    ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )
    n = sscanf(buffer,"%s %f %f %f %f",atname,&charge,&x,&y,&z);
  }

  /* store number of atoms in data structure */
  data->numatoms = numatoms;



  /* save a copy of atom charges in atomic_number array
   * of gamessdata */

  /* alocate temporary data structure */
  data->atomic_number = (int *)calloc(numatoms,sizeof(int));

  /* make sure memory was allocated properly */
  if (data->atomic_number == NULL) 
  {
    PRINTERR; 
    return FALSE;
  }

  
  /* move data over */ 
  for ( n = 0; n < numatoms; n++) 
  {
    *(data->atomic_number + n) = (int)(data->temporary + n)->charge;
  }


  /* at this point we create an empty scfenergy arrar */
  data->scfenergies = (double *)calloc(0,sizeof(double));


  /* of we are dealing with a single point energy run 
   * we read in the final SCF energy and store it;
   * TODO: for now also include HESSIAN runs here even though
   * this is not completely correct in the case of numerical
   * runs*/
  if ( (data->runtyp == ENERGY) || (data->runtyp == HESSIAN) )
  {
    status = fgets(buffer,sizeof(buffer),data->file);
     
    while ( strcmp(&word[0][0],"FINAL") || 
	    strcmp(&word[2][0],"ENERGY") ||
	    strcmp(&word[3][0],"IS"))
    {
      sscanf(buffer,"%s %s %s %s %lf",&word[0][0],&word[1][0],
	     &word[2][0],&word[3][0],&scfenergy);

      if ( fgets(buffer,sizeof(buffer),data->file) == NULL) break;
    }

    data->scfenergies = (double *)realloc(data->scfenergies,
             sizeof(double));
    *(data->scfenergies) = scfenergy;
    data->num_scfenergies++;
  }


  /* next we check for the Mulliken charges; if they are present
   * we set the have_mulliken charge flag, but don't bomb out */

  /* set have_mulliken initially to true */
  data->have_mulliken = TRUE;

  do 
  {
    status = fgets(buffer,sizeof(buffer),data->file);

    /* check for EOF */
    if (status == NULL) {
     
      data->have_mulliken = FALSE;
      break;
    }

    sscanf(buffer,"%s %s %s %s",&word[0][0],&word[1][0],
	&word[2][0],&word[3][0]);

  } while(strcmp(&word[0][0],"TOTAL") || 
          strcmp(&word[1][0],"MULLIKEN") ||
          strcmp(&word[2][0],"AND")   || 
	  strcmp(&word[3][0],"LOWDIN") );


  /* if Mulliken Charges are present read them, if not, rewind
   * the file and look for the ESP charges */
  if ( data->have_mulliken == TRUE ) 
  {
    /* reserve memory for dynamically allocated
     * Mulliken charge array */
    data->mulliken_charges = 
      (double *)calloc(numatoms,sizeof(double));

    /* check for success */
    if (data->mulliken_charges == NULL) 
    {
      PRINTERR; 
      return FALSE;
    }


    /* skip next line */
    eatline(data->file);

    for ( n = 0; n < numatoms; ++n) 
    {
      /* make sure we don't get nuked by a bad file */
      ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )

      sscanf(buffer,"%s %s %s %lf ",&word[0][0], &word[1][0], 
	  &word[2][0], (data->mulliken_charges+n));
    }
  }

  /* no Mulliken charges found; in this case we don't know where
   * we are in the file; hence rewind file */
  else 
  {
    rewind(data->file); 
  }


  /* next we check for the ESP charges; if they are present
   * we set the have_esp charge flag, but don't bomb out */

  /* set have_esd initially to true */
  data->have_esp = TRUE;

  do 
  {
    status = fgets(buffer,sizeof(buffer),data->file);

    /* check for EOF */
    if (status == NULL) {
     
      data->have_esp = FALSE;
      break;
    }

    sscanf(buffer,"%s %s %s",&word[0][0],&word[1][0],&word[2][0]);

  } while(strcmp(&word[0][0],"ATOM") || 
          strcmp(&word[1][0],"CHARGE") ||
          strcmp(&word[2][0],"E.S.D.") );


  /* if ESP Charges are present read them, if not, rewind
   * the file and return */
  if ( data->have_esp == TRUE ) 
  {
    /* reserve memory for dynamically allocated
     * Mulliken charge array */
    data->esp_charges = 
      (double *)calloc(numatoms,sizeof(double));

    /* check for success */
    if (data->esp_charges == NULL) 
    {
      PRINTERR; 
      return FALSE;
    }


    /* skip next line */
    eatline(data->file);

    for ( n = 0; n < numatoms; ++n) {

      /* make sure we don't get nuked by a bad file */
      ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )
      sscanf(buffer,"%s %lf ",&word[0][0], (data->esp_charges+n));
    }
  }

  /* no ESP charges found; in this case we don't know where
   * we are in the file; hence rewind file */
  else 
  {
    rewind(data->file); 
  }


  /* short message */ 
  printf("gamessplugin> Detected %d atoms in the input file. \n",
          data->numatoms);


  if ( data->runtyp == HESSIAN ) 
  {

#ifdef ANIMATE_MODE
    initialize_animated_mode(data);
#endif

    /* try reading the hessian matrix in internal and
     * cartesian coordinates as well as the internal
     * coordinates together with their associated
     * force constants */
    if ( get_int_coords(data) == TRUE ) 
    {
      data->have_internals = TRUE;
    }
    else 
    {
      printf("gamessplugin> \n");
      printf("gamessplugin> Could not determine the internal \n");
      printf("gamessplugin> coordinate and Hessian info!! \n");
      printf("gamessplugin> \n");
    }
    
 
    if ( get_cart_hessian(data) == TRUE ) 
    {
      data->have_cart_hessian = TRUE;
    }
    else 
    {
      printf("gamessplugin> \n");
      printf("gamessplugin> Could not determine the cartesian \n");
      printf("gamessplugin> Hessian matrix!! \n");
      printf("gamessplugin> \n");
    }


    /* read the wavenumbers, intensities of the normal modes 
     * as well as the modes themselves */
    if ( get_normal_modes(data) == FALSE ) 
    {
      printf("gamessplugin> \n");
      printf("gamessplugin> Could not scan the normal modes,\n");
      printf("gamessplugin> \n");

      have_normal_modes = FALSE;
    }

#ifdef ANIMATE_MODE
    /* generate an animated trajectory for mode i */ 
    if (data->have_cart_hessian && have_normal_modes) 
    {
      if ( animate_normal_mode(data, 62) == FALSE )
      {
       	printf("gamessplugin> \n");
	printf("gamessplugin> Error generating animated normal mode");
	printf("gamessplugin> \n \n");
      }
    }
#endif

  }

  /* done with this one */
  return TRUE; 
}



/******************************************************
 *
 * this function generates animated frames for normal
 * mode mode_to_animate
 *
 * *****************************************************/
int initialize_animated_mode(void *mydata)
{
  gamessdata *data = (gamessdata *)mydata;
  mode_data *animated_mode;

  /* allocate memory for a mode_data struct */
  animated_mode = (mode_data *)calloc(1,sizeof(mode_data)); 

  if ( animated_mode == NULL )
  {
    PRINTERR;
    return FALSE;
  }

  /* save pointer in gamessdata struct */
  data->animated_mode = animated_mode;

  /* set the number of frames per animated mode */
  animated_mode->mode_num_frames = 25;

  /* initialize the animated mode tracker */
  animated_mode->current_mode_frame = 0;

  /* set the scaling factor for the animated mode */
  animated_mode->mode_scaling = 0.01;

  /* reserve memory for animated mode trajectory */
  animated_mode->mode_frames = (double *)calloc((4 *
      animated_mode->mode_num_frames+3)*data->numatoms*3, 
      sizeof(double));

  if ( animated_mode->mode_frames == NULL )
  {
    PRINTERR;
    return FALSE;
  }

  return TRUE;
}



/************************************************************
 *
 * this function animates a given normal mode by means of
 * generating mod_num_frames frames away from the equilibrium
 * structure in a direction given by the hessiane 
 *
 ************************************************************/
int animate_normal_mode(void *mydata, unsigned int mode)
{
  gamessdata *data = (gamessdata *)mydata;
  mode_data *animated_mode = data->animated_mode;
  double *normal_modes = data->normal_modes;
  double scale = animated_mode->mode_scaling;
  unsigned int i = 0, k = 0; 
  int l = 0, m = 0;
  unsigned int natoms = data->numatoms;
  unsigned int num_frames = animated_mode->mode_num_frames;

  /* first sweep to max of interval */
  for ( k = 0; k < num_frames+1; ++k)
  {
    for ( i = 0; i < natoms; ++i)
    {
      *(animated_mode->mode_frames+(k*natoms*3)+(3*i)) = 
	  (data->temporary+i)->x * (1+( k*scale * 
	  (*(normal_modes+(mode*natoms*3)+(3*i)))));

      *(animated_mode->mode_frames+(k*natoms*3)+(3*i+1)) = 
	  (data->temporary+i)->y * (1+( k*scale * 
	  (*(normal_modes+(mode*natoms*3)+(3*i+1)))));

      *(animated_mode->mode_frames+(k*natoms*3)+(3*i+2)) = 
	  (data->temporary+i)->z * (1+( k*scale * 
	  (*(normal_modes+(mode*natoms*3)+(3*i+2)))));
    }
  }


  /* second sweep all the way back to min of interval */
  for ( l = 0; l < 2*num_frames+1; ++l)
  {
    for ( i = 0; i < natoms; ++i)
    {
      *(animated_mode->mode_frames+((l+k)*natoms*3)+(3*i)) = 
	  (data->temporary+i)->x * (1+((int)(num_frames-l)*scale * 
	  (*(normal_modes+(mode*natoms*3)+(3*i)))));

      *(animated_mode->mode_frames+((l+k)*natoms*3)+(3*i+1)) = 
	  (data->temporary+i)->y * (1+((int)(num_frames-l)*scale * 
	  (*(normal_modes+(mode*natoms*3)+(3*i+1)))));

      *(animated_mode->mode_frames+((l+k)*natoms*3)+(3*i+2)) = 
	  (data->temporary+i)->z * (1+((int)(num_frames-l)*scale * 
	  (*(normal_modes+(mode*natoms*3)+(3*i+2)))));
    }
  }


  /* third sweep back to the native starting structure */
  for ( m = 0; m < num_frames+1; ++m)
  {
    for ( i = 0; i < natoms; ++i)
    {
      *(animated_mode->mode_frames+((l+k+m)*natoms*3)+(3*i)) = 
	  (data->temporary+i)->x * (1+((int)(m-num_frames)*scale * 
	  (*(normal_modes+(mode*natoms*3)+(3*i)))));

      *(animated_mode->mode_frames+((l+k+m)*natoms*3)+(3*i+1)) = 
	  (data->temporary+i)->y * (1+((int)(m-num_frames)*scale * 
	  (*(normal_modes+(mode*natoms*3)+(3*i+1)))));

      *(animated_mode->mode_frames+((l+k+m)*natoms*3)+(3*i+2)) = 
	  (data->temporary+i)->z * (1+((int)(m-num_frames)*scale * 
	  (*(normal_modes+(mode*natoms*3)+(3*i+2)))));
    }
  }

  /* short message */
  printf("gamessplugin> Successfully animated mode %d \n", mode);

  return TRUE;
}



/***********************************************************
 *
 * this function reads in the wavenumbers and intensities of
 * the normal modes
 *
 * *********************************************************/
int get_normal_modes(void *mydata)
{
  gamessdata *data = (gamessdata *)mydata;
  char word[4][BUFSIZ];
  char buffer[BUFSIZ];
  unsigned int i = 0, k = 0, j = 0;
  unsigned int remaining_columns;
  double entry[6]; 
  char separator;
  char *item;
  unsigned int counter;


  /* initialize arrays */
  buffer[0] = '\0';
  memset(entry, 0, sizeof(entry));
  for ( i = 0; i < 4; ++i) word[i][0] = '\0';

    
  /* reserve memory for dynamically allocated data
   * arrays */
  item = (char *)calloc(BUFSIZ,sizeof(char));

  if ( item == NULL )
  {
    PRINTERR;
    return FALSE;
  }


  data->wavenumbers = 
    (double *)calloc(data->numatoms*3,sizeof(double));

  if ( data->wavenumbers == NULL ) 
  {
    PRINTERR;
    return FALSE;
  }


  data->intensities = 
    (double *)calloc(data->numatoms*3,sizeof(double));

  if ( data->intensities == NULL ) 
  {
    PRINTERR;
    return FALSE;
  }


  data->normal_modes = 
    (double *)calloc((data->numatoms*3)*(data->numatoms*3),
		     sizeof(double));

  if ( data->normal_modes == NULL ) 
  {
    PRINTERR;
    return FALSE;
  }


  /* look for FREQUENCYs and IR INTENSITIES */
  for ( i = 0; i < (unsigned int)(data->numatoms*3/5); ++i) 
  {
    do 
    {
      ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )
      sscanf(buffer,"%s ",&word[0][0]);

    } while(strcmp(&word[0][0],"FREQUENCY:"));


    /* scan the frequencies; 
     * this requires some care since there might
     * be imaginary modes present which leads to
     * an additional char per imaginary mode per
     * line */

    /* check for imaginary modes */
    if ( strchr(buffer,'I') != NULL ) 
    {
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
        if ( *item == 'I' ) 
	{
	  data->nimag++;
	}
	else 
	{
	  /* save only the first 5 modes - there NEVER should
	   * be more in any case, but just to make sure
	   * we don't overrun the array */
	  if ( counter < 5 ) 
	  {
	    *(data->wavenumbers+(i*5)+counter) = atof(item);
	    counter++;
	  }
	}
      } 
    }


    /* no imaginary mode, reading is straightforward */
    else 
    {
      sscanf(buffer,"%s %lf %lf %lf %lf %lf",&word[0][0],
	&entry[0],&entry[1],&entry[2],&entry[3],&entry[4]); 
    
      /* save 'em */
      for ( k = 0; k < 5; ++k) 
      {
	*(data->wavenumbers+(i*5)+k) = entry[k]; 
      }
    }


    /* skip next line */
    eatline(data->file);


    /* next line contains the IR INTENSITIES */
    ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )

    /* scan the IR INTENSITIES */
    sscanf(buffer,"%s %s %lf %lf %lf %lf %lf",&word[0][0],&word[1][0],
	&entry[0],&entry[1],&entry[2],&entry[3],&entry[4]);
 

    /* save 'em */
    for ( k = 0; k < 5; ++k) 
    {
      *(data->intensities+(i*5)+k) = entry[k]; 
    }

    /* skip the next line */
    eatline(data->file);

    /* read the following five modes */
    for ( k = 0; k < data->numatoms; ++k)
    {
      /* x */
      ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )
      sscanf(buffer,"%s %s %s %lf %lf %lf %lf %lf",&word[0][0], 
	  &word[1][0], &word[2][0], &entry[0], &entry[1], &entry[2],
	  &entry[3], &entry[4]);

      /* store 'em */
      for ( j = 0; j < 5; ++j)
      {
	*(data->normal_modes+(3*k)+((i*5+j)*3*data->numatoms)) = 
	    entry[j];
      }


      /* y */
      ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )
      sscanf(buffer,"%s %lf %lf %lf %lf %lf",&word[0][0], &entry[0],
	  &entry[1],&entry[2], &entry[3],&entry[4]);

      /* store 'em */
      for ( j = 0; j < 5; ++j)
      {
	*(data->normal_modes+(3*k+1)+((i*5+j)*3*data->numatoms)) = 
	    entry[j];
      }


      /* z */
      ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )
      sscanf(buffer,"%s %lf %lf %lf %lf %lf",&word[0][0], &entry[0],
	  &entry[1], &entry[2], &entry[3],&entry[4]);

      /* store 'em */
      for ( j = 0; j < 5; ++j)
      {
	*(data->normal_modes+(3*k+2)+((i*5+j)*3*data->numatoms)) = 
	    entry[j];
      }
    }
  }


  /* read the remaining columns */
  if ( ( remaining_columns = (data->numatoms*3)%5) != 0 ) { 

    /* move to next set of modes */
    do 
    {
      ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )
      sscanf(buffer,"%s ",&word[0][0]);

    } while(strcmp(&word[0][0],"FREQUENCY:"));


    /* scan the frequencies; 
     * this requires some care since there might
     * be imaginary modes present which leads to
     * an additional char per imaginary mode per
     * line */

    /* check for imaginary modes */
    if ( strchr(buffer,'I') != NULL ) 
    {
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
        if ( *item == 'I' ) 
	{
	  data->nimag++;
	}
	else 
	{
	  /* save only the first remaining columns modes - 
	   * there NEVER should
	   * be more in any case, but just to make sure
	   * we don't overrun the array */
	  if ( counter < remaining_columns ) 
	  {
	    *(data->wavenumbers+(i*5)+counter) = atof(item);
	    counter++;
	  }
	}
      } 
    }


    /* no imaginary mode, reading is straightforward */
    else 
    {
      sscanf(buffer,"%s %lf %lf %lf %lf %lf",&word[0][0],
	&entry[0],&entry[1],&entry[2],&entry[3],&entry[4]); 
 
      /* save 'em */
      for ( k = 0; k < remaining_columns; ++k) 
      {
	*(data->wavenumbers+(i*5)+k) = entry[k]; 
      }
    }


    /* skip next line */
    eatline(data->file);

    
    /* next line contains the IR INTENSITIES */
    ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )

    /* scan the IR INTENSITIES */
    sscanf(buffer,"%s %s %lf %lf %lf %lf %lf",&word[0][0],&word[1][0],
	&entry[0],&entry[1],&entry[2],&entry[3],&entry[4]);
 
    /* save 'em */
    for ( k = 0; k < remaining_columns; ++k) 
    {
      *(data->intensities+(i*5)+k) = entry[k];
    }


    /* skip the next line */
    eatline(data->file);

    /* read the following five modes */
    for ( k = 0; k < data->numatoms; ++k)
    {
      /* x */
      ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )
      sscanf(buffer,"%s %s %s %lf %lf %lf %lf %lf",&word[0][0], 
	  &word[1][0], &word[2][0], &entry[0], &entry[1], &entry[2],
	  &entry[3], &entry[4]);

      /* store 'em */
      for ( j = 0; j < remaining_columns; ++j)
      {
	*(data->normal_modes+(3*k)+((i*5+j)*3*data->numatoms)) = 
	    entry[j];
      }
      
      /* y */
      ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )
      sscanf(buffer,"%s %lf %lf %lf %lf %lf",&word[0][0], &entry[0],
	  &entry[1],&entry[2], &entry[3],&entry[4]);
      
      /* store 'em */
      for ( j = 0; j < remaining_columns; ++j)
      {
	*(data->normal_modes+(3*k+1)+((i*5+j)*3*data->numatoms)) = 
	    entry[j];
      }
      /* z */
      ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )
      sscanf(buffer,"%s %lf %lf %lf %lf %lf",&word[0][0], &entry[0],
	  &entry[1], &entry[2], &entry[3],&entry[4]);
      
      /* store 'em */
      for ( j = 0; j < remaining_columns; ++j)
      {
	*(data->normal_modes+(3*k+2)+((i*5+j)*3*data->numatoms)) = 
	    entry[j];
      }
    }
  }

  /* release memory that is not needed any more */
  free(item);

  /* print brief message */
  printf("gamessplugin> Successfully scanned normal modes \n");

  return TRUE;
}



/***********************************************************
 *
 * this function reads in the cartesian hessian matrix 
 *
 * *********************************************************/
int get_cart_hessian(void *mydata)
{

  gamessdata *data = (gamessdata *)mydata;
  char word[4][BUFSIZ];
  char buffer[BUFSIZ];
  char dummy; 
  unsigned int i,j,k;
  float entry[6]; 

  /* intitialize arrays */
  buffer[0] = '\0';
  memset(entry, 0, sizeof(entry));
  for ( i = 0; i < 4; ++i) word[i][0] = '\0';


  /* at this point we need to rewind the file, since
   * in case that there is no internal Hessian stuff the
   * previous call to get_int_coords scanned the file
   * until EOF */
  rewind(data->file);


  /* look for CARTESIAN FORCE CONSTANT MATRIX */
  do 
  {
    ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )
    sscanf(buffer,"%s %s %s %s",&word[0][0],&word[1][0],
	&word[2][0],&word[3][0]);

  } while(strcmp(&word[0][0],"CARTESIAN") || 
          strcmp(&word[1][0],"FORCE") ||
	  strcmp(&word[2][0],"CONSTANT") ||
	  strcmp(&word[3][0],"MATRIX"));


  /* skip next 5 lines */
  for ( i=0; i<5; ++i) eatline(data->file);


  /* reserve memory for array; 
   * NOTE: this is a lower triangular matrix, but for now
   * we save it in an square matrix of dim(3Nx3N) to 
   * facilitate element access */
  data->carthessian = 
    (double *)calloc((data->numatoms*3)*(data->numatoms*3),
		     sizeof(double));

  
  /* make sure memory was allocated properly */
  if (data->carthessian == NULL) 
  {
    PRINTERR;
    return FALSE;
  }


  /* start scanning; the cartesian hessian matrix is a lower
   * triangular matrix, organized in rows of 6 */

  /* read blocks with complete rows of 6 */
  for ( i = 0; i < (int)(data->numatoms/2) ; ++i) 
  {
    for ( j = 0; j < (data->numatoms*3)-(i*6); ++j) 
    {
      ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )
 
      if ( j%3 == 0 ) 
      {
	sscanf(buffer,"%s %s %c %f %f %f %f %f %f",
	    &word[0][0],&word[1][0],&dummy,&entry[0],&entry[1],
	    &entry[2],&entry[3],&entry[4],&entry[5]);
      }
      else 
      {
	sscanf(buffer,"%1s %f %f %f %f %f %f",
       	  &dummy,&entry[0],&entry[1],&entry[2],&entry[3],&entry[4],
	  &entry[5]);
      }


      /* save entries (lower triangular matrix) in a 
       * square matrix */
      for ( k = 0; k <= ( j<5 ? j : 5 ) ; ++k) {
	*(data->carthessian+((j+(i*6))*3*data->numatoms)+
	    (k+(i*6))) = entry[k];
      }
    }


    /* skip the three line separating the matrix entries */
    eatline(data->file);
    eatline(data->file);
    eatline(data->file);
    eatline(data->file);
  }



  /* short message */
  printf("gamessplugin> Scanned Hessian in CARTESIAN coordinates\n");

  return TRUE;
}
  
  
  
/***********************************************************
 *
 * this function reads in the hessian in internal coordinates
 * as well as the internal coordinates 
 *
 * *********************************************************/
int get_int_coords(void *mydata)
{

  gamessdata *data = (gamessdata *)mydata;
  char word[5][BUFSIZ];
  char buffer[BUFSIZ];
  long position;
  int first, second, third, fourth;
  double value;
  double hess[5];
  unsigned int i = 0, j = 0, k = 0, l = 0;
  int n, dummy, remaining_blocks;


  /* initialize arrays */
  buffer[0] = '\0';
  memset(hess, 0, sizeof(hess));
  for ( i = 0; i < 5; ++i) word[i][0] = '\0';


  /* initialize counters */
  data->nintcoords = 0;
  data->nbonds = 0;
  data->nangles = 0;
  data->ndiheds = 0;
  data->nimprops = 0;


  /* look for list of INTERNAL COORDINATES */
  do 
  {
    ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )
    sscanf(buffer,"%s %s ",&word[0][0],&word[1][0]);

  } while(strcmp(&word[0][0],"INTERNAL") || 
          strcmp(&word[1][0],"COORDINATES"));

  
  /* skip next 5 lines */
  for ( i = 0; i < 5; ++i) eatline(data->file);


  /* remember current position so we can jump back */
  position = ftell(data->file);


  /* scan the next line */
  ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )
  n = sscanf(buffer,"%s %s", &word[0][0], &word[1][0]); 
 

  /* read line by line */
  while ( n != -1) 
  {
    /* start counting the number of internal coordinates */
    data->nintcoords++;


    /* count the number of bonds, angles, dihedrals */
    if (!strcmp(&word[1][0],"STRETCH")) 
    {
      data->nbonds++;
    }
    else if (!strcmp(&word[1][0],"BEND")) 
    {
      data->nangles++;
    }
    else if (!strcmp(&word[1][0],"TORSION")) 
    {
      data->ndiheds++;
    }
    else if (!strcmp(&word[1][0],"PLA.BEND")) 
    {
      data->nimprops++;
    }

            
    /* scan next line */
    ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )
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
  for ( i = 0; i < data->nbonds; ++i) 
  {
    ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )
    sscanf(buffer,"%s %s %d %d %lf", &word[0][0], &word[1][0], 
	&first, &second, &value);

    /* save 'em */
    *(data->bonds+2*i) = first;
    *(data->bonds+2*i+1) = second;
    *(data->internal_coordinates+i) = value;
  }


  /* scan the BENDS */
  for ( j = 0; j < data->nangles; ++j) 
  {
    ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )
    sscanf(buffer,"%s %s %d %d %d %lf", &word[0][0], &word[1][0], 
	&first, &second, &third, &value);

    /* save 'em */
    *(data->angles+3*j) = first;
    *(data->angles+3*j+1) = second;
    *(data->angles+3*j+2) = third;
    *(data->internal_coordinates+i+j) = value;
  }


  /* scan the TORSIONS */
  for ( k = 0; k < data->ndiheds; ++k) 
  {
    ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )
    sscanf(buffer,"%s %s %d %d %d %d %lf", &word[0][0], &word[1][0],
	&first, &second, &third, &fourth, &value);

    /* save 'em */
    *(data->dihedrals+4*k) = first;
    *(data->dihedrals+4*k+1) = second;
    *(data->dihedrals+4*k+2) = third;
    *(data->dihedrals+4*k+3) = fourth;
    *(data->internal_coordinates+i+j+k) = value;
  }


  /* scan the IMPROPERS */
  for ( l = 0; l < data->nimprops; ++l) 
  {
    ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )
    sscanf(buffer,"%s %s %d %d %d %d %lf", &word[0][0], &word[1][0],
	&first, &second, &third, &fourth, &value);

    /* save 'em */
    *(data->impropers+4*l) = first;
    *(data->impropers+4*l+1) = second;
    *(data->impropers+4*l+2) = third;
    *(data->impropers+4*l+3) = fourth;
    *(data->internal_coordinates+i+j+k+l) = value;
  }


  /* short message */
  printf("gamessplugin> Scanned %d INTERNAL coordinates \n",
      data->nintcoords);
  printf("gamessplugin>    %d BONDS \n",data->nbonds);
  printf("gamessplugin>    %d ANGLES \n",data->nangles);
  printf("gamessplugin>    %d DIHEDRALS \n",data->ndiheds);
  printf("gamessplugin>    %d IMPROPERS \n",data->nimprops);


  /* next read in the hessian in internal coordinates;
   * we would expect the matrix immediately after the
   * internal coordinates in the output files;
   * we check this first */
  eatline(data->file);


  ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )
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
  while ( sscanf(buffer,"%d %lf %lf %lf %lf %lf", 
	         &dummy, &hess[0], &hess[1], &hess[2], 
		 &hess[3], &hess[4]) != 6 ) 
  {
    ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )
  }

  
  /* reserve memory for inthessian array */
  data->inthessian = 
    (double *)calloc((data->nintcoords)*(data->nintcoords),
		     sizeof(double));


  /* make sure memory was allocated properly */
  if (data->inthessian == NULL)  
  {
    PRINTERR;
    return FALSE;
  }


  /* start scanning; GAMESS organized the output of the
   * internal HESSIAN in rows of 5 */

  /* read blocks with complete rows of 5 */
  for ( i = 0; i < (int)(data->nintcoords/5); ++i) 
  {
    for ( j = 0; j < data->nintcoords; ++j) 
    {
      sscanf(buffer,"%d %lf %lf %lf %lf %lf", &dummy, &hess[0], 
	      &hess[1], &hess[2], &hess[3], &hess[4]);

      /* save entries */
      for ( k = 0; k < 5; ++k) 
      { 
        *(data->inthessian+(j*data->nintcoords)+(i*5)+k) = hess[k];
      }

      /* next line */
      ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )
    }

    /* skip the two lines separating the matrix entries 
     * and scan next line */
    eatline(data->file);
    eatline(data->file);

    ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )
  }

  
  /* read the remaining block with less then 5 rows
   * if present */
  remaining_blocks = data->nintcoords%5;
  
  if ( remaining_blocks != 0 ) 
  {
    for ( j = 0; j < data->nintcoords; ++j) 
    {
      sscanf(buffer,"%d %lf %lf %lf %lf %lf", &dummy, &hess[0], 
	      &hess[1], &hess[2], &hess[3], &hess[4]);

      for ( k = 0; k < remaining_blocks; ++k) 
      { 
        *(data->inthessian+(j*data->nintcoords)+(i*5)+k) = hess[k];
      }

      ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )
    }
  }


  /* short message */
  printf("gamessplugin> Scanned Hessian in INTERNAL coordinates\n");



  /* finally, dump the diagonal elements of the hessian into the
   * force constant arrays, after converting the units 
   * appropriately;
   * BONDS are in HARTREE/BOHR**2
   * ANGLES,DIHEDRALS,IMPROPERS are in HARTREE/RADIAN**2 */

  
  /* alocate dynamic arrays */
  data->bond_force_const = 
    (double *)calloc(data->nbonds,sizeof(double));
 
  if ( data->bond_force_const == NULL ) 
  {
    PRINTERR;
    return FALSE;
  }


  data->angle_force_const =
    (double *)calloc(data->nangles,sizeof(double));

  if ( data->angle_force_const == NULL ) 
  {
    PRINTERR;
    return FALSE;
  }


  data->dihedral_force_const =
    (double *)calloc(data->ndiheds,sizeof(double));

  if ( data->dihedral_force_const == NULL ) 
  {
    PRINTERR;
    return FALSE;
  }


  data->improper_force_const =
    (double *)calloc(data->nimprops,sizeof(double));

  if ( data->improper_force_const == NULL ) 
  {
    PRINTERR;
    return FALSE;
  }


  /* scan the bonds */
  for ( i = 0; i < data->nbonds; ++i) 
  {
    *(data->bond_force_const + i) = 
      *(data->inthessian+(i*data->nintcoords)+i) 
      * HARTREE_TO_KCAL / BOHR_TO_ANGS / BOHR_TO_ANGS;

    printf("%3d (BOND) %2d - %2d : %lf (CHARMM) %lf \n",i, 
	*(data->bonds+2*i), *(data->bonds+2*i+1),
	*(data->bond_force_const +i),
	*(data->bond_force_const +i)*0.5); 
  }
  

  /* scan the angles */
  for ( j = i; j < i+(data->nangles); ++j) 
  {
    *(data->angle_force_const + (j-i)) = 
      *(data->inthessian+(j*data->nintcoords)+j) 
      * HARTREE_TO_KCAL;
    
     printf("%3d (ANGLE) %2d - %2d - %2d : %lf (CHARMM) %lf \n",j,
	 *(data->angles+3*(j-i)), *(data->angles+3*(j-i)+1), 
	 *(data->angles+3*(j-i)+2), 
	 *(data->angle_force_const + (j-i)),
	 *(data->angle_force_const + (j-i))*0.5);
  }


  /* scan the dihedrals */
  for ( k = j; k < j+(data->ndiheds); ++k) 
  {
    *(data->dihedral_force_const + (k-j)) = 
      *(data->inthessian+(k*data->nintcoords)+k)
      * HARTREE_TO_KCAL;
    
     printf("%3d (DIHEDRAL) %2d - %2d - %2d - %2d : %lf \n",k,
	 *(data->dihedrals+4*(k-j)), *(data->dihedrals+4*(k-j)+1),
	 *(data->dihedrals+4*(k-j)+2), *(data->dihedrals+4*(k-j)+3),
	 *(data->dihedral_force_const + (k-j))); 
  }


  /* scan the impropers */
  for ( l = k; l < k+(data->nimprops); ++l) 
  {
    *(data->improper_force_const + (l-k)) = 
      *(data->inthessian+(l*data->nintcoords)+l)
      * HARTREE_TO_KCAL;
    
   printf("%3d (IMPROPERS) %2d - %2d - %2d - %2d : %lf \n",l,
      *(data->impropers+4*(l-k)), *(data->impropers+4*(l-k)+1),
      *(data->impropers+4*(l-k)+2), *(data->impropers+4*(l-k)+3),
      *(data->improper_force_const + (l-k)));
  }


  /* DONE */
  return TRUE;
}



/*******************************************************
 *
 * this function reads in the guess options
 *
 * ******************************************************/
int get_guess(void *mydata)
{
  gamessdata *data = (gamessdata *)mydata;
  char word[2][BUFSIZ];
  char buffer[BUFSIZ];
  unsigned int i = 0;

  
  /* rewind the file
   * TODO: check if this is neccessary at this point */
  rewind(data->file);


  /* initialize buffers */
  buffer[0] = '\0';
  for ( i = 0; i < 2; ++i) word[i][0] = '\0';


  /* parse for GUESS field */
  do
  {
    ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )
    sscanf(buffer,"%s %s",&word[0][0],&word[1][0]);

  } while(strcmp(&word[0][0],"GUESS") || 
          strcmp(&word[1][0],"OPTIONS"));


  /* next line contains all we need */
  eatline(data->file);
  ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )
  sscanf(buffer,"%s %s",&word[0][0],&word[1][0]);

  /* store it */
  strncpy(data->guess,(&word[1][0])+1,sizeof(data->guess));

  
  /* short info */
  printf("gamessplugin> Run was performed with GUESS = %s \n",
	  data->guess);

  return TRUE;
}



/*******************************************************
 *
 * this function reads in the basis set data 
 *
 * ******************************************************/
int get_basis(void *mydata)
{

  gamessdata *data = (gamessdata *)mydata;
  float *basis;
  int *basis_counter;
  int *atomic_shells;
  int *shell_primitives;
  char *orbital_symmetry;
  char buffer[BUFSIZ];
  char word[3][BUFSIZ];
  unsigned int i = 0; 
  int counter = 0, oldcounter = 0, nat = 0;

  
  /* initialize buffers */
  buffer[0] = '\0';
  for ( i = 0; i < 3; ++i ) word[i][0] = '\0';


  /* reserve space for pointer to array containing basis
   * info, i.e. contraction coeficients and expansion 
   * coefficients; need 2 entries per basis function, i.e.
   * exponent and contraction coefficient; also,
   * allocate space for the basis_counter (i.e. # of
   * basis functions per atoms as well as the array
   * holding the orbital symmetry information per primitive
   * Gaussian
   * Finally, initialize the arrays holding the number of 
   * shells per atom and the number of primitives per shell*/
  basis = (float *)calloc(2*MAXBASISFUNCTIONS,sizeof(float));

  /* make sure memory was allocated properly */
  if (basis == NULL) 
  {
    PRINTERR;
    return MOLFILE_ERROR;
  }


  basis_counter = (int *)calloc(MAXQMATOMS,sizeof(int));

  /* make sure memory was allocated properly */
  if (basis_counter == NULL) 
  {
    PRINTERR;	    
    return MOLFILE_ERROR;
  }


  orbital_symmetry = (char *)calloc(MAXBASISFUNCTIONS,
                     sizeof(char));
  
  /* make sure memory was allocated properly */
  if (orbital_symmetry == NULL) 
  {
    PRINTERR; 
    return MOLFILE_ERROR;
  }

  atomic_shells = (int *)calloc(MAXBASISFUNCTIONS, sizeof(int));

  /* make sure memory was allocated properly */
  if (atomic_shells == NULL) 
  {
    PRINTERR; 
    return MOLFILE_ERROR;
  }


  shell_primitives = (int *)calloc(MAXBASISFUNCTIONS, sizeof(int));

  /* make sure memory was allocated properly */
  if (shell_primitives == NULL) 
  {
    PRINTERR;	    
    return MOLFILE_ERROR;
  }


  /* store pointers in struct gamessdata */
  data->basis = basis;
  data->basis_counter = basis_counter;
  data->orbital_symmetry = orbital_symmetry;
  data->atomic_shells = atomic_shells;
  data->shell_primitives = shell_primitives;

  /* next parse through output file and look for the 
   * atomic basis set, i.e.
   *
   * ATOMIC BASIS SET
   *
   */

  rewind(data->file);
 
  /* initialize the per-atom basis function counter */
  oldcounter = 0;
  counter = 0;

  /* start reading the basis function field */
  do
  {
    ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )
    sscanf(buffer,"%s %s %s",&word[0][0],&word[1][0],&word[2][0]);

  } while(strcmp(&word[0][0],"ATOMIC") || 
          strcmp(&word[1][0],"BASIS") || 
          strcmp(&word[2][0],"SET"));


  /* skip the next 8 lines */
  for(i=0; i<7; i++)
  {
    eatline(data->file);
  }


  /* the number of atoms determines how often we have
   * to loop through the basis function array */
  nat = data->numatoms;
 

  /* now loop through the basis function section of the
   * output file and read in the exponents/coefficients */
  for(i=0; i<nat; i++)
  {
    /* here I go trough the basis set information on a 
     * per atom basis and put all the exps and coefficients
     * in the right place */
    counter = atomic_basis(oldcounter, data, basis, 
	 orbital_symmetry, atomic_shells,shell_primitives); 


    /* if atomic_basis() returns -1 we need to bomb
     * out */
    if ( counter == -1 ) return FALSE;


    /* store number of basis-functions in array
    * basis_counter */
    *basis_counter = counter;


    /* increase pointer to next entry */
    basis_counter++;


    /* update the pointer for the basis arrays */ 
    oldcounter += counter;


    /* increase atomic_shell pointer to hold value for
    * next atom */
    atomic_shells++;
  } 


  /* store the total number of basis functions */
  data->num_basis_funcs = oldcounter;


  /* short info */ 
  printf("gamessplugin> Parsed %d uncontracted basis functions. \n",
      data->num_basis_funcs);

  return TRUE;
}



/*******************************************************
 *
 * function renormalizing basis set coefficients of
 * primitives
 *
 *******************************************************/
float renorm_coefficient(float coefficient, float exponent, 
	char orb)
{
  /* normalization for S shells */
  if ( orb == 'S' )
  {
    return pow(exponent,1.5) *  pow((8/pow(MY_PI,3)),0.5)
      * coefficient;
  }

  /* normalization for P shells */
  else if ( orb == 'P' )
  {
    return pow(exponent,2.5) * pow((128/pow(MY_PI,3)),0.5) 
      * coefficient;
  }

  /* TODO: the renormalization factor for D and F shells
   * are taken from MOLDEN; check them! */

  /* normalization for D shells */
  else if ( orb == 'D' )
  {
    return pow(exponent,3.5) * pow((2048/9/pow(MY_PI,3)),0.5)
	* coefficient;
  }

  /* normalization for F shells */
  else if ( orb == 'F' )
  {
    return pow(exponent,4.5) * 1.875 * pow((512/pow(MY_PI,3)),0.5)
	* coefficient;
  }


  /* return 0 otherwise */
  return 0;
}



/******************************************************
 *
 * this function scans a single line in the basis set
 * data section 
 *
 * *****************************************************/
int atomic_basis(int oldcounter, void *mydata, float *basis, 
             char *orbital_symmetry, int *atomic_shells, 
             int *shell_primitives)
{
  gamessdata *data = (gamessdata *)mydata;
  char dummy[BUFSIZ];
  char buffer[BUFSIZ];
  char orbsym;
  float exponent = 0.0; 
  float contract_c1 = 0.0, contract_c2 = 0.0;
  int counter, success;
  static int prim_count = -1;
  int prev_shell = 0;
  int shell;


  /* initialize arrays */
  dummy[0] = '\0';
  buffer[0] = '\0';


  /* initialize basis function counter for current atom */
  counter = oldcounter;


  /* skip one line */
  eatline(data->file);


  /* loop through basis set info for next atom */
  do
  {
    ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )
    success = sscanf(buffer,"%d %c %s %f %f %f", &shell, \
	             &orbsym, dummy, \
	             &exponent, &contract_c1, &contract_c2); 

    
    /* make sure we have eiter S, L, P, or D shells;
     * otherwise we bomb out;
     * TODO: Add support for F SHELLS */
    if ( (success == 5  || success == 6 ) && 
         (orbsym != 'S' && orbsym != 'L' && orbsym != 'P' && 
	  orbsym != 'D' && orbsym != 'F' ))
    {
      printf("gamessplugin> WARNING ... %c shells are not supported \n", orbsym);
      return -1;
    }


    /* when this function is called prev_shell is
     * initialized to zero;
     * Next we compare it to the shell number just read in;
     * in case they are different we increase atomic_shell
     * by one and then update prev_shell with the new shell
     * value */
   if ( prev_shell != shell && success > 0 )
    {
      prim_count++;
      (*atomic_shells)++;
      prev_shell = shell;
    }
   

    /* store in basis array and increase the counter */ 
    switch(success)
    {
      case(5):
	{
	  /* store exponent */
	  *(basis+(2*counter)) = exponent;


	  /* renormalize and store coefficient */
	  *(basis+(2*counter+1)) = 
	      renorm_coefficient(contract_c1,exponent,orbsym);


	  /* store orbital symmetry */
	  *(orbital_symmetry+counter) = orbsym;


	  /* increase counter */
	  counter++;


	  /* increase counter for number of primitives
	   * for current shell */
	  (*(shell_primitives+prim_count))++;

	  break;
	}
      case(6):
	{ 
	  /* in the case of two contraction
	   * coefficients we expect "L" as orbital 
	   * symmetry, meaning that the first coefficient
	   * has S and the second P symmetry. Otherwise
	   * print a warning */
	  if ( orbsym != 'L' ) 
	  {
	    printf("gamessplugin>\n");
	    printf("gamessplugin> WARNING ... Non SP-Shell %c", 
		orbsym);
	    printf(" interpreted as SP-Shell\n");
	    printf("gamessplugin>\n");
	  }


	  /* line has two coefficient pairs,
	   * store the first one */
	  *(basis+(2*counter)) = exponent;

	  
	  /* renormalize and store coefficient */
	  *(basis+(2*counter+1)) = 
	      renorm_coefficient(contract_c1,exponent,'S');


	  /* the first coefficient has S symmetry */
	  *(orbital_symmetry+counter) = 'S';


	  /* increase counter */
	  counter++;


          /* increase counter for number of primitives
	   * for current shell */
	  (*(shell_primitives+prim_count))++;

	  
          /* same procedure for the second pair */
	  *(basis+(2*counter)) = exponent;


	  /* renormalize and store coefficient */
	  *(basis+(2*counter+1)) = 
	      renorm_coefficient(contract_c2,exponent,'P');


	  /* the second coefficient has P symmetry */
          *(orbital_symmetry+counter) = 'P';

	  counter++;
	  (*(shell_primitives+prim_count))++;

	  break;
	}
       case(-1):
	{
	  break; 
	}
       case(0):
	{
	  break;
	}
       case(3):
	{
	  break; 
	} 
       default:    /* this should never happen */
	{
        return MOLFILE_ERROR; 
	break;
	}
    }
  }  while((success!=0) && (success!=3));


  /* return only the counter difference */
  return (counter - oldcounter);
}


/********************************************************
 *
 * this function reads the number of A/B orbitals 
 *
 ********************************************************/
int get_num_orbitals(void *mydata)
{

  gamessdata *data = (gamessdata *)mydata; 
  char buffer[BUFSIZ];
  char word[7][BUFSIZ];
  unsigned int i;

  /* initialize buffers */
  buffer[0] = '\0';
  for ( i = 0; i < 7; ++i) word[i][0] = '\0';


  /* we start reading the file from the beginning */
  rewind(data->file);


  /* look for the orbital/charge/... info section */
  do
  {
    ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )
    sscanf(buffer,"%s %s %s %s",&word[0][0], &word[1][0],
	&word[2][0], &word[3][0]);

  } while(strcmp(&word[0][0],"TOTAL") || 
          strcmp(&word[1][0],"NUMBER") || 
          strcmp(&word[2][0],"OF")    || 
	  strcmp(&word[3][0],"BASIS"));

  
  /* go ahead reading the info */
  ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )
  sscanf(buffer,"%s %s %s %s %s %s %s %d", &word[0][0],
      &word[1][0], &word[2][0], &word[3][0], &word[4][0],
      &word[5][0], &word[6][0], &(data->num_gauss_basis_funcs));



  /* read the number of electrons */
  ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )
  sscanf(buffer,"%s %s %s %s %d", &word[0][0], &word[1][0],
      &word[2][0], &word[3][0], &(data->num_electrons));



  /* read the charge of the molecule */
  ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )
  sscanf(buffer,"%s %s %s %s %d", &word[0][0], &word[1][0],
      &word[2][0], &word[3][0], &(data->totalcharge));



  /* read the multiplicity of the molecule */
  ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )
  sscanf(buffer,"%s %s %s %d", &word[0][0], &word[1][0],
       &word[2][0], &(data->multiplicity));

  
  /* read number of A orbitals */
  ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )


  /* note the different number of items per line for A/B orbitals
   * due to "(ALPHA)" and "(BETA )" !! */
  sscanf(buffer,"%s %s %s %s %s %s %d", &word[0][0], &word[1][0],
      &word[2][0], &word[3][0], &word[4][0], &word[5][0],
      &(data->num_orbitals_A)); 
  
  
  
  /* read number of B orbitals */
  ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )
  sscanf(buffer,"%s %s %s %s %s %s %s %d", &word[0][0], &word[1][0],
      &word[2][0], &word[3][0], &word[4][0], &word[5][0], 
      &word[6][0], &(data->num_orbitals_B)); 


  
  /* short message */
  printf("gamessplugin> Number of Electrons: %d \n",
      data->num_electrons);

  printf("gamessplugin> Charge of Molecule : %d \n",
      data->totalcharge);

  printf("gamessplugin> Multiplicity of Wavefunction: %d \n",
      data->multiplicity);

  printf("gamessplugin> Number of A / B orbitals: %d / %d \n",\
      data->num_orbitals_A, data->num_orbitals_B);

  printf("gamessplugin> Number of gaussian basis functions: %d \n",\
      data->num_gauss_basis_funcs);

 
  return TRUE;
}


/*********************************************************
 *
 * this function reads the actual wavefunction, which is
 * punched at the end of the log file
 *
 **********************************************************/
int get_wavefunction(void *mydata)
{
  gamessdata *data = (gamessdata *)mydata; 
  float *wave_function;
  float *orbital_energy;
  char buffer[BUFSIZ];
  char word[5][BUFSIZ];
  char dummy1, dummy2, dummy3, dummy4;
  int orbital_counter = 0;
  unsigned int i = 0, j = 0, num_values = 0;
  int length = 0;


  /* initialize buffers */
  buffer[0] = '\0';
  for ( i = 0; i < 5; ++i ) word[i][0] = '\0';

  
  /* reserve space for arrays storing wavefunction and orbital
   * energies 
   * For the wavefunction I reserve twice the number of A orbitals
   * times MAXBASISFUNCTIONS;
   * I should check the gamess sources how much are printed out,
   * but this should cover me for now. 
   * Accordingly, for the energies I use the number of A orbitals */
  wave_function = (float *)calloc(data->num_gauss_basis_funcs * \
                  data->num_gauss_basis_funcs,sizeof(float)); 

  /* make sure memory was allocated properly */
  if (wave_function == NULL) 
  {
    PRINTERR;	    
    return FALSE;
  }
  
  orbital_energy = (float *)calloc(data->num_gauss_basis_funcs, \
                      sizeof(float));

  /* make sure memory was allocated properly */
  if (orbital_energy == NULL) 
  {
    PRINTERR; 
    return FALSE;
  }

  /* store the pointers in gamessdata */
  data->wave_function = wave_function;
  data->orbital_energy = orbital_energy;

  
  /* rewind the file before starting to look for
   * the wavefunction */

  rewind(data->file);
  
  do
  {
    ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )
    length = sscanf(buffer,"%s %s",&word[0][0],&word[1][0]);

  } while(strcmp(&word[0][0],"EIGENVECTORS") || length != 1);

  
  /* skip the next three lines */
  eatline(data->file);
  eatline(data->file);
  eatline(data->file);


  /* go ahead */
  ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )
  num_values = sscanf(buffer,"%s %s %s %s %s",&word[0][0],
      &word[1][0],&word[2][0],&word[3][0],&word[4][0]);


  while(strcmp(&word[0][0],"TOTAL"))
  {
    /* store the orbital energies in the appropriate arrays 
     * read them until we encounter an empty string */
    for(i=0; i<num_values; i++)
    {
      *(orbital_energy+i) = atof(&word[i][0]);
      orbital_counter++;
    }

    /* increase orbital energy pointer */
    orbital_energy = orbital_energy+5;

    /* skip the next line */
    eatline(data->file);

    /* now read in the wavefunction */
    /* markus: go ahead here, instead of just looping
     * over the lines actually store them :) */
    for(i=0; i< data->num_gauss_basis_funcs; i++)
    {
      ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )
      
      /* now read in the wavefunction coefficients for 5
       * orbitals at a time line by line */
      num_values = sscanf(buffer,"%s %s %s %s %s %s %s %s %s", 
	         &dummy1,&dummy2,&dummy3,&dummy4,&word[0][0], 
		 &word[1][0],&word[2][0],&word[3][0],&word[4][0]);


      /* each orbital has data->num_gauss_basis_funcs entries, 
       * hence we have to use this number as offset when storing 
       * them in groups of five */
      for(j=0 ; j<num_values-4; j++) {

	*(wave_function+(j*data->num_gauss_basis_funcs)+i) \
	        = atof(&word[j][0]);
      }
    }

    /* move wavefunction pointer to start of next five orbitals */
    wave_function = wave_function + 5*data->num_gauss_basis_funcs;

    /* skip next line */
    eatline(data->file);
    eatline(data->file);


    /* read next line to allow do/while loop to decide if
     * wavefunction data is over */
    ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )
    num_values = sscanf(buffer,"%s %s %s %s %s",&word[0][0],
	&word[1][0], &word[2][0], &word[3][0], &word[4][0]);
  } 

  /* store the number of orbitals read in */
  data->orbital_counter = orbital_counter;

  
  /* at this point we should have proper wavefunction information
   * in our arrays */
  data->got_wavefunction = TRUE;

  
  /* short message */ 
  printf("gamessplugin> Number of orbitals scanned: %d \n",\
             data->orbital_counter);


  return TRUE;
}



/*************************************************************
 *
 * this subroutine scans all the orbitals to find the HOMO 
 *
 *************************************************************/
int find_homo(void *mydata) {

  gamessdata *data = (gamessdata *)mydata; 
  float *orbital_energy;
  int orbital_counter;

  /* initialize the counter pointing to the HOMO */
  unsigned int homo_counter = 0;

  /* loop through the orbital energies and find the smallest
   * negative value, i.e. HOMO. Actually, for pathological
   * systems identifications of the HOMO in this way might
   * be flawed, since higher energy orbitals could be occupied
   * as well, but for now this should do */
  orbital_energy = data->orbital_energy;
  orbital_counter = data->orbital_counter;

  while ( *orbital_energy < 0 )
  {
    orbital_energy++;
    homo_counter++;
  }

  /* there is probably a lot that can go wrong here;
   * let's at least check if the HOMO is the last orbital
   * present in which case there is probably something
   * fishy going on */
  if ( homo_counter == orbital_counter)
  {
    printf("gamessplugin> WARNING: HOMO is last orbital");
  }

  
  /* store homo index */
  data->homo_index = homo_counter;

  
  /* short message */
  printf("gamessplugin> Identified orbital #%d (%f au) as HOMO \n",
	data->homo_index, 
	*(data->orbital_energy+(data->homo_index)-1)); 
 

  return TRUE;
}



/**********************************************************
 *
 * this function generates the computational grid for
 * the orbitals; as a proof of concept case we do the
 * HOMO first 
 *
 ***********************************************************/
int orbital_grid_driver(void *mydata)
{
  gamessdata *data = (gamessdata *)mydata; 

  /* in order to set up the grid properly, we
   * have to know the dimensions of the system.
   * Hence let's determine it ! */
  if (get_system_dimensions(data) == FALSE) return FALSE;
 

  /* finally, we're ready to evaluate the value
   * of the orbital at the computational grid
   * points */
  if (calculate_orbital(data) == FALSE ) return FALSE; 


  /* at this point we should have proper volumetric metadata/
     * data (otherwise it is the fault of generate_orbital_grid)
     * and can hence allow VMD to parse the apropriate arrays
     * in the metadata/data routines */
  data->have_volumetric = 1; 


  return TRUE;
}


/*********************************************************
 *
 * this function creates the computational grid given the
 * system dimensions 
 *
 **********************************************************/
int calculate_orbital(void *mydata) 
{
  gamessdata *data = (gamessdata *)mydata;
  float *orbital_grid;
  float *xyzdim;
  molfile_volumetric_t *vol_metadata;
  float *wave_function;
  float grid_x = 0.0, grid_y = 0.0, grid_z = 0.0;
  int numx, numy, numz, numxy, maxpoints;
  unsigned int nx = 0, ny = 0, nz = 0, i = 0;
  float *center;
  int success[3] = { 0, 0, 0};
  float grid_size = 1.0; 
  float buffer_region[3] = {3.0, 3.0, 3.0};
  int grid_bad[3] = { 1, 1, 1};
  float grid_dim[3]; /* gridsize in A */
  int shrink_grid_size = 0;
  int voxel[3];
  int tot_voxel;

  /* request memory for volume metadata */
  vol_metadata = (molfile_volumetric_t *)calloc(1, \
		sizeof(molfile_volumetric_t)); 

  /* make sure the call was succesful */
  if (vol_metadata == NULL) 
  {
     PRINTERR;
     return FALSE;
  }

  /* store pointer to volumetric metadata */
  data->vol = vol_metadata;


  /* retrieve the pointer for the temp arrays */
  /* temp_data = data->temporary; */


  /* retrieve array with # of primitives available
   * for each indivual atom of the system */
  /* basis_counter = data->basis_counter; */


  /* retrieve the pointer for the grid dimensions */
  xyzdim = data->system_dimensions;


  /* retrieve pointer for wavefunction */
  wave_function = data->wave_function;


  /* move wavefunction pointer to HOMO */
  wave_function += (data->num_gauss_basis_funcs*
      (data->homo_index - 1));


  /* retrieve orbital symmetry of shel primitives */
  /* orbital_symmetry = data->orbital_symmetry; */


  /* retrieve pointer to system center */
  center = data->system_center;


  /* now we scan the surface of the grid and see if there
   * are any values > 0.01; if not, we reduce the buffer
   * region by 1A and try again */
  printf("gamessplugin>\n");
  printf("gamessplugin> Optimizing size of orbital grid. \n");
 

  /* determine an optimal grid dimension and grid spacing */
  while ( grid_bad[0] || grid_bad[1] || grid_bad[2] )
  {
    /* examine the current grid and shrink it if possible */

    for ( i = 0 ; i < 3 ; ++i)
    {
      grid_dim[i] = (xyzdim[i+3] - xyzdim[i]); 

      /* pad with a buffer region subject to the constraint
       * the grid_dim[i] has an integer number of voxels */
        
      grid_dim[i] = grid_dim[i] - fmod(grid_dim[i],2*grid_size) 
		    + 2.0 * buffer_region[i]; 

      vol_metadata->origin[i] = center[i] - grid_dim[i]/2.0; 
    }

    /* make results available as volumetric metadata:
     * the {xyz}axis point in the xyz direction with
     * lengths determined by the size of the grid in
     * each direction;
     * the {xyz}size is given by the size of the grid
     * plus a buffer region BUFFER_REGION
     * in the respective direction divided by the
     * GRIDSIZE (i.e. number of voxels) */
      
    vol_metadata->xaxis[0] = grid_dim[0];
    vol_metadata->yaxis[1] = grid_dim[1];
    vol_metadata->zaxis[2] = grid_dim[2];
    vol_metadata->xsize = (int)(grid_dim[0] / grid_size)+1;
    vol_metadata->ysize = (int)(grid_dim[1] / grid_size)+1;
    vol_metadata->zsize = (int)(grid_dim[2] / grid_size)+1; 

    /* don't have color information */
    vol_metadata->has_color = 0;
    
    /* default name for data set */
    sprintf(vol_metadata->dataname,"wavefunction data"); 

     
    /* calculate the number of grid points in each dimension
     * as well as the total number of points needed
     * in order to dimension the orbital_grid
     * array properly */
    numx = vol_metadata->xsize;
    numy = vol_metadata->ysize;
    numz = vol_metadata->zsize;


    /* scan xy surfaces */
    nz = 0;

    if ( !success[0] )
    {
      for ( nx = 0 ; nx < numx ; ++nx)
      {
	if ( !success[0] ) 
	{
	  for ( ny = 0 ; ny < numy ; ++ny)
  	  {
  	    /* calculate the xyz coordinate of the current
 	     * grid point */
	    grid_x = vol_metadata->origin[0] + nx * grid_size;
	    grid_y = vol_metadata->origin[1] + ny * grid_size;
	    grid_z = vol_metadata->origin[2] + nz * grid_size;

	
	    /* calculate the value of the wavefunction of the
	     * selected orbital at the current grid point */
	    if ( fabs(orbital_at_grid_xyz(data, wave_function, 
		  grid_size, grid_x, grid_y, grid_z)) > 0.01 )
	    {
	      success[0] = 1;
	      break;
	    }
	  } 
	}
      }
    } 
    	
    /* scan back side of xy plane */
    nz = numz - 1;

    if ( ! success[0] )
    {
      for ( nx = 0 ; nx < numx ; ++nx)
      {
	if ( !success[0] ) 
	{
	  for ( ny = 0 ; ny < numy ; ++ny)
	  {
	    /* calculate the xyz coordinate of the current
	    * grid point */
	    grid_x = vol_metadata->origin[0] + nx * grid_size;
	    grid_y = vol_metadata->origin[1] + ny * grid_size;
	    grid_z = vol_metadata->origin[2] + nz * grid_size;
	  
	
	    /* calculate the value of the wavefunction of the
	     * selected orbital at the current grid point */
	    if ( fabs(orbital_at_grid_xyz(data, wave_function, 
		  grid_size, grid_x, grid_y, grid_z)) > 0.01 )
	    {
	      success[0] = 1;
	      break;
	    }
	  } 
	}
      }
    }


    /* adjust buffer region in z direction */
    if ( success[0] )
    {
      grid_bad[2] = 0;

    }
    else
    {
      buffer_region[2] -= 0.5;
    }



    /* scan xz surfaces */
    ny = 0;

    if ( !success[1] )
    {
      for ( nx = 0 ; nx < numx ; ++nx)
      {
	if ( !success[1] ) 
	{
	  for ( nz = 0 ; nz < numz ; ++nz)
	  {
	    /* calculate the xyz coordinate of the current
	    * grid point */
	    grid_x = vol_metadata->origin[0] + nx * grid_size;
	    grid_y = vol_metadata->origin[1] + ny * grid_size;
	    grid_z = vol_metadata->origin[2] + nz * grid_size;
	  
	
	    /* calculate the value of the wavefunction of the
	     * selected orbital at the current grid point */
	    if ( fabs(orbital_at_grid_xyz(data, wave_function, 
		  grid_size, grid_x, grid_y, grid_z)) > 0.01 )
	    {
	      success[1] = 1;
	      break;
	    }
	  } 
	}
      }
    }


    /* scan back side of xz surface */	    
    ny = numy - 1; 
    
    if ( !success[1] )
    {
      for ( nx = 0 ; nx < numx ; ++nx)
      {
	if ( !success[1] ) 
	{
	  for ( nz = 0 ; nz < numz ; ++nz)
	  {
	    /* calculate the xyz coordinate of the current
	    * grid point */
	    grid_x = vol_metadata->origin[0] + nx * grid_size;
	    grid_y = vol_metadata->origin[1] + ny * grid_size;
	    grid_z = vol_metadata->origin[2] + nz * grid_size;
	  
	
	    /* calculate the value of the wavefunction of the
	     * selected orbital at the current grid point */
	    if ( fabs(orbital_at_grid_xyz(data, wave_function, 
		  grid_size, grid_x, grid_y, grid_z)) > 0.01 )
	    {
	      success[1] = 1;
	      break;
	    }
	  } 
	}
      }
    }

    /* adjust buffer region in y direction */
    if ( success[1] )
    {
      grid_bad[1] = 0;
    }
    else
    {
      buffer_region[1] -= 0.5;
    }


    /* scan yz surfaces */
    nx = 0;

    if ( !success[2] )
    {
      for ( ny = 0 ; ny < numy ; ++ny)
      {
	if ( !success[2] ) 
	{
	  for ( nz = 0 ; nz < numz ; ++nz)
	  {
	    /* calculate the xyz coordinate of the current
	    * grid point */
	    grid_x = vol_metadata->origin[0] + nx * grid_size;
	    grid_y = vol_metadata->origin[1] + ny * grid_size;
	    grid_z = vol_metadata->origin[2] + nz * grid_size;
	  
	
	    /* calculate the value of the wavefunction of the
	     * selected orbital at the current grid point */
	    if ( fabs(orbital_at_grid_xyz(data, wave_function, 
		  grid_size, grid_x, grid_y, grid_z)) > 0.01 )
	    {
	      success[2] = 1;
	      break;
	    }
	  } 
	}
      }
    }


    /* scan back side of yz surface */	    
    nx = numx - 1; 
    
    if ( !success[2] )
    {
      for ( ny = 0 ; ny < numy ; ++ny)
      {
	if ( !success[2] ) 
	{
	  for ( nz = 0 ; nz < numz ; ++nz)
	  {
	    /* calculate the xyz coordinate of the current
	    * grid point */
	    grid_x = vol_metadata->origin[0] + nx * grid_size;
	    grid_y = vol_metadata->origin[1] + ny * grid_size;
	    grid_z = vol_metadata->origin[2] + nz * grid_size;
	  
	
	    /* calculate the value of the wavefunction of the
	     * selected orbital at the current grid point */
	    if ( fabs(orbital_at_grid_xyz(data, wave_function, 
		  grid_size, grid_x, grid_y, grid_z)) > 0.01 )
	    {
	      success[2] = 1;
	      break;
	    }
	  } 
	}
      }
    }

    /* adjust buffer region in x direction */
    if ( success[2] )
    {
      grid_bad[0] = 0;
    }
    else
    {
      buffer_region[0] -= 0.5;
    }

  }

  /* at this point we have optimized the dimensions of
   * the grid and we can now chose an optimal low value for
   * the grid size depending such that we either have 
   * less than MAX_GRIDPOINTS or a value of 0.1 A */
  shrink_grid_size = 1;

  while ( shrink_grid_size ) 
  {
    for ( i = 0 ; i < 3 ; ++i)
    {
      grid_dim[i] = (xyzdim[i+3] - xyzdim[i]);

      voxel[i] = (int)((grid_dim[i] + 2.0 * buffer_region[i]) /
	  grid_size);
    }

    tot_voxel = voxel[0] * voxel[1] * voxel[2];


    if (tot_voxel > MAX_GRIDPOINTS)
    {
      shrink_grid_size = 0;

      /* we need to bump up the gridsize again
       * since at this point we've already gone
       * beyond MAX_GRIDPOINTS */
      grid_size += 0.1;
    }
    else
    {
      if ( grid_size > 0.1 )
      {
	grid_size -= 0.1;
      }
      else
      {
	shrink_grid_size = 0;
      }
    }
  }

  for ( i = 0 ; i < 3 ; ++i)
  {
    /* pad with a buffer region subject to the constraint
     * the grid_dim[i] has an integer number of voxels */
        
    grid_dim[i] = grid_dim[i] - fmod(grid_dim[i],2*grid_size) 
		    + 2.0 * buffer_region[i]; 

    vol_metadata->origin[i] = center[i] - grid_dim[i]/2.0; 
  }

  /* make results available as volumetric metadata:
   * the {xyz}axis point in the xyz direction with
   * lengths determined by the size of the grid in
   * each direction;
   * the {xyz}size is given by the size of the grid
   * plus a buffer region BUFFER_REGION
   * in the respective direction divided by the
   * GRIDSIZE (i.e. number of voxels) */
      
  vol_metadata->xaxis[0] = grid_dim[0];
  vol_metadata->yaxis[1] = grid_dim[1];
  vol_metadata->zaxis[2] = grid_dim[2];
  vol_metadata->xsize = (int)(grid_dim[0] / grid_size)+1;
  vol_metadata->ysize = (int)(grid_dim[1] / grid_size)+1;
  vol_metadata->zsize = (int)(grid_dim[2] / grid_size)+1; 


  /* calculate the number of grid points in each dimension
   * as well as the total number of points needed
   * in order to dimension the orbital_grid
   * array properly */
  numx = vol_metadata->xsize;
  numy = vol_metadata->ysize;
  numz = vol_metadata->zsize;
  maxpoints = numx * numy * numz;


  /* useful shorthand notation needed for storage
   * of orbital values in an array */
  numxy = numx*numy;


  /* store number of gridpoints */
  data->num_gridpoints = maxpoints;


  /* now we can reserve the space for the grid array */
  orbital_grid = (float *)calloc(maxpoints,sizeof(float));


  /* make sure we have enough memory */
  if (orbital_grid== NULL) 
  {
    PRINTERR; 
    return FALSE;
  }

  /* store the pointer to the grid */
  data->orbital_grid = orbital_grid;


  /* let's give the user a warning, since the calculation
   * could take a while, otherwise they might think the
   * system is borked */
  printf("gamessplugin> Calculating orbital grid. \n");
  printf("              Please be patient ....... \n");


  /* now we start one heck of a large loop going through the
   * whole grid calculating the value of the orbital at each
   * point */

  /* start looping and calculating ........... */
  for ( nx = 0 ; nx < numx ; ++nx)
  {
    for ( ny = 0 ; ny < numy ; ++ny)
    {
      for ( nz = 0 ; nz < numz ; ++nz)
      {

	/* calculate the xyz coordinate of the current
	 * grid point */
	grid_x = vol_metadata->origin[0] + nx * grid_size;
	grid_y = vol_metadata->origin[1] + ny * grid_size;
	grid_z = vol_metadata->origin[2] + nz * grid_size;

	
	/* calculate the value of the wavefunction of the
	 * selected orbital at the current grid point */
	*(orbital_grid + nx + ny*numx + nz*numxy) = 
	    orbital_at_grid_xyz(data, wave_function, grid_size, 
	       	grid_x, grid_y, grid_z);
      }
    } 
  }

  /* everything went fine */ 
  return TRUE;
}



/*********************************************************
 *
 * this function provides the system dimensions 
 *
 *********************************************************/
int get_system_dimensions(void *mydata) {

  gamessdata *data = (gamessdata *)mydata; 
  gamess_temp *temp_data;
  unsigned int i = 0;

  /* array holding geometric center of system */
  float *center;   
 
  /* array containing system dimensions in the form
   * [x_min,y_min,....z_max] */
  float *xyzdim;   
  
  /* request memory for center array */
  center = (float *)calloc(3,sizeof(float));

  /* make sure memory could be allocated */
  if (center == NULL)
  {
    PRINTERR;
    return FALSE;
  }


  /* request memory for xyzdim array */
  xyzdim = (float *)calloc(6,sizeof(float));

  /* make sure the call was succesful */
  if (xyzdim == NULL) 
  {
    PRINTERR;
    return FALSE;
  }


  /* retrieve pointer of gamess_temp struct */
  temp_data = data->temporary;

  /* store pointer to volumetric metadata */
  /* data->vol = vol_meta_data; */

  /* store pointer to center of system */
  data->system_center = center;


  /* determine the origin of the system; this is needed
   * to properly initialize the volumetric metadata */
  for ( i = 0 ; i != data->numatoms ; ++i) 
  {
    center[0] = ( center[0] + (temp_data+i)->x ); 
    center[1] = ( center[1] + (temp_data+i)->y ); 
    center[2] = ( center[2] + (temp_data+i)->z ); 
  }

  for ( i = 0 ; i < 3; ++i) 
  {
    center[i] = center[i] / data->numatoms;
  }


  /* determine the the minimum and maximum {xyz}
   * coordintes of the system; having that we
   * can dimension the orbital grid */

  /* set initial values of temp values to the coordinates
   * of the first atom */
  xyzdim[0] = xyzdim[3] = temp_data->x;
  xyzdim[1] = xyzdim[4] = temp_data->y;
  xyzdim[2] = xyzdim[5] = temp_data->z;  

  /* now loop over the rest of the atoms to check if there's
   * something larger/smaller for the maximum and minimum
   * respectively */
  for( i = 0 ; i != data->numatoms; ++i)
  {
    if ( (temp_data+i)->x < xyzdim[0] ) 
       xyzdim[0] = (temp_data+i)->x;
    if ( (temp_data+i)->y < xyzdim[1] )
       xyzdim[1] = (temp_data+i)->y;
    if ( (temp_data+i)->z < xyzdim[2] ) 
       xyzdim[2] = (temp_data+i)->z;
    if ( (temp_data+i)->x > xyzdim[3] ) 
       xyzdim[3] = (temp_data+i)->x;
    if ( (temp_data+i)->y > xyzdim[4] ) 
       xyzdim[4] = (temp_data+i)->y;
    if ( (temp_data+i)->z > xyzdim[5] )
       xyzdim[5] = (temp_data+i)->z;
  }

  /* store final results in array */
  data->system_dimensions = xyzdim;

  /* done :) */
  return TRUE;
}


/**********************************************************
 *
 * this subroutine reads the number of procs and the amount
 * of memory requested
 *
 **********************************************************/
static int get_proc_mem(void *mydata) {

  gamessdata *data = (gamessdata *)mydata;
  char word[3][BUFSIZ];
  char buffer[BUFSIZ];
  unsigned int nproc;
  unsigned int i;


  /* initialize arrays */
  buffer[0] = '\0';
  for ( i = 0; i < 3; ++i) word[i][0] = '\0';


  /* first, we need to rewind the file */
  rewind(data->file);

  
  /* scan for the number of processors */
  do 
  {
    ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )
    sscanf(buffer,"%s %d %s",&word[0][0],&nproc,&word[1][0]);

  } while( strcmp(&word[0][0],"Initiating") || 
           strcmp(&word[1][0],"compute") );

  
  /* store the number of processors */
  data->nproc = nproc;

  
  /* scan for the amount of memory requested */
  do 
  {
    ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )
    sscanf(buffer,"%s %s",&word[0][0],&word[1][0]);

  } while( strcmp(&word[0][0],"$SYSTEM") || 
           strcmp(&word[1][0],"OPTIONS") );

   
  /* skip next line */
  eatline(data->file);


  /* next line contains the amount of memory requested */
  ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )
  sscanf(buffer,"%s %s %s",&word[0][0],&word[1][0],&word[2][0]);

  /* store it */
  strncpy(data->memory,&word[2][0],sizeof(data->memory));


  /* short message */
  printf("gamessplugin> GAMESS used %d compute processes \n",nproc);
  printf("gamessplugin> GAMESS used %s words of memory \n",
      data->memory);


  /* done here */
  return TRUE;
}



/**********************************************************
 *
 * this subroutine checks if the provided files is
 * actually a GAMESS file;
 *
 **********************************************************/
static int have_gamess(void *mydata) 
{
  gamessdata *data = (gamessdata *)mydata;
  char word[3][BUFSIZ];
  char buffer[BUFSIZ];
  int day, year;
  char month[BUFSIZ], rev[BUFSIZ];
  unsigned int i = 0;
 

  /* initialize arrays */
  buffer[0] = '\0';
  for ( i = 0; i < 3; ++i ) word[i][0] = '\0';


  /* check if the file is GAMESS format 
   * for now I just read line by line until 
   * the word GAMESS appears                */
  do
  {
    ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )
    sscanf(buffer,"%s %s %s",&word[0][0],&word[1][0],&word[2][0]);

  } while( strcmp(&word[1][0],"GAMESS") || 
           strcmp(&word[2][0],"VERSION") );


  /* extract the version number if possible; otherwise
   * return empty string */
  if ( strstr(buffer,"=") != NULL ) 
  {
    strncpy(data->version_string,strstr(buffer,"=")+2,16); 
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
  else
  { 
    data->version = 1;
  }


  /* short messsage to stdout */
  printf("gamessplugin> Detected GAMESS format :)\n");
  printf("gamessplugin> GAMESS version = %s \n", 
      data->version_string);

  return TRUE;
}



/**********************************************************
 *
 * this subroutine extracts the GBASIS
 *
 **********************************************************/
static int get_gbasis(void *mydata) {

  gamessdata *data = (gamessdata *)mydata;
  char word[4][BUFSIZ];
  char buffer[BUFSIZ];
  char diffuse[BUFSIZ];
  char polarization[BUFSIZ];
  unsigned int i = 0;


  /* initialize arrays */
  buffer[0] = '\0';
  diffuse[0] = '\0';
  polarization[0] = '\0';
  for ( i = 0; i < 4; ++i) word[i][0] = '\0';


  /* to be safe let's rewind the file */
  rewind(data->file);


  /* start scanning */
  do 
  {
    ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )
    sscanf(buffer,"%s %s",&word[0][0],&word[1][0]);

  } while((strcmp(&word[0][0],"BASIS")) || 
          (strcmp(&word[1][0],"OPTIONS")));


  /* skip next line */
  eatline(data->file);


  /* the first string in the current line contains the
   * GBASIS used; copy it over into the gbasis variable
   * of gamessdata */
  ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )
  sscanf(buffer,"%s %s %s",&word[0][0],&word[1][0],&word[2][0]);
 
  strncpy(data->gbasis,(&word[0][0])+7,sizeof(data->gbasis));

  printf("gamessplugin> Run was performed with GBASIS = %s \n",
          data->gbasis);


  /* in case we're using a pople style basis set, i.e. 
   * GBASIS=N311,N31,N21 or STO we also scan for the number 
   * of gaussians, as well as p,d,f and diffuse functions
   * and use this info to assemble a "basis set string" */
  if ( !strncmp(data->gbasis,"N311",sizeof(data->gbasis)) ||
       !strncmp(data->gbasis,"N31",sizeof(data->gbasis)) ||
       !strncmp(data->gbasis,"N21",sizeof(data->gbasis)) ||
       !strncmp(data->gbasis,"STO",sizeof(data->gbasis)) ) 
  {
    /* word[2] read previously should still contain the
     * number of gaussians used */
    data->ngauss = atoi(&word[2][0]);


    /* the next line gives us the d,f and diffuse sp
     * functions */
    ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )
    sscanf(buffer,"%s %d %s %d %s %s",&word[0][0],&data->ndfunc,
	&word[1][0],&data->nffunc,&word[2][0],&word[3][0]);

    /* convert GAMESS' .TRUE./.FALSE. for DIFFSP into 1/0 */
    if ( !strncmp(&word[3][0],"T",sizeof(&word[3][0])) ) 
	data->diffsp = TRUE;


    /* the next line gives us the p and diffuse s
     * functions */
    ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )
    sscanf(buffer,"%s %d %s %s",&word[0][0],&data->npfunc,
	&word[1][0],&word[2][0]);

    /* convert GAMESS' .TRUE./.FALSE. for DIFFSP into 1/0 */
    if ( !strncmp(&word[2][0],"T",sizeof(&word[3][0])) ) 
	data->diffs = TRUE;


    /* now we need some logic to assemble this info into
     * some nice looking string a la "6-31G*" */

    /* create the diffuse function string */
    if ( data->diffs && data->diffsp )
    {
    	strncpy(diffuse,"++",sizeof(diffuse));
    }
    else if ( data->diffsp ) 
    {
    	strncpy(diffuse,"+",sizeof(diffuse));
    }
    else 
    {
    	strncpy(diffuse,"",sizeof(diffuse));
    }


    /* create the polarization function string */
    if ( ( data->npfunc > 0 ) && ( data->ndfunc > 0) &&
	 ( data->nffunc > 0 )) 
    {
#if 1
        /* Windows doesn't provide snprintf in MSVC 6 */
	sprintf(polarization,
	    "(%dp,%dd,%df)", data->npfunc, data->ndfunc,
	    data->nffunc);
#else
	snprintf(polarization,sizeof(polarization),
	    "(%dp,%dd,%df)", data->npfunc, data->ndfunc,
	    data->nffunc);
#endif
    }
    else if ( ( data->npfunc > 0 ) && ( data->ndfunc > 0) ) 
    {
#if 1
        /* Windows doesn't provide snprintf in MSVC 6 */
	sprintf(polarization,
	    "(%dp,%dd)", data->npfunc, data->ndfunc);
#else
	snprintf(polarization,sizeof(polarization),
	    "(%dp,%dd)", data->npfunc, data->ndfunc);
#endif
    }
    else if ( ( data->npfunc > 0 ) )  
    {
#if 1
        /* Windows doesn't provide snprintf in MSVC 6 */
        sprintf(polarization,
            "(%dp)", data->npfunc);
#else
        snprintf(polarization,sizeof(polarization),
            "(%dp)", data->npfunc);
#endif
    }
    else if ( ( data->ndfunc > 0) ) 
    {
#if 1
        /* Windows doesn't provide snprintf in MSVC 6 */
        sprintf(polarization,
            "(%dd)", data->ndfunc);
#else
        snprintf(polarization,sizeof(polarization),
            "(%dd)", data->ndfunc);
#endif
    } 
    else 
    {
        strncpy(polarization,"",sizeof(polarization));
    } 
   

    /* assemble the bits */ 
    if ( !strcmp(data->gbasis,"STO") ) 
    {
#if 1
      /* Windows doesn't provide snprintf in MSVC 6 */
      sprintf(data->basis_string,
	  "STO-%dG%s%s", data->ngauss, diffuse, polarization);
#else
      snprintf(data->basis_string,sizeof(data->basis_string),
	  "STO-%dG%s%s", data->ngauss, diffuse, polarization);
#endif
     }
     else 
     {
#if 1
      /* Windows doesn't provide snprintf in MSVC 6 */
      sprintf(data->basis_string,
	  "%d-%s%sG%s", data->ngauss, (data->gbasis+1), diffuse, 
	  polarization);
#else
      snprintf(data->basis_string,sizeof(data->basis_string),
	  "%d-%s%sG%s", data->ngauss, (data->gbasis+1), diffuse, 
	  polarization);
#endif
     }      
  }

  /* for non pople style basis sets we just use the GBASIS
   * for the basis string;
   * TODO: make the basis_string more comprehensive for non
   *       pople-style basis sets */
  else 
  {
    strncpy(data->basis_string,data->gbasis,
	sizeof(data->basis_string));
  }


  /* short message */
  printf("gamessplugin> Run was performed with BASIS = %s \n",
      data->basis_string);


  return TRUE;
}



/**********************************************************
 *
 * this subroutine extracts run title 
 *
 **********************************************************/
static int get_runtitle(void *mydata) 
{
  gamessdata *data = (gamessdata *)mydata;
  char word[2][BUFSIZ];
  char buffer[BUFSIZ];


  /* initialize arrays */
  word[0][0] = '\0';
  word[1][0] = '\0';
  buffer[0] = '\0';

  
  /* look for RUN TITLE section */
  do 
  {
    ERR_FALSE( fgets(buffer,sizeof(buffer),data->file))
    sscanf(buffer,"%s %s",&word[0][0],&word[1][0]);

  } while( (strcmp(&word[0][0],"RUN")) || 
           (strcmp(&word[1][0],"TITLE")));


  /* the second to next line has what we want */
  eatline(data->file);
  ERR_FALSE( fgets(buffer,sizeof(buffer),data->file))
  strncpy(data->runtitle,chop_string_nl(buffer),sizeof(data->runtitle));

  return TRUE;
} 



/**********************************************************
 *
 * this subroutine checks the is we support wavefunction
 * and orbital stuff for the current gbasis
 *
 **********************************************************/
static int have_supported_gbasis(void *mydata) {

  gamessdata *data = (gamessdata *)mydata;

  /* check for supported gbasis, otherwise don't even attempt
   * to parse the orbital data/wavefunction */
  if ( !strcmp(data->gbasis,"MINI") ||
       !strcmp(data->gbasis,"MIDI") ||
       !strcmp(data->gbasis,"STO") ||
       !strcmp(data->gbasis,"N21") ||
       !strcmp(data->gbasis,"N31") ||
       !strcmp(data->gbasis,"N311") ||
       !strcmp(data->gbasis,"DZV") ||
       !strcmp(data->gbasis,"DH") ||
       !strcmp(data->gbasis,"TZV") ||
       !strcmp(data->gbasis,"MC") ) return TRUE;
       

  /* I guess we dont' support it then */
  return FALSE;
}



/**********************************************************
 *
 * this subroutine checks the RUNTYP
 *
 **********************************************************/
static int check_contrl(void *mydata) {

  gamessdata *data = (gamessdata *)mydata;
  char word[2][BUFSIZ];
  char buffer[BUFSIZ];


  /* initialize buffer */
  word[0][0] = '\0';
  word[1][0] = '\0';
  buffer[0] = '\0';


  /* start scanning; currently we support
   * RUNTYP = ENERGY, OPTIMIZE, SADPOINT, HESSIAN */
  do 
  {
    ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )
    sscanf(buffer,"%s %s",&word[0][0],&word[1][0]);

  } while( strcmp(&word[0][0],"\044CONTRL") || 
	   strcmp(&word[1][0],"OPTIONS") );

	   
  /* skip next line */
  eatline(data->file);


  /* current line contains RUNTYP info; scan it */
  ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )
  sscanf(buffer,"%s %s",&word[0][0],&word[1][0]);


  /* check for supported RUNTYPs */
  if (!strcmp(&word[1][0],"RUNTYP=ENERGY")) 
  {
    printf("gamessplugin> File generated via %s \n",&word[1][0]);
    data->runtyp = ENERGY;
    strncpy(data->runtyp_string,"Single point",
	    sizeof(data->runtyp_string));
  }
  else if (!strcmp(&word[1][0],"RUNTYP=OPTIMIZE")) 
  {
    printf("gamessplugin> File generated via %s \n",&word[1][0]);
    data->runtyp = OPTIMIZE;
    strncpy(data->runtyp_string,"Geometry optimization",
	    sizeof(data->runtyp_string));
  }
  else if (!strcmp(&word[1][0],"RUNTYP=SADPOINT")) 
  {
    printf("gamessplugin> File generated via %s \n",&word[1][0]);
    data->runtyp = SADPOINT;
    strncpy(data->runtyp_string,"Saddle point search",
	    sizeof(data->runtyp_string));
  }
  else if (!strcmp(&word[1][0],"RUNTYP=HESSIAN")) 
  {
    printf("gamessplugin> File generated via %s \n",&word[1][0]);
    data->runtyp = HESSIAN;
    strncpy(data->runtyp_string,"Hessian calculation",
	    sizeof(data->runtyp_string));
  }
  else 
  {
    printf("gamessplugin> The %s is currently not supported \n",
	    &word[1][0]);
    return FALSE;
  }

  
  /* reserve memory for scfenergy array */
  data->scfenergies = (double *)calloc(1,sizeof(double));


  /* determine SCFTYP */
  if (!strcmp(&word[0][0],"SCFTYP=RHF")) 
  {
    printf("gamessplugin> Type of wavefunction used %s \n",
	&word[0][0]);
    data->scftyp = RHF;
    strncpy(data->scftyp_string,"RHF",sizeof(data->scftyp_string));
  }
  else if (!strcmp(&word[0][0],"SCFTYP=UHF")) 
  {
    printf("gamessplugin> Type of wavefunction used %s \n",
	&word[0][0]);
    data->scftyp = UHF;
    strncpy(data->scftyp_string,"UHF",sizeof(data->scftyp_string));
  }
  else if (!strcmp(&word[0][0],"SCFTYP=ROHF")) 
  {
    printf("gamessplugin> Type of wavefunction used %s \n",
	&word[0][0]);
    data->scftyp = ROHF;
    strncpy(data->scftyp_string,"ROHF",sizeof(data->scftyp_string));
  }
  else if (!strcmp(&word[0][0],"SCFTYP=GVB")) 
  {
    printf("gamessplugin> Type of wavefunction used %s \n",
	&word[0][0]);
    data->scftyp = GVB;
    strncpy(data->scftyp_string,"GVB",sizeof(data->scftyp_string));
  }
  else if (!strcmp(&word[0][0],"SCFTYP=MCSCF")) 
  {
    printf("gamessplugin> Type of wavefunction used %s \n",
	&word[0][0]);
    data->scftyp = MCSCF;
    strncpy(data->scftyp_string,"MCSCF",sizeof(data->scftyp_string));
  }
  else 
  {
    /* if we don't find a supported SCFTYP we bomb out; this
     * might be a little drastic */
    printf("gamessplugin> %s is currently not supported \n",
	    &word[0][0]);
    strncpy(data->scftyp_string,"\0",sizeof(data->scftyp_string));
    return FALSE;
  }


  eatline(data->file);

  /* next we determine the coord type (i.e. zmatrix, unique ..)
   * as provided with the COORD keyword */
  ERR_FALSE( fgets(buffer,sizeof(buffer),data->file) )
  strcpy(data->geometry,
      chop_string_all(strstr(buffer,"COORD =")+7)); 


  /* brief message */
  printf("gamessplugin> Coordinate type %s used \n",data->geometry);


  /* at this point everything should be in good shape */
  return TRUE;
}




/********************************************************
 *
 * THIS IS THE END OF THE GAMESS SPECIFIC SUBROUTINES
 *
 ********************************************************/




/*********************************************************
 *
 * this function prints the current date and time 
 *
 *********************************************************/
void 
get_time(char *mytime)
{
  time_t current_time;
  struct tm *curtime;

  /* get the current time */

  time(& current_time);

  /* figure out the local time */

  curtime = localtime(& current_time);

  /* print it */

  (void) strftime(mytime,BUFSIZ, 
		  "programm start: %a, %0d %b %0H:%0M:%0S", 
		  curtime);
  
  return;
}



/********************************************************
 *
 * this function removes trailing spaces/newlines off
 * a character array
 *
 ********************************************************/
char* chop_string_all( char* the_string )
{
  int i = 0;

  while ( (*(the_string+i) != '\n') && (*(the_string+i) != ' ') && 
          (*(the_string+i) != '\0')) 
  {
    ++i;
  }

  *(the_string+i) = '\0';

  return the_string;
}



/********************************************************
 *
 * this function removes trailing newlines off
 * a character array
 *
 ********************************************************/
char* chop_string_nl( char* the_string )
{
  int i = 0;

  while ( (*(the_string+i) != '\n') && (*(the_string+i) != '\0')) 
  {
    ++i;
  }
  *(the_string+i) = '\0';

  return the_string;
}



/*********************************************************
 *
 * this function calculates the value of the wavefunction
 * corresponding to a particular orbital at grid point
 * grid_x, grid_y, grid_z
 *
 *********************************************************/
float orbital_at_grid_xyz(void *mydata, float *wave_function,
    float grid_size, float grid_x, float grid_y, float grid_z)
{
  gamessdata *data = (gamessdata *)mydata;
  gamess_temp *temp_data;
  int *basis_counter, position;
  float xtemp = 0.0, ytemp = 0.0, ztemp = 0.0;
  float xtemp2 = 0.0, ytemp2 = 0.0, ztemp2 = 0.0;
  char *orbital_symmetry;
  int at;
  int prim, shell;
  int orbital_pointer = 0;
  float dist2 = 0; 
  float value, value_grid;
  int *temp_shell;
  int *temp_prims;
  int sym_counter;
  float temp_p_value, temp_s_value, temp_d_value, temp_f_value;
  int have_s, have_p, have_d, have_f;
  float temp;

    
  /* initialize value of orbital at gridpoint */
  value_grid = 0.0;
  value = 0.0;


  /* retrieve orbital symmetry of shel primitives */
  orbital_symmetry = data->orbital_symmetry;

	
  /* retrieve array with # of primitives available
  * for each indivual atom of the system */
  basis_counter = data->basis_counter;


  /* retrieve the pointer for the temp arrays */
  temp_data = data->temporary;

  /* initialize the orbital pointer */
  orbital_pointer = 0; 


  /* initialize the primitive counter */
  position = 0;
  sym_counter = 0;


  /* store shells and primitive counter in
  * temporary arrays */
  temp_shell = data->atomic_shells;
  temp_prims = data->shell_primitives;

	
  /* loop over all the QM atoms */
  for ( at = 0 ; at < data->numatoms ; ++at) 
  {
    /* at this point we can calculate the
     * distance of the current grid point
     * from the center of atom at */
    xtemp =  ( grid_x - (temp_data+at)->x );
    ytemp =  ( grid_y - (temp_data+at)->y );
    ztemp =  ( grid_z - (temp_data+at)->z );

    xtemp2 =  pow(xtemp,2.0);
    ytemp2 =  pow(ytemp,2.0);
    ztemp2 =  pow(ztemp,2.0);

	  
    dist2 = xtemp2 + ytemp2 + ztemp2;

    /* loop over the primitives belonging to
     * this atom;
     * the number of shells per atom is encoded
     * in the array atomic_shells;
     * the number of primitives per shell is 
     * stored in the array shell_primitives */
    for ( shell = 0; shell < *temp_shell; ++shell)
    {

      /* initialize temp arrays */
      temp_p_value = 0.0; 
      temp_s_value = 0.0;
      temp_d_value = 0.0;
      temp_f_value = 0.0;


      /* flags keepting track of what orbital
       * symmetries we parsed in the current 
       * shell */
      have_s = FALSE;
      have_p = FALSE;
      have_d = FALSE;
      have_f = FALSE;


      /* start looping over the shell primitives */
      for ( prim = 0; prim < *temp_prims;  ++prim)
      {
	temp = *(data->basis+position+1) * 
	  exp(*(data->basis + position)*-1.0*dist2);

	
	/* compute-routine for S-shells */
	if ( *(orbital_symmetry+sym_counter) == 'S' )
	{
	  /* indicate presence of S shell primitives */ 
	  have_s = TRUE;

	  /* add primitive to s shell */
	  temp_s_value = temp_s_value + temp;

	  /* increase the counter keeping track of the total
	   * number of primitives */
	  position += 2;
	}
	   
	/* compute-routine for L-shells */
	else if ( *(orbital_symmetry+sym_counter) == 'P' )
	{
	  /* indicate presence of P shell primitives */
	  have_p = TRUE;
	       
	  /* add primitive to p shell */
	  temp_p_value = temp_p_value + temp;


	  /* increase the counter keeping track of the total
	   * number of primitives */
	  position += 2;
	} 
	       
	        
	/* compute-routine for D-shells */
	else if ( *(orbital_symmetry+sym_counter) == 'D' )
	{
	  /* indicate presence of P shell primitives */
	  have_d = TRUE;
	       
	  /* add primitive to p shell */
	  temp_d_value = temp_d_value + temp;


	  /* increase the counter keeping track of the total
           * number of primitives */
	  position += 2;
	} 

	       
	/* compute-routine for F-shells */
	else if ( *(orbital_symmetry+sym_counter) == 'F' )
	{
	  /* indicate presence of P shell primitives */
	  have_f = TRUE;
	       
	  /* add primitive to p shell */
	  temp_f_value = temp_f_value + temp;


	  /* increase the counter keeping track of the total
	   * number of primitives */
	  position += 2;
	}


	/* at this point we have fallen through all the
      	 * shell types we know and we better quit */
	else
	{
	  printf("gamessplugin> WARNING ... ");
	  printf("Encountered unknown shell type %d \n", 
	    *(orbital_symmetry+sym_counter));
	  return FALSE;
	}


	/* update counter that keeps track of the orbital
	 * symmetry of the current shell */
	sym_counter++;	 
      }

	
      /* multiply with the appropriate wavefunction
      * coefficient */

      /* S shells */
      if ( have_s ) 
      {
	/* s */
	value = *(wave_function + orbital_pointer) * temp_s_value;
	orbital_pointer++;
      }
	  

      /* p shells */
      if ( have_p )
      {
	/* px */
	value = *(wave_function + orbital_pointer) * temp_p_value * 
	    xtemp;
	orbital_pointer++;

	/* py */
	value = *(wave_function + orbital_pointer) * temp_p_value * 
	    ytemp;
	orbital_pointer++;
	     
	/* pz */ 
	value = *(wave_function + orbital_pointer) * temp_p_value * 
	    ztemp;
	orbital_pointer++;
      }

            
      /* d shells */
      if ( have_d )
      {
	/* d_xx */
	value = *(wave_function + orbital_pointer) * temp_d_value * 
	    xtemp * xtemp;
	orbital_pointer++;

	/* d_yy */
	value = *(wave_function + orbital_pointer) * temp_d_value * 
	    ytemp * ytemp;
	orbital_pointer++;

	/* d_zz */
	value = *(wave_function + orbital_pointer) * temp_d_value * 
	    ztemp * ztemp;
	orbital_pointer++;

	/* d_xy */
	value = *(wave_function + orbital_pointer) * temp_d_value * 
	    xtemp * ytemp;
	orbital_pointer++;

	/* d_xz */
	value = *(wave_function + orbital_pointer) * temp_d_value * 
	    xtemp * ztemp;
	orbital_pointer++;

	/* d_yz */
	value = *(wave_function + orbital_pointer) * temp_d_value * 
	    ytemp * ztemp;
	orbital_pointer++;
      }


      /* f shells */
      if ( have_f )
      {
	/* f_xxx */
	value = *(wave_function + orbital_pointer) * temp_d_value * 
	    xtemp * xtemp * xtemp;
	orbital_pointer++;

	/* f_yyy */
	value = *(wave_function + orbital_pointer) * temp_d_value * 
	  ytemp * ytemp * ytemp;
	orbital_pointer++;

	/* f_zzz */
	value = *(wave_function + orbital_pointer) * temp_d_value * 
	    ztemp * ztemp * ztemp;
	orbital_pointer++;

	/* f_xxy */
	value = *(wave_function + orbital_pointer) * temp_d_value * 
	    xtemp * xtemp * ytemp;
	orbital_pointer++;

	/* f_xxz */
	value = *(wave_function + orbital_pointer) * temp_d_value * 
	    xtemp * xtemp * ztemp;
	orbital_pointer++;

	/* f_yyx */
	value = *(wave_function + orbital_pointer) * temp_d_value * 
	    ytemp * ytemp * xtemp;
	orbital_pointer++;

	/* f_yyz */
	value = *(wave_function + orbital_pointer) * temp_d_value * 
	    ytemp * ytemp * ztemp;
	orbital_pointer++;

	/* f_zzx */
	value = *(wave_function + orbital_pointer) * temp_d_value * 
	    ztemp * ztemp * xtemp;
	orbital_pointer++;

	/* f_zzy */
	value = *(wave_function + orbital_pointer) * temp_d_value * 
	    ztemp * ztemp * ytemp;
	orbital_pointer++;

	/* f_xyz */
	value = *(wave_function + orbital_pointer) * temp_d_value * 
	    xtemp * ytemp * ztemp;
	orbital_pointer++;
      }


      /* move primitives pointer to next shell */
      ++temp_prims;


      /* store value in temporary array, and 
      * reset value */
      value_grid = value_grid + value;
      value = 0.0; 
    } 


    /* move shell pointer to next atom */
    temp_shell++;
  }

  /* return value at grid point */
  return value_grid;
}


