/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_mol2plugin
#define STATIC_PLUGIN 1

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
 *      $RCSfile: mol2plugin.C,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.20 $       $Date: 2006/03/30 02:41:28 $
 *
 ***************************************************************************/

/*
 * mol2 file reader
 * More information on this format is available at
 *   http://www.tripos.com/data/support/mol2.pdf
 *   http://www.tripos.com/mol2/
 *
 *   DOCK mol2 page: 
 *     http://www.csb.yale.edu/userguides/datamanip/dock/DOCK_4.0.1/html/Manual.41.html
 *
 * This plugin currently reads the following record types:
 *  MOLECULE
 *  ATOM
 *  BOND
 *
 */

#include "molfile_plugin.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#if defined(_AIX)
#include <strings.h>
#endif

#define LINESIZE 256

typedef struct {
  FILE *file;
  molfile_atom_t *atomlist;
  int natoms, nbonds, optflags, coords_read;
  int *from, *to;
  float *bondorder;
} mol2data;

// Open the file and create the mol2 struct used to pass data to the other
// functions.
static void *open_mol2_read(const char *path, const char *filetype, 
    int *natoms) {
  FILE *fd;
  mol2data *mol2;
  char line[LINESIZE]; 
  int match, nbonds, optflags;

  fd = fopen(path, "r");
  if (!fd)
    return NULL;
  
  // Find and read the MOLECULE record
  do {
    fgets(line, LINESIZE, fd);
    if ( ferror(fd) || feof(fd) ) {
      fprintf(stderr, "mol2plugin: No molecule record found in file.\n");
      return NULL;
    }
  } while ( strncmp(line, "@<TRIPOS>MOLECULE", 17) );

  fgets(line, LINESIZE, fd);  // Read and ignore the mol_name
  fgets(line, LINESIZE, fd);  // Read the molecule info
  match = sscanf(line, " %d %d", natoms, &nbonds);
  if (match == 1) {
    nbonds = 0;
  }
  else if (match != 2) {
    fprintf(stderr, "mol2plugin: Cannot determine the number of atoms.\n");
    return NULL;
  }
  fgets(line, LINESIZE, fd);  // Read and ignore the mol_type
  fgets(line, LINESIZE, fd);  // Read the charge_type
  if ( strncmp(line, "NO_CHARGES", 10) == 0 ) {
    optflags = MOLFILE_NOOPTIONS;
  }
  else {
    optflags = MOLFILE_CHARGE;
  }

  // Allocate and initialize the mol2 structure
  mol2 = new mol2data;
  mol2->file = fd;
  mol2->natoms = *natoms;
  mol2->nbonds = nbonds;
  mol2->optflags = optflags;
  mol2->coords_read = 0;
  mol2->from = NULL;
  mol2->to = NULL;
  mol2->bondorder = NULL;

  return mol2;
}

// Read atom information, but not coordinates.
static int read_mol2(void *v, int *optflags, molfile_atom_t *atoms) {
  mol2data *mol2 = (mol2data *)v;
  char line[LINESIZE]; 
  int i, match;
  molfile_atom_t *atom;

  *optflags = mol2->optflags;

  // Find and read the ATOM record
  rewind(mol2->file);
  do {
    fgets(line, LINESIZE, mol2->file);
    if ( ferror(mol2->file) || feof(mol2->file) ) {
      fprintf(stderr, "mol2plugin: No atom record found in file.\n");
      return MOLFILE_ERROR;
    }
  } while ( strncmp(line, "@<TRIPOS>ATOM", 13) );

  // Read the atoms
  for (i = 0; i < mol2->natoms; i++) {
    atom = atoms+i;

    fgets(line, LINESIZE, mol2->file);
    if ( ferror(mol2->file) || feof(mol2->file) ) {
      fprintf(stderr, "mol2plugin: Error occurred reading atom record.\n");
      return MOLFILE_ERROR;
    }

    match = sscanf(line, " %*d %s %*f %*f %*f %*s %d %s %f", 
      atom->name, &atom->resid, atom->resname, &atom->charge);

    // The last three records are optional for mol2 files, supply values if
    // any are missing. Note that these cases are meant to fall through.
    switch (match) {
      case 0: 
        fprintf(stderr, "mol2plugin: Improperly formatted atom record.\n");
        return MOLFILE_ERROR;

      case 1:
        atom->resid = 0;

      case 2:
        sprintf(atom->resname, "%d", atom->resid);

      case 3:
        atom->charge = 0.0;

      default:
        break;
    }

    strcpy(atom->type, atom->name);
    // Leave these blank when not provided by the file.
    atom->chain[0] = '\0';
    atom->segid[0] = '\0';
  }

  return MOLFILE_SUCCESS;
}



// Create arrays of one-based bond indicies.
static int read_mol2_bonds(void *v, int *nbonds, int **fromptr, int **toptr, float **bondorderptr) {
  mol2data *mol2 = (mol2data *)v;
  char line[LINESIZE], bond_type[16]; 
  int i, match, bond_from, bond_to, bond_index, current_nbonds;
  float curr_order;

  if (mol2->nbonds == 0) {
    *nbonds = 0;
    *fromptr = NULL;
    *toptr = NULL;
    return MOLFILE_SUCCESS;
  }

  // Allocate memory for the from and to arrays. This will be freed in
  // close_mol2_read
  mol2->from = new int[mol2->nbonds];
  mol2->to = new int[mol2->nbonds];
  mol2->bondorder=new float[mol2->nbonds];
  current_nbonds = mol2->nbonds;

  // Find and read the BOND record
  rewind(mol2->file);
  do {
    fgets(line, LINESIZE, mol2->file);
    if ( ferror(mol2->file) || feof(mol2->file) ) {
      fprintf(stderr, "mol2plugin: No atom record found in file.\n");
      return MOLFILE_ERROR;
    }
  } while ( strncmp(line, "@<TRIPOS>BOND", 13) );

  // Read the bonds
  bond_index = 0;
  for (i = 0; i < mol2->nbonds; i++) {
    fgets(line, LINESIZE, mol2->file);
    if ( ferror(mol2->file) || feof(mol2->file) ) {
      fprintf(stderr, "mol2plugin: Error occurred reading atom record.\n");
      return MOLFILE_ERROR;
    }

    //Move on if the next line is a header
    if (strncmp(line, "@", 1) == 0) {
      //Then the bonds are over
      break;
    }

    match = sscanf(line, " %*d %d %d %s", &bond_from, &bond_to, bond_type);
    if (match < 3) {
      fprintf(stderr, "mol2plugin: Improperly formatted bond record.\n");
      continue;
    }
    if ( strncmp(bond_type, "nc", 2) == 0 ) {
      // Not an actual bond, don't add it to the list
      current_nbonds--;
    }
    else {
      // Add the bond to the list
      curr_order=strtod(bond_type,NULL);
      if (curr_order<1.0 || curr_order>4.0) curr_order=1;
//      fprintf(stdout,"mol2plugin: Bond from %d to %d of order %f\n", bond_from, bond_to, curr_order);
      mol2->from[bond_index] = bond_from;
      mol2->to[bond_index] = bond_to;
      mol2->bondorder[bond_index]=curr_order;
      bond_index++;
    }
  }
  if (bond_index > 0) {
    *nbonds = current_nbonds;
    *fromptr = mol2->from;
    *toptr = mol2->to;
    *bondorderptr = mol2->bondorder; 
  } else {
    printf("mol2plugin) WARNING: no bonds defined in mol2 file\n");
    *nbonds = 0;
    *fromptr = NULL;
    *toptr = NULL;
    *bondorderptr = NULL; 
  }
    
//printf("End of read_mol2_bonds\n");
  return MOLFILE_SUCCESS;
}


// Read atom coordinates
static int read_mol2_timestep(void *v, int natoms, molfile_timestep_t *ts) {
  mol2data *mol2 = (mol2data *)v;
  char line[LINESIZE];
  int i, match;
  float x, y, z;

  // Since the file is rewound when coordinates are read, EOF shouldn't
  // happen. Instead, use a flag to indicate that the single timestep has
  // been read
  if (mol2->coords_read) {
    return MOLFILE_EOF;
  }

  // Find and read the ATOM record
  rewind(mol2->file);
  do {
    fgets(line, LINESIZE, mol2->file);
    if ( ferror(mol2->file) || feof(mol2->file) ) {
      fprintf(stderr, "mol2plugin: No atom record found in file.\n");
      return MOLFILE_ERROR;
    }
  } while ( strncmp(line, "@<TRIPOS>ATOM", 13) );

  // Read the atoms
  for (i = 0; i < mol2->natoms; i++) {
    fgets(line, LINESIZE, mol2->file);
    if ( ferror(mol2->file) || feof(mol2->file) ) {
      fprintf(stderr, "mol2plugin: Error occurred reading atom coordinates.\n");
      return MOLFILE_ERROR;
    }


    match = sscanf(line, " %*d %*s %f %f %f", &x, &y, &z);
    if (match < 3) {
      fprintf(stderr, "mol2plugin: Improperly formatted atom coordinates.\n");
      return MOLFILE_ERROR;
    }

    if (ts) {
      ts->coords[3*i  ] = x;
      ts->coords[3*i+1] = y;
      ts->coords[3*i+2] = z;
    }
  }

  mol2->coords_read = 1;
  return MOLFILE_SUCCESS;
}


static void *open_mol2_write(const char *filename, const char *filetype, 
                           int natoms) {
  FILE *fd;
  mol2data *data;

  fd = fopen(filename, "w");
  if (!fd) { 
    fprintf(stderr, "Error) Unable to open mol2 file %s for writing\n",
            filename);
    return NULL;
  }
  
  data = (mol2data *)malloc(sizeof(mol2data));
  data->natoms = natoms;
  data->file = fd;
//  data->file_name = strdup(filename);
  return data;
}


static int write_mol2_structure(void *mydata, int optflags, 
                               const molfile_atom_t *atoms) {
  mol2data *data = (mol2data *)mydata;
  data->atomlist = (molfile_atom_t *)malloc(data->natoms*sizeof(molfile_atom_t));
  memcpy(data->atomlist, atoms, data->natoms*sizeof(molfile_atom_t));
  return MOLFILE_SUCCESS;
}

void getmol2ff(char* outputtype, const char* psftype) {
//fprintf(stdout,"Doing ff typing on %s\n",psftype);
  if (strncmp(psftype,"H",1)==0) {
    //It's a hydrogen
    strncpy(outputtype, "H   ",4);
    return;
  } else if (strncmp(psftype,"C",1)==0) {
    //It's a carbon... probably
    if (strncmp(psftype,"C ",2)==0 || strncmp(psftype,"CA ",3)==0 || strncmp(psftype,"CPH",3)==0 || strncmp(psftype,"CPT",3)==0 || strncmp(psftype,"CC ",3)==0 || strncmp(psftype,"CD ",3)==0 || strncmp(psftype,"CN1",3)==0 || strncmp(psftype,"CN2",3)==0 || strncmp(psftype,"CN3",3)==0 || strncmp(psftype,"CN4",3)==0 || strncmp(psftype,"CN5",3)==0 || strncmp(psftype,"CNA",3)==0) {
	  strncpy(outputtype, "C.2 ",4);
	  return;
    } else {
	  strncpy(outputtype, "C.3 ",4);
	  return;
    }  
  } else if (strncmp(psftype,"N",1)==0) {
     //It"s probably nitrogen
     if (strncmp(psftype,"NR",2)==0 || strncmp(psftype,"NH1",3)==0 || strncmp(psftype,"NH2",3)==0 || strncmp(psftype,"NC2",3)==0 || strncmp(psftype,"NY",2)==0 || (strncmp(psftype,"NN",2)==0 && strncmp(psftype,"NN6",3)!=0)) {
       strncpy(outputtype, "N.am",4);
       return;
       } else {
       strncpy(outputtype, "N.3 ",4);
       return;
       }
  } else if (strncmp(psftype,"O",1)==0) {
     //Probably an oxygen
     if (strncmp(psftype,"OH1",3)==0 || strncmp(psftype,"OS",2)==0 || strncmp(psftype,"OT ",3)==0 || strncmp(psftype,"ON4",3)==0 || strncmp(psftype,"ON5",3)==0 || strncmp(psftype,"ON6",3)==0) {
        strncpy(outputtype, "O.3 ",4);
	return;
     } else {
        strncpy(outputtype, "O.2 ",4);
	return;
     } 
  } else if (strncmp(psftype,"S",1)==0) {
     strncpy(outputtype, "S.3 ",4);
     return;
  } else if (strncmp(psftype,"P",1)==0) {
     strncpy(outputtype, "P.3 ",4);
     return;
  } else {
     strncpy(outputtype, "X.  ",4);
     return;
  }
}






static int write_mol2_timestep(void *mydata, const molfile_timestep_t *ts) {
  mol2data *data = (mol2data *)mydata; 
  const molfile_atom_t *atom;
  const float *pos;
  float chrgsq;
  int i;

  // try to guess whether we have charge information.
  chrgsq=0.0;
  atom = data->atomlist;
  for (i = 0; i < data->natoms; i++) {
      chrgsq += atom->charge*atom->charge;
      ++atom;
  }

  //print header block
  fprintf(data->file, "@<TRIPOS>MOLECULE\n");
  fprintf(data->file, "generated by VMD\n");
  fprintf(data->file, "%4d %4d 1 0 0\n", data->natoms, data->nbonds);
  fprintf(data->file, "SMALL\n");
  // educated guess
  if (chrgsq > 0.0001) {
      fprintf(data->file, "USER_CHARGES\n");
  } else {
      fprintf(data->file, "NO_CHARGES\n");
  }
  fprintf(data->file, "****\n");
  fprintf(data->file, "Energy = 0\n\n");
  
  //print atoms block
  fprintf(data->file, "@<TRIPOS>ATOM\n");
  atom = data->atomlist;
  pos = ts->coords;
  char mol2fftype[5];
  for (i = 0; i < data->natoms; i++) {
    getmol2ff(mol2fftype, atom->type);
    fprintf(data->file, "%7d %-4s      %8.4f  %8.4f  %8.4f %4s %4d  %3s        %8.6f\n",
            i+1, atom->name, pos[0], pos[1], pos[2], mol2fftype, 
            atom->resid, atom->resname, atom->charge);
    ++atom; 
    pos += 3;
  }

  //print bond info

  int l=1; //number of bond record
  printf("mol2plugin) numbonds: %d\n", data->nbonds);
  if (data->nbonds>0) fprintf(data->file, "@<TRIPOS>BOND\n");
  for (i=0; i<data->nbonds; i++) {
    // For mol2, only write bonds for fromptr[i]<toptr[i]
    // bondorder is either 1, 2, 3 or a textual representation: am,ar,du,un,nc
    // we don't have the info for the text, so we truncate to integer.
    fprintf(data->file, "%5d %5d %5d %2d\n", l ,data->from[i], data->to[i],
            (int)data->bondorder[i]);
    l++;
  } 

  // Print out substructure info to keep some programs sane
  fprintf(data->file,"\n@<TRIPOS>SUBSTRUCTURE\n");
  fprintf(data->file,"1 ****        1 TEMP                        ");
  fprintf(data->file,"0 ****  **** 0 ROOT\n");

  return MOLFILE_SUCCESS;
}

static int write_bonds(void *v, int nbonds, int *fromptr, int *toptr, float *bondorderptr) {
  mol2data *data = (mol2data *)v;
  data->from = new int[nbonds];
  data->to = new int[nbonds];
  data->bondorder = new float[nbonds];
  //set the pointers for use later
  for (int i=0;i<nbonds;i++) {
	  data->from[i]=fromptr[i];
	  data->to[i]=toptr[i];
	  data->bondorder[i]=bondorderptr[i];
  }
  data->nbonds = nbonds;
  return MOLFILE_SUCCESS;
}


static void close_mol2_write(void *mydata) {
  mol2data *data = (mol2data *)mydata;
  if (data) {
    if (data->file) fclose(data->file);
    if (data->from != NULL) free(data->from);
    if (data->to != NULL)   free(data->to);
    if (data->bondorder != NULL)   free(data->bondorder);
    delete data;
  }
}

//
// Free the memory used by the mol2 structure
static void close_mol2_read(void *v) {
  mol2data *mol2 = (mol2data *)v;
  if (mol2) {
    if (mol2->file) fclose(mol2->file);
    if (mol2->from != NULL) free(mol2->from);
    if (mol2->to != NULL)   free(mol2->to);
    if (mol2->bondorder != NULL)   free(mol2->bondorder);
    delete mol2;
  }
}


static int read_bonds(void *v, int *nbonds, int **fromptr, int **toptr, float **bondorderptr) {
  mol2data *mol2 = (mol2data *)v;

  /* now read bond data */
//  *nbonds = start_psf_bonds(psf->fp);

  if (mol2->nbonds > 0) {
    mol2->from = (int *) malloc(*nbonds*sizeof(int));
    mol2->to = (int *) malloc(*nbonds*sizeof(int));
    mol2->bondorder = (float *) malloc(*nbonds*sizeof(float));
    if ((read_mol2_bonds(mol2, nbonds, &(mol2->from), &(mol2->to), &(mol2->bondorder))) != MOLFILE_SUCCESS) {
      fclose(mol2->file);
      mol2->file = NULL;
      return MOLFILE_ERROR;
    }
    *fromptr = mol2->from;
    *toptr = mol2->to;
    *bondorderptr = mol2->bondorder; 
  } else {
    printf("mol2plugin) WARNING: no bonds defined in mol2 file.\n");
    *fromptr=NULL;
    *toptr=NULL;
    *bondorderptr=NULL;
  }
  return MOLFILE_SUCCESS;
}

static molfile_plugin_t mol2plugin = {
  vmdplugin_ABIVERSION,
  MOLFILE_PLUGIN_TYPE,                         /* type */
  "mol2",                                      /* short name */
  "MDL mol2",                                  /* pretty name */
  "Peter Freddolino and Eamon Caddigan",       /* author */
  0,                                           /* major version */
  8,                                           /* minor version */
  VMDPLUGIN_THREADSAFE,                        /* is reentrant */
  "mol2",
  open_mol2_read,		/*open file read*/
  read_mol2,			/*read structure*/
  read_bonds,			/* read bond list*/
  read_mol2_timestep,		/* read next timestep*/
  close_mol2_read,		/* close_file_read */
  open_mol2_write,		/* open file for writing */
  write_mol2_structure,		/* write structure */
  write_mol2_timestep,		/* write tempestep*/
  close_mol2_write,		/*close file for writing*/
  0,                            /* read_volumetric_metadata */
  0,                            /* read_volumetric_data */
  0,                            /* read_rawgraphics */
  0,				/*read molecule metadata */
  write_bonds			/* write bonds */
};

VMDPLUGIN_EXTERN int VMDPLUGIN_init() {
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_EXTERN int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)&mol2plugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_EXTERN int VMDPLUGIN_fini() {
  return VMDPLUGIN_SUCCESS;
}

