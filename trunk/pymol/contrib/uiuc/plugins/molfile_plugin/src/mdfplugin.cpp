/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_mdfplugin
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
 *      $RCSfile: mdfplugin.C,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.15 $       $Date: 2009/04/29 15:45:31 $
 *
 ***************************************************************************/

/*
 * Molecular data file (.mdf) reader
 * Insight II, Discover, etc. structure and bond information. This plugin
 * reads only the topology section, ignoring the optional symmertry and
 * atomset sections.
 *
 * Format specification can be found at:
 * http://instinct.v24.uthscsa.edu/~hincklab/html/soft_packs/msi_docs/insight980/formats980/File_Formats_1998.html#484257
 *
 * TODO: The current code reads the file *four* times -- once on open, once
 * to read the structure, and twice to read the bonds. Perhaps these could
 * be consolidated, e.g. by counting the bonds and populating the hash
 * tables during open or read_structure.
 *
 */

#include "molfile_plugin.h"

#define VMDPLUGIN_STATIC
#include "hash.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#if defined(_AIX)
#include <strings.h>
#endif

#define LINESIZE 256
#define NAMESIZE 32

typedef struct {
  FILE *file;
  int natoms, nmols, *from, *to;
  long mol_data_location;
} mdfdata;

// Read a line of atom data and store the values in the atom structure
// Return 1 on success, 0 on error
static int read_mdf_structure_line(molfile_atom_t *atom, const char *line) {
  // Read pertinent structure information from the line
  if ( sscanf(line, "%[^:]:%s %s %*s %*s %*d %*s %f %*d %*d %*d %f",
              atom->resname, atom->name, atom->type, 
              &atom->charge, &atom->occupancy) != 5 ) {
    return 0;
  }

  // Get the resid from the resname
  if ( sscanf(atom->resname, "%*[^_]_%d", &atom->resid) != 1 ) {
    return 0;
  }

  // Provide defaults for missing values
  atom->chain[0] = '\0';
  atom->segid[0] = '\0';

  return 1;
}

// Read the atom info from src and copy the connectivity record to dest.
// Convert each record to resname_resnumber:atom form
// Return 1 on success, 0 on error
static int get_mdf_bonds(char *dest, const char *src) {
  char resinfo[NAMESIZE], bond_records[LINESIZE], *curr, *next, *tmp;

  // Get the connectivity records
  if ( sscanf(src, "%[^:]:%*s %*s %*s %*s %*d %*s %*f %*d %*d %*d %*f %*f %256c", resinfo, bond_records) != 2 ) {
    return 0;
  }

  // Append the bonds to the destination string, converting then to the
  // correct format along the way.
  dest[0] = '\0';
  for ( curr = bond_records; (next = strchr(curr, ' ')) != NULL;
        curr = next + 1 ) {
    *next = '\0';

    // Prepend the resname and resid to the destination atom name if it's
    // not already present.
    if ( strchr(curr, ':') == NULL ) {
      strcat(dest, resinfo);
      strcat(dest, ":");
    }

    // Remove cell/sympop/bondorder information from the bond
    if ( ((tmp = strchr(curr, '%')) != NULL) ||
         ((tmp = strchr(curr, '#')) != NULL) ||
         ((tmp = strchr(curr, '/')) != NULL) ||
         ((tmp = strchr(curr, '\n')) != NULL) ) {
      *tmp = '\0';
    }
    strcat(dest, curr);
    strcat(dest, " ");
  }

  return 1;
}

// Return the number of bond records on a line
static int count_mdf_bonds(const char *line) {
  char bond_records[LINESIZE];
  int bonds = 0;
  char *tmp;

  if ( !get_mdf_bonds(bond_records, line) ) {
    return -1;
  }
  
  for ( tmp = bond_records; (tmp = strchr(tmp, ' ')) != NULL;
        tmp++ ) {
    bonds++;
  }

  return bonds;
}

// Open the file and create the mdf struct used to pass data to the other
// functions.
static void *open_mdf_read(const char *path, const char *filetype, 
    int *natoms) {
  FILE *fd;
  mdfdata *mdf;
  long mol_data_location;
  char line[LINESIZE]; 
  int nmols = 0;

  fd = fopen(path, "r");
  if (!fd)
    return NULL;
  
  // Find the first molecule record
  do {
    fgets(line, LINESIZE, fd);
    if ( ferror(fd) || feof(fd) ) {
      fprintf(stderr, "mdfplugin) No molecule record found in file.\n");
      return NULL;
    }
  } while ( strncmp(line, "@molecule", 9) );

  // Remember the location of the beginning of the molecule data
  mol_data_location = ftell(fd);

  // Count the atoms in each molecule
  while ( line[0] != '#' ) {
    fgets(line, LINESIZE, fd);

    // Count atoms until a new molecule or the end of the section is reached
    while ( (line[0] != '@') && (line[0] != '#') ) {
      // Ignore blank and comment lines
      if ( !isspace(line[0]) && (line[0] != '!') )
        *natoms = *natoms + 1;
      fgets(line, LINESIZE, fd);
      if ( ferror(fd) || feof(fd) ) {
        fprintf(stderr, "mdfplugin) Error while counting atoms.\n");
        return NULL;
      }
    }
    nmols++;
  }

  // Allocate and initialize the mdf structure
  mdf = new mdfdata;
  mdf->file = fd;
  mdf->natoms = *natoms;
  mdf->nmols = nmols;
  mdf->from = NULL;
  mdf->to = NULL;
  mdf->mol_data_location = mol_data_location; 

  return mdf;
}

// Read the atom information for each molecule, but not bonds.
// XXX - this ignores the column records, which may cause the atom records
// to be read incorrectly.
static int read_mdf_structure(void *v, int *optflags, molfile_atom_t *atoms) {
  mdfdata *mdf = (mdfdata *)v;
  char line[LINESIZE];
  int mol_num;
  molfile_atom_t *atom = atoms;

  *optflags = MOLFILE_OCCUPANCY | MOLFILE_CHARGE;

  // Seek to the first molecule record
  fseek(mdf->file, mdf->mol_data_location, SEEK_SET);
  line[0] = '\0';

  // Read the atom structure for each molecule
  mol_num = 0;
  while ( line[0] != '#' ) {
    fgets(line, LINESIZE, mdf->file);

    // Read atom structure for the current molecule
    while ( (line[0] != '@') && (line[0] != '#') ) {
      // Ignore blank and comment lines
      if ( !isspace(line[0]) && (line[0] != '!') ) {
        if ( !read_mdf_structure_line(atom, line) ) {
          fprintf(stderr, "mdfplugin) Improperly formatted atom record encountered while reading structure.\n");
          return MOLFILE_ERROR;
        }

        // XXX - use the chain name to identify different molecules
        sprintf(atom->chain, "%d", mol_num);

        atom++;
      }

      fgets(line, LINESIZE, mdf->file);
      if ( ferror(mdf->file) || feof(mdf->file) ) {
        fprintf(stderr, "mdfplugin) File error while reading structure.\n");
        return MOLFILE_ERROR;
      }
    }
    mol_num++;
  }

  return MOLFILE_SUCCESS;
}

// Create arrays of one-based bond indicies.
static int read_mdf_bonds(void *v, int *nbonds, int **from_data, int **to_data, 
                          float **bondorderptr, int **bondtype, 
                          int *nbondtypes, char ***bondtypename) {
  mdfdata *mdf = (mdfdata *)v;
  int mol, atom, bond_count, *fromptr, *toptr, tmp_to;
  char *curr, *next, line[LINESIZE], bond_records[LINESIZE];
  char (*atomnames)[NAMESIZE]; // Dynamic array of cstrings
  hash_t *hasharray;           // Array of hash tables

  // Allocate and initialize the hash table for each molecule.
  hasharray = new hash_t[mdf->nmols];
  for (mol = 0; mol < mdf->nmols; mol++) {
    hash_init(&hasharray[mol], 256);
  }
  atomnames = new char[mdf->natoms][NAMESIZE];

  // Populate the hash table; key: atom name; value: one-based atom index.
  // Count the bonds, each bond is counted twice.
  fseek(mdf->file, mdf->mol_data_location, SEEK_SET);
  line[0] = '\0';
  atom = 1;
  mol = 0;
  bond_count = 0;
  while ( line[0] != '#' ) {
    fgets(line, LINESIZE, mdf->file);

    // Read the atom names
    while ( (line[0] != '@') && (line[0] != '#') ) {
      // Ignore blank and comment lines
      if ( !isspace(line[0]) && (line[0] != '!') ) {
        if ( sscanf(line, "%s %*s", atomnames[atom-1]) != 1 ) {
          fprintf(stderr, "mdfplugin) Improperly formatted atom record encountered while reading bonds.\n");
          return MOLFILE_ERROR;
        }
        if ( hash_insert(&hasharray[mol], atomnames[atom-1], atom) != HASH_FAIL ) {
          fprintf(stderr, "mdfplugin) Could not add atom to hash table.\n");
          return MOLFILE_ERROR;
        }

        bond_count += count_mdf_bonds(line);
        atom++;
      }

      fgets(line, LINESIZE, mdf->file);
      if ( ferror(mdf->file) || feof(mdf->file) ) {
        fprintf(stderr, "mdfplugin) File error while reading bonds.\n");
        return MOLFILE_ERROR;
      }
    }

    mol++;
  }

  bond_count /= 2;
  mdf->from = new int[bond_count];
  mdf->to = new int[bond_count];
  fromptr = mdf->from;
  toptr = mdf->to;

  // Read the molecules, storing the bond-indicies in fromptr and toprt
  fseek(mdf->file, mdf->mol_data_location, SEEK_SET);
  line[0] = '\0';
  atom = 1;
  mol = 0;
  while ( line[0] != '#' ) {
    fgets(line, LINESIZE, mdf->file);

    // Read the bonds
    while ( (line[0] != '@') && (line[0] != '#') ) {
      // Ignore blank and comment lines
      if ( !isspace(line[0]) && (line[0] != '!') ) {
        if ( !get_mdf_bonds(bond_records, line) ) {
          fprintf(stderr, "mdfplugin) Error reading bonds from atom data.\n");
          return MOLFILE_ERROR;
        }

        // Read each bond in the line
        for ( curr = bond_records; (next = strchr(curr, ' ')) != NULL; 
              curr = next+1 ) {
          *next = '\0';
          tmp_to = hash_lookup(&hasharray[mol], curr);
          if (tmp_to == HASH_FAIL) {
            fprintf(stderr, "mdfplugin) Could not find atom in hash table.\n");
            return MOLFILE_ERROR;
          }
          else if (tmp_to > atom) {
            // Only count bonds to atoms greater than the current one, since
            // each bond is listed twice
            *fromptr = atom;
            *toptr = tmp_to;
            fromptr++;
            toptr++;
          }
        }

        atom++;
      }

      fgets(line, LINESIZE, mdf->file);
      if ( ferror(mdf->file) || feof(mdf->file) ) {
        fprintf(stderr, "mdfplugin) File error while reading bonds.\n");
        return MOLFILE_ERROR;
      }
    }

    mol++;
  }

  for (mol = 0; mol < mdf->nmols; mol++) {
    hash_destroy(&hasharray[mol]);
  }
  delete [] hasharray;
  delete [] atomnames;

  *nbonds = bond_count;
  *from_data = mdf->from;
  *to_data = mdf->to;
  *bondorderptr = NULL; // not implemented yet
  *bondtype = NULL;
  *nbondtypes = 0;
  *bondtypename = NULL;

  return MOLFILE_SUCCESS;
}

// Free the memory used by the mdf structure
static void close_mdf_read(void *v) {
  mdfdata *mdf = (mdfdata *)v;
  if (mdf) {
    if (mdf->file) fclose(mdf->file);
    if (mdf->from) delete [] mdf->from;
    if (mdf->to)   delete [] mdf->to;
    delete mdf;
  }
}

// Plugin Initialization
static molfile_plugin_t plugin;

VMDPLUGIN_EXTERN int VMDPLUGIN_init(void) { 
  memset(&plugin, 0, sizeof(molfile_plugin_t));
  plugin.abiversion = vmdplugin_ABIVERSION;
  plugin.type = MOLFILE_PLUGIN_TYPE;
  plugin.name = "mdf";
  plugin.prettyname = "InsightII MDF";
  plugin.author = "Eamon Caddigan";
  plugin.majorv = 0;
  plugin.minorv = 4;
  plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
  plugin.filename_extension = "mdf";
  plugin.open_file_read = open_mdf_read;
  plugin.read_structure = read_mdf_structure;
  plugin.read_bonds = read_mdf_bonds;
  plugin.close_file_read = close_mdf_read;
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_EXTERN int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_EXTERN int VMDPLUGIN_fini(void) { return VMDPLUGIN_SUCCESS; }

