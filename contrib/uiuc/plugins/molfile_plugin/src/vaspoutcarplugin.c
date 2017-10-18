/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_vaspoutcarplugin
#define STATIC_PLUGIN 1

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: vaspoutcarplugin.c,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.9 $       $Date: 2014/10/10 14:41:01 $
 *
 ***************************************************************************/

/*
 *  VASP plugins for VMD
 *  Sung Sakong, Dept. of Phys., Univsity Duisburg-Essen
 *  
 *  VASP manual   
 *  http://cms.mpi.univie.ac.at/vasp/
 * 
 *  LINUX
 *  g++ -O2 -Wall -I. -I$VMDBASEDIR/plugins/include -c vaspoutcarplugin.c
 *  ld -shared -o vaspoutcarplugin.so vaspoutcarplugin.o
 *
 *  MACOSX
 *  c++ -O2 -Wall -I. -I$VMDBASEDIR/plugins/include -c vaspoutcarplugin.c
 *  c++ -bundle -o vaspoutcarplugin.so vaspoutcarplugin.o
 *
 *  Install
 *  copy vaspoutcarplugin.so $VMDBASEDIR/plugins/$ARCH/molfile
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "molfile_plugin.h"
#include "vaspplugin.h"
#include "periodic_table.h"


static void *open_vaspoutcar_read(const char *filename, const char *filetype, int *natoms)
{
  vasp_plugindata_t *data;
  char lineptr[LINESIZE];

  /* Verify that input is OK */
  if (!filename || !natoms) return NULL;

  /* Start with undefined value; set it after successful read */
  *natoms = MOLFILE_NUMATOMS_UNKNOWN;

  data = vasp_plugindata_malloc();
  if (!data) return NULL;

  data->file = fopen(filename, "rb");
  if (!data->file) {
    vasp_plugindata_free(data);
    return NULL;
  }
  
  data->filename = strdup(filename);

  /* Catch total number of atoms */
  data->numatoms = 0;
  while (fgets(lineptr, LINESIZE, data->file) && data->numatoms == 0) {
   if (strstr(lineptr, "NIONS =") != NULL) {
      sscanf(lineptr, " %*[ a-zA-Z] = %*d %*[ a-zA-Z] = %d", &data->numatoms);
      break;
    }
  }

  if (data->numatoms <= 0) {
    vasp_plugindata_free(data);
    fprintf(stderr, "\n\nVASP OUTCAR read) ERROR: file '%s' does not contain the number of atoms.\n", filename);
    return NULL;
  }

  *natoms = data->numatoms;

  /* Catch the lattice vectors */
  while (fgets(lineptr, LINESIZE, data->file)) {
     if (strstr(lineptr, "direct lattice vectors") != NULL) {
       int i;
       for (i = 0; i < 3; ++i) {
	 fgets(lineptr, LINESIZE, data->file);
	 if (3 != sscanf(lineptr, "%f %f %f", &data->cell[i][0], &data->cell[i][1], &data->cell[i][2])) {
           vasp_plugindata_free(data);
           fprintf(stderr, "\n\nVASP OUTCAR read) ERROR: file '%s' does not contain lattice vectors.\n", filename);
           return NULL;
	 }
       }
       break;
     }
  }
  vasp_buildrotmat(data);

  rewind(data->file);

  return data;
}


static int read_vaspoutcar_structure(void *mydata, int *optflags, molfile_atom_t *atoms)
{
  vasp_plugindata_t *data = (vasp_plugindata_t *)mydata;
  FILE *potcar;
  char lineptr[LINESIZE], potcarfile[1000];
  float atommass[MAXATOMTYPES];
  int atomcount, typecount, i;

  if (!data || !optflags || !atoms) return MOLFILE_ERROR;

  *optflags = MOLFILE_MASS; /* we set atom mass from the PTE. */
  *optflags |= MOLFILE_ATOMICNUMBER | MOLFILE_RADIUS;

  typecount = 0;
  atomcount = 0;
  while (fgets(lineptr, LINESIZE, data->file) && atomcount < data->numatoms) {
    if (strstr(lineptr, "POMASS") != NULL) sscanf(lineptr, " POMASS = %f;", &atommass[typecount++]);

    if (strstr(lineptr, "ions per type =") != NULL) {
      char const *token = strtok(lineptr, "=");
      for (i = 0; i < typecount; ++i) {
        token = strtok(NULL, " ");
        atomcount += data->eachatom[i] = atoi(token);
      }
    }
  }

  if (atomcount != data->numatoms) {
    fprintf(stderr, "\n\nVASP OUTCAR read) ERROR: file '%s' does not have number of each atom.\n", data->filename);
    return MOLFILE_ERROR;
  }

  /* Read POTCAR file to determine atom types.
   * Each atom type section in POTCAR starts with a line
   * that contains the name of the element (H, He, C etc.).
   * Otherwise try to find similar mass in periodic table.
   */
  if (strstr(data->filename, "OUTCAR")) {
    strcpy(potcarfile, data->filename);
    strcpy(strstr(potcarfile, "OUTCAR"), "POTCAR");
    potcar = fopen(potcarfile, "r");
  } else {
    potcar = NULL;
  }

  for (atomcount = i = 0; atomcount < data->numatoms; ++i) {
    char const *label;
    float mass, radius;
    int k, idx = 0;

    if (potcar) {
      /* Obtain atom types from POTCAR file */
      char atomtype[5] = "X";
      if (fgets(lineptr, LINESIZE, potcar)) sscanf(lineptr, "%*s %4[^_. 0-9]", atomtype);
      idx = get_pte_idx(atomtype);

      /* Skip lines in potcar file until next element */
      while (fgets(lineptr, LINESIZE, potcar)) if (strstr(lineptr, "End of Dataset")) break;
    }
    
    if (idx == 0) {
        /* Try to find atom type by browsing through periodic table */
        idx = sizeof(pte_mass)/sizeof(pte_mass[0]);
	do idx--;
	while (idx > 0 && fabs(get_pte_mass(idx) - atommass[i]) > 0.01);
    }

    label = get_pte_label(idx);
    mass = (idx ? get_pte_mass(idx) : atommass[i]);
    radius = get_pte_vdw_radius(idx);
    for (k = 0; k < data->eachatom[i]; ++k, ++atomcount) {
      molfile_atom_t *atom = &(atoms[atomcount]);

      /* Required settings */
      strncpy(atom->name, label, sizeof(atom->name));
      strncpy(atom->type, atom->name, sizeof(atom->type));
      atom->resname[0] = '\0';
      atom->resid = 1;
      atom->segid[0]='\0';
      atom->chain[0]='\0';

      /* Optional flags (as defined in *optflags) */
      atom->mass = mass;
      atom->radius = radius;
      atom->atomicnumber = idx;
    }
  }
  if (potcar) fclose(potcar);

  if (atomcount != data->numatoms) {
    fprintf(stderr, "\n\nVASP OUTCAR read) ERROR: file '%s' does contain list of atom names.\n", data->filename);
    return MOLFILE_ERROR;
  }

  atomcount = 0;
  while (fgets(lineptr, LINESIZE, data->file) && atomcount == 0) {
    if (strstr(lineptr, "position of ions in cartesian coordinates") != NULL) {
      for (i = 0; i < data->numatoms; ++i, ++atomcount) {
	float coord;
        fgets(lineptr, LINESIZE, data->file);
        if (3 != sscanf(lineptr, "%f %f %f", &coord, &coord, &coord)) {
	  fprintf(stderr, "\n\nVASP OUTCAR read) missing type or coordinate(s) in file '%s' for atom '%d'\n", data->filename, i+1);
	  return MOLFILE_ERROR;
	}
      }
    }
  }

  if (atomcount != data->numatoms) {
    fprintf(stderr, "\n\nVASP OUTCAR read) ERROR: file '%s' does contain list of atom names.\n", data->filename);
    return MOLFILE_ERROR;
  }
 
  rewind(data->file);

  return MOLFILE_SUCCESS;
}


static int read_vaspoutcar_timestep(void *mydata, int natoms, molfile_timestep_t *ts)
{
  vasp_plugindata_t *data = (vasp_plugindata_t *)mydata;
  char lineptr[LINESIZE];
  int atomcount;

  /* Save coords only if we're given a timestep pointer,
   * otherwise assume that VMD wants us to skip past it.
   */
  if (!data || !ts) return MOLFILE_EOF;

  atomcount = 0;
  while (fgets(lineptr, LINESIZE, data->file) && atomcount == 0) {
    if (strstr(lineptr, "TOTAL-FORCE") != NULL) {
      int i;
      fgets(lineptr, LINESIZE, data->file);
      for (i = 0; i < data->numatoms; ++i, ++atomcount) {
	float x, y, z;
	fgets(lineptr, LINESIZE, data->file);
	if (3 != sscanf(lineptr, "%f %f %f", &x, &y, &z)) return MOLFILE_EOF;
      
	ts->coords[3*i  ] = data->rotmat[0][0]*x+data->rotmat[0][1]*y+data->rotmat[0][2]*z;
	ts->coords[3*i+1] = data->rotmat[1][0]*x+data->rotmat[1][1]*y+data->rotmat[1][2]*z;
	ts->coords[3*i+2] = data->rotmat[2][0]*x+data->rotmat[2][1]*y+data->rotmat[2][2]*z;
      }
    }
  }
  if (atomcount != data->numatoms) return MOLFILE_EOF;

  vasp_timestep_unitcell(ts, data);

  return MOLFILE_SUCCESS;
}


static void close_vaspoutcar_read(void *mydata)
{
  vasp_plugindata_t *data = (vasp_plugindata_t *)mydata;
  vasp_plugindata_free(data);
}


/* registration stuff */
static molfile_plugin_t plugin;

int VMDPLUGIN_init() {
  memset(&plugin, 0, sizeof(molfile_plugin_t));
  plugin.abiversion = vmdplugin_ABIVERSION;
  plugin.type = MOLFILE_PLUGIN_TYPE;
  plugin.name = "OUTCAR";
  plugin.prettyname = "VASP_OUTCAR";
  plugin.author = "Sung Sakong";
  plugin.majorv = 0;
  plugin.minorv = 7;
  plugin.is_reentrant = VMDPLUGIN_THREADUNSAFE;
  plugin.filename_extension = "OUTCAR";
  plugin.open_file_read = open_vaspoutcar_read;
  plugin.read_structure = read_vaspoutcar_structure;
  plugin.read_next_timestep = read_vaspoutcar_timestep;
  plugin.close_file_read = close_vaspoutcar_read;
  return VMDPLUGIN_SUCCESS;
}

int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

int VMDPLUGIN_fini() {
  return VMDPLUGIN_SUCCESS;
}
