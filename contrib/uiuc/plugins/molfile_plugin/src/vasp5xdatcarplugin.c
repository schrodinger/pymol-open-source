/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_vasp5xdatcarplugin
#define STATIC_PLUGIN 1

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: vasp5xdatcarplugin.c,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.3 $       $Date: 2014/10/10 14:41:01 $
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
 *  g++ -O2 -Wall -fPIC -I. -I$VMDBASEDIR/plugins/include -c vasp5xdatcarplugin.c
 *  ld -shared -o vasp5xdatcarplugin.so vasp5xdatcarplugin.o
 *
 *  MACOSX
 *  c++ -O2 -Wall -I. -I$VMDBASEDIR/plugins/include -c vasp5xdatcarplugin.c
 *  c++ -bundle -o vasp5xdatcarplugin.so vasp5xdatcarplugin.o
 *
 *  Install
 *  copy vasp5xdatcarplugin.so $VMDBASEDIR/plugins/$ARCH/molfile
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "molfile_plugin.h"
#include "vaspplugin.h"
#include "periodic_table.h"


static void *open_vasp5xdatcar_read(const char *filename, const char *filetype, int *natoms)
{
  vasp_plugindata_t *data;
  char lineptr[LINESIZE];
  int i;

  /* Verify that input is OK */
  if (!filename || !natoms) return NULL;

  /* Start with undefined value; set it after successful read */
  *natoms = MOLFILE_NUMATOMS_UNKNOWN;

  data = vasp_plugindata_malloc();
  if (!data) return NULL;

  /* VASP5 is assumed in default */
  data->version = 5;
  data->file = fopen(filename, "rb");
  if (!data->file) {
    vasp_plugindata_free(data);
    return NULL;
  }

  data->filename = strdup(filename);

  /* Ignore rest of header up to the line with atom numbers */
  for (i = 0; i < 5; ++i) fgets(lineptr, LINESIZE, data->file);

  /* Read title line */
  fgets(lineptr, LINESIZE, data->file);
  data->titleline = strdup(lineptr);

  /* Read the number of atoms per atom type */
  data->numatoms = 0;
  fgets(lineptr, LINESIZE, data->file);
  for (i = 0; i < MAXATOMTYPES; ++i) {
    char const *token = (i == 0 ? strtok(lineptr, " ") : strtok(NULL, " "));
    int const n = (token ? atoi(token) : -1);
    
    if (n <= 0) break;
    
    data->eachatom[i] = n;
    data->numatoms += n;
  }


  if (data->numatoms == 0) {
    vasp_plugindata_free(data);
    fprintf(stderr, "\n\nVASP5 XDATCAR read) ERROR: file '%s' does not have list of atom numbers.\n", filename);
    return NULL;
  }

  *natoms = data->numatoms;
  rewind(data->file);

  return data;
}


static int read_vasp5xdatcar_structure(void *mydata, int *optflags, molfile_atom_t *atoms)
{
  vasp_plugindata_t *data = (vasp_plugindata_t *)mydata;
  FILE *potcar = NULL;
  int atomcount, i;
  char lineptr[LINESIZE], potcarfile[1000], *cp;
  float lc;
 
  if (!data || !optflags || !atoms) return MOLFILE_ERROR;

  *optflags = MOLFILE_MASS; /* we set atom mass from the PTE. */
  *optflags |= MOLFILE_ATOMICNUMBER | MOLFILE_RADIUS; 

  strcpy(potcarfile, data->filename);
  cp = strstr(potcarfile, "XDATCAR");

  if (cp) {
    strcpy(cp, "POTCAR");
    potcar = fopen(potcarfile, "r");
  }

  /* Read POTCAR file to determine atom types.
   * Each atom type section in POTCAR starts with a line
   * that contains the name of the element (H, He, C etc.).
   * Otherwise try the title line instead.
   */
  for (atomcount = i = 0; atomcount < data->numatoms; ++i) {
    int idx, j;
    char const *label;
    float mass, radius;

    if (potcar) {
       char atomtype[5] = "X";
       /* Obtain atom types from POTCAR file */
       if (fgets(lineptr, LINESIZE, potcar)) sscanf(lineptr, "%*s %4[^_. 0-9]", atomtype);
       idx = get_pte_idx(atomtype);
       /* Skip lines in potcar file until next element */
       while (fgets(lineptr, LINESIZE, potcar)) if (strstr(lineptr, "End of Dataset")) break;
    } else {
       /* Try to obtain atom types from title line */
       char const *token = (i == 0 ? strtok(data->titleline, " ") : strtok(NULL, " "));
       idx = get_pte_idx(token);
    }

    label = get_pte_label(idx);
    mass = get_pte_mass(idx);
    radius = get_pte_vdw_radius(idx);
    for (j = 0; j < data->eachatom[i]; ++j, ++atomcount) {
      molfile_atom_t *const atom = &(atoms[atomcount]);

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
    fprintf(stderr, "\n\nVASP5 XDATCAR read) ERROR: file '%s' doesn't seem to have list of atoms.\n", data->filename);
    return MOLFILE_ERROR;
  }

  for (i = 0; i < 2; ++i) fgets(lineptr, LINESIZE, data->file);
  sscanf(lineptr, "%f", &lc);
   fprintf(stderr, "%f\n", lc);

  for (i = 0; i < 3; ++i) {
    float x, y, z;
    fgets(lineptr, LINESIZE, data->file);
    sscanf(lineptr, "%f %f %f", &x, &y, &z);
    data->cell[i][0] = x*lc;
    data->cell[i][1] = y*lc;
    data->cell[i][2] = z*lc;
  }
  vasp_buildrotmat(data);

 /* Ignore header until X,Y,Z-coordinates */
 for (i = 0; i < 3; ++i) fgets(lineptr, LINESIZE, data->file);

 /* Check whether all coordinates are present in the file */
 for (i = 0; i < data->numatoms; ++i) {
   float coord;
   fgets(lineptr, LINESIZE, data->file);
   if (3 != sscanf(lineptr, "%f %f %f", &coord, &coord, &coord)) {
     fprintf(stderr, "\n\nVASP5 XDATCAR read) ERROR: structure is missing type or coordinate(s) in file '%s' for atom '%d'\n", data->filename, i+1);
     return MOLFILE_ERROR;
   }
 }

 rewind(data->file);

 /* Ignore header until X,Y,Z-coordinates */
 for (i = 0; i < 8; ++i) fgets(lineptr, LINESIZE, data->file);

 return MOLFILE_SUCCESS;
}


static int read_vasp5xdatcar_timestep(void *mydata, int natoms, molfile_timestep_t *ts)
{
  int i;
  char lineptr[LINESIZE];
  
  vasp_plugindata_t *data = (vasp_plugindata_t *)mydata;

  /* Save coords only if we're given a timestep pointer,
   * otherwise assume that VMD wants us to skip past it.
   */
  if (!ts || !data) return MOLFILE_EOF;

  for (i = 0; i < data->numatoms; ++i) {
    float x, y, z, rotx, roty, rotz;
    fgets(lineptr, LINESIZE, data->file);
    if (3 != sscanf(lineptr, "%f %f %f", &x, &y, &z)) {
      fprintf(stderr, "VASP5 XDATCAR read) missing type or coordinate(s) in file '%s' for atom '%d'\n", data->filename, i+1);
      return MOLFILE_EOF;
    }

    rotx = x*data->cell[0][0]+y*data->cell[1][0]+z*data->cell[2][0];
    roty = x*data->cell[0][1]+y*data->cell[1][1]+z*data->cell[2][1];
    rotz = x*data->cell[0][2]+y*data->cell[1][2]+z*data->cell[2][2];

    ts->coords[3*i  ] = data->rotmat[0][0]*rotx+data->rotmat[0][1]*roty+data->rotmat[0][2]*rotz;
    ts->coords[3*i+1] = data->rotmat[1][0]*rotx+data->rotmat[1][1]*roty+data->rotmat[1][2]*rotz;
    ts->coords[3*i+2] = data->rotmat[2][0]*rotx+data->rotmat[2][1]*roty+data->rotmat[2][2]*rotz;
  }

  vasp_timestep_unitcell(ts, data);

  /* VMD keeps calling for a next timestep, until we reach End-Of-File here */
  if (fgets(lineptr, LINESIZE, data->file) == NULL) return MOLFILE_EOF;

  return MOLFILE_SUCCESS;
}


static void close_vasp5xdatcar_read(void *mydata)
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
  plugin.name = "VASP_XDATCAR5";
  plugin.prettyname = "VASP_XDATCAR5";
  plugin.author = "Sung Sakong";
  plugin.majorv = 0;
  plugin.minorv = 7;
  plugin.is_reentrant = VMDPLUGIN_THREADUNSAFE;
  plugin.filename_extension = "VASP_XDATCAR5";
  plugin.open_file_read = open_vasp5xdatcar_read;
  plugin.read_structure = read_vasp5xdatcar_structure;
  plugin.read_next_timestep = read_vasp5xdatcar_timestep;
  plugin.close_file_read = close_vasp5xdatcar_read;
  return VMDPLUGIN_SUCCESS;
}

int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

int VMDPLUGIN_fini() {
  return VMDPLUGIN_SUCCESS;
}
