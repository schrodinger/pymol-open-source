/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_vaspxdatcarplugin
#define STATIC_PLUGIN 1

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: vaspxdatcarplugin.c,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.9 $       $Date: 2009/06/22 21:37:42 $
 *
 ***************************************************************************/

/*
 *  VASP plugins for VMD
 *  Sung Sakong, Dept. of Phys., Univsity Duisburg-Essen
 *  
 *  VASP manual   
 *  http:

 * 
 *  LINUX
 *  gcc -O2 -Wall -I. -I$VMDBASEDIR/plugins/include -c vaspxdatcarplugin.c
 *  ld -shared -o vaspxdatcarplugin.so vaspxdatcarplugin.o
 *
 *  MACOSX
 *  c++ -O2 -Wall -I. -I$VMDBASEDIR/plugins/include -c vaspxdatcarplugin.c
 *  c++ -bundle -o vaspxdatcarplugin.so vaspxdatcarplugin.o
 *
 *  Install
 *  copy vaspxdatcarplugin.so $VMDBASEDIR/plugins/$ARCH/molfile
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "molfile_plugin.h"
#include "vaspplugin.h"
#include "periodic_table.h"


static void *open_vaspxdatcar_read(const char *filename, const char *filetype, int *natoms)
{
  FILE *poscar;
  vasp_plugindata_t *data;
  char lineptr[LINESIZE], poscarfile[1000];
  float lc;
  int i;
  
  /* Verify that input is OK */
  if (!filename || !natoms) return NULL;

  /* Start with undefined value; set it after successful read */
  *natoms = MOLFILE_NUMATOMS_UNKNOWN;

  /* Use POSCAR or CONTCAR file for cell data and number of atoms */
  if (strstr(filename, "XDATCAR") == NULL) {
    fprintf(stderr, "\n\nVASP XDATCAR read) ERROR: file name '%s' does not contain 'XDATCAR'.\n", filename);
    return NULL;
  }
  strcpy(poscarfile, filename);
  strcpy(strstr(poscarfile, "XDATCAR"), "POSCAR");
  poscar = fopen(poscarfile, "r");
  if (!poscar) {
    strcpy(poscarfile, filename);
    strcpy(strstr(poscarfile, "XDATCAR"), "CONTCAR");
    poscar = fopen(poscarfile, "r");
    if (!poscar) {
      fprintf(stderr, "\n\nVASP XDATCAR read) ERROR: corresponding POSCAR or CONTCAR file not found.\n");
      return NULL;
    }
  }
  fprintf(stderr, "\n\nVASP XDATCAR read) determining lattice vectors and number of atoms from file '%s'.\n", poscarfile);

  data = vasp_plugindata_malloc();
  if (!data) return NULL;

  /* VASP4 is assumed in default */
  data->version = 4;
  data->file = fopen(filename, "rb");
  if (!data->file) {
    vasp_plugindata_free(data);
    return NULL;
  }

  data->filename = strdup(filename);

  fgets(lineptr, LINESIZE, poscar);
  data->titleline = strdup(lineptr);

  fgets(lineptr, LINESIZE, poscar);
  lc = atof(strtok(lineptr, " "));
  
  for (i = 0; i < 3; ++i) {
    float x, y, z;
    fgets(lineptr, LINESIZE, poscar);
    if (3 != sscanf(lineptr, "%f %f %f", &x, &y, &z)) {
      vasp_plugindata_free(data);
      fprintf(stderr, "\n\nVASP XDATCAR read) ERROR: POSCAR file '%s' does not have lattice vectors.\n", poscarfile);
      return NULL;
    }
    data->cell[i][0] = x*lc;
    data->cell[i][1] = y*lc;
    data->cell[i][2] = z*lc;
  }
  vasp_buildrotmat(data);

  /* Read the number of atoms per atom type */
  data->numatoms = 0;
  fgets(lineptr, LINESIZE, poscar);
  for (i = 0; i < MAXATOMTYPES; ++i) {
    char const *tmplineptr = strdup(lineptr);
    char const *token = (i == 0 ? strtok(lineptr, " ") : strtok(NULL, " "));
    int const n = (token ? atoi(token) : -1);

    /* if fails to read number of atoms, then assume VASP5 */
    if (i == 0 && n <= 0) {
      data->version = 5;
      data->titleline =  strdup(tmplineptr);
      fgets(lineptr, LINESIZE, poscar);
      break;
    }else if (n <= 0) break;

    data->eachatom[i] = n;
    data->numatoms += n;
  }

  if (data->version == 5) {
    data->numatoms = 0;
    for (i = 0; i < MAXATOMTYPES; ++i) {
      char const *token = (i == 0 ? strtok(lineptr, " ") : strtok(NULL, " "));
      int const n = (token ? atoi(token) : -1);
      
      if (n <= 0) break;
      
      data->eachatom[i] = n;
      data->numatoms += n;
    }
  }

  fclose(poscar);

  if (data->numatoms == 0) {
    vasp_plugindata_free(data);
    fprintf(stderr, "\n\nVASP XDATCAR read) ERROR: POSCAR file '%s' does not have the list of atom numbers.\n", poscarfile);
    return NULL;
  }

  *natoms = data->numatoms;

  return data;
}


static int read_vaspxdatcar_structure(void *mydata, int *optflags, molfile_atom_t *atoms)
{
  vasp_plugindata_t *data = (vasp_plugindata_t *)mydata;
  FILE *potcar;
  char lineptr[LINESIZE], potcarfile[1000];
  int atomcount, i;

  /* Verify that input is OK */
  if (!data || !optflags || !atoms) return MOLFILE_ERROR;

  *optflags = MOLFILE_MASS; /* we set atom mass from the PTE. */
  *optflags |= MOLFILE_ATOMICNUMBER | MOLFILE_RADIUS; 

  /* Read POTCAR file to determine atom types.
   * Each atom type section in POTCAR starts with a line
   * that contains the name of the element (H, He, C etc.).
   * Otherwise, try the title line instead.
   */
  strcpy(potcarfile, data->filename);
  strcpy(strstr(potcarfile, "XDATCAR"), "POTCAR");
  potcar = fopen(potcarfile, "r");
  if (potcar) fprintf(stderr, "\n\nVASP XDATCAR read) using file '%s' for determining atom types.\n", potcarfile);

  for (atomcount = i = 0; atomcount < data->numatoms; ++i) {
    int j;
    int idx;
    char const *label;
    float mass, radius;

    if (potcar) {
       /* Obtain atom types from POTCAR file */
       char atomtype[5] = "X";
       if (fgets(lineptr, LINESIZE, potcar)) sscanf(lineptr, "%*s %4s", atomtype);
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
    /* This can generally not happen; if it happens, there's a bug in the for-loop(i) above */
    fprintf(stderr, "\n\nVASP XDATCAR read) ERROR: problem occurred when setting the atom types.\n");
    return MOLFILE_ERROR;
  }

  /* Ignore header until X,Y,Z-coordinates */
  for (i = 0; i < 4; ++i) fgets(lineptr, LINESIZE, data->file);

 /* Determine VASP4 and VASP5 */
  if (tolower(lineptr[0]) == 'd'){
    data->version = 5;
    fgets(lineptr, LINESIZE, data->file);
  }
  else {
      data->version = 4;
      for (i = 0; i < 2; ++i)  fgets(lineptr, LINESIZE, data->file);
    }

  /* Check whether all coordinates are present in the file */
  for (i = 0; i < data->numatoms && fgets(lineptr, LINESIZE, data->file); ++i) {
    float coord;
    if (3 != sscanf(lineptr, "%f %f %f", &coord, &coord, &coord)) break;
  }
  if (i != data->numatoms) {
    fprintf(stderr, "\n\nVASP XDATCAR read) ERROR: file '%s' does not contain all coordinates of the atoms.\n", data->filename);
    return MOLFILE_ERROR;
  }

  /* Set file pointer to the line of the atoms' coordinates */
  rewind(data->file);
  for (i = 0; i < 10 - data->version; ++i) fgets(lineptr, LINESIZE, data->file);

  return MOLFILE_SUCCESS;
}


static int read_vaspxdatcar_timestep(void *mydata, int natoms, molfile_timestep_t *ts)
{
  vasp_plugindata_t *data = (vasp_plugindata_t *)mydata;
  char lineptr[LINESIZE];
  int i;

  /* only save coords if we're given a timestep pointer, */
  /* otherwise assume that VMD wants us to skip past it. */
  if (!data || !ts) return MOLFILE_EOF;

  for (i = 0; i < data->numatoms && fgets(lineptr, LINESIZE, data->file); ++i) {
    float x, y, z, rotx, roty, rotz;
    if (3 != sscanf(lineptr, "%f %f %f", &x, &y, &z)) break;

    rotx = x*data->cell[0][0] + y*data->cell[1][0] + z*data->cell[2][0];
    roty = x*data->cell[0][1] + y*data->cell[1][1] + z*data->cell[2][1];
    rotz = x*data->cell[0][2] + y*data->cell[1][2] + z*data->cell[2][2];

    ts->coords[3*i  ] = data->rotmat[0][0]*rotx + data->rotmat[0][1]*roty + data->rotmat[0][2]*rotz;
    ts->coords[3*i+1] = data->rotmat[1][0]*rotx + data->rotmat[1][1]*roty + data->rotmat[1][2]*rotz;
    ts->coords[3*i+2] = data->rotmat[2][0]*rotx + data->rotmat[2][1]*roty + data->rotmat[2][2]*rotz;
  }
  if (i != data->numatoms) return MOLFILE_EOF;

  /* Skip the empty line after coordinates */
  fgets(lineptr, LINESIZE, data->file);

  vasp_timestep_unitcell(ts, data);

  return MOLFILE_SUCCESS;
}


static void close_vaspxdatcar_read(void *mydata)
{
  vasp_plugindata_t *data = (vasp_plugindata_t *)mydata;
  vasp_plugindata_free(data);
}


/* registration stuff */
static molfile_plugin_t vaspxdatcarplugin = {
  vmdplugin_ABIVERSION,
  MOLFILE_PLUGIN_TYPE,         /* type */
  "XDATCAR",                   /* name */
  "VASP_XDATCAR",              /* pretty name */
  "Sung Sakong",               /* author */
  0,                           /* major version */
  6,                           /* minor version */
  VMDPLUGIN_THREADUNSAFE,      /* is reentrant */

  "XDATCAR",                   /* filename_extension */
  open_vaspxdatcar_read,       /* open_file_read */
  read_vaspxdatcar_structure,  /* read_structure */
  0,                           /* read_bonds */
  read_vaspxdatcar_timestep,   /* read_next_timestep */
  close_vaspxdatcar_read,      /* close_file_read */
  0,                           /* open_file_write */
  0,                           /* write_structure */
  0,                           /* write_timestep */
  0,                           /* close_file_write */
  0,                           /* read_volumetric_metadata */
  0,                           /* read_volumetric_data */
  0,                         
  0,                         
  0                          /* read_rawgraphics */
                             /* read_molecule_metadata */
                             /* write_bonds */
};

int VMDPLUGIN_init() {
  return VMDPLUGIN_SUCCESS;
}

int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)&vaspxdatcarplugin);
  return VMDPLUGIN_SUCCESS;
}

int VMDPLUGIN_fini() {
  return VMDPLUGIN_SUCCESS;
}
