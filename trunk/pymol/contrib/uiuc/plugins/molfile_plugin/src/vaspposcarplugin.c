/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_vaspposcarplugin
#define STATIC_PLUGIN 1

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: vaspposcarplugin.c,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.6 $       $Date: 2009/01/29 14:57:00 $
 *
 ***************************************************************************/

/*
 *  VASP plugins for VMD
 *  Sung Sakong, Dept. of Theo. Chem., Univsity of Ulm 
 *  
 *  VASP manual   
 *  http:

 * 
 *  LINUX
 *  g++ -O2 -Wall -I. -I$VMDBASEDIR/plugins/include -c vaspposcarplugin.c
 *  ld -shared -o vaspposcarplugin.so vaspposcarplugin.o
 *
 *  MACOSX
 *  c++ -O2 -Wall -I. -I$VMDBASEDIR/plugins/include -c vaspposcarplugin.c
 *  c++ -bundle -o vaspposcarplugin.so vaspposcarplugin.o
 *
 *  Install
 *  copy vaspxatcarplugin.so $VMDBASEDIR/plugins/$ARCH/molfile
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "molfile_plugin.h"
#include "vaspplugin.h"
#include "periodic_table.h"


static void *open_vaspposcar_read(const char *filename, const char *filetype, int *natoms)
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

  data->file = fopen(filename, "rb");
  if (!data->file) {
    vasp_plugindata_free(data);
    return NULL;
  }

  data->filename = strdup(filename);

  /* Read title line */
  fgets(lineptr, LINESIZE, data->file);
  data->titleline = strdup(lineptr);

  /* Ignore rest of header up to the line with atom numbers */
  for (i = 0; i < 5; ++i) fgets(lineptr, LINESIZE, data->file);

  /* Read the number of atoms per atom type */
  data->numatoms = 0;
  for (i = 0; i < MAXATOMTYPES; ++i) {
    char const *token = (i == 0 ? strtok(lineptr, " ") : strtok(NULL, " "));
    int const n = (token ? atoi(token) : -1);

    if (n <= 0) break;

    data->eachatom[i] = n;
    data->numatoms += n;
  }

  if (data->numatoms == 0) {
    vasp_plugindata_free(data);
    fprintf(stderr, "\n\nVASP POSCAR read) ERROR: file '%s' does not have list of atom numbers.\n", filename);
    return NULL;
  }

  *natoms = data->numatoms;
  rewind(data->file);

  return data;
}


static int read_vaspposcar_structure(void *mydata, int *optflags, molfile_atom_t *atoms)
{
  vasp_plugindata_t *data = (vasp_plugindata_t *)mydata;
  FILE *potcar = NULL;
  int atomcount, i;
  char lineptr[LINESIZE], potcarfile[1000], *cp;
 
  if (!data || !optflags || !atoms) return MOLFILE_ERROR;

  *optflags = MOLFILE_MASS; /* we set atom mass from the PTE. */
  *optflags |= MOLFILE_ATOMICNUMBER | MOLFILE_RADIUS; 

  /* This plugin can read both POSCAR and CONTCAR files */
  strcpy(potcarfile, data->filename);
  cp = strstr(potcarfile, "POSCAR");
  if (!cp) cp = strstr(potcarfile, "CONTCAR");

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
    fprintf(stderr, "\n\nVASP POSCAR read) ERROR: file '%s' doesn't seem to have list of atoms.\n", data->filename);
    return MOLFILE_ERROR;
  }

 /* Ignore header until X,Y,Z-coordinates */
 for (i = 0; i < 7; ++i) fgets(lineptr, LINESIZE, data->file);

 /* Ignore selective tag-line, starting with either 's' or 'S'. */
 if (tolower(lineptr[0]) == 's') fgets(lineptr, LINESIZE, data->file);

 /* Check whether all coordinates are present in the file */
 for (i = 0; i < data->numatoms; ++i) {
   float coord;
   fgets(lineptr, LINESIZE, data->file);
   if (3 != sscanf(lineptr, "%f %f %f", &coord, &coord, &coord)) {
     fprintf(stderr, "\n\nVASP POSCAR read) ERROR: structure is missing type or coordinate(s) in file '%s' for atom '%d'\n", data->filename, i+1);
     return MOLFILE_ERROR;
   }
 }

 rewind(data->file);

 return MOLFILE_SUCCESS;
}


static int read_vaspposcar_timestep(void *mydata, int natoms, molfile_timestep_t *ts)
{
  int i, direct;
  char lineptr[LINESIZE];
  float lc;
  
  vasp_plugindata_t *data = (vasp_plugindata_t *)mydata;

  /* Save coords only if we're given a timestep pointer,
   * otherwise assume that VMD wants us to skip past it.
   */
  if (!ts || !data) return MOLFILE_EOF;

  /* VMD keeps calling for a next timestep, until we reach End-Of-File here */
  if (fgets(lineptr, LINESIZE, data->file) == NULL) return MOLFILE_EOF;

  fgets(lineptr, LINESIZE, data->file);
  sscanf(lineptr, "%f", &lc);

  for (i = 0; i < 3; ++i) {
    float x, y, z;
    fgets(lineptr, LINESIZE, data->file);
    sscanf(lineptr, "%f %f %f", &x, &y, &z);
    data->cell[i][0] = x*lc;
    data->cell[i][1] = y*lc;
    data->cell[i][2] = z*lc;
  }
  vasp_buildrotmat(data);

  /* Skip numbers of atom types */
  for (i = 0; i < 2; ++i) fgets(lineptr, LINESIZE, data->file);

  /* Skip selective tag-line, starting with 's' or 'S'. */
  if (tolower(lineptr[0]) == 's') fgets(lineptr, LINESIZE, data->file);

  /* Detect direct coordinates tag, starting with 'd' or 'D'. */
  direct = (tolower(lineptr[0]) == 'd' ? 1 : 0);

  for (i = 0; i < data->numatoms; ++i) {
    float x, y, z, rotx, roty, rotz;
    fgets(lineptr, LINESIZE, data->file);
    if (3 != sscanf(lineptr, "%f %f %f", &x, &y, &z)) {
      fprintf(stderr, "VASP POSCAR read) missing type or coordinate(s) in file '%s' for atom '%d'\n", data->filename, i+1);
      return MOLFILE_EOF;
    }

    if (direct) {
       rotx = x*data->cell[0][0]+y*data->cell[1][0]+z*data->cell[2][0];
       roty = x*data->cell[0][1]+y*data->cell[1][1]+z*data->cell[2][1];
       rotz = x*data->cell[0][2]+y*data->cell[1][2]+z*data->cell[2][2];
    } else {
       rotx = x*lc;
       roty = y*lc;
       rotz = z*lc;
    }
    ts->coords[3*i  ] = data->rotmat[0][0]*rotx+data->rotmat[0][1]*roty+data->rotmat[0][2]*rotz;
    ts->coords[3*i+1] = data->rotmat[1][0]*rotx+data->rotmat[1][1]*roty+data->rotmat[1][2]*rotz;
    ts->coords[3*i+2] = data->rotmat[2][0]*rotx+data->rotmat[2][1]*roty+data->rotmat[2][2]*rotz;
  }

  vasp_timestep_unitcell(ts, data);

  /* POSCAR type files have only one single timestep.
   * Therefore proceed till end of file after reading the coordinates.
   */
  fseek(data->file, 0, SEEK_END);

  return MOLFILE_SUCCESS;
}


static void close_vaspposcar_read(void *mydata)
{
  vasp_plugindata_t *data = (vasp_plugindata_t *)mydata;

  vasp_plugindata_free(data);
}


static void *open_vaspposcar_write(const char *filename, const char *filetype, int natoms)
{
  vasp_plugindata_t *data;

  data = vasp_plugindata_malloc();
  if (!data) return NULL;

  data->file = fopen(filename, "w");
  if (!data->file) {
    vasp_plugindata_free(data);
    fprintf(stderr, "VASP POSCAR write) ERROR: Unable to open vaspposcar file '%s' for writing\n", filename);
    return NULL;
  }

  data->filename = strdup(filename);
  data->numatoms = natoms;

  return data;
}


static int write_vaspposcar_structure(void *mydata, int optflags, const molfile_atom_t *atoms)
{
  vasp_plugindata_t *data = (vasp_plugindata_t *)mydata;

  if (!data || !atoms) return MOLFILE_ERROR;

  data->atomlist = (molfile_atom_t *)malloc(data->numatoms*sizeof(molfile_atom_t));
  if (!data->atomlist) return MOLFILE_ERROR;

  memcpy(data->atomlist, atoms, data->numatoms*sizeof(molfile_atom_t));

  return MOLFILE_SUCCESS;
}


static int write_vaspposcar_timestep(void *mydata, const molfile_timestep_t *ts)
{
  vasp_plugindata_t *data = (vasp_plugindata_t *)mydata; 
  molfile_atom_t const *atom;
  int i, maxtype, eachatom[MAXATOMTYPES];
  float x1, x2, y2, x3, y3, z3;
  char tmptype[LINESIZE] = "";

  if (!data || !ts) {
    fprintf(stderr, "VASP POSCAR write) ERROR: Wrong input for writing POSCAR file\n");
    return MOLFILE_ERROR;
  }

  /* restore unit cell information in Cartesian */
  x1 = ts->A;
  x2 = ts->B*cos(ts->gamma*M_PI/180.0);
  y2 = ts->B*sin(ts->gamma*M_PI/180.0);
  x3 = ts->C*cos(ts->beta*M_PI/180.0);
  y3 = (ts->B*ts->C*cos(ts->alpha*M_PI/180.0)-x2*x3)/y2;
  z3 = sqrt(ts->C*ts->C - x3*x3 - y3*y3);

  maxtype = -1;
  atom = data->atomlist;
  for (i = 0; i < data->numatoms && maxtype < MAXATOMTYPES - 1; ++i, ++atom) {
    if (strcmp(tmptype, atom->type) != 0) {
      fprintf(data->file, "%-2s  ", atom->type);
      eachatom[++maxtype] = 0;
    }
    eachatom[maxtype]++;
    strncpy(tmptype, atom->type, sizeof(atom->type));
  }

  fprintf(data->file, "\n%20.12f\n", 1.0);  
  fprintf(data->file, "%20.12f  %20.12f  %20.12f\n", x1, 0.0, 0.0);
  fprintf(data->file, "%20.12f  %20.12f  %20.12f\n", x2, y2,  0.0);
  fprintf(data->file, "%20.12f  %20.12f  %20.12f\n", x3, y3,  z3);

  for (i = 0; i <= maxtype; ++i) fprintf(data->file, " %d ", eachatom[i]);

  fprintf(data->file, "\nCartesian\n");

  for (i = 0; i < data->numatoms; ++i) {
    float const *pos = ts->coords + 3*i;
    fprintf(data->file, "%20.12f %20.12f %20.12f \n", pos[0], pos[1], pos[2]);
  }

  return MOLFILE_SUCCESS;
}


static void close_vaspposcar_write(void *mydata)
{
  vasp_plugindata_t *data = (vasp_plugindata_t *)mydata;
  vasp_plugindata_free(data);
}


/* registration stuff */
static molfile_plugin_t vaspposcarplugin = {
  vmdplugin_ABIVERSION,
  MOLFILE_PLUGIN_TYPE,         /* type */
  "POSCAR",                    /* name */
  "VASP_POSCAR",               /* pretty name */
  "Sung Sakong",               /* author */
  0,                           /* major version */
  3,                           /* minor version */
  VMDPLUGIN_THREADUNSAFE,      /* is not reentrant */

  "POSCAR",                    /* filename_extension */
  open_vaspposcar_read,        /* open_file_read */
  read_vaspposcar_structure,   /* read_structure */
  0,                           /* read_bonds */
  read_vaspposcar_timestep,    /* read_next_timestep */
  close_vaspposcar_read,       /* close_file_read */
  open_vaspposcar_write,       /* open_file_write */
  write_vaspposcar_structure,  /* write_structure */
  write_vaspposcar_timestep,   /* write_timestep */
  close_vaspposcar_write,      /* close_file_write */
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
  (*cb)(v, (vmdplugin_t *)(void *)&vaspposcarplugin);
  return VMDPLUGIN_SUCCESS;
}

int VMDPLUGIN_fini() {
  return VMDPLUGIN_SUCCESS;
}
