/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_vaspxmlplugin
#define STATIC_PLUGIN 1

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: vaspxmlplugin.c,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.5 $       $Date: 2008/08/05 20:13:49 $
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
 *  g++ -O2 -Wall -I. -I$VMDBASEDIR/plugins/include -c vaspxmlplugin.c
 *  ld -shared -o vaspxmlplugin.so vaspxmlplugin.o
 *
 *  MACOSX
 *  c++ -O2 -Wall -I. -I$VMDBASEDIR/plugins/include -c vaspxmlplugin.c
 *  c++ -bundle -o vaspxmlplugin.so vaspxmlplugin.o
 *
 *  Install
 *  copy vaspxmlplugin.so $VMDBASEDIR/plugins/$ARCH/molfile
 */

 /* Be aware that the XML file has 2 extra frames in addition to the
  * ionic time steps. There is a 'initial position' frame, then the
  * frames of each time step, and in the end once again a 'final position'
  * frame. Hence, the xml file will load two frames more than the ones
  * in XDATCAR or OUTCAR; the first two and last two are thus identical.
  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "molfile_plugin.h"
#include "vaspplugin.h"
#include "periodic_table.h"


static void *open_vaspxml_read(const char *filename, const char *filetype, int *natoms)
{
  vasp_plugindata_t *data;
  char lineptr[LINESIZE];
  int cellcoords, finished;

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

  /* Scan xml file */
  data->numatoms = cellcoords = finished = 0;
  while (fgets(lineptr, LINESIZE, data->file) && !finished) {

    if (strstr(lineptr, "SYSTEM") != NULL && data->titleline == NULL) {
       /* Extract title line */
       char *begin = strstr(lineptr, ">") + 1;
       char *end = strstr(lineptr, "</i>");
       if (end) end = '\0';
       if (begin) data->titleline = strdup(begin);

    } else if (strstr(lineptr, "atominfo") != NULL && data->numatoms == 0) {
       /* Extract number of atoms */
       fgets(lineptr, LINESIZE, data->file);
       sscanf(lineptr, " <atoms> %d </atoms>", &data->numatoms);

    } else if (strstr(lineptr, "crystal") != NULL && cellcoords == 0) {
       /* Extract lattice vectors */
       int i;
       fgets(lineptr, LINESIZE, data->file);
       for (i = 0; i < 3 && fgets(lineptr, LINESIZE, data->file); ++i) cellcoords += sscanf(lineptr, " <v> %f %f %f </v>", &data->cell[i][0], &data->cell[i][1], &data->cell[i][2]);
    }

    finished = data->titleline != NULL && data->numatoms != 0 && cellcoords != 0;
  }

  if (data->numatoms <= 0) {
     vasp_plugindata_free(data);
     fprintf(stderr, "\n\nVASP xml read) ERROR: file '%s' does not contain the number of atoms.\n", filename);
     return NULL;
  }

  if (cellcoords != 9) {
     vasp_plugindata_free(data);
     fprintf(stderr, "\n\nVASP xml read) ERROR: file '%s' does not contain lattice vectors.\n", filename);
     return NULL;
  }

  vasp_buildrotmat(data);

  *natoms = data->numatoms;
  rewind(data->file);

  return data;
}


static int read_vaspxml_structure(void *mydata, int *optflags, molfile_atom_t *atoms)
{
  vasp_plugindata_t *data = (vasp_plugindata_t *)mydata;
  int atomcount, coordscount, finished;
  char lineptr[LINESIZE];

  /* Verify that input is OK */
  if (!data || !optflags || !atoms) return MOLFILE_ERROR;

  *optflags = MOLFILE_MASS; /* we set atom mass from the PTE. */
  *optflags |= MOLFILE_ATOMICNUMBER | MOLFILE_RADIUS; 

  /* Scan xml file */
  atomcount = coordscount = finished = 0;
  while (fgets(lineptr, LINESIZE, data->file) && !finished) {

    /* Extract atom types */
    if (strstr(lineptr, "atomtype") != NULL && atomcount == 0) {
      int i;
      fgets(lineptr, LINESIZE, data->file);
      for (i = 0; i < data->numatoms; ++i, ++atomcount) {
        molfile_atom_t *atom = &(atoms[i]);
        char atomtype[5];
        int idx;
        fgets(lineptr, LINESIZE, data->file);
	if (1 != sscanf(lineptr, " <rc><c> %4s </c>", atomtype)) break;

        idx = get_pte_idx(atomtype);

        /* Required settings */
        strncpy(atom->name, get_pte_label(idx), sizeof(atom->name));
        strncpy(atom->type, atom->name, sizeof(atom->type));
	atom->resname[0] = '\0';
        atom->resid = 1;
        atom->segid[0]  ='\0';
        atom->chain[0] = '\0';

        /* Optional flags (as defined in *optflags) */
        atom->mass = get_pte_mass(idx);
        atom->radius = get_pte_vdw_radius(idx);
        atom->atomicnumber = idx;
      }

    /* Verify presence of coordinates for all atoms */
    } else if (strstr(lineptr, "positions") != NULL && coordscount == 0) {
        int i;
        for (i = 0; i < data->numatoms && fgets(lineptr, LINESIZE, data->file); ++i) {
           float x, y, z;
	   if (3 != sscanf(lineptr, " <v> %f %f %f <\v>", &x, &y, &z)) break;
        }
	coordscount = 3 * i;
    }

    finished = atomcount != 0 && coordscount != 0;
  }

  if (atomcount != data->numatoms) {
    fprintf(stderr, "\n\nVASP xml read) ERROR: file '%s' does not have list of atom names.\n", data->filename);
    return MOLFILE_ERROR;
  }

  if (coordscount != 3 * data->numatoms) {
    fprintf(stderr, "\n\nVASP xml read)  file '%s' does not contain coordinates of all atoms.\n", data->filename);
    return MOLFILE_ERROR;
  }

  rewind(data->file);

  return MOLFILE_SUCCESS;
}


static int read_vaspxml_timestep(void *mydata, int natoms, molfile_timestep_t *ts)
{
  vasp_plugindata_t *data = (vasp_plugindata_t *)mydata;
  char lineptr[LINESIZE];

  /* only save coords if we're given a timestep pointer, */
  /* otherwise assume that VMD wants us to skip past it. */
  if (!data || !ts) return MOLFILE_EOF;

  /* Scan xml file */
  while (fgets(lineptr, LINESIZE, data->file)) {

    /* Extract coordinates of all atoms */
    if (strstr(lineptr, "positions") != NULL) {
      int i;
      for (i = 0; i < data->numatoms && fgets(lineptr, LINESIZE, data->file); ++i) {
        float x, y, z, rotx, roty, rotz;
	if (3 != sscanf(lineptr, " <v> %f %f %f </v>", &x, &y, &z)) return MOLFILE_EOF;

        rotx = x*data->cell[0][0] + y*data->cell[1][0] + z*data->cell[2][0];
        roty = x*data->cell[0][1] + y*data->cell[1][1] + z*data->cell[2][1];
        rotz = x*data->cell[0][2] + y*data->cell[1][2] + z*data->cell[2][2];

        ts->coords[3*i  ] = data->rotmat[0][0]*rotx + data->rotmat[0][1]*roty + data->rotmat[0][2]*rotz;
        ts->coords[3*i+1] = data->rotmat[1][0]*rotx + data->rotmat[1][1]*roty + data->rotmat[1][2]*rotz;
        ts->coords[3*i+2] = data->rotmat[2][0]*rotx + data->rotmat[2][1]*roty + data->rotmat[2][2]*rotz;
      }
      vasp_timestep_unitcell(ts, data);
      return MOLFILE_SUCCESS;
    }
  }

  return MOLFILE_EOF;
}


static void close_vaspxml_read(void *mydata)
{
  vasp_plugindata_t *data = (vasp_plugindata_t *)mydata;
  vasp_plugindata_free(data);
}


/* registration stuff */
static molfile_plugin_t vaspxmlplugin = {
  vmdplugin_ABIVERSION,
  MOLFILE_PLUGIN_TYPE,     /* type */
  "xml",                   /* name */
  "VASP_xml",              /* pretty name */
  "Sung Sakong",           /* author */
  0,                       /* major version */
  3,                       /* minor version */
  VMDPLUGIN_THREADSAFE,    /* is reentrant */

  "xml",                   /* filename_extension */
  open_vaspxml_read,       /* open_file_read */
  read_vaspxml_structure,  /* read_structure */
  0,                       /* read_bonds */
  read_vaspxml_timestep,   /* read_next_timestep */
  close_vaspxml_read,      /* close_file_read */
  0,                       /* open_file_write */
  0,                       /* write_structure */
  0,                       /* write_timestep */
  0,                       /* close_file_write */
  0,                       /* read_volumetric_metadata */
  0,                       /* read_volumetric_data */
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
  (*cb)(v, (vmdplugin_t *)(void *)&vaspxmlplugin);
  return VMDPLUGIN_SUCCESS;
}

int VMDPLUGIN_fini() {
  return VMDPLUGIN_SUCCESS;
}
