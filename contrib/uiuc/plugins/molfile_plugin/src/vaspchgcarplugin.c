/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_vaspchgcarplugin
#define STATIC_PLUGIN 1

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: vaspchgcarplugin.c,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.10 $       $Date: 2009/02/24 23:33:21 $
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
 *  g++ -O2 -Wall -I. -I$VMDBASEDIR/plugins/include -c vaspchgcarplugin.c
 *  ld -shared -o vaspchgcarplugin.so vaspchgcarplugin.o
 *
 *  MACOSX
 *  c++ -O2 -Wall -I. -I$VMDBASEDIR/plugins/include -c vaspchgcarplugin.c
 *  c++ -bundle -o vaspchgcarplugin.so vaspchgcarplugin.o
 *
 *  Install
 *  copy vaspchgcarplugin.so $VMDBASEDIR/plugins/$ARCH/molfile
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "molfile_plugin.h"
#include "vaspplugin.h"


static void *open_vaspchgcar_read(const char *filename, const char *filetype, int *natoms)
{
  vasp_plugindata_t *data;
  char lineptr[LINESIZE];
  float lc;
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

  /* Read system title */
  fgets(lineptr, LINESIZE, data->file);
  data->titleline = strdup(lineptr);

  /* Read lattice constant */
  fgets(lineptr, LINESIZE, data->file);
  lc = atof(strtok(lineptr, " "));

  /* Read unit cell lattice vectors and multiply by lattice constant */
  for(i = 0; i < 3; ++i) {
    float x, y, z;
    fgets(lineptr, LINESIZE, data->file);
    sscanf(lineptr, "%f %f %f", &x, &y, &z);
    data->cell[i][0] = x*lc;
    data->cell[i][1] = y*lc;
    data->cell[i][2] = z*lc;
  }

  /* Build rotation matrix */
  vasp_buildrotmat(data);

  /* Count number of atoms */
  fgets(lineptr, LINESIZE, data->file);
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
    fprintf(stderr, "\n\nVASP CHGCAR read) ERROR: file '%s' does not contain list of atom numbers.\n", filename);
    return NULL;
  }

  /* Skip lines up to the grid numbers */
  for (i = 0; i < data->numatoms + 2; ++i) fgets(lineptr, LINESIZE, data->file);

  *natoms = data->numatoms;

  return data;
}


static int read_vaspchgcar_metadata(void *mydata, int *nvolsets, molfile_volumetric_t **metadata)
{
  vasp_plugindata_t *data = (vasp_plugindata_t *)mydata;
  char lineptr[LINESIZE];
  int gridx, gridy, gridz, i;
  char const spintext[4][20] = { "spin up+down", "spin up-down", "spin up", "spin down" };

  /* Verify that input is OK */
  if (!data || !nvolsets || !metadata) return MOLFILE_ERROR;

  /* Read the grid size */
  fgets(lineptr, LINESIZE, data->file);
  if (3 != sscanf(lineptr, "%d %d %d", &gridx, &gridy, &gridz)) {
     fprintf(stderr, "\n\nVASP CHGCAR read) ERROR: file '%s' does not contain grid dimensions.\n", data->filename);
     return MOLFILE_ERROR;
  }

  fprintf(stderr, "\n\nVASP CHGCAR read) found grid data block...\n");

  /* Initialize the volume set list with 4 entries:
   * spin up+down : always present
   * spin up-down / spin up /spin down : only there for spin-polarized calculations
   *                (the latter remain empty for non-spin-polarized calculations)
   */
  data->nvolsets = 4;
  data->vol = (molfile_volumetric_t *)malloc(data->nvolsets * sizeof(molfile_volumetric_t));
  if (!data->vol) {
     fprintf(stderr, "\n\nVASP CHGCAR read) ERROR: Cannot allocate space for volume data.\n");
     return MOLFILE_ERROR;
  }

  for (i = 0; i < data->nvolsets; ++i) {
    molfile_volumetric_t *const set = &(data->vol[i]); /* get a handle to the current volume set meta data */
    int k;

    set->has_color = 0;

    /* put volume data name */
    sprintf(set->dataname, "Charge density (%s)", spintext[i]);

    set->origin[0] = set->origin[1] = set->origin[2] = 0;
    set->xsize = gridx + 1;
    set->ysize = gridy + 1;
    set->zsize = gridz + 1;

    /* Rotate unit cell vectors */
    for (k = 0; k < 3; ++k) {
      set->xaxis[k] = data->rotmat[k][0] * data->cell[0][0]
		+ data->rotmat[k][1] * data->cell[0][1]
		+ data->rotmat[k][2] * data->cell[0][2];
      
      set->yaxis[k] = data->rotmat[k][0] * data->cell[1][0] 
		+ data->rotmat[k][1] * data->cell[1][1]
		+ data->rotmat[k][2] * data->cell[1][2];
      
      set->zaxis[k] = data->rotmat[k][0] * data->cell[2][0] 
		+ data->rotmat[k][1] * data->cell[2][1]
		+ data->rotmat[k][2] * data->cell[2][2];
    }
  }

  *nvolsets = data->nvolsets;
  *metadata = data->vol;  

  return MOLFILE_SUCCESS;
}


static int read_vaspchgcar_data(void *mydata, int set, float *datablock, float *colorblock)
{
  vasp_plugindata_t *data = (vasp_plugindata_t *)mydata;
  char lineptr[LINESIZE];
  int chargedensity, error, iset, n;
  float volume;

  /* Verify that input is OK */
  if (!data || !datablock) return MOLFILE_ERROR;
  if (set >= data->nvolsets) return MOLFILE_ERROR;

  if (strstr(data->filename, "LOCPOT") == NULL && strstr(data->filename, "ELFCAR") == NULL) {
    chargedensity = 1;
    fprintf(stderr, "\nVASP CHGCAR read) Charge density is assumed. Each value will be divided by unit cell volume.\n");
  } else {
    if (set == 1) {
      fprintf(stderr, "\n\nVASP CHGCAR read) ERROR: ELF or local potential do not include spin difference information.\n");
      return MOLFILE_ERROR;
    }
    chargedensity = 0;
    fprintf(stderr, "\nVASP CHGCAR read) ELF or local potential is assumed.\n");
  }

  volume = fabs(
            data->cell[0][0]*(data->cell[1][1]*data->cell[2][2]-data->cell[1][2]*data->cell[2][1])
	  + data->cell[0][1]*(data->cell[1][2]*data->cell[2][0]-data->cell[1][0]*data->cell[2][2])
	  + data->cell[0][2]*(data->cell[1][0]*data->cell[2][1]-data->cell[2][0]*data->cell[1][1])
               );

  /* Set file pointer to beginning of file and then skip header up to density data */
  rewind(data->file);
  for (n = 0; n < data->numatoms + 9; ++n) fgets(lineptr, LINESIZE, data->file);

  for(error = iset = 0; iset <= set && iset < 2 && !error; ++iset) {
    char const *dataname = data->vol[iset].dataname;
    int const xsize = data->vol[iset].xsize; 
    int const ysize = data->vol[iset].ysize;
    int const zsize = data->vol[iset].zsize;
    int const numberOfDatapoints = (xsize - 1) * (ysize - 1) * (zsize - 1);
    int ix, iy, iz;

    fprintf(stderr, "\nVASP CHGCAR read) Patience! Reading volume set %d (%d points): %s\n", iset + 1, numberOfDatapoints, dataname);

    for (n = iz = 0; iz < zsize; ++iz) {
      for (iy = 0; iy < ysize; ++iy) {
        for (ix = 0; ix < xsize; ++ix, ++n) {
          float value;
	  if (ix == xsize - 1) value = datablock[n - ix];
	  else if (iy == ysize - 1) value = datablock[n - iy*xsize];
	  else if (iz == zsize - 1) value = datablock[n - iz*ysize*xsize];
	  else {
	    if (1 != fscanf(data->file, "%f", &value)) return MOLFILE_ERROR;
	    if (chargedensity) value /= volume;
	  }

	  /* for set == 2: spin-up   = 0.5 * set0 + 0.5 * set1
	   * for set == 3: spin-down = 0.5 * set0 - 0.5 * set1 */
	  if (iset == set) datablock[n] = value;
	  else if (set >= 2 && iset == 0) datablock[n] = 0.5 * value;
	  else if (set == 2 && iset == 1) datablock[n] += 0.5 * value;
	  else if (set == 3 && iset == 1) datablock[n] -= 0.5 * value;
        }
      }
    }

    /* Skip paw-augmentation part 
     * augmentation parts are found only in CHGCAR not in PARCHG
     * I have no good idea for that */
    for (iy = 0; iy < data->numatoms; ++iy) {
      int numaug;
      if (1 != fscanf(data->file, "%*s %*s %*d %d", &numaug)) error = 1;
      
      for (ix = 0; ix < numaug && !error; ++ix) {
        float val;
        if (1 != fscanf(data->file, "%f", &val)) error = 2;
      }
      fgets(lineptr, LINESIZE, data->file);
    }

    /* After the charge density data there is float numbers (number of atoms)
     * and a line with three grid integers, all of which should be skipped. */
    if(iset==0){
      for (iy = 0; iy < data->numatoms; ++iy) {
        float val;
        if (1 != fscanf(data->file, "%f", &val)) error = 2;
      }
      for (iy = 0; iy < 3; ++iy) {
	int ival;
	if (1 != fscanf(data->file, "%d", &ival)) error = 2;
      }
    }

    fprintf(stderr, "\nVASP CHGCAR read) %s finished.\n", dataname);
  }

  if (error) fprintf(stderr, "\nVASP CHGCAR read) PAW-augmentation part is incomplete, but it is ignored anyway.\n");

  return MOLFILE_SUCCESS;
}


static void close_vaspchgcar_read(void *mydata)
{
  vasp_plugindata_t *data = (vasp_plugindata_t *)mydata;
  vasp_plugindata_free(data);
}


static molfile_plugin_t vaspchgcarplugin = {
  vmdplugin_ABIVERSION,
  MOLFILE_PLUGIN_TYPE,       /* type */ 	
  "CHGCAR",                  /* name */
  "VASP_CHGCAR",             /* pretty name */
  "Sung Sakong",             /* author */
  0,                         /* major version */
  5,                         /* minor version */
  VMDPLUGIN_THREADUNSAFE,    /* is reentrant */

  "CHGCAR",                  /* filename_extension */
  open_vaspchgcar_read,      /* open_file_read */          
  0,                         /* read_structure */
  0,                         /* read_bonds */
  0,                         /* read_next_timestep */
  close_vaspchgcar_read,     /* close_file_read */
  0,                         /* open_file_write */
  0,                         /* write_structure */
  0,                         /* write_timestep */
  0,                         /* close_file_write */
  read_vaspchgcar_metadata,  /* read_volumetric_metadata */
  read_vaspchgcar_data,      /* read_volumetric_data */
  0,                         
  0,                         
  0                          /* read_rawgraphics */
                             /* read_molecule_metadata */
                             /* write_bonds */
};

int VMDPLUGIN_init(void) {
  return VMDPLUGIN_SUCCESS;
}

int VMDPLUGIN_fini(void) {
  return VMDPLUGIN_SUCCESS;
}

int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)(void *)&vaspchgcarplugin);
  return VMDPLUGIN_SUCCESS;
}
