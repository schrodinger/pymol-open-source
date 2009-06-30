/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: vaspplugin.h,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.5 $       $Date: 2009/06/22 19:45:49 $
 *
 ***************************************************************************/


#ifndef _VASPPLUGIN_H_
#define _VASPPLUGIN_H_

#include <stdio.h>
#include <math.h>
#include "molfile_plugin.h"

#define LINESIZE 1024
#define MAXATOMTYPES 100

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

typedef struct {
  FILE *file;
  char *filename;
  char *titleline;              /* Comment line provides system information */
  int version;                  /* VASP version control */
  int numatoms;                 /* total number of atoms */
  int eachatom[MAXATOMTYPES];   /* number of atoms per atom type */
  molfile_atom_t *atomlist;

  float cell[3][3];             /* lattice vectors of the unit cell */
  float rotmat[3][3];           /* rotation matrix, stored for periodic display hack */

  /* volumetric variables for charge density data */
  int nvolsets;                 /* number of volumetric datasets */
  molfile_volumetric_t *vol;    /* volume set metadata */
} vasp_plugindata_t;


/* allocate memory for plugin data and NULLify the pointers */
static vasp_plugindata_t *vasp_plugindata_malloc()
{
  vasp_plugindata_t *data = (vasp_plugindata_t *)malloc(sizeof(vasp_plugindata_t));

  if (!data) {
    fprintf(stderr, "\n\nVASP plugin) ERROR: cannot allocate memory for plugin data.\n");
    return NULL;
  }

  data->file = NULL;
  data->filename = NULL;
  data->titleline = NULL;
  data->atomlist = NULL;
  data->vol = NULL;

  return data;
}


/* free up the plugin data */
static void vasp_plugindata_free(vasp_plugindata_t *data)
{
  if (!data) return;

  if (data->file) fclose(data->file);
  if (data->filename) free(data->filename);
  if (data->titleline) free(data->titleline);
  if (data->atomlist) free(data->atomlist);
  if (data->vol) free(data->vol);
  free(data);
  data = NULL;
}


/* calculate and store rotation matrix to realign everything later. */
static void vasp_buildrotmat(vasp_plugindata_t *data)
{
  float const *const a = data->cell[0];
  float const *const b = data->cell[1];

  /* rotate first around y and z to align a along the x-axis... */
  const double len   = sqrt(a[0]*a[0] + a[1]*a[1]);
  const double phi   = atan2((double) a[2], (double) len);
  const double theta = atan2((double) a[1], (double) a[0]);

  const double cph = cos(phi);
  const double cth = cos(theta);
  const double sph = sin(phi);
  const double sth = sin(theta);

  /* ...then rotate around x to put b into the xy-plane. */
  const double psi = atan2(-sph*cth*b[0] - sph*sth*b[1] + cph*b[2],-sth*b[0] + cth*b[1]);
  const double cps = cos(psi);
  const double sps = sin(psi);

  data->rotmat[0][0] =  cph*cth;
  data->rotmat[0][1] =  cph*sth;
  data->rotmat[0][2] =  sph;
  data->rotmat[1][0] = -sth*cps - sph*cth*sps;
  data->rotmat[1][1] =  cth*cps - sph*sth*sps;
  data->rotmat[1][2] =  cph*sps; 
  data->rotmat[2][0] =  sth*sps - sph*cth*cps;
  data->rotmat[2][1] = -cth*sps - sph*sth*cps; 
  data->rotmat[2][2] =  cph*cps;
}


static void vasp_timestep_unitcell(molfile_timestep_t *ts, vasp_plugindata_t *data)
{
  if (!ts || !data) return;

  ts->A = sqrt(data->cell[0][0]*data->cell[0][0]+data->cell[0][1]*data->cell[0][1]+data->cell[0][2]*data->cell[0][2]);
  ts->B = sqrt(data->cell[1][0]*data->cell[1][0]+data->cell[1][1]*data->cell[1][1]+data->cell[1][2]*data->cell[1][2]);
  ts->C = sqrt(data->cell[2][0]*data->cell[2][0]+data->cell[2][1]*data->cell[2][1]+data->cell[2][2]*data->cell[2][2]);

  ts->gamma = acos((data->cell[0][0]*data->cell[1][0]+data->cell[0][1]*data->cell[1][1]+data->cell[0][2]*data->cell[1][2])/(ts->A*ts->B))*180.0/M_PI;
  ts->beta  = acos((data->cell[0][0]*data->cell[2][0]+data->cell[0][1]*data->cell[2][1]+data->cell[0][2]*data->cell[2][2])/(ts->A*ts->C))*180.0/M_PI;
  ts->alpha = acos((data->cell[1][0]*data->cell[2][0]+data->cell[1][1]*data->cell[2][1]+data->cell[1][2]*data->cell[2][2])/(ts->B*ts->C))*180.0/M_PI;
}

#endif /* _VASPPLUGIN_H_ */
