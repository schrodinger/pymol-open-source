/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_netcdfplugin
#define STATIC_PLUGIN 1

/***************************************************************************
 *cr
 *cr            (C) Copyright 1995-2016 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: netcdfplugin.c,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.29 $       $Date: 2016/11/28 05:01:54 $
 *
 ***************************************************************************/

/*
 * NetCDF based trajectories, used by AMBER 9, MMTK, etc.
 */

#include <netcdf.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "molfile_plugin.h"

#define CDF_TYPE_UNKNOWN 0
#define CDF_TYPE_AMBER   1
#define CDF_TYPE_MMTK    2

#define CDF_SUCCESS      0
#define CDF_ERR         -1

typedef struct {
  int trajectorytype;
  int step_numberdimid;
  size_t step_numberdim;
  int minor_step_numberdimid;
  size_t minor_step_numberdim;
  int atom_numberdimid;
  size_t atom_numberdim;
  int xyzdimid;
  size_t xyzdim;
  int box_size_lengthdimid;
  size_t box_size_lengthdim;
  int description_lengthdimid;
  size_t description_lengthdim;
  char *description;
  int description_id;
  int step_id;
  int time_id;
  int box_size_id;
  int configuration_id;
  int has_box;
  char *comment;
} mmtkdata;

typedef struct {
  int is_restart;
  int is_trajectory;
  int has_box;
  int atomdimid;
  size_t atomdim;
  int spatialdimid;
  size_t spatialdim;
  int framedimid;
  size_t framedim;
  char *conventionversion;
  char *title;
  char *application;
  char *program;
  char *programversion;
  int spatial_id;                  /* "xyz" */
  int cell_spatial_id;             /* "abc" */
  int cell_angular_id;             /* "alpha, beta, gamma" */
  int time_id;                     /* frame time in picoseconds */
  int coordinates_id;              /* coords in angstroms */
  char *coordinates_units;         /* coordinates units */
  float coordinates_scalefactor;   /* coordinates scaling factor */
  int cell_lengths_id;             /* cell lengths in angstroms */
  char *cell_lengths_units;        /* cell lengths units */
  float cell_lengths_scalefactor;  /* cell lengths scaling factor */
  int cell_angles_id;              /* cell angles in degrees */
  char *cell_angles_units;         /* cell angles units */
  float cell_angles_scalefactor;   /* cell angles scaling factor */
  int velocities_id;               /* velocities in angstroms/picosecond */
} amberdata;


typedef struct {
  /* sub-format independent data */
  int ncid;
  int type;
  int natoms; 
  int curframe;
  char *conventions;

  /* stuff used by AMBER */
  amberdata amber;

  /* stuff used by MMTK */
  mmtkdata mmtk;

} cdfdata;


static void close_cdf_read(void *mydata) {
  cdfdata *cdf = (cdfdata *)mydata;

  nc_close(cdf->ncid);

  /* AMBER stuff */
  if (cdf->amber.title)
    free(cdf->amber.title);

  if (cdf->amber.application)
    free(cdf->amber.application);

  if (cdf->amber.program)
    free(cdf->amber.program);

  if (cdf->amber.programversion)
    free(cdf->amber.programversion);

  if (cdf->amber.conventionversion)
    free(cdf->amber.conventionversion);

  if (cdf->amber.coordinates_units)
    free(cdf->amber.coordinates_units);

  if (cdf->amber.cell_lengths_units)
    free(cdf->amber.cell_lengths_units);

  /* MMTK stuff */
  if (cdf->mmtk.comment)
    free(cdf->mmtk.comment);

  /* format independent stuff */
  if (cdf->conventions)
    free(cdf->conventions);

  free(cdf);
}



static int open_amber_cdf_read(cdfdata *cdf) {
  int rc;
  size_t len; 
  amberdata *amber = &cdf->amber;

  /* check if this is a restart file or not */
  if (!strcmp(cdf->conventions, "AMBERRESTART")) {
    amber->is_restart = 1;
  } else {
    amber->is_trajectory = 1;
  }

  /* global attrib: "ConventionVersion" -- required */
  rc = nc_inq_attlen(cdf->ncid, NC_GLOBAL, "ConventionVersion", &len);
  if (rc == NC_NOERR && len > 0) {
    amber->conventionversion = (char *) malloc((len+1) * sizeof(char));
    nc_get_att_text(cdf->ncid, NC_GLOBAL, "ConventionVersion", amber->conventionversion);
    amber->conventionversion[len] = '\0';
    printf("netcdfplugin) %s follows AMBER conventions version '%s'\n", 
           (amber->is_restart) ? "restart file" : "trajectory",
           amber->conventionversion);
  } else {
    return CDF_ERR;
  }

  /* at this point we know that this is an AMBER trajectory or restart file */
  cdf->type = CDF_TYPE_AMBER;

  /* initialize default scaling factors so they are always set to a sane */
  /* value even if problems occur later */
  amber->coordinates_scalefactor = 1.0;
  amber->cell_lengths_scalefactor = 1.0;
  amber->cell_angles_scalefactor = 1.0;

  /* global attrib: "program" -- required */
  rc = nc_inq_attlen(cdf->ncid, NC_GLOBAL, "program", &len);
  if (rc == NC_NOERR && len > 0) {
    amber->program = (char *) malloc((len+1) * sizeof(char));
    nc_get_att_text(cdf->ncid, NC_GLOBAL, "program", amber->program);
    amber->program[len] = '\0';
    printf("netcdfplugin) AMBER: program '%s'\n", amber->program);
  } else {
    printf("netcdfplugin) AMBER: Missing required 'program' global attribute, corrupt file?\n");
  }


  /* global attrib: "programVersion" -- required */
  rc = nc_inq_attlen(cdf->ncid, NC_GLOBAL, "programVersion", &len);
  if (rc == NC_NOERR && len > 0) {
    amber->programversion = (char *) malloc((len+1) * sizeof(char));
    nc_get_att_text(cdf->ncid, NC_GLOBAL, "programVersion", amber->programversion);
    amber->programversion[len] = '\0';
    printf("netcdfplugin) AMBER: program version '%s'\n", amber->programversion);
  } else {
    printf("netcdfplugin) AMBER: Missing required 'programVersion' global attribute, corrupt file?\n");
  }


  /* global attrib: "title" -- optional */
  rc = nc_inq_attlen(cdf->ncid, NC_GLOBAL, "title", &len);
  if (rc == NC_NOERR && len > 0) {
    amber->title = (char *) malloc((len+1) * sizeof(char));
    nc_get_att_text(cdf->ncid, NC_GLOBAL, "title", amber->title);
    amber->title[len] = '\0';
    printf("netcdfplugin) AMBER: title '%s'\n", amber->title);
  } 


  /* global attrib: "application" -- optional */
  rc = nc_inq_attlen(cdf->ncid, NC_GLOBAL, "application", &len);
  if (rc == NC_NOERR && len > 0) {
    amber->application = (char *) malloc((len+1) * sizeof(char));
    nc_get_att_text(cdf->ncid, NC_GLOBAL, "application", amber->application);
    amber->application[len] = '\0';
    printf("netcdfplugin) AMBER: application '%s'\n", amber->application);
  } 


/* XXX lots of additional error checking is needed below... */

  /* read in spatial dimension */
  rc = nc_inq_dimid(cdf->ncid, "spatial", &amber->spatialdimid);
  if (rc == NC_NOERR) {    
    rc = nc_inq_dimlen(cdf->ncid, amber->spatialdimid, &amber->spatialdim);
    if (rc == NC_NOERR) {
      printf("netcdfplugin) AMBER: spatial dimension: %ld\n", (long)amber->spatialdim);
    } else {
      printf("netcdfplugin) AMBER: Missing spatial dimension, corrupt file?\n");
      printf("netcdfplugin) AMBER: Fixing by guessing spatialdim as '3'\n");
      amber->spatialdim = 3;
    }
  } else {
    printf("netcdfplugin) AMBER: Missing spatial dimension, corrupt file?\n");
    printf("netcdfplugin) AMBER: Fixing by guessing spatialdim as '3'\n");
    amber->spatialdim = 3;
  }
 
  /* read in atom dimension */
  rc = nc_inq_dimid(cdf->ncid, "atom", &amber->atomdimid);
  if (rc == NC_NOERR) {    
    rc = nc_inq_dimlen(cdf->ncid, amber->atomdimid, &amber->atomdim);
    if (rc == NC_NOERR) {
      printf("netcdfplugin) AMBER: atom dimension: %ld\n", (long)amber->atomdim);
      cdf->natoms = amber->atomdim; /* copy to format independent part */
    } else  {
      printf("netcdfplugin) AMBER: missing atom dimension, aborting\n");
      return CDF_ERR;
    }
  } else {
    printf("netcdfplugin) AMBER: missing atom dimension, aborting\n");
    return CDF_ERR;
  }
 
  /* if this is a trajectory, read in frame dimension */
  if (amber->is_trajectory) {
    rc = nc_inq_dimid(cdf->ncid, "frame", &amber->framedimid);
    if (rc == NC_NOERR) {    
      rc = nc_inq_dimlen(cdf->ncid, amber->framedimid, &amber->framedim);
      if (rc == NC_NOERR) {
        printf("netcdfplugin) AMBER: frame dimension: %ld\n", (long)amber->framedim);
      } else {
        printf("netcdfplugin) AMBER: missing frame dimension, aborting\n");
        return CDF_ERR;
      }
    } else {
      printf("netcdfplugin) AMBER: missing frame dimension, aborting\n");
      return CDF_ERR;
    }
  }


  /* 
   * get ID values for all of the variables we're interested in 
   */
#if 0
  /* VMD can live without the various human readable label variables. */
  rc = nc_inq_varid(cdf->ncid, "spatial", &amber->spatial_id);
  if (rc != NC_NOERR)
    return CDF_ERR;

  rc = nc_inq_varid(cdf->ncid, "cell_spatial", &amber->cell_spatial_id);
  if (rc != NC_NOERR)
    return CDF_ERR;

  rc = nc_inq_varid(cdf->ncid, "cell_angular", &amber->cell_angular_id);
  if (rc != NC_NOERR)
    return CDF_ERR;
#endif

  /* VMD requires coordinates at a minimum */
  rc = nc_inq_varid(cdf->ncid, "coordinates", &amber->coordinates_id);
  if (rc != NC_NOERR) {
    printf("netcdfplugin) AMBER: no coordinates variable, nothing to load\n");
    return CDF_ERR;
  }

  /* Coordinate units */
  rc = nc_inq_attlen(cdf->ncid, amber->coordinates_id, "units", &len);
  if (rc == NC_NOERR && len > 0) {
    amber->coordinates_units = (char *) malloc((len+1) * sizeof(char));
    nc_get_att_text(cdf->ncid, amber->coordinates_id, "units", amber->coordinates_units);
    amber->coordinates_units[len] = '\0';
    printf("netcdfplugin) AMBER: coordinates units: '%s'\n", amber->coordinates_units);
  } else {
    printf("netcdfplugin) AMBER: no coordinates units attribute, Angstroms assumed\n");
  }

  /* Coordinate scaling factor to get to Angstroms */
  if (nc_get_att_float(cdf->ncid, amber->coordinates_id, "scale_factor", &amber->coordinates_scalefactor) != NC_NOERR) {
    printf("netcdfplugin) AMBER: no coordinates scalefactor attribute, 1.0 assumed\n");
  }
  printf("netcdfplugin) AMBER: coordinates scalefactor: %f\n", amber->coordinates_scalefactor);

#if 0
  /* we don't need velocities at this time */
  rc = nc_inq_varid(cdf->ncid, "velocities", &amber->velocities_id);
  if (rc != NC_NOERR) {
    printf("netcdfplugin) AMBER: missing velocities variable, aborting\n");
    return CDF_ERR;
  }
#endif

  /* optional periodic cell info */
  rc = nc_inq_varid(cdf->ncid, "cell_lengths", &amber->cell_lengths_id);
  if (rc == NC_NOERR) {
    rc = nc_inq_varid(cdf->ncid, "cell_angles", &amber->cell_angles_id);
    if (rc == NC_NOERR) {
      printf("netcdfplugin) AMBER trajectory contains periodic cell information\n");
      amber->has_box = 1;

      /* Cell lengths units */
      rc = nc_inq_attlen(cdf->ncid, amber->cell_lengths_id, "units", &len);
      if (rc == NC_NOERR && len > 0) {
        amber->cell_lengths_units = (char *) malloc((len+1) * sizeof(char));
        nc_get_att_text(cdf->ncid, amber->cell_lengths_id, "units", amber->cell_lengths_units);
        amber->cell_lengths_units[len] = '\0';
        printf("netcdfplugin) AMBER: cell lengths units: '%s'\n", amber->cell_lengths_units);
      } else {
        printf("netcdfplugin) AMBER: no cell lengths units attribute, Angstroms assumed\n");
      }

      /* Cell lengths scaling factor to get to Angstroms */
      if (nc_get_att_float(cdf->ncid, amber->cell_lengths_id, "scale_factor", &amber->cell_lengths_scalefactor) != NC_NOERR) {
        printf("netcdfplugin) AMBER: no cell lengths scalefactor attribute, 1.0 assumed\n");
      }
      printf("netcdfplugin) AMBER: cell lengths scalefactor: %f\n", amber->cell_lengths_scalefactor);

      /* Cell angles units */
      rc = nc_inq_attlen(cdf->ncid, amber->cell_angles_id, "units", &len);
      if (rc == NC_NOERR && len > 0) {
        amber->cell_angles_units = (char *) malloc((len+1) * sizeof(char));
        nc_get_att_text(cdf->ncid, amber->cell_angles_id, "units", amber->cell_angles_units);
        amber->cell_angles_units[len] = '\0';
        printf("netcdfplugin) AMBER: cell angles units: '%s'\n", amber->cell_angles_units);
      } else {
        printf("netcdfplugin) AMBER: no cell angles units attribute, Degrees assumed\n");
      }

      /* Cell angles scaling factor to get to degrees */
      if (nc_get_att_float(cdf->ncid, amber->cell_angles_id, "scale_factor", &amber->cell_angles_scalefactor) != NC_NOERR) {
        printf("netcdfplugin) AMBER: no cell angles scalefactor attribute, 1.0 assumed\n");
      }
      printf("netcdfplugin) AMBER: cell angles scalefactor: %f\n", amber->cell_angles_scalefactor);
    }
  }

  return CDF_SUCCESS;
}


static int open_mmtk_cdf_read(cdfdata *cdf, int conventionsknown) {
  int rc;
  size_t len; 
  mmtkdata *mmtk = &cdf->mmtk;

  /* If conventions specify MMTK then we're safe to continue */
  /* and we know what we're dealing with */
  if (conventionsknown) {
    cdf->type = CDF_TYPE_MMTK;
  }

  /* global attrib: "trajectory_type" (new format) */
  rc = nc_get_att_int(cdf->ncid, NC_GLOBAL, "trajectory_type", &mmtk->trajectorytype);
  if (rc == NC_NOERR) {
    printf("netcdfplugin) MMTK trajectory type: %d\n", mmtk->trajectorytype);
  } else {
    printf("netcdfplugin) Assuming MMTK trajectory type: %d\n", mmtk->trajectorytype);
    mmtk->trajectorytype = 0;
  }

  /* read in spatial dimension */
  rc = nc_inq_dimid(cdf->ncid, "xyz", &mmtk->xyzdimid);
  if (rc == NC_NOERR) {
    rc = nc_inq_dimlen(cdf->ncid, mmtk->xyzdimid, &mmtk->xyzdim);
    if (rc == NC_NOERR)
      printf("netcdfplugin) MMTK: xyz dimension: %ld\n", (long)mmtk->xyzdim);
    else 
      return CDF_ERR;
  } else {
    return CDF_ERR;
  }

  /* read in atom dimension */
  rc = nc_inq_dimid(cdf->ncid, "atom_number", &mmtk->atom_numberdimid); 
  if (rc == NC_NOERR) {
    rc = nc_inq_dimlen(cdf->ncid, mmtk->atom_numberdimid, &mmtk->atom_numberdim);
    if (rc == NC_NOERR) {
      printf("netcdfplugin) MMTK: atom_number dimension: %ld\n", (long)mmtk->atom_numberdim);
      cdf->natoms = mmtk->atom_numberdim; /* copy to format independent part */
    } else {
      return CDF_ERR;
    }
  } else {
    return CDF_ERR;
  }


  /* read in frame dimension */
  rc = nc_inq_dimid(cdf->ncid, "step_number", &mmtk->step_numberdimid);
  if (rc == NC_NOERR) {
    rc = nc_inq_dimlen(cdf->ncid, mmtk->step_numberdimid, &mmtk->step_numberdim);
    if (rc == NC_NOERR)
      printf("netcdfplugin) MMTK: step_number dimension: %ld\n", (long)mmtk->step_numberdim);
    else 
      return CDF_ERR;
  } else {
    return CDF_ERR;
  }


  /* read in minor step number dimension */
  rc = nc_inq_dimid(cdf->ncid, "minor_step_number", &mmtk->minor_step_numberdimid);
  if (rc == NC_NOERR) {
    rc = nc_inq_dimlen(cdf->ncid, mmtk->minor_step_numberdimid, &mmtk->minor_step_numberdim);
    if (rc == NC_NOERR)
      printf("netcdfplugin) MMTK: minor_step_number dimension: %ld\n", (long)mmtk->minor_step_numberdim);
    else 
      return CDF_ERR;
  } else if (rc == NC_EBADDIM) {
    printf("netcdfplugin) MMTK: no minor_step_number dimension\n");
    mmtk->minor_step_numberdim = 0;
  } else {
    return CDF_ERR;
  }


  /* read in description_length dimension */
  rc = nc_inq_dimid(cdf->ncid, "description_length", &mmtk->description_lengthdimid); 
  if (rc == NC_NOERR) {
    rc = nc_inq_dimlen(cdf->ncid, mmtk->description_lengthdimid, &mmtk->description_lengthdim);
    if (rc == NC_NOERR)
      printf("netcdfplugin) MMTK: description_length dimension: %ld\n", (long)mmtk->description_lengthdim);
    else
      return CDF_ERR;
  } else {
    return CDF_ERR;
  }


  /* get ID values for all of the variables we're interested in */
  rc = nc_inq_varid(cdf->ncid, "configuration", &mmtk->configuration_id);
  if (rc != NC_NOERR)
    return CDF_ERR;

  rc = nc_inq_varid(cdf->ncid, "description", &mmtk->description_id);
  if (rc != NC_NOERR)
    return CDF_ERR;

  /* check for PBC */
  rc = nc_inq_varid(cdf->ncid, "box_size", &mmtk->box_size_id);
  if (rc == NC_NOERR) {
    mmtk->has_box = 1;
    printf("netcdfplugin) MMTK: system has periodic boundary conditions\n");
  }
  else if (rc == NC_ENOTVAR)
    mmtk->has_box = 0;
  else
    return CDF_ERR;


  /* global attrib: "comment" -- optional */
  rc = nc_inq_attlen(cdf->ncid, NC_GLOBAL, "comment", &len);
  if (rc == NC_NOERR && len > 0) {
    mmtk->comment = (char *) malloc((len+1) * sizeof(char));
    nc_get_att_text(cdf->ncid, NC_GLOBAL, "comment", mmtk->comment);
    mmtk->comment[len] = '\0';
    printf("netcdfplugin) MMTK: comment '%s'\n", mmtk->comment);
  } 

  /* at this point we know that this is an MMTK trajectory */
  if (!conventionsknown) {
    printf("netcdfplugin) File is an old format MMTK trajectory without conventions\n");    
    cdf->type = CDF_TYPE_MMTK;
  }

  return CDF_SUCCESS;
}

 
static void *open_cdf_read(const char *filename, const char *filetype, 
                           int *natoms) {
  int ncid, rc;
  size_t len;
  cdfdata *cdf;
 
  rc = nc_open(filename, NC_NOWRITE, &ncid);
  if (rc != NC_NOERR) return NULL;

  cdf = (cdfdata *) malloc(sizeof(cdfdata));
  memset(cdf, 0, sizeof(cdfdata));

  cdf->ncid = ncid;
  cdf->type = CDF_TYPE_UNKNOWN;

  /* Determine what NetCDF conventions apply to this data, if any */
  rc = nc_inq_attlen(cdf->ncid, NC_GLOBAL, "Conventions", &len);
  if (rc == NC_NOERR && len > 0) {
    cdf->conventions = (char *) malloc((len+1) * sizeof(char));
    nc_get_att_text(cdf->ncid, NC_GLOBAL, "Conventions", cdf->conventions);
    cdf->conventions[len] = '\0';
    printf("netcdfplugin) conventions: '%s'\n", cdf->conventions);
  } 

  if (cdf->conventions != NULL) {
    /* Check if this is a file generated by AMBER */
    if (strstr(cdf->conventions, "AMBER") != NULL) {
      if (!open_amber_cdf_read(cdf)) {
        *natoms = cdf->natoms;
        return cdf;
      }
    } 

    /* Check if this is a file generated by MMTK */
    if (strstr(cdf->conventions, "MMTK") != NULL) {
      if (!open_mmtk_cdf_read(cdf, 1)) {
        *natoms = cdf->natoms;
        return cdf;
      }
    } 
  } 

  printf("netcdfplugin) Missing or unrecognized conventions attribute\n");
  printf("netcdfplugin) checking for old format MMTK NetCDF file...\n");

  /* If no conventions are specified, then maybe it's from MMTK */
  if (!open_mmtk_cdf_read(cdf, 0)) {
    *natoms = cdf->natoms;
    return cdf;
  } 

  /* if no conventions are recognized, then we free everything */
  /* and return failure                                        */
  close_cdf_read(cdf);

  return NULL; 
}

/* A very basic bracket counter. It assumes that the expression
   is syntactically correct. */
static char *find_closing_bracket(char *s) {
  int count = 1;
  while (*s && count > 0) {
    if (*s == '(' || *s == '[')
      count++;
    if (*s == ')' || *s == ']')
      count--;
    s++;
  }
  return s;
}

/* Simple string replacement routine for fixing atom names. */
static void atom_name_replace(char *name, char *substring, char letter) {
  char *s = strstr(name, substring);
  if (s != NULL) {
    *s = letter;
    strcpy(s+1, s+strlen(substring));
  }
}

static void atom_name_remove_underscores(char *name) {
  char *s = name;
  while (1) {
    s = strchr(s, '_');
    if (s == NULL)
      break;
    strcpy(s, s+1);
  }
}

/* Set chainid, resname, and resnum for a range of atoms
   and fix atom names. */
static void set_atom_attributes(molfile_atom_t *atoms, int natoms,
				char **atom_pointers, char chain_id,
				char *resname, int resnum,
				char *start, char *end,
				int name_correction_type) {
  int i;
  for (i=0; i<natoms; i++)
    if (atom_pointers[i] > start && atom_pointers[i] < end) {
      molfile_atom_t *atom = atoms + i;
      atom->chain[0] = chain_id;      
      atom->chain[1] = '\0';      
      strcpy(atom->resname, resname);
      atom->resid = resnum;
      if (name_correction_type == 1 /* proteins */) {
	atom_name_replace(atom->name, "_alpha", 'A');
	atom_name_replace(atom->name, "_beta", 'B');
	atom_name_replace(atom->name, "_gamma", 'G');
	atom_name_replace(atom->name, "_delta", 'D');
	atom_name_replace(atom->name, "_epsilon", 'E');
	atom_name_replace(atom->name, "_zeta", 'Z');
	atom_name_replace(atom->name, "_eta", 'H');
	atom_name_remove_underscores(atom->name);
      }
      else if (name_correction_type == 2 /* nucleic acids */) {
	if (strcmp(atom->name, "O_1") == 0)
	  strcpy(atom->name, "O1P");
	else if (strcmp(atom->name, "O_2") == 0)
	  strcpy(atom->name, "O2P");
	else if (strcmp(atom->name, "C_1") == 0)
	  strcpy(atom->name, "C1'");
	else if (strcmp(atom->name, "C_2") == 0)
	  strcpy(atom->name, "C2'");
	else if (strcmp(atom->name, "C_3") == 0)
	  strcpy(atom->name, "C3'");
	else if (strcmp(atom->name, "O_3") == 0)
	  strcpy(atom->name, "O3'");
	else if (strcmp(atom->name, "C_4") == 0)
	  strcpy(atom->name, "C4'");
	else if (strcmp(atom->name, "O_4") == 0)
	  strcpy(atom->name, "O4'");
	else if (strcmp(atom->name, "C_5") == 0)
	  strcpy(atom->name, "C5'");
	else if (strcmp(atom->name, "O_5") == 0)
	  strcpy(atom->name, "O5'");
	else
	  atom_name_remove_underscores(atom->name);
      }
    }
}

/* Get structure from an MMTK trajectory file */
static int read_mmtk_cdf_structure(void *mydata, int *optflags,
                                   molfile_atom_t *atoms) {
  int i, rc;
  molfile_atom_t *atom;
  cdfdata *cdf = (cdfdata *) mydata;
  mmtkdata *mmtk = &cdf->mmtk;
  size_t start[3], count[3];
  char *dstr;
  char **atom_pointers;
  int resnum;
  char resname[8];

  *optflags = MOLFILE_NOOPTIONS;

  mmtk->description = (char *) malloc((mmtk->description_lengthdim + 1) * sizeof(char));
  if (mmtk->description == NULL) 
    return MOLFILE_ERROR;

  start[0] = cdf->curframe; /* frame */
  count[0] = mmtk->description_lengthdim;

  rc = nc_get_vara_text(cdf->ncid, mmtk->description_id,
                        start, count, mmtk->description);
  if (rc != NC_NOERR)
    return MOLFILE_ERROR;

  /* initialize all atoms with name "X" to start with */
  /* indicating unknown atom types etc..              */
  for (i=0; i<cdf->natoms; i++) {
    atom = atoms + i;
    strncpy(atom->name, "X", sizeof(atom->name)-1);
    atom->name[sizeof(atom->name)-1] = '\0';
    strncpy(atom->type, atom->name, sizeof(atom->type)-1);
    atom->type[sizeof(atom->type)-1] = '\0';
    atom->resname[0] = '\0';
    atom->resid = 1;
    atom->chain[0] = '\0';
    atom->segid[0] = '\0';
  }

  /* Allocate a pointer array that will hold each atom's location in
     the description string. This will be used in a second pass through
     the description string in which residue names and indices will
     be assigned. */
  atom_pointers = (char **) malloc(cdf->natoms * sizeof(char *));
  if (atom_pointers == NULL)
    return MOLFILE_ERROR;

  /* First pass: look only at atoms */
  dstr = mmtk->description;
  while (dstr < (mmtk->description + mmtk->description_lengthdim)) {
    char *atomstr;
    atomstr = strstr(dstr, "A('");

    if (atomstr != NULL) {
      char name[1024];
      char *nmstart = NULL;
      char *nmend = NULL;
      char *indstart = NULL;
      char *endp = NULL;
      int index, len;

      endp = strchr(atomstr, ')');
      nmstart = strchr(atomstr, '\'');
      if (nmstart != NULL)
        nmend = strchr(nmstart+1, '\'');
      indstart = strchr(atomstr, ',');
      if (endp == NULL || nmstart == NULL || nmend == NULL || indstart == NULL) {
        printf("netcdfplugin) mmtk_read_structure(): unable to parse atom tag\n");
        break; /* something went wrong */
      }

      len = nmend - nmstart - 1;
      if (len > sizeof(name)) {
        printf("netcdfplugin) mmtk_read_structure(): bad length: %d\n", len);
        break; /* something went wrong */
      }
      memcpy(name, nmstart+1, len); 
      name[len] = '\0';

      index = -1;
      sscanf(indstart, ",%d)", &index);
      atom_pointers[index] = atomstr;

      if (index >= 0 && index < cdf->natoms) {
        atom = atoms + index;
        strncpy(atom->name, name, sizeof(atom->name)-1);
        atom->name[sizeof(atom->name)-1] = '\0';
        strncpy(atom->type, atom->name, sizeof(atom->type)-1);
        atom->type[sizeof(atom->type)-1] = '\0';
      }

      dstr = atomstr+1;
    } else {
      break; /* no more atom records found */
    }
  }

  /* Second pass: peptide chains */
  dstr = mmtk->description;
  while (dstr < (mmtk->description + mmtk->description_lengthdim)) {
    char *peptide, *pend;
    char *group, *gend;
    char *nmstart, *nmend;
    char chain_id = 'A';
    char *s;

    peptide = strstr(dstr, "S('");
    if (peptide == NULL)
      break;
    pend = find_closing_bracket(peptide+2);

    resnum = 1;
    group = peptide;
    while (1) {
      group = strstr(group, "G('");
      if (group == NULL || group >= pend)
	break;
      gend = find_closing_bracket(group+2);
      nmstart = strchr(group, '\'') + 1;
      nmend = strchr(nmstart, '\'');
      while (nmend > nmstart && isdigit(*(nmend-1)))
	nmend--;
      if (nmend-nmstart > 7)
	nmend = nmstart+7;
      strncpy(resname, nmstart, nmend-nmstart);
      resname[nmend-nmstart] = '\0';
      s = resname;
      while (*s) {
	*s = toupper(*s);
	s++;
      }
      set_atom_attributes(atoms, cdf->natoms, atom_pointers,
			  chain_id, resname, resnum, group, gend, 1);
      group = gend;
      resnum++;
    }

    if (chain_id == 'Z')
      chain_id = 'A';
    else
	chain_id++;
    dstr = pend;
  }

  /* Third pass: nucleic acid chains */
  dstr = mmtk->description;
  while (dstr < (mmtk->description + mmtk->description_lengthdim)) {
    char *nacid, *nend;
    char *group, *gend;
    char *nmstart, *nmend;
    char chain_id = 'a';
    char *s;

    nacid = strstr(dstr, "N('");
    if (nacid == NULL)
      break;
    nend = find_closing_bracket(nacid+2);

    resnum = 1;
    group = nacid;
    while (1) {
      group = strstr(group, "G('");
      if (group == NULL || group >= nend)
	break;
      gend = find_closing_bracket(group+2);
      nmstart = strchr(group, '\'') + 1;
      nmend = strchr(nmstart, '\'');
      while (nmend > nmstart && isdigit(*(nmend-1)))
	nmend--;
      if (nmend > nmstart && nmend[-1] == '_')
	nmend--;
      if (nmend-nmstart > 7)
	nmend = nmstart+7;
      strncpy(resname, nmstart, nmend-nmstart);
      resname[nmend-nmstart] = '\0';
      s = resname;
      while (*s) {
	*s = toupper(*s);
	s++;
      }
      if (resname[0] == 'R' || resname[0] == 'D') {
	switch (resname[1]) {
	case 'A':
	  strcpy(resname, "ADE");
	  break;
	case 'C':
	  strcpy(resname, "CYT");
	  break;
	case 'G':
	  strcpy(resname, "GUA");
	  break;
	case 'T':
	  strcpy(resname, "THY");
	  break;
	case 'U':
	  strcpy(resname, "URA");
	  break;
	}
      }
      set_atom_attributes(atoms, cdf->natoms, atom_pointers,
			  chain_id, resname, resnum, group, gend, 2);
      group = gend;
      resnum++;
    }

    if (chain_id == 'z')
      chain_id = 'a';
    else
	chain_id++;
    dstr = nend;
  }

  /* Fourth pass: non-chain molecules */
  resnum = 1;
  dstr = mmtk->description;
  while (dstr < (mmtk->description + mmtk->description_lengthdim)) {
    char *molecule, *mend;
    char *nmstart, *nmend;

    molecule = strstr(dstr, "M('");
    if (molecule == NULL)
      break;
    mend = find_closing_bracket(molecule+2);
    nmstart = strchr(molecule, '\'') + 1;
    nmend = strchr(nmstart, '\'');
    if (strncmp(nmstart, "water", 5) == 0)
      strcpy(resname, "HOH");
    else {
      if (nmend-nmstart > 7)
	nmend = nmstart+7;
      strncpy(resname, nmstart, nmend-nmstart);
      resname[nmend-nmstart] = '\0';
    }
    set_atom_attributes(atoms, cdf->natoms, atom_pointers,
			'_', resname, resnum, molecule, mend, 0);
    resnum++;
    dstr = mend;
  }

  free(atom_pointers);

  return MOLFILE_SUCCESS;
}


static int read_cdf_structure(void *mydata, int *optflags,
                                   molfile_atom_t *atoms) {
  cdfdata *cdf = (cdfdata *)mydata;

  switch (cdf->type) {
    case CDF_TYPE_AMBER:
      return MOLFILE_NOSTRUCTUREDATA; /* not an error, just no data */

    case CDF_TYPE_MMTK:
      return read_mmtk_cdf_structure(mydata, optflags, atoms);
  }

  return MOLFILE_ERROR;
}


static int read_amber_cdf_timestep(void *mydata, int natoms, molfile_timestep_t *ts) {
  size_t start[3], count[3];
  cdfdata *cdf = (cdfdata *)mydata;
  amberdata *amber = &cdf->amber;
  int rc;

  /* Read in the atom coordinates and unit cell information */
  /* only save coords if we're given a valid ts pointer     */ 
  /* otherwise VMD wants us to skip it.                     */
  if (ts != NULL) {
    if (amber->is_restart) {
      // restart files only contain one frame, so return if we get called again
      if (cdf->curframe > 0)
        return MOLFILE_ERROR;

      start[0] = 0;             /* atom */
      start[1] = 0;             /* spatial */
      start[2] = 0;

      count[0] = amber->atomdim;
      count[1] = amber->spatialdim;
      count[2] = 0;

      rc = nc_get_vara_float(cdf->ncid, amber->coordinates_id, 
                             start, count, ts->coords);
      if (rc != NC_NOERR) { 
        printf("netcdfplugin) AMBER: failed to parse restart file coordinates!\n");
      }
    } else {
      start[0] = cdf->curframe; /* frame */
      start[1] = 0;             /* atom */
      start[2] = 0;             /* spatial */

      count[0] = 1;
      count[1] = amber->atomdim;
      count[2] = amber->spatialdim;

      /* parse trajectory timestep */
      rc = nc_get_vara_float(cdf->ncid, amber->coordinates_id, 
                             start, count, ts->coords);
    }

    if (rc != NC_NOERR) 
      return MOLFILE_ERROR;

    /* apply coordinate scaling factor if not 1.0 */
    if (amber->coordinates_scalefactor != 1.0) {
      int i;
      float s = amber->coordinates_scalefactor;
      for (i=0; i<natoms*3; i++) {
        ts->coords[i] *= s;
      }
    }

    /* Read the PBC box info. */
    if (amber->has_box) {
      size_t start[3], count[3];
      double lengths[3];
      double angles[3];

      if (amber->is_restart) {
        start[0] = 0;             /* spatial */
        start[1] = 0;
        start[2] = 0;
        count[0] = amber->spatialdim;
        count[1] = 0;
        count[2] = 0;
      } else {
        start[0] = cdf->curframe; /* frame */
        start[1] = 0;             /* spatial */
        start[2] = 0;
        count[0] = 1;
        count[1] = amber->spatialdim;
        count[2] = 0;
      }

      rc = nc_get_vara_double(cdf->ncid, amber->cell_lengths_id, 
                              start, count, lengths);
      if (rc != NC_NOERR) 
        return MOLFILE_ERROR;

      rc = nc_get_vara_double(cdf->ncid, amber->cell_angles_id, 
                              start, count, angles);
      if (rc != NC_NOERR) 
        return MOLFILE_ERROR;

      ts->A = lengths[0] * amber->cell_lengths_scalefactor;
      ts->B = lengths[1] * amber->cell_lengths_scalefactor;
      ts->C = lengths[2] * amber->cell_lengths_scalefactor;

      ts->alpha = angles[0] * amber->cell_angles_scalefactor;
      ts->beta  = angles[1] * amber->cell_angles_scalefactor;
      ts->gamma = angles[2] * amber->cell_angles_scalefactor;
    }
  }

  cdf->curframe++;
  return MOLFILE_SUCCESS;
}


static int read_mmtk_cdf_timestep(void *mydata, int natoms, molfile_timestep_t *ts) {
  cdfdata *cdf = (cdfdata *)mydata;
  mmtkdata *mmtk = &cdf->mmtk;
  int rc;

  /* Read in the atom coordinates and unit cell information */
  /* only save coords if we're given a valid ts pointer     */ 
  /* otherwise VMD wants us to skip it.                     */
  if (ts != NULL) {
    size_t start[4], count[4];
    int i;

    if (mmtk->minor_step_numberdim == 0) {
      start[0] = cdf->curframe; /* step */
      start[1] = 0;             /* atom */
      start[2] = 0;             /* spatial */
      start[3] = 0;             /* minor step */
    }
    else {
      start[0] = cdf->curframe/mmtk->minor_step_numberdim;   /* step */
      start[1] = 0;             /* atom */
      start[2] = 0;             /* spatial */
      start[3] = cdf->curframe % mmtk->minor_step_numberdim; /* minor step */
    }

    count[0] = 1;
    count[1] = mmtk->atom_numberdim;
    count[2] = mmtk->xyzdim;
    count[3] = 1;             /* only want one minor step, regardless */

    rc = nc_get_vara_float(cdf->ncid, mmtk->configuration_id, 
                           start, count, ts->coords);
    if (rc != NC_NOERR) 
      return MOLFILE_ERROR;

    /* check for allocated but not yet used frame */
    if (ts->coords[0] == NC_FILL_FLOAT)
      return MOLFILE_ERROR;

    /* scale coordinates from nanometers to angstroms */
    for (i=0; i<(3 * mmtk->atom_numberdim); i++) {
      ts->coords[i] *= 10.0f;
    }

    /* Read the PBC box info. */
    if (mmtk->has_box) {
      float lengths[3];

      if (mmtk->minor_step_numberdim == 0) {
	start[0] = cdf->curframe; /* step */
	start[1] = 0;             /* box_size */
	start[2] = 0;             /* minor step */
      }
      else {
	start[0] = cdf->curframe/mmtk->minor_step_numberdim;   /* step */
	start[1] = 0;             /* box_size */
	start[2] = cdf->curframe % mmtk->minor_step_numberdim; /* minor step */
      }

      count[0] = 1;
      count[1] = 3;
      count[2] = 1;

      rc = nc_get_vara_float(cdf->ncid, mmtk->box_size_id,
                             start, count, lengths);
      if (rc != NC_NOERR) 
        return MOLFILE_ERROR;

      ts->A = 10.*lengths[0];
      ts->B = 10.*lengths[1];
      ts->C = 10.*lengths[2];

      ts->alpha = 90.;
      ts->beta  = 90.;
      ts->gamma = 90.;
    }
  }

  cdf->curframe++;
  return MOLFILE_SUCCESS;
}



static int read_cdf_timestep(void *mydata, int natoms, molfile_timestep_t *ts) {
  cdfdata *cdf = (cdfdata *)mydata;

  switch (cdf->type) {
    case CDF_TYPE_AMBER: 
      return read_amber_cdf_timestep(mydata, natoms, ts); 

    case CDF_TYPE_MMTK:
      return read_mmtk_cdf_timestep(mydata, natoms, ts); 
  }

  return MOLFILE_ERROR;
}


static molfile_plugin_t plugin;

VMDPLUGIN_API int VMDPLUGIN_init(void) {
  memset(&plugin, 0, sizeof(molfile_plugin_t));
  plugin.abiversion = vmdplugin_ABIVERSION;
  plugin.type = MOLFILE_PLUGIN_TYPE;
  plugin.name = "netcdf";
  plugin.prettyname = "NetCDF (AMBER, MMTK)";
  plugin.author = "Konrad Hinsen, John Stone";
  plugin.majorv = 1;
  plugin.minorv = 1;
  plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
  plugin.filename_extension = "nc,ncrst";
  plugin.open_file_read = open_cdf_read;
  plugin.read_structure = read_cdf_structure;
  plugin.read_next_timestep = read_cdf_timestep;
  plugin.close_file_read = close_cdf_read;
  return VMDPLUGIN_SUCCESS; 
}

VMDPLUGIN_API int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_fini(void) { return VMDPLUGIN_SUCCESS; }

