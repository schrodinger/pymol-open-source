/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_moldenplugin
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
 *      $RCSfile: moldenplugin.c,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.12 $       $Date: 2006/02/23 19:36:45 $
 *
 ***************************************************************************/

/* This is a plugin that will read input from a MOLDEN
** generated output file 
** some more details will go here soon 
** NOTE: The current version of the plugin relies
** on the fact that the [Atom] field comes before
** the [Geometries] field */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "molfile_plugin.h"



typedef struct {
  FILE *file;
  int trajectory;
  int numatoms;
  char *file_name;
  molfile_atom_t *atomlist;
} moldendata;

/* this will skip one line at a time */

static void eatline(FILE * fd)
{
  char readbuf[1025];
  fgets(readbuf, 1024, fd);
}

 
static void *open_molden_read(const char *filename, const char *filetype, int *natoms) {

  FILE *fd;
  moldendata *data;
  char moldentest1[7];
  char moldentest2[7];
  char buffer[1024];
  char tester[20];
  int i;
  char *s;

  fd = fopen(filename, "rb");
  if (!fd) return NULL;
  
  data = (moldendata *)malloc(sizeof(moldendata));
  data->file = fd;
  data->file_name = strdup(filename);

/* check if the file is MOLDEN format */

  fscanf(data->file, "%s %s",moldentest1,moldentest2);
  if (!strcmp(moldentest1,"[Molden") && \
      !strcmp(moldentest2,"Format]"))
  {
    printf("Detected MOLDEN file format!\n");
  }
  else
  {
    printf("The file does not seem to be in MOLDEN format!\n");
    return NULL;
  }

/* Unfortunately the molden file format has two possibilities:
** either there is a [ATOMS] section which provides the atom
** name, charges and coordinates or there is only a 
** XYZ style [GEOMETRIES] section, which still provides atom
** type and geometry information, hence I have to check which
** case I have */
  
/* check if there is an [ATOMS] section */

  do
  {
    i=fscanf(data->file, "%s",tester); 
    if (!strcmp(tester,"[Atoms]"))
    {
/* start counting the atoms; 
** read until I hit the first line that starts with a "["
** bracket */

      eatline(fd); 
      (*natoms)=0; 
      s=fgets(buffer,1024,fd);

/* Here I assume that the [Atoms] section goes
** on until another section starts, i.e. ther is
** a "[" or I encounter EOF */

      while ((buffer[0]!='[') && (s != NULL))
      {
       	(*natoms)++;     
	s=fgets(buffer,1024,fd);
      }
      data->numatoms=*natoms;
      rewind(fd);
      data->trajectory = 0;
      return data;
    }
    else if (!strcmp(tester,"[GEOMETRIES]"))
    {
      printf("Found [Geometry] section ...\n");
      data->trajectory = 1;	    
       
/* In this case I am lucky because the first line
** of the XYZ type [GEOMETRIES] input contains the
** number of atoms, i.e. skip to this line and the
** read this entry */

      eatline(fd);

      i=fscanf(data->file, "%d",natoms);
      if (i!=1)
      {
	printf("The [GEOMTRIES] output does not have \n");
	printf("the number of atoms in line number one !! \n");
      }

      data->numatoms=*natoms;
     
/* skip the next two lines, so the file can be parsed in the 
** structure section */      

      eatline(fd);
      eatline(fd);

      return data; 
    }
   } while (i>0);
  
  return NULL;
}

static int read_molden_structure(void *mydata, int *optflags, 
    molfile_atom_t *atoms) 
 {
  int i;
  char atname[1024];
  char buffer[1024];
  char geotest[11];
  int num,charge;
  float x,y,z;
  molfile_atom_t *atom;
  moldendata *data = (moldendata *)mydata;

  *optflags = MOLFILE_NOOPTIONS; /* no optional data */

/* here I have two possibilities, either there is an
** [Atoms] section (i.e. data->trajectory=0) and I can
** read there structure information right there,
** or I have to extract it from the [GEOMETRIES]
** output (i.e. data->trajectory=1) */

  if(data->trajectory==0)
  { 

/* Skip the first three lines */

    eatline(data->file);
    eatline(data->file);

/* Now read in the atom types, names, charges as well
** as x,y,z coordinates */

    for(i=0;i<data->numatoms;i++) 
    {
      atom = atoms+i;
      fgets(buffer,1024,data->file);    
      sscanf(buffer,"%s %d %d %f %f %f",atname,&num,&charge,\
	  &x,&y,&z);
      strncpy(atom->name,atname,sizeof(atom->name)); 
      strncpy(atom->type, atom->name, sizeof(atom->type));
      atom->resname[0] = '\0';
      atom->resid = 1;
      atom->chain[0] = '\0';
      atom->segid[0] = '\0';
    }


/* finally and important skip the the beginning of the xyz
** section */

    do 
    {
      fscanf(data->file, "%s",geotest);   
    } while (strcmp(geotest,"[GEOMETRIES]")!=0);
  
     printf("Found Geometry Section\n");

  /* advance to the beginning of the XYZ list */

    eatline(data->file);
    eatline(data->file);
    eatline(data->file); 

/* time to go back */

    return MOLFILE_SUCCESS; 
  }
  else if(data->trajectory==1)
  {
    
/* in the read section I already forwarded to the correct
** location in the file hence I can start reading right
** away */


    for(i=0;i<data->numatoms;i++) 
    {
      atom = atoms+i;
      fgets(buffer,1024,data->file);    
      sscanf(buffer,"%s %f %f %f",atname,&x,&y,&z);
      strncpy(atom->name,atname,sizeof(atom->name)); 
      strncpy(atom->type, atom->name, sizeof(atom->type));
      atom->resname[0] = '\0';
      atom->resid = 1;
      atom->chain[0] = '\0';
      atom->segid[0] = '\0';
    }

/* now rewind the file and go back to the [GEOMETRIES]
** section */

    rewind(data->file);
    
    do 
    {
      fscanf(data->file, "%s",geotest);   
    } while (strcmp(geotest,"[GEOMETRIES]")!=0);
  
     printf("Found Geometry Section\n");

/* advance to the beginning of the XYZ list */

     eatline(data->file);
     eatline(data->file);
     eatline(data->file);

/* time to go back */

     return MOLFILE_SUCCESS;

   }

  printf("Sorry, could not obtain structure information \n");
  printf("from either the [Atoms] or [GEOMETRIES] section! \n");
  printf("Please check your MOLDEN output file! \n"); 
  return MOLFILE_ERROR; 
}


static int read_next_timestep(void *mydata, int natoms, molfile_timestep_t *ts) {

  int i;
  char *k;
  char atname[1024];
  char buffer[1024];
  float x, y, z;
  
  moldendata *data = (moldendata *)mydata;
  
/* read the coordinates */

  for(i=0;i<data->numatoms;i++) 
  {
      k=fgets(buffer,1024,data->file);
      sscanf(buffer,"%s %f %f %f",atname,&x,&y,&z);

/* save coordinates only if given a timestep pointer
** otherwise assume that VMD would like to skip past
** it */

/* if I don't check for EOF file I get stuck in
** an infinite loop */

      if ( k==NULL)
      {
	return MOLFILE_ERROR;
      }
      
      if (ts!=NULL) 
      {
       ts->coords[3*i  ] = x;
       ts->coords[3*i+1] = y;
       ts->coords[3*i+2] = z;     
      } 
  }     

/* skip the two comment lines, which separate
** the individual coordinate entries 
** (for now, cause contains useful information
** like energies) */

  eatline(data->file);
  eatline(data->file); 

/* done and go back */

  return MOLFILE_SUCCESS;
}
  

static void close_molden_read(void *mydata) {
  moldendata *data = (moldendata *)mydata;
  fclose(data->file);
  free(data->file_name);
  free(data);
}


/* registration stuff */
static molfile_plugin_t moldenplugin = {
  vmdplugin_ABIVERSION,
  MOLFILE_PLUGIN_TYPE,                         /* type */
  "molden",                                    /* short name */
  "Molden",                                    /* pretty name */
  "Markus Dittrich",                           /* author */
  0,                                           /* major version */
  1,                                           /* minor version */
  VMDPLUGIN_THREADSAFE,                        /* is reentrant */
  "molden",
  open_molden_read,
  read_molden_structure,
  0,
  read_next_timestep, 
  close_molden_read,
  0,
  0,
  0, 
  0,
};

VMDPLUGIN_API int VMDPLUGIN_init() {
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)&moldenplugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_fini() {
  return VMDPLUGIN_SUCCESS;
}

