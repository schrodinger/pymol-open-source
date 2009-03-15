/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_xbgfplugin
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
 *      $RCSfile: xbgfplugin.C,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.32 $       $Date: 2009/02/20 22:28:42 $
 *
 ***************************************************************************/

#include "molfile_plugin.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#if defined(_AIX)
#include <strings.h>
#endif

#define LINESIZE 256
#define MAXBONDS 16

typedef struct {
  FILE *file;
  molfile_atom_t *atomlist;
  molfile_metadata_t *meta;
  int natoms, nbonds, optflags, coords_read;
  int *from, *to;
  float *bondorder;
} xbgfdata;


// Open the file and create the xbgf struct used to pass data to the other
// functions.
static void *open_xbgf_read(const char *path, const char *filetype, 
    int *natoms) {
  FILE *fd;
  xbgfdata *bgf;
  char line[LINESIZE]; 
  int nbonds, optflags;
  int numat=0;
  nbonds=0;
  int nbline; //Number of bonds in current line

  // Allocate and initialize the bgf structure
  bgf = (xbgfdata *) malloc(sizeof(xbgfdata));
  memset(bgf, 0, sizeof(xbgfdata));

  bgf->meta = (molfile_metadata_t *) malloc(sizeof(molfile_metadata_t));
  memset(bgf->meta, 0, sizeof(molfile_metadata_t));

  bgf->meta->remarklen = 0;
  bgf->meta->remarks = NULL;

  if ((fd = fopen(path, "r")) == NULL)
    return NULL;

  do {
    fgets(line, LINESIZE, fd);
    if ( ferror(fd) || feof(fd) ) {
      printf("xbgfplugin) Improperly terminated bgf file\n");
      return NULL;
    }

    if ((strncmp(line, "ATOM", 4) == 0) || (strncmp(line, "HETATM", 6)==0)) 
      numat++;

    if (strncmp(line,"CONECT",6)==0) {
      nbline=(strlen(line)-1)/6; 
      nbline -= 2;
      nbonds += nbline;
    }

    // Read the remarks
    if (strncmp(line, "REMARK", 4)==0 || strncmp(line, "LEWIS", 4)==0 || 
        strncmp(line, "VDW", 3)==0) {
      int len=strlen(line);
      int newlen = len + bgf->meta->remarklen;
      char *newstr=(char*) realloc(bgf->meta->remarks, newlen + 1);
      if (newstr != NULL) {
        bgf->meta->remarks = newstr;
        bgf->meta->remarks[bgf->meta->remarklen] = '\0';
        memcpy(bgf->meta->remarks + bgf->meta->remarklen, line, len);
        bgf->meta->remarks[newlen] = '\0';
        bgf->meta->remarklen = newlen;
      }
    }

    optflags = MOLFILE_INSERTION | MOLFILE_CHARGE | MOLFILE_BFACTOR | MOLFILE_OCCUPANCY | MOLFILE_ATOMICNUMBER;
  } while ( strncmp(line, "END", 3) );
    
  *natoms = numat;
  rewind(fd);

  bgf->file = fd;
  bgf->natoms = *natoms;
  bgf->nbonds = nbonds;

  bgf->optflags = optflags;
  bgf->coords_read = 0;
  bgf->from = NULL;
  bgf->to = NULL;
  bgf->bondorder = NULL;
                                                                                
  return bgf;
}


static void adjust_xbgf_field_string(char *field) {
  int i, len;

  len = strlen(field);
  while (len > 0 && field[len-1] == ' ') {
    field[len-1] = '\0';
    len--;
  }

  while (len > 0 && field[0] == ' ') {
    for (i=0; i < len; i++)
      field[i] = field[i+1];
    len--;
  }
}

static void get_xbgf_coordinates(const char *record, 
                                float *x, float *y, float *z) {
  char numstr[50]; /* store all fields in one array to save memset calls */
  memset(numstr, 0, sizeof(numstr));
  if (x != NULL) {
    strncpy(numstr, record + 32, 10);
    *x = (float) atof(numstr);
  }

  if (y != NULL) {
    strncpy(numstr+10, record + 42, 10);
    *y = (float) atof(numstr+10);
  }

  if (z != NULL) {
    strncpy(numstr+20, record + 52, 10);
    *z = (float) atof(numstr+20);
  }
}


static void get_xbgf_fields(const char *record, char *name, char *resname, 
                           char *chain, char* segname, float *occupancy, 
                           float *bfactor, int *elementindex,
                           int *resid, char *type, float *charge,
                           float *x, float *y, float *z, char* insert) {
  char tempresid[6];
  char tempcharge[8];
  char tempbeta[7];
  char tempocc[7];
  char tempelem[4];
  strcpy(insert, " ");

  /* get atom name */
  strncpy(name, record + 14, 5);
  name[5] = '\0';
  adjust_xbgf_field_string(name); /* remove spaces from the name */

  /* get residue name */
  strncpy(resname, record + 20, 4);
  resname[4] = '\0';
  adjust_xbgf_field_string(resname); /* remove spaces from the resname */

  /* set segname */
  strncpy(segname, record + 101, 4);
  segname[4]='\0';
  adjust_xbgf_field_string(segname); /* remove spaces from the segname */

  /* get chain name */
  chain[0] = record[25];
  chain[1] = '\0';

  /* get residue id number */
  strncpy(tempresid, record + 27, 5);
  tempresid[5] = '\0';
  adjust_xbgf_field_string(tempresid); /* remove spaces from the resid */
  *resid=atoi(tempresid);

  /* get force field type */
  strncpy(type, record+63, 5);
  type[5]='\0';
  adjust_xbgf_field_string(type); /* remove spaces */

  /* get charge*/
  strncpy(tempcharge, record + 74, 7);
  tempcharge[7] = '\0';
  adjust_xbgf_field_string(tempcharge); /* remove spaces from the charge */
  *charge=atof(tempcharge);

  /* Get B factor, occupancy, and element */
  strncpy(tempbeta, record + 83, 6);
  tempbeta[6] = '\0';
  adjust_xbgf_field_string(tempbeta); /* remove spaces from the beta field */
  *bfactor=atof(tempbeta);
  
  strncpy(tempocc, record + 90, 6);
  tempocc[6] = '\0';
  adjust_xbgf_field_string(tempocc); /* remove spaces from the occupancy field */
  *occupancy=atof(tempocc);
  
  strncpy(tempelem, record + 97, 3);
  tempelem[3] = '\0';
  adjust_xbgf_field_string(tempelem); /* remove spaces from the element */
  *elementindex=atoi(tempelem);

  /* get x, y, and z coordinates */
  get_xbgf_coordinates(record, x, y, z);
}  


// Read atom information, but not coordinates.
static int read_xbgf_structure(void *v, int *optflags, molfile_atom_t *atoms) {
  xbgfdata *bgf = (xbgfdata *)v;
  char line[LINESIZE]; 
  molfile_atom_t *atom;
  int natoms=0;

  //optflags = MOLFILE_INSERTION | MOLFILE_CHARGE | MOLFILE_BFACTOR | MOLFILE_OCCUPANCY | MOLFILE_ATOMICNUMBER;
  *optflags = bgf->optflags;

  // Find and read the ATOM record
  rewind(bgf->file);
  do {
    fgets(line, LINESIZE, bgf->file);
    if ( ferror(bgf->file) || feof(bgf->file) ) {
      printf("xbgfplugin) FORMAT ATOM record found in file.\n");
      return MOLFILE_ERROR;
    }
  } while ( strncmp(line, "FORMAT ATOM", 11) );

  // Read the atoms
  do {
    fgets(line, LINESIZE, bgf->file);
    if ( ferror(bgf->file) || feof(bgf->file) ) {
      printf("xbgfplugin) Error occurred reading atom record.\n");
      return MOLFILE_ERROR;
    }

    if (strncmp(line, "ATOM", 4) && strncmp(line, "HETATM", 6)) continue;
    atom=atoms+natoms;
    natoms++;
    get_xbgf_fields(line, atom->name, atom->resname, atom->chain, atom->segid, 
                   &atom->occupancy, &atom->bfactor, &atom->atomicnumber, 
                   &atom->resid, atom->type, &atom->charge, 
                   NULL, NULL, NULL, atom->insertion);
  } while (strncmp(line, "END", 3));

  bgf->natoms = natoms;

  return MOLFILE_SUCCESS;
}


// Read atom coordinates
static int read_xbgf_timestep(void *v, int natoms, molfile_timestep_t *ts) {
  xbgfdata *bgf = (xbgfdata *)v;
  char line[LINESIZE];
  int i;
  float x, y, z;

  // Since the file is rewound when coordinates are read, EOF shouldn't
  // happen. Instead, use a flag to indicate that the single timestep has
  // been read
  if (bgf->coords_read) {
    return MOLFILE_EOF;
  }

  // Find and read the ATOM record
  rewind(bgf->file);
  do {
    fgets(line, LINESIZE, bgf->file);
    if ( ferror(bgf->file) || feof(bgf->file) ) {
      printf("xbgfplugin) No FORMAT ATOM record found in file.\n");
      return MOLFILE_ERROR;
    }
  } while ( strncmp(line, "FORMAT ATOM", 11) );

  // Read the atoms
  for (i = 0; i < bgf->natoms; i++) {
    fgets(line, LINESIZE, bgf->file);
    if ( ferror(bgf->file) || feof(bgf->file) ) {
      printf("xbgfplugin) Error occurred reading atom coordinates.\n");
      return MOLFILE_ERROR;
    }

    // skip comments and blank lines
    if (strncmp(line,"ATOM",4)!=0 && strncmp(line,"HETATM",6)!=0) continue;

    get_xbgf_coordinates(line, &x, &y, &z);

    if (ts) {
      ts->coords[3*i  ] = x;
      ts->coords[3*i+1] = y;
      ts->coords[3*i+2] = z;
    }
  }

  bgf->coords_read = 1;
  return MOLFILE_SUCCESS;
}


static void *open_xbgf_write(const char *filename, const char *filetype, 
                           int natoms) {
  FILE *fd;
  xbgfdata *data;

  if ((fd = fopen(filename, "w")) == NULL) { 
    printf("xbgfplugin) Error, unable to open xbgf file %s for writing\n",
            filename);
    return NULL;
  }
  
  data = (xbgfdata *) malloc(sizeof(xbgfdata));
  memset(data, 0, sizeof(xbgfdata));
  data->natoms = natoms;
  data->file = fd;
  data->nbonds = 0;
  return data;
}


static int write_xbgf_structure(void *mydata, int optflags, 
                               const molfile_atom_t *atoms) {
  fflush(stdout);
  xbgfdata *data = (xbgfdata *)mydata;
  data->atomlist = (molfile_atom_t *)malloc(data->natoms*sizeof(molfile_atom_t));
  memcpy(data->atomlist, atoms, data->natoms*sizeof(molfile_atom_t));
  return MOLFILE_SUCCESS;
}


static int read_xbgf_bonds_aux(void *v, int *nbonds, int **fromptr, int **toptr, float **bondorderptr) {
  xbgfdata *bgf = (xbgfdata *)v;
  char line[LINESIZE]; 
  char nextline[LINESIZE]; 
  if (bgf->nbonds == 0) {
    *nbonds = 0;
    *fromptr = NULL;
    *toptr = NULL;
    *bondorderptr = NULL;
    return MOLFILE_SUCCESS;
  }

  // Find and read the BOND record
  rewind(bgf->file);
  do {
    fgets(line, LINESIZE, bgf->file);
    if ( ferror(bgf->file) || feof(bgf->file) ) {
      printf("xbgfplugin) No bond record found in file.\n");
      return MOLFILE_ERROR;
    }
  } while ( strncmp(line, "FORMAT CONECT", 13) != 0 );

  // Read the bonds
  int j; //From atom
  int k; //To atom
  bool conline=false; //true if line after the conect line is an order line
  char currbond[7]="xxxxxx"; //Stores current bond field
  char currcon[7]="xxxxxx"; //Stores current ORDER field
  char* bondptr; //pointer to current position in bond line
  char* conptr; //pointer to current position in order line
  int bonds[MAXBONDS]; //Stores bonds of current atom
  float orders[MAXBONDS]; //Stores bond orders of current atom
  int numbonds; //Stores number of bonds of current atom
  int numords; //Stores number of bond order records of current atom
  float bo; //current bond order
  int i=0; //Number of the current bond
  int numfields=0; //number of fields in the current line
  fgets(line, LINESIZE, bgf->file);
  while (1) {
    // bondptr=NULL;
    //conptr=NULL;
    conline=false;
    if (strncmp(line,"END",3)==0)
      break;

    fgets(nextline, LINESIZE, bgf->file);
    if ( ferror(bgf->file) || feof(bgf->file) ) {
      printf("xbgfplugin) Error occurred reading bond record.\n");
      return MOLFILE_ERROR;
    }

    if (strncmp(nextline,"ORDER",5)==0) 
      conline=true;

    if (strncmp(line,"CONECT",6)==0) {
      numfields=(strlen(line)-1)/6;
      bondptr=&line[0];
      numfields--;
      bondptr += 6;
      numbonds=0;
      numords=0;
      strncpy(currbond,bondptr,6);
      j=atoi(currbond);
      numfields--;
      bondptr += 6;

      while ((numfields > 0) && (numbonds <= MAXBONDS)) {
        strncpy(currbond,bondptr,6);
        numfields--;
        bondptr += 6;
        bonds[numbonds]=atoi(currbond);
        numbonds++;
      }

      if (conline) {
        numfields=(strlen(line)-1)/6;
        conptr=&nextline[0];
        numfields -= 2;
        conptr += 12;
        numords=0;
        while ((numfields > 0) && (numords < numbonds)) {
          strncpy(currcon,conptr,6);
          numfields--;
          conptr+=6;
          bo=atof(currcon);
          orders[numords]=bo;
          numords++;
        }
      }

      for (int l=0;l<numbonds;l++) {
        k=bonds[l];
        if (j<k) {
          bgf->from[i]=j;
          bgf->to[i]=k;
          if (conline) {
            bgf->bondorder[i]=orders[l];
          } else {
            bgf->bondorder[i]=1.0;
          }
          i++;
        }
      }
        
      if (conline) {
        fgets(line, LINESIZE, bgf->file);
      } else {
        strncpy(line,nextline,LINESIZE);
      }
    } else {
      strncpy(line,nextline,LINESIZE);
    }
  }

  *nbonds = i;
  *fromptr = bgf->from;
  *toptr = bgf->to;
  *bondorderptr = bgf->bondorder; // not implemented yet

  return MOLFILE_SUCCESS;
}


static int read_xbgf_bonds(void *v, int *nbonds, int **fromptr, int **toptr, 
                           float **bondorderptr, int **bondtype,
                           int *nbondtypes, char ***bondtypename) {
  xbgfdata *bgf = (xbgfdata *)v;

  /* now read bond data */
  *nbonds=bgf->nbonds;
  if (bgf->nbonds > 0) {
    bgf->from = (int *) malloc(*nbonds*sizeof(int));
    bgf->to = (int *) malloc(*nbonds*sizeof(int));
    bgf->bondorder = (float *) malloc(*nbonds*sizeof(float));

    if ((read_xbgf_bonds_aux(bgf, nbonds, &(bgf->from), &(bgf->to), &(bgf->bondorder))) != MOLFILE_SUCCESS) {
      fclose(bgf->file);
      bgf->file = NULL;
      return MOLFILE_ERROR;
    }
    *fromptr = bgf->from;
    *toptr = bgf->to;
    *bondorderptr = bgf->bondorder; // not implemented yet
    *bondtype = NULL;
    *nbondtypes = 0;
    *bondtypename = NULL;
  } else {
    printf("xbgfplugin) WARNING: no bonds defined in xbgf file.\n");
    *fromptr = NULL;
    *toptr = NULL;
    *bondorderptr = NULL;
    *bondtype = NULL;
    *nbondtypes = 0;
    *bondtypename = NULL;
  }
  return MOLFILE_SUCCESS;
}


static int write_xbgf_timestep(void *mydata, const molfile_timestep_t *ts) {
  fflush(stdout);
  xbgfdata *data = (xbgfdata *)mydata; 
  const molfile_atom_t *atom;
  const float *pos;
  int i;

  //print header block
  fprintf(data->file, "BIOGRF  332\n");
  fprintf(data->file, "REMARK NATOM %4i\n", data->natoms);
  fprintf(data->file, "FORCEFIELD DREIDING\n");
  fprintf(data->file, "FORMAT ATOM   (a6,1x,i6,1x,a5,1x,a4,1x,a1,1x,i5,3f10.5,1x,a5,i3,i2,1x,f8.5,1x,f6.3,1x,f6.3,1x,i3,1x,a4)\n");

  //print atoms block
  atom = data->atomlist;
  pos = ts->coords;
  int numbonds=0;
  int lp=0;
  for (i = 0; i < data->natoms; i++) {
    fprintf(data->file, "%-6s %6i %5s %4s %1s %5i%10.5f%10.5f%10.5f %-5s%3i%2i %8.5f %6.3f %6.3f %3i %4s\n", "ATOM", i+1, atom->name, atom->resname, atom->chain, atom->resid, pos[0], pos[1], pos[2], atom->type, numbonds, lp, atom->charge, atom->bfactor, atom->occupancy, atom->atomicnumber, atom->segid);
    ++atom; 
    pos += 3;
  }

  //write the connectivity data
  fprintf(data->file,"FORMAT CONECT (a6,14i6) \nFORMAT ORDER (a6,i6,13f6.3)\n");
    
  //iterate through the bond arrays and write them all
  int* bonds=(int *)malloc((data->natoms+1) * sizeof(int) * MAXBONDS);
  float* orders=(float *)malloc((data->natoms+1)*sizeof(float) * MAXBONDS);
  int* numcons=(int *)malloc((data->natoms+1)*sizeof(int));
  for (i=0;i<data->natoms+1;i++) {
    numcons[i]=0;
  }

  int j,k; //indices for atoms being bonded
  float o; //bond order
  for (i=0;i<data->nbonds;i++) {
    j=data->from[i];
    k=data->to[i];

    if (data->bondorder != NULL)
      o=data->bondorder[i];
    else 
      o=1.0f;

    numcons[j]++;
    numcons[k]++;
    if (numcons[j]>MAXBONDS) {
      printf("xbgfplugin) Warning: Bond overflow. Not all bonds were written\n");
      numcons[j]--;
      numcons[k]--;
      continue;
    }
       
    if (numcons[k]>MAXBONDS) {
      printf("xbgfplugin) Warning: Bond overflow. Not all bonds were written\n");
      numcons[k]--;
      numcons[j]--;
      continue;
    }
    bonds[6*j+numcons[j]-1]=k;
    bonds[6*k+numcons[k]-1]=j;
    orders[6*j+numcons[j]-1]=o;
    orders[6*k+numcons[k]-1]=o;
  }

  for (i=1;i<=data->natoms;i++) {
    fprintf(data->file,"CONECT%6i",i);
    for (j=0;j<numcons[i];j++) {
      fprintf(data->file,"%6i",bonds[6*i+j]);
    }
    fprintf(data->file,"\nORDER %6i",i);
    for (j=0;j<numcons[i];j++) {
      fprintf(data->file,"%6.3f",orders[6*i+j]);
    }
    fprintf(data->file,"\n");
  }

  if (bonds != NULL) {
    free(bonds);
    bonds = NULL;
  }
  if (orders != NULL) {
    free(orders);
    orders = NULL;
  }
  if (numcons != NULL) {
    free(numcons);
    numcons = NULL;
  }

  fprintf(data->file,"END\n");
  return MOLFILE_SUCCESS;
}

static int write_xbgf_bonds(void *v, int nbonds, int *fromptr, int *toptr, 
                            float *bondorderptr,  int *bondtype, 
                            int nbondtypes, char **bondtypename) {
  xbgfdata *data = (xbgfdata *)v;
  data->from = (int*) malloc (nbonds * sizeof(int));
  data->to = (int*) malloc (nbonds * sizeof(int));
  data->nbonds = nbonds;
  fflush(stdout);


  //set the pointers for use later
  for (int i=0;i<nbonds;i++) {
    data->from[i]=fromptr[i];
    data->to[i]=toptr[i];
  }

  if (bondorderptr != NULL) {
    data->bondorder = (float*) malloc (nbonds * sizeof(float));
    for (int i=0;i<nbonds;i++) {
      data->bondorder[i]=bondorderptr[i];
    }
  }


  return MOLFILE_SUCCESS;
}

static void close_xbgf_write(void *mydata) {
  xbgfdata *data = (xbgfdata *)mydata;
  if (data) {
    fclose(data->file);

    if (data->atomlist != NULL) free(data->atomlist);
    data->atomlist = NULL;
    if (data->from != NULL) free(data->from);
    data->from = NULL;
    if (data->to != NULL) free(data->to);
    data->to = NULL;
    if (data->bondorder != NULL) free(data->bondorder);
    data->bondorder = NULL;
    free(data);
  }
}

//
// Free the memory used by the bgf structure
static void close_xbgf_read(void *v) {
  xbgfdata *bgf = (xbgfdata *)v;
  if (bgf) {
    if (bgf->file) fclose(bgf->file);
    if (bgf->from != NULL) free(bgf->from);
    if (bgf->to != NULL)   free(bgf->to);
    if (bgf->bondorder != NULL) free(bgf->bondorder);

    if (bgf->meta->remarks != NULL)
      free(bgf->meta->remarks);
    if (bgf->meta != NULL) 
      free(bgf->meta);
    free(bgf);
  }
  bgf=NULL;
}


static int read_xbgf_molecule_metadata(void *v, molfile_metadata_t **metadata) {
  xbgfdata *bgf = (xbgfdata *)v; 
  *metadata = bgf->meta;
  return MOLFILE_SUCCESS;
}


static molfile_plugin_t plugin;

VMDPLUGIN_EXTERN int VMDPLUGIN_init() {
  memset(&plugin, 0, sizeof(molfile_plugin_t));
  plugin.abiversion = vmdplugin_ABIVERSION;
  plugin.type = MOLFILE_PLUGIN_TYPE;
  plugin.name = "xbgf";
  plugin.prettyname = "Internal Paratool Format";
  plugin.author = "Peter Freddolino ";
  plugin.majorv = 0;
  plugin.minorv = 13;
  plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
  plugin.filename_extension = "xbgf";
  plugin.open_file_read = open_xbgf_read;
  plugin.read_structure = read_xbgf_structure;
  plugin.read_bonds = read_xbgf_bonds;
  plugin.read_next_timestep = read_xbgf_timestep;
  plugin.close_file_read = close_xbgf_read;
  plugin.open_file_write = open_xbgf_write;
  plugin.write_structure = write_xbgf_structure;
  plugin.write_timestep = write_xbgf_timestep;
  plugin.close_file_write = close_xbgf_write;
  plugin.read_molecule_metadata = read_xbgf_molecule_metadata;
  plugin.write_bonds = write_xbgf_bonds;
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_EXTERN int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)(void *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_EXTERN int VMDPLUGIN_fini() {
  return VMDPLUGIN_SUCCESS;
}

