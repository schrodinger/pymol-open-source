/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_bgfplugin
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
 *      $RCSfile: bgfplugin.C,v $
 *      $Author: petefred $       $Locker:  $             $State: Exp $
 *      $Revision: 1.18 $       $Date: 2006/03/22 18:11:56 $
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

typedef struct {
  FILE *file;
  molfile_atom_t *atomlist;
  int natoms, nbonds, optflags, coords_read;
  int *from, *to;
  float *bondorder;
} bgfdata;

// Open the file and create the bgf struct used to pass data to the other
// functions.
static void *open_bgf_read(const char *path, const char *filetype, 
    int *natoms) {
  FILE *fd;
  bgfdata *bgf;
  char line[LINESIZE]; 
  int nbonds, optflags;
  int numat=0;
  nbonds=0;
  int nbline; //Number of bonds in current line

  if ((fd = fopen(path, "r")) == NULL)
    return NULL;

  do {
    fgets(line, LINESIZE, fd);
    if ( ferror(fd) || feof(fd) ) {
      printf("bgfplugin) Improperly terminated bgf file\n");
      return NULL;
    }

    if ((strncmp(line, "ATOM", 4) == 0) || (strncmp(line, "HETATM", 6)==0)) 
      numat++;

    if (strncmp(line,"CONECT",6)==0) {
      nbline=(strlen(line)-1)/6; 
      nbline -= 2;
      nbonds += nbline;
    }

  } while ( strncmp(line, "END", 3) );
    
  optflags = MOLFILE_INSERTION | MOLFILE_CHARGE; 
  *natoms = numat;
  rewind(fd);

  // Allocate and initialize the bgf structure
  bgf = new bgfdata;
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


static void adjust_bgf_field_string(char *field) {
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


static void get_bgf_coordinates(const char *record, 
                                float *x, float *y, float *z) {
  char numstr[50]; /* store all fields in one array to save memset calls */
  memset(numstr, 0, sizeof(numstr));
  if (x != NULL) {
    strncpy(numstr, record + 31, 10);
    *x = (float) atof(numstr);
  }

  if (y != NULL) {
    strncpy(numstr+10, record + 41, 10);
    *y = (float) atof(numstr+10);
  }

  if (z != NULL) {
    strncpy(numstr+20, record + 51, 10);
    *z = (float) atof(numstr+20);
  }
}


static void get_bgf_fields(const char *record, char *name, char *resname, 
                           char *chain, char* segname,
                           int *resid, char *type, float *charge,
                           float *x, float *y, float *z) {
  char tempresid[6];
  char tempcharge[9];

  /* get atom name */
  strncpy(name, record + 13, 5);
  name[5] = '\0';
  adjust_bgf_field_string(name); /* remove spaces from the name */

  /* get residue name */
  strncpy(resname, record + 19, 4);
  resname[4] = '\0';
  adjust_bgf_field_string(resname); /* remove spaces from the resname */

  /* set segname */
  segname[0]='\0';

  /* get chain name */
  chain[0] = record[23];
  chain[1] = '\0';

  /* get residue id number */
  strncpy(tempresid, record + 26, 5);
  tempresid[5] = '\0';
  adjust_bgf_field_string(tempresid); /* remove spaces from the resid */
  *resid=atoi(tempresid);

  /* get force field type */
  strncpy(type, record+61, 5);
  type[5]='\0';
  adjust_bgf_field_string(type);

  /* get charge*/
  strncpy(tempcharge, record + 72, 8);
  tempcharge[8] = '\0';
  adjust_bgf_field_string(tempcharge); /* remove spaces from the charge */
  *charge=atof(tempcharge);

  /* get x, y, and z coordinates */
  get_bgf_coordinates(record, x, y, z);
}  


// Read atom information, but not coordinates.
static int read_bgf_structure(void *v, int *optflags, molfile_atom_t *atoms) {
  bgfdata *bgf = (bgfdata *)v;
  char line[LINESIZE]; 
  molfile_atom_t *atom;
  int natoms=0;

  *optflags = bgf->optflags;

  // Find and read the ATOM record
  rewind(bgf->file);
  do {
    fgets(line, LINESIZE, bgf->file);
    if ( ferror(bgf->file) || feof(bgf->file) ) {
      printf("bgfplugin) FORMAT ATOM record not found in file.\n");
      return MOLFILE_ERROR;
    }
  } while ( strncmp(line, "FORMAT ATOM", 11) );

  // Read the atoms
  do {
    fgets(line, LINESIZE, bgf->file);
    if ( ferror(bgf->file) || feof(bgf->file) ) {
      printf("bgfplugin) Error occurred reading atom record.\n");
      return MOLFILE_ERROR;
    }

    if (strncmp(line, "ATOM", 4) && strncmp(line, "HETATM", 6)) 
      continue;

    atom=atoms+natoms;
    natoms++;

    get_bgf_fields(line, atom->name, atom->resname, atom->chain, 
                   atom->segid, &atom->resid, atom->type, &atom->charge, 
                   NULL, NULL, NULL);
  } while (strncmp(line, "END", 3));

  bgf->natoms = natoms;

  return MOLFILE_SUCCESS;
}


// Read atom coordinates
static int read_bgf_timestep(void *v, int natoms, molfile_timestep_t *ts) {
  bgfdata *bgf = (bgfdata *)v;
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
      printf("bgfplugin) No FORMAT ATOM record found in file.\n");
      return MOLFILE_ERROR;
    }
  } while ( strncmp(line, "FORMAT ATOM", 11) );

  // Read the atoms
  for (i = 0; i < bgf->natoms; i++) {
    fgets(line, LINESIZE, bgf->file);
    if ( ferror(bgf->file) || feof(bgf->file) ) {
      printf("bgfplugin) Error occurred reading atom coordinates.\n");
      return MOLFILE_ERROR;
    }

    // skip comments and blank lines
    if (strncmp(line,"ATOM",4)!=0 && strncmp(line,"HETATM",6)!=0) continue;

    get_bgf_coordinates(line, &x, &y, &z);

    if (ts) {
      ts->coords[3*i  ] = x;
      ts->coords[3*i+1] = y;
      ts->coords[3*i+2] = z;
    }
  }

  bgf->coords_read = 1;
  return MOLFILE_SUCCESS;
}


static void *open_bgf_write(const char *filename, const char *filetype, 
                           int natoms) {
  FILE *fd;
  bgfdata *data;

  if ((fd = fopen(filename, "w")) == NULL) {
    printf("Error) Unable to open bgf file %s for writing\n", filename);
    return NULL;
  }
  
  data = (bgfdata *)malloc(sizeof(bgfdata));
  data->natoms = natoms;
  data->file = fd;
  return data;
}


static int write_bgf_structure(void *mydata, int optflags, 
                               const molfile_atom_t *atoms) {
  bgfdata *data = (bgfdata *)mydata;
  data->atomlist = (molfile_atom_t *)malloc(data->natoms*sizeof(molfile_atom_t));
  memcpy(data->atomlist, atoms, data->natoms*sizeof(molfile_atom_t));
  return MOLFILE_SUCCESS;
}

void getatomfield(char* atomfield, const char* resname) {
  if ((strncmp(resname,"ALA",3) == 0) || (strncmp(resname,"ASP",3) == 0) || (strncmp(resname,"ARG",3) == 0) || (strncmp(resname,"ASN",3) == 0) || (strncmp(resname,"CYS",3) == 0) || (strncmp(resname,"GLN",3) == 0) || (strncmp(resname,"GLU",3) == 0) || (strncmp(resname,"GLY",3) == 0) || (strncmp(resname,"HIS",3) == 0) || (strncmp(resname,"ILE",3) == 0) || (strncmp(resname,"LEU",3) == 0) || (strncmp(resname,"LYS",3) == 0) || (strncmp(resname,"MET",3) == 0) || (strncmp(resname,"PHE",3) == 0) || (strncmp(resname,"PRO",3) == 0) || (strncmp(resname,"SER",3) == 0) || (strncmp(resname,"THR",3) == 0) || (strncmp(resname,"TRP",3) == 0) || (strncmp(resname,"TYR",3) == 0) || (strncmp(resname,"VAL",3) == 0) || (strncmp(resname,"ADE",3) == 0) || (strncmp(resname,"THY",3) == 0) || (strncmp(resname,"GUA",3) == 0) || (strncmp(resname,"CYT",3) == 0) || (strncmp(resname,"URA",3) == 0) || (strncmp(resname,"HSD",3) == 0) || (strncmp(resname,"HSE",3) == 0) || (strncmp(resname,"HSP",3) == 0)) {
    strncpy(atomfield, "ATOM  \0", 7);
  } else {
    strncpy(atomfield, "HETATM\0", 7);
  }
}
      
static void getdreiidff(char* outputtype, const char* psftype, int& numbonds, int& lp) {
  //Note that while this function isn't used yet, it actually IS important, and will be enabled in the future pending dicussion with some bgf users
  if (strncmp(psftype,"H",1)==0) {
    //It's a hydrogen
    //FIXME: Doesn't properly identify acidic hydrogens yet
    strncpy(outputtype, "H_  ",4);
    numbonds=1;
    lp=0;
    return;
  } else if (strncmp(psftype,"C",1)==0) {
    //It's a carbon... probably
    if (strncmp(psftype,"C ",2)==0 || strncmp(psftype,"CA ",3)==0 || strncmp(psftype,"CPH",3)==0 || strncmp(psftype,"CPT",3)==0 || strncmp(psftype,"CC ",3)==0 || strncmp(psftype,"CD ",3)==0 || strncmp(psftype,"CN1",3)==0 || strncmp(psftype,"CN2",3)==0 || strncmp(psftype,"CN3",3)==0 || strncmp(psftype,"CN4",3)==0 || strncmp(psftype,"CN5",3)==0 || strncmp(psftype,"CNA",3)==0) {
      strncpy(outputtype, "C_2 ",4);
      numbonds=3;
      lp=0;
      return; 
    } else {
      strncpy(outputtype, "C_3 ",4);
      numbonds=4;
      lp=0;
      return; 
    }  
  } else if (strncmp(psftype,"N",1)==0) {
    //It"s probably nitrogen
    if (strncmp(psftype,"NR",2)==0 || strncmp(psftype,"NH1",3)==0 || strncmp(psftype,"NH2",3)==0 || strncmp(psftype,"NC2",3)==0 || strncmp(psftype,"NY",2)==0 || (strncmp(psftype,"NN",2)==0 && strncmp(psftype,"NN6",3)!=0)) {
      strncpy(outputtype, "N_R",4);
      numbonds=3;
      lp=0;
      return;
    } else {
      strncpy(outputtype, "N_3 ",4);
      numbonds=3;
      lp=1;
      return;
    }
  } else if (strncmp(psftype,"O",1)==0) {
    //Probably an oxygen
    if (strncmp(psftype,"OH1",3)==0 || strncmp(psftype,"OS",2)==0 || strncmp(psftype,"OT ",3)==0 || strncmp(psftype,"ON4",3)==0 || strncmp(psftype,"ON5",3)==0 || strncmp(psftype,"ON6",3)==0) {
      strncpy(outputtype, "O_3 ",4);
      numbonds=2;
      lp=2;
      return;
   } else {
      strncpy(outputtype, "O_2 ",4);
      numbonds=1;
      lp=2;
      return;
    }
  } else if (strncmp(psftype,"S",1)==0) {
    strncpy(outputtype, "S_3 ",4);
    numbonds=2;
    lp=2;
    return;
  } else if (strncmp(psftype,"P",1)==0) {
    strncpy(outputtype, "P_3 ",4);
    numbonds=6;
    lp=0;
    return;
  } else {
    strncpy(outputtype, "X_  ",4);
    numbonds=0;
    lp=0;
    return;
  }
}


static int read_bgf_bonds(void *v, int *nbonds, int **fromptr, int **toptr, float **bondorderptr) {
  bgfdata *bgf = (bgfdata *)v;
  char line[LINESIZE]; 
  char nextline[LINESIZE]; 
  if (bgf->nbonds == 0) {
    *nbonds = 0;
    *fromptr = NULL;
    *toptr = NULL;
    return MOLFILE_SUCCESS;
  }

  // Allocate memory for the from and to arrays. This will be freed in
  // close_mol2_read
  bgf->from = new int[bgf->nbonds];
  bgf->to = new int[bgf->nbonds];
  bgf->bondorder = new float[bgf->nbonds];

  // Find and read the BOND record
  rewind(bgf->file);
  do {
    if ( ferror(bgf->file) || feof(bgf->file) ) {
      printf("bgfplugin) No bond record found in file.\n");
      return MOLFILE_ERROR;
    }
    fgets(line, LINESIZE, bgf->file);
  } while ( strncmp(line, "FORMAT CONECT", 13) != 0 );

  // Read the bonds
  int j; //From atom
  int k; //To atom
  bool conline=false; //true if line after the conect line is an order line
  char currbond[7]="xxxxxx"; //Stores current bond field
  char currcon[7]="xxxxxx"; //Stores current ORDER field
  char* bondptr; //pointer to current position in bond line
  char* conptr; //pointer to current position in order line
  int bonds[8]; //Stores bonds of current atom
  float orders[8]; //Stores bond orders of current atom
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

    if (strncmp(line,"END", 3)==0) 
      break;

    fgets(nextline, LINESIZE, bgf->file);
    if ( ferror(bgf->file) || feof(bgf->file) ) {
      printf("bgfplugin) Error occurred reading bond record.\n");
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

      while ((numfields > 0) && (numbonds < 8)) {
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
  *bondorderptr = bgf->bondorder; 

  return MOLFILE_SUCCESS;
}


static int read_bonds(void *v, int *nbonds, int **fromptr, int **toptr, float **bondorderptr) {
  bgfdata *bgf = (bgfdata *)v;

  *nbonds=bgf->nbonds;
  if (bgf->nbonds > 0) {
    bgf->from = (int *) malloc(*nbonds*sizeof(int));
    bgf->to = (int *) malloc(*nbonds*sizeof(int));
    bgf->bondorder = (float *) malloc(*nbonds*sizeof(float));

    if ((read_bgf_bonds(bgf, nbonds, &(bgf->from), &(bgf->to), &(bgf->bondorder))) != MOLFILE_SUCCESS) {
      fclose(bgf->file);
      bgf->file = NULL;
      return MOLFILE_ERROR;
    }

    *fromptr = bgf->from;
    *toptr = bgf->to;
    *bondorderptr = bgf->bondorder; 
  } else {
    printf("bgfplugin) WARNING: no bonds defined in bgf file.\n");
    *fromptr = NULL;
    *toptr = NULL;
    *bondorderptr = NULL;
  }

  return MOLFILE_SUCCESS;
}


static int write_bgf_timestep(void *mydata, const molfile_timestep_t *ts) {
  bgfdata *data = (bgfdata *)mydata; 
  const molfile_atom_t *atom;
  const float *pos;
  int i;

  //print header block
  fprintf(data->file, "BIOGRF  332\n");
  fprintf(data->file, "REMARK NATOM %4i\n", data->natoms);
  fprintf(data->file, "FORCEFIELD DREIDING\n");
  fprintf(data->file, "FORMAT ATOM   (a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,a5,i3,i2,1x,f8.5,i2,i4,f10.5)\n");

  atom = data->atomlist;
  pos = ts->coords;
  int numbonds=0;
  int lp=0;
  char atomfield[7];
  for (i = 0; i < data->natoms; i++) {
    getatomfield(&atomfield[0], atom->resname);
    fprintf(data->file, "%-6s %5i %5s %3.3s %1s %5i%10.5f%10.5f%10.5f %-5s%3i%2i %8.5f%2i%4i\n", atomfield, i+1, atom->name, atom->resname, atom->chain, atom->resid, pos[0], pos[1], pos[2], atom->type, numbonds, lp, atom->charge, 0, 0);
    ++atom; 
    pos += 3;
  }

  //write the connectivity data
  fprintf(data->file,"FORMAT CONECT (a6,14i6) \nFORMAT ORDER (a6,i6,13f6.3)\n");

  //iterate through the bond arrays and write them all
  int *bonds    =(int *)  malloc((data->natoms+1) * sizeof(int) * 6);
  float *orders =(float *)malloc((data->natoms+1) * sizeof(float) * 6);
  int *numcons  =(int *)  malloc((data->natoms+1) * sizeof(int));

  for (i=0;i<data->natoms+1;i++) {
    numcons[i]=0;
  }

  int j, k;         //indices for atoms being bonded
  float o;          //bond order
  bool printorder;  //flag to print bond orders
  int l;            //dummy iterator in bond order loop
  for (i=0;i<data->nbonds;i++) {
    j=data->from[i];
    k=data->to[i];
    o=data->bondorder[i];
    numcons[j]++;
    numcons[k]++;
    if (numcons[j]>6) {
      printf("bgfplugin) Warning: Bond overflow. Not all bonds were written\n");
      numcons[j]--;
      numcons[k]--;
      continue;
    }
       
    if (numcons[k]>6) {
      printf("bgfplugin) Warning: Bond overflow. Not all bonds were written\n");
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
    fprintf(data->file,"\n");
    printorder = false;
    for (l=0;l<numcons[i];l++) {
      if (orders[6*i+l] != 1.0) {
        printorder = true;
      }
    }
    if (printorder) {
      fprintf(data->file,"ORDER %6i",i);
      for (j=0;j<numcons[i];j++) {
        fprintf(data->file,"%6i",int(orders[6*i+j]));
      }
      fprintf(data->file,"\n");
    }
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

static int write_bonds(void *v, int nbonds, int *fromptr, int *toptr, float *bondorderptr) {
  bgfdata *data = (bgfdata *)v;
  data->from = new int[nbonds];
  data->to = new int[nbonds];
  data->bondorder = new float[nbonds];

  //set the pointers for use later
  for (int i=0;i<nbonds;i++) {
    data->from[i]=fromptr[i];
    data->to[i]=toptr[i];
    data->bondorder[i]=bondorderptr[i];
  }

  data->nbonds = nbonds;
  return MOLFILE_SUCCESS;
}

static void close_bgf_write(void *mydata) {
  bgfdata *data = (bgfdata *)mydata;
  if (data) {
    if (data->file != NULL) fclose(data->file);

    if (data->atomlist != NULL) free(data->atomlist);
    if (data->from != NULL) free(data->from);
    if (data->to != NULL) free(data->to);
    if (data->bondorder != NULL) free(data->bondorder);
    free(data);
  }
}

//
// Free the memory used by the bgf structure
static void close_bgf_read(void *v) {
  bgfdata *bgf = (bgfdata *)v;
  if (bgf) {
    if (bgf->file != NULL) fclose(bgf->file);
    if (bgf->from != NULL) free(bgf->from);
    if (bgf->to != NULL)   free(bgf->to);
    if (bgf->bondorder != NULL)   free(bgf->bondorder);
    delete bgf;
  }
}


static molfile_plugin_t bgfplugin = {
  vmdplugin_ABIVERSION,
  MOLFILE_PLUGIN_TYPE,                      
  "bgf",                                    
  "MSI Biograf Format",
  "Peter Freddolino ",    
  0,                                        
  10,                                        
  VMDPLUGIN_THREADSAFE,                     
  "bgf",
  open_bgf_read,
  read_bgf_structure,
  read_bonds,
  read_bgf_timestep,
  close_bgf_read,
  open_bgf_write,
  write_bgf_structure,
  write_bgf_timestep,
  close_bgf_write,
  0,                            
  0,                            
  0,
  0,
  write_bonds  
};

VMDPLUGIN_EXTERN int VMDPLUGIN_init() {
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_EXTERN int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)&bgfplugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_EXTERN int VMDPLUGIN_fini() {
  return VMDPLUGIN_SUCCESS;
}


