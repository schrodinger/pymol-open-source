/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: ReadPARM7.h,v $
 *      $Author: johns $        $Locker:  $                $State: Exp $
 *      $Revision: 1.28 $      $Date: 2009/03/06 03:24:53 $
 *
 ***************************************************************************
 * DESCRIPTION:
 * NOTE:: Significant modifications were made to the VMD version of 
 *        Bill Ross's original code in order to make it easy to hook 
 *        into VMD plugin structures.  
 *        Further modifications were made to the VMD code to 
 *        read amber 7 parm files, courtesy of Brian Bennion
 * Here is what has changed:
 *     Functions became Class Methods, data became instance variables
 *     The Code to check for compressed files before opening was disabled
 *     Methods get_parm7_atom, get_parm7_bond, get_hydrogen_bond,
 *     get_parm7_natoms, get_parm7_nbonds, get_parm7_boxInfo were added in 
 *     order to convert from prm.c parlance to VMD conventions.
 ***************************************************************************/

/*
 * COPYRIGHT 1992, REGENTS OF THE UNIVERSITY OF CALIFORNIA
 *
 *  prm.c - read information from an amber PARM topology file:
 *      atom/residue/bond/charge info, plus force field data.
 *      This file and the accompanying prm.h may be distributed
 *      provided this notice is retained unmodified and provided
 *      that any modifications to the rest of the file are noted
 *      in comments.
 *
 *      Bill Ross, UCSF 1994
 */

#ifndef READPARM7_H
#define READPARM7_H

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <string.h>
#include "molfile_plugin.h"  // needed for molfile return codes etc

#if defined(WIN32) || defined(WIN64)
#define strcasecmp stricmp
#endif

#if 0 
#define _REAL   double
#define DBLFMT  "%lf"
#else
#define _REAL   float
#define DBLFMT  "%f"
#endif


typedef struct parm {
  char title[85];
  char version[85];
  int   IfBox, Nmxrs, IfCap,
        Natom,  Ntypes,  Nbonds, Nbonh,  Mbona,  Ntheth,  Mtheta, 
        Nphih,  Mphia,  Nhparm, Nparm, Nnb, Nres,Mptra,
        Nbona,  Ntheta,  Nphia,  Numbnd,  Numang,  Nptra,Jparm,
        Natyp,  Nphb, Nat3, Ntype2d, Nttyp, Nspm, Iptres, Nspsol,
        Ipatm, Natcap,Ifpert,Nbper,Ngper,Ndper,Mbper,Mgper,Mdper,
        Numextra;
  _REAL Box[3], Cutcap, Xcap, Ycap, Zcap;
} parmstruct;

static int read_parm7_flag(FILE *file, const char *flag, const char *format) {
  char buf[1024];
    
  /* read the %FLAG text */
  fscanf(file, "%s\n", buf);
  if (strcmp("%FLAG", buf)) {
    printf("AMBER 7 parm read error, at flag section %s,\n", flag);
    printf("        expected %%FLAG but got %s\n", buf);
    return 0; /* read of flag data failed */
  }

  /* read field name specifier */
  fscanf(file, "%s\n", buf);
  if (flag != NULL) {
    if (strcmp(flag, buf)) {
      printf("AMBER 7 parm read error at flag section %s,\n", flag);
      printf("      expected flag field %s but got %s\n", flag, buf);
      return 0; /* read of flag data failed */
    }
  }

  /* read format string */
  fscanf(file, "%s\n", buf);
  if (format != NULL) {
    if (strcmp(format, buf)) {
      if (!strcmp(flag, "TITLE") && !strcmp(format, "%FORMAT(20a4)")) {
        if (!strcmp(buf, "%FORMAT(a80)"))
          return 1; /* accept a80 as substitute for 20a4 */
      }
      printf("AMBER 7 parm read error at flag section %s,\n", flag);
      printf("      expected format %s but got %s\n", format, buf);
      return 0; /* read of flag data failed */
    }
  }

  return 1; /* read of flag data succeeded */
}

/*
 *  open_parm7_file() - fopen regular or popen compressed file for reading
 *  Return FILE handle on success.
 *  set as_pipe to 1 if opened with popen, or 0 if opened with fopen.
 */

static FILE *open_parm7_file(const char *name, int *as_pipe) {
  struct stat buf;
  char cbuf[8192];
  int length;
  int &compressed = *as_pipe;
  FILE *fp;

  length = strlen(name);
  compressed = 0;  // Just to start
  strcpy(cbuf, name);

  /*
   *  if file doesn't exist, maybe it has been compressed/decompressed
   */
  if (stat(cbuf, &buf) == -1) {
    switch (errno) {
      case ENOENT:
        if (!compressed) {
          strcat(cbuf, ".Z");
          if (stat(cbuf, &buf) == -1) {
            printf("%s, %s: does not exist\n", name, cbuf);
            return(NULL);
          }
          compressed++;

          // Don't modify the filename
          //strcat(name, ".Z"); /* TODO: add protection */
        } else {
          cbuf[length-2] = '\0';
          if (stat(cbuf, &buf) == -1) {
            printf("%s, %s: does not exist\n", name, cbuf);
            return(NULL);
          }
          compressed = 0;
        }
        break;

      default:
        return(NULL);
    }
  }

  /*
   *  open the file
   */
#if defined(_MSC_VER)
  if (compressed) {
    /* NO "zcat" on Win32 */
    printf("Cannot load compressed PARM files on Windows.\n");
    return NULL;
  }
#else
  if (compressed) {
    char pcmd[120];
    sprintf(pcmd, "zcat %s", cbuf);
    if ((fp = popen(pcmd, "r")) == NULL) {
      perror(pcmd);
      return NULL;
    }
  }
#endif
  else {
    if ((fp = fopen(cbuf, "r")) == NULL) {
      perror(cbuf);
      return NULL;
    }
  }

  return(fp);
}


static int parse_parm7_atoms(const char *fmt, 
    int natoms, molfile_atom_t *atoms, FILE *file) {
  char buf[85];

  if (strcasecmp(fmt, "%FORMAT(20a4)"))
    return 0;

  int j=0;
  for (int i=0; i<natoms; i++) {
    molfile_atom_t *atom = atoms+i;
    if (!(i%20)) {
      j=0;
      fgets(buf, 85, file);
    }
    strncpy(atom->name, buf+4*j, 4);
    atom->name[4]='\0';
    j++;
  }
  return 1;
}


static int parse_parm7_charge(const char *fmt, 
    int natoms, molfile_atom_t *atoms, FILE *file) {
  if (strcasecmp(fmt, "%FORMAT(5E16.8)")) 
    return 0;

  for (int i=0; i<natoms; i++) {
    double q=0;
    if (fscanf(file, " %lf", &q) != 1) {
      fprintf(stderr, "PARM7: error reading charge at index %d\n", i);
      return 0;
    }
    atoms[i].charge = 0.0548778 * (float)q; /* convert to elementary charge units */
  }

  return 1;
}


static int parse_parm7_mass(const char *fmt,
    int natoms, molfile_atom_t *atoms, FILE *file) {
  if (strcasecmp(fmt, "%FORMAT(5E16.8)")) return 0;
  for (int i=0; i<natoms; i++) {
    double m=0;
    if (fscanf(file, " %lf", &m) != 1) {
      fprintf(stderr, "PARM7: error reading mass at index %d\n", i);
      return 0;
    }
    atoms[i].mass = (float)m;
  }
  return 1;
}


static int parse_parm7_atype(const char *fmt,
    int natoms, molfile_atom_t *atoms, FILE *file) {
  if (strcasecmp(fmt, "%FORMAT(20a4)")) return 0;
  char buf[85];
  int j=0;
  for (int i=0; i<natoms; i++) {
    molfile_atom_t *atom = atoms+i;
    if (!(i%20)) {
      j=0;
      fgets(buf, 85, file);
    }
    strncpy(atom->type, buf+4*j, 4);
    atom->type[4]='\0';
    j++;
  }
  return 1;
}


static int parse_parm7_resnames(const char *fmt,
    int nres, char *resnames, FILE *file) {
  if (strcasecmp(fmt, "%FORMAT(20a4)")) return 0;
  char buf[85];
  int j=0;
  for (int i=0; i<nres; i++) {
    if (!(i%20)) {
      j=0;
      fgets(buf, 85, file);
    }
    strncpy(resnames, buf+4*j, 4);
    resnames += 4;
    j++;
  }

  return 1;
}

static int parse_parm7_respointers(const char *fmt, int natoms, 
    molfile_atom_t *atoms, int nres, const char *resnames, FILE *file) {
  if (strcasecmp(fmt, "%FORMAT(10I8)")) return 0;
  int cur, next;
  fscanf(file, " %d", &cur);
  for (int i=1; i<nres; i++) {
    if (fscanf(file, " %d", &next) != 1) {
      fprintf(stderr, "PARM7: error reading respointer records at residue %d\n", i);
      return 0;
    }
    while (cur < next) {
      if (cur > natoms) {
        fprintf(stderr, "invalid atom index: %d\n", cur);
        return 0;
      }
      strncpy(atoms[cur-1].resname, resnames, 4);
      atoms[cur-1].resname[4] = '\0';
      atoms[cur-1].resid = i;
      cur++;
    }
    resnames += 4;
  }
  // store the last residue name
  while (cur <= natoms) {
    strncpy(atoms[cur-1].resname, resnames, 4);
    atoms[cur-1].resname[4] = '\0';
    atoms[cur-1].resid = nres;
    cur++;
  }

  return 1;
}


static int parse_parm7_bonds(const char *fmt,
    int nbonds, int *from, int *to, FILE *file) {
  if (strcasecmp(fmt, "%FORMAT(10I8)")) 
    return 0;

  int a, b, tmp;
  for (int i=0; i<nbonds; i++) {
    if (fscanf(file, " %d %d %d", &a, &b, &tmp) != 3) {
      fprintf(stderr, "PARM7: error reading bond number %d\n", i);
      return 0;
    }
    from[i] = a/3 + 1;
    to[i]   = b/3 + 1;
  }

  return 1;
}


/*
 *  close_parm7_file() - close fopened or popened file
 */
static void close_parm7_file(FILE *fileptr, int popn) {
#if defined(_MSC_VER)
  if (popn) {
    printf("pclose() no such function on win32!\n");
  } else {
   if (fclose(fileptr) == -1)
     perror("fclose");
  }
#else
  if (popn) {
    if (pclose(fileptr) == -1)
      perror("pclose");
  } else {
    if (fclose(fileptr) == -1)
      perror("fclose");
  }
#endif
}

static const char *parm7 = "%8d%8d%8d%8d%8d%8d%8d%8d%8d%8d\n";

static parmstruct *read_parm7_header(FILE *file) {
  char sdum[512]; 
  parmstruct *prm;
  prm = new parmstruct;

  /* READ VERSION */
  fgets(sdum, 512, file);

  /* READ TITLE */
  if (!read_parm7_flag(file, "TITLE", "%FORMAT(20a4)")) {
    delete prm;
    return NULL;
  }

  // read the title string itself, and handle empty lines
#if 0
  // XXX this code fails with some AMBER 9 test files
  fscanf(file, "%s\n", prm->title);
#else
  // XXX this hack causes AMBER 9 prmtop files to load
  fgets(prm->title, sizeof(prm->title), file);
#endif

  if (strstr(prm->title, "%FLAG") == NULL) {
    // Got a title string
    if (!read_parm7_flag(file, "POINTERS", "%FORMAT(10I8)")) {
      delete prm;
      return NULL;
    }
  } else {
    // NO title string, use a special method to pick up next flag
#if 0
    fscanf(file,"%s\n", sdum);
    if (strcmp("POINTERS", sdum)) {
      printf("AMBER 7 parm read error at flag section POINTERS\n");
      printf("      expected flag field POINTERS but got %s\n", sdum);
#else
    if (strstr(prm->title, "POINTERS") == NULL) {
      printf("AMBER 7 parm read error at flag section POINTERS\n");
      printf("      expected flag field POINTERS but got %s\n", prm->title);
#endif
      delete prm;
      return NULL;
    }
#if 0
    fscanf(file,"%s\n", sdum);
    if (strcasecmp("%FORMAT(10I8)", sdum)) {
#else
    fgets(sdum, sizeof(sdum), file);
    if ((strstr(sdum, "%FORMAT(10I8)") == NULL) &&
        (strstr(sdum, "%FORMAT(10i8)") == NULL)) {
#endif
      printf("AMBER 7 parm read error at flag section POINTERS,\n");
      printf("      expected format %%FORMAT(10I8) but got %s\n", sdum);
      delete prm;
      return NULL;
    }
  }

  /* READ POINTERS (CONTROL INTEGERS) */
  fscanf(file,parm7,
         &prm->Natom,  &prm->Ntypes, &prm->Nbonh, &prm->Nbona,
         &prm->Ntheth, &prm->Ntheta, &prm->Nphih, &prm->Nphia,
         &prm->Jparm,  &prm->Nparm);
  fscanf(file, parm7,  
         &prm->Nnb,   &prm->Nres,   &prm->Mbona,  &prm->Mtheta,
         &prm->Mphia, &prm->Numbnd, &prm->Numang, &prm->Mptra,
         &prm->Natyp, &prm->Nphb);
  fscanf(file, parm7,  &prm->Ifpert, &prm->Nbper,  &prm->Ngper,
         &prm->Ndper, &prm->Mbper,  &prm->Mgper, &prm->Mdper,
         &prm->IfBox, &prm->Nmxrs,  &prm->IfCap);

  fscanf(file,"%8d",&prm->Numextra); //BB
  prm->Nptra=prm->Mptra; //BB new to amber 7 files...

  prm->Nat3 = 3 * prm->Natom;
  prm->Ntype2d = prm->Ntypes * prm->Ntypes;
  prm->Nttyp = prm->Ntypes*(prm->Ntypes+1)/2;

  return prm;
}


#endif
