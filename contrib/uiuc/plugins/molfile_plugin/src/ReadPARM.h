/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: ReadPARM.h,v $
 *      $Author: johns $        $Locker:  $                $State: Exp $
 *      $Revision: 1.13 $      $Date: 2006/07/26 22:03:05 $
 *
 ***************************************************************************
 * DESCRIPTION:
 * NOTE:: Significant were made to the VMD version of
 *        Bill Ross's original code in order to make it easy to hook
 *        into VMD plugin structures.
 * Here is what was changed:
 *     Functions became Class Methods, data became instance variables
 *     The Code to check for compressed files before opening was disabled
 *     Methods get_parm_atom, get_parm_bond, get_hydrogen_bond,
 *     get_parm_natoms, get_parm_nbonds, get_parm_boxInfo were added in 
 *     order to convert from prm.c parlance to VMD conventions.
 *     RCS Information headers and footers were added.
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

#ifndef READPARM_H
#define READPARM_H

// XXX enable the new AMBER reading code, deals with packed 
// fortran 12I6 integer formats.
#define USENEWCODE 1

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <string.h>
#include "molfile_plugin.h"  // needed for molfile return codes etc

#if 0 
#define _REAL		double
#define DBLFMT          "%lf"
#else
#define _REAL		float
#define DBLFMT          "%f"
#endif

typedef struct parm {
	char	ititl[81];
	int 	IfBox, Nmxrs, IfCap,
		 Natom,  Ntypes,  Nbonh,  Mbona,  Ntheth,  Mtheta, 
		 Nphih,  Mphia,  Nhparm, Nparm, Nnb, Nres,
		 Nbona,  Ntheta,  Nphia,  Numbnd,  Numang,  Nptra,
		 Natyp,  Nphb, Nat3, Ntype2d, Nttyp, Nspm, Iptres, Nspsol,
		 Ipatm, Natcap;
	char 	*AtomNames, *ResNames, *AtomSym, *AtomTree;
	_REAL	*Charges, *Masses, *Rk, *Req, *Tk, *Teq, *Pk, *Pn, *Phase,
		 *Solty, *Cn1, *Cn2, *HB12, *HB6;
	_REAL	Box[3], Cutcap, Xcap, Ycap, Zcap;
	int 	*Iac, *Iblo, *Cno, *Ipres, *ExclAt, *TreeJoin, 
		*AtomRes, *BondHAt1, *BondHAt2, *BondHNum, *BondAt1, *BondAt2, 
		*BondNum, *AngleHAt1, *AngleHAt2, *AngleHAt3, *AngleHNum, 
		*AngleAt1, *AngleAt2, *AngleAt3, *AngleNum, *DihHAt1, 
		*DihHAt2, *DihHAt3, *DihHAt4, *DihHNum, *DihAt1, *DihAt2, 
		*DihAt3, *DihAt4, *DihNum, *Boundary;
} parmstruct;

// put the class in an anonymous namespace to give it internal linkage
namespace {
  class ReadPARM {
  public:
    ReadPARM() {popn = 0;}
    ~ReadPARM(void) {}
 
    int popn;
    parmstruct *prm;
    FILE *open_parm_file(const char *name);
    void close_parm_file(FILE *fileptr);
    char *get(int size);
    int preadln(FILE *file, char *string);
    int readparm(FILE *file);
    void get_parm_atom(int, char *, char *, char *, char *, int *, float *,
                       float *);

    void get_parm_bond(int, int fromAtom[], int toAtom[]);
    void get_hydrogen_bond(int, int fromAtom[], int toAtom[]);
    int get_parm_natoms();
    int get_parm_nbonds();
    int get_parm_boxInfo();
    int read_fortran_12I6(FILE *fp, int *data, int count);  
  };
}


/*	fortran formats 
 *	 9118 FORMAT(12I6)
 *	 9128 FORMAT(5E16.8)
static char	*f9118 = (char *) "%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d\n";
 */

static int readtoeoln(FILE *f) {
  int c;

  /* skip to eoln */
  while((c = getc(f)) != '\n') {
    if (c == EOF) 
      return -1;
  }

  return 0;
}  

/***********************************************************************
 			    	open_parm_file()
************************************************************************/

/*
 *  open_parm_file() - fopen regular or popen compressed file for reading
 */

FILE *ReadPARM::open_parm_file(const char *name) {
	struct stat	buf;
	char		cbuf[120];
	int		length, compressed;
	FILE		*fp;

	length = strlen(name);
	compressed = 0;  // Just to start
//	compressed = iscompressed(name);
	strcpy(cbuf, name);

	/*
	 *  if file doesn't exist, maybe it has been compressed/decompressed
	 */

	if (stat(cbuf, &buf) == -1) {
		switch (errno) {
		case ENOENT:	{
			if (!compressed) {
				strcat(cbuf, ".Z");
				if (stat(cbuf, &buf) == -1) {
					printf("%s, %s: does not exist\n", 
						name, cbuf);
					return(NULL);
				}
				compressed++;
                                // Don't modify the filename
				//strcat(name, ".Z"); /* TODO: add protection */
			} else {
				cbuf[length-2] = '\0';
				if (stat(cbuf, &buf) == -1) {
					printf("%s, %s: does not exist\n", 
							name, cbuf);
					return(NULL);
				}
				compressed = 0;
			}
			break;
		}
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
		popn = 1;

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

/***********************************************************************
 			    close_parm_file   
************************************************************************/

/*
 *  close_parm_file() - close fopened or popened file
 */
void ReadPARM::close_parm_file(FILE *fileptr) {
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


/***********************************************************************
                                      GET()
************************************************************************/
char *ReadPARM::get(int size) {
	char	*ptr;

#ifdef DEBUG
	printf("malloc %d\n", size);
	fflush(stdout);
#endif
	if (size ==0)
		return((char *) NULL);

	if ((ptr = (char *) malloc((unsigned)size)) == NULL) {
		printf("malloc %d", size);
		fflush(stdout);
		perror("malloc err:");
        fprintf(stderr, "Exiting due to ReadPARM memory allocation error.\n");
	}
	return(ptr);
}

/***********************************************************************
 			  	PREADLN()
************************************************************************/
int ReadPARM::preadln(FILE *file, char *string) {
  int i, j;

  for (i=0; i<81; i++) {
    if ((j = getc(file)) == EOF) {
      printf("Error: unexpected EOF in Parm file\n");
      return MOLFILE_ERROR;
    }

    string[i] = (char) j;
    if (string[i] == '\n') {
      break;
    }
  }

  if (i == 80  &&  string[i] != '\n') {
    printf("Error: line too long in Parm file:\n%.80s", string);
    return MOLFILE_ERROR;
  }

  return 0; /* success */
}

/***************************************************************************
		      		READPARM()
****************************************************************************/
/*
 * readparm() - instantiate a given parmstruct
 */
int ReadPARM::readparm(FILE *file) {
        //
        // XXX This code leaks memory every time it's run!  prm is allocated
        // but its data is never freed, even when ReadPARM itself is deleted.
        // I've added exception code so we at least have a chance of 
        // recovering gracefully from a memory allocation error.

	_REAL 		*H;
	int		i, res, ifpert;
        int *buffer; // used for reading fortran integer blocks

	prm = (parmstruct *) get(sizeof(parmstruct));
        if (prm == NULL) {
          return MOLFILE_ERROR; /* allocation failure */
        }

	/* READ TITLE */
	if (preadln(file, prm->ititl) != 0) {
          return MOLFILE_ERROR;
        }

	/* READ CONTROL INTEGERS */
#if !defined(USENEWCODE)
	fscanf(file, f9118, 
		&prm->Natom,  &prm->Ntypes, &prm->Nbonh, &prm->Mbona, 
		&prm->Ntheth, &prm->Mtheta, &prm->Nphih, &prm->Mphia, 
		&prm->Nhparm, &prm->Nparm,  &prm->Nnb,   &prm->Nres);

	fscanf(file, f9118, 
		&prm->Nbona,  &prm->Ntheta, &prm->Nphia, &prm->Numbnd, 
		&prm->Numang, &prm->Nptra,  &prm->Natyp, &prm->Nphb, 
		&ifpert,      &idum,        &idum,       &idum);

	fscanf(file, " %d %d %d %d %d %d", 
		&idum, &idum,&idum,&prm->IfBox,&prm->Nmxrs,&prm->IfCap);
#else
        buffer = new int[30];
        if (!read_fortran_12I6(file,buffer,30)) { 
          return MOLFILE_ERROR;
        }
        prm->Natom = buffer[0];
        prm->Ntypes = buffer[1];
        prm->Nbonh = buffer[2];
        prm->Mbona = buffer[3];
        prm->Ntheth = buffer[4];
        prm->Mtheta = buffer[5];
        prm->Nphih = buffer[6];
        prm->Mphia = buffer[7];
        prm->Nhparm = buffer[8];
        prm->Nparm = buffer[9];
        prm->Nnb = buffer[10];
        prm->Nres = buffer[11];
        prm->Nbona = buffer[12];
        prm->Ntheta = buffer[13];
        prm->Nphia = buffer[14];
        prm->Numbnd = buffer[15];
        prm->Numang = buffer[16];
        prm->Nptra = buffer[17];
        prm->Natyp = buffer[18];
        prm->Nphb = buffer[19];
        ifpert = buffer[20];
        // items 21 through 26 are ignored currently.
        prm->IfBox = buffer[27];
        prm->Nmxrs = buffer[28];
        prm->IfCap = buffer[29];
        delete [] buffer;
#endif
        readtoeoln(file);


       if (ifpert) {
         printf("not equipped to read perturbation prmtop\n");
         free(prm);
         return MOLFILE_ERROR;
       }


	/* ALLOCATE MEMORY */
	prm->Nat3 = 3 * prm->Natom;
	prm->Ntype2d = prm->Ntypes * prm->Ntypes;
	prm->Nttyp = prm->Ntypes*(prm->Ntypes+1)/2;

	/*
	 * get most of the indirect stuff; some extra allowed for char arrays
	 */
	prm->AtomNames = (char *) get(4*prm->Natom+81);
	prm->Charges = (_REAL *) get(sizeof(_REAL)*prm->Natom);
	prm->Masses = (_REAL *) get(sizeof(_REAL)*prm->Natom);
	prm->Iac = (int *) get(sizeof(int)*prm->Natom);
	prm->Iblo = (int *) get(sizeof(int)*prm->Natom);
	prm->Cno = (int *) get(sizeof(int)* prm->Ntype2d);
	prm->ResNames = (char *) get(4* prm->Nres+81);
	prm->Ipres = (int *) get(sizeof(int)*( prm->Nres+1));
        prm->Rk = (_REAL *) get(sizeof(_REAL)* prm->Numbnd);
        prm->Req = (_REAL *) get(sizeof(_REAL)* prm->Numbnd);
        prm->Tk = (_REAL *) get(sizeof(_REAL)* prm->Numang);
        prm->Teq = (_REAL *) get(sizeof(_REAL)* prm->Numang);
        prm->Pk = (_REAL *) get(sizeof(_REAL)* prm->Nptra);
        prm->Pn = (_REAL *) get(sizeof(_REAL)* prm->Nptra);
        prm->Phase = (_REAL *) get(sizeof(_REAL)* prm->Nptra);
        prm->Solty = (_REAL *) get(sizeof(_REAL)* prm->Natyp);
        prm->Cn1 = (_REAL *) get(sizeof(_REAL)* prm->Nttyp);
        prm->Cn2 = (_REAL *) get(sizeof(_REAL)* prm->Nttyp);
	prm->BondHAt1 = (int *) get(sizeof(int)* prm->Nbonh);
	prm->BondHAt2 = (int *) get(sizeof(int)* prm->Nbonh);
	prm->BondHNum = (int *) get(sizeof(int)* prm->Nbonh);
	prm->BondAt1 = (int *) get(sizeof(int)* prm->Nbona);
	prm->BondAt2 = (int *) get(sizeof(int)* prm->Nbona);
	prm->BondNum = (int *) get(sizeof(int)* prm->Nbona);
	prm->AngleHAt1 = (int *) get(sizeof(int)* prm->Ntheth);
	prm->AngleHAt2 = (int *) get(sizeof(int)* prm->Ntheth);
	prm->AngleHAt3 = (int *) get(sizeof(int)* prm->Ntheth);
	prm->AngleHNum = (int *) get(sizeof(int)* prm->Ntheth);
	prm->AngleAt1 = (int *) get(sizeof(int)* prm->Ntheta);
	prm->AngleAt2 = (int *) get(sizeof(int)*prm->Ntheta);
	prm->AngleAt3 = (int *) get(sizeof(int)*prm->Ntheta);
	prm->AngleNum = (int *) get(sizeof(int)*prm->Ntheta);
	prm->DihHAt1 = (int *) get(sizeof(int)*prm->Nphih);
	prm->DihHAt2 = (int *) get(sizeof(int)*prm->Nphih);
	prm->DihHAt3 = (int *) get(sizeof(int)*prm->Nphih);
	prm->DihHAt4 = (int *) get(sizeof(int)*prm->Nphih);
	prm->DihHNum = (int *) get(sizeof(int)*prm->Nphih);
	prm->DihAt1 = (int *) get(sizeof(int)*prm->Nphia);
	prm->DihAt2 = (int *) get(sizeof(int)*prm->Nphia);
	prm->DihAt3 = (int *) get(sizeof(int)*prm->Nphia);
	prm->DihAt4 = (int *) get(sizeof(int)*prm->Nphia);
	prm->DihNum = (int *) get(sizeof(int)*prm->Nphia);
	prm->ExclAt = (int *) get(sizeof(int)*prm->Nnb);
	prm->HB12 = (_REAL *) get(sizeof(_REAL)*prm->Nphb);
	prm->HB6 = (_REAL *) get(sizeof(_REAL)*prm->Nphb);
	prm->AtomSym = (char *) get(4*prm->Natom+81);
	prm->AtomTree = (char *) get(4*prm->Natom+81);
	prm->TreeJoin = (int *) get(sizeof(int)*prm->Natom);
	prm->AtomRes = (int *) get(sizeof(int)*prm->Natom);

	/* 
	 * READ ATOM NAMES -IH(M04)
	 */
	for (i=0; i<(prm->Natom/20 + (prm->Natom%20 ? 1 : 0)); i++) {
		preadln(file, &prm->AtomNames[i*80]);
        }

	/* 
	 * READ ATOM CHARGES -X(L15)
	 *	(pre-multiplied by an energy factor of 18.2223 == sqrt(332)
	 *	 for faster force field calculations)
	 */
	for (i=0; i<prm->Natom; i++) {
          fscanf(file, " " DBLFMT, &prm->Charges[i]);
          prm->Charges[i] *= 0.0548778; /* convert back to elementary charge units */
        }
        readtoeoln(file);

	/* 
	 * READ ATOM MASSES -X(L20)
	 */
	for (i=0; i<prm->Natom; i++)
	  fscanf(file, " " DBLFMT, &prm->Masses[i]);
        readtoeoln(file);

	/* 
	 * READ ATOM L-J TYPES -IX(I04)
	 */
#if !defined(USENEWCODE)
	for (i=0; i<prm->Natom; i++) 
          fscanf(file, " %d", &prm->Iac[i]);
#else
        if (!read_fortran_12I6(file, prm->Iac, prm->Natom)) {
          return MOLFILE_ERROR;
        }
#endif
        readtoeoln(file);

	/* 
	 * READ ATOM INDEX TO 1st IN EXCLUDED ATOM LIST "NATEX" -IX(I08)
	 */
#if !defined(USENEWCODE)
	for (i=0; i<prm->Natom; i++)
          fscanf(file, " %d", &prm->Iblo[i]);
#else
        if (!read_fortran_12I6(file, prm->Iblo, prm->Natom)) {
          return MOLFILE_ERROR;
        }
#endif
        readtoeoln(file);

	/* 
	 * READ TYPE INDEX TO N-B TYPE -IX(I06)
	 */
#if !defined(USENEWCODE)
	for (i=0; i<prm->Ntype2d; i++)
          fscanf(file, " %d", &prm->Cno[i]);
#else
        if (!read_fortran_12I6(file, prm->Cno, prm->Ntype2d)) {
          return MOLFILE_ERROR;
        } 
#endif
        readtoeoln(file);

	/* 
	 * READ RES NAMES (4 chars each, 4th blank) -IH(M02)
	 */
	for (i=0; i<(prm->Nres/20 + (prm->Nres%20 ? 1 : 0)); i++)
          preadln(file, &prm->ResNames[i*80]);

	/* 
	 * READ RES POINTERS TO 1st ATOM 		-IX(I02)
	 */
#if !defined(USENEWCODE)
	for (i=0; i<prm->Nres; i++) 
          fscanf(file, " %d", &prm->Ipres[i]);
#else
        if (!read_fortran_12I6(file, prm->Ipres, prm->Nres)) {
          return MOLFILE_ERROR;
        }
#endif
	prm->Ipres[prm->Nres] = prm->Natom + 1;
        readtoeoln(file);

	/* 
	 * READ BOND FORCE CONSTANTS 			-RK()
	 */
	for (i=0; i< prm->Numbnd; i++) 
          fscanf(file, " " DBLFMT, &prm->Rk[i]);
        readtoeoln(file);

	/* 
	 * READ BOND LENGTH OF MINIMUM ENERGY  		-REQ()
	 */
	for (i=0; i< prm->Numbnd; i++) 
          fscanf(file, " " DBLFMT, &prm->Req[i]);
        readtoeoln(file);

	/* 
	 * READ BOND ANGLE FORCE CONSTANTS (following Rk nomen) -TK()
	 */
	for (i=0; i< prm->Numang; i++) 
          fscanf(file, " " DBLFMT, &prm->Tk[i]);
        readtoeoln(file);

	/* 
	 * READ BOND ANGLE OF MINIMUM ENERGY (following Req nomen) -TEQ()
	 */
	for (i=0; i< prm->Numang; i++) 
          fscanf(file, " " DBLFMT, &prm->Teq[i]);
        readtoeoln(file);

	/* 
	 * READ DIHEDRAL PEAK MAGNITUDE 		-PK()
	 */
	for (i=0; i< prm->Nptra; i++) 
          fscanf(file, " " DBLFMT, &prm->Pk[i]);
        readtoeoln(file);

	/* 
	 * READ DIHEDRAL PERIODICITY 			-PN()
	 */
	for (i=0; i< prm->Nptra; i++) 
          fscanf(file, " " DBLFMT, &prm->Pn[i]);
        readtoeoln(file);

	/* 
	 * READ DIHEDRAL PHASE  			-PHASE()
	 */
	for (i=0; i< prm->Nptra; i++) 
          fscanf(file, " " DBLFMT, &prm->Phase[i]);
        readtoeoln(file);

	/* 
	 * ?? "RESERVED" 				-SOLTY()
	 */
	for (i=0; i< prm->Natyp; i++) 
          fscanf(file, " " DBLFMT, &prm->Solty[i]);
        readtoeoln(file);

	/* 
	 * READ L-J R**12 FOR ALL PAIRS OF ATOM TYPES  	-CN1()
	 *	(SHOULD BE 0 WHERE H-BONDS)
	 */
	for (i=0; i< prm->Nttyp; i++) 
          fscanf(file, " " DBLFMT, &prm->Cn1[i]);
        readtoeoln(file);

	/* 
	 * READ L-J R**6 FOR ALL PAIRS OF ATOM TYPES 	-CN2()
	 *	(SHOULD BE 0 WHERE H-BONDS)
	 */
	for (i=0; i< prm->Nttyp; i++) 
          fscanf(file, " " DBLFMT, &prm->Cn2[i]);
        readtoeoln(file);

	/* 
	 * READ COVALENT BOND W/ HYDROGEN (3*(atnum-1)): 
	 *	IBH = ATOM1 		-IX(I12)
	 *	JBH = ATOM2 		-IX(I14)
	 *	ICBH = BOND ARRAY PTR	-IX(I16)
	 */
#if !defined(USENEWCODE)
	for (i=0; i<prm->Nbonh; i++) 
          fscanf(file, " %d %d %d", 
                 &prm->BondHAt1[i], &prm->BondHAt2[i], &prm->BondHNum[i]);
#else
        buffer = new int[3*prm->Nbonh];
        if (!read_fortran_12I6(file, buffer, 3*prm->Nbonh)) { 
          return MOLFILE_ERROR;
        }
        for (i=0; i<prm->Nbonh; i++) { 
          prm->BondHAt1[i] = buffer[3*i];
          prm->BondHAt2[i] = buffer[3*i+1];
          prm->BondHNum[i] = buffer[3*i+2];
        }
        delete [] buffer;
#endif
        readtoeoln(file);

	/* 
	 * READ COVALENT BOND W/OUT HYDROGEN (3*(atnum-1)):
	 *	IB = ATOM1		-IX(I18)
	 *	JB = ATOM2		-IX(I20)
	 *	ICB = BOND ARRAY PTR	-IX(I22)
	 */
#if !defined(USENEWCODE)
	for (i=0; i<prm->Nbona; i++)
          fscanf(file, " %d %d %d", 
                 &prm->BondAt1[i], &prm->BondAt2[i], &prm->BondNum[i]);
#else
        buffer = new int[3*prm->Nbona];
        if (!read_fortran_12I6(file, buffer, 3*prm->Nbona)) { 
          return MOLFILE_ERROR;
        }
        for (i=0; i<prm->Nbona; i++) { 
          prm->BondAt1[i] = buffer[3*i];
          prm->BondAt2[i] = buffer[3*i+1];
          prm->BondNum[i] = buffer[3*i+2];
        }
        delete [] buffer;
#endif
        readtoeoln(file);

	/* 
	 * READ ANGLE W/ HYDROGEN: 
	 *	ITH = ATOM1			-IX(I24)
	 *	JTH = ATOM2			-IX(I26)
	 *	KTH = ATOM3			-IX(I28)
	 *	ICTH = ANGLE ARRAY PTR		-IX(I30)
	 */
#if !defined(USENEWCODE)
	for (i=0; i<prm->Ntheth; i++)
          fscanf(file, " %d %d %d %d", 
                 &prm->AngleHAt1[i], &prm->AngleHAt2[i], 
                 &prm->AngleHAt3[i], &prm->AngleHNum[i]);
#else
        buffer = new int[4*prm->Ntheth];
        if (!read_fortran_12I6(file, buffer, 4*prm->Ntheth)) { 
          return MOLFILE_ERROR;
        }
        for (i=0; i<prm->Ntheth; i++) { 
          prm->AngleHAt1[i] = buffer[4*i];
          prm->AngleHAt2[i] = buffer[4*i+1];
          prm->AngleHAt3[i] = buffer[4*i+2];
          prm->AngleHNum[i] = buffer[4*i+3];
        }
        delete [] buffer;
#endif
        readtoeoln(file);

	/* 
	 * READ ANGLE W/OUT HYDROGEN: 
	 *	IT = ATOM1			-IX(I32)
	 *	JT = ATOM2			-IX(I34)
	 *	KT = ATOM3			-IX(I36)
	 *	ICT = ANGLE ARRAY PTR		-IX(I38)
	 */
#if !defined(USENEWCODE)
	for (i=0; i<prm->Ntheta; i++)
          fscanf(file, " %d %d %d %d", 
                 &prm->AngleAt1[i], &prm->AngleAt2[i], 
                 &prm->AngleAt3[i], &prm->AngleNum[i]);
#else
        buffer = new int[4*prm->Ntheta];
        if (!read_fortran_12I6(file, buffer, 4*prm->Ntheta)) {
          return MOLFILE_ERROR;
        }
        for (i=0; i<prm->Ntheta; i++) { 
          prm->AngleAt1[i] = buffer[4*i];
          prm->AngleAt2[i] = buffer[4*i+1];
          prm->AngleAt3[i] = buffer[4*i+2];
          prm->AngleNum[i] = buffer[4*i+3];
        }
        delete [] buffer;
#endif
        readtoeoln(file);

	/* 
	 * READ DIHEDRAL W/ HYDROGEN: 
	 *	ITH = ATOM1			-IX(40)
	 *	JTH = ATOM2			-IX(42)
	 *	KTH = ATOM3			-IX(44)
	 *	LTH = ATOM4			-IX(46)
	 *	ICTH = DIHEDRAL ARRAY PTR	-IX(48)
	 */
#if !defined(USENEWCODE)
	for (i=0; i<prm->Nphih; i++)
          fscanf(file, " %d %d %d %d %d", 
                 &prm->DihHAt1[i], &prm->DihHAt2[i], &prm->DihHAt3[i], 
                 &prm->DihHAt4[i], &prm->DihHNum[i]);
#else
        buffer = new int[5*prm->Nphih];
        if (!read_fortran_12I6(file, buffer, 5*prm->Nphih)) {
          return MOLFILE_ERROR;
        }
        for (i=0; i<prm->Nphih; i++) { 
          prm->DihHAt1[i] = buffer[5*i];
          prm->DihHAt2[i] = buffer[5*i+1];
          prm->DihHAt3[i] = buffer[5*i+2];
          prm->DihHAt4[i] = buffer[5*i+3];
          prm->DihHNum[i] = buffer[5*i+4];
        }
        delete [] buffer;
#endif
        readtoeoln(file);

	/* 
	 * READ DIHEDRAL W/OUT HYDROGEN: 
	 *	IT = ATOM1
	 *	JT = ATOM2
	 *	KT = ATOM3
	 *	LT = ATOM4
	 *	ICT = DIHEDRAL ARRAY PTR
	 */
#if !defined(USENEWCODE)
	for (i=0; i<prm->Nphia; i++) {
          fscanf(file, " %d %d %d %d %d", 
                 &prm->DihAt1[i], &prm->DihAt2[i], &prm->DihAt3[i], 
                 &prm->DihAt4[i], &prm->DihNum[i]);
	}
#else
        buffer = new int[5*prm->Nphia];
        if (!read_fortran_12I6(file, buffer, 5*prm->Nphia)) {
          return MOLFILE_ERROR;
        }
        for (i=0; i<prm->Nphia; i++) {
          prm->DihAt1[i] = buffer[5*i];
          prm->DihAt2[i] = buffer[5*i+1];
          prm->DihAt3[i] = buffer[5*i+2];
          prm->DihAt4[i] = buffer[5*i+3];
          prm->DihNum[i] = buffer[5*i+4];
        }
        delete [] buffer;
#endif
        readtoeoln(file);

	/*
	 * READ EXCLUDED ATOM LIST	-IX(I10)
	 */
#if !defined(USENEWCODE)
	for (i=0; i<prm->Nnb; i++)
		fscanf(file, " %d", &prm->ExclAt[i]);
#else
        if (!read_fortran_12I6(file, prm->ExclAt, prm->Nnb)) {
           return MOLFILE_ERROR;
        }
#endif
        readtoeoln(file);

	/*
	 * READ H-BOND R**12 TERM FOR ALL N-B TYPES	-ASOL()
	 */
	for (i=0; i<prm->Nphb; i++) 
          fscanf(file, " " DBLFMT, &prm->HB12[i]);
        readtoeoln(file);

	/*
	 * READ H-BOND R**6 TERM FOR ALL N-B TYPES	-BSOL()
	 */
	for (i=0; i<prm->Nphb; i++) 
          fscanf(file, " " DBLFMT, &prm->HB6[i]);
        readtoeoln(file);

	/*
	 * READ H-BOND CUTOFF (NOT USED) ??		-HBCUT()
	 */
	H = (_REAL *) get(prm->Nphb * sizeof(_REAL));
	for (i=0; i<prm->Nphb; i++) 
          fscanf(file, " " DBLFMT, &H[i]);
	free((char *)H);
        readtoeoln(file);

	/*
	 * READ ATOM SYMBOLS (FOR ANALYSIS PROGS)	-IH(M06)
	 */
	for (i=0; i<(prm->Natom/20 + (prm->Natom%20 ? 1 : 0)); i++)
          preadln(file, &prm->AtomSym[i*80]);

	/*
	 * READ TREE SYMBOLS (FOR ANALYSIS PROGS)	-IH(M08)
	 */
	for (i=0; i<(prm->Natom/20 + (prm->Natom%20 ? 1 : 0)); i++)
          preadln(file, &prm->AtomTree[i*80]);
      
	/*
	 * READ TREE JOIN INFO (FOR ANALYSIS PROGS)	-IX(I64)
	 */
#if !defined(USENEWCODE)
	for (i=0; i<prm->Natom; i++)
          fscanf(file, " %d", &prm->TreeJoin[i]);
#else
        if (!read_fortran_12I6(file, prm->TreeJoin, prm->Natom)) {
          return MOLFILE_ERROR;
        }
#endif
        readtoeoln(file);

	/*
	 * READ PER-ATOM RES NUMBER			-IX(I66)
	 *	NOTE: this appears to be something entirely different
	 *	NOTE: overwriting this with correct PER-ATOM RES NUMBERs
	 */
#if !defined(USENEWCODE)
	for (i=0; i<prm->Natom; i++)
          fscanf(file, " %d", &prm->AtomRes[i]);
#else
        if (!read_fortran_12I6(file, prm->AtomRes, prm->Natom)) {
          return MOLFILE_ERROR;
        }
#endif
        res = 0;
	for (i=0; i<prm->Natom; i++) {
          if (i+1 == prm->Ipres[res+1])	/* atom is 1st of next res */
            res++;
          prm->AtomRes[i] = res;
	}
      
	/*
	 * BOUNDARY CONDITION STUFF
	 */
	if (!prm->IfBox) {
		prm->Nspm = 1;
		prm->Boundary = (int *) get(sizeof(int)*prm->Nspm);
		prm->Boundary[0] = prm->Natom;
	} else {
                readtoeoln(file);
#if !defined(USENEWCODE)
		fscanf(file, " %d %d %d", &prm->Iptres, &prm->Nspm, 
								&prm->Nspsol);
#else
                buffer = new int[3];
                if (!read_fortran_12I6(file, buffer, 3)) { 
                  return MOLFILE_ERROR;
                }
                prm->Iptres = buffer[0];
                prm->Nspm = buffer[1];
                prm->Nspsol = buffer[2];
                delete [] buffer;
#endif
                readtoeoln(file);
		prm->Boundary = (int *) get(sizeof(int)*prm->Nspm);
#if !defined(USENEWCODE)
		for (i=0; i<prm->Nspm; i++)
			fscanf(file, " %d", &prm->Boundary[i]);
#else
                if (!read_fortran_12I6(file, prm->Boundary, prm->Nspm)) {
                  return MOLFILE_ERROR;
                }
#endif
                readtoeoln(file);
		fscanf(file, " " DBLFMT " " DBLFMT " " DBLFMT, 
				&prm->Box[0], &prm->Box[1], &prm->Box[2]);
                readtoeoln(file);
		if (prm->Iptres)
			prm->Ipatm = prm->Ipres[prm->Iptres] - 1; 
      		/* IF(IPTRES.GT.0) IPTATM = IX(I02+IPTRES-1+1)-1 */
	}

	/*
	 * ----- LOAD THE CAP INFORMATION IF NEEDED -----
	 */
	if (prm->IfCap) {
		/* if (prm->IfBox) 
			skipeoln(file); */
		fscanf(file, " %d " DBLFMT " " DBLFMT " " DBLFMT " " DBLFMT,
				&prm->Natcap, &prm->Cutcap, 
				&prm->Xcap, &prm->Ycap, &prm->Zcap);
	}

  return MOLFILE_SUCCESS;
}


void ReadPARM::get_parm_atom(int i, char *name, char *atype, char *resname,
char *segname, int *resid, float *q, float *m) {
  int nres = prm->Nres;

  int j,k;
  int flag = 0;
  char *blank = (char *) " ";

  *q = prm->Charges[i];
  *m = prm->Masses[i];

  for (k = 0; k < 4; k++) {
    if (prm->AtomNames[i*4+k] == *blank)
       name[k] = '\0';
    else
       name[k] = prm->AtomNames[i*4+k];
  }
  name[k] = '\0';

  for (k = 0; k < 4; k++) {
    if ((prm->AtomSym[i*4+k]) == *blank)
       atype[k] = '\0';
    else
       atype[k] = prm->AtomSym[i*4+k];
  }
  atype[k] = '\0';

  for (j = 0; j < nres-1; j++)
    if (((i+1) >= prm->Ipres[j]) && ((i+1) < prm->Ipres[j+1])) {
      *resid = j;
      resname[0] = prm->ResNames[j*4];
      resname[1] = prm->ResNames[j*4+1];
      resname[2] = prm->ResNames[j*4+2];
      resname[3] = '\0';
      flag = 1;
    }
  if (flag == 0) {
     *resid = j;
     resname[0] = prm->ResNames[j*4];
     resname[1] = prm->ResNames[j*4+1];
     resname[2] = prm->ResNames[j*4+2];
     resname[3] = '\0';
     flag = 1;
  }

  segname[0] = '\0'; 

}

void ReadPARM::get_parm_bond (int i, int fromAtom[], int toAtom[]) {

  if (i < prm->Nbona) {
     fromAtom[i] = (int) ((prm->BondAt1[i])/3 + 1);
     toAtom[i] =   (int) ((prm->BondAt2[i])/3 + 1);
  }
  else get_hydrogen_bond (i, fromAtom, toAtom);
}

void ReadPARM::get_hydrogen_bond(int i, int fromAtom[], int toAtom[]) {
  fromAtom[i] = (int) ((prm->BondHAt1[i-prm->Nbona])/3 + 1);
  toAtom[i] = (int) ((prm->BondHAt2[i-prm->Nbona])/3 + 1);
}


int ReadPARM::get_parm_natoms() {return prm->Natom;}

int ReadPARM::get_parm_nbonds() {return (prm->Nbona + prm->Nbonh);}

int ReadPARM::get_parm_boxInfo() {return(prm->IfBox);}


// Read FORTRAN 12I6 format data (no space between adjacent data is assumed)
// One needs to read the whole data block into the buffer here
// fp - file pointer.
// data - buffer to hold the whole block. Should be allocated before the call
// count - number of ints in the block.
int ReadPARM::read_fortran_12I6(FILE *fp, int *data, int count) {
  int i, j;
  char buf[7];

  for (i=0; i<count; ++i) { 
    for (j=0; j<6; ++j) { 
      buf[j]=getc(fp);
      if (buf[j]=='\n' || buf[j]=='\0' || buf[j]==EOF)
        return 0;
    }
    buf[6] = '\0';

    if (sscanf(buf,"%d",data+i) != 1)
      return 0;

    if (i%12==11 && i<count-1)
      readtoeoln(fp);
  }

  return 1;
}


#endif

