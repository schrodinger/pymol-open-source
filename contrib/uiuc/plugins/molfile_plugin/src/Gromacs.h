/***************************************************************************
 *cr
 *cr            (C) Copyright 1995-2006 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/
/***************************************************************************
 * RCS INFORMATION:
 *      $RCSfile: Gromacs.h,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.22 $       $Date: 2006/02/24 00:57:10 $
 ***************************************************************************/

/*
 * GROMACS file format reader for VMD
 *
 * This code provides a high level I/O library for reading
 * and writing the following file formats:
 *	gro	GROMACS format or trajectory
 *	g96	GROMOS-96 format or trajectory
 *	trj	Trajectory - x, v and f (binary, full precision)
 *	trr	Trajectory - x, v and f (binary, full precision, portable)
 *	xtc	Trajectory - x only (compressed, portable, any precision)
 *      top
 * Currently supported: gro trj trr g96 [xtc]
 *
 * TODO list
 *   o  velocities are ignored because VMD doesn't use them, but some other 
 *      program might ...
 *   o  gro_rec() assumes positions in .gro files are nanometers and
 *      converts to angstroms, whereas they really could be any unit
 */

#ifndef GROMACS_H
#define GROMACS_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#if defined(_AIX)
#include <strings.h>
#endif

#include "endianswap.h"

#if defined(WIN32) || defined(WIN64)
#define strcasecmp stricmp
#endif

#ifndef M_PI_2
#define M_PI_2 1.57079632679489661922
#endif

// Error codes for mdio_errno
#define MDIO_SUCCESS		0
#define MDIO_BADFORMAT		1
#define MDIO_EOF		2
#define MDIO_BADPARAMS		3
#define MDIO_IOERROR		4
#define MDIO_BADPRECISION	5
#define MDIO_BADMALLOC		6
#define MDIO_CANTOPEN		7
#define MDIO_BADEXTENSION	8
#define MDIO_UNKNOWNFMT		9
#define MDIO_CANTCLOSE		10
#define MDIO_WRONGFORMAT	11
#define MDIO_SIZEERROR		12
#define MDIO_UNKNOWNERROR	1000

#define MDIO_READ	0
#define MDIO_WRITE	1

#define MDIO_MAX_ERRVAL		11

// Format extensions
const char *mdio_fmtexts[] = {
    "",
    ".gro",
    ".trr",
    ".g96",
    ".trj",
    ".xtc",
    NULL
};


static int mdio_errcode;	// Last error code

#define TRX_MAGIC	1993	// Magic number for .trX files
#define XTC_MAGIC	1995	// Magic number for .xtc files
#define MAX_GRO_LINE	500	// Maximum line length of .gro files
#define MAX_G96_LINE	500	// Maximum line length of .g96 files
#define MAX_TRX_TITLE	80	// Maximum length of a title in .trX
#define MAX_MDIO_TITLE	80	// Maximum supported title length
#define ANGS_PER_NM	10	// Unit conversion factor


// All the supported file types and their respective extensions
#define MDFMT_GRO		1
#define MDFMT_TRR		2
#define MDFMT_G96		3
#define MDFMT_TRJ		4
#define MDFMT_XTC		5


// A structure to hold .trX file format header information. This
// is an optional member of the md_file structure that is used
// when .trX files are being dealt with.
typedef struct {
	int version;		// File version number
	char title[MAX_TRX_TITLE + 1];	// File title
	int ir_size;
	int e_size;
	int box_size;
	int vir_size;
	int pres_size;
	int top_size;
	int sym_size;
	int x_size;		// Positions of atoms
	int v_size;		// Velocities of atoms
	int f_size;
	int natoms;		// Number of atoms in the system
	int step;
	int nre;
	float t;
	float lambda;
} trx_hdr;


// A generic i/o structure that contains information about the
// file itself and the input/output state
typedef struct {
	FILE *	f;	// Pointer to the file
	int	fmt;	// The file format
	int	prec;	// Real number precision
	int	rev;	// Reverse endiannism?
	trx_hdr * trx;	// Trx files require a great deal more
			// header data to be stored.
} md_file;


// A format-independent structure to hold header data from files
typedef struct {
	char title[MAX_MDIO_TITLE + 1];
	int natoms;
	float timeval;
} md_header;


// A format-independent structure to hold unit cell data
typedef struct {
  float A, B, C, alpha, beta, gamma;
} md_box;


// Timestep information
typedef struct {
	float *pos;	// Position array (3 * natoms)
	//float *vel;	// Velocity array ** (VMD doesn't use this) **
	//float *f;	// Force array ** (VMD doesn't use this) **
	//float *box;	// Computational box ** (VMD doesn't use this) **
	int natoms;	// Number of atoms
	int step;	// Simulation step
	float time;	// Time of simulation
  md_box *box;
} md_ts;


// Atom information
typedef struct {
	char resid[7];		// Residue index number
	char resname[7];	// Residue name
	int atomnum;		// Atom index number
	char atomname[7];	// Atom name
	float pos[3];		// Position array (3 * natoms)
	//float vel[3];	// Velocity array ** (VMD doesn't use this) **
} md_atom;


// Open a molecular dynamics file. The second parameter specifies
// the format of the file. If it is zero, the format is determined
// from the file extension. the third argument (if given) decides
// whether to read (==0) or to write (!= 0).
// using a default argument set to read for backward compatibility.
static md_file *mdio_open(const char *, const int, const int=MDIO_READ);

// Closes a molecular dynamics file.
static int mdio_close(md_file *);


// Format-independent file I/O routines
static int mdio_header(md_file *, md_header *);
static int mdio_timestep(md_file *, md_ts *);


// .gro file functions
static int gro_header(md_file *, char *, int, float *, int *, int = 1);
static int gro_rec(md_file *, md_atom *);
static int gro_timestep(md_file *, md_ts *);


// .trX file functions
static int trx_header(md_file *, int = 0);
static int trx_int(md_file *, int *);
static int trx_real(md_file *, float *);

static int trx_rvector(md_file *, float *);
static int trx_string(md_file *, char *, int);
static int trx_timestep(md_file *, md_ts *);

// .g96 file functions
static int g96_header(md_file *, char *, int, float *);
static int g96_timestep(md_file *, md_ts *);
static int g96_rec(md_file *, md_atom *);
static int g96_countatoms(md_file *);


// .xtc file functions
static int xtc_int(md_file *, int *);
static int xtc_float(md_file *, float *);
/* 
static int xtc_receivebits(int *, int);
static void xtc_receiveints(int *, int, int, const unsigned *, int *);
*/
static int xtc_timestep(md_file *, md_ts *);
static int xtc_3dfcoord(md_file *, float *, int *, float *);


// Error reporting functions
static int mdio_errno(void);
static const char *mdio_errmsg(int);
static int mdio_seterror(int);


// Miscellaneous functions
static int strip_white(char *);
static int mdio_readline(md_file *, char *, int, int = 1);
static int mdio_tsfree(md_ts *, int = 0);
static int mdio_readbox(md_box *, float *, float *, float *);



static int xtc_receivebits(int *, int);

// Error descriptions for mdio_errno
static const char *mdio_errdescs[] = {
	"no error",
	"file does not match format",
	"unexpected end-of-file reached",
	"function called with bad parameters",
	"file i/o error",
	"unsupported precision",
	"out of memory",
	"cannot open file",
	"bad file extension",
	"unknown file format",
	"cannot close file",
	"wrong file format for this function",
	"binary i/o error: sizeof(int) != 4",
	NULL
};

/*! \fn static inline bool host_is_little_endian(void)
 * detect endiannes of host machine. returns true on little endian machines. */
static inline int host_is_little_endian(void) 
{
  const union { char c[4]; unsigned int i; } 
  fixed = { { 0x10 , 0x20 , 0x40 , 0x80 } };
  const unsigned int i = 0x80402010U;
        
  if (fixed.i == i) {
    return 1;
  }
  return 0;
}



// Open a molecular dynamics file. The second parameter specifies
// the format of the file. If it is zero, the format is determined
// from the file extension.
md_file *mdio_open(const char *fn, const int fmt, const int rw) {
	md_file *mf;

	if (!fn) {
		mdio_seterror(MDIO_BADPARAMS);
		return NULL;
	}

	// Allocate memory
	mf = (md_file *) malloc(sizeof(md_file));
	if (!mf) {
		mdio_seterror(MDIO_BADMALLOC);
		return NULL;
	}

	// Zero out the structure
	memset(mf, 0, sizeof(md_file));

	// Determine the file type from the extension
	if (!fmt) {
		char *p;
		int n;

		// Seek to the extension part of the filename
		for (p = (char *) &fn[strlen(fn) - 1]; *p != '.' && p > fn; p--);
		if (p == fn) {
			free(mf);
			mdio_seterror(MDIO_BADEXTENSION);
			return NULL;
		}

		// Check the extension against known extensions
		for (n = 1; mdio_fmtexts[n]; n++)
			if (!strcasecmp(p, mdio_fmtexts[n])) break;

		// If !mdio_fmtexts[n], we failed (unknown ext)
		if (!mdio_fmtexts[n]) {
			free(mf);
			mdio_seterror(MDIO_UNKNOWNFMT);
			return NULL;
		}

		// All set
		mf->fmt = n;
	}
	else {
		mf->fmt = fmt;
	}

	// Differentiate between binary and ascii files. Also,
	// .trX files need a header information structure allocated.
	switch (mf->fmt) {
    case MDFMT_GRO:
	case MDFMT_G96: /* fallthrough */
        if (rw) 
            mf->f = fopen(fn, "wt");
        else
            mf->f = fopen(fn, "rt");

		break;
	case MDFMT_TRR:
	case MDFMT_TRJ: /* fallthrough */
		// Allocate the trx header data struct
		mf->trx = (trx_hdr *) malloc(sizeof(trx_hdr));
		if (!mf->trx) {
			free(mf);
			mdio_seterror(MDIO_BADMALLOC);
			return NULL;
		}
		memset(mf->trx, 0, sizeof(trx_hdr));
	case MDFMT_XTC:  /* fallthrough */
		// Finally, open the file
        if (rw)
            mf->f = fopen(fn, "wb");
        else
            mf->f = fopen(fn, "rb");

		break;
	default:
		free(mf);
		mdio_seterror(MDIO_UNKNOWNFMT);
		return NULL;
	}

	// Check for opening error
	if (!mf->f) {
		if (mf->trx) free(mf->trx);
		free(mf);
		mdio_seterror(MDIO_CANTOPEN);
		return NULL;
	}

	// File is opened, we're all set!
	mdio_seterror(MDIO_SUCCESS);
	return mf;
}


// Closes a molecular dynamics file.
static int mdio_close(md_file *mf) {
	if (!mf) return mdio_seterror(MDIO_BADPARAMS);

	if (fclose(mf->f) == EOF) return mdio_seterror(MDIO_CANTCLOSE);

	// Free the dynamically allocated memory
	if (mf->trx) free(mf->trx);
	free(mf);
	return mdio_seterror(MDIO_SUCCESS);
}


// Returns the last error code reported by any of the mdio functions
static int mdio_errno(void) {
	return mdio_errcode;
}


// Returns a textual message regarding an mdio error code
static const char *mdio_errmsg(int n) {
	if (n < 0 || n > MDIO_MAX_ERRVAL) return (char *) "unknown error";
	else return mdio_errdescs[n];
}


// Sets the error code and returns an appropriate return value
// for the calling function to return to its parent
static int mdio_seterror(int code) {
	mdio_errcode = code;
	return code ? -1 : 0;
}


// Reads a line from the text file, strips leading/trailing whitespace
// and newline, checks for errors, and returns the number of characters
// in the string on success or -1 on error.
static int mdio_readline(md_file *mf, char *buf, int n, int strip) {
	if (!buf || n < 1 || !mf) return mdio_seterror(MDIO_BADPARAMS);

	// Read the line
	fgets(buf, n, mf->f);

	// End of file reached?
	if (feof(mf->f)) return mdio_seterror(MDIO_EOF);

	// File I/O error?
	if (ferror(mf->f)) return mdio_seterror(MDIO_IOERROR);

	// Strip whitespace
	if (strip) strip_white(buf);

	return strlen(buf);
}


// Strips leading and trailing whitespace from a string. Tabs,
// spaces, newlines and carriage returns are stripped. Example:
// "\n   hello\t \r" becomes "hello".
static int strip_white(char *buf) {
	int i, j, k;

	// Protect against NULL pointer
	if (!buf) return -1;
	if (!strlen(buf)) return -1;

	// Kill trailing whitespace first
	for (i = strlen(buf) - 1;
	     buf[i] == ' ' || buf[i] == '\t' ||
	     buf[i] == '\n' || buf[i] == '\r';
	     i--)
		buf[i] = 0;

	// Skip past leading whitespace
	for (i = 0; buf[i] == ' ' || buf[i] == '\t' ||
	     buf[i] == '\n' || buf[i] == '\r'; i++);
	if (i) {
		k = 0;
		for (j = i; buf[j]; j++)
			buf[k++] = buf[j];
		buf[k] = 0;
	}

	return strlen(buf);
}


// Frees the memory allocated in a ts structure. The holderror
// parameter defaults to zero. Programs that are calling this
// function because of an error reported by another function should
// set holderror so that mdio_tsfree() does not overwrite the error
// code with mdio_seterror().
static int mdio_tsfree(md_ts *ts, int holderror) {
	if (!ts) {
		if (holderror) return -1;
		else return mdio_seterror(MDIO_BADPARAMS);
	}

	if (ts->pos && ts->natoms > 0) free(ts->pos);

  if (ts->box) free(ts->box);

	if (holderror) return -1;
	else return mdio_seterror(MDIO_SUCCESS);
}


// Converts box basis vectors to A, B, C, alpha, beta, and gamma.  
// Stores values in md_box struct, which should be allocated before calling
// this function.
static int mdio_readbox(md_box *box, float *x, float *y, float *z) {
  float A, B, C;

  if (!box) {
    return mdio_seterror(MDIO_BADPARAMS);
  }

  // A, B, C are the lengths of the x, y, z vectors, respectively
  A = sqrt( x[0]*x[0] + x[1]*x[1] + x[2]*x[2] ) * ANGS_PER_NM;
  B = sqrt( y[0]*y[0] + y[1]*y[1] + y[2]*y[2] ) * ANGS_PER_NM;
  C = sqrt( z[0]*z[0] + z[1]*z[1] + z[2]*z[2] ) * ANGS_PER_NM;
  if ((A<=0) || (B<=0) || (C<=0)) {
    /* Use zero-length box size and set angles to 90. */
    box->A = box->B = box->C = 0;
    box->alpha = box->beta = box->gamma = 90;
  } else {
    box->A = A;
    box->B = B;
    box->C = C;
  
    // gamma, beta, alpha are the angles between the x & y, x & z, y & z
    // vectors, respectively
    box->gamma = acos( (x[0]*y[0]+x[1]*y[1]+x[2]*y[2])/(A*B) ) * 90.0/M_PI_2;
    box->beta = acos( (x[0]*z[0]+x[1]*z[1]+x[2]*z[2])/(A*C) ) * 90.0/M_PI_2;
    box->alpha = acos( (y[0]*z[0]+y[1]*z[1]+y[2]*z[2])/(B*C) ) * 90.0/M_PI_2; 
  }
  return mdio_seterror(MDIO_SUCCESS);
}


// Reads the header of a file (format independent)
static int mdio_header(md_file *mf, md_header *mdh) {
	int n;
	if (!mf || !mdh) return mdio_seterror(MDIO_BADPARAMS);
	if (!mf->f) return mdio_seterror(MDIO_BADPARAMS);

	switch (mf->fmt) {
	case MDFMT_GRO:
		if (gro_header(mf, mdh->title, MAX_MDIO_TITLE,
		&mdh->timeval, &mdh->natoms, 1) < 0)
			return -1;
		return 0;

	case MDFMT_TRR: 
	case MDFMT_TRJ: /* fallthrough */
		if (trx_header(mf, 1) < 0) return -1;
		mdh->natoms = mf->trx->natoms;
		mdh->timeval = (float) mf->trx->t;
		strncpy(mdh->title, mf->trx->title, MAX_MDIO_TITLE);
		return 0;

	case MDFMT_G96:
		if (g96_header(mf, mdh->title, MAX_MDIO_TITLE,
		&mdh->timeval) < 0) return -1;
		mdh->natoms = -1;
		return 0;

	case MDFMT_XTC:
		memset(mdh, 0, sizeof(md_header));
		// Check magic number
		if (xtc_int(mf, &n) < 0) return -1;
		if (n != XTC_MAGIC) return mdio_seterror(MDIO_BADFORMAT);

		// Get number of atoms
		if (xtc_int(mf, &n) < 0) return -1;
		mdh->natoms = n;
		rewind(mf->f);
		return 0;

	default:
		return mdio_seterror(MDIO_UNKNOWNFMT);
	}
}


// Reads in a timestep from a file (format independent)
static int mdio_timestep(md_file *mf, md_ts *ts) {
	if (!mf || !ts) return mdio_seterror(MDIO_BADPARAMS);
	if (!mf->f) return mdio_seterror(MDIO_BADPARAMS);

	switch (mf->fmt) {
	case MDFMT_GRO:
		return gro_timestep(mf, ts);

	case MDFMT_TRR:
	case MDFMT_TRJ: /* fallthrough */
		return trx_timestep(mf, ts);

	case MDFMT_G96:
		return g96_timestep(mf, ts);

	case MDFMT_XTC:
		return xtc_timestep(mf, ts);

	default:
		return mdio_seterror(MDIO_UNKNOWNFMT);
	}
}



static int g96_header(md_file *mf, char *title, int titlelen, float *timeval) {
	char buf[MAX_G96_LINE + 1];
	char *p;

	// Check parameters
	if (!mf) return mdio_seterror(MDIO_BADPARAMS);

	// The header consists of blocks. The title block
	// is mandatory, and a TIMESTEP block is optional.
	// Example:
	//
	// TITLE
	// Generated by trjconv :  t=  90.00000
	// more title info
	// .
	// .
	// .
	// END
	// .
	// .
	// .

	if (mdio_readline(mf, buf, MAX_G96_LINE + 1) < 0) return -1;
	if (strcasecmp(buf, "TITLE")) return mdio_seterror(MDIO_BADFORMAT);

	// Read in the title itself
	if (mdio_readline(mf, buf, MAX_G96_LINE + 1) < 0) return -1;

        // The timevalue can be included in the title string
        // after a "t=" prefix.
        if ((p = (char *) strstr(buf, "t="))) {
                char *q = p;
                *(q--) = 0;

                // Skip the `t=' and strip whitespace from
                // the resulting strings
                p += 2;
                strip_white(p);
                strip_white(buf);

                // Grab the timevalue from the title string
                if (timeval) *timeval = (float) atof(p);
        }
        else {
                // No timevalue - just copy the string and strip
                // any leading/trailing whitespace
                if (timeval) *timeval = 0;
                strip_white(buf);
        }

	// Copy the title string
	if (title && titlelen) strncpy(title, buf, titlelen);

	// Now ignore subsequent title lines and get the END string
	while (strcasecmp(buf, "END"))
		if (mdio_readline(mf, buf, MAX_G96_LINE + 1) < 0) return -1;

	// Done!
	return mdio_seterror(MDIO_SUCCESS);
}


// Used to determine the number of atoms in a g96 file, because
// VMD needs to know this for some reason.
static int g96_countatoms(md_file *mf) {
	char buf[MAX_G96_LINE + 1];
	int natoms;
	int n;
	long fpos;
	float lastf;

	if (!mf) return mdio_seterror(MDIO_BADPARAMS);

	fpos = ftell(mf->f);

	natoms = 0;
	for (;;) {
		if (mdio_readline(mf, buf, MAX_G96_LINE + 1, 0) < 0)
			break;
		n = sscanf(buf, "%*6c%*6c%*6c%*6c %*f %*f %f", &lastf);
		if (n == 1) natoms++;
		else {
			strip_white(buf);
			if (!strcasecmp(buf, "END")) break;
		}
	}

	fseek(mf->f, fpos, SEEK_SET);
	return natoms;
}


// Reads an atom line from the G96 file
static int g96_rec(md_file *mf, md_atom *ma) {
	char buf[MAX_G96_LINE + 1];
	char atomnum[7];
	int n;

	// Check parameters
	if (!mf || !ma) return mdio_seterror(MDIO_BADPARAMS);

	// Read in a line, assuming it is an atom line
	do {
		if (mdio_readline(mf, buf, MAX_G96_LINE + 1, 0) < 0) return -1;
	} while (buf[0] == '#' || strlen(buf) == 0);

	n = sscanf(buf, "%6c%6c%6c%6c %f %f %f",
		ma->resid, ma->resname, ma->atomname, atomnum,
		&ma->pos[0], &ma->pos[1], &ma->pos[2]);
	if (n == 7) {
		atomnum[6] = 0;
		ma->resid[6] = 0;
		ma->resname[6] = 0;
		ma->atomname[6] = 0;

		strip_white(atomnum);
		strip_white(ma->resid);
		strip_white(ma->resname);
		strip_white(ma->atomname);

		ma->atomnum = atoi(atomnum);

		ma->pos[0] *= ANGS_PER_NM;
		ma->pos[1] *= ANGS_PER_NM;
		ma->pos[2] *= ANGS_PER_NM;

		return 0;
	}

	return mdio_seterror(MDIO_BADFORMAT);
}


// Reads a timestep from a G96 file and stores the data in
// the generic md_ts structure. Returns 0 on success or a
// negative number on error and sets mdio_errcode.
static int g96_timestep(md_file *mf, md_ts *ts) {
	char		buf[MAX_G96_LINE + 1];
	char		stripbuf[MAX_G96_LINE + 1];
	float		pos[3], x[3], y[3], z[3], *currAtom;
	long		fpos;
	int		n, i, boxItems;

	// Check parameters
	if (!mf || !ts) return mdio_seterror(MDIO_BADPARAMS);

  // Allocate data space for the timestep, using the number of atoms
  // determined by open_g96_read().
	ts->pos = (float *) malloc(sizeof(float) * 3 * ts->natoms);
	if (!ts->pos) {
		return mdio_seterror(MDIO_BADMALLOC);
	}
  currAtom = ts->pos;

	// The timesteps follow the header in a fixed block
	// format:
	//
	// TIMESTEP
	//         <step number> <time value>
	// END
	// POSITIONRED
	//     <x float> <y float> <z float>
	//     .         .         .
	//     .         .         .
	//     .         .         .
	// END
	// VELOCITYRED
	//     <x float> <y float> <z float>
	//     .         .         .
	//     .         .         .
	//     .         .         .
	// END
	// BOX
	//     <x float> <y float> <z float>
	// END
	//
	// -----
	//
	// The TIMESTEP, VELOCITY and BOX blocks are optional.
	// Floats are written in 15.9 precision.
	//
	// Reference: GROMACS 2.0 user manual
	//            http://rugmd4.chem.rug.nl/~gmx/online2.0/g96.html

	// First, look for an (optional) title block and skip it
	if (mdio_readline(mf, buf, MAX_G96_LINE + 1) < 0) return -1;

  if (!strcasecmp(buf, "TITLE")) {
    // skip over the text until we reach 'END'
    while (strcasecmp(buf, "END")) {
      if (mdio_readline(mf, buf, MAX_G96_LINE + 1) < 0) return -1;
    }

    // Read in the next line
    if (mdio_readline(mf, buf, MAX_G96_LINE + 1) < 0) return -1;
  }

	// Next, look for a timestep block
	if (!strcasecmp(buf, "TIMESTEP")) {
		// Read in the value line
		if (mdio_readline(mf, buf, MAX_G96_LINE + 1) < 0) return -1;

		// Extract the time value and the timestep index
		n = sscanf(buf, "%d %f", &ts->step, &ts->time);
		if (n != 2) return mdio_seterror(MDIO_BADFORMAT);

		// Read the "END" line
		if (mdio_readline(mf, buf, MAX_G96_LINE + 1) < 0) return -1;
		if (strcasecmp(buf, "END"))
			return mdio_seterror(MDIO_BADFORMAT);

		// Read in the next line
		if (mdio_readline(mf, buf, MAX_G96_LINE + 1) < 0) return -1;
	}
	else {
		// No timestep specified -- set to zero
		ts->step = 0;
		ts->time = 0;
	}

	// At this point a POSITION or POSITIONRED block
	// is REQUIRED by the format
	if (!strcasecmp(buf, "POSITIONRED")) {

    // So now we read in some atoms
    i = 0;
		while (i < ts->natoms) {
			// Read in an atom
			if (mdio_readline(mf, buf, MAX_G96_LINE + 1) < 0)
				return -1;
 
      // We shouldn't reach the end yet
      if (!strcasecmp(buf, "END"))
        return mdio_seterror(MDIO_BADFORMAT);

			// Get the x,y,z coordinates
			n = sscanf(buf, "%f %f %f", &pos[0], &pos[1], &pos[2]);
      
      // Ignore improperly formatted lines
			if (n == 3) {
				pos[0] *= ANGS_PER_NM;
				pos[1] *= ANGS_PER_NM;
				pos[2] *= ANGS_PER_NM;

				// Copy the atom data into the array
				memcpy(currAtom, pos, sizeof(float) * 3);
        currAtom += 3;
        i++;
			}
		}
	}
	else if (!strcasecmp(buf, "POSITION") || !strcasecmp(buf, "REFPOSITION")) {
		/*
		char resnum[7];
		char resname[7];
		char atomname[7];
		char atomnum[7];
		*/

		// So now we read in some atoms
    i = 0;
		while (i < ts->natoms) {
			// Read in the first line
			if (mdio_readline(mf, buf, MAX_G96_LINE + 1, 0) < 0)
				return -1;
 
      // We shouldn't reach the end yet
      strcpy(stripbuf, buf);
      strip_white(stripbuf); 
      if (!strcasecmp(stripbuf, "END"))
        return mdio_seterror(MDIO_BADFORMAT);

			// Get the x,y,z coordinates and name data
			n = sscanf(buf, "%*6c%*6c%*6c%*6c %f %f %f",
				&pos[0], &pos[1], &pos[2]);

      // Ignore improperly formatted lines
			if (n == 3) {
				pos[0] *= ANGS_PER_NM;
				pos[1] *= ANGS_PER_NM;
				pos[2] *= ANGS_PER_NM;

				// Copy the atom data into the linked list item
				memcpy(currAtom, pos, sizeof(float) * 3);
				currAtom += 3;
        i++;
			}
		}
	}
	else {
		return mdio_seterror(MDIO_BADFORMAT);
	}

  // Read the END keyword
  if (mdio_readline(mf, buf, MAX_G96_LINE + 1) < 0)
    return -1;
  if (strcasecmp(buf, "END"))
    return mdio_seterror(MDIO_BADFORMAT);

	// ... another problem: there may or may not be a VELOCITY
	// block or a BOX block, so we need to read one line beyond
	// the POSITION block to determine this. If neither VEL. nor
	// BOX are present we've read a line too far and infringed
	// on the next timestep, so we need to keep track of the
	// position now for a possible fseek() later to backtrack.
	fpos = ftell(mf->f);

	// Now we must read in the velocities and the box, if present
	if (mdio_readline(mf, buf, MAX_G96_LINE + 1) < 0) {
    // It's okay if we end the file here; any other errors need to be
    // reported.
    if (mdio_errcode == MDIO_EOF) 
      return mdio_seterror(MDIO_SUCCESS);
    else 
      return -1;
  }

	// Is there a velocity block present ?
	if (!strcasecmp(buf, "VELOCITY") || !strcasecmp(buf, "VELOCITYRED")) {
		// Ignore all the coordinates - VMD doesn't use them
		for (;;) {
			if (mdio_readline(mf, buf, MAX_G96_LINE + 1) < 0)
				return -1;
			if (!strcasecmp(buf, "END")) break;
		}

		// Again, record our position because we may need
		// to fseek here later if we read too far.
		fpos = ftell(mf->f);

		// Go ahead and read the next line.
		if (mdio_readline(mf, buf, MAX_G96_LINE + 1) < 0) return -1;
	}

	// Is there a box present ?
	if (!strcasecmp(buf, "BOX")) {
		if (mdio_readline(mf, buf, MAX_G96_LINE + 1) < 0) return -1;
    boxItems = sscanf(buf, " %f %f %f %f %f %f %f %f %f", 
               &x[0], &y[1], &z[2], &x[1], &x[2], &y[0], &y[2], &z[0], &z[1]);
    if (boxItems == 3) {
      x[1] = x[2] = 0;
      y[0] = y[2] = 0;
      z[0] = z[1] = 0;
    }
    else if (boxItems != 9) 
      return mdio_seterror(MDIO_BADFORMAT);

    // Allocate the box and convert the vectors.
    ts->box = (md_box *) malloc(sizeof(md_box));
    if (mdio_readbox(ts->box, x, y, z) < 0) {
      free(ts->box);
      ts->box = NULL;
      return mdio_seterror(MDIO_BADFORMAT);
    }

		if (mdio_readline(mf, buf, MAX_G96_LINE + 1) < 0) {
      free(ts->box);
      ts->box = NULL;
      return -1;
    }
		if (strcasecmp(buf, "END")) {
      free(ts->box);
      ts->box = NULL;
			return mdio_seterror(MDIO_BADFORMAT);
    }
	}
	else {
		// We have read too far, so fseek back to the
		// last known safe position so we don't return
		// with the file pointer set infringing on the
		// next timestep data.
		fseek(mf->f, fpos, SEEK_SET);
	}

	// We're done!
	return mdio_seterror(MDIO_SUCCESS);
}


// Attempts to read header data from a GROMACS structure file
// The GROMACS header format is as follows (fixed, 2 lines ASCII):
// <title> [ n= <timevalue> ]
//     <num atoms>
static int gro_header(md_file *mf, char *title, int titlelen, float *timeval,
               int *natoms, int rewind) {
	char buf[MAX_GRO_LINE + 1];
	long fpos;
	char *p;

	// Check parameters
	if (!mf) return mdio_seterror(MDIO_BADPARAMS);

	// Get the current file position for rewinding later
	fpos = ftell(mf->f);

	// The header consists of 2 lines - get the first line
	if (mdio_readline(mf, buf, MAX_GRO_LINE + 1) < 0) return -1;

	// The timevalue can be included in the title string
	// after a "t=" prefix.
	if ((p = (char *) strstr(buf, "t="))) {
		char *q = p;
		*(q--) = 0;

		// Skip the `t=' and strip whitespace from
		// the resulting strings
		p += 2;
		strip_white(p);
		strip_white(buf);

		// Grab the timevalue from the title string
		if (timeval) *timeval = (float) atof(p);
	}
	else {
		// No timevalue - just copy the string
		if (timeval) *timeval = 0;
	}

	// Copy the title string
	if (title && titlelen) strncpy(title, buf, titlelen);

	// Get the second line and grab the number of atoms
	if (mdio_readline(mf, buf, MAX_GRO_LINE + 1) < 0) return -1;

	// Store the number of atoms
	if (natoms) if (!(*natoms = atoi(buf)))
		return mdio_seterror(MDIO_BADFORMAT);

	// Now we rewind the file so that subsequent calls to
	// gro_timestep() will succeed. gro_timestep() requires
	// the header to be at the current file pointer.
	if (rewind) fseek(mf->f, fpos, SEEK_SET);

	// Done!
	return 0;
}


// Reads one atom record from a GROMACS file. Returns GMX_SUCCESS
// on success or a negative number on error.
//
// Record format (one line, fixed):
//    rrrrrRRRRRaaaaaAAAAA <x pos> <y pos> <z pos> <x vel> <y vel> <z vel>
//
//    r = residue number
//    R = residue name
//    a = atom name
//    A = atom number
//
static int gro_rec(md_file *mf, md_atom *ma) {
	char	buf[MAX_GRO_LINE + 1];
	char	atomnum[6];
	int	n;

	if (!mf) return mdio_seterror(MDIO_BADPARAMS);

	do {
		if (mdio_readline(mf, buf, MAX_GRO_LINE + 1, 0) < 0) return -1;
	} while (buf[0] == '#' || !strlen(buf));

	// Read in the fields
	n = sscanf(buf, "%5c%5c%5c%5c%f %f %f", ma->resid,
		ma->resname, ma->atomname, atomnum, ma->pos,
		&ma->pos[1], &ma->pos[2]);
	if (n != 7) return mdio_seterror(MDIO_BADFORMAT);

	// Null terminate the strings
	ma->resname[5] = 0;
	ma->atomname[5] = 0;
	ma->resid[5] = 0;
	atomnum[5] = 0;

	// Convert strings to numbers
	strip_white(atomnum);
	ma->atomnum = atoi(atomnum);

	// Convert nanometers to angstroms
	ma->pos[0] *= ANGS_PER_NM;
	ma->pos[1] *= ANGS_PER_NM;
	ma->pos[2] *= ANGS_PER_NM;

	// Strip leading and trailing whitespace
	strip_white(ma->atomname);
	strip_white(ma->resname);
	strip_white(ma->resid);

	return 0;
}


// Reads in a timestep from a .gro file. Ignores the data
// not needed for a timestep, so is a little faster than
// calling gro_rec() for each atom. Also reads in the
// header block.
//
static int gro_timestep(md_file *mf, md_ts *ts) {
	char buf[MAX_GRO_LINE + 1];
	long coord;
	int i, n, boxItems;
  float x[3], y[3], z[3];

	if (!mf || !ts) return mdio_seterror(MDIO_BADPARAMS);

	if (gro_header(mf, NULL, 0, &ts->time, &ts->natoms, 0) < 0)
		return -1;
	ts->pos = (float *) malloc(3 * sizeof(float) * ts->natoms);
	if (!ts->pos)
		return mdio_seterror(MDIO_BADMALLOC);

	coord = 0;
	for (i = 0; i < ts->natoms; i++) {
		if (mdio_readline(mf, buf, MAX_GRO_LINE + 1, 0) < 0) {
			free(ts->pos);
			return -1;
		}
	
		n = sscanf(buf, "%*5c%*5c%*5c%*5c%f %f %f",
			&ts->pos[coord], &ts->pos[coord + 1],
			&ts->pos[coord + 2]);

		ts->pos[coord] *= ANGS_PER_NM;
		ts->pos[coord + 1] *= ANGS_PER_NM;
		ts->pos[coord + 2] *= ANGS_PER_NM;

		if (n != 3) return mdio_seterror(MDIO_BADFORMAT);
		coord += 3;
	}

	// Read the box, stored as three vectors representing its edges
	if (mdio_readline(mf, buf, MAX_GRO_LINE + 1, 0) < 0) {
		free(ts->pos);
		return -1;
	}
  boxItems = sscanf(buf, " %f %f %f %f %f %f %f %f %f", 
             &x[0], &y[1], &z[2], &x[1], &x[2], &y[0], &y[2], &z[0], &z[1]);
  // File may only include three scalars for the box information -- if
  // that's the case, the box is orthoganal.
  if (boxItems == 3) {
    x[1] = x[2] = 0;
    y[0] = y[2] = 0;
    z[0] = z[1] = 0;
  }
  else if (boxItems != 9) {
    free(ts->pos);
    return -1;
  }

  // Allocate the box and convert the vectors.
  ts->box = (md_box *) malloc(sizeof(md_box));
  if (mdio_readbox(ts->box, x, y, z) < 0) {
    free(ts->pos);
    free(ts->box);
    ts->box = NULL;
    return -1;
  }

	return 0;
}


// Attempts to read header data from a .trX trajectory file
//
// The .trX header format is as follows:
//
//	4 bytes		- magic number (0x07C9)
//	...
//
static int trx_header(md_file *mf, int rewind) {
	int magic;
	trx_hdr *hdr;
	long fpos;

	if (!mf) return mdio_seterror(MDIO_BADPARAMS);

	// In case we need to rewind
	fpos = ftell(mf->f);

	// We need to store some data to the trX header data
	// structure inside the md_file structure
	hdr = mf->trx;
	if (!mf->trx) return mdio_seterror(MDIO_BADPARAMS);

	// Read the magic number
	if (trx_int(mf, &magic) < 0) return -1;
	if (magic != TRX_MAGIC) {
		// Try reverse endianism
		swap4_aligned(&magic, 1);
		if (magic != TRX_MAGIC) return mdio_seterror(MDIO_BADFORMAT);

		// Enable byte swapping (actually works, too!)
		mf->rev = 1;
	}

	// Read the version number. 
        // XXX. this is not the version number, but the storage size
	// of the following XDR encoded string.
	// the 'title' string is in fact the version identifier.
	// since VMD does not use any of that, it does no harm,
	// but is should still be fixed occasionally. AK 2005/01/08.

	if(mf->fmt!=MDFMT_TRJ) {
		// It appears that TRJ files either don't contain a version
		// number or don't have a length-delimiter on the string,
		// whereas TRR files do contain both.  Thus, with TRJ, we just
		// assume that the version number is the string length and 
		// just hope for the best. -- WLD 2006/07/09
		if (trx_int(mf, &hdr->version) < 0) return -1;
	}

	// Read in the title string
	if (trx_string(mf, hdr->title, MAX_TRX_TITLE) < 0)
		return -1;

	// Read in some size data
	if (trx_int(mf, &hdr->ir_size) < 0) return -1;
	if (trx_int(mf, &hdr->e_size) < 0) return -1;
	if (trx_int(mf, &hdr->box_size) < 0) return -1;
	if (trx_int(mf, &hdr->vir_size) < 0) return -1;
	if (trx_int(mf, &hdr->pres_size) < 0) return -1;
	if (trx_int(mf, &hdr->top_size) < 0) return -1;
	if (trx_int(mf, &hdr->sym_size) < 0) return -1;
	if (trx_int(mf, &hdr->x_size) < 0) return -1;
	if (trx_int(mf, &hdr->v_size) < 0) return -1;
	if (trx_int(mf, &hdr->f_size) < 0) return -1;
	if (trx_int(mf, &hdr->natoms) < 0) return -1;
	if (trx_int(mf, &hdr->step) < 0) return -1;
	if (trx_int(mf, &hdr->nre) < 0) return -1;

	// Make sure there are atoms...
	if (!hdr->natoms) return mdio_seterror(MDIO_BADFORMAT);

	// Try to determine precision (float? double?)
	if (hdr->x_size) mf->prec = hdr->x_size / (hdr->natoms * 3);
	else if (hdr->v_size) mf->prec = hdr->v_size / (hdr->natoms * 3);
	else if (hdr->f_size) mf->prec = hdr->f_size / (hdr->natoms * 3);
	else return mdio_seterror(MDIO_BADPRECISION);

	if (mf->prec != sizeof(float) && mf->prec != sizeof(double)) {
		// We have no data types this size! The
		// file must've been generated on another
		// platform
		return mdio_seterror(MDIO_BADPRECISION);
	}

	// Read in timestep and lambda
	if (trx_real(mf, &hdr->t) < 0) return -1;
	if (trx_real(mf, &hdr->lambda) < 0) return -1;

	// Rewind if necessary
	if (rewind) fseek(mf->f, fpos, SEEK_SET);

	return 0;
}


// Reads in an integer and stores it in y. Returns GMX_SUCCESS
// on success or a negative number on error.
static int trx_int(md_file *mf, int *y) {
	if (!mf) return mdio_seterror(MDIO_BADPARAMS);

        // sanity check.
        if (sizeof(int) != 4) return mdio_seterror(MDIO_SIZEERROR);

	if (y) {
		if (fread(y, 4, 1, mf->f) != 1)
			return mdio_seterror(MDIO_IOERROR);
		if (mf->rev) swap4_aligned(y, 1);
	}
	else if (fseek(mf->f, 4, SEEK_CUR) != 0)
		return mdio_seterror(MDIO_IOERROR);

	return mdio_seterror(MDIO_SUCCESS);
}


// Reads in either a float or a double, depending on the
// precision, and stores that number in y. Returns
// GMX_SUCCESS on success or a negative number on error.
static int trx_real(md_file *mf, float *y) {
	double x;

	if (!mf) return mdio_seterror(MDIO_BADPARAMS);

	switch (mf->prec) {
		case sizeof(float):
			if (!y) {
				if (fseek(mf->f, mf->prec, SEEK_CUR) != 0)
					return mdio_seterror(MDIO_IOERROR);
			} else {
				if (fread(y, mf->prec, 1, mf->f) != 1)
					return mdio_seterror(MDIO_IOERROR);
				if (mf->rev) swap4_aligned(y, 1);
			}
			return mdio_seterror(MDIO_SUCCESS);

		case sizeof(double):
			if (!y) {
				if (fseek(mf->f, mf->prec, SEEK_CUR) != 0)
					return mdio_seterror(MDIO_IOERROR);
			} else {
				if (fread(&x, mf->prec, 1, mf->f) != 1)
					return mdio_seterror(MDIO_IOERROR);
				if (mf->rev) swap8_aligned(&x, 1);
				*y = (float) x;
			}
			return mdio_seterror(MDIO_SUCCESS);

		default:
			return mdio_seterror(MDIO_BADPRECISION);
	}

}


// Reads in a real-valued vector (taking precision into account).
// Stores the vector in vec, and returns GMX_SUCCESS on success
// or a negative number on error.
static int trx_rvector(md_file *mf, float *vec) {
	if (!mf) return mdio_seterror(MDIO_BADPARAMS);

	if (!vec) {
		if (trx_real(mf, NULL) < 0) return -1;
		if (trx_real(mf, NULL) < 0) return -1;
		if (trx_real(mf, NULL) < 0) return -1;
		return mdio_seterror(MDIO_SUCCESS);
	} else {
		if (trx_real(mf, &vec[0]) < 0) return -1;
		if (trx_real(mf, &vec[1]) < 0) return -1;
		if (trx_real(mf, &vec[2]) < 0) return -1;
		return mdio_seterror(MDIO_SUCCESS);
	}
}


// Reads in a string by first reading an integer containing the
// string's length, then reading in the string itself and storing
// it in str. If the length is greater than max, it is truncated
// and the rest of the string is skipped in the file. Returns the
// length of the string on success or a negative number on error.
static int trx_string(md_file *mf, char *str, int max) {
	int size;
  size_t ssize;

	if (!mf) return mdio_seterror(MDIO_BADPARAMS);

	if (trx_int(mf, &size) < 0) return -1;
  ssize = (size_t)size;

	if (str && size <= max) {
		if (fread(str, 1, size, mf->f) != ssize)
			return mdio_seterror(MDIO_IOERROR);
		str[size] = 0;
		return size;
	} else if (str) {
		if (fread(str, 1, max, mf->f) != ssize)
			return mdio_seterror(MDIO_IOERROR);
		if (fseek(mf->f, size - max, SEEK_CUR) != 0)
			return mdio_seterror(MDIO_IOERROR);
		str[max] = 0;
		return max;
	} else {
		if (fseek(mf->f, size, SEEK_CUR) != 0)
			return mdio_seterror(MDIO_IOERROR);
		return 0;
	}
}


// Reads in a timestep frame from the .trX file and returns the
// data in a timestep structure. Returns NULL on error.
static int trx_timestep(md_file *mf, md_ts *ts) {
	int i;
  float x[3], y[3], z[3];
	trx_hdr *hdr;

	if (!mf || !ts) return mdio_seterror(MDIO_BADPARAMS);
	if (mf->fmt != MDFMT_TRJ && mf->fmt != MDFMT_TRR)
		return mdio_seterror(MDIO_WRONGFORMAT);

	// Read the header
	if (trx_header(mf) < 0) return -1;

	// We need some data from the trX header
	hdr = mf->trx;
	if (!hdr) return mdio_seterror(MDIO_BADPARAMS);

	if (hdr->box_size) { // XXX need to check value of box_size!!
		if (trx_rvector(mf, x) < 0) return -1;
		if (trx_rvector(mf, y) < 0) return -1;
		if (trx_rvector(mf, z) < 0) return -1;
    // Allocate the box and convert the vectors.
    ts->box = (md_box *) malloc(sizeof(md_box));
    if (mdio_readbox(ts->box, x, y, z) < 0) {
      free(ts->box);
      ts->box = NULL;
      return -1;
    }
	}

	if (hdr->vir_size) {
		if (trx_rvector(mf, NULL) < 0) return -1;
		if (trx_rvector(mf, NULL) < 0) return -1;
		if (trx_rvector(mf, NULL) < 0) return -1;
	}

        if (hdr->pres_size) {
                if (trx_rvector(mf, NULL) < 0) return -1;
                if (trx_rvector(mf, NULL) < 0) return -1;
                if (trx_rvector(mf, NULL) < 0) return -1;
        }

        if (hdr->x_size) {
                ts->pos = (float *) malloc(sizeof(float) * 3 * hdr->natoms);
                if (!ts->pos) return mdio_seterror(MDIO_BADMALLOC);

		ts->natoms = hdr->natoms;

		for (i = 0; i < hdr->natoms; i++) {
			if (trx_rvector(mf, &ts->pos[i * 3]) < 0) {
				mdio_tsfree(ts, 1);
				return -1;
			}
			ts->pos[i * 3] *= ANGS_PER_NM;
			ts->pos[i * 3 + 1] *= ANGS_PER_NM;
			ts->pos[i * 3 + 2] *= ANGS_PER_NM;
		}
        }

        if (hdr->v_size) {
		for (i = 0; i < hdr->natoms; i++) {
			if (trx_rvector(mf, NULL) < 0) {
				mdio_tsfree(ts, 1);
				return -1;
			}
		}
        }

        if (hdr->f_size) {
		for (i = 0; i < hdr->natoms; i++) {
			if (trx_rvector(mf, NULL) < 0) {
				mdio_tsfree(ts, 1);
				return -1;
			}
		}
        }

	return mdio_seterror(MDIO_SUCCESS);
}


// writes an int in big endian. Returns GMX_SUCCESS
// on success or a negative number on error.
static int put_trx_int(md_file *mf, int y) {
      if (!mf) return mdio_seterror(MDIO_BADPARAMS);

      // sanity check.
      if (sizeof(int) != 4) return mdio_seterror(MDIO_SIZEERROR);

      if (mf->rev) swap4_aligned(&y, 1);
      if (fwrite(&y, 4, 1, mf->f) != 1)
    return mdio_seterror(MDIO_IOERROR);

  return mdio_seterror(MDIO_SUCCESS);
}

// writes a real in big-endian. Returns GMX_SUCCESS
// on success or a negative number on error.
static int put_trx_real(md_file *mf, float y) {
      if (!mf) return mdio_seterror(MDIO_BADPARAMS);

      if (mf->rev) swap4_aligned(&y, 1);
      if (fwrite(&y, 4, 1, mf->f) != 1)
        return mdio_seterror(MDIO_IOERROR);

      return mdio_seterror(MDIO_SUCCESS);
}


// writes an xdr encoded string. Returns GMX_SUCCESS
// on success or a negative number on error.
static int put_trx_string(md_file *mf, const char *s) {
	if (!mf || !s) return mdio_seterror(MDIO_BADPARAMS);
        
        // write: size of object, string length, string data
        size_t len = strlen(s);
        if ( put_trx_int(mf, len+1)
             || put_trx_int(mf, len)
             || (fwrite(s, len, 1, mf->f) != 1))
          return mdio_seterror(MDIO_IOERROR);

	return mdio_seterror(MDIO_SUCCESS);
}


// xtc_int() - reads an integer from an xtc file
static int xtc_int(md_file *mf, int *i) {
	unsigned char c[4];

	if (!mf) return mdio_seterror(MDIO_BADPARAMS);
        // sanity check.
        if (sizeof(int) != 4) return mdio_seterror(MDIO_SIZEERROR);

	if (fread(c, 1, 4, mf->f) != 4) {
		if (feof(mf->f)) return mdio_seterror(MDIO_EOF);
		else if (ferror(mf->f)) return mdio_seterror(MDIO_IOERROR);
		else return mdio_seterror(MDIO_UNKNOWNERROR);
	}

	if (i) *i = c[3] + (c[2] << 8) + (c[1] << 16) + (c[0] << 24);
	return mdio_seterror(MDIO_SUCCESS);
}


// xtc_float() - reads a float from an xtc file
static int xtc_float(md_file *mf, float *f) {
	unsigned char c[4];
	int i;

	if (!mf) return mdio_seterror(MDIO_BADPARAMS);

	if (fread(c, 1, 4, mf->f) != 4) {
		if (feof(mf->f)) return mdio_seterror(MDIO_EOF);
		else if (ferror(mf->f)) return mdio_seterror(MDIO_IOERROR);
		else return mdio_seterror(MDIO_UNKNOWNERROR);
	}

	if (f) {
		// By reading the number in as an integer and then
		// copying it to a floating point number we can
		// ensure proper endianness
		i = c[3] + (c[2] << 8) + (c[1] << 16) + (c[0] << 24);
		memcpy(f, &i, 4);
	}
	return mdio_seterror(MDIO_SUCCESS);
}


// xtc_data() - reads a specific amount of data from an xtc
// file using the xdr format.
static int xtc_data(md_file *mf, char *buf, int len) {
	if (!mf || len < 1) return mdio_seterror(MDIO_BADPARAMS);
  size_t slen = (size_t)len;
	if (buf) {
		if (fread(buf, 1, slen, mf->f) != slen) {
			if (feof(mf->f)) return mdio_seterror(MDIO_EOF);
			if (ferror(mf->f)) return mdio_seterror(MDIO_IOERROR);
			else return mdio_seterror(MDIO_UNKNOWNERROR);
		}
		if (len % 4) {
			if (fseek(mf->f, 4 - (len % 4), SEEK_CUR)) {
				if (feof(mf->f)) return mdio_seterror(MDIO_EOF);
				if (ferror(mf->f)) return mdio_seterror(MDIO_IOERROR);
				else return mdio_seterror(MDIO_UNKNOWNERROR);
			}
		}
	}
	else {
		int newlen;
		newlen = len;
		if (len % 4) newlen += (4 - (len % 4));
		if (fseek(mf->f, newlen, SEEK_CUR)) {
			if (feof(mf->f)) return mdio_seterror(MDIO_EOF);
			if (ferror(mf->f)) return mdio_seterror(MDIO_IOERROR);
			else return mdio_seterror(MDIO_UNKNOWNERROR);
		}
	}
	return len;
}


// xtc_timestep() - reads a timestep from an .xtc file.
static int xtc_timestep(md_file *mf, md_ts *ts) {
	int n;
	float f, x[3], y[3], z[3];

	int size = 0; // explicitly initialized to zero.
	float precision;

	if (!mf || !ts) return mdio_seterror(MDIO_BADPARAMS);
	if (!mf->f) return mdio_seterror(MDIO_BADPARAMS);
	if (mf->fmt != MDFMT_XTC) return mdio_seterror(MDIO_WRONGFORMAT);

	// Check magic number
	if (xtc_int(mf, &n) < 0) return -1;
	if (n != XTC_MAGIC) return mdio_seterror(MDIO_BADFORMAT);

	// Get number of atoms
	if (xtc_int(mf, &n) < 0) return -1;
	ts->natoms = n;

	// Get the simulation step
	if (xtc_int(mf, &n) < 0) return -1;
	ts->step = n;

	// Get the time value
	if (xtc_float(mf, &f) < 0) return -1;
	ts->time = f;

	// Read the basis vectors of the box
  if ( (xtc_float(mf, &x[0]) < 0) ||
       (xtc_float(mf, &x[1]) < 0) ||
       (xtc_float(mf, &x[2]) < 0) ||
       (xtc_float(mf, &y[0]) < 0) ||
       (xtc_float(mf, &y[1]) < 0) ||
       (xtc_float(mf, &y[2]) < 0) ||
       (xtc_float(mf, &z[0]) < 0) ||
       (xtc_float(mf, &z[1]) < 0) ||
       (xtc_float(mf, &z[2]) < 0) )
    return -1;
  // Allocate the box and convert the vectors.
  ts->box = (md_box *) malloc(sizeof(md_box));
  if (mdio_readbox(ts->box, x, y, z) < 0) {
    free(ts->box);
    ts->box = NULL;
    return -1;
  }

	ts->pos = (float *) malloc(sizeof(float) * 3 * ts->natoms);
	if (!ts->pos) return mdio_seterror(MDIO_BADMALLOC);
	n = xtc_3dfcoord(mf, ts->pos, &size, &precision);
	if (n < 0) return -1;

	/* Now we're left with the job of scaling... */
	for (n = 0; n < ts->natoms * 3; n++)
		ts->pos[n] *= ANGS_PER_NM;

	return mdio_seterror(MDIO_SUCCESS);
}


///////////////////////////////////////////////////////////////////////
// This algorithm is an implementation of the 3dfcoord algorithm
// written by Frans van Hoesel (hoesel@chem.rug.nl) as part of the
// Europort project in 1995.
///////////////////////////////////////////////////////////////////////

// integer table used in decompression
static int xtc_magicints[] = {
	0, 0, 0, 0, 0, 0, 0, 0, 0,8, 10, 12, 16, 20, 25, 32, 40, 50, 64,
	80, 101, 128, 161, 203, 256, 322, 406, 512, 645, 812, 1024, 1290,
	1625, 2048, 2580, 3250, 4096, 5060, 6501, 8192, 10321, 13003, 16384,
	20642, 26007, 32768, 41285, 52015, 65536, 82570, 104031, 131072,
	165140, 208063, 262144, 330280, 416127, 524287, 660561, 832255,
	1048576, 1321122, 1664510, 2097152, 2642245, 3329021, 4194304,
	5284491, 6658042, 8388607, 10568983, 13316085, 16777216 };

#define FIRSTIDX 9
/* note that magicints[FIRSTIDX-1] == 0 */
#define LASTIDX (sizeof(xtc_magicints) / sizeof(*xtc_magicints))


// returns the number of bits in the binary expansion of
// the given integer.
static int xtc_sizeofint(int size) {
	unsigned int num = 1;
  unsigned int ssize = (unsigned int)size;
	int nbits = 0;

	while (ssize >= num && nbits < 32) {
		nbits++;
		num <<= 1;
	}
	return nbits;
}

// calculates the number of bits a set of integers, when compressed,
// will take up.
static int xtc_sizeofints(int nints, unsigned int *sizes) {
	int i;
  unsigned int num;
	unsigned int nbytes, nbits, bytes[32], bytecnt, tmp;
	nbytes = 1;
	bytes[0] = 1;
	nbits = 0;
	for (i=0; i < nints; i++) {	
		tmp = 0;
		for (bytecnt = 0; bytecnt < nbytes; bytecnt++) {
			tmp = bytes[bytecnt] * sizes[i] + tmp;
			bytes[bytecnt] = tmp & 0xff;
			tmp >>= 8;
		}
		while (tmp != 0) {
			bytes[bytecnt++] = tmp & 0xff;
			tmp >>= 8;
		}
		nbytes = bytecnt;
	}
	num = 1;
	nbytes--;
	while (bytes[nbytes] >= num) {
		nbits++;
		num *= 2;
	}
	return nbits + nbytes * 8;
}

// reads bits from a buffer.    
static int xtc_receivebits(int *buf, int nbits) {
	int cnt, num; 
	unsigned int lastbits, lastbyte;
	unsigned char * cbuf;
	int mask = (1 << nbits) -1;

	cbuf = ((unsigned char *)buf) + 3 * sizeof(*buf);
	cnt = buf[0];
	lastbits = (unsigned int) buf[1];
	lastbyte = (unsigned int) buf[2];

	num = 0;
	while (nbits >= 8) {
		lastbyte = ( lastbyte << 8 ) | cbuf[cnt++];
		num |=  (lastbyte >> lastbits) << (nbits - 8);
		nbits -=8;
	}
	if (nbits > 0) {
		if (lastbits < (unsigned int)nbits) {
			lastbits += 8;
			lastbyte = (lastbyte << 8) | cbuf[cnt++];
		}
		lastbits -= nbits;
		num |= (lastbyte >> lastbits) & ((1 << nbits) -1);
	}
	num &= mask;
	buf[0] = cnt;
	buf[1] = lastbits;
	buf[2] = lastbyte;
	return num; 
}

// decompresses small integers from the buffer
static void xtc_receiveints(int *buf, const int nints, int nbits,
			unsigned int *sizes, int *nums) {
	int bytes[32];
	int i, j, nbytes, p, num;

	bytes[1] = bytes[2] = bytes[3] = 0;
	nbytes = 0;
	while (nbits > 8) {
		bytes[nbytes++] = xtc_receivebits(buf, 8);
		nbits -= 8;
	}
	if (nbits > 0) {
		bytes[nbytes++] = xtc_receivebits(buf, nbits);
	}
	for (i = nints-1; i > 0; i--) {
		num = 0;
		for (j = nbytes-1; j >=0; j--) {
			num = (num << 8) | bytes[j];
			p = num / sizes[i];
			bytes[j] = p;
			num = num - p * sizes[i];
		}
		nums[i] = num;
	}
	nums[0] = bytes[0] | (bytes[1] << 8) | (bytes[2] << 16) | (bytes[3] << 24);
}

// function that actually reads and writes compressed coordinates    
static int xtc_3dfcoord(md_file *mf, float *fp, int *size, float *precision) {
	static int *ip = NULL;
	static int oldsize;
	static int *buf;

	int minint[3], maxint[3], *lip;
	int smallidx;
	unsigned sizeint[3], sizesmall[3], bitsizeint[3], size3;
	int flag, k;
	int small, smaller, i, is_smaller, run;
	float *lfp;
	int tmp, *thiscoord,  prevcoord[3];

	int bufsize, lsize;
	unsigned int bitsize;
	float inv_precision;


	if (xtc_int(mf, &lsize) < 0) return -1;

	if (*size != 0 && lsize != *size) return mdio_seterror(MDIO_BADFORMAT);

	*size = lsize;
	size3 = *size * 3;
	if (*size <= 9) {
		for (i = 0; i < *size; i++) {
			if (xtc_float(mf, fp + (3 * i)) < 0) return -1;
			if (xtc_float(mf, fp + (3 * i) + 1) < 0) return -1;
			if (xtc_float(mf, fp + (3 * i) + 2) < 0) return -1;
		}
		return *size;
	}
	xtc_float(mf, precision);
	if (ip == NULL) {
		ip = (int *)malloc(size3 * sizeof(*ip));
		if (ip == NULL) return mdio_seterror(MDIO_BADMALLOC);
		bufsize = (int) (size3 * 1.2);
		buf = (int *)malloc(bufsize * sizeof(*buf));
		if (buf == NULL) return mdio_seterror(MDIO_BADMALLOC);
		oldsize = *size;
	} else if (*size > oldsize) {
		ip = (int *)realloc(ip, size3 * sizeof(*ip));
		if (ip == NULL) return mdio_seterror(MDIO_BADMALLOC);
		bufsize = (int) (size3 * 1.2);
		buf = (int *)realloc(buf, bufsize * sizeof(*buf));
		if (buf == NULL) return mdio_seterror(MDIO_BADMALLOC);
		oldsize = *size;
	}
	buf[0] = buf[1] = buf[2] = 0;

	xtc_int(mf, &(minint[0]));
	xtc_int(mf, &(minint[1]));
	xtc_int(mf, &(minint[2]));

	xtc_int(mf, &(maxint[0]));
	xtc_int(mf, &(maxint[1]));
	xtc_int(mf, &(maxint[2]));
		
	sizeint[0] = maxint[0] - minint[0]+1;
	sizeint[1] = maxint[1] - minint[1]+1;
	sizeint[2] = maxint[2] - minint[2]+1;
	
	/* check if one of the sizes is to big to be multiplied */
	if ((sizeint[0] | sizeint[1] | sizeint[2] ) > 0xffffff) {
		bitsizeint[0] = xtc_sizeofint(sizeint[0]);
		bitsizeint[1] = xtc_sizeofint(sizeint[1]);
		bitsizeint[2] = xtc_sizeofint(sizeint[2]);
		bitsize = 0; /* flag the use of large sizes */
	} else {
		bitsize = xtc_sizeofints(3, sizeint);
	}

	xtc_int(mf, &smallidx);
	smaller = xtc_magicints[FIRSTIDX > smallidx - 1 ? FIRSTIDX : smallidx - 1] / 2;
	small = xtc_magicints[smallidx] / 2;
	sizesmall[0] = sizesmall[1] = sizesmall[2] = xtc_magicints[smallidx] ;

	/* buf[0] holds the length in bytes */

	if (xtc_int(mf, &(buf[0])) < 0) return -1;

	if (xtc_data(mf, (char *) &buf[3], (int) buf[0]) < 0) return -1;

	buf[0] = buf[1] = buf[2] = 0;

	lfp = fp;
	inv_precision = 1.0f / (*precision);
	run = 0;
	i = 0;
	lip = ip;
	while (i < lsize) {
		thiscoord = (int *)(lip) + i * 3;

		if (bitsize == 0) {
			thiscoord[0] = xtc_receivebits(buf, bitsizeint[0]);
			thiscoord[1] = xtc_receivebits(buf, bitsizeint[1]);
			thiscoord[2] = xtc_receivebits(buf, bitsizeint[2]);
		} else {
			xtc_receiveints(buf, 3, bitsize, sizeint, thiscoord);
		}

		i++;
		thiscoord[0] += minint[0];
		thiscoord[1] += minint[1];
		thiscoord[2] += minint[2];

		prevcoord[0] = thiscoord[0];
		prevcoord[1] = thiscoord[1];
		prevcoord[2] = thiscoord[2];
 

		flag = xtc_receivebits(buf, 1);
		is_smaller = 0;
		if (flag == 1) {
			run = xtc_receivebits(buf, 5);
			is_smaller = run % 3;
			run -= is_smaller;
			is_smaller--;
		}
		if (run > 0) {
			thiscoord += 3;
			for (k = 0; k < run; k+=3) {
				xtc_receiveints(buf, 3, smallidx, sizesmall, thiscoord);
				i++;
				thiscoord[0] += prevcoord[0] - small;
				thiscoord[1] += prevcoord[1] - small;
				thiscoord[2] += prevcoord[2] - small;
				if (k == 0) {
					/* interchange first with second atom for better
					 * compression of water molecules
					 */
					tmp = thiscoord[0]; thiscoord[0] = prevcoord[0];
					prevcoord[0] = tmp;
					tmp = thiscoord[1]; thiscoord[1] = prevcoord[1];
					prevcoord[1] = tmp;
					tmp = thiscoord[2]; thiscoord[2] = prevcoord[2];
					prevcoord[2] = tmp;
					*lfp++ = prevcoord[0] * inv_precision;
					*lfp++ = prevcoord[1] * inv_precision;
					*lfp++ = prevcoord[2] * inv_precision;
				} else {
					prevcoord[0] = thiscoord[0];
					prevcoord[1] = thiscoord[1];
					prevcoord[2] = thiscoord[2];
				}
				*lfp++ = thiscoord[0] * inv_precision;
				*lfp++ = thiscoord[1] * inv_precision;
				*lfp++ = thiscoord[2] * inv_precision;
			}
		} else {
			*lfp++ = thiscoord[0] * inv_precision;
			*lfp++ = thiscoord[1] * inv_precision;
			*lfp++ = thiscoord[2] * inv_precision;		
		}
		smallidx += is_smaller;
		if (is_smaller < 0) {
			small = smaller;
			if (smallidx > FIRSTIDX) {
				smaller = xtc_magicints[smallidx - 1] /2;
			} else {
				smaller = 0;
			}
		} else if (is_smaller > 0) {
			smaller = small;
			small = xtc_magicints[smallidx] / 2;
		}
		sizesmall[0] = sizesmall[1] = sizesmall[2] = xtc_magicints[smallidx] ;
	}
	return 1;
}
#endif
