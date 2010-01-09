

/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* Copyright (c) Schrodinger, LLC. 
D* -------------------------------------------------------------------
E* It is unlawful to modify or remove this copyright notice.
F* -------------------------------------------------------------------
G* Please see the accompanying LICENSE file for further information. 
H* -------------------------------------------------------------------
I* Additional authors of this source file include:
-* 
-* 
-*
Z* -------------------------------------------------------------------
*/
#ifndef _H_Raw
#define _H_Raw

#include"PyMOLGlobals.h"


/* interface class for device-independent, raw binary I/O */

typedef struct {
  PyMOLGlobals *G;
  int mode;
  FILE *f;
  char *bufVLA;
  int swap;
  int header[4];
} CRaw;


/* object types */

#define cRaw_EOF                     0
#define cRaw_AtomInfo1               1
#define cRaw_Coords1                 2
#define cRaw_Spheroid1               3
#define cRaw_SpheroidNormals1        4
#define cRaw_SpheroidInfo1           5
#define cRaw_Bonds1                  6


/* FYI: file format;

4 byte magic for determining endianess 0x04030201

followed by any number of records with the following header...

4 bytes (record length, less header) 
4 bytes (record type)
4 bytes (program version)
4 bytes (serial number)
...

*/

CRaw *RawOpenRead(PyMOLGlobals * G, char *fname);
CRaw *RawOpenWrite(PyMOLGlobals * G, char *fname);
CRaw *RawOpenAppend(PyMOLGlobals * G, char *fname);
void RawFree(CRaw * I);

int RawWrite(CRaw * I, int type, unsigned int size, int serial, char *bytes);

int RawGetNext(CRaw * I, int *size, int *version);

char *RawRead(CRaw * I, int *type, unsigned int *size, int *serial);
char *RawReadPtr(CRaw * I, int type, int *size);
char *RawReadVLA(CRaw * I, int type, unsigned int rec_size, int grow_factor,
                 int auto_zero);
int RawReadInto(CRaw * I, int type, unsigned int size, char *buffer);
int RawReadSkip(CRaw * I);

#endif
