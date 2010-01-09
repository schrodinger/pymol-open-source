

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

#ifndef _H_Queue
#define _H_Queue

typedef struct {
  char *ptr;
  unsigned int inp, out;
  unsigned int mask, size;
} CQueue;

CQueue *QueueNew(PyMOLGlobals * G, unsigned int mask);
void QueueFree(CQueue * I);
int QueueStrCheck(CQueue * I);

void QueueStrIn(CQueue * I, char *c);
int QueueStrOut(CQueue * I, char *c);

#endif
