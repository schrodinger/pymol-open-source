

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
#ifndef _H_Util
#define _H_Util

#include "os_predef.h"

#include "PyMOLGlobals.h"
#include "Base.h"
#include <string>

void UtilZeroMem(void *ptr, ov_size howMuch);
void UtilCopyMem(void *dst, const void *src, ov_size howMuch);
void *UtilArrayCalloc(unsigned int *dim, ov_size ndim, ov_size atom_size);
char *UtilConcat(char *where,const char *what);
void UtilNConcat(char *dst, const char *str, ov_size n);
void UtilConcatVLA(char **vla, ov_size * cc, const char *str);
void UtilNPadVLA(char **vla, ov_size * cc, const char *str, ov_size len);
void UtilFillVLA(char **vla, ov_size * cc, char what, ov_size len);
void UtilNCopy(char *dst, const char *src, ov_size n);        /* up to N-1 chars */
void UtilNCopyToLower(char *dst, const char *src, ov_size n); /* up to N-1 chars */
void UtilCleanStr(char *s);
std::string UtilCleanStdStr(const std::string& s);
void UtilStripANSIEscapes(char *s);
void UtilStripANSIEscapes(std::string& str);
int UtilCountStringVLA(char *vla);

double UtilGetSecondsEpoch();
double UtilGetSeconds(PyMOLGlobals * G);
int UtilInit(PyMOLGlobals * G);
void UtilFree(PyMOLGlobals * G);

typedef int UtilOrderFn(const void *array, int l, int r);
void UtilSortIndex(int n, void *array, int *x, UtilOrderFn * fOrdered);

int UtilSemiSortFloatIndex(int n, float *array, int *x, int forward);
int UtilSemiSortFloatIndexWithNBins(int n, int nbins, float *array,int *x, int forward);
int UtilSemiSortFloatIndexWithNBinsImpl(int *start1, int n, int nbins, float *array, int *destx, int forward);

void UtilApplySortedIndices(int n, int *x, int rec_size, void *src, void *dst);

void UtilSortInPlace(PyMOLGlobals * G, void *array, int nItem, unsigned int itemSize,
                     UtilOrderFn * fOrdered);

void UtilExpandArrayElements(void *src, void *dst, int n_entries, int old_rec_size,
                             int new_rec_size);
typedef int UtilOrderFnGlobals(PyMOLGlobals * G, const void *array, int l, int r);
void UtilSortIndexGlobals(PyMOLGlobals * G, int n, const void *array, int *x,
                          UtilOrderFnGlobals * fOrdered);

int UtilShouldWePrintQuantity(int quantity);

#endif
