#ifdef _PYMOL_NUMPY

#define NUMERIC
#include "_glumodule.c"

#else
#ifdef WIN32
#include<windows.h>
#endif
#include<Python.h>
DL_EXPORT(void)
init_glu_num(void) {};
#endif
 