#ifdef _PYMOL_NUMPY

#define NUMERIC
#include "_openglmodule.c"

#else

#ifdef WIN32
#include<windows.h>
#endif

#include<Python.h>

DL_EXPORT(void)
init_opengl_num(void) {}

#endif
