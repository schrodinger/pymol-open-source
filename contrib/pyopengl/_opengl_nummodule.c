#ifdef _PYMOL_WIN32
#include"os_predef.h"
#endif
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

/* for distutils compatibility on WIN32 */
DL_EXPORT(void)
init_opengl_nummodule(void) {init_opengl_num();}

#endif

