#ifdef _PYMOL_WIN32
#include"os_predef.h"
#endif
#if 0 
/* we no longer support numpy with this ancient version of PyOpenGL */
def _PYMOL_NUMPY
*/

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

