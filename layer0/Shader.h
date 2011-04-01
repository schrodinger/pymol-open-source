

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
#if 0
#ifndef _H_Shader
#define _H_Shader

#include"PyMOLGlobals.h"

int ShaderInit(PyMOLGlobals *G);

int ShaderEnable(PyMOLGlobals *G);
int ShaderDisable(PyMOLGlobals *G);



/* getGLVersion -- determine user's GL version
 * 
 * RETURNS
 * writes to major and minor
 */
void getGLVersion(PyMOLGlobals * G, int * major, int * minor);

/* getGLSLVersion -- determine user's GLSL version
 *
 * RETURNS
 * writes to major and minor 
 */
void getGLSLVersion(PyMOLGlobals * G, int * major, int * minor);


/* allocate memory for the shader */
CShaderMgr* ShaderMgr_Alloc();

/* intialize the shader object and 
 * initialize shaders on this machine */
void ShaderMgr_Init(CShaderMgr * s, PyMOLGlobals * G);

/* try to init the shaders on this machine */
int ShaderMgr_InitShaders(CShaderMgr * s);

/* clean up */
void ShaderMgr_Free(CShaderMgr * s);

#endif
#endif
