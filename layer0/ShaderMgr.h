/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* copyright Schrodinger, LLC.
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
#ifndef _SHADER_MGR_H
#define _SHADER_MGR_H

#include "os_gl.h"
#include "PyMOLGlobals.h"
#include "Executive_pre.h"
#include "OVContext.h"
#include "Rep.h"

#ifndef GL_FRAGMENT_PROGRAM_ARB
#define GL_FRAGMENT_PROGRAM_ARB                         0x8804
#endif

/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
#ifdef __WIN32
PFNGLTEXIMAGE3DPROC getTexImage3D();
static PFNGLTEXIMAGE3DPROC glTexImage3D;
static PFNGLACTIVETEXTUREPROC glActiveTexture = -100l;
static PFNGLGENPROGRAMSARBPROC glGenProgramsARB;
static PFNGLBINDPROGRAMARBPROC glBindProgramARB;
static PFNGLDELETEPROGRAMSARBPROC glDeleteProgramsARB;
static PFNGLPROGRAMSTRINGARBPROC glProgramStringARB;
static PFNGLPROGRAMENVPARAMETER4FARBPROC glProgramEnvParameter4fARB;
static PFNGLGETPROGRAMIVARBPROC glGetProgramivARB;
static PFNGLGETPROGRAMIVARBPROC glGetProgramiv;
static PFNGLATTACHSHADERPROC glAttachShader;
static PFNGLGETSHADERINFOLOGPROC glGetShaderInfoLog;
static PFNGLGETSHADERIVPROC glGetShaderiv;
static PFNGLCOMPILESHADERPROC glCompileShader;
static PFNGLSHADERSOURCEPROC glShaderSource;
static PFNGLCREATESHADERPROC glCreateShader;
static PFNGLCREATEPROGRAMPROC glCreateProgram;
static PFNGLDELETEPROGRAMPROC glDeleteProgram;
static PFNGLDELETESHADERPROC glDeleteShader;
static PFNGLUSEPROGRAMPROC glUseProgram;
static PFNGLGETPROGRAMINFOLOGPROC glGetProgramInfoLog;
static PFNGLLINKPROGRAMPROC glLinkProgram;
static PFNGLUNIFORM1IPROC glUniform1i;
static PFNGLUNIFORM1FPROC glUniform1f;
static PFNGLGETUNIFORMLOCATIONPROC glGetUniformLocation;
#endif

/* END PROPRIETARY CODE SEGMENT */
typedef struct _CShaderPrg {
  PyMOLGlobals * G;
  /* for retrieving from the ShaderMgr */
  char * name;
  /* openGL assigned id */
  int id;
  /* openGL fragment and vertex shader ids */
  GLuint vid;
  GLuint fid;
  /* fragment and vertex source */
  char * f;
  char * v;
  struct _CShaderPrg * prev;
  struct _CShaderPrg * next;
  int uniform_set ;
} CShaderPrg;

struct _CShaderMgr { 
  PyMOLGlobals * G;
  CShaderPrg * programs;
  int ShadersPresent;
  GLuint *vbos_to_free;
  int number_of_vbos_to_free;
  CShaderPrg *current_shader;
  int stereo_flag; /* -1 left; 0 = off; 1 = right */
  int is_picking;

  /* These lookups and arrays are used to dynamically change the shaders
     based on preprocessor-like statements within the shaders.  These need
     string lookups. */
  OVLexicon *ShaderLex;
  OVOneToOne *ShaderLexLookup;
  char **shader_replacement_strings;
  int *shader_include_values;
  char **shader_update_when_include_filename;
  char **shader_update_when_include;
  char **shader_update_replace_with;
  short print_warnings;
  int reload_bits;
} _CShaderMgr;

#include "OVreturns.h"

/* ShaderMgrInit -- Called from PyMOL.c, allocates and attaches to G->ShaderMgr
 */
OVstatus ShaderMgrInit(PyMOLGlobals * G);

/* ShaderMgrConfig -- Called from PyMOL.c, configures the global ShaderMgr;
 */
void ShaderMgrConfig(PyMOLGlobals * G);

/* ============================================================================
 * CShaderMgr class -- simple ShaderMgr for PyMOL
 * ============================================================================*/

/* New */
CShaderMgr * CShaderMgr_New(PyMOLGlobals * G);
CShaderPrg * CShaderPrg_NewARB(PyMOLGlobals * G, const char * name, const char * v, const char * f);
/* Delete */
void CShaderMgr_Delete(CShaderMgr * I);
void CShaderMgrFree(PyMOLGlobals *G);

/* AddShader -- Adds to the global shader library */ 
int CShaderMgr_AddShaderPrg(CShaderMgr * I, CShaderPrg * s);

/* RemoveShader -- Removes shader program by name */
int CShaderMgr_RemoveShaderPrg(CShaderMgr * I, const char * name);

/* GetShaderPrg -- gets a ptr to the installed shader */
CShaderPrg * CShaderMgr_GetShaderPrg(CShaderMgr * I, const char * name);
CShaderPrg * CShaderMgr_GetShaderPrg_NoSet(CShaderMgr * I, const char * name);
int CShaderMgr_ShaderPrgExists(CShaderMgr * I, const char * name);

/* runtime check for shaders */
int CShaderMgr_ShadersPresent(CShaderMgr * I);

/* ============================================================================
 * CShaderPrg class -- a simple facade that wraps a shader and shader program
 * into a simple API
 * ============================================================================*/

/* New -- from shader source strings */
CShaderPrg * CShaderPrg_New(PyMOLGlobals * G, const char * name, const char * v, const char * f);

/* New -- from files on disk */
CShaderPrg * CShaderPrg_NewFromFile(PyMOLGlobals * G, const char * name, const char * vFile, const char * fFile);

/* Delete */
void CShaderPrg_Delete(CShaderPrg * I);

/* Enable */
int CShaderPrg_Enable(CShaderPrg * I);

/* Disable */
int CShaderPrg_Disable(CShaderPrg * I);
int CShaderPrg_DisableARB(CShaderPrg *I);

/* Link and IsLinked */
int CShaderPrg_Link(CShaderPrg * I);
int CShaderPrg_IsLinked(CShaderPrg * I);

/* accessors/mutators/uniform setters */
int CShaderPrg_Set1i(CShaderPrg * I, const char * name, int i);
int CShaderPrg_Set1f(CShaderPrg * I, const char * name, float f);
int CShaderPrg_Set2f(CShaderPrg * I, const char * name, float f1, float f2);
int CShaderPrg_Set3f(CShaderPrg * I, const char * name, float f1, float f2, float f3);
int CShaderPrg_Set4fv(CShaderPrg * p, const char * name, float *f);
int CShaderPrg_Set3fv(CShaderPrg * p, const char * name, float *f);
int CShaderPrg_Set4f(CShaderPrg * p, const char * name, float f1, float f2, float f3, float f4);
int CShaderPrg_SetMat3f(CShaderPrg * p, const char * name, GLfloat * m);
int CShaderPrg_SetMat4f(CShaderPrg * p, const char * name, GLfloat * m);
//int CShaderPrg_SetTexture2D(CShaderPrg * p, const char * name, GLuint );

int CShaderPrg_SetMat3fc(CShaderPrg * p, const char * name, GLfloat * m);
int CShaderPrg_SetMat4fc(CShaderPrg * p, const char * name, GLfloat * m);

int CShaderPrg_GetAttribLocation(CShaderPrg * p, const char * name);
void CShaderPrg_SetAttrib4fLocation(CShaderPrg * p, const char * name, float f1, float f2, float f3, float f4);
void CShaderPrg_SetAttrib1fLocation(CShaderPrg * p, const char * name, float f1);

int CShaderPrg_SetLightingEnabled(CShaderPrg * p, int);

void CShaderMgr_AddVBOsToFree(CShaderMgr * I, GLuint *vboid, int nvbos);
void CShaderMgr_AddVBOToFree(CShaderMgr * I, GLuint vboid);
void CShaderMgr_FreeAllVBOs(CShaderMgr * I);

CShaderPrg *CShaderPrg_Enable_DefaultShader(PyMOLGlobals * G);
CShaderPrg *CShaderPrg_Enable_DefaultShaderWithSettings(PyMOLGlobals * G, CSetting * set1, CSetting * set2);
CShaderPrg *CShaderPrg_Enable_DefaultScreenShader(PyMOLGlobals * G);

CShaderPrg *CShaderPrg_Enable_CylinderShader(PyMOLGlobals * G);
CShaderPrg *CShaderPrg_Enable_DefaultSphereShader(PyMOLGlobals * G);

CShaderPrg *CShaderPrg_Enable_SphereShader(PyMOLGlobals * G, char *name);
CShaderPrg *CShaderPrg_Enable_SphereShaderARB(PyMOLGlobals * G);
CShaderPrg *CShaderPrg_Enable_RampShader(PyMOLGlobals * G);

CShaderPrg *CShaderPrg_Get_ScreenShader(PyMOLGlobals * G);
CShaderPrg *CShaderPrg_Get_DefaultShader(PyMOLGlobals * G);
CShaderPrg *CShaderPrg_Get_DefaultScreenShader(PyMOLGlobals * G);
CShaderPrg *CShaderPrg_Get_CylinderShader(PyMOLGlobals * G);
CShaderPrg *CShaderPrg_Get_CylinderShader_NoSet(PyMOLGlobals * G);
CShaderPrg *CShaderPrg_Get_DefaultSphereShader(PyMOLGlobals * G);
CShaderPrg *CShaderPrg_Get_RampShader(PyMOLGlobals * G);

CShaderPrg *CShaderPrg_Get_Current_Shader(PyMOLGlobals * G);

CShaderPrg *CShaderPrg_Get_IndicatorShader(PyMOLGlobals * G);
CShaderPrg *CShaderPrg_Enable_IndicatorShader(PyMOLGlobals * G);

CShaderPrg *CShaderPrg_Enable_BackgroundShader(PyMOLGlobals * G);
CShaderPrg *CShaderPrg_Get_BackgroundShader(PyMOLGlobals * G);

CShaderPrg *CShaderPrg_Enable_ScreenShader(PyMOLGlobals * G);
CShaderPrg *CShaderPrg_Enable_LabelScreenShader(PyMOLGlobals * G);
CShaderPrg *CShaderPrg_Enable_LabelShader(PyMOLGlobals * G);
CShaderPrg *CShaderPrg_Get_LabelShader(PyMOLGlobals * G);
CShaderPrg *CShaderPrg_Get_LabelScreenShader(PyMOLGlobals * G);

void ShaderMgrResetUniformSet(PyMOLGlobals * G);

void CShaderPrg_SetIsPicking(PyMOLGlobals * G, int is_picking);
void CShaderPrg_SetFogUniforms(PyMOLGlobals * G, CShaderPrg * shaderPrg);
void CShaderPrg_Update_Shaders_For_Background(PyMOLGlobals *G);
int CShaderPrg_Reload(PyMOLGlobals * G, char *name, char *v, char *f);
void CShaderMgr_Reload_Default_Shader(PyMOLGlobals *G);
void CShaderMgr_Reload_Cylinder_Shader(PyMOLGlobals *G);
void CShaderMgr_Reload_Sphere_Shader(PyMOLGlobals *G);
void CShaderMgr_Reload_Indicator_Shader(PyMOLGlobals *G);
void CShaderMgr_Reload_Shader_Variables(PyMOLGlobals * G);
void CShaderPrg_Reload_CallComputeColorForLight(PyMOLGlobals * G, char *name);
void CShaderPrg_Reload_All_Shaders_For_CallComputeColorForLight(PyMOLGlobals * G);
void CShaderPrg_Reload_All_Shaders(PyMOLGlobals * G);
void CShaderMgr_Set_Reload_Bits(PyMOLGlobals *G, int bits);
void CShaderMgr_Check_Reload(PyMOLGlobals *G);

/* for reload_bits */
//RELOAD_ALL_SHADERS calls CShaderPrg_Reload_All_Shaders
#define RELOAD_ALL_SHADERS 1
//RELOAD_SHADERS_FOR_LIGHTING calls CShaderPrg_Reload_All_Shaders_For_CallComputeColorForLight
#define RELOAD_SHADERS_FOR_LIGHTING 2
//RELOAD_SHADERS_UPDATE_FOR_BACKGROUND calls CShaderPrg_Update_Shaders_For_Background
#define RELOAD_SHADERS_UPDATE_FOR_BACKGROUND 4
//RELOAD_SHADERS_CYLINDER calls CShaderMgr_Reload_Cylinder_Shader
#define RELOAD_SHADERS_CYLINDER 8

#define VERTEX_POS_SIZE    3
#define VERTEX_COLOR_SIZE  4

#define VERTEX_POS    0
#define VERTEX_NORMAL 1
#define VERTEX_COLOR  2

#define CYLINDER_ORIGIN 0
#define CYLINDER_AXIS   1
#define CYLINDER_COLOR  2
#define CYLINDER_COLOR2 3
#endif

