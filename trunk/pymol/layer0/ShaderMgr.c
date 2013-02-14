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

#include "os_python.h"
#include <string.h>
#include "ShaderMgr.h"
#include "OOMac.h"
#include "ListMacros.h"
#include "PyMOLOptions.h"
#include "os_gl.h"
#include "Feedback.h"
#include "MemoryDebug.h"
#include "Setting.h"
#include "Scene.h"
#include "Color.h"

#define MAX_LOG_LEN 1024

#include "ShaderText.h"

void getGLVersion(PyMOLGlobals * G, int *major, int* minor);
void getGLSLVersion(PyMOLGlobals * G, int* major, int* minor);

static void disableShaders(PyMOLGlobals * G);

#ifdef WIN32
/* REMOVE US */
PFNGLTEXIMAGE3DPROC getTexImage3D(){
  static PFNGLTEXIMAGE3DPROC my_glTexImage3D = NULL;
  if (!my_glTexImage3D)
    my_glTexImage3D = (PFNGLTEXIMAGE3DPROC) wglGetProcAddress("glTexImage3D");
  return my_glTexImage3D;
}

PFNGLACTIVETEXTUREPROC getActiveTexture() {
  static PFNGLACTIVETEXTUREPROC my_glActiveTexture = NULL;
  if (!my_glActiveTexture)
    my_glActiveTexture = (PFNGLACTIVETEXTUREPROC) wglGetProcAddress("glActiveTexture");
  return my_glActiveTexture;
}
#endif

/*
 * Use this to turn off shaders if the renderer cannot use them.
 */
void disableShaders(PyMOLGlobals * G) {
    /* Auto-disable shader-based rendering */
    SettingSetGlobal_b(G, cSetting_use_shaders, 0);
    SettingSetGlobal_i(G, cSetting_sphere_mode, 0);
}

/* ============================================================================
 * ShaderMgrInit is called from PyMOL.c during start up; it just allocates
 * the global ShaderMgr
 */
void ShaderMgrInit(PyMOLGlobals * G) {
  G->ShaderMgr = CShaderMgr_New(G);
}

/* ShaderMgrConfig -- Called from PyMOL.c, configures the global ShaderMgr
 * This needs to be called once the OpenGL context has been created, it is 
 * called from MainInit() for PyMol, and from PyMOL_ConfigureShadersGL() for
 * other programs (i.e., JyMOL, AxPyMOL, etc.).
 */
void ShaderMgrConfig(PyMOLGlobals * G) {
  int major, minor;
  char buf[50];
  CShaderPrg *defaultShader, *volumeShader, *sphereShader,
    *cylinderShader, *cylinderNoFFShader, *spheredirectShader;

  int hasShaders = 0, hasCylinderShader = 0;
  int ok = 0;
  GLenum err; 

  ok = (G && G->HaveGUI); /* && G->ValidContext); */

#ifndef _PYMOL_PURE_OPENGL_ES
  if (ok) {
    err = glewInit();
  } else {
    return;
  }
#else
  err = -1;
#endif
#ifndef ANDROID
  if (GLEW_OK==err) {
    if (GLEW_VERSION_2_0) {
      FeedbackAdd(G, " Detected OpenGL version 2.0 or greater. Shaders available.\n");
    }
    else { 
      FeedbackAdd(G, " Detected OpenGL version prior to 2.0. Shaders and volumes unavailable.\n");
      disableShaders(G);
      return;
    }
  } 
  else {
    /* print info on glew error? */
    FeedbackAdd(G, " There was an error intializing GLEW.  Basic graphics, including\n shaders and volumes may be unavailable.\n");
    disableShaders(G);
#ifndef _PYMOL_PURE_OPENGL_ES
    fprintf(stderr, " GLEW-Error: %s\n", glewGetErrorString(err));
#endif
    return;
  }
#endif

  /* First try to load configuration from files in $PYMOL_PATH, if they don't exist 
     load from const char * */
  //  FeedbackEnable(G, FB_ShaderMgr, FB_Everything);
#if defined(PURE_OPENGL_ES_2)
  /* these are the shaders for the IPhone (Compiled but not used yet) */
  PRINTFB(G, FB_ShaderMgr, FB_Results) "reading in default_es.vs and default_es.fs\n" ENDFB(G);
  defaultShader = CShaderPrg_NewFromFile(G, "default", "default_es.vs", "default_es.fs");
#elif defined(OPENGL_ES_2)
  /* these are the shaders that uses ES2 */
  PRINTFB(G, FB_ShaderMgr, FB_Debugging) "reading in default_es2.vs and default_es2.fs\n" ENDFB(G);
  defaultShader = CShaderPrg_NewFromFile(G, "default", "default_es2.vs", "default_es2.fs");
#else
  /* these are the original default shaders */
  PRINTFB(G, FB_ShaderMgr, FB_Results) "reading in default.vs and default.fs\n" ENDFB(G);
  defaultShader = CShaderPrg_NewFromFile(G, "default", "default.vs", "default.fs");
#endif
  if (!defaultShader){
    PRINTFB(G, FB_ShaderMgr, FB_Results)
      " PyMOLShader_NewFromFile-Warning: default shader files not found, loading from memory.\n" ENDFB(G);
    defaultShader = CShaderPrg_New(G, "default", default_vs, default_fs);
  }

#if defined(OPENGL_ES_2)
  if (defaultShader){
    GLenum err ;
    glBindAttribLocation(defaultShader->id, VERTEX_POS, "a_Vertex");
    if ((err = glGetError())){
        PRINTFB(G, FB_ShaderMgr, FB_Warnings) "GLERROR: a_Vertex: %d\n", err ENDFB(G);
    }
    glBindAttribLocation(defaultShader->id, VERTEX_NORMAL, "a_Normal");
    if ((err = glGetError())){
        PRINTFB(G, FB_ShaderMgr, FB_Warnings) "GLERROR: a_Normal: %d\n", err ENDFB(G);
    }
    glBindAttribLocation(defaultShader->id, VERTEX_COLOR, "a_Color");
    if ((err = glGetError())){
        PRINTFB(G, FB_ShaderMgr, FB_Warnings) "GLERROR: a_Color: %d\n", err ENDFB(G);
    }
    CShaderPrg_Link(defaultShader);	  
  }
#endif

  hasShaders = (defaultShader!=0);

  CShaderMgr_AddShaderPrg(G->ShaderMgr, defaultShader);

#ifndef _PYMOL_PURE_OPENGL_ES
  volumeShader = CShaderPrg_NewFromFile(G, "volume", "volume.vs", "volume.fs");
  if (!volumeShader){
    PRINTFB(G, FB_ShaderMgr, FB_Warnings)
      " PyMOLShader_NewFromFile-Warning: volume shader files not found, loading from memory.\n" ENDFB(G);
    volumeShader = CShaderPrg_New(G, "volume", volume_vs, volume_fs);
  }
  hasShaders &= (volumeShader!=0);
  CShaderMgr_AddShaderPrg(G->ShaderMgr, volumeShader);

  sphereShader = CShaderPrg_NewFromFile(G, "sphere", "sphere.vs", "sphere.fs");
  if (!sphereShader){
    PRINTFB(G, FB_ShaderMgr, FB_Warnings)
      " PyMOLShader_NewFromFile-Warning: sphere shader files not found, loading from memory.\n" ENDFB(G);
    sphereShader = CShaderPrg_New(G, "sphere", sphere_vs, sphere_fs);
  }
  hasShaders &= (sphereShader!=0);
  CShaderMgr_AddShaderPrg(G->ShaderMgr, sphereShader);

  spheredirectShader = CShaderPrg_NewFromFile(G, "spheredirect", "spheredirect.vs", 0);
  if (!spheredirectShader){
    PRINTFB(G, FB_ShaderMgr, FB_Warnings)
      " PyMOLShader_NewFromFile-Warning: sphere shader files not found, loading from memory.\n" ENDFB(G);
    spheredirectShader = CShaderPrg_New(G, "spheredirect", spheredirect_vs, 0);
  }
  if (sphereShader){
    spheredirectShader->fid = sphereShader->fid;
    glAttachShader(spheredirectShader->id, spheredirectShader->fid);
    CShaderPrg_Link(spheredirectShader);
  }
  hasShaders &= (spheredirectShader!=0);
  CShaderMgr_AddShaderPrg(G->ShaderMgr, spheredirectShader);
  

  cylinderShader = CShaderPrg_NewFromFile(G, "cylinder", "cylinder.vs", "cylinder.fs");
  if (!cylinderShader){
    PRINTFB(G, FB_ShaderMgr, FB_Warnings)
      " PyMOLShader_NewFromFile-Warning: cylinder shader files not found, loading from memory.\n" ENDFB(G);
    cylinderShader = CShaderPrg_New(G, "cylinder", cylinder_vs, cylinder_fs);
  }
  hasCylinderShader &= (cylinderShader!=0);
  CShaderMgr_AddShaderPrg(G->ShaderMgr, cylinderShader);

  cylinderNoFFShader = CShaderPrg_NewFromFile(G, "cylinder_no_ff", 0, "cylinder_no_ff.fs");
  if (!cylinderNoFFShader){
    PRINTFB(G, FB_ShaderMgr, FB_Warnings)
      " PyMOLShader_NewFromFile-Warning: cylinder_no_ff shader files not found, loading from memory.\n" ENDFB(G);
    cylinderNoFFShader = CShaderPrg_New(G, "cylinder_no_ff", 0, cylinder_no_ff_fs);
  }
  if (cylinderShader){
    cylinderNoFFShader->vid = cylinderShader->vid;
    glAttachShader(cylinderNoFFShader->id, cylinderNoFFShader->vid);
    CShaderPrg_Link(cylinderNoFFShader);
  }
  hasShaders &= (cylinderNoFFShader!=0);
  CShaderMgr_AddShaderPrg(G->ShaderMgr, cylinderNoFFShader);

  /* report GLSL version */
  if (G && G->Option && !G->Option->quiet) {
    getGLSLVersion(G, &major, &minor);
    sprintf(buf, " Detected GLSL version %d.%d.\n", major, minor);
    FeedbackAdd(G, buf);
  }
#endif

  G->ShaderMgr->ShadersPresent = hasShaders;

#ifdef PURE_OPENGL_ES_2
  CShaderPrg_Enable_DefaultShader(G);
#endif
  if (hasShaders) {
    SettingSetGlobal_b(G, cSetting_use_shaders, 1);
  } else {
      disableShaders(G);
  }
}

/* getGLVersion -- determine user's GL version
 * PARAMS
 * major, return value for major
 * minor, return value for minor
 *
 * RETURNS
 * nothing; writes to major and minor
 */
void getGLVersion(PyMOLGlobals * G, int *major, int* minor) {
  /* query the version string */
  const char* verstr = (const char*) glGetString(GL_VERSION);
  /* attempt to store the values into major and minor */
  if ((verstr==NULL) || (sscanf(verstr,"%d.%d", major, minor) != 2)) {
    *major = *minor = 0;
    /* Use PyMOL FB system, instead of fprintf */
    PRINTFD(G, FB_ObjectVolume) 
      "Invalid GL_VERSION format.\n" ENDFD;
  }
}


/* getGLSLVersion -- determine user's GLSL version
 * PARAMS
 * major, rval for major
 * minor, rval for minor
 */
#ifndef _PYMOL_PURE_OPENGL_ES
void getGLSLVersion(PyMOLGlobals * G, int* major, int* minor) {
  int gl_major, gl_minor;
  *major = *minor = 0;

  /* grab the GL version */
  getGLVersion(G, &gl_major, &gl_minor);

  /* GL version 1 */
  if (1==gl_major) {
    const char* extstr = (const char*) glGetString(GL_EXTENSIONS);
    if ((extstr!=NULL)  &&
        (strstr(extstr, "GL_ARB_shading_language_100")!=NULL)){
      *major = 1;
      *minor = 0;
    }
  }
  /* GL > version 1 */
  else if (gl_major>=2) {
    const char* verstr = (const char*) glGetString(GL_SHADING_LANGUAGE_VERSION);

    if ((verstr==NULL) || (sscanf(verstr, "%d.%d", major, minor)!=2)){
      *major = *minor = 0;

      if (G && G->Option && !G->Option->quiet) {
        PRINTFD(G, FB_ObjectVolume) 
          "Invalid GL_SHADING_LANGUAGE_VERSION format.\n" ENDFD;
      }
    }
  }
}
#endif



/* ============================================================================
 * CShaderMgr -- Simple Shader Manager class
 */
CShaderMgr * CShaderMgr_New(PyMOLGlobals * G)
{
  /* init/alloc the new ShaderMgr, now called 'I' */
  OOAlloc(G, CShaderMgr);

  if (!G) { /* error out */
    return NULL;
  }
   
  if (!I) {
    /* error out */
    if (G && G->Option && !G->Option->quiet) {
    PRINTFB(G, FB_ShaderMgr, FB_Errors)
      " CShaderMgr_New-Error: Failed to create the shader manager.  Shader disabled.\n" ENDFB(G)
    }
    return NULL;
  }

  I->G = G;
  I->current_shader = 0;
  DListInit(I->programs, prev, next, CShaderPrg);
  I->ShadersPresent = 0;
  I->vbos_to_free = 0;
  I->number_of_vbos_to_free = 0;
  I->stereo_flag = 0;
  return I;
}

void CShaderMgr_Delete(CShaderMgr * I)
{
  CShaderPrg * ptr, *target;
  if (!I) { /* error out */
    return;
  }

  if (I->programs) {
    ptr = I->programs;
    while (ptr != I->programs) {
      target = ptr;
      ptr = ptr->next;
      DListRemove(target, prev, next);
      DListElemFree(target);
      target=NULL;
    }
  }
  OOFreeP(I);
}

int CShaderMgr_AddShaderPrg(CShaderMgr * I, CShaderPrg * s)
{
  if (!I || !s)
    return 0;

  DListInsert(I->programs, s, prev, next);
  return 1;
}

int CShaderMgr_RemoveShaderPrg(CShaderMgr * I, const char * name)
{
  CShaderPrg * p = NULL;
  DListIterate(I->programs, p, next) 
    {
      if (p && strcmp(p->name,name)==0) break;
    }
  DListRemove(p, prev, next);
  return 1;
}

CShaderPrg * CShaderMgr_GetShaderPrg(CShaderMgr * I, const char * name)
{
  CShaderPrg * p = NULL, *ret = NULL;
  DListIterate(I->programs, p, next) 
    {
      if (p && strcmp(p->name,name)==0){
	ret = p;
	break;
      }
    }

  I->current_shader = ret;
  return ret;
}

int CShaderMgr_ShadersPresent(CShaderMgr * I)
{
  return I->ShadersPresent;
}

char * CShaderMgr_ReadShaderFromDisk(PyMOLGlobals * G, const char * fileName) {
  FILE* f;
  long size;
  char* buffer = NULL, *p, *pymol_path, *shader_path, *fullFile;
  size_t res;

  PRINTFB(G, FB_ShaderMgr, FB_Debugging)
    "CShaderMgr_ReadShaderFromDisk: fileName='%s'\n", fileName
    ENDFB(G);
  /* check the input */
  if (!strlen(fileName)) {
    PRINTFB(G, FB_ShaderMgr, FB_Errors)
      " PyMOLShader_NewFromFile-Error: empty filename, cannot create shader. " ENDFB(G);
    return NULL;
  }
  
  pymol_path = getenv("PYMOL_PATH");
  if (!pymol_path){
    PRINTFB(G, FB_ShaderMgr, FB_Warnings)
      " PyMOLShader_NewFromFile-Warning: PYMOL_PATH not set, cannot read shader config files from disk\n" ENDFB(G);
    return NULL;
  }
  /* make this a setting */
  shader_path = "/data/shaders/";
  fullFile = malloc( sizeof(char) * (strlen(pymol_path)+strlen(shader_path)+strlen(fileName)+1));
  fullFile = strcpy(fullFile, pymol_path);
  fullFile = strcat(fullFile, shader_path);
  fullFile = strcat(fullFile, fileName);

  /* read the file from disk */
  f = fopen(fullFile, "rb");

  if (!f) {
    PRINTFB(G, FB_ShaderMgr, FB_Errors)
      " PyMOLShader_NewFromFile-Error: Unable to open file '%s' PYMOL_PATH='%s'\n", fullFile, pymol_path ENDFB(G);
    return NULL;
  } else {
    PRINTFB(G, FB_ShaderMgr, FB_Blather)
      " PyMOLShader_NewFromFile: Loading shader from '%s'.\n", fullFile ENDFB(G);
  }
  
  /* rewind and get file size */
  fseek(f, 0, SEEK_END);
  size = ftell(f);
  fseek(f, 0, SEEK_SET);

  buffer = (char*) mmalloc(size+255);
  ErrChkPtr(G,buffer);
  p = buffer;
  fseek(f, 0, SEEK_SET);
  res = fread(p, size, 1, f);
  /* error reading shader */
  if(1!=res) {
    PRINTFB(G, FB_ShaderMgr, FB_Errors)
      " PyMOLShader_NewFromFile-Error: size(%ld)!=res(%ld)\n", size,res ENDFB(G);
    return NULL;
  }

  p[size] = 0;
  fclose(f);

  free(fullFile);
  return buffer;
}


#ifdef _PYMOL_OPENGL_SHADERS

#ifndef GL_FRAGMENT_PROGRAM_ARB
#define GL_FRAGMENT_PROGRAM_ARB                         0x8804
#endif

static GLboolean ProgramStringIsNative(PyMOLGlobals * G,
                                       GLenum target, GLenum format,
                                       GLsizei len, const GLvoid * string)
{
  GLint errorPos, isNative;
  glProgramStringARB(target, format, len, string);
  glGetIntegerv(GL_PROGRAM_ERROR_POSITION_ARB, &errorPos);
  glGetProgramivARB(target, GL_PROGRAM_UNDER_NATIVE_LIMITS_ARB, &isNative);
  if((errorPos == -1) && (isNative == 1))
    return GL_TRUE;
  else if(errorPos >= 0) {
    if(Feedback(G, FB_OpenGL, FB_Errors)) {
      printf("OpenGL-Error: ARB shader error at char %d\n---->%s\n", errorPos,
             ((char *) string) + errorPos);
    }
  }
  return GL_FALSE;
}
#endif

/* ============================================================================
 * CShaderPrg -- Simple Shader class
 */
CShaderPrg * CShaderPrg_NewARB(PyMOLGlobals * G, const char * name, const char * v, const char * f)
{
  /* if v == f == NULL, read 'name.vs' and 'name.fs' from disk */

  /* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
#ifdef WIN32
  if(!(glGenProgramsARB && glBindProgramARB &&
       glDeleteProgramsARB && glProgramStringARB && glProgramEnvParameter4fARB)) {
    glGenProgramsARB = (PFNGLGENPROGRAMSARBPROC) wglGetProcAddress("glGenProgramsARB");
    glBindProgramARB = (PFNGLBINDPROGRAMARBPROC) wglGetProcAddress("glBindProgramARB");
    glDeleteProgramsARB =
      (PFNGLDELETEPROGRAMSARBPROC) wglGetProcAddress("glDeleteProgramsARB");
    glProgramStringARB =
      (PFNGLPROGRAMSTRINGARBPROC) wglGetProcAddress("glProgramStringARB");
    glProgramEnvParameter4fARB =
      (PFNGLPROGRAMENVPARAMETER4FARBPROC) wglGetProcAddress("glProgramEnvParameter4fARB");
    glGetProgramivARB = (PFNGLGETPROGRAMIVARBPROC) wglGetProcAddress("glGetProgramivARB");
  }

  if(glGenProgramsARB && glBindProgramARB &&
     glDeleteProgramsARB && glProgramStringARB && glProgramEnvParameter4fARB)
#endif
  {
    /* END PROPRIETARY CODE SEGMENT */

#ifdef _PYMOL_OPENGL_SHADERS
    int ok = true;
    GLuint programs[2];
    glGenProgramsARB(2, programs);
    
    /* load the vertex program */
    glBindProgramARB(GL_VERTEX_PROGRAM_ARB, programs[0]);
    
    ok = ok && (ProgramStringIsNative(G, GL_VERTEX_PROGRAM_ARB,
				      GL_PROGRAM_FORMAT_ASCII_ARB, strlen(v), v));
    
    if(Feedback(G, FB_OpenGL, FB_Debugging))
      PyMOLCheckOpenGLErr("loading vertex program");
    
    /* load the fragment program */
    glBindProgramARB(GL_FRAGMENT_PROGRAM_ARB, programs[1]);
    
    ok = ok && (ProgramStringIsNative(G, GL_FRAGMENT_PROGRAM_ARB,
				      GL_PROGRAM_FORMAT_ASCII_ARB, strlen(f), f));
    
    if(Feedback(G, FB_OpenGL, FB_Debugging))
      PyMOLCheckOpenGLErr("loading fragment program");
    if(ok) {
      CShaderPrg * I = NULL;
      
      DListElemAlloc(G, I, CShaderPrg);
      DListElemInit(I, prev, next);

      I->G = G;
      I->name = strdup(name);

      I->vid = programs[0];
      I->fid = programs[1];

      CShaderMgr_AddShaderPrg(G->ShaderMgr, I);

      return I;
    } else {
      glDeleteProgramsARB(2, programs);
    }
#endif
  }
  return NULL;
}

/* ============================================================================
 * CShaderPrg -- Simple Shader class
 */
CShaderPrg * CShaderPrg_New(PyMOLGlobals * G, const char * name, const char * v, const char * f)
{
  int status, howLong;
  char infoLog[MAX_LOG_LEN];
  /* if v == f == NULL, read 'name.vs' and 'name.fs' from disk */
  CShaderPrg * I = NULL;
  DListElemAlloc(G, I, CShaderPrg);
  DListElemInit(I, prev, next);

  I->G = G;
  I->name = strdup(name);

  I->id = glCreateProgram();
  PRINTFB(G, FB_ShaderMgr, FB_Debugging)
    "Created program with id: %d\n", I->id ENDFB(G);

  /* VERTEX shader setup */
  /* CShaderPrg_InitShader(I, GL_VERTEX_SHADER); */
  if (v){
    I->v = strdup(v);
    I->vid = glCreateShader(GL_VERTEX_SHADER);
    PRINTFB(G, FB_ShaderMgr, FB_Debugging)
      "Created vertex shader with id: %d\n", I->vid ENDFB(G);
    glShaderSource(I->vid, 1, (const GLchar**) &I->v, NULL);
    glCompileShader((GLuint) I->vid);
    /* verify compilation */
    glGetShaderiv(I->vid, GL_COMPILE_STATUS, &status);
    if (!status) {
      if (G && G->Option && !G->Option->quiet) {
	PRINTFB(G, FB_ShaderMgr, FB_Errors) " CShaderPrg_New-Error: vertex shader compilation failed name='%s'; log follows.\n", I->name ENDFB(G);
	glGetShaderInfoLog(I->vid, MAX_LOG_LEN, &howLong, infoLog);
	PRINTFB(G, FB_ShaderMgr, FB_Errors)
	  "infoLog=%s\n", infoLog ENDFB(G);
      }
      return 0;
    }
    PRINTFB(G, FB_ShaderMgr, FB_Debugging)
      "CShaderPrg_New-Message: vertex shader compiled.\n" ENDFB(G);
    glAttachShader(I->id, I->vid);
  }

  if (f){
    /* FRAGMENT source setup */
    /* CShaderPrg_InitShader(I, GL_FRAGMENT_SHADER); */
    I->f = strdup(f);
    I->fid = glCreateShader(GL_FRAGMENT_SHADER);
    PRINTFB(G, FB_ShaderMgr, FB_Debugging)
      "Created fragment shader with id: %d\n", I->fid 
      ENDFB(G);
    
    glShaderSource(I->fid, 1, (const GLchar **) &I->f, NULL);
    glCompileShader((GLuint) I->fid);
    /* verify compilation */
    glGetShaderiv(I->fid, GL_COMPILE_STATUS, &status);
    if (!status) {
      if (G && G->Option && !G->Option->quiet) {
	PRINTFB(G, FB_ShaderMgr, FB_Errors)
	  " CShaderPrg-Error: fragment shader compilation failed name='%s'; log follows.\n", I->name
	  ENDFB(G);
	glGetShaderInfoLog(I->fid, MAX_LOG_LEN, &howLong, infoLog);
	PRINTFB(G, FB_ShaderMgr, FB_Errors)
	  "infoLog=%s\n", infoLog ENDFB(G);
      }
      return 0;
    }    
    glAttachShader(I->id, I->fid);
  }

  if (v && f){
    /* Link the new program */
    if (!CShaderPrg_Link(I)){
      CShaderPrg_Delete(I);
      return NULL;
    }
  }
  return I;
}


CShaderPrg * CShaderPrg_NewFromFile(PyMOLGlobals * G, const char * name, const char * vFile, const char * fFile)
{
  char *vFileStr = NULL, *fFileStr = NULL;
  if (vFile){
    vFileStr = CShaderMgr_ReadShaderFromDisk(G, vFile);
    if (!vFileStr){
      return (NULL);
    }
  }
  if (fFile){
    fFileStr = CShaderMgr_ReadShaderFromDisk(G, fFile);
    if (!fFileStr){
      return (NULL);
    }
  }
  return CShaderPrg_New(G, name, 
			vFileStr, 
			fFileStr);
}


void CShaderPrg_Delete(CShaderPrg * I)
{
  glDeleteShader(I->vid);
  glDeleteShader(I->fid);
  glDeleteProgram(I->id);
  OOFreeP(I->v);
  OOFreeP(I->f);
  OOFreeP(I->name);
  I->next = I->prev = NULL;
  DListElemFree(I);
}


int CShaderPrg_Enable(CShaderPrg * I)
{
  int howLong, ok;
  char infoLog[MAX_LOG_LEN];
  PyMOLGlobals * G = I->G;

  if (!I) return 0;

  /* linked? */
  ok = CShaderPrg_IsLinked(I);
  /* no, so give it a shot */
  if (!ok) {
    ok = CShaderPrg_Link(I);    
  }
  /* did that work? */
  if (!ok) {
    if (G && G->Option && !G->Option->quiet) {
      PRINTFB(G, FB_ShaderMgr, FB_Errors)
	"CShaderPrg_Enable-Error: Cannot enable the shader program; linking failed.  Shaders disabled.  Log follows.\n"
	ENDFB(G);
      glGetProgramInfoLog(I->id, MAX_LOG_LEN, &howLong, infoLog);
      PRINTFB(G, FB_ShaderMgr, FB_Errors)
	"%s\n", infoLog
	ENDFB(G);
    }
    return 0;
  }
  /* if so, use the program */
  glUseProgram(I->id);
  return 1;
}

int CShaderPrg_Disable(CShaderPrg * p)
{
  glUseProgram(0);
  return 1;
}

int CShaderPrg_DisableARB(CShaderPrg * p)
{
  glDisable(GL_FRAGMENT_PROGRAM_ARB);
#ifndef _PYMOL_PURE_OPENGL_ES
  glDisable(GL_VERTEX_PROGRAM_ARB);
#endif
  return 1;
}

 int CShaderPrg_IsLinked(CShaderPrg * I )
 {
   int status;
   glGetProgramiv(I->id, GL_LINK_STATUS, &status);
   return status==GL_TRUE;
 }

int CShaderPrg_Link(CShaderPrg * I)
{
  int howLong;
  char infoLog[MAX_LOG_LEN];
  PyMOLGlobals * G = I->G;

  glLinkProgram(I->id);
  
  if (!CShaderPrg_IsLinked(I)) {
    if (G && G->Option && !G->Option->quiet) {
      GLint maxVarFloats[1];
#ifndef _PYMOL_PURE_OPENGL_ES
      glGetIntegerv(GL_MAX_VARYING_FLOATS, maxVarFloats);
#endif
    PRINTFB(G, FB_ShaderMgr, FB_Errors)
	" CShaderPrg_Link-Error: Shader program failed to link name='%s'; GL_MAX_VARYING_FLOATS=%d log follows.\n", I->name, maxVarFloats[0]
	ENDFB(G);
      glGetProgramInfoLog(I->id, MAX_LOG_LEN, &howLong, infoLog);
      PRINTFB(G, FB_ShaderMgr, FB_Errors)
	"%s\n", infoLog
	ENDFB(G);
    }
    return 0;
  }
  return 1;      
}


/* accessors/mutators/uniform setters */
int CShaderPrg_Set1i(CShaderPrg * p, const char * name, int i)
{
  if (p && p->id) {
    GLint loc = glGetUniformLocation(p->id, name);
    if (loc < 0)
      return 0;
    glUniform1i(loc, i);
  }
  return 1;
}

int CShaderPrg_Set3f(CShaderPrg * p, const char * name, float f1, float f2, float f3)
{
  if (p && p->id) {
    GLint loc = glGetUniformLocation(p->id, name);
    if (loc < 0)
      return 0;
    glUniform3f(loc, f1, f2, f3);
  }
  return 1;
}

int CShaderPrg_SetMat3f(CShaderPrg * p, const char * name, float* m) {
  if (p && p->id) {
    GLint loc = glGetUniformLocation(p->id, name);
    if (loc < 0)
      return 0;
    glUniformMatrix3fv(loc, 1, GL_TRUE, m);
  }
  return 1;
}

int CShaderPrg_Set4f(CShaderPrg * p, const char * name, float f1, float f2, float f3, float f4)
{
  if (p && p->id) {
    GLint loc = glGetUniformLocation(p->id, name);
    if (loc < 0)
      return 0;
    glUniform4f(loc, f1, f2, f3, f4);
  }
  return 1;
}

int CShaderPrg_Set1f(CShaderPrg * p, const char * name, float f)
{
  if (p && p->id) {
    GLint loc = glGetUniformLocation(p->id, name);
    if (loc < 0)
      return 0;
    glUniform1f(loc, f);
  }
  return 1;
}

int CShaderPrg_GetAttribLocation(CShaderPrg * p, const char * name)
{
  GLint loc = -1;

  if (p && p->id) {
    loc = glGetAttribLocation(p->id, name);
    if (loc < 0)
      return -1;
  }
  return loc;
}

void CShaderPrg_SetAttrib4fLocation(CShaderPrg * p, const char * name, float f1, float f2, float f3, float f4){
  if (p){
    int attr = CShaderPrg_GetAttribLocation(p, name);
    if (attr>=0){
      glVertexAttrib4f(attr, f1, f2, f3, f4);
    }
  }
}

void CShaderMgr_FreeAllVBOs(CShaderMgr * I){
  GLuint *vboids = I->vbos_to_free, nvbos = I->number_of_vbos_to_free;
  I->vbos_to_free = 0;
  I->number_of_vbos_to_free = 0;
  if (I && vboids){
    GLuint i, nvbo=0 ;
    for (i=0; i<nvbos; i++){
      if (glIsBuffer(vboids[i])){
	vboids[nvbo++] = vboids[i];
      } else {
	PRINTFB(I->G, FB_ShaderMgr, FB_Warnings) "WARNING: CShaderMgr_FreeAllVBOs() buffer is not a VBO i=%d vboids[i]=%d\n", i, vboids[i] ENDFB(I->G);
      }
    }
    if (nvbo){
      glDeleteBuffers(nvbo, vboids);
    }
    VLAFree(vboids);
  }
}

#define STEP_FOR_BUF 100

void CShaderMgr_AddVBOToFree(CShaderMgr * I, GLuint vboid){
  if (I && I->vbos_to_free){
    int nvbostofree = I->number_of_vbos_to_free++;
    int exp = STEP_FOR_BUF*(1+((nvbostofree+1)/STEP_FOR_BUF));
    VLACheck(I->vbos_to_free, GLuint, exp);
    I->vbos_to_free[nvbostofree] = vboid;
  } else {
    I->vbos_to_free = VLAlloc(GLuint, STEP_FOR_BUF);
    I->vbos_to_free[0] = vboid;
    I->number_of_vbos_to_free = 1;
  }
}

void CShaderPrg_Set_Specular_Values(PyMOLGlobals * G, CShaderPrg * shaderPrg){
  float spec_value = SettingGet(G, cSetting_specular);
  float settingSpecReflect, settingSpecDirect, settingSpecDirectPower, settingSpecPower;
  int spec_count = SettingGet(G, cSetting_spec_count);

  settingSpecPower = SettingGet(G, cSetting_spec_power);

  if(settingSpecPower < 0.0F) {
    settingSpecPower = SettingGet(G, cSetting_shininess);
  }

  CShaderPrg_Set1f(shaderPrg, "shininess", settingSpecPower);

  if (spec_count < 0){
    spec_count = SettingGet(G, cSetting_light_count);
  }
  if(spec_value == 1.0F)
    spec_value = SettingGet(G, cSetting_specular_intensity);

  settingSpecReflect = SettingGet(G, cSetting_spec_reflect);
  settingSpecReflect = SceneGetSpecularValue(G, settingSpecReflect, 10);
  settingSpecDirect = SettingGet(G, cSetting_spec_direct);
  settingSpecDirectPower = SettingGet(G, cSetting_spec_direct_power);

  if(settingSpecReflect < 0.0F)
    settingSpecReflect = spec_value;
  if(settingSpecDirect < 0.0F)
    settingSpecDirect = spec_value;
  if(settingSpecDirectPower < 0.0F)
    settingSpecDirectPower = settingSpecPower;


  if(settingSpecReflect > 1.0F)
    settingSpecReflect = 1.0F;
  if(SettingGet(G, cSetting_specular) < R_SMALL4) {
    settingSpecReflect = 0.0F;
  }
  CShaderPrg_Set1f(shaderPrg, "spec_value_0", settingSpecDirect);
  CShaderPrg_Set1f(shaderPrg, "shininess_0", settingSpecDirectPower);
  CShaderPrg_Set1f(shaderPrg, "spec_value", settingSpecReflect);
  CShaderPrg_Set1i(shaderPrg, "spec_count", spec_count);

}

void CShaderPrg_Set_AnaglyphMode(PyMOLGlobals * G, CShaderPrg * shaderPrg, int mode) {
  extern float anaglyphR_constants[6][9];
  extern float anaglyphL_constants[6][9];
  /** Coefficients from: http://3dtv.at/Knowhow/AnaglyphComparison_en.aspx */
  /** anaglyph[R|L]_constants are found in Scene.c b/c of ray tracing */
  CShaderPrg_SetMat3f(shaderPrg, "matR", anaglyphR_constants[mode]);
  CShaderPrg_SetMat3f(shaderPrg, "matL", anaglyphL_constants[mode]);
  CShaderPrg_Set1f(shaderPrg, "gamma", SettingGet(G, cSetting_gamma));
}

CShaderPrg *CShaderPrg_Enable_DefaultShader(PyMOLGlobals * G){
  float fog_enabled, *fog_color_top, *fog_color_bottom;
  int bg_gradient, stereo, stereo_mode;
  CShaderPrg * shaderPrg = CShaderMgr_GetShaderPrg(G->ShaderMgr, "default");
  CShaderPrg_Enable(shaderPrg);
  fog_enabled = SettingGet(G, cSetting_depth_cue) ? 1.0 : 0.0;
  bg_gradient = SettingGet(G, cSetting_bg_gradient);
  if (bg_gradient){
    fog_color_top = SettingGetfv(G, cSetting_bg_rgb_top);
    fog_color_bottom = SettingGetfv(G, cSetting_bg_rgb_bottom);
  } else {
    fog_color_top = SettingGetfv(G, cSetting_bg_rgb);
    fog_color_bottom = fog_color_top;
  }

  stereo = SettingGetGlobal_i(G, cSetting_stereo);
  stereo_mode = SettingGetGlobal_i(G, cSetting_stereo_mode);
  CShaderPrg_Set1i(shaderPrg, "stereo_flag", G->ShaderMgr->stereo_flag);
  if (stereo && stereo_mode==cStereo_anaglyph)
    CShaderPrg_Set_AnaglyphMode(G, shaderPrg, SettingGetGlobal_i(G, cSetting_anaglyph_mode));

  CShaderPrg_Set1i(shaderPrg, "bg_gradient", bg_gradient);
  CShaderPrg_Set3f(shaderPrg, "fog_color_top", fog_color_top[0], fog_color_top[1], fog_color_top[2]);
  CShaderPrg_Set3f(shaderPrg, "fog_color_bottom", fog_color_bottom[0], fog_color_bottom[1], fog_color_bottom[2]);
  CShaderPrg_Set1f(shaderPrg, "fog_enabled", fog_enabled);
  CShaderPrg_Set1i(shaderPrg, "lighting_enabled", 1); // lighting on by default
  CShaderPrg_Set1i(shaderPrg, "two_sided_lighting_enabled", SceneGetTwoSidedLighting(G));
  CShaderPrg_Set1i(shaderPrg, "light_count", SettingGet(G, cSetting_light_count));
  CShaderPrg_Set1f(shaderPrg, "ambient_occlusion_scale", 0.f);
  CShaderPrg_Set1i(shaderPrg, "accessibility_mode", SettingGetGlobal_i(G, cSetting_ambient_occlusion_mode) / 4);

  {
    int interior_color = SettingGetGlobal_i(G, cSetting_ray_interior_color);
    float *color, inter[] = { 0.f, 0.f, 0.f }, threshold = 0.f;
    if (interior_color < 0){
      threshold = .22f;  // this is hardcoded for now, need to figure out exactly what Ray.c does
                         // to cull the backfacing polygons
    }
    CShaderPrg_Set1f(shaderPrg, "interior_color_threshold", threshold);
    if (interior_color < 0){
      color = inter;
    } else {
      ColorGetEncoded(G, interior_color, inter);
      color = inter;
    }
    CShaderPrg_Set4f(shaderPrg, "interior_color", color[0], color[1], color[2], 1.f);
  }
  CShaderPrg_Set1i(shaderPrg, "use_interior_color_threshold", 0);

  CShaderPrg_Set_Specular_Values(G, shaderPrg);

  return (shaderPrg);
}
CShaderPrg *CShaderPrg_Enable_CylinderShader(PyMOLGlobals * G){
  int fog_enabled, bg_gradient;
  float *fog_color_top, *fog_color_bottom;
  int ortho, stereo, stereo_mode;
  int width, height;
  CShaderPrg *shaderPrg;
  float *m;

  SceneGetWidthHeight(G, &width, &height);
  m = SceneGetMatrix(G);
  if (SettingGetGlobal_i(G, cSetting_cylinder_shader_ff_workaround)){
    shaderPrg = CShaderMgr_GetShaderPrg(G->ShaderMgr, "cylinder_no_ff");
  } else {
    shaderPrg = CShaderMgr_GetShaderPrg(G->ShaderMgr, "cylinder");
  }
  CShaderPrg_Enable(shaderPrg);
  CShaderPrg_Set1f(shaderPrg, "uni_radius", 0.f);
  fog_enabled = SettingGet(G, cSetting_depth_cue) ? 1.0 : 0.0;
  bg_gradient = SettingGet(G, cSetting_bg_gradient);
  if (bg_gradient){
    fog_color_top = SettingGetfv(G, cSetting_bg_rgb_top);
    fog_color_bottom = SettingGetfv(G, cSetting_bg_rgb_bottom);
  } else {
    fog_color_top = SettingGetfv(G, cSetting_bg_rgb);
    fog_color_bottom = fog_color_top;
  }
  stereo = SettingGetGlobal_i(G, cSetting_stereo);
  stereo_mode = SettingGetGlobal_i(G, cSetting_stereo_mode);
  CShaderPrg_Set1i(shaderPrg, "stereo_flag", G->ShaderMgr->stereo_flag);
  if (stereo && stereo_mode==cStereo_anaglyph)
    CShaderPrg_Set_AnaglyphMode(G, shaderPrg, SettingGetGlobal_i(G, cSetting_anaglyph_mode));

  CShaderPrg_Set1i(shaderPrg, "bg_gradient", bg_gradient);
  CShaderPrg_Set3f(shaderPrg, "fog_color_top", fog_color_top[0], fog_color_top[1], fog_color_top[2]);
  CShaderPrg_Set3f(shaderPrg, "fog_color_bottom", fog_color_bottom[0], fog_color_bottom[1], fog_color_bottom[2]);
  CShaderPrg_Set1f(shaderPrg, "fog_enabled", fog_enabled);
  CShaderPrg_Set1f(shaderPrg, "inv_height", 1.0/height);
  ortho = SettingGet(G, cSetting_ortho);
  CShaderPrg_Set1f(shaderPrg, "ortho", ortho ? 1.0 : 0.0);
  CShaderPrg_Set1f(shaderPrg, "flat_caps", 0.0);
  CShaderPrg_Set1i(shaderPrg, "filter_front_facing", SettingGet(G, cSetting_cylinders_shader_filter_faces));
  CShaderPrg_Set1i(shaderPrg, "two_sided_lighting_enabled", SceneGetTwoSidedLighting(G));
  CShaderPrg_Set1i(shaderPrg, "light_count", SettingGet(G, cSetting_light_count));
  CShaderPrg_Set1i(shaderPrg, "filter_front_facing", SettingGet(G, cSetting_cylinders_shader_filter_faces));
  {
    float smooth_half_bonds = (SettingGetGlobal_i(G, cSetting_smooth_half_bonds)) ? .2f : 0.f;
    CShaderPrg_Set1f(shaderPrg, "half_bond", smooth_half_bonds);
  }
  CShaderPrg_Set_Specular_Values(G, shaderPrg);
  return shaderPrg;
}

CShaderPrg *CShaderPrg_Enable_SphereShader(PyMOLGlobals * G, char *name){
  int fog_enabled, bg_gradient;
  float *fog_color_top, *fog_color_bottom;
  int ortho, stereo, stereo_mode;
  CShaderPrg *shaderPrg;
  int width, height;
  SceneGetWidthHeight(G, &width, &height);
  shaderPrg = CShaderMgr_GetShaderPrg(G->ShaderMgr, name);
  CShaderPrg_Enable(shaderPrg);
  CShaderPrg_Set1i(shaderPrg, "lighting_enabled", 1);
  CShaderPrg_Set1f(shaderPrg, "sphere_size_scale", 1.f);
  fog_enabled = SettingGet(G, cSetting_depth_cue) ? 1.0 : 0.0;
  bg_gradient = SettingGet(G, cSetting_bg_gradient);
  if (bg_gradient){
    fog_color_top = SettingGetfv(G, cSetting_bg_rgb_top);
    fog_color_bottom = SettingGetfv(G, cSetting_bg_rgb_bottom);
  } else {
    fog_color_top = SettingGetfv(G, cSetting_bg_rgb);
    fog_color_bottom = fog_color_top;
  }

  stereo = SettingGetGlobal_i(G, cSetting_stereo);
  stereo_mode = SettingGetGlobal_i(G, cSetting_stereo_mode);
  CShaderPrg_Set1i(shaderPrg, "stereo_flag", G->ShaderMgr->stereo_flag);
  if (stereo && stereo_mode==cStereo_anaglyph)
    CShaderPrg_Set_AnaglyphMode(G, shaderPrg, SettingGetGlobal_i(G, cSetting_anaglyph_mode));

  CShaderPrg_Set1i(shaderPrg, "bg_gradient", bg_gradient);
  CShaderPrg_Set3f(shaderPrg, "fog_color_top", fog_color_top[0], fog_color_top[1], fog_color_top[2]);
  CShaderPrg_Set3f(shaderPrg, "fog_color_bottom", fog_color_bottom[0], fog_color_bottom[1], fog_color_bottom[2]);
  CShaderPrg_Set1f(shaderPrg, "fog_enabled", fog_enabled);
  CShaderPrg_Set1f(shaderPrg, "inv_height", 1.0/height);
  ortho = SettingGet(G, cSetting_ortho);
  CShaderPrg_Set1f(shaderPrg, "ortho", ortho ? 1.0 : 0.0);
  CShaderPrg_Set1i(shaderPrg, "light_count", SettingGet(G, cSetting_light_count));

  {
    float adj;
    float fov = SettingGetGlobal_f(G, cSetting_field_of_view);
    /* Polynomial fitting for adjustment values relative to the field of view */
    if (fov <= 90.f){
      adj = 1.0027+0.000111*fov+0.000098*fov*fov;
    } else {
      adj = 2.02082 - 0.033935*fov + 0.00037854*fov*fov;
    }
    CShaderPrg_Set1f(shaderPrg, "horizontal_adjustment", adj);
    CShaderPrg_Set1f(shaderPrg, "vertical_adjustment", adj);
  }
  CShaderPrg_Set_Specular_Values(G, shaderPrg);
  return (shaderPrg);
}

CShaderPrg *CShaderPrg_Enable_SphereShaderARB(PyMOLGlobals * G){
  CShaderPrg *shaderPrg = NULL;
#ifdef _PYMOL_OPENGL_SHADERS
  /* load the vertex program */
  shaderPrg = CShaderMgr_GetShaderPrg(G->ShaderMgr, "sphere_arb");
  glBindProgramARB(GL_VERTEX_PROGRAM_ARB, shaderPrg->vid);
  
  /* load the fragment program */
  glBindProgramARB(GL_FRAGMENT_PROGRAM_ARB, shaderPrg->fid);
  
  /* load some safe initial values  */
  
  glProgramEnvParameter4fARB(GL_VERTEX_PROGRAM_ARB, 0, 0.0F, 0.0F, 1.0, 0.0F);
  glProgramEnvParameter4fARB(GL_FRAGMENT_PROGRAM_ARB, 0, 0.5F, 2.0F, 0.0F, 0.0F);
  
  glEnable(GL_VERTEX_PROGRAM_ARB);
  glEnable(GL_FRAGMENT_PROGRAM_ARB);
  
#endif
  return shaderPrg;
}
