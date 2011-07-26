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

/* TODO
 * (1) Debugging; Feeding back
 * (2) printfs needs to be FeedbackAdd or PRINTFD
 */

#include <string.h>
#include "ShaderMgr.h"
#include "OOMac.h"
#include "ListMacros.h"
#include "PyMOLOptions.h"
#include "os_gl.h"
#include "Feedback.h"

#define MAX_LOG_LEN 1024

const char* default_vs = "varying vec3 N, L0, H0, L1, H1;"
"varying vec4 D0, A0, D1, A1;"
"varying float fog;"
"void main()"
"{"
"  N = normalize(gl_NormalMatrix * gl_Normal);"
"  vec3 eye_pos = vec3(gl_ModelViewMatrix * gl_Vertex);"
"  vec3 aux = vec3(gl_LightSource[0].position - vec4(eye_pos, 1.0));"
"  L0 = normalize(aux);"
"  H0 = normalize(gl_LightSource[0].halfVector.xyz);"
"  aux = vec3(gl_LightSource[1].position - vec4(eye_pos, 1.0));"
"  L1 = normalize(aux);"
"  H1 = normalize(gl_LightSource[1].halfVector.xyz);"
"  A0 = gl_LightSource[0].ambient * gl_Color;"
"  D0 = gl_LightSource[0].diffuse * gl_Color;"
"  A1 = gl_LightSource[1].ambient * gl_Color;"
"  D1 = gl_LightSource[1].diffuse * gl_Color;"
"  gl_FogFragCoord = -eye_pos.z;"
"  fog = (gl_Fog.end - gl_FogFragCoord) * gl_Fog.scale;"
"  gl_Position = ftransform();"
"}";

const char* default_fs = "varying vec3 N, L0, H0, L1, H1;"
"varying vec4 D0, A0, D1, A1;"
""
"uniform float fog_r;"
"uniform float fog_g;"
"uniform float fog_b;"
"uniform float fog_enabled;"
""
"varying float fog;"
""
"void main()"
"{"
"  vec3 n, h;"
"  float NdotL, NdotH;"
"  vec4 color = A0 + A1;"
"  n = normalize(N);"
"  NdotL = max(dot(n, normalize(L0)), 0.0);"
"  float shininess = gl_FrontMaterial.shininess;"
"  if (NdotL > 0.0) {"
"      color += D0 * NdotL;"
"      h = normalize(H0);"
"      NdotH = max(dot(n, h), 0.0);"
"      color += gl_LightSource[0].specular * pow(NdotH, shininess);"
"  }"
"  NdotL = max(dot(n, normalize(L1)), 0.0);"
"  if (NdotL > 0.0) {"
"      color += D1 * NdotL;"
"      h = normalize(H1);"
"      NdotH = max(dot(n, h), 0.0);"
"      color += gl_LightSource[1].specular * pow(NdotH, shininess);"
"  }"
""
"  vec3 fog_color = vec3(fog_r, fog_g, fog_b);"
"  float cfog = mix(1.0, clamp(fog, 0.0, 1.0), fog_enabled);"
""
"  gl_FragColor = vec4(mix(fog_color, color.rgb, cfog), color.a);"
"}";

const char *volume_vs = "varying float fog;"
""
"void main()"
"{"
"  vec4 vertex = gl_ModelViewMatrix * gl_Vertex;"
"  gl_TexCoord[0] = gl_MultiTexCoord0;"
"  gl_ClipVertex = vertex;"
"  gl_Position = ftransform();"
"  gl_FogFragCoord = -vertex.z;"
"  fog = (gl_Fog.end - gl_FogFragCoord) * gl_Fog.scale;"
"}";
const char *volume_fs = "uniform sampler3D volumeTex;"
"uniform sampler1D colorTex;"
""
"uniform float fog_r;"
"uniform float fog_g;"
"uniform float fog_b;"
"uniform float fog_enabled;"
"varying float fog;"
""
"void main()"
"{"
"  vec3 fog_color = vec3(fog_r, fog_g, fog_b);"
"  vec4 color = texture1D(colorTex, texture3D(volumeTex, gl_TexCoord[0].xyz).x);"
"  if (color.a == 0.0) discard;"
"  float cfog = mix(1.0, clamp(fog, 0.0, 1.0), fog_enabled);"
"  gl_FragColor = vec4(vec3(mix(fog_color, color.rgb, cfog)), color.a);"
"}";




void getGLVersion(PyMOLGlobals * G, int *major, int* minor);
void getGLSLVersion(PyMOLGlobals * G, int* major, int* minor);

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
  CShaderPrg *p;
  int hasShaders = 0;
  int ok = 0;
  GLenum err; 

  ok = (G && G->HaveGUI); /* && G->ValidContext); */

  if (ok) {
    err = glewInit();
  } else {
    return;
  }

  if (GLEW_OK==err) {
    if (GLEW_VERSION_2_0) {
      FeedbackAdd(G, " Detected OpenGL version 2.0 or greater.  Shaders available.\n");
    }
    else { 
      FeedbackAdd(G, " Detected OpenGL version prior to 2.0.  Shaders and volumes unavailable.\n");
      return;
    }
  } 
  else {
    /* print info on glew error? */
    FeedbackAdd(G, " There was an error intializing GLEW.  Basic graphics, including\n shaders and volumes may be unavailable.\n");
    /*fprintf(stderr, " GLEW-Error: %s\n", glewGetErrorString(err));*/
    return;
  }

  /* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
  /* USE GLEW instead */
/* #ifdef __WIN32 */
#if 0
      printf("glActiveTexture is: %x\n", glActiveTexture);

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
    glGetProgramiv = (PFNGLGETPROGRAMIVPROC) wglGetProcAddress("glGetProgramiv");
    glAttachShader = (PFNGLATTACHSHADERPROC) wglGetProcAddress("glAttachShader");
    glGetShaderInfoLog = (PFNGLGETSHADERINFOLOGPROC) wglGetProcAddress("glGetShaderInfoLog");
    glGetShaderiv = (PFNGLGETSHADERIVPROC) wglGetProcAddress("glGetShaderiv");
    glCompileShader = (PFNGLCOMPILESHADERPROC) wglGetProcAddress("glCompileShader");
    glShaderSource = (PFNGLSHADERSOURCEPROC) wglGetProcAddress("glShaderSource");
    glCreateShader = (PFNGLCREATESHADERPROC) wglGetProcAddress("glCreateShader");
    glCreateProgram = (PFNGLCREATEPROGRAMPROC) wglGetProcAddress("glCreateProgram");
    glDeleteProgram  = (PFNGLDELETEPROGRAMPROC) wglGetProcAddress("glDeleteProgram");
    glDeleteShader = (PFNGLDELETESHADERPROC) wglGetProcAddress("glDeleteShader");
    glUseProgram = (PFNGLUSEPROGRAMPROC) wglGetProcAddress("glUseProgram");
    glGetProgramInfoLog = (PFNGLGETPROGRAMINFOLOGPROC) wglGetProcAddress("glGetProgramInfoLog");
    glLinkProgram = (PFNGLLINKPROGRAMPROC) wglGetProcAddress("glLinkProgram");
    glUniform1i = (PFNGLUNIFORM1IPROC) wglGetProcAddress("glUniform1i");
    glGetUniformLocation = (PFNGLGETUNIFORMLOCATIONPROC) wglGetProcAddress("glGetUniformLocation");
    glUniform1f = (PFNGLUNIFORM1FPROC) wglGetProcAddress("glUniform1f");
    glTexImage3D = getTexImage3D(); /*(PFNGLTEXIMAGE3DPROC) wglGetProcAddress("glTexImage3D");*/
    glActiveTexture = getActiveTexture(); /*(PFNGLACTIVETEXTUREPROC) wglGetProcAddress("glActiveTexture");*/

    if (! glActiveTexture) {
      FeedbackAdd(G, " Could not bind the glActiveTexture function.  No volume shaders.\n");
      abort();
    } 
    /* else  */
    /*   printf("glActiveTexture is: %x\n", glActiveTexture); */

    if (! glTexImage3D) {
      FeedbackAdd(G, " Could not bind the glTexImage3D function.  No volume shaders.\n");
      abort();
    } 
    /* else */
    /*   printf("glTexImage3D is: %x\n", glTexImage3D); */
  }
#endif

  /* First try to load configuration from files in $PYMOL_PATH, if they don't exist 
     load from const char * */
  //  FeedbackEnable(G, FB_ShaderMgr, FB_Everything);
  p = CShaderPrg_NewFromFile(G, "default", "default.vs", "default.fs");
  if (!p){
    PRINTFB(G, FB_ShaderMgr, FB_Warnings)
      " PyMOLShader_NewFromFile-Warning: default shader files not found, loading from memory.\n" ENDFB(G);
    p = CShaderPrg_New(G, "default", default_vs, default_fs);
  }
  hasShaders = (p!=0);

  CShaderMgr_AddShaderPrg(G->ShaderMgr, p);
  p = CShaderPrg_NewFromFile(G, "volume", "volume.vs", "volume.fs");
  if (!p){
    PRINTFB(G, FB_ShaderMgr, FB_Warnings)
      " PyMOLShader_NewFromFile-Warning: volume shader files not found, loading from memory.\n" ENDFB(G);
    p = CShaderPrg_New(G, "volume", volume_vs, volume_fs);
  }
  hasShaders &= (p!=0);
  CShaderMgr_AddShaderPrg(G->ShaderMgr, p);
  /* report GLSL version */
  if (G && G->Option && !G->Option->quiet) {
    getGLSLVersion(G, &major, &minor);
    sprintf(buf, " Detected GLSL version %d.%d.\n", major, minor);
    FeedbackAdd(G, buf);
  }

  G->ShaderMgr->ShadersPresent = hasShaders;

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
      " CShaderMgr_New-Error: Failed to create the shader manager.  Shaders disabled.\n" ENDFB(G)
    }
    return NULL;
  }


  I->G = G;
  DListInit(I->programs, prev, next, CShaderPrg);
  I->ShadersPresent = 0;
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
  CShaderPrg * p = NULL;
  DListIterate(I->programs, p, next) 
    {
      if (p && strcmp(p->name,name)==0) break;
    }

  return p;
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
    printf("size(%d)!=res(%d)\n", size,res);
    return NULL;
  }

  p[size] = 0;
  fclose(f);

  free(fullFile);
  return buffer;
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

  I->name = strdup(name);
  I->id = glCreateProgram();
  PRINTFB(G, FB_ShaderMgr, FB_Debugging)
    "Created program with id: %d\n", I->id ENDFB(G);

  /* VERTEX shader setup */
  /* CShaderPrg_InitShader(I, GL_VERTEX_SHADER); */
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
      PRINTFB(G, FB_ShaderMgr, FB_Errors) " CShaderPrg_New-Error: vertex shader compilation failed; log follows.\n" ENDFB(G);
      glGetShaderInfoLog(I->vid, MAX_LOG_LEN, &howLong, infoLog);
      PRINTFB(G, FB_ShaderMgr, FB_Errors)
	"infoLog=%s\n", infoLog ENDFB(G);
    }
    return 0;
  }
  PRINTFB(G, FB_ShaderMgr, FB_Debugging)
    "CShaderPrg_New-Message: vertex shader compiled.\n" ENDFB(G);
  glAttachShader(I->id, I->vid);

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
	" CShaderPrg-Error: fragment shader compilation failed; log follows.\n"
	ENDFB(G);
      glGetShaderInfoLog(I->fid, MAX_LOG_LEN, &howLong, infoLog);
      PRINTFB(G, FB_ShaderMgr, FB_Errors)
	"infoLog=%s\n", infoLog ENDFB(G);
    }
    return 0;
  }    
  glAttachShader(I->id, I->fid);
  
  /* Link the new program */
  CShaderPrg_Link(I);

  return I;
}


CShaderPrg * CShaderPrg_NewFromFile(PyMOLGlobals * G, const char * name, const char * vFile, const char * fFile)
{
  char *vFileStr, *fFileStr;
  vFileStr = CShaderMgr_ReadShaderFromDisk(G, vFile);
  if (!vFileStr){
    return (NULL);
  }
  fFileStr = CShaderMgr_ReadShaderFromDisk(G, fFile);
  if (!fFileStr){
    return (NULL);
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
      PRINTFB(G, FB_ShaderMgr, FB_Errors)
	" CShaderPrg_Link-Error: Shader program failed to link; log follows.\n"
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

