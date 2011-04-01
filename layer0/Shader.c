
/* -------------------------------------------------------------------
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
#include"os_gl.h"
#include <stdio.h>
#include <string.h>

#include "OOMac.h"
#include "Shader.h"
#include "PyMOLGlobals.h"
#include "PyMOLOptions.h"
#include "Feedback.h"

enum ShaderStatus { FAIL=-1, NOT_READY, READY };
enum ShaderReason { HW=-1, SW, OTHER };

struct _CShaderMgr {
  PyMOLGlobals * G;
  enum ShaderStatus status;
  enum ShaderReason reason;
  const char* v;
  const char* f;
};

/* allocate memory for the shader */
CShaderMgr* ShaderMgr_Alloc(PyMOLGlobals * G) {
  /* Singleton, don't reinit */
  if (G->ShaderMgr) return G->ShaderMgr;

  OOAlloc(G, CShaderMgr);
  return I;
}

/* intialize the shader object and 
 * initialize shaders on this machine */
void ShaderMgr_Init(CShaderMgr * s, PyMOLGlobals * G)
{
  int ok;
  CShaderMgr * I = s;
  I->G = G;
  I->status = NOT_READY;
  I->reason = OTHER;
  I->v = NULL;
  I->f = NULL;

  ok = ShaderMgr_InitShaders(I);

  if (ok)
    I->status = READY;
  else
    I->status = FAIL;
}

int ShaderMgr_InitShaders(CShaderMgr * s)
{
  CShaderMgr * I = s;
  return ShaderInit(I->G);
}


const char *vertex_shader_src =
  "varying vec3 N, L0, H0, L1, H1;\n"
  "varying vec4 D0, A0, D1, A1;\n"
  "void main()\n"
  "{\n"
  "  N = normalize(gl_NormalMatrix * gl_Normal);\n"
  "  vec3 eye_pos = vec3(gl_ModelViewMatrix * gl_Vertex);\n"
  "  vec3 aux = vec3(gl_LightSource[0].position - vec4(eye_pos, 1.0));\n"
  "  L0 = normalize(aux);\n"
  "  H0 = normalize(gl_LightSource[0].halfVector.xyz);\n"
  "  aux = vec3(gl_LightSource[1].position - vec4(eye_pos, 1.0));\n"
  "  L1 = normalize(aux);\n"
  "  H1 = normalize(gl_LightSource[1].halfVector.xyz);\n"
  "  A0 = gl_LightSource[0].ambient * gl_Color;\n"
  "  D0 = gl_LightSource[0].diffuse * gl_Color;\n"
  "  A1 = gl_LightSource[1].ambient * gl_Color;\n"
  "  D1 = gl_LightSource[1].diffuse * gl_Color;\n"
  "  gl_Position = ftransform();\n"
  "}";

const char *fragment_shader_src = 
  "varying vec3 N, L0, H0, L1, H1;\n"
  "varying vec4 D0, A0, D1, A1;\n"
  "void main()\n"
  "{\n"
  "  vec3 n, h;\n"
  "  float NdotL, NdotH;\n"
  "  vec4 color = A0 + A1;\n"
  "  n = normalize(N);\n"
  "  NdotL = max(dot(n, normalize(L0)), 0.0);\n"
  "  float shininess = gl_FrontMaterial.shininess;\n"
  "  if (NdotL > 0.0) {\n"
  "      color += D0 * NdotL;\n"
  "      h = normalize(H0);\n"
  "      NdotH = max(dot(n, h), 0.0);\n"
  "      color += gl_LightSource[0].specular * pow(NdotH, shininess);\n"
  "  }\n"
  "  NdotL = max(dot(n, normalize(L1)), 0.0);\n"
  "  if (NdotL > 0.0) {\n"
  "      color += D1 * NdotL;\n"
  "      h = normalize(H1);\n"
  "      NdotH = max(dot(n, h), 0.0);\n"
  "      color += gl_LightSource[1].specular * pow(NdotH, shininess);\n"
  "  }\n"
  "  gl_FragColor = color;\n"
  "}";

GLuint shader_program = 0;

/* Initializes all shaders. This should be done once at the time of
 * GL context initialization. 
 */
int ShaderInit(PyMOLGlobals *G)
{
  GLuint vert_shader, frag_shader;
  GLuint status;

  char buf[50];
  int major, minor;

  int MAX_SHADER_LEN = 1024;
  int howLong;
  char infoLog[MAX_SHADER_LEN];

  /* report GLSL version */
  if (G && G->Option && !G->Option->quiet) {
    getGLSLVersion(G, &major, &minor);
    sprintf(buf, "Detected GLSL version %d.%d\n", major, minor);
    FeedbackAdd(G, buf);
  }

  // initalize the ShaderManager and load the proper shaders given the version info

  /* create the shader objects */
  vert_shader = glCreateShader(GL_VERTEX_SHADER);
  frag_shader = glCreateShader(GL_FRAGMENT_SHADER);

  /* assign the objects their source code */
  glShaderSource(vert_shader, 1, &vertex_shader_src, NULL);
  glShaderSource(frag_shader, 1, &fragment_shader_src, NULL);

  /* compile the vertex shader */
  glCompileShader(vert_shader);
  /* look for vertex compilation errors */
  glGetShaderiv(vert_shader, GL_COMPILE_STATUS, &status);
  if (!status) {
    if (G && G->Option && !G->Option->quiet) {
      FeedbackAdd(G, " GLSL: (ERROR) vertex shader compilation failed; shaders disabled; log follows.\n");
      glGetShaderInfoLog(vert_shader, MAX_SHADER_LEN, &howLong, infoLog);
      FeedbackAdd(G, infoLog);
    }
    return 0;
  } else if (G && G->Option && !G->Option->quiet) {
    FeedbackAdd(G, " GLSL: shader compilation successful (vertex).\n");
  }
  /* compile the fragment shaders */
  glCompileShader(frag_shader);
  /* look for fragment shader compilation errors */
  glGetShaderiv(frag_shader, GL_COMPILE_STATUS, &status);
  if (!status) {
    if (G && G->Option && !G->Option->quiet) {
      FeedbackAdd(G, "GLSL: (ERROR) fragment shader compilation failed; shaders disabled.\n");
    }
    return 0;
  } else if (G && G->Option && !G->Option->quiet) {
    /* alert the user */
    FeedbackAdd(G, " GLSL: shader compilation successful (fragment).\n");
  }

  /* create the shader program */
  shader_program = glCreateProgram();

  glAttachShader(shader_program, vert_shader);
  glAttachShader(shader_program, frag_shader);

  /* link */
  glLinkProgram(shader_program);

  if (G && G->Option && !G->Option->quiet) {
    char s[50];
    sprintf(s, " GLSL: (DEBUG) shader_program = %d\n\n", shader_program);
    FeedbackAdd(G, s);
  }
  
  return 1;
}

int ShaderEnable(PyMOLGlobals *G)
{
#if 0
  CShaderMgr * I;
  int ok;

  if (G && G->ShaderMgr)
    I = G->ShaderMgr;
  else 
    if (G)
      FeedbackAdd(G, "  GLSL: (ERROR) PyMOL has not initialized its ShaderManager.\n");

  if (!I->status)
    ok = ShaderMgr_InitShaders(I);
    
  if (ok) {
#endif
    glUseProgram(shader_program);
    return 1;
#if 0
  }
  else {
    FeedbackAdd(G, "  GLSL: (ERROR) Could not initialize openGL shaders on this machine.\n");
    return 0;
  }
  return 0;
#endif
}

int ShaderDisable(PyMOLGlobals *G) 
{
  glUseProgram(0);
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
      "Invalid GL_VERSION format.\n" ENDFD(G);
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
	  "Invalid GL_SHADING_LANGUAGE_VERSION format.\n" ENDFD(G);
      }
    }
  }
}

#endif
