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

#include "os_gl.h"
#include "os_python.h"
#include <string.h>
#include <iostream>
#include "ShaderMgr.h"
#include "OOMac.h"
#include "ListMacros.h"
#include "PyMOLOptions.h"
#include "Feedback.h"
#include "MemoryDebug.h"
#include "Setting.h"
#include "Scene.h"
#include "Color.h"
#include "Vector.h"
#include "Util.h"
#include "Util2.h"
#include "Texture.h"
#include "FileStream.h"
#include "Matrix.h"
#include "Parse.h"

#ifndef _PYMOL_NO_AA_SHADERS
#endif

#include "CGO.h"
#ifdef _WEBGL
#include "Matrix.h"
#include "WebPyMOLLibrary.h"
#endif
#define MAX_LOG_LEN 1024

#include <algorithm>
#include <sstream>
#include <stack>
#include <vector>
#include <functional>
/*  Texture Usage:

    0 - for not-PURE_OPENGL_ES_2: ObjectVolume: volumeTex
        for WEBGL: SCHRODINGER logo
    1 - for not-PURE_OPENGL_ES_2: Volume: either colorTex1D or colorTex2D
        for _PYMOL_PRECOMPUTED_LIGHTING: Lighting Texture (ShaderMgr->lightingTexture)
    2 - FXAA - color_texture
        SMAA1 - colorTex
        SMAA3 - colorTex
    3 - SMAA3 - blendTex
        Label/Indicator Shader : textureMap (both PURE_OPENGL_ES_2 and non-PURE_OPENGL_ES_2
    4 - Background Texture: bgTextureMap
    5 - OIT - 2nd pass : accumTex
        Volumes - carvemask
    6 - OIT - 2nd pass : revealageTex
        SMAA2 - edgesTex
    7 - SMAA2 - areaTex
        OIT Copy - colorTex
    8 - SMAA2 - searchTex

 */

#ifndef _DEAD_CODE_DIE
#define SUPPRESS_GEOMETRY_SHADER_ERRORS
#endif

#define CONNECTOR_GS_NUM_VERTICES  31

#define SCENEGETIMAGESIZE SceneGetWidthHeight

#include "ShaderText.h"

#define WARNING_IF_GLERROR(msg) {		\
  GLenum err; \
  if ((err = glGetError())){ \
    PRINTFB(G, FB_ShaderMgr, FB_Warnings) "GLERROR 0x%04x: %s\n", err, msg ENDFB(G); \
  } \
}

static void glShaderSource1String(GLuint shad, const std::string &strobj){
  const GLchar *str = (const GLchar *)strobj.c_str();
  glShaderSource(shad, 1, (const GLchar **)&str, nullptr);
}

bool CShaderPrg::reload(){
  // skip programs with empty file names, assume their code is managed
  // outside of the reload logic (like ARB shaders).
  if (is_valid || vertfile.empty())
    return true;

  std::string gs, vs, fs;
  CShaderMgr *I = G->ShaderMgr;
  GLint status;

  if (!geomfile.empty())
    gs = I->GetShaderSource(geomfile);

  vs = I->GetShaderSource(vertfile);
  fs = I->GetShaderSource(fragfile);

  WARNING_IF_GLERROR("CShaderPrg::reload begin");

  PRINTFB(G, FB_ShaderMgr, FB_Blather)
    "Loading shader named: %s\n", name.c_str()
    ENDFB(G);

  if (!id) {
    id = glCreateProgram();
  }

#ifndef PURE_OPENGL_ES_2
  if (!gs.empty() && SettingGetGlobal_b(G, cSetting_use_geometry_shaders)) {
    if (!gid) {
      gid = glCreateShader(GL_GEOMETRY_SHADER);

      GLenum err;
      if ((err=glGetError()) || !gid) {
        PRINTFB(G, FB_ShaderMgr, FB_Errors)
          " Error: geometry shader creation failed. name=%s err=0x%x\n", name.c_str(), err ENDFB(G);
        return false;
      }

      glAttachShader(id, gid);
    }

    glShaderSource1String(gid, gs);
    glCompileShader((GLuint) gid);
    glGetShaderiv(gid, GL_COMPILE_STATUS, &status);

    if (!status) {
#ifndef SUPPRESS_GEOMETRY_SHADER_ERRORS
      ErrorMsgWithShaderInfoLog(gid, "geometry shader compilation failed.");
#endif
      glDetachShader(id, gid);
      glDeleteShader(gid);
      gid = 0;
      return false;
    }

    glProgramParameteriEXT(id, GL_GEOMETRY_INPUT_TYPE_EXT, gsInput);
    glProgramParameteriEXT(id, GL_GEOMETRY_OUTPUT_TYPE_EXT, gsOutput);
    glProgramParameteriEXT(id, GL_GEOMETRY_VERTICES_OUT_EXT, ngsVertsOut);

    PRINTFB(G, FB_ShaderMgr, FB_Debugging)
      " ShaderPrg-Debug: geometry shader compiled.\n" ENDFB(G);
  } else if (gid) {
    // for manually switching off geometry shaders (set use_geometry_shaders, off)
    glDetachShader(id, gid);
    glDeleteShader(gid);
    gid = 0;
  }

  WARNING_IF_GLERROR("CShaderPrg::reload after geometry shader");
#endif

  // vertex shader
  {
    if (!vid) {
      vid = glCreateShader(GL_VERTEX_SHADER);
      glAttachShader(id, vid);
    }

    glShaderSource1String(vid, vs);
    glCompileShader((GLuint) vid);
    glGetShaderiv(vid, GL_COMPILE_STATUS, &status);

    if (!status) {
      ErrorMsgWithShaderInfoLog(vid, "vertex shader compilation failed.");
      return false;
    }
  }

  // fragment shader
  {
    if (!fid) {
      fid = glCreateShader(GL_FRAGMENT_SHADER);
      glAttachShader(id, fid);
    }

    glShaderSource1String(fid, fs);
    glCompileShader((GLuint) fid);
    glGetShaderiv(fid, GL_COMPILE_STATUS, &status);

    if (!status) {
      ErrorMsgWithShaderInfoLog(fid, "fragment shader compilation failed.");
      return false;
    }
  }

  uniforms.clear();
  uniform_set = 0;

  // it is valid to bind unused names, and to bind multiple names to the same index
  if (!name.compare(0, 8, "cylinder")){
    glBindAttribLocation(id, CYLINDER_VERTEX1, "attr_vertex1");
    glBindAttribLocation(id, CYLINDER_VERTEX2, "attr_vertex2");
    glBindAttribLocation(id, CYLINDER_COLOR, "a_Color");
    glBindAttribLocation(id, CYLINDER_COLOR2, "a_Color2");
    glBindAttribLocation(id, CYLINDER_RADIUS, "attr_radius");
    glBindAttribLocation(id, CYLINDER_CAP, "a_cap");
  } else {
    glBindAttribLocation(id, VERTEX_POS, "a_Vertex");
    glBindAttribLocation(id, VERTEX_COLOR, "a_Color");
    glBindAttribLocation(id, VERTEX_NORMAL, "a_Normal");
    glBindAttribLocation(id, 0, "attr_worldpos");
  }
  WARNING_IF_GLERROR("after glBindAttribLocation");

  is_linked = false;
  is_valid = true;

  return true;
}

#define MASK_SHADERS_PRESENT_GEOMETRY 0x2;
#define MASK_SHADERS_PRESENT_SMAA 0x4;

static void getGLVersion(PyMOLGlobals * G, int *major, int* minor);
static void getGLSLVersion(PyMOLGlobals * G, int* major, int* minor);

static void disableShaders(PyMOLGlobals * G);

#ifdef WIN32
/* REMOVE US */
PFNGLTEXIMAGE3DPROC getTexImage3D(){
  static PFNGLTEXIMAGE3DPROC my_glTexImage3D = NULL;
  if (!my_glTexImage3D)
    my_glTexImage3D = (PFNGLTEXIMAGE3DPROC) wglGetProcAddress("glTexImage3D");
  return my_glTexImage3D;
}
#endif

/**
 * Use this to turn off shaders if the renderer cannot use them.
 */
void disableShaders(PyMOLGlobals * G) {
    /* Auto-disable shader-based rendering */
    SettingSetGlobal_b(G, cSetting_use_shaders, false);
}

static void disableGeometryShaders(PyMOLGlobals * G) {
  SettingSetGlobal_b(G, cSetting_use_geometry_shaders, false);
  if(G->ShaderMgr)
    G->ShaderMgr->SetPreprocVar("use_geometry_shaders", false);

  if (G->Option && !G->Option->quiet)
    PRINTFB(G, FB_ShaderMgr, FB_Warnings)
      " Geometry shaders not available\n" ENDFB(G);
}

/**
 * Replace strings from a list of pairs.
 *
 * src: string to modify
 * replaceStrings: map of strings to replace (as consecutive elements in an
 *                 array like {from1, to1, from2, to2, ..., ""}
 * returns: new string
 */
static std::string stringReplaceAll(const std::string &src, const std::string * replaceStrings) {
  std::string dest = src;
  for (int i = 0; !replaceStrings[i].empty(); i += 2) {
    int slen1 = replaceStrings[i].length();
    int slen2 = replaceStrings[i + 1].length();
    for (size_t pl = 0;
        (pl = dest.find(replaceStrings[i], pl)) != std::string::npos;
        pl += slen2) {
      dest.replace(pl, slen1, replaceStrings[i + 1]);
    }
  }
  return dest;
}

/**
 * Reload "CallComputeColorForLight" shader replacement string
 */
void CShaderMgr::Reload_CallComputeColorForLight(){
  if ((reload_bits & RELOAD_CALLCOMPUTELIGHTING)) {
    reload_bits &= ~RELOAD_CALLCOMPUTELIGHTING;
  } else {
    return;
  }

  if (SettingGetGlobal_b(G, cSetting_precomputed_lighting)) {
    Generate_LightingTexture();
    return;
  }

  int light_count = SettingGetGlobal_i(G, cSetting_light_count);
  int spec_count = SettingGetGlobal_i(G, cSetting_spec_count);
  std::ostringstream accstr;

  std::string rawtemplate = GetShaderSource("call_compute_color_for_light.fs");

  std::string lightstrings[] = {
    "`light`", "0",
    "`postfix`", "_0",
    ""
  };

  accstr << stringReplaceAll(rawtemplate, lightstrings);

  if (light_count > 8){
    PRINTFB(G, FB_ShaderMgr, FB_Details)
      " ShaderMgr-Detail: using 8 lights (use precomputed_lighting for light_count > 8)\n"
      ENDFB(G);
    light_count = 8;
  }

  // no postfix for 1..light_count
  lightstrings[3] = "";

  for (int i=1; i<light_count; i++){
    std::ostringstream lstr;
    lstr << i;
    lightstrings[1] = lstr.str(); // std::to_string(i)

    if (i == spec_count + 1) {
      // no specular for [spec_count + 1 .. light_count]
      lightstrings[3] = " * 0.0";
    }

    accstr << stringReplaceAll(rawtemplate, lightstrings);
  }

  SetShaderSource("CallComputeColorForLight", accstr.str());
}

void CShaderMgr::Invalidate_All_Shaders(){
  for (auto& prog : programs) {
    prog.second->Invalidate();
  }
}

void CShaderMgr::Reload_All_Shaders(){
  Reload_Shader_Variables();
  Reload_CallComputeColorForLight();

  if (SettingGetGlobal_i(G, cSetting_transparency_mode) == 3) {
    Reload_Derivatives("NO_ORDER_TRANSP");
  }

  for (auto& prog : programs) {
    if (prog.second->derivative.empty())
      prog.second->reload();
  }
}

// bitmasks for preprocessor parsing
#define IFDEF    1   // #ifdef or #ifndef
#define IFNDEF   2   // #ifndef
#define ELSE     4   // #else
#define ENDIF    8   // #endif
#define INCLUDE 16   // #include
#define LOOKUP  32   // #ifdef or #ifndef or #include

// preprocessor directive (like '#ifdef') -> bitmask
static std::map<std::string, short> preprocmap;

// filename -> contents (static filesystem)
static std::map<std::string, const char *> shader_cache_raw;

// preproc variable -> NULL terminated list of filenames ("used by")
std::map<std::string, const char **> ifdef_deps;

// filename -> NULL terminated list of filenames ("included by")
std::map<std::string, const char **> include_deps;

/**
 * Return a pointer to the next whitespace character or to the end of the string
 */
static const char * nextwhitespace(const char * p) {
  for (;; p++) {
    switch (*p) {
      case ' ': case '\0': case '\n': case '\r': case '\t':
        return p;
    }
  }
}

/**
 * Return a pointer to the next line beginning or to the end of the string.
 * Skips blank lines.
 */
static const char * nextline(const char * p) {
  for (;; p++) {
    switch (*p) {
      case '\0': case '\n': case '\r':
        goto switch2;
    }
  }
  for (;; p++) {
switch2:
    switch (*p) {
      case ' ': case '\n': case '\r': case '\t':
        break;
      default:
        return p;
    }
  }
}

/**
 * Get the processed shader file contents with all `#ifdef` and `#include`
 * preprocessors processed.
 *
 * Note: There must be a single whitespace character between preprocessor
 * directive and argument.
 *
 * Valid:
 *
 *     #ifdef foo
 *
 * Invalid:
 *
 *     # ifdef foo
 *     #ifdef  foo
 *
 * @param filename file name of the shader file inside `$PYMOL_DATA/shaders`
 */
std::string CShaderMgr::GetShaderSource(const std::string &filename)
{
  // processed cache
  auto it = shader_cache_processed.find(filename);
  if (it != shader_cache_processed.end()) {
    return it->second;
  }

  std::string buffer;
  const char *pl = nullptr, *newpl, *tpl;
  std::ostringstream newbuffer;

  /* "if_depth" counts the level of nesting, and "true_depth" how far the
   * if conditions were actually true. So if the current block is true, then
   * if_depth == true_depth, otherwise if_depth > true_depth.
   */
  int if_depth = 0, true_depth = 0;

#ifndef _PYMOL_IOS
  /* read the file from disk */
  if (SettingGetGlobal_b(G, cSetting_shaders_from_disk)) {
    const char * pymol_data = getenv("PYMOL_DATA");

    if (pymol_data && pymol_data[0]) {
      std::string path(pymol_data);
      path.append(PATH_SEP).append("shaders").append(PATH_SEP).append(filename);

      try {
        buffer = pymol::file_get_contents(path);
        pl = buffer.c_str();
      } catch (...) {
        PRINTFB(G, FB_ShaderMgr, FB_Warnings)
          " Warning: shaders_from_dist=on, but unable to open file '%s'\n",
          path.c_str() ENDFB(G);
      }
    } else {
      PRINTFB(G, FB_ShaderMgr, FB_Warnings)
        " Warning: shaders_from_dist=on, but PYMOL_DATA not set\n" ENDFB(G);
    }
  }
#endif

  if (!pl) {
    pl = shader_cache_raw[filename];
    if (!pl) {
      PRINTFB(G, FB_ShaderMgr, FB_Errors)
        " GetShaderSource-Error: No such file: '%s'\n", filename.c_str() ENDFB(G);
      return "";
    }
  }

  /* Now we need to read through the shader and do processing if necessary */
  for (; *pl; pl = newpl) {
    int preproc = 0;

    // only do preprocessor lookup if line starts with a hash
    if (pl[0] == '#') {
      // next white space
      tpl = nextwhitespace(pl);

      // copy of first word
      std::string tmp_str(pl, tpl - pl);

      // lookup word in preprocmap
      std::map<std::string, short>::const_iterator
        preprocit = preprocmap.find(tmp_str);
      if (preprocit != preprocmap.end()) {
        preproc = preprocit->second;

        if (preproc & LOOKUP) { // #ifdef or #ifndef or #include
          if (if_depth == true_depth) {
            // copy of second word
            tpl++;
            tmp_str = std::string(tpl, nextwhitespace(tpl) - tpl);

            if (preproc & IFDEF) { // #ifdef or #ifndef
              bool if_value = false;

              // lookup for boolean shader preprocessor values
              auto item = preproc_vars.find(tmp_str);
              if (item != preproc_vars.end())
                if_value = item->second;

              if (preproc & IFNDEF)
                if_value = !if_value; // #ifndef

              if (if_value)
                true_depth++;

            } else if (preproc & INCLUDE) { //#include
              tmp_str = std::string(tpl, nextwhitespace(tpl) - tpl);
              newbuffer << GetShaderSource(tmp_str);
            }
          }

          if (preproc & IFDEF)
            if_depth++;

        } else if (preproc & ENDIF){ // #endif
          if (if_depth-- == true_depth)
            true_depth--;
        } else if (preproc & ELSE){ // #else
          if (if_depth == true_depth)
            true_depth--;
          else if (if_depth == true_depth + 1)
            true_depth++;
        }
      }
    }

    newpl = nextline(pl);

    // add to the output buffer if this is a regular active line
    if (!preproc && if_depth == true_depth) {
      newbuffer.write(pl, newpl - pl);
    }
  }

  std::string result = newbuffer.str();
  shader_cache_processed[filename] = result;
  return result;
}

#define FREE_AND_REPLACE_WITH(var, with) if (var) free(var);  var = with;

void CShaderMgr::Reload_Shader_Variables() {
  if ((reload_bits & RELOAD_VARIABLES)) {
    reload_bits &= ~RELOAD_VARIABLES;
  } else {
    return;
  }

  int bg_image_mode = SettingGetGlobal_i(G, cSetting_bg_image_mode);  
  int bg_gradient = SettingGetGlobal_b(G, cSetting_bg_gradient);  
  int bg_image_mode_solid;
  int stereo, stereo_mode;
  const char * bg_image_filename = SettingGet_s(G, nullptr, nullptr, cSetting_bg_image_filename);
  short bg_image = bg_image_filename && bg_image_filename[0];
  bg_image_mode_solid = !(bg_gradient || bg_image || OrthoBackgroundDataIsSet(*G->Ortho));

  SetPreprocVar("bg_image_mode_solid", bg_image_mode_solid);
  if (!bg_image_mode_solid) {
    SetPreprocVar("bg_image_mode_1_or_3", (bg_image_mode == 1 || bg_image_mode == 3));
    SetPreprocVar("bg_image_mode_2_or_3", (bg_image_mode == 2 || bg_image_mode == 3));
  }

#ifdef _PYMOL_IP_EXTRAS
  SetPreprocVar("volume_mode", SettingGetGlobal_i(G, cSetting_volume_mode));
#endif

  SetPreprocVar("ortho", SettingGetGlobal_i(G, cSetting_ortho));
  SetPreprocVar("depth_cue", SettingGetGlobal_b(G, cSetting_depth_cue)
      && SettingGetGlobal_b(G, cSetting_fog) != 0.0F);

#ifndef PURE_OPENGL_ES_2
  SetPreprocVar("use_geometry_shaders", SettingGetGlobal_b(G, cSetting_use_geometry_shaders));
#endif

  SetPreprocVar("line_smooth", SettingGetGlobal_b(G, cSetting_line_smooth));

  stereo = SettingGetGlobal_i(G, cSetting_stereo);
  stereo_mode = SettingGetGlobal_i(G, cSetting_stereo_mode);

#ifdef _PYMOL_OPENVR
  SetPreprocVar("openvr_enabled", stereo && stereo_mode == cStereo_openvr);
#endif

  SetPreprocVar("ANAGLYPH", stereo && stereo_mode == cStereo_anaglyph);
  SetPreprocVar("ray_trace_mode_3", SettingGetGlobal_i(G, cSetting_ray_trace_mode) == 3);
  SetPreprocVar("transparency_mode_3", SettingGetGlobal_i(G, cSetting_transparency_mode)==3);

#ifndef _PYMOL_NO_AA_SHADERS
#endif

  SetPreprocVar("precomputed_lighting", SettingGetGlobal_b(G, cSetting_precomputed_lighting));
  SetPreprocVar("ray_transparency_oblique", SettingGetGlobal_f(G, cSetting_ray_transparency_oblique) > R_SMALL4);

  int chromadepth = SettingGetGlobal_i(G, cSetting_chromadepth);
  SetPreprocVar("chromadepth", chromadepth != 0);
  SetPreprocVar("chromadepth_postlighting", chromadepth == 2);
}

/* ============================================================================
 * ShaderMgrInit is called from PyMOL.c during start up; it just allocates
 * the global ShaderMgr
 */
bool ShaderMgrInit(PyMOLGlobals * G) {
  // initialize some globals (do this only once)
  if (preprocmap.empty()) {
    preprocmap["#ifdef"] = LOOKUP | IFDEF;
    preprocmap["#ifndef"] = LOOKUP | IFDEF | IFNDEF;
    preprocmap["#else"] = ELSE;
    preprocmap["#endif"] = ENDIF;
    preprocmap["#include"] = LOOKUP | INCLUDE;

    // make #include dependency map from flat array
    for (const char ** ptr = _include_deps; *ptr; ++ptr) {
      include_deps[ptr[0]] = ptr + 1;
      while (*(++ptr)) {}
    }

    // make #ifdef dependency map from flat array
    for (const char ** ptr = _ifdef_deps; *ptr; ++ptr) {
      ifdef_deps[ptr[0]] = ptr + 1;
      while (*(++ptr)) {}
    }

    // make shader file cache from flat array
    for (const char ** ptr = _shader_cache_raw; *ptr; ptr += 2) {
      shader_cache_raw[ptr[0]] = *(ptr + 1);
    }
  }

  G->ShaderMgr = new CShaderMgr(G);

  if(!G->ShaderMgr)
    return false;

  return true;
}

/**
 * Print the given message as ShaderMgr-Error, followed by the shader info log.
 */
void CShaderPrg::ErrorMsgWithShaderInfoLog(const GLuint sid, const char * msg) {
  if (!G->Option || G->Option->quiet)
    return;

  GLint infoLogLength = 0;
  glGetShaderiv(sid, GL_INFO_LOG_LENGTH, &infoLogLength);
  std::vector<GLchar> infoLog(infoLogLength);
  glGetShaderInfoLog(sid, infoLogLength, nullptr, infoLog.data());

  PRINTFB(G, FB_ShaderPrg, FB_Errors) " ShaderPrg-Error: %s; name='%s'\n",
    msg, name.c_str() ENDFB(G);

  PRINTFB(G, FB_ShaderPrg, FB_Errors) " ShaderPrg-Error-InfoLog:\n%s\n",
    infoLog.data() ENDFB(G);
}

/* ShaderMgrConfig -- Called from PyMOL.c, configures the global ShaderMgr
 * This needs to be called once the OpenGL context has been created, it is 
 * called from MainInit() for PyMol, and from PyMOL_ConfigureShadersGL() for
 * other programs (i.e., JyMOL, AxPyMOL, etc.).
 */
void CShaderMgr::Config() {
  if (!G || !G->HaveGUI) /* && G->ValidContext); */
    return;

  glGetFloatv(GL_ALIASED_LINE_WIDTH_RANGE, line_width_range);

#ifndef PURE_OPENGL_ES_2
  GLenum err = glewInit();

  if (GLEW_OK==err) {
    GLint gl_major = 0, gl_minor = 0;
    getGLVersion(G, &gl_major, &gl_minor);

    PRINTFB(G, FB_ShaderMgr, FB_Details)
      " Detected OpenGL version %d.%d.", gl_major, gl_minor ENDFB(G);

    if (GLEW_VERSION_2_0) {
      G->Feedback->add(" Shaders available.\n");
    }
    else { 
      G->Feedback->add(" Shaders and volumes unavailable.\n");
      disableShaders(G);
      return;
    }
  } 
  else {
    /* print info on glew error? */
    G->Feedback->add(" There was an error intializing GLEW.  Basic graphics, including\n shaders and volumes may be unavailable.\n");
    disableShaders(G);
    fprintf(stderr, " GLEW-Error: %s\n", glewGetErrorString(err));
    return;
  }
#endif

  // static preprocessor values
  preproc_vars["GLEW_VERSION_3_0"] = GLEW_VERSION_3_0 ? true : false;
  if (TM3_IS_ONEBUF){
    preproc_vars["ONE_DRAW_BUFFER"] = true;
  }
#ifdef PURE_OPENGL_ES_2
  preproc_vars["PURE_OPENGL_ES_2"] = 1;
  preproc_vars["PYMOL_WEBGL"] = 1;
  preproc_vars["PYMOL_WEBGL_IOS"] = 1;
#else
  preproc_vars["gl_VertexID_enabled"] = GLEW_EXT_gpu_shader4;
#endif

  // shaders
#define make_program(name, ...) programs[name] = new CShaderPrg(G, name, __VA_ARGS__)

  make_program("bg", "bg.vs", "bg.fs");
  make_program("indicator", "indicator.vs", "indicator.fs");
  make_program("label", "label.vs", "label.fs");

#ifndef PURE_OPENGL_ES_2
  make_program("volume", "volume.vs", "volume.fs");
#endif

  make_program("default", "default.vs", "default.fs");
  make_program("surface", "surface.vs", "surface.fs");

  make_program("line", "line.vs", "line.fs");

  make_program("screen", "screen.vs", "screen.fs");

  if (GLEW_EXT_geometry_shader4 && GLEW_EXT_gpu_shader4){
    make_program("connector", "connector.vs", "connector.fs",
		 "connector.gs", GL_POINTS, GL_TRIANGLE_STRIP, CONNECTOR_GS_NUM_VERTICES);
  } else {
    make_program("connector", "connector.vs", "connector.fs");
  }

  if (GET_FRAGDEPTH_SUPPORT()) {
    make_program("cylinder", "cylinder.vs", "cylinder.fs");
    make_program("sphere", "sphere.vs", "sphere.fs");
  }  

  make_program("ramp", "ramp.vs", "ramp.fs");
  programs["ramp"]->uniformLocations[RAMP_OFFSETPT] = "offsetPt";

  make_program("oit", "oit.vs", "oit.fs");
  make_program("copy", "copy.vs", "copy.fs");
  make_program("trilines", "trilines.vs", "trilines.fs");

  Reload_Shader_Variables();
  Reload_CallComputeColorForLight();

  // shaders availability test
  ok_assert(1, programs["default"]->reload());

#ifndef PURE_OPENGL_ES_2
  // geometry shaders availability test
  if (programs["connector"]->reload() && programs["connector"]->gid) {
    shaders_present |= MASK_SHADERS_PRESENT_GEOMETRY;
  } else {
    disableGeometryShaders(G);
  }
#else
  SettingSetGlobal_b(G, cSetting_use_geometry_shaders, 0);
#endif

#define check_program(name, setting, value) {   \
    if (!programs[name]->reload()) {            \
      SettingSetGlobal_i(G, setting, value);    \
      programs.erase(name);                     \
    }}                                          \
  // other shader compilation tests
  if (GET_FRAGDEPTH_SUPPORT()) {
  check_program("cylinder", cSetting_render_as_cylinders, 0);
  check_program("sphere", cSetting_sphere_mode, 0);
  }

#ifndef _PYMOL_NO_AA_SHADERS
#endif

  // get filename -> shader program dependencies
  for (auto& prog : programs) {
    RegisterDependantFileNames(prog.second);
  }

  // make transparency_mode_3 shader derivatives
  MakeDerivatives("_t", "NO_ORDER_TRANSP");

#ifndef PURE_OPENGL_ES_2
  /* report GLSL version */
  if (G && G->Option && !G->Option->quiet) {
    char buf[50];
    int major, minor;
    getGLSLVersion(G, &major, &minor);
    sprintf(buf, " Detected GLSL version %d.%d.\n", major, minor);
    G->Feedback->add(buf);
  }
#endif
  shaders_present |= 0x1;
  SettingSetGlobal_b(G, cSetting_use_shaders, true);

  is_configured = true;
  return;
ok_except1:
  disableShaders(G);
  G->ShaderMgr->shaders_present = 0;
  is_configured = true;
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
  if (!verstr || sscanf(verstr, "%d.%d", major, minor) != 2) {
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
#ifndef PURE_OPENGL_ES_2
void getGLSLVersion(PyMOLGlobals * G, int* major, int* minor) {
  int gl_major, gl_minor;
  *major = *minor = 0;

  /* grab the GL version */
  getGLVersion(G, &gl_major, &gl_minor);

  /* GL version 1 */
  if (1==gl_major) {
    const char* extstr = (const char*) glGetString(GL_EXTENSIONS);
    if (extstr && strstr(extstr, "GL_ARB_shading_language_100")) {
      *major = 1;
      *minor = 0;
    }
  }
  /* GL > version 1 */
  else if (gl_major>=2) {
    const char* verstr = (const char*) glGetString(GL_SHADING_LANGUAGE_VERSION);

    if (!verstr || sscanf(verstr, "%d.%d", major, minor) != 2) {
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
CShaderMgr::CShaderMgr(PyMOLGlobals * G_)
{
  G = G_;
  current_shader = nullptr;
  shaders_present = 0;
  stereo_flag = 0;
  stereo_blend = 0;
#ifdef _PYMOL_LIB
  print_warnings = 0;
#else
  print_warnings = 1;
#endif

  lightingTexture = 0;
  is_picking = 0;
  reload_bits = RELOAD_ALL_SHADERS;

#ifndef _WEBGL
  vbos_to_free.reserve(256);
#endif
}

CShaderMgr::~CShaderMgr() {
  for (auto& prog : programs) {
    delete prog.second;
  }
  programs.clear();

  shader_cache_processed.clear();

  freeGPUBuffer(offscreen_rt);

  FreeAllVBOs();
}

int CShaderMgr::AddShaderPrg(CShaderPrg * s) {
  if (!s)
    return 0;
  const std::string& name = s->name;
  if (programs.find(name)!=programs.end()){
    delete programs[name];
  }
  programs[name] = s;
  return 1;
}

int CShaderMgr::RemoveShaderPrg(const std::string& name) {
  if (programs.find(name) != programs.end()){
    delete programs[name];
  }
  return 1;
}

/**
 * Lookup a shader program by name and set it as the `current_shader` of the
 * shader manager. If `pass` is provided and is less than zero, and we are
 * in transparency_mode 3, then look up the NO_ORDER_TRANSP derivative.h
 */
CShaderPrg * CShaderMgr::GetShaderPrg(std::string name, short set_current_shader, RenderPass pass) {
  if (pass == RenderPass::Transparent && SettingGetGlobal_i(G, cSetting_transparency_mode) == 3) {
    name += "_t";
  }

  auto it = programs.find(name);
  if (it == programs.end())
    return nullptr;

  if (set_current_shader)
    current_shader = it->second;

  return it->second;
}

int CShaderMgr::ShaderPrgExists(const char * name){
  return (programs.find(name) != programs.end());
}

int CShaderMgr::ShadersPresent() {
  return shaders_present;
}

int CShaderMgr::GeometryShadersPresent() {
  return shaders_present & MASK_SHADERS_PRESENT_GEOMETRY;
}

/**
 * glDeleteBuffers for vbos_to_free
 */
void CShaderMgr::FreeAllVBOs() {
#ifndef _WEBGL
  freeAllGPUBuffers();

  LOCK_GUARD_MUTEX(lock, vbos_to_free_mutex);

  if (vbos_to_free.empty())
    return;

  glDeleteBuffers(vbos_to_free.size(), &vbos_to_free[0]);

  vbos_to_free.clear();
#endif
}

void CShaderMgr::AddVBOsToFree(GLuint *vboid, int nvbos){
  int i;
  for (i=0; i<nvbos; i++){
    if (vboid[i]>0)
        AddVBOToFree(vboid[i]);
  }
}

/**
 * thread-safe deferred glDeleteBuffers(1, &vboid)
 */
void CShaderMgr::AddVBOToFree(GLuint vboid){
#ifdef _WEBGL // No threads, immediately delete
  if (glIsBuffer(vboid)) {
    glDeleteBuffers(1, &vboid);
  } else {
    PRINTFB(G, FB_ShaderMgr, FB_Warnings) "WARNING: CShaderMgr_AddVBOToFree() buffer is not a VBO %d", vboid ENDFB(G);
  }
#else
  LOCK_GUARD_MUTEX(lock, vbos_to_free_mutex);
  vbos_to_free.push_back(vboid);
#endif
}

CShaderPrg *CShaderMgr::Enable_DefaultShaderWithSettings(
    const CSetting *set1,
    const CSetting *set2, RenderPass pass) {
  CShaderPrg * shaderPrg = Get_DefaultShader(pass);
  return Setup_DefaultShader(shaderPrg, set1, set2);
}

CShaderPrg *CShaderMgr::Enable_DefaultShader(RenderPass pass){
  CShaderPrg * shaderPrg = Get_DefaultShader(pass);
  return Setup_DefaultShader(shaderPrg, nullptr, nullptr);
}

CShaderPrg *CShaderMgr::Enable_LineShader(RenderPass pass){
  CShaderPrg * shaderPrg = Get_LineShader(pass);
  return Setup_DefaultShader(shaderPrg, nullptr, nullptr);
}

CShaderPrg *CShaderMgr::Enable_SurfaceShader(RenderPass pass){
  CShaderPrg * shaderPrg = Get_SurfaceShader(pass);
  return Setup_DefaultShader(shaderPrg, nullptr, nullptr);
}

CShaderPrg *CShaderMgr::Enable_ConnectorShader(RenderPass pass){
  CShaderPrg * shaderPrg = Get_ConnectorShader(pass);
  if (!shaderPrg)
    return nullptr;
  shaderPrg = Setup_DefaultShader(shaderPrg, nullptr, nullptr);
  shaderPrg->SetLightingEnabled(0);
  {
    float front, back;
    front = SceneGetCurrentFrontSafe(G);
    back = SceneGetCurrentBackSafe(G);
    shaderPrg->Set1f("front", front);
    shaderPrg->Set1f("clipRange", back - front);
  }

  int width, height;
  SceneGetWidthHeightStereo(G, &width, &height);
  shaderPrg->Set2f("screenSize", width, height);

  {
    float v_scale = SceneGetScreenVertexScale(G, nullptr);
    shaderPrg->Set1f("screenOriginVertexScale", v_scale/2.f);
  }

  return shaderPrg;
}

CShaderPrg *CShaderMgr::Setup_DefaultShader(CShaderPrg * shaderPrg,
    const CSetting *set1,
    const CSetting *set2) {
  if (!shaderPrg){
    current_shader = nullptr; 
    return shaderPrg;
  }
  shaderPrg->Enable();
  shaderPrg->SetBgUniforms();
  shaderPrg->Set_Stereo_And_AnaglyphMode();

  bool two_sided_lighting_enabled = SceneGetTwoSidedLightingSettings(G, set1, set2);

  shaderPrg->SetLightingEnabled(1); // lighting on by default
  shaderPrg->Set1i("two_sided_lighting_enabled", two_sided_lighting_enabled);
  shaderPrg->Set1f("ambient_occlusion_scale", 0.f);
  shaderPrg->Set1i("accessibility_mode", SettingGetGlobal_i(G, cSetting_ambient_occlusion_mode) / 4);
  shaderPrg->Set1f("accessibility_mode_on", SettingGetGlobal_i(G, cSetting_ambient_occlusion_mode) ? 1.f : 0.f);

  // interior_color
  {
    int interior_color = SettingGet_i(G, set1, set2, cSetting_ray_interior_color);
    if (interior_color == cColorDefault || two_sided_lighting_enabled) {
      shaderPrg->Set1i("use_interior_color", 0);
    } else {
      float inter[] = { 0.f, 0.f, 0.f };
      ColorGetEncoded(G, interior_color, inter);
      shaderPrg->Set1i("use_interior_color", 1);
      shaderPrg->Set4f("interior_color", inter[0], inter[1], inter[2], 1.f);
    }
  }

  shaderPrg->Set_Specular_Values();
  shaderPrg->Set_Matrices();
  return (shaderPrg);
}
CShaderPrg *CShaderMgr::Enable_CylinderShader(RenderPass pass){
  return Enable_CylinderShader("cylinder", pass);
}

CShaderPrg *CShaderMgr::Enable_CylinderShader(const char *shader_name, RenderPass pass){
  int width, height;
  CShaderPrg *shaderPrg;

  SceneGetWidthHeightStereo(G, &width, &height);
  shaderPrg = GetShaderPrg(shader_name, 1, pass);
  if (!shaderPrg)
      return nullptr;
  shaderPrg->Enable();

  shaderPrg->SetLightingEnabled(1); // lighting on by default

  shaderPrg->Set1f("uni_radius", 0.f);

  shaderPrg->Set_Stereo_And_AnaglyphMode();

  shaderPrg->Set1f("inv_height", 1.0/height);
  shaderPrg->Set1i("no_flat_caps", 1);
  {
    float smooth_half_bonds = (SettingGetGlobal_i(G, cSetting_smooth_half_bonds)) ? .2f : 0.f;
    shaderPrg->Set1f("half_bond", smooth_half_bonds);
  }
  shaderPrg->Set_Specular_Values();
  shaderPrg->Set_Matrices();

  shaderPrg->SetBgUniforms();

  // always enable backface culling for cylinders
  glCullFace(GL_BACK);
  glEnable(GL_CULL_FACE);
  return shaderPrg;
}

CShaderPrg *CShaderMgr::Get_DefaultSphereShader(RenderPass pass){
  return GetShaderPrg("sphere", 1, pass);
}


CShaderPrg *CShaderMgr::Enable_DefaultSphereShader(RenderPass pass) {
  CShaderPrg *shaderPrg = Get_DefaultSphereShader(pass);
  if (!shaderPrg) return nullptr;
  shaderPrg->Enable();
  shaderPrg->SetLightingEnabled(1);
  shaderPrg->Set1f("sphere_size_scale", 1.f);

  shaderPrg->Set_Stereo_And_AnaglyphMode();

  shaderPrg->Set_Specular_Values();
  shaderPrg->Set_Matrices();

  shaderPrg->SetBgUniforms();

  return (shaderPrg);
}

CShaderPrg *CShaderMgr::Get_ConnectorShader(RenderPass pass){
  return GetShaderPrg("connector", 1, pass);
}

CShaderPrg *CShaderMgr::Get_DefaultShader(RenderPass pass){
  return GetShaderPrg("default", 1, pass);
}

CShaderPrg *CShaderMgr::Get_LineShader(RenderPass pass){
  return GetShaderPrg("line", 1, pass);
}

CShaderPrg *CShaderMgr::Get_SurfaceShader(RenderPass pass){
  return GetShaderPrg("surface", 1, pass);
}

CShaderPrg *CShaderMgr::Get_CylinderShader(RenderPass pass, short set_current_shader) {
  return GetShaderPrg("cylinder", set_current_shader, pass);
}

CShaderPrg *CShaderMgr::Get_CylinderNewShader(RenderPass pass, short set_current_shader) {
  return GetShaderPrg("cylinder_new", set_current_shader, pass);
}

CShaderPrg *CShaderMgr::Get_Current_Shader(){
  return current_shader;
}

CShaderPrg *CShaderMgr::Get_BackgroundShader(){
  return GetShaderPrg("bg");
}

CShaderPrg *CShaderMgr::Enable_BackgroundShader(){
  CShaderPrg * shaderPrg = Get_BackgroundShader();
  if (!shaderPrg) return shaderPrg;
  shaderPrg->Enable();

  glDisable(GL_DEPTH_TEST);
  shaderPrg->SetBgUniforms();

  return shaderPrg;
}

CShaderPrg *CShaderMgr::Enable_TriLinesShader() {
  CShaderPrg * shaderPrg = GetShaderPrg("trilines");
  if (!shaderPrg) return shaderPrg;
  shaderPrg->Enable();
  shaderPrg->SetBgUniforms();
  shaderPrg->Set_Stereo_And_AnaglyphMode();
  shaderPrg->Set_Matrices();
  {
    int width, height;
    SceneGetWidthHeightStereo(G, &width, &height);

      shaderPrg->Set2f("inv_dimensions", 1.f/width, 1.f/height);
  }
  return shaderPrg;
}

#ifndef _PYMOL_NO_AA_SHADERS
#endif

CShaderPrg *CShaderMgr::Enable_OITShader() {
  CShaderPrg * shaderPrg = GetShaderPrg("oit");
  if (!shaderPrg) return shaderPrg;
  shaderPrg->Enable();

  constexpr GLuint accumTexUnit = 5;
  constexpr GLuint revealageTexUnit = 6;
  oit_pp->activateRTAsTexture(OIT_PostProcess::OITRT::ACCUM, accumTexUnit);
  oit_pp->activateRTAsTexture(OIT_PostProcess::OITRT::REVEALAGE, revealageTexUnit);
  shaderPrg->Set1i("accumTex", accumTexUnit);
  shaderPrg->Set1i("revealageTex", revealageTexUnit);
  shaderPrg->Set1f("isRight", stereo_flag > 0 ? 1. : 0);
  glEnable(GL_BLEND);
  glBlendFuncSeparate(
      GL_SRC_ALPHA,   GL_ONE_MINUS_SRC_ALPHA,
      GL_ONE,         GL_ONE_MINUS_SRC_ALPHA);

  glDisable(GL_DEPTH_TEST);
#ifndef PURE_OPENGL_ES_2
  glDisable(GL_ALPHA_TEST);
#endif
  return shaderPrg;
}

CShaderPrg *CShaderMgr::Enable_OITCopyShader() {
  CShaderPrg * shaderPrg = GetShaderPrg("copy");
  if (!shaderPrg) return shaderPrg;
  shaderPrg->Enable();

  constexpr GLuint colorTexUnit = 7;
  activateOffscreenTexture(colorTexUnit);
  shaderPrg->Set1i("colorTex", colorTexUnit);
  if (G->ShaderMgr->stereo_blend){
    // for full-screen stereo
    glEnable(GL_BLEND);
    glBlendFunc(GL_ONE, GL_ONE);
  } else {
    glDisable(GL_BLEND);
  }
  glDisable(GL_DEPTH_TEST);
#ifndef PURE_OPENGL_ES_2
  glDisable(GL_ALPHA_TEST);
#endif
  return shaderPrg;
}

CShaderPrg *CShaderMgr::Enable_LabelShader(RenderPass pass){
  CShaderPrg *shaderPrg;
  shaderPrg = Get_LabelShader(pass);  
  if (!shaderPrg)
      return nullptr;
  shaderPrg->Enable();
  return Setup_LabelShader(shaderPrg);
}

CShaderPrg *CShaderMgr::Enable_ScreenShader(){
  CShaderPrg *shaderPrg;
  shaderPrg = Get_ScreenShader();  
  if (!shaderPrg)
      return nullptr;
  shaderPrg->Enable();

  int ortho_width, ortho_height;
  std::tie(ortho_width, ortho_height) = OrthoGetSize(*G->Ortho);
  shaderPrg->Set2f("t2PixelSize", 2.f / ortho_width, 2.f / ortho_height);

  return Setup_LabelShader(shaderPrg);
}

CShaderPrg *CShaderMgr::Enable_RampShader(){
  CShaderPrg *shaderPrg;
  shaderPrg = Get_RampShader();  
  if (!shaderPrg)
      return nullptr;
  shaderPrg->Enable();
  return Setup_LabelShader(shaderPrg);
}

CShaderPrg *CShaderMgr::Setup_LabelShader(CShaderPrg *shaderPrg) {
  int width = 0, height = 0;

  shaderPrg->Set_Matrices();

  glActiveTexture(GL_TEXTURE3);
  glBindTexture(GL_TEXTURE_2D, TextureGetTextTextureID(G));
  if (!(shaderPrg->uniform_set & 8)){
    shaderPrg->uniform_set |= 8;
    shaderPrg->Set1i("textureMap", 3);
  }

  SceneGetWidthHeightStereo(G, &width, &height);

  if (width)
    shaderPrg->Set2f("screenSize", width, height);

  shaderPrg->SetBgUniforms();

  {
    float v_scale = SceneGetScreenVertexScale(G, nullptr);
    shaderPrg->Set1f("screenOriginVertexScale", v_scale/2.f);
  }
  {
    float front, back;
    front = SceneGetCurrentFrontSafe(G);
    back = SceneGetCurrentBackSafe(G);
    shaderPrg->Set1f("front", front);
    shaderPrg->Set1f("clipRange", back - front);
  }

  return shaderPrg;
}

CShaderPrg *CShaderMgr::Get_LabelShader(RenderPass pass){
  return GetShaderPrg("label", 1, pass);
}

CShaderPrg *CShaderMgr::Get_ScreenShader() {
  if (is_picking)
    return nullptr;
  return GetShaderPrg("screen");
}

CShaderPrg *CShaderMgr::Get_RampShader() {
  return GetShaderPrg("ramp");
}

CShaderPrg *CShaderMgr::Get_IndicatorShader() {
  return GetShaderPrg("indicator");
}

CShaderPrg *CShaderMgr::Enable_IndicatorShader() {
  CShaderPrg * shaderPrg = Get_IndicatorShader();
  if (!shaderPrg) return shaderPrg;
  shaderPrg->Enable();

  shaderPrg->Set_Stereo_And_AnaglyphMode();
  shaderPrg->Set_Matrices();

  glActiveTexture(GL_TEXTURE3);
  glBindTexture(GL_TEXTURE_2D, TextureGetTextTextureID(G));
  if (!(shaderPrg->uniform_set & 8)){
    shaderPrg->Set1i("textureMap", 3);
    shaderPrg->uniform_set |= 8;
  }
#ifdef PURE_OPENGL_ES_2
  shaderPrg->SetMat4fc("g_ModelViewMatrix", SceneGetModelViewMatrix(G));
  shaderPrg->SetMat4fc("g_ProjectionMatrix", SceneGetProjectionMatrix(G));
#endif

  return (shaderPrg);
}


void CShaderMgr::ResetUniformSet() {
  for (auto & prog : programs) {
    prog.second->uniform_set = 0;
  }
}

void CShaderMgr::SetIsPicking(int is_picking) {
  this->is_picking = is_picking;
}
int CShaderMgr::GetIsPicking() {
  return is_picking;
}

#define LIGHTINGTEXTUREWIDTH 64

/*
 * Lighting setting indices are not contiguous, so we need a mapping array
 */
int light_setting_indices[] = {
  cSetting_light, cSetting_light2, cSetting_light3, cSetting_light4,
  cSetting_light5, cSetting_light6, cSetting_light7, cSetting_light8,
  cSetting_light9
};

/**
 * Generate and upload a precomputed vec2(ambient, specular) lighting texture.
 *
 * Must be equivalent to "ComputeLighting" in "compute_color_for_light.fs"
 */
void CShaderMgr::Generate_LightingTexture() {
  const int light_max = 10;
  int light_count = SettingGetGlobal_i(G, cSetting_light_count);
  int spec_count = SettingGetGlobal_i(G, cSetting_spec_count);
  float ambient = SettingGetGlobal_f(G, cSetting_ambient);
  float direct = SettingGetGlobal_f(G, cSetting_direct);
  float reflect = SettingGetGlobal_f(G, cSetting_reflect) * SceneGetReflectScaleValue(G, light_max);
  float shininess, spec_value;
  float shininess_0, spec_value_0;
  float diffuse, spec, shine;
  float power, power_0 = SettingGetGlobal_f(G, cSetting_power);
  float reflect_power = SettingGetGlobal_f(G, cSetting_reflect_power);

  float light_positions[light_max][3] = {{0.F, 0.F, 1.F}};

  // (ambient, specular) 2D texture
  unsigned char texture_AS[LIGHTINGTEXTUREWIDTH][LIGHTINGTEXTUREWIDTH][2];

  SceneGetAdjustedLightValues(G,
      &spec_value,
      &shininess,
      &spec_value_0,
      &shininess_0,
      light_max);

  if (light_count < 2) {
    light_count = 1;
    direct += reflect;
  } else if (light_count > light_max) {
    light_count = light_max;
  }

  if(spec_count < 0) {
    spec_count = light_count - 1;
  }

  for (int i = 1; i < light_count; ++i) {
    const float * setting = SettingGetGlobal_3fv(G, light_setting_indices[i - 1]);
    copy3f(setting, light_positions[i]);
    normalize3f(light_positions[i]);
    invert3f(light_positions[i]);
  }

  glGenTextures(1, &lightingTexture);

  glActiveTexture(GL_TEXTURE1);
  glBindTexture(GL_TEXTURE_CUBE_MAP, lightingTexture);
  glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

#ifndef PURE_OPENGL_ES_2
  glEnable(GL_TEXTURE_CUBE_MAP_SEAMLESS);
#else
  glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
#endif

  float normal[3];
  const float vz = LIGHTINGTEXTUREWIDTH / 2;

  for (int face = 0; face < 6; ++face) {
    for (int y = 0; y < LIGHTINGTEXTUREWIDTH; ++y) {
      for (int x = 0; x < LIGHTINGTEXTUREWIDTH; ++x) {

        float vx =   x + .5f - vz;
        float vy = -(y + .5f - vz);

        switch (face) {
          case 0: set3f(normal,  vz,  vy, -vx); break;
          case 1: set3f(normal, -vz,  vy,  vx); break;
          case 2: set3f(normal,  vx,  vz, -vy); break;
          case 3: set3f(normal,  vx, -vz,  vy); break;
          case 4: set3f(normal,  vx,  vy,  vz); break;
          case 5: set3f(normal, -vx,  vy, -vz); break;
        }

        normalize3f(normal);

        float ambient_sum = ambient;
        float specular_sum = 0.F;

        for (int i = 0; i < light_count; ++i) {
          if (i == 0) {
            diffuse = direct;
            spec = spec_value_0;
            shine = shininess_0;
            power = power_0;
          } else {
            diffuse = reflect;
            spec = spec_value;
            shine = shininess;
            power = reflect_power;
          }

          // light direction (normalized)
          const float * L = light_positions[i];

          // cosine of angle between normal and light
          float NdotL = dot_product3f(normal, L);

          // normal points away from light
          if (NdotL <= 0.0)
            continue;

          // power/reflect_power, was ray trace only until 1.7.7
          NdotL = pow(NdotL, power);

          // diffuse
          ambient_sum += NdotL * diffuse;

          // specular
          if (i <= spec_count) {
            // H = normalize(L + vec3(0., 0., 1.));
            float H[] = {0., 0., 1.};
            add3f(L, H, H);
            normalize3f(H);

            float NdotH = std::max(dot_product3f(normal, H), 0.f);
            specular_sum += spec * pow(NdotH, shine);
          }
        }

        texture_AS[y][x][0] = pymol_roundf(255.F * std::min(1.F, ambient_sum));
        texture_AS[y][x][1] = pymol_roundf(255.F * std::min(1.F, specular_sum));
      }
    }

    glTexImage2D(
        GL_TEXTURE_CUBE_MAP_POSITIVE_X + face,
        /* level */ 0,
        /* internalformat */ GL_LUMINANCE_ALPHA,
        /* width  */ LIGHTINGTEXTUREWIDTH,
        /* height */ LIGHTINGTEXTUREWIDTH,
        /* border */ 0,
        /* format */ GL_LUMINANCE_ALPHA,
        /* type */ GL_UNSIGNED_BYTE,
        /* data */ (void*) texture_AS);
  }
}

void CShaderMgr::Set_Reload_Bits(int bits){
  reload_bits |= bits;
}

void CShaderMgr::Check_Reload() {
  if(!SettingGetGlobal_b(G, cSetting_use_shaders)) {
    return;
  }

  if (reload_bits){
    if (reload_bits == RELOAD_ALL_SHADERS) {
      for (auto& prog : programs)
        prog.second->is_valid = false;
      shader_cache_processed.clear();
    }

    Reload_All_Shaders();
    reload_bits = 0;
  }
}

GLfloat *CShaderMgr::GetLineWidthRange() {
  return line_width_range;
}

#ifndef _PYMOL_NO_AA_SHADERS
#endif

/**
 * Register filename -> shader dependencies for `shader`
 */
void CShaderMgr::RegisterDependantFileNames(CShaderPrg * shader) {
  shader_deps[shader->vertfile].push_back(shader->name);
  shader_deps[shader->fragfile].push_back(shader->name);
  if (!shader->geomfile.empty())
    shader_deps[shader->geomfile].push_back(shader->name);
}

/**
 * Recursive function to insert `filename` and all the files where
 * `filename` is included into the given output vector.
 */
void CShaderMgr::CollectDependantFileNames(const std::string &filename,
    std::vector<std::string> &filenames) {
  auto it = include_deps.find(filename);
  if (it != include_deps.end()) {
    for (const char ** filenameptr = it->second;
        *filenameptr; ++filenameptr) {
      CollectDependantFileNames(*filenameptr, filenames);
    }
  }
  filenames.push_back(filename);
}

/**
 * Make derived shaders for all shaders that depend on `variable`
 */
void CShaderMgr::MakeDerivatives(const std::string &suffix, const std::string &variable) {
  std::set<std::string> shadernames;
  std::vector<std::string> filenames;

  // variable -> files
  for (const char ** filenameptr = ifdef_deps[variable];
      *filenameptr; ++filenameptr) {
    CollectDependantFileNames(*filenameptr, filenames);
  }

  // files -> shaders
  for (auto& filename : filenames) {
    auto &vec = shader_deps[filename];
    for (auto& n_it : vec) {
      shadernames.insert(n_it);
    }
  }

  // create shader derivatives
  for (const auto& shadername : shadernames) {
    auto shader = programs[shadername]->DerivativeCopy(shadername + suffix, variable);
    programs[shader->name] = shader;

    // register dependency
    RegisterDependantFileNames(shader);
  }
}

/**
 * Reload the derivative shaders for `variable`
 */
void CShaderMgr::Reload_Derivatives(const std::string &variable, bool value) {
  SetPreprocVar(variable, value, false);

  for (auto& prog : programs) {
    if (prog.second->derivative == variable)
      prog.second->reload();
  }

  SetPreprocVar(variable, !value, false);
}

/**
 * Removes `filename` and all it's parents from the shader source cache,
 * and if `invshaders` is true, also clear the `is_valid` flag for all
 * shader infos that depend on `filename`.
 */
void CShaderMgr::ShaderSourceInvalidate(const char * filename, bool invshaders) {
  // recursion for includes
  auto it = include_deps.find(filename);
  if (it != include_deps.end()) {
    for (const char ** filenameptr = it->second;
        *filenameptr; ++filenameptr) {
      ShaderSourceInvalidate(*filenameptr, invshaders);
    }
  }

  // invalidate shaders
  if (invshaders) {
    auto &vec = shader_deps[filename];
    for (const auto& shadername : vec) {
      programs[shadername]->is_valid = false;
    }
  }

  // invalidate source file
  auto repl_it = shader_cache_processed.find(filename);
  if (repl_it != shader_cache_processed.end()) {
    shader_cache_processed.erase(repl_it);
  }
}

/**
 * Set the value for the `#ifdef` variable `key` and if the value has changed,
 * then invalidate all its dependant shader source files.
 */
void CShaderMgr::SetPreprocVar(const std::string &key, bool value, bool invshaders) {
  auto &ref = preproc_vars[key];
  if (ref != value) {
    for (const char ** filenameptr = ifdef_deps[key];
        *filenameptr; ++filenameptr) {
      ShaderSourceInvalidate(*filenameptr, invshaders);
    }
    ref = value;
  }
}

/**
 * Insert `filename` -> `contents` (processed source) into the shader source
 * cache and invalidate its parents
 */
void CShaderMgr::SetShaderSource(const char * filename, const std::string &contents) {
  ShaderSourceInvalidate(filename);
  shader_cache_processed[filename] = contents;
}

void CShaderMgr::bindGPUBuffer(size_t hashid) {
  auto search = _gpu_object_map.find(hashid);
  if (search != _gpu_object_map.end())
    search->second->bind();
}

void CShaderMgr::freeGPUBuffer(size_t hashid) {
  if (!hashid)
    return;
  LOCK_GUARD_MUTEX(lock, gpu_objects_to_free_mutex);
  _gpu_objects_to_free_vector.push_back(hashid);

#ifdef _WEBGL
  freeAllGPUBuffers(); // immediate free on web
#endif
}

void CShaderMgr::freeGPUBuffers(std::vector<size_t> &&hashids) {
  LOCK_GUARD_MUTEX(lock, gpu_objects_to_free_mutex);
  _gpu_objects_to_free_vector.insert(_gpu_objects_to_free_vector.end(),
                                     hashids.begin(), hashids.end());
#ifdef _WEBGL
  freeAllGPUBuffers(); // immediate free on web
#endif
}

void CShaderMgr::freeGPUBuffers(size_t * arr, size_t len) {
  for (unsigned int i = 0; i < len; ++i)
    freeGPUBuffer(arr[i]);
}

void CShaderMgr::freeAllGPUBuffers() {
  LOCK_GUARD_MUTEX(lock, gpu_objects_to_free_mutex);
  for (auto hashid : _gpu_objects_to_free_vector) {
    auto search = _gpu_object_map.find(hashid);
    if (search != _gpu_object_map.end()) {
      if (search->second)
        delete search->second;
      _gpu_object_map.erase(search);
    }
  }
  _gpu_objects_to_free_vector.clear();
}

int CShaderMgr::GetAttributeUID(const char * name)
{
  auto uloc = attribute_uids_by_name.find(name);
  if (uloc != attribute_uids_by_name.end())
    return uloc->second;

  int uid = attribute_uids_by_name.size() + 1;
  attribute_uids_by_name[name] = uid;
  attribute_uids[uid] = name;
  return uid;
}
const char *CShaderMgr::GetAttributeName(int uid)
{
  auto uloc = attribute_uids.find(uid);
  if (uloc == attribute_uids.end())
    return nullptr;
  return attribute_uids[uid].c_str();
}

// SceneRenderBindToOffscreen
void CShaderMgr::bindOffscreen(int width, int height, GridInfo *grid) {
  using namespace tex;
  renderTarget_t::shape_type req_size(width, height);
  renderTarget_t* rt = nullptr;

#ifndef _PYMOL_NO_AA_SHADERS
#endif

  // Doesn't exist, create
  if (!offscreen_rt) {
    CGOFree(G->Scene->offscreenCGO);
    rt = newGPUBuffer<renderTarget_t>(req_size);
    rt->layout({ { 4, rt_layout_t::UBYTE } });
    offscreen_rt = rt->get_hash_id();
  } else {
    rt = getGPUBuffer<renderTarget_t>(offscreen_rt);

    // resize
    if (req_size != rt->size()) {
      rt->resize(req_size);
#ifndef _PYMOL_NO_AA_SHADERS
#endif
    }
  }

  if (rt)
    rt->bind(!stereo_blend);
  glEnable(GL_BLEND);

  SceneInitializeViewport(G, 1);
  if (grid->active) {
    grid->cur_view[0] = grid->cur_view[1] = 0;
    grid->cur_view[2] = req_size.x;
    grid->cur_view[3] = req_size.y;
  }
}

// SceneRenderBindToOffscreenOIT
void CShaderMgr::bindOffscreenOIT(int width, int height, int drawbuf) {
  using namespace tex;

  renderTarget_t::shape_type req_size(width, height);

  if(!oit_pp || oit_pp->size() != req_size) {
    oit_pp = pymol::make_unique<OIT_PostProcess>(
        width, height, getGPUBuffer<renderTarget_t>(offscreen_rt)->_rbo);
  } else {
    if (!TM3_IS_ONEBUF) {
      drawbuf = 1;
    }
    oit_pp->bindFBORBO(drawbuf - 1);
  }
}

void CShaderMgr::activateOffscreenTexture(GLuint textureIdx) {
  glActiveTexture(GL_TEXTURE0 + textureIdx);
  auto t = getGPUBuffer<renderTarget_t>(offscreen_rt);
  if (t->_textures[0])
    t->_textures[0]->bind();
}

void CShaderMgr::Disable_Current_Shader()
{
  if(current_shader){
    current_shader->Disable();
  }
}
