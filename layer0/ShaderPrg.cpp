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
#include "ShaderMgr.h"
#include "OOMac.h"
#include "PyMOLOptions.h"
#include "Feedback.h"
#include "MemoryDebug.h"
#include "Setting.h"
#include "Scene.h"
#include "Color.h"
#include "Vector.h"
#include "Texture.h"
#include "Matrix.h"
#ifdef _WEBGL
#include "WebPyMOLLibrary.h"
#endif

const float mat3identity[] = { 1., 0., 0., 0., 1., 0., 0., 0., 1. };

// -----------------------------------------------------------------------------
// Uniforms and Attributes


int CShaderPrg::Enable() {
  if (!id) return 0;

  /* linked? */
  if (!IsLinked() && !Link()) {
    // Link prints error message
    return 0;
  }
  /* if so, use the program */
  glUseProgram(id);
  // uniform 
  Set1i("isPicking", SettingGetGlobal_b(G, cSetting_pick_shading) ? 1 : G->ShaderMgr->is_picking);
  return 1;
}

int CShaderPrg::Disable() {
  glUseProgram(0);
  G->ShaderMgr->current_shader = nullptr;
  glBindTexture(GL_TEXTURE_2D, 0);
  glActiveTexture(GL_TEXTURE0);
  return 1;
}

int CShaderPrg::IsLinked() {
  GLint status = GL_FALSE;
  if (is_linked) {
    glGetProgramiv(id, GL_LINK_STATUS, &status);
  }
  return status==GL_TRUE;
}

int CShaderPrg::Link() {
  int howLong;

  glLinkProgram(id);
  is_linked = true;
  
  if (!IsLinked()) {
    if (G && G->Option && !G->Option->quiet) {
      GLint maxVarFloats;
      int infoLogLength = 0;

#ifndef PURE_OPENGL_ES_2
      glGetIntegerv(GL_MAX_VARYING_FLOATS, &maxVarFloats);
#endif
      PRINTFB(G, FB_ShaderMgr, FB_Errors)
	" CShaderPrg_Link-Error: Shader program failed to link name='%s'; GL_MAX_VARYING_FLOATS=%d log follows.\n", name.c_str(), maxVarFloats
	ENDFB(G);
      glGetProgramiv(id, GL_INFO_LOG_LENGTH, &infoLogLength);
      if (!glGetError() && infoLogLength>0){
	char *infoLog = pymol::malloc<char>(infoLogLength);
	glGetProgramInfoLog(id, infoLogLength, &howLong, infoLog);
	PRINTFB(G, FB_ShaderMgr, FB_Errors)
	  "%s\n", infoLog
	  ENDFB(G);
	FreeP(infoLog);
      }
    }
    return 0;
  }
  uniforms.clear();
  return 1;      
}

/**
 * Get a uniform location id by name. Caches ids for faster lookup (avoiding
 * expensive `glGetUniformLocation` calls).
 */
int CShaderPrg::GetUniformLocation(const char * name)
{
  if (!id)
    return -1;

  auto uloc = uniforms.find(name);
  if (uloc != uniforms.end())
    return uloc->second;

  GLint loc = glGetUniformLocation(id, name);
  uniforms[name] = loc;
  return loc;
}

/* accessors/mutators/uniform setters */
int CShaderPrg::Set1i(const char * name, int i)
{
  GLint loc = GetUniformLocation(name);
  if (loc < 0)
    return 0;
  glUniform1i(loc, i);
  return 1;
}

int CShaderPrg::Set3f(const char * name, float f1, float f2, float f3)
{
  GLint loc = GetUniformLocation(name);
  if (loc < 0)
    return 0;
  glUniform3f(loc, f1, f2, f3);
  return 1;
}

int CShaderPrg::Set2f(const char * name, float f1, float f2)
{
  GLint loc = GetUniformLocation(name);
  if (loc < 0)
    return 0;
  glUniform2f(loc, f1, f2);
  return 1;
}

/**
 * Set column major 3x3 matrix
 */
int CShaderPrg::SetMat3fc(const char * name, const GLfloat * m){
  GLint loc = GetUniformLocation(name);
  if (loc < 0)
    return 0;
  glUniformMatrix3fv(loc, 1, GL_FALSE, m);
  return 1;
}

/**
 * Set column major 4x4 matrix
 */
int CShaderPrg::SetMat4fc(const char * name, const GLfloat * m){
  GLint loc = GetUniformLocation(name);
  if (loc < 0)
    return 0;
  glUniformMatrix4fv(loc, 1, GL_FALSE, m);
  return 1;
}

int CShaderPrg::Set4fv(const char * name, const float *f){
  return (Set4f(name, f[0], f[1], f[2], f[3]));
}
int CShaderPrg::Set3fv(const char * name, const float *f){
  return (Set3f(name, f[0], f[1], f[2]));
}

int CShaderPrg::Set4f(const char * name, float f1, float f2, float f3, float f4)
{
  GLint loc = GetUniformLocation(name);
  if (loc < 0)
    return 0;
  glUniform4f(loc, f1, f2, f3, f4);
  return 1;
}

int CShaderPrg::Set1f(const char * name, float f)
{
  GLint loc = GetUniformLocation(name);
  if (loc < 0)
    return 0;
  glUniform1f(loc, f);
  return 1;
}

int CShaderPrg::GetAttribLocation(const char * name)
{
  GLint loc = -1;

  if (id && name) {
    loc = glGetAttribLocation(id, name);
    if (loc < 0)
      return -1;
  }
  return loc;
}

void CShaderPrg::SetAttrib4fLocation(const char * name, float f1, float f2, float f3, float f4){
  if (id){
    int attr = GetAttribLocation(name);
    if (attr>=0){
      glVertexAttrib4f(attr, f1, f2, f3, f4);
    }
  }
}

void CShaderPrg::SetAttrib1fLocation(const char * name, float f1){
  if (id){
    int attr = GetAttribLocation(name);
    if (attr>=0){
      glVertexAttrib1f(attr, f1);
    }
  }
}


void CShaderPrg::Invalidate() {
  if (!id)
    return;
  if (gid){
    glDetachShader(id, gid);
    glDeleteShader(gid);
    gid = 0;
  }
  if (vid){
    glDetachShader(id, vid);
    glDeleteShader(vid);
    vid = 0;
  }
  if (fid){
    glDetachShader(id, fid);
    glDeleteShader(fid);
    fid = 0;
  }
  glDeleteProgram(id);
  id = 0;
}


void CShaderPrg::Set_Matrices() {
  if (!(uniform_set & 2) &&
      SettingGetGlobal_b(G, cSetting_precomputed_lighting)) {
    Set1i("lightingTex", 1);
    uniform_set |= 2;
  }

  const float * mvm = SceneGetModelViewMatrix(G);

  // normalmatrix = transpose(inverse(mvm))
  // Since we only support orthogonal normal matrices (e.g. with our sphere
  // shader which uses a uniform radius), scaling by the length^2 of one row
  // is sufficient.
  float normalmatrix[9];
  copy44f33f(mvm, normalmatrix);
  float len_sq = lengthsq3f(normalmatrix);
  for (int i = 0; i < 9; ++i)
    normalmatrix[i] /= len_sq;

  SetMat3fc("g_NormalMatrix", normalmatrix);
  SetMat4fc("g_ModelViewMatrix", mvm);
  SetMat4fc("g_ProjectionMatrix", SceneGetProjectionMatrix(G));
}

int CShaderPrg::SetLightingEnabled(int lighting_enabled){
  return Set1i("lighting_enabled", lighting_enabled); 
}

void CShaderPrg::Set_AnaglyphMode(int mode) {
  extern float anaglyphR_constants[6][9];
  extern float anaglyphL_constants[6][9];
  /** Coefficients from: http://3dtv.at/Knowhow/AnaglyphComparison_en.aspx */
  /** anaglyph[R|L]_constants are found in Scene.c b/c of ray tracing */
  SetMat3fc("matL", G->ShaderMgr->stereo_flag < 0 ?
                          anaglyphL_constants[mode] :
                          anaglyphR_constants[mode]);
  Set1f("gamma", SettingGetGlobal_f(G, cSetting_gamma));
}

void CShaderPrg::Set_Stereo_And_AnaglyphMode() {
  int stereo, stereo_mode;

  stereo = SettingGetGlobal_i(G, cSetting_stereo);
  stereo_mode = SettingGetGlobal_i(G, cSetting_stereo_mode);

  if (stereo && stereo_mode==cStereo_anaglyph){
    Set_AnaglyphMode(SettingGetGlobal_i(G, cSetting_anaglyph_mode));
  } else {
    SetMat3fc("matL", (GLfloat*)mat3identity);
    Set1f("gamma", 1.0);
  }
  // for using one draw buffer, switch on which_pass between color and weights
  if (TM3_IS_ONEBUF)
    Set1f("which_pass", G->ShaderMgr->stereo_draw_buffer_pass ? 1.f : 0.f);
}

void CShaderPrg::Set_Specular_Values() {
  auto trans_oblique = SettingGet<float>(G, cSetting_ray_transparency_oblique);
  if (trans_oblique > R_SMALL4) {
    Set1f("trans_oblique", trans_oblique);
    Set1f("oblique_power", SettingGetGlobal_f(G, cSetting_ray_transparency_oblique_power));
  }

  if (SettingGetGlobal_b(G, cSetting_precomputed_lighting)) {
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_CUBE_MAP, G->ShaderMgr->lightingTexture);
    return;
  }

  SceneProgramLighting(G, this);

  float settingSpecReflect, settingSpecDirect, settingSpecDirectPower, settingSpecPower;

  SceneGetAdjustedLightValues(G,
      &settingSpecReflect,
      &settingSpecPower,
      &settingSpecDirect,
      &settingSpecDirectPower);

  Set1f("spec_value_0", settingSpecDirect);
  Set1f("shininess_0", settingSpecDirectPower);
  Set1f("spec_value", settingSpecReflect);
  Set1f("shininess", settingSpecPower);
}

void CShaderPrg::SetBgUniforms() {
  auto bg_image_tilesize = SettingGet<const float*>(G, cSetting_bg_image_tilesize);
  int bg_width, bg_height;
  int scene_width, scene_height;

  Set3fv("bgSolidColor", ColorGet(G, SettingGet_color(G, nullptr, nullptr, cSetting_bg_rgb)));

  SceneGetWidthHeight(G, &scene_width, &scene_height);
  std::tie(bg_width, bg_height) = OrthoGetBackgroundSize(*G->Ortho);

  Set2f("tiledSize", bg_image_tilesize[0]/(float)scene_width, bg_image_tilesize[1]/(float)scene_height);
  Set2f("tileSize", 1.f/(float)bg_image_tilesize[0], 1.f/(float)bg_image_tilesize[1]);
  Set2f("viewImageSize", bg_width/(float)scene_width, bg_height/(float)scene_height);

  glActiveTexture(GL_TEXTURE4);
  glBindTexture(GL_TEXTURE_2D, OrthoGetBackgroundTextureID(G));
  if (!(uniform_set & 4)){
    Set1i("bgTextureMap", 4);
    uniform_set |= 4;
  }

  SceneSetFogUniforms(G, this);

  if (SettingGet<bool>(G, cSetting_chromadepth) &&
      !SettingGet<bool>(G, cSetting_orthoscopic)) {
    Set2f("clippingplanes",
        SceneGetCurrentFrontSafe(G),
        SceneGetCurrentBackSafe(G));
  }
}
