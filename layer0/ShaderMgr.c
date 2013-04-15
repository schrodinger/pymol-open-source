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
#include "ListMacros.h"
#include "PyMOLOptions.h"
#include "Feedback.h"
#include "MemoryDebug.h"
#include "Setting.h"
#include "Scene.h"
#include "Color.h"
#include "Vector.h"
#include "Util.h"

#define MAX_LOG_LEN 1024

#include "ShaderText.h"

#if defined(OPENGL_ES_2)
#define DEFAULT_VS_FILENAME "default_es2.vs"
#define DEFAULT_FS_FILENAME "default_es2.fs"
#else
#define DEFAULT_VS_FILENAME "default.vs"
#define DEFAULT_FS_FILENAME "default.fs"
#endif

#define DEFAULTSCREEN_VS_FILENAME "defaultscreen.vs"
#define DEFAULTSCREEN_FS_FILENAME "defaultscreen.fs"

#define CYLINDER_VS_FILENAME "cylinder.vs"
#define CYLINDER_FS_FILENAME "cylinder.fs"
#define SPHERE_VS_FILENAME "sphere.vs"
#define SPHERE_FS_FILENAME "sphere.fs"

#define INDICATOR_VS_FILENAME "indicator.vs"
#define INDICATOR_FS_FILENAME "indicator.fs"


void getGLVersion(PyMOLGlobals * G, int *major, int* minor);
void getGLSLVersion(PyMOLGlobals * G, int* major, int* minor);

static void disableShaders(PyMOLGlobals * G);


void CShaderPrg_SetFogUniforms(PyMOLGlobals * G, CShaderPrg * shaderPrg){
  int bg_width, bg_height;
  int scene_width, scene_height;
  int ortho_width, ortho_height;
  int bg_gradient = SettingGet_b(G, NULL, NULL, cSetting_bg_gradient);

  CShaderPrg_Set1f(shaderPrg, "fogIsSolidColor", bg_gradient ? 0.f : 1.f);
  CShaderPrg_Set3fv(shaderPrg, "fogSolidColor", ColorGet(G, SettingGet_color(G, NULL, NULL, cSetting_bg_rgb)));

  SceneGetWidthHeight(G, &scene_width, &scene_height);

  OrthoGetBackgroundSize(G, &bg_width, &bg_height);
  OrthoGetSize(G, &ortho_width, &ortho_height);
  CShaderPrg_Set1f(shaderPrg, "isStretched", bg_gradient ? 1.f : 0.f);
  CShaderPrg_Set2f(shaderPrg, "viewImageSize", bg_width/(float)scene_width, bg_height/(float)scene_height);
  CShaderPrg_Set2f(shaderPrg, "pixelSize", 1.f/(float)scene_width, 1.f/(float)scene_height);
  CShaderPrg_Set2f(shaderPrg, "tPixelSize", 1.f/(float)ortho_width, 1.f/(float)ortho_height);
  CShaderPrg_Set2f(shaderPrg, "t2PixelSize", 2.f/(float)ortho_width, 2.f/(float)ortho_height);
  {
    float hpixelx = floor(scene_width / 2.f)/(float)scene_width, 
      hpixely = floor(scene_height / 2.f)/(float)scene_height;
    CShaderPrg_Set2f(shaderPrg, "halfPixel", hpixelx, hpixely);
  }
}

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

int SHADERLEX_LOOKUP(PyMOLGlobals * G, char *strarg){
  CShaderMgr *I = G->ShaderMgr;
  OVreturn_word result, result2;
  if(!OVreturn_IS_OK((result = OVLexicon_BorrowFromCString(I->ShaderLex, strarg))))
    return -1;
  result2 = OVOneToOne_GetForward(I->ShaderLexLookup, result.word);
  return result2.word;
}

void CShaderMgr_Free_Shader_Arrays(CShaderMgr *I){
  int i, sz = VLAGetSize(I->shader_replacement_strings);
  for (i=0; i<sz;i++){
    if (I->shader_replacement_strings[i]){
      VLAFreeP(I->shader_replacement_strings[i]);
      I->shader_replacement_strings[i] = 0;
    }
    I->shader_include_values[i] = 0;
  }
}

#define MIN_CHAR(x,y) (!x ? y : !y ? x : ((x < y) ? x : y))

char *CShaderPrg_ReadFromFile_Or_Use_String(PyMOLGlobals * G, char *name, char *fileName, char *fallback_str);

char *CShaderPrg_ReadFromFile_Or_Use_String_Replace_Strings(PyMOLGlobals * G, char *name, char *fileName, char *fallback_str, char **replaceStrings);

char *CShaderPrg_ReadFromFile_Or_Use_String(PyMOLGlobals * G, char *name, char *fileName, char *fallback_str){
  return CShaderPrg_ReadFromFile_Or_Use_String_Replace_Strings(G, name, fileName, fallback_str, NULL);
}

void CShaderPrg_ReplaceStringsInPlace(PyMOLGlobals *G, char *dest_line, char **replaceStrings){
  int i;
  OrthoLineType tmp_line;
  int slen, rlen;
  char *rstr;
  if (replaceStrings){
    i = 0;
    while (replaceStrings[i]){
      slen = strlen(replaceStrings[i]);
      rlen = strlen(replaceStrings[i+1]);
      while(rstr=strstr(dest_line, replaceStrings[i])){
	strcpy(tmp_line, rstr + slen);
	strcpy(rstr, replaceStrings[i+1]);
	strcpy(rstr+rlen, tmp_line);
      }
      i+=2;
    }
  }
}

void CShaderPrg_Reload_CallComputeColorForLight(PyMOLGlobals * G, char *name){
  CShaderMgr *I = G->ShaderMgr;
  int light_count = SettingGetGlobal_i(G, cSetting_light_count);
  char **reparr = Alloc(char*, 5);
  char *accstr, *tmpstr ;
  int tmpstrlen, accstrlen, i, idx;
  reparr[0] = "`light`";
  reparr[1] = "0";
  reparr[2] = "`postfix`";
  reparr[3] = "_0";
  reparr[4] = 0 ;
  accstr = CShaderPrg_ReadFromFile_Or_Use_String_Replace_Strings(G, name, "call_compute_color_for_light.fs", (char*)call_compute_color_for_light_fs, reparr);

  reparr[3] = "";
  reparr[1] = Alloc(char, 5);

  if (light_count > 8){
    PRINTFB(G, FB_Setting, FB_Warnings)
      "CShaderPrg-Error: light_count cannot be higher than 8, setting light_count to 8\n"
      ENDFB(G);
    SettingSet_i(G->Setting, cSetting_light_count, 8);
    light_count = 8;
  }
  for (i=1; i<light_count; i++){
    sprintf(reparr[1], "%d", i);
    tmpstr = CShaderPrg_ReadFromFile_Or_Use_String_Replace_Strings(G, name, "call_compute_color_for_light.fs", (char*)call_compute_color_for_light_fs, reparr);
    tmpstrlen = strlen(tmpstr);
    accstrlen = strlen(accstr);
    VLASize(accstr, char, tmpstrlen + accstrlen);
    strcpy(accstr + accstrlen-1, tmpstr);
    VLAFreeP(tmpstr);    
  }
  FreeP(reparr[1]);
  FreeP(reparr);
  idx = SHADERLEX_LOOKUP(G, "CallComputeColorForLight");
  if (I->shader_replacement_strings[idx]){
    VLAFreeP(I->shader_replacement_strings[idx]);
  }
  I->shader_replacement_strings[idx] = accstr;
}
void CShaderPrg_Reload_All_Shaders_For_CallComputeColorForLight(PyMOLGlobals * G){
  CShaderMgr_Reload_Shader_Variables(G);
  CShaderMgr_Reload_Default_Shader(G);
  CShaderMgr_Reload_Cylinder_Shader(G);
  CShaderMgr_Reload_Sphere_Shader(G);
}

void CShaderPrg_Reload_All_Shaders(PyMOLGlobals * G){
  CShaderMgr_Reload_Shader_Variables(G);
  CShaderMgr_Reload_Default_Shader(G);
  CShaderMgr_Reload_Cylinder_Shader(G);
  CShaderMgr_Reload_Sphere_Shader(G);
  CShaderMgr_Reload_Indicator_Shader(G);
}

char *CShaderPrg_ReadFromFile_Or_Use_String_Replace_Strings(PyMOLGlobals * G, char *name, char *fileName, char *fallback_str, char **replaceStrings){
  CShaderMgr *I = G->ShaderMgr;
  FILE* f = NULL;
  long size;
  char* buffer = NULL, *p, *pymol_path, *shader_path, *fullFile = NULL, *pl, *newpl, *tpl;
  size_t res;
  char *newbuffer;
  int newbuffersize;
  short allocated = 0;
  int i, len, tlen;
  OrthoLineType tmp_line, tmp_str;
  short *ifdefstack = VLAlloc(short, 10), current_include = 1;
  int ifdefstacksize = 1;
  ifdefstack[0] = 1;

  pymol_path = getenv("PYMOL_PATH");
  if (!pymol_path){
    if (I->print_warnings){
      PRINTFB(G, FB_ShaderMgr, FB_Warnings)
	" CShaderPrg_ReadFromFile_Or_Use_String: PYMOL_PATH not set, cannot read shader config files from disk\n" ENDFB(G);
    }
  } else {
    shader_path = "/data/shaders/";
    fullFile = malloc( sizeof(char) * (strlen(pymol_path)+strlen(shader_path)+strlen(fileName)+1));
    fullFile = strcpy(fullFile, pymol_path);
    fullFile = strcat(fullFile, shader_path);
    fullFile = strcat(fullFile, fileName);
    /* read the file from disk */
    f = fopen(fullFile, "rb");
  }
  if (!f) {
    if (I->print_warnings){
      PRINTFB(G, FB_ShaderMgr, FB_Errors)
	" CShaderPrg_ReadFromFile_Or_Use_String-Error: Unable to open file '%s' loading from memory\n", fullFile ENDFB(G);
    }
    buffer = fallback_str;
    res = strlen(buffer) -1;
  } else {
    fseek(f, 0, SEEK_END);
    size = ftell(f);
    fseek(f, 0, SEEK_SET);
    buffer = (char*) Alloc(char, size+255);
    allocated = 1;
    ErrChkPtr(G,buffer);
    p = buffer;
    fseek(f, 0, SEEK_SET);
    res = fread(p, size, 1, f);
    if (res) res = size;
  }
  newbuffer = VLAlloc(char, 1000);
  newbuffer[0] = 0;
  newbuffersize = 1;
  /* Now we need to read through the shader and do processing if necessary */
  pl = buffer;
  i = 1;

  while ((pl - buffer) < res){
    short pass_line = 0;
    newpl = strchr(pl, '\n');
    len = newpl - pl + 1;
    strncpy(tmp_line, pl, len);
    tmp_line[len] = 0;

    tpl = (char*)MIN_CHAR(strchr(pl, ' '), strchr(pl, '\n'));
    if (tpl <= newpl){ // && tlen < len){
      short ifl = 0, ifdefl = 0, ifdefnot = 0, elsel = 0, endifl = 0, includel = 0, lookup = 0;
      tlen = tpl - pl;
      strncpy(tmp_str, pl, tlen);
      tmp_str[tlen] = 0;
      if (!strcmp(tmp_str, "#if")){
	lookup = ifl = 1;	
      } else if (!strcmp(tmp_str, "#ifdef")){
	lookup = ifdefl = 1;
      } else if (!strcmp(tmp_str, "#ifndef")){
	lookup = ifdefl = ifdefnot = 1;
      } else if (!strcmp(tmp_str, "#else")){
	pass_line = elsel = 1;
      } else if (!strcmp(tmp_str, "#endif")){
	pass_line = endifl = 1;
      } else if (!strcmp(tmp_str, "#include")){
	lookup = includel = 1;
      }
      if (lookup){
	int off;
	char *tpl2 = (char*)MIN_CHAR(strchr(tpl + 1, '\n'), strchr(tpl + 1, ' '));
	int t2len = tpl2 - tpl - 1, is_name;
	pass_line = 1;
	strncpy(tmp_str, tpl + 1, t2len);
	tmp_str[t2len] = 0;
	off = SHADERLEX_LOOKUP(G, tmp_str);
	if (ifl){
	  /*	  char *op = (char*)MIN_CHAR(strchr(tpl2 + 1, '\n'), strchr(tpl2 + 1, ' '));
		  char *compval = (char*)MIN_CHAR(strchr(op + 1, '\n'), strchr(op + 1, ' '));
		  printf("op='%s' compval='%s'\n", op, compval); */
	} else {
	  is_name = !strcmp(tmp_str, name);
	  if (off >= 0 || is_name){
	    if (ifdefl){
	      int ifr;
	      if (is_name){
		ifr = 1;
	      } else {
		ifr = I->shader_include_values[off];
	      }
	      if (ifdefnot) ifr = !ifr;
	      VLACheck(ifdefstack, short, ifdefstacksize+1);
	      ifdefstack[ifdefstacksize++] = ifr;
	      current_include = ifr;
	    } else if (includel){
	      if (I->shader_update_when_include[off]){
		I->shader_replacement_strings[off] = CShaderPrg_ReadFromFile_Or_Use_String(G, name, I->shader_update_when_include_filename[off], (char*)I->shader_update_when_include[off]);
	      }
	      {
		int slen = strlen(I->shader_replacement_strings[off]);	    
		VLACheck(newbuffer, char, newbuffersize + slen);
		strcpy(&newbuffer[newbuffersize-1], I->shader_replacement_strings[off]);
		newbuffer[newbuffersize + slen-1] = 0;
        newbuffersize += slen;
	      }
	    }
	  } else {
	    /* Lookup doesn't exist, fails check */
	    VLACheck(ifdefstack, short, ifdefstacksize+1);
	    ifdefstack[ifdefstacksize++] = 0;
	    current_include = 0;
	  }
	}
      }
      if (endifl){
	int pl;
	ifdefstacksize--;
	pl = ifdefstacksize - 1;
	current_include = (pl >= 0) ? ifdefstack[pl] : 1;
	pass_line = 1;
      } else if (elsel){
	current_include = !current_include;
	pass_line = 1;
      }
    } 
    if (!pass_line && current_include){
      if (replaceStrings){
	CShaderPrg_ReplaceStringsInPlace(G, tmp_line, replaceStrings);
	len = strlen(tmp_line);
      }
      VLACheck(newbuffer, char, newbuffersize + len);
      strcpy(&newbuffer[newbuffersize-1], tmp_line);
      newbuffer[newbuffersize + len -1] = 0;
      newbuffersize += len;
    }
    pl = newpl + 1;
    i++;
  }
  if (allocated){
    FreeP(buffer);
  }
  VLAFreeP(ifdefstack);
  if (fullFile)
    free(fullFile);
  if (f)
    fclose(f);
  return newbuffer;
}

#define FREE_AND_REPLACE_WITH(var, with) if (var) free(var);  var = with;

void CShaderMgr_Reload_Shader_Variables(PyMOLGlobals * G){
  CShaderMgr *I = G->ShaderMgr;
  int bg_gradient = SettingGetGlobal_b(G, cSetting_bg_gradient);  
  int bg_image_mode_solid;
  int stereo, stereo_mode;
  bg_image_mode_solid = !bg_gradient;
  CShaderMgr_Free_Shader_Arrays(I);

  I->shader_include_values[SHADERLEX_LOOKUP(G, "bg_image_mode_solid")] = !bg_gradient;
  I->shader_include_values[SHADERLEX_LOOKUP(G, "bg_image_mode_stretched")] = bg_gradient;
  I->shader_include_values[SHADERLEX_LOOKUP(G, "cylinder_shader_ff_workaround")] = SettingGetGlobal_b(G, cSetting_cylinder_shader_ff_workaround);

  stereo = SettingGetGlobal_i(G, cSetting_stereo);
  stereo_mode = SettingGetGlobal_i(G, cSetting_stereo_mode);

  I->shader_include_values[SHADERLEX_LOOKUP(G, "ANAGLYPH")] = (stereo && stereo_mode==cStereo_anaglyph) ? 1 : 0;

  I->shader_replacement_strings[SHADERLEX_LOOKUP(G, "ComputeFogColor")] = CShaderPrg_ReadFromFile_Or_Use_String(G, "ComputeFogColor", "compute_fog_color.fs", (char*)compute_fog_color_fs);

  {
    int ComputeColorForLightOffset = SHADERLEX_LOOKUP(G, "ComputeColorForLight");
    FREE_AND_REPLACE_WITH(I->shader_update_when_include_filename[ComputeColorForLightOffset], strdup("compute_color_for_light.fs"));
    I->shader_update_when_include[ComputeColorForLightOffset] = (char*)compute_color_for_light_fs;
  }

  {
    int anaglyphHeaderOffset = SHADERLEX_LOOKUP(G, "ANAGLYPH_HEADER");
    FREE_AND_REPLACE_WITH(I->shader_update_when_include_filename[anaglyphHeaderOffset], strdup("anaglyph_header.fs"));
    I->shader_update_when_include[anaglyphHeaderOffset] = (char*)anaglyph_header_fs;
  }
  {
    int anaglyphOffset = SHADERLEX_LOOKUP(G, "ANAGLYPH_BODY");
    FREE_AND_REPLACE_WITH(I->shader_update_when_include_filename[anaglyphOffset], strdup("anaglyph.fs"));
    I->shader_update_when_include[anaglyphOffset] = (char*)anaglyph_fs;
  }

}


/* ============================================================================
 * ShaderMgrInit is called from PyMOL.c during start up; it just allocates
 * the global ShaderMgr
 */
OVstatus ShaderMgrInit(PyMOLGlobals * G) {
  OVreturn_word result;
  CShaderMgr *I = G->ShaderMgr = CShaderMgr_New(G);
  OVContext *C = G->Context;
  I->reload_bits = 0;
  G->ShaderMgr->is_picking = 0;

  I->ShaderLex = OVLexicon_New(C->heap);
  I->ShaderLexLookup = OVOneToOne_New(C->heap);
  
#define SHADERLEX(ARG, OFFSET)							\
  if(!OVreturn_IS_OK( (result= OVLexicon_GetFromCString(I->ShaderLex,#ARG))))  \
    return_OVstatus_FAILURE;						\
  if(!OVreturn_IS_OK( OVOneToOne_Set(I->ShaderLexLookup, result.word, OFFSET)))  \
    return_OVstatus_FAILURE;

  SHADERLEX(ComputeFogColor, 0);
  /* 1-3 reserved for incentive */
  SHADERLEX(bg_image_mode_stretched, 4);
  SHADERLEX(bg_image_mode_solid, 5);
  SHADERLEX(default_vs, 6);
  SHADERLEX(default_fs, 7);
  SHADERLEX(bg_vs, 8);
  SHADERLEX(bg_fs, 9);
  SHADERLEX(cylinder_vs, 10);
  SHADERLEX(cylinder_fs, 11);
  SHADERLEX(cylinder_shader_ff_workaround, 12);
  SHADERLEX(label_vs, 13);
  SHADERLEX(label_fs, 14);
  SHADERLEX(sphere_vs, 15);
  SHADERLEX(sphere_fs, 16);
  SHADERLEX(volume_vs, 17);
  SHADERLEX(volume_fs, 18);
  SHADERLEX(ComputeColorForLight, 19);
  SHADERLEX(CallComputeColorForLight, 20);
  /* 21 reserved for incentive */
  SHADERLEX(ANAGLYPH, 22);
  SHADERLEX(ANAGLYPH_HEADER, 23);
  SHADERLEX(ANAGLYPH_BODY, 24);
  SHADERLEX(indicator_vs, 25);
  SHADERLEX(indicator_fs, 26);
  SHADERLEX(labelscreen_vs, 27);
  SHADERLEX(labelscreen_fs, 28);
  SHADERLEX(defaultscreen_vs, 29);
  SHADERLEX(defaultscreen_fs, 30);
  SHADERLEX(screen_vs, 31);
  SHADERLEX(screen_fs, 32);
  SHADERLEX(ramp_vs, 33);
  SHADERLEX(ramp_fs, 34);
  {
    int nlexvals = 35;
    I->shader_replacement_strings = VLACalloc(char*, nlexvals);
    I->shader_include_values = VLACalloc(int, nlexvals);
    I->shader_update_when_include_filename = VLACalloc(char*, nlexvals);
    I->shader_update_when_include = VLACalloc(char*, nlexvals);
  }
  return_OVstatus_SUCCESS;
}

void CShaderPrg_BindAttribLocations(PyMOLGlobals * G, char *name){
#ifdef OPENGL_ES_2
  CShaderPrg *I = CShaderMgr_GetShaderPrg_NoSet(G->ShaderMgr, name);
  if (I){
    GLenum err ;
    glBindAttribLocation(I->id, VERTEX_POS, "a_Vertex");
    if (err = glGetError()){
      PRINTFB(G, FB_ShaderMgr, FB_Warnings) "GLERROR: a_Vertex: %d\n", err ENDFB(G);
    }
    glBindAttribLocation(I->id, VERTEX_NORMAL, "a_Normal");
    if (err = glGetError()){
      PRINTFB(G, FB_ShaderMgr, FB_Warnings) "GLERROR: a_Normal: %d\n", err ENDFB(G);
    }
    glBindAttribLocation(I->id, VERTEX_COLOR, "a_Color");
    if (err = glGetError()){
      PRINTFB(G, FB_ShaderMgr, FB_Warnings) "GLERROR: a_Color: %d\n", err ENDFB(G);
    }
    CShaderPrg_Link(I);
  }
#endif
}

void CShaderPrg_BindCylinderAttribLocations(PyMOLGlobals * G){
#ifdef OPENGL_ES_2
  CShaderPrg *I = CShaderPrg_Get_CylinderShader_NoSet(G);
  if (I){
    GLenum err ;
    glBindAttribLocation(I->id, CYLINDER_ORIGIN, "attr_origin");
    if (err = glGetError()){
      PRINTFB(G, FB_ShaderMgr, FB_Warnings) "GLERROR: attr_origin: %d\n", err ENDFB(G);
    }
    glBindAttribLocation(I->id, CYLINDER_AXIS, "attr_axis");
    if (err = glGetError()){
      PRINTFB(G, FB_ShaderMgr, FB_Warnings) "GLERROR: attr_axis: %d\n", err ENDFB(G);
    }
    glBindAttribLocation(I->id, CYLINDER_COLOR, "attr_color");
    if (err = glGetError()){
      PRINTFB(G, FB_ShaderMgr, FB_Warnings) "GLERROR: attr_color: %d\n", err ENDFB(G);
    }
    glBindAttribLocation(I->id, CYLINDER_COLOR2, "attr_color2");
    if (err = glGetError()){
      PRINTFB(G, FB_ShaderMgr, FB_Warnings) "GLERROR: attr_color2: %d\n", err ENDFB(G);
    }
    CShaderPrg_Link(I);	  
  }
#endif
}
void CShaderMgr_Reload_Sphere_Shader(PyMOLGlobals *G){
  CShaderMgr *I = G->ShaderMgr;
  char *vs, *fs;
  int vs_pl, fs_pl;
  CShaderPrg_Reload_CallComputeColorForLight(G, "sphere");
  vs_pl = SHADERLEX_LOOKUP(G, "sphere_vs");
  fs_pl = SHADERLEX_LOOKUP(G, "sphere_fs");
  vs = CShaderPrg_ReadFromFile_Or_Use_String(G, "sphere", SPHERE_VS_FILENAME, (char*)sphere_vs);
  fs = CShaderPrg_ReadFromFile_Or_Use_String(G, "sphere", SPHERE_FS_FILENAME, (char*)sphere_fs);
  if (I->shader_replacement_strings[vs_pl])
    VLAFreeP(I->shader_replacement_strings[vs_pl]);    
  if (I->shader_replacement_strings[fs_pl])
    VLAFreeP(I->shader_replacement_strings[fs_pl]);    
  I->shader_replacement_strings[vs_pl] = vs;
  I->shader_replacement_strings[fs_pl] = fs;
  CShaderPrg_Reload(G, "sphere", vs, fs);  
}

void CShaderMgr_Reload_Indicator_Shader(PyMOLGlobals *G){
  CShaderMgr *I = G->ShaderMgr;
  char *vs, *fs;
  int vs_pl, fs_pl;
  CShaderPrg_Reload_CallComputeColorForLight(G, "indicator");
  vs_pl = SHADERLEX_LOOKUP(G, "indicator_vs");
  fs_pl = SHADERLEX_LOOKUP(G, "indicator_fs");
  vs = CShaderPrg_ReadFromFile_Or_Use_String(G, "indicator", INDICATOR_VS_FILENAME, (char*)indicator_vs);
  fs = CShaderPrg_ReadFromFile_Or_Use_String(G, "indicator", INDICATOR_FS_FILENAME, (char*)indicator_fs);
  if (I->shader_replacement_strings[vs_pl])
    VLAFreeP(I->shader_replacement_strings[vs_pl]);    
  if (I->shader_replacement_strings[fs_pl])
    VLAFreeP(I->shader_replacement_strings[fs_pl]);    
  I->shader_replacement_strings[vs_pl] = vs;
  I->shader_replacement_strings[fs_pl] = fs;
  CShaderPrg_Reload(G, "indicator", vs, fs);  
}

void CShaderMgr_Reload_Default_Shader(PyMOLGlobals *G){
  CShaderMgr *I = G->ShaderMgr;
  char *vs, *fs;
  int vs_pl, fs_pl;
  CShaderPrg_Reload_CallComputeColorForLight(G, "default");
  vs_pl = SHADERLEX_LOOKUP(G, "default_vs");
  fs_pl = SHADERLEX_LOOKUP(G, "default_fs");
  vs = CShaderPrg_ReadFromFile_Or_Use_String(G, "default", DEFAULT_VS_FILENAME, (char*)default_vs);
  fs = CShaderPrg_ReadFromFile_Or_Use_String(G, "default", DEFAULT_FS_FILENAME, (char*)default_fs);
  if (I->shader_replacement_strings[vs_pl])
    VLAFreeP(I->shader_replacement_strings[vs_pl]);    
  if (I->shader_replacement_strings[fs_pl])
    VLAFreeP(I->shader_replacement_strings[fs_pl]);    
  I->shader_replacement_strings[vs_pl] = vs;
  I->shader_replacement_strings[fs_pl] = fs;
  if (CShaderPrg_Reload(G, "default", vs, fs))
      CShaderPrg_BindAttribLocations(G, "default");

  CShaderPrg_Reload_CallComputeColorForLight(G, "defaultscreen");
  vs_pl = SHADERLEX_LOOKUP(G, "defaultscreen_vs");
  fs_pl = SHADERLEX_LOOKUP(G, "defaultscreen_fs");
  vs = CShaderPrg_ReadFromFile_Or_Use_String(G, "defaultscreen", "defaultscreen.vs", (char*)defaultscreen_vs);
  fs = CShaderPrg_ReadFromFile_Or_Use_String(G, "defaultscreen", "defaultscreen.fs", (char*)defaultscreen_fs);
  if (I->shader_replacement_strings[vs_pl])
    VLAFreeP(I->shader_replacement_strings[vs_pl]);    
  if (I->shader_replacement_strings[fs_pl])
    VLAFreeP(I->shader_replacement_strings[fs_pl]);    
  I->shader_replacement_strings[vs_pl] = vs;
  I->shader_replacement_strings[fs_pl] = fs;
  if (CShaderPrg_Reload(G, "defaultscreen", vs, fs))
      CShaderPrg_BindAttribLocations(G, "defaultscreen");
}

void CShaderMgr_Reload_Cylinder_Shader(PyMOLGlobals *G){
  char *vs, *fs;
  int vs_pl, fs_pl;
  CShaderMgr *I = G->ShaderMgr;
  CShaderPrg_Reload_CallComputeColorForLight(G, "cylinder");
  vs_pl = SHADERLEX_LOOKUP(G, "cylinder_vs");
  fs_pl = SHADERLEX_LOOKUP(G, "cylinder_fs");
  vs = CShaderPrg_ReadFromFile_Or_Use_String(G, "cylinder", "cylinder.vs", (char*)cylinder_vs);
  fs = CShaderPrg_ReadFromFile_Or_Use_String(G, "cylinder", "cylinder.fs", (char*)cylinder_fs);
  if (I->shader_replacement_strings[vs_pl])
    VLAFreeP(I->shader_replacement_strings[vs_pl]);    
  if (I->shader_replacement_strings[fs_pl])
    VLAFreeP(I->shader_replacement_strings[fs_pl]);    
  I->shader_replacement_strings[vs_pl] = vs;
  I->shader_replacement_strings[fs_pl] = fs;
  CShaderPrg_Reload(G, "cylinder", vs, fs);
  CShaderPrg_BindCylinderAttribLocations(G);
}

void CShaderPrg_Update_Shaders_For_Background(PyMOLGlobals * G) {
  char *vs, *fs;
  CShaderMgr *I = G->ShaderMgr;
  CShaderMgr_Reload_Shader_Variables(G);
  if (!I)
    return;
  CShaderMgr_Reload_Default_Shader(G);

  vs = CShaderPrg_ReadFromFile_Or_Use_String(G, "bg", "bg.vs", (char*)bg_vs);
  fs = CShaderPrg_ReadFromFile_Or_Use_String(G, "bg", "bg.fs", (char*)bg_fs);
  I->shader_replacement_strings[SHADERLEX_LOOKUP(G, "bg_vs")] = vs;
  I->shader_replacement_strings[SHADERLEX_LOOKUP(G, "bg_fs")] = fs;
  CShaderPrg_Reload(G, "bg", vs, fs);  

  vs = CShaderPrg_ReadFromFile_Or_Use_String(G, "label", "label.vs", (char*)label_vs);
  fs = CShaderPrg_ReadFromFile_Or_Use_String(G, "label", "label.fs", (char*)label_fs);
  I->shader_replacement_strings[SHADERLEX_LOOKUP(G, "label_vs")] = vs;
  I->shader_replacement_strings[SHADERLEX_LOOKUP(G, "label_fs")] = fs;
  CShaderPrg_Reload(G, "label", vs, fs);  

  vs = CShaderPrg_ReadFromFile_Or_Use_String(G, "labelscreen", "labelscreen.vs", (char*)labelscreen_vs);
  fs = CShaderPrg_ReadFromFile_Or_Use_String(G, "labelscreen", "labelscreen.fs", (char*)labelscreen_fs);
  I->shader_replacement_strings[SHADERLEX_LOOKUP(G, "labelscreen_vs")] = vs;
  I->shader_replacement_strings[SHADERLEX_LOOKUP(G, "labelscreen_fs")] = fs;
  CShaderPrg_Reload(G, "labelscreen", vs, fs);

  CShaderMgr_Reload_Sphere_Shader(G);
  CShaderMgr_Reload_Cylinder_Shader(G);

  vs = CShaderPrg_ReadFromFile_Or_Use_String(G, "volume", "volume.vs", (char*)volume_vs);
  fs = CShaderPrg_ReadFromFile_Or_Use_String(G, "volume", "volume.fs", (char*)volume_fs);
  I->shader_replacement_strings[SHADERLEX_LOOKUP(G, "volume_vs")] = vs;
  I->shader_replacement_strings[SHADERLEX_LOOKUP(G, "volume_fs")] = fs;
  CShaderPrg_Reload(G, "volume", vs, fs);  

  vs = CShaderPrg_ReadFromFile_Or_Use_String(G, "indicator", "indicator.vs", (char*)indicator_vs);
  fs = CShaderPrg_ReadFromFile_Or_Use_String(G, "indicator", "indicator.fs", (char*)indicator_fs);
  I->shader_replacement_strings[SHADERLEX_LOOKUP(G, "indicator_vs")] = vs;
  I->shader_replacement_strings[SHADERLEX_LOOKUP(G, "indicator_fs")] = fs;
  CShaderPrg_Reload(G, "indicator", vs, fs);  
}

int CShaderPrg_Reload(PyMOLGlobals * G, char *name, char *v, char *f){
  int status, howLong;
  char infoLog[MAX_LOG_LEN];
  CShaderPrg * I = CShaderMgr_GetShaderPrg_NoSet(G->ShaderMgr, name);
  if (!I){
    CShaderMgr *SM = G->ShaderMgr;
    if (SM && SM->ShadersPresent){
      PRINTFB(G, FB_ShaderMgr, FB_Warnings)
	" CShaderPrg_Reload: cannot find shader '%s'\n", name ENDFB(G);
    }
    return 0;
  } else if(v) {
    if (I->v)
      free(I->v);
    I->v = strdup(v);
    glShaderSource(I->vid, 1, (const GLchar**) &I->v, NULL);
    glCompileShader((GLuint) I->vid);
    glGetShaderiv(I->vid, GL_COMPILE_STATUS, &status);
    if (!status) {
      if (G && G->Option && !G->Option->quiet) {
	PRINTFB(G, FB_ShaderMgr, FB_Errors) " CShaderPrg_Reload-Error: vertex shader compilation failed name='%s'; log follows.\n", I->name ENDFB(G);
	glGetShaderInfoLog(I->vid, MAX_LOG_LEN, &howLong, infoLog);
	PRINTFB(G, FB_ShaderMgr, FB_Errors)
	  "infoLog=%s\n", infoLog ENDFB(G);
      }
      return 0;
    }
    PRINTFB(G, FB_ShaderMgr, FB_Debugging)
      "CShaderPrg_Reload-Message: vertex shader compiled.\n" ENDFB(G);
  }
  if (f){
    if (I->f)
      free(I->f);
    I->f = strdup(f);
    glShaderSource(I->fid, 1, (const GLchar**) &I->f, NULL);
    glCompileShader((GLuint) I->fid);
    glGetShaderiv(I->fid, GL_COMPILE_STATUS, &status);
    if (!status) {
      if (G && G->Option && !G->Option->quiet) {
	PRINTFB(G, FB_ShaderMgr, FB_Errors) " CShaderPrg_Reload-Error: vertex shader compilation failed name='%s'; log follows.\n", I->name ENDFB(G);
	glGetShaderInfoLog(I->fid, MAX_LOG_LEN, &howLong, infoLog);
	PRINTFB(G, FB_ShaderMgr, FB_Errors)
	  "infoLog=%s\n", infoLog ENDFB(G);
      }
      return 0;
    }
    PRINTFB(G, FB_ShaderMgr, FB_Debugging)
      "CShaderPrg_Reload-Message: vertex shader compiled.\n" ENDFB(G);
  }
  if (v && f){
    // Link the new program 
    if (!CShaderPrg_Link(I)){
//        CShaderPrg_Delete(I);
        return 0;
    }
  }
  I->uniform_set = 0;
  return 1;
}

/* ShaderMgrConfig -- Called from PyMOL.c, configures the global ShaderMgr
 * This needs to be called once the OpenGL context has been created, it is 
 * called from MainInit() for PyMol, and from PyMOL_ConfigureShadersGL() for
 * other programs (i.e., JyMOL, AxPyMOL, etc.).
 */
void ShaderMgrConfig(PyMOLGlobals * G) {
  int major, minor;
  char buf[50];
  CShaderPrg *defaultShader, *volumeShader, *sphereShader, *defaultScreenShader,
    *cylinderShader, *spheredirectShader, *labelShader, *labelScreenShader, *indicatorShader, *bgShader, *screenShader, *rampShader;
  int hasShaders = 0, hasCylinderShader = 0;
  int ok = 0;
  GLenum err; 
  CShaderMgr *I = G->ShaderMgr;
  ok = (G && G->HaveGUI); /* && G->ValidContext); */

  if (ok) {
    err = glewInit();
  } else {
    return;
  }
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
    fprintf(stderr, " GLEW-Error: %s\n", glewGetErrorString(err));
    return;
  }

  CShaderMgr_Reload_Shader_Variables(G);
  /* First try to load configuration from files in $PYMOL_PATH, if they don't exist 
     load from const char * */
  //  FeedbackEnable(G, FB_ShaderMgr, FB_Everything);

  PRINTFB(G, FB_ShaderMgr, FB_Debugging) "reading in %s and %s\n", DEFAULT_VS_FILENAME, DEFAULT_FS_FILENAME ENDFB(G);
  CShaderPrg_Reload_CallComputeColorForLight(G, "default");
  I->shader_replacement_strings[SHADERLEX_LOOKUP(G, "default_vs")] = CShaderPrg_ReadFromFile_Or_Use_String(G, "default", DEFAULT_VS_FILENAME, (char*)default_vs);
  I->shader_replacement_strings[SHADERLEX_LOOKUP(G, "default_fs")] = CShaderPrg_ReadFromFile_Or_Use_String(G, "default", DEFAULT_FS_FILENAME, (char*)default_fs);
  defaultShader = CShaderPrg_New(G, "default", I->shader_replacement_strings[SHADERLEX_LOOKUP(G, "default_vs")],
                                               I->shader_replacement_strings[SHADERLEX_LOOKUP(G, "default_fs")]);
  if (!defaultShader){
    PRINTFB(G, FB_ShaderMgr, FB_Results)
      " PyMOLShader_NewFromFile-Warning: default shader files not found, loading from memory.\n" ENDFB(G);
    defaultShader = CShaderPrg_New(G, "default", default_vs, default_fs);
  }

  PRINTFB(G, FB_ShaderMgr, FB_Debugging) "reading in %s and %s\n", DEFAULTSCREEN_VS_FILENAME, DEFAULTSCREEN_FS_FILENAME ENDFB(G);
  CShaderPrg_Reload_CallComputeColorForLight(G, "defaultscreen");
  I->shader_replacement_strings[SHADERLEX_LOOKUP(G, "defaultscreen_vs")] = CShaderPrg_ReadFromFile_Or_Use_String(G, "defaultscreen", DEFAULTSCREEN_VS_FILENAME, (char*)defaultscreen_vs);
  I->shader_replacement_strings[SHADERLEX_LOOKUP(G, "defaultscreen_fs")] = CShaderPrg_ReadFromFile_Or_Use_String(G, "defaultscreen", DEFAULTSCREEN_FS_FILENAME, (char*)defaultscreen_fs);
  defaultScreenShader = CShaderPrg_New(G, "defaultscreen", I->shader_replacement_strings[SHADERLEX_LOOKUP(G, "defaultscreen_vs")],
				       I->shader_replacement_strings[SHADERLEX_LOOKUP(G, "defaultscreen_fs")]);
  if (!defaultScreenShader){
    PRINTFB(G, FB_ShaderMgr, FB_Results)
      " PyMOLShader_NewFromFile-Warning: defaultscreen shader files not found, loading from memory.\n" ENDFB(G);
    defaultScreenShader = CShaderPrg_New(G, "defaultscreen", defaultscreen_vs, defaultscreen_fs);
  }

  hasShaders = (defaultShader!=0);

  PRINTFB(G, FB_ShaderMgr, FB_Debugging) "reading in label.vs and label.fs\n" ENDFB(G);
  I->shader_replacement_strings[SHADERLEX_LOOKUP(G, "label_vs")] = CShaderPrg_ReadFromFile_Or_Use_String(G, "label", "label.vs", (char*)label_vs);
  I->shader_replacement_strings[SHADERLEX_LOOKUP(G, "label_fs")] = CShaderPrg_ReadFromFile_Or_Use_String(G, "label", "label.fs", (char*)label_fs);
  labelShader = CShaderPrg_New(G, "label", I->shader_replacement_strings[SHADERLEX_LOOKUP(G, "label_vs")],
                                           I->shader_replacement_strings[SHADERLEX_LOOKUP(G, "label_fs")]);
  if (labelShader){
      CShaderPrg_Link(labelShader);
      CShaderMgr_AddShaderPrg(G->ShaderMgr, labelShader);
  }

  PRINTFB(G, FB_ShaderMgr, FB_Debugging) "reading in labelscreen.vs and labelscreen.fs\n" ENDFB(G);
  I->shader_replacement_strings[SHADERLEX_LOOKUP(G, "labelscreen_vs")] = CShaderPrg_ReadFromFile_Or_Use_String(G, "labelscreen", "labelscreen.vs", (char*)labelscreen_vs);
  I->shader_replacement_strings[SHADERLEX_LOOKUP(G, "labelscreen_fs")] = CShaderPrg_ReadFromFile_Or_Use_String(G, "labelscreen", "labelscreen.fs", (char*)labelscreen_fs);
  labelScreenShader = CShaderPrg_New(G, "labelscreen", I->shader_replacement_strings[SHADERLEX_LOOKUP(G, "labelscreen_vs")],
				     I->shader_replacement_strings[SHADERLEX_LOOKUP(G, "labelscreen_fs")]);
  if (labelScreenShader){
      CShaderPrg_Link(labelScreenShader);
      CShaderMgr_AddShaderPrg(G->ShaderMgr, labelScreenShader);
  }

  PRINTFB(G, FB_ShaderMgr, FB_Debugging) "reading in screen.vs and screen.fs\n" ENDFB(G);
  I->shader_replacement_strings[SHADERLEX_LOOKUP(G, "screen_vs")] = CShaderPrg_ReadFromFile_Or_Use_String(G, "screen", "screen.vs", (char*)screen_vs);
  I->shader_replacement_strings[SHADERLEX_LOOKUP(G, "screen_fs")] = CShaderPrg_ReadFromFile_Or_Use_String(G, "screen", "screen.fs", (char*)screen_fs);
  screenShader = CShaderPrg_New(G, "screen", I->shader_replacement_strings[SHADERLEX_LOOKUP(G, "screen_vs")],
				     I->shader_replacement_strings[SHADERLEX_LOOKUP(G, "screen_fs")]);
  if (screenShader){
      CShaderPrg_Link(screenShader);
      CShaderMgr_AddShaderPrg(G->ShaderMgr, screenShader);
  }

  PRINTFB(G, FB_ShaderMgr, FB_Debugging) "reading in ramp.vs and ramp.fs\n" ENDFB(G);
  I->shader_replacement_strings[SHADERLEX_LOOKUP(G, "ramp_vs")] = CShaderPrg_ReadFromFile_Or_Use_String(G, "ramp", "ramp.vs", (char*)ramp_vs);
  I->shader_replacement_strings[SHADERLEX_LOOKUP(G, "ramp_fs")] = CShaderPrg_ReadFromFile_Or_Use_String(G, "ramp", "ramp.fs", (char*)ramp_fs);
  rampShader = CShaderPrg_New(G, "ramp", I->shader_replacement_strings[SHADERLEX_LOOKUP(G, "ramp_vs")],
				     I->shader_replacement_strings[SHADERLEX_LOOKUP(G, "ramp_fs")]);
  if (rampShader){
      CShaderPrg_Link(rampShader);
      CShaderMgr_AddShaderPrg(G->ShaderMgr, rampShader);
  }


  {
    char *vs, *fs;
    vs = CShaderPrg_ReadFromFile_Or_Use_String(G, "indicator", "indicator.vs", (char*)indicator_vs);
    fs = CShaderPrg_ReadFromFile_Or_Use_String(G, "indicator", "indicator.fs", (char*)indicator_fs);
    indicatorShader = CShaderPrg_New(G, "indicator", vs, fs);
  }
  if (indicatorShader && defaultShader){
    CShaderPrg_Link(indicatorShader);
    CShaderMgr_AddShaderPrg(G->ShaderMgr, indicatorShader);
    if (indicatorShader){
      GLenum err ;
      glBindAttribLocation(indicatorShader->id, VERTEX_POS, "a_Vertex");
      if (err = glGetError()){
	PRINTFB(G, FB_ShaderMgr, FB_Warnings) "GLERROR: a_Vertex: %d\n", err ENDFB(G);
      }
      glBindAttribLocation(indicatorShader->id, VERTEX_COLOR, "a_Color");
      if (err = glGetError()){
	PRINTFB(G, FB_ShaderMgr, FB_Warnings) "GLERROR: a_Color: %d\n", err ENDFB(G);
      }
      CShaderPrg_Link(indicatorShader);	
      CShaderMgr_AddShaderPrg(G->ShaderMgr, indicatorShader);
    }
  }

  PRINTFB(G, FB_ShaderMgr, FB_Debugging) "reading in bg.vs and bg.fs\n" ENDFB(G);
  I->shader_replacement_strings[SHADERLEX_LOOKUP(G, "bg_vs")] = CShaderPrg_ReadFromFile_Or_Use_String(G, "bg", "bg.vs", (char*)bg_vs);
  I->shader_replacement_strings[SHADERLEX_LOOKUP(G, "bg_fs")] = CShaderPrg_ReadFromFile_Or_Use_String(G, "bg", "bg.fs", (char*)bg_fs);
  bgShader = CShaderPrg_New(G, "bg", I->shader_replacement_strings[SHADERLEX_LOOKUP(G, "bg_vs")],
                                     I->shader_replacement_strings[SHADERLEX_LOOKUP(G, "bg_fs")]);
  if (bgShader){
    CShaderPrg_Link(bgShader);
    CShaderMgr_AddShaderPrg(G->ShaderMgr, bgShader);
  }

  CShaderPrg_BindAttribLocations(G, "default");
  CShaderPrg_BindAttribLocations(G, "defaultscreen");

  hasShaders = (defaultShader!=0);

  CShaderMgr_AddShaderPrg(G->ShaderMgr, defaultShader);
  CShaderMgr_AddShaderPrg(G->ShaderMgr, defaultScreenShader);

  PRINTFB(G, FB_ShaderMgr, FB_Debugging) "reading in volume.vs and volume.fs\n" ENDFB(G);
  I->shader_replacement_strings[SHADERLEX_LOOKUP(G, "volume_vs")] = CShaderPrg_ReadFromFile_Or_Use_String(G, "volume", "volume.vs", (char*)volume_vs);
  I->shader_replacement_strings[SHADERLEX_LOOKUP(G, "volume_fs")] = CShaderPrg_ReadFromFile_Or_Use_String(G, "volume", "volume.fs", (char*)volume_fs);
  volumeShader = CShaderPrg_New(G, "volume", I->shader_replacement_strings[SHADERLEX_LOOKUP(G, "volume_vs")],
                                             I->shader_replacement_strings[SHADERLEX_LOOKUP(G, "volume_fs")]);
  hasShaders &= (volumeShader!=0);

  if (volumeShader){
    CShaderMgr_AddShaderPrg(G->ShaderMgr, volumeShader);
  }

  CShaderPrg_Reload_CallComputeColorForLight(G, "sphere");
  I->shader_replacement_strings[SHADERLEX_LOOKUP(G, "sphere_vs")] = CShaderPrg_ReadFromFile_Or_Use_String(G, "sphere", SPHERE_VS_FILENAME, (char*)sphere_vs);
  I->shader_replacement_strings[SHADERLEX_LOOKUP(G, "sphere_fs")] = CShaderPrg_ReadFromFile_Or_Use_String(G, "sphere", SPHERE_FS_FILENAME, (char*)sphere_fs);
  sphereShader = CShaderPrg_New(G, "sphere", I->shader_replacement_strings[SHADERLEX_LOOKUP(G, "sphere_vs")],
                                             I->shader_replacement_strings[SHADERLEX_LOOKUP(G, "sphere_fs")]);
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
  
  CShaderPrg_Reload_CallComputeColorForLight(G, "cylinder");
  I->shader_replacement_strings[SHADERLEX_LOOKUP(G, "cylinder_vs")] = CShaderPrg_ReadFromFile_Or_Use_String(G, "cylinder", CYLINDER_VS_FILENAME, (char*)cylinder_vs);
  I->shader_replacement_strings[SHADERLEX_LOOKUP(G, "cylinder_fs")] = CShaderPrg_ReadFromFile_Or_Use_String(G, "cylinder", CYLINDER_FS_FILENAME, (char*)cylinder_fs);
  cylinderShader = CShaderPrg_New(G, "cylinder", I->shader_replacement_strings[SHADERLEX_LOOKUP(G, "cylinder_vs")],
                                                 I->shader_replacement_strings[SHADERLEX_LOOKUP(G, "cylinder_fs")]);
  hasCylinderShader &= (cylinderShader!=0);
  CShaderMgr_AddShaderPrg(G->ShaderMgr, cylinderShader);

  /* report GLSL version */
  if (G && G->Option && !G->Option->quiet) {
    getGLSLVersion(G, &major, &minor);
    sprintf(buf, " Detected GLSL version %d.%d.\n", major, minor);
    FeedbackAdd(G, buf);
  }

  G->ShaderMgr->ShadersPresent = hasShaders;

  CShaderPrg_Reload_All_Shaders(G);

  if (hasShaders) {
    SettingSetGlobal_b(G, cSetting_use_shaders, 1);
  } else {
      disableShaders(G);
  }
  I->print_warnings = 1;
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
  I->print_warnings = 1;
  return I;
}

void CShaderMgrFree(PyMOLGlobals *G){
  CShaderMgr_Delete(G->ShaderMgr);
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
  OVLexicon_DEL_AUTO_NULL(I->ShaderLex);
  OVOneToOne_Del(I->ShaderLexLookup);

  CShaderMgr_Free_Shader_Arrays(I);
  VLAFreeP(I->shader_replacement_strings);
  VLAFreeP(I->shader_include_values);

  {
    int i, sz = VLAGetSize(I->shader_update_when_include_filename);
    for (i=0; i<sz;i++){
      if (I->shader_update_when_include_filename[i]){
	free(I->shader_update_when_include_filename[i]);
	I->shader_update_when_include_filename[i] = 0;
	I->shader_update_when_include[i] = 0;
      }
    }
  }

  VLAFreeP(I->shader_update_when_include_filename);
  VLAFreeP(I->shader_update_when_include);
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

CShaderPrg * CShaderMgr_GetShaderPrgImpl(CShaderMgr * I, const char * name, short set_current_shader);

CShaderPrg * CShaderMgr_GetShaderPrg_NoSet(CShaderMgr * I, const char * name){
  return CShaderMgr_GetShaderPrgImpl(I, name, 0);
}

CShaderPrg * CShaderMgr_GetShaderPrg(CShaderMgr * I, const char * name){
  return CShaderMgr_GetShaderPrgImpl(I, name, 1);
}

CShaderPrg * CShaderMgr_GetShaderPrgImpl(CShaderMgr * I, const char * name, short set_current_shader)
{
  CShaderPrg * p = NULL, *ret = NULL;
  DListIterate(I->programs, p, next) 
    {
      if (p && strcmp(p->name,name)==0){
	ret = p;
	break;
      }
    }
  if (set_current_shader){
    I->current_shader = ret;
  }
  return ret;
}

int CShaderMgr_ShaderPrgExists(CShaderMgr * I, const char * name){
  CShaderPrg * p = NULL, *ret = NULL;
  DListIterate(I->programs, p, next) 
    {
      if (p && strcmp(p->name,name)==0){
	ret = p;
	break;
      }
    }
  return ret!=NULL;
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
	PRINTFB(G, FB_ShaderMgr, FB_Errors)
	  "shader: %s\n", I->v ENDFB(G);
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
  I->uniform_set = 0;
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
  if (p)
    p->G->ShaderMgr->current_shader = 0;
  glBindTexture(GL_TEXTURE_2D, 0);
  glActiveTexture(GL_TEXTURE0);
  return 1;
}

int CShaderPrg_DisableARB(CShaderPrg * p)
{
  glDisable(GL_FRAGMENT_PROGRAM_ARB);
  glDisable(GL_VERTEX_PROGRAM_ARB);
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
  PyMOLGlobals * G = I->G;
  int howLong;
  char infoLog[MAX_LOG_LEN];

  glLinkProgram(I->id);
  
  if (!CShaderPrg_IsLinked(I)) {
    if (G && G->Option && !G->Option->quiet) {
      GLint maxVarFloats;
      glGetIntegerv(GL_MAX_VARYING_FLOATS, &maxVarFloats);
    PRINTFB(G, FB_ShaderMgr, FB_Errors)
	" CShaderPrg_Link-Error: Shader program failed to link name='%s'; GL_MAX_VARYING_FLOATS=%d log follows.\n", I->name, maxVarFloats
	ENDFB(G);
      glGetProgramiv(I->id, GL_INFO_LOG_LENGTH, &howLong);
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

int CShaderPrg_Set2f(CShaderPrg * p, const char * name, float f1, float f2)
{
  if (p && p->id) {
    GLint loc = glGetUniformLocation(p->id, name);
    if (loc < 0)
      return 0;
    glUniform2f(loc, f1, f2);
  }
  return 1;
}

int CShaderPrg_SetMat3f_Impl(CShaderPrg * p, const char * name, GLfloat* m, GLboolean transpose) {
  if (p && p->id) {
    GLint loc = glGetUniformLocation(p->id, name);
    if (loc < 0)
      return 0;
      glUniformMatrix3fv(loc, 1, transpose, m);
  }
  return 1;
}

int CShaderPrg_SetMat4f_Impl(CShaderPrg * p, const char * name, GLfloat* m, GLboolean transpose) {
  if (p && p->id) {
    GLint loc = glGetUniformLocation(p->id, name);
    if (loc < 0)
      return 0;
    glUniformMatrix4fv(loc, 1, transpose, m);
  }
  return 1;
}

int CShaderPrg_SetMat3f(CShaderPrg * p, const char * name, GLfloat * m){
  return (CShaderPrg_SetMat3f_Impl(p, name, m, GL_TRUE));
}

int CShaderPrg_SetMat4f(CShaderPrg * p, const char * name, GLfloat * m){
  return (CShaderPrg_SetMat4f_Impl(p, name, m, GL_TRUE));
}

/*
int CShaderPrg_SetTexture2D(CShaderPrg * p, const char * name, GLuint i){
  if (p && p->id) {
    GLint loc = glGetUniformLocation(p->id, name);
    if (loc < 0)
      return 0;
    glUniform1i(loc, 1, i);
  }
  return 1;
}
*/
int CShaderPrg_SetMat3fc(CShaderPrg * p, const char * name, GLfloat * m){
  return (CShaderPrg_SetMat3f_Impl(p, name, m, GL_FALSE));
}

int CShaderPrg_SetMat4fc(CShaderPrg * p, const char * name, GLfloat * m){
  return (CShaderPrg_SetMat4f_Impl(p, name, m, GL_FALSE));
}

int CShaderPrg_Set4fv(CShaderPrg * p, const char * name, float *f){
  return (CShaderPrg_Set4f(p, name, f[0], f[1], f[2], f[3]));
}
int CShaderPrg_Set3fv(CShaderPrg * p, const char * name, float *f){
  return (CShaderPrg_Set3f(p, name, f[0], f[1], f[2]));
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

void CShaderPrg_SetAttrib1fLocation(CShaderPrg * p, const char * name, float f1){
  if (p){
    int attr = CShaderPrg_GetAttribLocation(p, name);
    if (attr>=0){
      glVertexAttrib1f(attr, f1);
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

void CShaderMgr_AddVBOsToFree(CShaderMgr * I, GLuint *vboid, int nvbos){
  int i;
  for (i=0; i<nvbos; i++){
    if (vboid[i]>0)
        CShaderMgr_AddVBOToFree(I, vboid[i]);
  }
}

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
  float spec_value = SettingGetGlobal_f(G, cSetting_specular);
  float settingSpecReflect, settingSpecDirect, settingSpecDirectPower, settingSpecPower;
  int spec_count = SettingGetGlobal_i(G, cSetting_spec_count);

  settingSpecPower = SettingGetGlobal_f(G, cSetting_spec_power);

  if(settingSpecPower < 0.0F) {
    settingSpecPower = SettingGetGlobal_f(G, cSetting_shininess);
  }

  CShaderPrg_Set1f(shaderPrg, "shininess", settingSpecPower);

  if (spec_count < 0){
    spec_count = SettingGetGlobal_i(G, cSetting_light_count);
  }
  if(spec_value == 1.0F)
    spec_value = SettingGetGlobal_f(G, cSetting_specular_intensity);

  settingSpecReflect = SettingGetGlobal_f(G, cSetting_spec_reflect);
  settingSpecReflect = SceneGetSpecularValue(G, settingSpecReflect, 10);
  settingSpecDirect = SettingGetGlobal_f(G, cSetting_spec_direct);
  settingSpecDirectPower = SettingGetGlobal_f(G, cSetting_spec_direct_power);

  if(settingSpecReflect < 0.0F)
    settingSpecReflect = spec_value;
  if(settingSpecDirect < 0.0F)
    settingSpecDirect = spec_value;
  if(settingSpecDirectPower < 0.0F)
    settingSpecDirectPower = settingSpecPower;


  if(settingSpecReflect > 1.0F)
    settingSpecReflect = 1.0F;
  if(SettingGetGlobal_f(G, cSetting_specular) < R_SMALL4) {
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
  CShaderPrg_Set1f(shaderPrg, "gamma", SettingGetGlobal_f(G, cSetting_gamma));
}

void CShaderPrg_Set_Stereo_And_AnaglyphMode(PyMOLGlobals * G, CShaderPrg * shaderPrg) {
  int stereo, stereo_mode;

  stereo = SettingGetGlobal_i(G, cSetting_stereo);
  stereo_mode = SettingGetGlobal_i(G, cSetting_stereo_mode);

  if (stereo && stereo_mode==cStereo_anaglyph){
    CShaderPrg_Set1f(shaderPrg, "stereo_flag_l", G->ShaderMgr->stereo_flag < 0 ? 1.f : 0.f);
    CShaderPrg_Set1f(shaderPrg, "stereo_flag_r", G->ShaderMgr->stereo_flag > 0 ? 1.f : 0.f);
    CShaderPrg_Set_AnaglyphMode(G, shaderPrg, SettingGetGlobal_i(G, cSetting_anaglyph_mode));
  } else {
    CShaderPrg_Set1f(shaderPrg, "stereo_flag", G->ShaderMgr->stereo_flag==0 ? 1.f : 0.f);
  }
}

int CShaderPrg_SetLightingEnabled(CShaderPrg *shaderPrg, int lighting_enabled){
  return CShaderPrg_Set1i(shaderPrg, "lighting_enabled", lighting_enabled); 
}

CShaderPrg *CShaderPrg_Enable_DefaultShaderImpl(PyMOLGlobals * G, CShaderPrg * shaderPrg, CSetting *set1, CSetting *set2);

CShaderPrg *CShaderPrg_Enable_DefaultScreenShader(PyMOLGlobals * G){
  CShaderPrg * shaderPrg = CShaderPrg_Get_DefaultScreenShader(G);
  return CShaderPrg_Enable_DefaultShaderImpl(G, shaderPrg, NULL, NULL);
}

CShaderPrg *CShaderPrg_Enable_DefaultShaderWithSettings(PyMOLGlobals * G, CSetting *set1, CSetting *set2){
  CShaderPrg * shaderPrg = CShaderPrg_Get_DefaultShader(G);
  return CShaderPrg_Enable_DefaultShaderImpl(G, shaderPrg, set1, set2);
}

CShaderPrg *CShaderPrg_Enable_DefaultShader(PyMOLGlobals * G){
  CShaderPrg * shaderPrg = CShaderPrg_Get_DefaultShader(G);
  return CShaderPrg_Enable_DefaultShaderImpl(G, shaderPrg, NULL, NULL);
}

CShaderPrg *CShaderPrg_Enable_DefaultShaderImpl(PyMOLGlobals * G, CShaderPrg * shaderPrg, CSetting *set1, CSetting *set2){
  float fog_enabled, *fog_color_top, *fog_color_bottom;
  int bg_gradient;
  if (!shaderPrg){
    G->ShaderMgr->current_shader = NULL; 
    return shaderPrg;
  }
  CShaderPrg_Enable(shaderPrg);
  fog_enabled = SettingGetGlobal_b(G, cSetting_depth_cue) ? 1.0 : 0.0;
  bg_gradient = SettingGetGlobal_b(G, cSetting_bg_gradient);
  if (bg_gradient){
    fog_color_top = ColorGet(G, SettingGet_color(G, NULL, NULL, cSetting_bg_rgb_top));
    fog_color_bottom = ColorGet(G, SettingGet_color(G, NULL, NULL, cSetting_bg_rgb_bottom));
  } else {
    fog_color_top = ColorGet(G, SettingGet_color(G, NULL, NULL, cSetting_bg_rgb));
    fog_color_bottom = fog_color_top;
  }

  CShaderPrg_SetFogUniforms(G, shaderPrg);


  glActiveTexture(GL_TEXTURE4);
  glBindTexture(GL_TEXTURE_2D, OrthoGetBackgroundTextureID(G));
  if (!(shaderPrg->uniform_set & 8)){
    CShaderPrg_Set1i(shaderPrg, "bgTextureMap", 4);
    shaderPrg->uniform_set |= 8;
  }

  CShaderPrg_Set_Stereo_And_AnaglyphMode(G, shaderPrg);

  CShaderPrg_Set1i(shaderPrg, "bg_gradient", bg_gradient);
  CShaderPrg_Set3f(shaderPrg, "fog_color_top", fog_color_top[0], fog_color_top[1], fog_color_top[2]);
  CShaderPrg_Set3f(shaderPrg, "fog_color_bottom", fog_color_bottom[0], fog_color_bottom[1], fog_color_bottom[2]);
  CShaderPrg_Set1f(shaderPrg, "fog_enabled", fog_enabled);
  
  CShaderPrg_SetLightingEnabled(shaderPrg, 1); // lighting on by default
  CShaderPrg_Set1i(shaderPrg, "two_sided_lighting_enabled", SceneGetTwoSidedLightingSettings(G, set1, set2));
  CShaderPrg_Set1i(shaderPrg, "light_count", SettingGetGlobal_i(G, cSetting_light_count));
  CShaderPrg_Set1f(shaderPrg, "ambient_occlusion_scale", 0.f);
  CShaderPrg_Set1i(shaderPrg, "accessibility_mode", SettingGetGlobal_i(G, cSetting_ambient_occlusion_mode) / 4);
  {
    int interior_color = SettingGet_i(G, set1, set2, cSetting_ray_interior_color);
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
  int ortho;
  int width, height;
  CShaderPrg *shaderPrg;
  float *m;

  SceneGetWidthHeight(G, &width, &height);
  m = SceneGetMatrix(G);
  shaderPrg = CShaderPrg_Get_CylinderShader(G);
  if (!shaderPrg)
      return NULL;
  CShaderPrg_Enable(shaderPrg);
  CShaderPrg_Set1f(shaderPrg, "uni_radius", 0.f);
  fog_enabled = SettingGetGlobal_b(G, cSetting_depth_cue) ? 1.0 : 0.0;
  bg_gradient = SettingGetGlobal_b(G, cSetting_bg_gradient);
  if (bg_gradient){
    fog_color_top = ColorGet(G, SettingGet_color(G, NULL, NULL, cSetting_bg_rgb_top));
    fog_color_bottom = ColorGet(G, SettingGet_color(G, NULL, NULL, cSetting_bg_rgb_bottom));
  } else {
    fog_color_top = ColorGet(G, SettingGet_color(G, NULL, NULL, cSetting_bg_rgb));
    fog_color_bottom = fog_color_top;
  }

  CShaderPrg_Set_Stereo_And_AnaglyphMode(G, shaderPrg);

  CShaderPrg_Set1i(shaderPrg, "bg_gradient", bg_gradient);
  CShaderPrg_Set3f(shaderPrg, "fog_color_top", fog_color_top[0], fog_color_top[1], fog_color_top[2]);
  CShaderPrg_Set3f(shaderPrg, "fog_color_bottom", fog_color_bottom[0], fog_color_bottom[1], fog_color_bottom[2]);
  CShaderPrg_Set1f(shaderPrg, "fog_enabled", fog_enabled);
  CShaderPrg_Set1f(shaderPrg, "inv_height", 1.0/height);
  ortho = SettingGetGlobal_b(G, cSetting_ortho);
  CShaderPrg_Set1f(shaderPrg, "ortho", ortho ? 1.0 : 0.0);
  CShaderPrg_Set1f(shaderPrg, "no_flat_caps", 1.0);
  CShaderPrg_Set1i(shaderPrg, "filter_front_facing", SettingGetGlobal_b(G, cSetting_cylinders_shader_filter_faces));
  CShaderPrg_Set1i(shaderPrg, "two_sided_lighting_enabled", SceneGetTwoSidedLighting(G));
  CShaderPrg_Set1i(shaderPrg, "light_count", SettingGetGlobal_i(G, cSetting_light_count));
  CShaderPrg_Set1i(shaderPrg, "filter_front_facing", SettingGetGlobal_b(G, cSetting_cylinders_shader_filter_faces));
  {
    float smooth_half_bonds = (SettingGetGlobal_i(G, cSetting_smooth_half_bonds)) ? .2f : 0.f;
    CShaderPrg_Set1f(shaderPrg, "half_bond", smooth_half_bonds);
  }
  CShaderPrg_Set_Specular_Values(G, shaderPrg);

  CShaderPrg_SetFogUniforms(G, shaderPrg);

  CShaderPrg_Set1f(shaderPrg, "fog_enabled", SettingGetGlobal_b(G, cSetting_depth_cue) ? 1.f : 0.f);
  glActiveTexture(GL_TEXTURE4);
  glBindTexture(GL_TEXTURE_2D, OrthoGetBackgroundTextureID(G));

  if (!(shaderPrg->uniform_set & 4)){
    CShaderPrg_Set1i(shaderPrg, "bgTextureMap", 4);
    shaderPrg->uniform_set |= 4;
  }
  {
    float fog[4];
    SceneSetFog(G, fog);
  }

  return shaderPrg;
}

CShaderPrg *CShaderPrg_Enable_DefaultSphereShader(PyMOLGlobals * G){
  return CShaderPrg_Enable_SphereShader(G, "sphere");
}

CShaderPrg *CShaderPrg_Get_DefaultSphereShader(PyMOLGlobals * G){
  return CShaderMgr_GetShaderPrg(G->ShaderMgr, "sphere");
}


CShaderPrg *CShaderPrg_Enable_SphereShader(PyMOLGlobals * G, char *name){
  int fog_enabled, bg_gradient;
  int ortho;
  CShaderPrg *shaderPrg;
  int width, height;
  SceneGetWidthHeight(G, &width, &height);
  shaderPrg = CShaderMgr_GetShaderPrg(G->ShaderMgr, name);
  if (!shaderPrg) return NULL;
  CShaderPrg_Enable(shaderPrg);
  CShaderPrg_SetLightingEnabled(shaderPrg, 1);
  CShaderPrg_Set1f(shaderPrg, "sphere_size_scale", 1.f);
  fog_enabled = SettingGetGlobal_b(G, cSetting_depth_cue) ? 1.0 : 0.0;
  bg_gradient = SettingGetGlobal_b(G, cSetting_bg_gradient);

  CShaderPrg_Set_Stereo_And_AnaglyphMode(G, shaderPrg);

  CShaderPrg_Set1i(shaderPrg, "bg_gradient", bg_gradient);
  CShaderPrg_Set1f(shaderPrg, "inv_height", 1.0/height);
  ortho = SettingGetGlobal_b(G, cSetting_ortho);
  CShaderPrg_Set1f(shaderPrg, "ortho", ortho ? 1.0 : 0.0);
  CShaderPrg_Set1i(shaderPrg, "light_count", SettingGetGlobal_i(G, cSetting_light_count));

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
  CShaderPrg_Set1f(shaderPrg, "fog_enabled", fog_enabled);

  CShaderPrg_SetFogUniforms(G, shaderPrg);

  glActiveTexture(GL_TEXTURE4);
  glBindTexture(GL_TEXTURE_2D, OrthoGetBackgroundTextureID(G));
  if (!(shaderPrg->uniform_set & 4)){
    CShaderPrg_Set1i(shaderPrg, "bgTextureMap", 4);
    shaderPrg->uniform_set |= 4;
  }
  {
    float fog[4];
    SceneSetFog(G, fog);
  }


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

CShaderPrg *CShaderPrg_Get_DefaultShader(PyMOLGlobals * G){
  return CShaderMgr_GetShaderPrg(G->ShaderMgr, "default");
}

CShaderPrg *CShaderPrg_Get_DefaultScreenShader(PyMOLGlobals * G){
  if (G->ShaderMgr->is_picking)
    return NULL;
  return CShaderMgr_GetShaderPrg(G->ShaderMgr, "defaultscreen");
}

CShaderPrg *CShaderPrg_Get_CylinderShaderImpl(PyMOLGlobals * G, short set_current_shader){
  CShaderPrg *shaderPrg;
  shaderPrg = CShaderMgr_GetShaderPrgImpl(G->ShaderMgr, "cylinder", set_current_shader);
  return shaderPrg;
}

CShaderPrg *CShaderPrg_Get_CylinderShader(PyMOLGlobals * G){
  return CShaderPrg_Get_CylinderShaderImpl(G, 1);
}

CShaderPrg *CShaderPrg_Get_CylinderShader_NoSet(PyMOLGlobals * G){
  return CShaderPrg_Get_CylinderShaderImpl(G, 0);
}

CShaderPrg *CShaderPrg_Get_Current_Shader(PyMOLGlobals * G){
  return G->ShaderMgr->current_shader;
}

CShaderPrg *CShaderPrg_Get_BackgroundShader(PyMOLGlobals * G){
  {
    return CShaderMgr_GetShaderPrg(G->ShaderMgr, "bg");
  }
}

CShaderPrg *CShaderPrg_Enable_BackgroundShader(PyMOLGlobals * G){
  CShaderPrg * shaderPrg = CShaderPrg_Get_BackgroundShader(G);
  if (!shaderPrg) return shaderPrg;
  CShaderPrg_Enable(shaderPrg);

  glBindTexture(GL_TEXTURE_2D, 0);
  glActiveTexture(GL_TEXTURE4);
  glBindTexture(GL_TEXTURE_2D, OrthoGetBackgroundTextureID(G));

  glDisable(GL_DEPTH_TEST);
  CShaderPrg_SetFogUniforms(G, shaderPrg);

  if (!(shaderPrg->uniform_set & 8)){
    CShaderPrg_Set1i(shaderPrg, "bgTextureMap", 4);
    shaderPrg->uniform_set |= 8;
  }
  return shaderPrg;
}

#include "Texture.h"
CShaderPrg *CShaderPrg_Enable_LabelShaderImpl(PyMOLGlobals * G, CShaderPrg *shaderPrg);

CShaderPrg *CShaderPrg_Enable_LabelScreenShader(PyMOLGlobals * G){
  CShaderPrg *shaderPrg;
  shaderPrg = CShaderPrg_Get_LabelScreenShader(G);  
  if (!shaderPrg)
      return NULL;
  CShaderPrg_Enable(shaderPrg);
  return CShaderPrg_Enable_LabelShaderImpl(G, shaderPrg);
}

CShaderPrg *CShaderPrg_Enable_LabelShader(PyMOLGlobals * G){
  CShaderPrg *shaderPrg;
  shaderPrg = CShaderPrg_Get_LabelShader(G);  
  if (!shaderPrg)
      return NULL;
  CShaderPrg_Enable(shaderPrg);
  return CShaderPrg_Enable_LabelShaderImpl(G, shaderPrg);
}

CShaderPrg *CShaderPrg_Enable_ScreenShader(PyMOLGlobals * G){
  CShaderPrg *shaderPrg;
  shaderPrg = CShaderPrg_Get_ScreenShader(G);  
  if (!shaderPrg)
      return NULL;
  CShaderPrg_Enable(shaderPrg);
  return CShaderPrg_Enable_LabelShaderImpl(G, shaderPrg);
}

CShaderPrg *CShaderPrg_Enable_RampShader(PyMOLGlobals * G){
  CShaderPrg *shaderPrg;
  shaderPrg = CShaderPrg_Get_RampShader(G);  
  if (!shaderPrg)
      return NULL;
  CShaderPrg_Enable(shaderPrg);
  return CShaderPrg_Enable_LabelShaderImpl(G, shaderPrg);
}

CShaderPrg *CShaderPrg_Enable_LabelShaderImpl(PyMOLGlobals * G, CShaderPrg *shaderPrg){
  glActiveTexture(GL_TEXTURE3);
  glBindTexture(GL_TEXTURE_2D, TextureGetTextTextureID(G));
  CShaderPrg_Set1i(shaderPrg, "textureMap", 3);
  if (!(shaderPrg->uniform_set & 8)){
    int width, height;
    SceneGetImageSizeFast(G, &width, &height);
    CShaderPrg_Set2f(shaderPrg, "screenSize", width, height);
    CShaderPrg_Set1f(shaderPrg, "aspectRatioAdjustment", 1.f);
    shaderPrg->uniform_set |= 8;
  }
  if (SceneIsGridModeActive(G)){
    int width, height, gwidth, gheight;
    float aspectRatioAdjustment;
    SceneGetImageSizeFast(G, &width, &height);
    SceneGetImageSizeFastAdjustForGrid(G, &gwidth, &gheight);
    aspectRatioAdjustment = (width/(float)height) / (gwidth/(float)gheight);
    CShaderPrg_Set1f(shaderPrg, "aspectRatioAdjustment", aspectRatioAdjustment);
  }
  CShaderPrg_Set1f(shaderPrg, "isPicking", G->ShaderMgr->is_picking ? 1.f : 0.f );

  CShaderPrg_SetFogUniforms(G, shaderPrg);

  CShaderPrg_Set1f(shaderPrg, "fog_enabled", SettingGetGlobal_b(G, cSetting_depth_cue) ? 1.f : 0.f);
  glActiveTexture(GL_TEXTURE4);
  glBindTexture(GL_TEXTURE_2D, OrthoGetBackgroundTextureID(G));
  if (!(shaderPrg->uniform_set & 4)){
    CShaderPrg_Set1i(shaderPrg, "bgTextureMap", 4);
    shaderPrg->uniform_set |= 4;
  }
  {
    float fog[4];
    SceneSetFog(G, fog);
  }
  {
    float origin[3], v_scale;
    SceneOriginGet(G, origin);
    v_scale = SceneGetScreenVertexScale(G, origin);
    CShaderPrg_Set1f(shaderPrg, "screenOriginVertexScale", v_scale/2.f);
  }
  return shaderPrg;
}

CShaderPrg *CShaderPrg_Get_LabelShader(PyMOLGlobals * G){
  return CShaderMgr_GetShaderPrg(G->ShaderMgr, "label");
}

CShaderPrg *CShaderPrg_Get_LabelScreenShader(PyMOLGlobals * G){
  if (G->ShaderMgr->is_picking)
    return NULL;
  return CShaderMgr_GetShaderPrg(G->ShaderMgr, "labelscreen");
}

CShaderPrg *CShaderPrg_Get_ScreenShader(PyMOLGlobals * G){
  if (G->ShaderMgr->is_picking)
    return NULL;
  return CShaderMgr_GetShaderPrg(G->ShaderMgr, "screen");
}

CShaderPrg *CShaderPrg_Get_RampShader(PyMOLGlobals * G){
  return CShaderMgr_GetShaderPrg(G->ShaderMgr, "ramp");
}

CShaderPrg *CShaderPrg_Get_IndicatorShader(PyMOLGlobals * G){
  return CShaderMgr_GetShaderPrg(G->ShaderMgr, "indicator");
}

CShaderPrg *CShaderPrg_Enable_IndicatorShader(PyMOLGlobals * G){
  CShaderPrg * shaderPrg = CShaderPrg_Get_IndicatorShader(G);
  if (!shaderPrg) return shaderPrg;
  CShaderPrg_Enable(shaderPrg);

  CShaderPrg_Set_Stereo_And_AnaglyphMode(G, shaderPrg);

  glActiveTexture(GL_TEXTURE3);
  glBindTexture(GL_TEXTURE_2D, TextureGetTextTextureID(G));
  if (!(shaderPrg->uniform_set & 8)){
    CShaderPrg_Set1i(shaderPrg, "textureMap", 3);
    shaderPrg->uniform_set |= 8;
  }
  return (shaderPrg);
}


void ShaderMgrResetUniformSet(PyMOLGlobals * G){
  CShaderPrg * p = NULL;
  DListIterate(G->ShaderMgr->programs, p, next) 
    {
      p->uniform_set = 0;
    }
}

void CShaderPrg_SetIsPicking(PyMOLGlobals * G, int is_picking){
  G->ShaderMgr->is_picking = is_picking;
}
int CShaderPrg_GetIsPicking(PyMOLGlobals * G){
  return G->ShaderMgr->is_picking;
}

void CShaderMgr_Set_Reload_Bits(PyMOLGlobals *G, int bits){
  CShaderMgr *I = G->ShaderMgr;
  if (I)
    I->reload_bits |= bits;
}

void CShaderMgr_Check_Reload(PyMOLGlobals *G){
  CShaderMgr *I = G->ShaderMgr;
  if (I->reload_bits){
    if (I->reload_bits & RELOAD_ALL_SHADERS){
      CShaderPrg_Reload_All_Shaders(G);
    } else {
      if (I->reload_bits & RELOAD_SHADERS_FOR_LIGHTING){
	CShaderPrg_Reload_All_Shaders_For_CallComputeColorForLight(G);
      }
      if (I->reload_bits & RELOAD_SHADERS_UPDATE_FOR_BACKGROUND){
	CShaderPrg_Update_Shaders_For_Background(G);
      }
      if (I->reload_bits & RELOAD_SHADERS_CYLINDER){
	CShaderMgr_Reload_Shader_Variables(G);
	CShaderMgr_Reload_Cylinder_Shader(G);
      }
    }
    I->reload_bits = 0;
  }
}

