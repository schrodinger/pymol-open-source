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
#include "Rep.h"
#include "GenericBuffer.h"
#include "SceneDef.h"
#include "ShaderPreprocessor.h"
#include "PostProcess.h"
#include <map>
#include <set>
#include <string>
#include <unordered_map>

#ifdef _WEBGL
#define GET_FRAGDEPTH_SUPPORT() webpymol_get_fragdepth_support()
#else
#define GET_FRAGDEPTH_SUPPORT() 1
#endif

#if !defined(_WEBGL)
#include <mutex>
#define LOCK_GUARD_MUTEX(name, var) std::lock_guard<std::mutex> name(var)
#else
#define LOCK_GUARD_MUTEX(name, var)
#endif

struct GLFramebufferConfig {
  std::size_t framebuffer{};
  GLenum drawBuffer{};
};

class CShaderPrg {
public:
  struct GeometryShaderParams
  {
    std::string file = "";
    GLenum input{};
    GLenum output{};
    int numVertsOut{};
    GLuint id{};
  };

  struct TessellationShaderParams
  {
    std::string controlFile = "";
    std::string evaluationFile = "";
    GLuint controlID{};
    GLuint evaluationID{};
  };

  const std::string name, vertfile, fragfile;
  pymol::copyable_ptr<GeometryShaderParams> geomParams;
  pymol::copyable_ptr<TessellationShaderParams> tessParams;

  std::map<int, std::string> uniformLocations;

  std::string derivative;
  bool is_valid;
  bool is_linked;

  CShaderPrg(PyMOLGlobals* G_, const std::string& name,
      const std::string& vertfile, const std::string& fragfile,
      pymol::copyable_ptr<GeometryShaderParams> geomParams = nullptr,
      pymol::copyable_ptr<TessellationShaderParams> tessParams = nullptr)
      : name(name)
      , geomParams(std::move(geomParams))
      , tessParams(std::move(tessParams))
      , vertfile(vertfile)
      , fragfile(fragfile)
      , is_valid(false)
      , is_linked(false)
      , G(G_)
  {
  }

  bool reload();

  /*
   * Create a derivative copy. This will reload with ShaderMgr::Reload_Derivatives.
   */
  CShaderPrg* DerivativeCopy(
      const std::string& name, const std::string& variable)
  {
    CShaderPrg* copy =
        new CShaderPrg(G, name, vertfile, fragfile, geomParams, tessParams);
    copy->derivative = variable;
    return copy;
  }

/* Enable */
  int Enable();

/* Disable */
  int Disable();

/* Link and IsLinked */
  int Link();
  int IsLinked();

/* accessors/mutators/uniform setters */
  int Set1i(const char * name, int i);
  int Set1f(const char * name, float f);
  int Set2f(const char * name, float f1, float f2);
  int Set3f(const char * name, float f1, float f2, float f3);
  int Set4fv(const char * name, const float *f);
  int Set3fv(const char * name, const float *f);
  int Set4f(const char * name, float f1, float f2, float f3, float f4);
  int SetMat3fc(const char * name, const GLfloat * m);
  int SetMat4fc(const char * name, const GLfloat * m);

  int GetUniformLocation(const char * name);
  int GetAttribLocation(const char * name);
  void SetAttrib4fLocation(const char * name, float f1, float f2, float f3, float f4);
  void SetAttrib1fLocation(const char * name, float f1);

  int SetLightingEnabled(int);
  void SetBgUniforms();
  void Set_Stereo_And_AnaglyphMode();
  void Set_AnaglyphMode(int mode);
  void Set_Matrices();
  void Set_Specular_Values();

  void Invalidate();
  void ErrorMsgWithShaderInfoLog(const GLuint sid, const char * msg);

public:
  PyMOLGlobals * G;

  /* openGL assigned id */
  int id{};

  /* openGL fragment and vertex shader ids */
  GLuint vid{};
  GLuint fid{};

  std::map<std::string, int> uniforms;
  std::map<std::string, int> attributes;

  int uniform_set{}; // bitmask
};

/* ============================================================================
 * CShaderMgr class -- simple ShaderMgr for PyMOL
 * ============================================================================*/
class CShaderMgr {
public:
  CShaderMgr(PyMOLGlobals * G);
  ~CShaderMgr();

  void Config();
  void Set_Reload_Bits(int bits);
  void Check_Reload();
  GLfloat *GetLineWidthRange();

  template <typename T, typename... TArgs>
  T * newGPUBuffer(TArgs&&... args) {
    std::hash<gpuBuffer_t *> hash_fn;
    T * buffer = new T(std::forward<TArgs>(args)...);
    auto bufptr = dynamic_cast<gpuBuffer_t *>(buffer);
    const size_t hashid = hash_fn(bufptr);
    buffer->set_hash_id(hashid);
    _gpu_object_map[hashid] = buffer;
    return buffer;
  }

  template <typename T>
  T * getGPUBuffer(size_t hashid) {
    auto search = _gpu_object_map.find(hashid);
    if (search != _gpu_object_map.end())
      return dynamic_cast<T*>(search->second);
    else
      return nullptr;
  }

  void bindGPUBuffer(size_t hashid);
  void freeGPUBuffer(size_t handle);
  void freeGPUBuffers(std::vector<size_t> && handles);
  void freeGPUBuffers(size_t * arr, size_t len);

#ifndef _PYMOL_NO_AA_SHADERS
  bool Reload_SMAA_Shaders();
  int SMAAShadersPresent();
  CShaderPrg *Enable_FXAAShader();
  CShaderPrg *Enable_SMAA1Shader();
  CShaderPrg *Enable_SMAA2Shader();
  CShaderPrg *Enable_SMAA3Shader();
  CShaderPrg *Enable_SMAAShader(const char *shaderName);
#endif
  void Reload_Shader_Variables();

/* AddShader -- Adds to the global shader library */ 
  int AddShaderPrg(CShaderPrg * s);

/* RemoveShader -- Removes shader program by name */
  int RemoveShaderPrg(const std::string& name);

/* GetShaderPrg -- gets a ptr to the installed shader */
  CShaderPrg * GetShaderPrg(std::string name, short set_current_shader = 1, RenderPass pass = RenderPass::Antialias);

  int ShaderPrgExists(const char * name);

/* runtime check for shaders */
  int ShadersPresent();
  int GeometryShadersPresent();

  void AddVBOsToFree(GLuint *vboid, int nvbos);
  void AddVBOToFree(GLuint vboid);
  void FreeAllVBOs();

  CShaderPrg *Enable_DefaultShader(RenderPass pass);
  CShaderPrg *Enable_LineShader(RenderPass pass);
  CShaderPrg *Enable_SurfaceShader(RenderPass pass);
  CShaderPrg *Enable_DefaultShaderWithSettings(const CSetting * set1, const CSetting * set2, RenderPass pass);
  CShaderPrg *Enable_CylinderShader(const char *, RenderPass pass);
  CShaderPrg *Enable_CylinderShader(RenderPass pass);
  CShaderPrg *Enable_DefaultSphereShader(RenderPass pass);
  CShaderPrg *Enable_RampShader();
  CShaderPrg *Enable_ConnectorShader(RenderPass pass);
  CShaderPrg *Enable_TriLinesShader();
  CShaderPrg *Enable_ScreenShader();
  CShaderPrg *Enable_LabelShader(RenderPass pass);
  CShaderPrg *Enable_OITShader();
  CShaderPrg *Enable_OITCopyShader();
  CShaderPrg *Enable_IndicatorShader();
  CShaderPrg *Enable_BackgroundShader();
  CShaderPrg* Enable_BezierShader();

  void Disable_Current_Shader();

  CShaderPrg *Get_ScreenShader();
  CShaderPrg *Get_ConnectorShader(RenderPass pass);
  CShaderPrg *Get_DefaultShader(RenderPass pass);
  CShaderPrg *Get_LineShader(RenderPass pass);
  CShaderPrg *Get_SurfaceShader(RenderPass pass);
  CShaderPrg *Get_CylinderShader(RenderPass pass, short set_current_shader=1);
  CShaderPrg *Get_CylinderNewShader(RenderPass pass, short set_current_shader=1);
  CShaderPrg *Get_DefaultSphereShader(RenderPass pass);
  CShaderPrg *Get_RampShader();
  CShaderPrg *Get_Current_Shader();
  CShaderPrg *Get_IndicatorShader();
  CShaderPrg *Get_BackgroundShader();
  CShaderPrg *Get_LabelShader(RenderPass pass);
  CShaderPrg* Get_BezierShader();

  void Reload_CallComputeColorForLight();
  void Reload_All_Shaders();

  void Invalidate_All_Shaders();

  void ResetUniformSet();

  void CShaderPrg_SetIsPicking(int is_picking);

  void Generate_LightingTexture();

  std::string GetShaderSource(const std::string &filename);

  CShaderPrg *Setup_LabelShader(CShaderPrg *shaderPrg);
  CShaderPrg *Setup_DefaultShader(CShaderPrg *shaderPrg, const CSetting *set1, const CSetting *set2);

  void SetIsPicking(int is_picking);
  int GetIsPicking();

  void SetPreprocVar(const std::string &key, bool value, bool invshaders = true);

  /**
   * @return a map of shader filenames to their raw source code
   */
  static std::map<std::string, const char*>* GetRawShaderCache();

  inline static std::size_t OpenGLDefaultFramebufferID = 0;

private:
  void freeAllGPUBuffers();
  void RegisterDependantFileNames(CShaderPrg * shader);
  void CollectDependantFileNames(const std::string &filename, std::vector<std::string> &filenames);
  void MakeDerivatives(const std::string &suffix, const std::string &variable);
  void Reload_Derivatives(const std::string &variable, bool value = true);
  void ShaderSourceInvalidate(const char * filename, bool invshaders = true);
  void SetShaderSource(const char * filename, const std::string &contents);

public:
  PyMOLGlobals * G;
  int shaders_present;

#ifndef _WEBGL
  // for deleting buffers in the correct thread
  std::vector<GLuint> vbos_to_free;
  std::mutex vbos_to_free_mutex;
  std::mutex gpu_objects_to_free_mutex;
#endif

  CShaderPrg *current_shader;
  int is_picking;
  GLuint lightingTexture;

  /* These lookups and arrays are used to dynamically change the shaders
     based on preprocessor-like statements within the shaders.  These need
     string lookups. */

private:
  ShaderPreprocessor m_shaderPreprocessor;

  std::unordered_map<size_t, gpuBuffer_t*> _gpu_object_map;
  std::vector<size_t> _gpu_objects_to_free_vector;
public:
  std::map<std::string, CShaderPrg*> programs;

  std::map<int, std::string> attribute_uids;
  std::map<const std::string, int> attribute_uids_by_name;
  int GetAttributeUID(const char * name);
  const char *GetAttributeName(int);

  short print_warnings;
  int reload_bits;
  GLfloat line_width_range[2];
  short stereo_flag; /* -1 left; 0 = off; 1 = right */
  short stereo_blend;  /* 0 - no blend, 1 - blend  : for right eye stereo in full-screen e.g., anaglyph */
  bool stereo_draw_buffer_pass;

  // for glDrawBuffer (e.g. GL_BACK, unless we're using QOpenGLWidget)
  GLFramebufferConfig defaultBackbuffer{}; // Default framebuffer configuration
  GLFramebufferConfig currentFBConfig{}; // Current framebuffer configuration
  GLFramebufferConfig topLevelConfig{};  // Final render target in draw loop
                                         // (currently only known by Ortho)

private:
  bool is_configured { false };
public:
  bool IsConfigured(){ return is_configured; }
  // filename -> used by shader
  std::map<std::string, std::vector<std::string> > shader_deps;

  // Post process render targets
  std::size_t offscreen_rt { 0 }; //Texture before postprocessing;
  std::size_t offscreen_ortho_rt {};
#ifndef _PYMOL_NO_AA_SHADERS
  std::unique_ptr<PostProcess> smaa_pp;
#endif
  std::unique_ptr<PostProcess> oit_pp;

  void bindOffscreen(int width, int height, GridInfo * grid);
  void bindOffscreenOIT(int width, int height, int drawbuf = 0);

  /**
   * @brief Binds an offscreen render target for UI overlay rendering
   * @param width width of the offscreen render target
   * @param height height of the offscreen render target
   * @param clear whether to clear the offscreen render target
   */
  void bindOffscreenOrtho(int width, int height, bool clear = true);

  /**
   * @brief Reads pixel data from framebuffer
   * @param rect the rectangle to read from
   * @param srcConfig the framebuffer configuration to read from
   * @return pixel data in RGBA format
   */
  std::vector<unsigned char> readPixelsFrom(
      PyMOLGlobals* G, const Rect2D& rect, const GLFramebufferConfig& srcConfig);

  /**
   * @brief Draws pixel data to framebuffer
   * @param rect the rectangle to draw to
   * @param srcPixels the pixel data (RGBA) to draw
   * @param dstConfig the framebuffer configuration to draw to
   */
  void drawPixelsTo(PyMOLGlobals* G, const Rect2D& rect,
      const std::byte* srcPixels,
      const GLFramebufferConfig& dstConfig);

  /**
   * Activates/Binds offscreen render target.
   * @param textureIdx offset of texture unit to assign (0 for GL_TEXTURE0, 1
   * for GL_TEXTURE1, etc...)
   * @note: indices should preferably be passed in as enum or named variable for
   * clarity
   */
  void activateOffscreenTexture(GLuint textureIdx);

  /**
   * Sets the current draw buffer for rendering
   * @param drawBuffer Buffer to draw to (e.g. GL_BACK)
   * @note: This currently assumes the default framebuffer. Prefer
   * GLFramebufferConfig coverload.
   */
  void setDrawBuffer(GLenum drawBuffer);

  /**
   * Sets the current framebuffer or draw buffer for rendering
   * @param drawBuffer Framebuffer to draw to
   */
  void setDrawBuffer(const GLFramebufferConfig& config);
};

bool ShaderMgrInit(PyMOLGlobals * G);

/* for reload_bits */
enum {
  RELOAD_VARIABLES              = 0x01,
  RELOAD_CALLCOMPUTELIGHTING    = 0x02,
  RELOAD_ALL_SHADERS            = 0xff,
};

#define VERTEX_POS_SIZE    3
#define VERTEX_COLOR_SIZE  4

#define VERTEX_POS    0
#define VERTEX_NORMAL 1
#define VERTEX_COLOR  2

#define CYLINDER_VERTEX1 0
#define CYLINDER_VERTEX2 1
#define CYLINDER_COLOR   2
#define CYLINDER_COLOR2  3
#define CYLINDER_RADIUS  4
#define CYLINDER_CAP     5

#define RAMP_OFFSETPT 0

#endif

// vi:sw=2:expandtab
