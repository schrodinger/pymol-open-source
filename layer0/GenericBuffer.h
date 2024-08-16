#pragma once
// -----------------------------------------------------------------------------
#include "GraphicsUtil.h"
#include "Vector.h"
#include <vector>
#include <tuple>
#include <map>
#include <array>
#include <type_traits>
#include <cstdlib>
#include <string.h>
#include <glm/vec2.hpp>

enum class VertexFormat {
  // 8 bit
  Byte,
  Byte2,
  Byte3,
  Byte4,
  ByteNorm,
  Byte2Norm,
  Byte3Norm,
  Byte4Norm,
  UByte,
  UByte2,
  UByte3,
  UByte4,
  UByteNorm,
  UByte2Norm,
  UByte3Norm,
  UByte4Norm,

  //Single Precision
  Float,
  Float2,
  Float3,
  Float4,

  // 32 bit
  Int,
  Int2,
  Int3,
  Int4,
  UInt,
  UInt2,
  UInt3,
  UInt4,
};

enum class VertexFormatBaseType { Byte, UByte, Float, Int, UInt };

/**
 * @param format Vertex format
 * @return The base type of the vertex format
*/
VertexFormatBaseType GetVertexFormatBaseType(VertexFormat format);

/**
 * @param format Vertex format
 * @return The GL type of the vertex format
*/
GLenum VertexFormatToGLType(VertexFormat format);

/**
 * @param format Vertex format
 * @return The number of components in the vertex format
*/
GLint VertexFormatToGLSize(VertexFormat format);

/**
 * @param format Vertex format
 * @return Whether the vertex format is a normalized type
*/
bool VertexFormatIsNormalized(VertexFormat format);

/**
 * @param format Vertex format
 * @return Whether the vertex format is a normalized type as a GLboolean
*/
GLboolean VertexFormatToGLNormalized(VertexFormat format);

/**
 * @param format Vertex format
 * @return The size of the vertex format in bytes
*/
std::size_t GetSizeOfVertexFormat(VertexFormat format);


// -----------------------------------------------------------------------------
// DESCRIPTORS
// -----------------------------------------------------------------------------
// Describes a single array held in the vbo
struct BufferDesc {
  BufferDesc(const char* _attr_name, VertexFormat _format,
      std::size_t _data_size = 0, const void* _data_ptr = nullptr,
      std::uint32_t offset = 0)
      : attr_name(_attr_name)
      , m_format(_format)
      , data_size(_data_size)
      , data_ptr(_data_ptr)
      , offset(offset)
  {
  }

  const char* attr_name{nullptr};
  VertexFormat m_format{VertexFormat::Float};
  std::size_t data_size{};
  const void* data_ptr{nullptr};
  std::uint32_t offset{};
};

using BufferDataDesc = std::vector< BufferDesc >;

/* different types of AttribOp */
enum attrib_op_type {
  NO_COPY = 0,
  FLOAT_TO_FLOAT,
  FLOAT2_TO_FLOAT2,
  FLOAT3_TO_FLOAT3,
  FLOAT4_TO_FLOAT4,
  FLOAT3_TO_UB3,
  FLOAT1_TO_UB_4TH,
  UB3_TO_UB3,
  UINT_INT_TO_PICK_DATA,
  UB1_INTERP_TO_CAP,
  FLOAT1_TO_INTERP,
  UB4_TO_UB4,
  PICK_DATA_TO_PICK_DATA,
  CYL_CAP_TO_CAP,
  FLOAT1_INTERP_TO_CAP,
  UB1_TO_INTERP,
  CYL_CAPS_ARE_ROUND,
  CYL_CAPS_ARE_FLAT,
  CYL_CAPS_ARE_CUSTOM,
  FLOAT4_TO_UB4
};

struct AttribDesc;

typedef void (*AttribOpFuncDataFunctionPtr)(void *varData, const float * pc, void *globalData, int idx);

/* AttribOpFuncData : This structure holds information a callback that sets/post-processes
                      data for a particular attribute.  Currently, a list of these functions
                      are attached to the AttribOp so that when vertices are created (i.e., incr_vertices > 0)
                      then for each vertex, this function is called for the particular attribute attribName.

   funcDataConversion - pointer to function that sets/post-processes the attribute data
   funcDataGlobalArg - pointer to global structure that can be used in each call to the callback
   attribName - attribute name this function is processing. (this calling function specifies the name, and the
                CGOConvertToShader() sets attrib to the associated AttribDesc.
 */
struct AttribOpFuncData {
  void (*funcDataConversion)(void *varData, const float * pc, void *globalData, int idx); // if set, should be called on every output value for this attribute
  void *funcDataGlobalArg;
  const char *attribName;
  AttribDesc *attrib;
  AttribOpFuncDataFunctionPtr _funcDataConversion;
  AttribOpFuncData(AttribOpFuncDataFunctionPtr _funcDataConversion,
                   void *_funcDataGlobalArg,
                   const char *_attribName)
  : funcDataConversion(_funcDataConversion), funcDataGlobalArg(_funcDataGlobalArg), attribName(_attribName), attrib(nullptr){}
};

using AttribOpFuncDataDesc = std::vector< AttribOpFuncData >;

/**
 * defines an operation that copies and (optionally) creates new vertices in 
 * a VBO operation for a particular CGO operation (op).
 *
 * op - the CGO operation
 * order - the order for this operation to be executed for the given CGO operation
 * offset - the offset into the CGO operation to copy
 * conv_type - type of copy (can be general or specific, see above, e.g. FLOAT3_TO_FLOAT3, UB1_TO_INTERP)
 * incr_vertices - the number of vertices (if any) that are generated for the VBO after this operation 
 *                 is executed.
 *
 */
struct AttribOp {
  AttribOp(unsigned short _op, size_t _order, size_t _conv_type, size_t _offset, size_t _incr_vertices=0, int _copyFromAttr=-1)
    : op(_op)
    , order(_order)
    , offset(_offset)
    , conv_type(_conv_type)
    , incr_vertices(_incr_vertices)
    , copyFromAttr(_copyFromAttr)
    {}
  unsigned short op { 0 };
  size_t order { 0 };
  size_t offset  { 0 };
  size_t conv_type { 0 };
  size_t incr_vertices { 0 };
  int copyFromAttr { -1 };
  struct AttribDesc *desc { 0 };
  struct AttribDesc *copyAttribDesc { 0 };
  std::vector<AttribOpFuncData> funcDataConversions;
};
using AttribDataOp = std::vector< AttribOp >;

/**
 * defines an attribute that is used in a shader.  this description has all of the necessary information
 * for our "optimize" function to generate either an array for input into the VBO or a call to the
 * related glVertexAttrib() call when this attribute has the same value throughout the CGO.
 *
 * attr_name - the name of this attribute inside the shaders
 * order - order of attribute used in VBO
 * attrOps - all AttribOp for this particular attribute.  This allows the user to define how this
 *           attribute gets populated from the primitive CGO's one or many CGO operations.
 * default_value - pointer to the default value of this attribute (optional, needs to be the same
 *                 size of the attribute's type)
 * repeat_value/repeat_value_length - specified if the attribute has repeating values
 *                                    repeat_value - a pointer to the type and data for repeat values
 *                                    repeat_value_length - number of repeat values
 * type_size - size of type for this attribute (e.g., GL_FLOAT, GL_UNSIGNED_BYTE)
 * type_dim - number of primitives (i.e., type_size) for each vertex of this attribute
 * data_norm - whether this attribute is normalized when passed to the VBO (GL_TRUE or GL_FALSE)
 *
 */
struct AttribDesc {
  AttribDesc(
      const char* _attr_name, VertexFormat _format, AttribDataOp _attrOps = {})
      : attr_name(_attr_name)
      , m_format(_format)
      , attrOps(_attrOps)
  {
  }
  const char * attr_name { nullptr };
  VertexFormat m_format { VertexFormat::Float };
  int order { 0 };
  AttribDataOp attrOps { };
  unsigned char *default_value { nullptr };
  unsigned char *repeat_value { nullptr };
  int repeat_value_length { 0 };
};
using AttribDataDesc = std::vector< AttribDesc >;

class gpuBuffer_t {
  friend class CShaderMgr;
public:
  virtual ~gpuBuffer_t() {};
  virtual size_t get_hash_id() { return _hashid; }
  virtual void bind() const = 0;
protected:
  virtual void set_hash_id(size_t id) { _hashid = id; }
private:
  size_t _hashid { 0 };
};

enum class buffer_layout {
  SEPARATE,   // multiple vbos
  SEQUENTIAL, // single vbo
  INTERLEAVED // single vbo
};

// -----------------------------------------------------------------------------
/* Vertexbuffer rules:
 * -----------------------------------------------------------------------------
 * - If the buffer data is interleaved then buffer sub data functionality cannot
 *   be used.
 * - The same order of buffer data must be maintained when uploading and binding
 *
 *-----------------------------------------------------------------------
 * USAGE_PATTERN:
 * SEPARATE:
 *   vbo1 [ data1 ]
 *   vbo2 [ data2 ]
 *   ...
 *   vboN [ dataN ]
 * SEQUENTIAL:
 *   vbo [ data1 | data2 | ... | dataN ]
 * INTERLEAVED:
 *   vbo [ data1[0], data2[0], ..., dataN[0] | ... | data1[M], data2[M], ..., dataN[M] ]
 */
class GenericBuffer : public gpuBuffer_t
{
  friend class CShaderMgr;

public:

  GenericBuffer(buffer_layout layout = buffer_layout::SEPARATE, GLenum usage = GL_STATIC_DRAW);
  GenericBuffer(const GenericBuffer&) = delete;
  GenericBuffer& operator=(const GenericBuffer&) = delete;
  GenericBuffer(GenericBuffer&&) = delete;
  GenericBuffer& operator=(GenericBuffer&&) = delete;
  ~GenericBuffer();

  /**
   * Conditionally generates a GPU buffer for the given data descriptor
   * @param desc The buffer data descriptor
   * @return Whether the buffer data was successfully buffered
   * @note The supplied data ptr in the struct can
   * be zero, in which case if the default usage is STATIC_DRAW then no
   * opengl buffer will be generated for that, else it is assumed that the
   * data will be supplied at a later point because it's dynamic draw.
   */
  bool bufferData(BufferDataDesc&& desc);

  /**
   * Generates a GPU buffer for the given data descriptor
   * @param desc The buffer data descriptor
   * @param data The data to buffer
   * @param len The length of the data
   * @param stride The stride of the data
   * @note assumes the data is interleaved
   */
  bool bufferData(
      BufferDataDesc&& desc, const void* data, size_t len, size_t stride);

  // -----------------------------------------------------------------------------

  /**
   * Updates (a portion of) the buffer data
   * @param offset The offset to start updating the buffer data
   * @param size The size of the data to update
   * @param data The data to update
   * @param index The index of the buffer data to update
   * @note This function assumes that the data layout hasn't changed
  */
  void bufferSubData(size_t offset, size_t size, void* data, size_t index = 0);

  /**
   * Replaces the whole interleaved buffer data
   * @param data The data to replace the buffer data with
   * @param len The length of the data
   * @param data The data to replace the buffer data with
   * @note This function assumes that the data layout hasn't changed
   */
  void bufferReplaceData(size_t offset, size_t len, const void* data);

protected:

  /**
   * Generates GPU buffer(s) for the given data descriptor
   * @return Whether the buffer data was successfully buffered
   */
  bool evaluate();

  // USAGE PATTERNS

  /**
   * Generates a separate buffer for each data descriptor
   * @return Whether the buffer data was successfully buffered
   */
  bool sepBufferData();

  /**
   * Generates a single sequential buffer for all data descriptors
   * @return Whether the buffer data was successfully buffered
   */
  bool seqBufferData();

  /**
   * Generates a single interleaved buffer for all data descriptors
   * @return Whether the buffer data was successfully buffered
   */
  bool interleaveBufferData();

  /**
   * Generates an OpenGL buffer with given size and data
   * @param id The OpenGL buffer ID
   * @param size The size of the buffer
   * @param ptr The data to buffer
   * @return Whether the buffer was successfully generated
   */
  bool genBuffer(GLuint& id, size_t size, const void* ptr);

  /**
   * @return OpenGL buffer type
   */
  virtual GLenum bufferType() const = 0;

protected:
  bool m_status{false};
  bool m_interleaved{false};
  GLuint m_interleavedID{0};
  const GLenum m_buffer_usage{GL_STATIC_DRAW};
  const buffer_layout m_layout{ buffer_layout::SEPARATE };
  size_t m_stride{0};
  BufferDataDesc m_desc;
  std::vector<GLuint> desc_glIDs; // m_desc's gl buffer IDs
};

/**
 * Vertex buffer specialization
 */
class VertexBuffer : public GenericBuffer {
  void bind_attrib(GLuint prg, const BufferDesc& d, GLuint glID);

public:
  VertexBuffer(buffer_layout layout = buffer_layout::SEPARATE,
      GLenum usage = GL_STATIC_DRAW);

  void bind() const;

  void bind(GLuint prg, int index = -1);

  void unbind();

  void maskAttributes(std::vector<GLint> attrib_locs);
  void maskAttribute(GLint attrib_loc);

  GLenum bufferType() const override;

private:
  // m_locs is only for interleaved data
  std::vector<GLint> m_locs;
  std::vector<GLint> m_attribmask;
};

/**
 * Index buffer specialization
 */
class IndexBuffer : public GenericBuffer {
public:
  using GenericBuffer::GenericBuffer;

  void bind() const;
  void unbind();
  GLenum bufferType() const override;
};

// Forward Decls
class frameBuffer_t;
class renderBuffer_t;

/***********************************************************************
 * RENDERBUFFER
 ***********************************************************************/
namespace rbo {
  enum storage {
    DEPTH16 = 0,
    DEPTH24,
    COUNT
  };

  void unbind();
};

class renderBuffer_t : public gpuBuffer_t {
  friend class frameBuffer_t;
  friend class CShaderMgr;
public:
  renderBuffer_t(int width, int height, rbo::storage storage) :
    _width(width), _height(height), _storage(storage) {
    genBuffer();
  }
  ~renderBuffer_t() {
    freeBuffer();
  }

  void bind() const;
  void unbind() const;

private:
  void genBuffer();
  void freeBuffer();

protected:
  uint32_t     _id;
  int          _dim[2];
  int          _width;
  int          _height;
  rbo::storage _storage;
};

/***********************************************************************
 * TEXTURE
 ***********************************************************************/
namespace tex {
  enum class dim : int {
    D1 = 0,
    D2,
    D3,
    COUNT
  };
  enum class format : int {
    R = (int)dim::COUNT,
    RG,
    RGB,
    RGBA,
    COUNT
  };
  enum class data_type : int {
    UBYTE = (int)format::COUNT,
    FLOAT,
    HALF_FLOAT,
    COUNT
  };
  enum class filter : int {
    NEAREST = (int)data_type::COUNT,
    LINEAR,
    NEAREST_MIP_NEAREST,
    NEAREST_MIP_LINEAR,
    LINEAR_MIP_NEAREST,
    LINEAR_MIP_LINEAR,
    COUNT
  };
  enum class wrap : int {
    REPEAT = (int)filter::COUNT,
    CLAMP,
    MIRROR_REPEAT,
    CLAMP_TO_EDGE,
    CLAMP_TO_BORDER,
    MIRROR_CLAMP_TO_EDGE,
    COUNT
  };
  enum class env_name : int {
    ENV_MODE = (int)wrap::COUNT,
    COUNT
  };
  enum class env_param : int {
    REPLACE = (int)env_name::COUNT,
    COUNT
  };
  const uint32_t max_params = (int)env_param::COUNT;

  void env(tex::env_name, tex::env_param);
};

class textureBuffer_t : public gpuBuffer_t {
  friend class frameBuffer_t;
public:
  // Generates a 1D texture
  textureBuffer_t(tex::format format, tex::data_type type,
                  tex::filter mag, tex::filter min,
                  tex::wrap wrap_s) :
    _dim(tex::dim::D1), _format(format), _type(type),
    _sampling({(int)mag, (int)min, (int)wrap_s, 0, 0})
    {
      genBuffer();
    };
  // Generates a 2D texture
  textureBuffer_t(tex::format format, tex::data_type type,
                  tex::filter mag, tex::filter min,
                  tex::wrap wrap_s, tex::wrap wrap_t) :
    _dim(tex::dim::D2), _format(format), _type(type),
    _sampling({(int)mag, (int)min, (int)wrap_s, (int)wrap_t, 0})
    {
      genBuffer();
    };
  // Generates a 3D texture
  textureBuffer_t(tex::format format, tex::data_type type,
                  tex::filter mag, tex::filter min,
                  tex::wrap wrap_s, tex::wrap wrap_t,
                  tex::wrap wrap_r) :
    _dim(tex::dim::D3), _format(format), _type(type),
    _sampling({(int)mag, (int)min, (int)wrap_s, (int)wrap_t, (int)wrap_r})
    {
      genBuffer();
    };
  ~textureBuffer_t() {
    freeBuffer();
  }

  void bind() const;
  /**
   * Binds the texture to a specific texture unit
   * @param textureUnit The texture unit to bind to (0, 1, 2, etc)
   */
  void bindToTextureUnit(std::uint8_t textureUnit) const;
  void unbind() const;

  void texture_data_1D(int width, const void * data);
  void texture_data_2D(int width, int height, const void * data);
  void texture_data_3D(int width, int height, int depth, const void * data);

  /**
   * Specifies a 2D texture subimage
   * @param xoffset The x offset of the subimage
   * @param yoffset The y offset of the subimage
   * @param width The width of the subimage
   * @param height The height of the subimage
   * @param data The data to upload
  */
  void texture_subdata_2D(int xoffset, int yoffset, int width, int height, const void* data);
private:
  void genBuffer();
  void freeBuffer();

private:
  const tex::dim           _dim;
  const tex::format        _format;
  const tex::data_type     _type;
  const std::array<int, 5> _sampling;
  uint32_t                 _id     { 0 };
  int                      _width  { 0 };
  int                      _height { 0 };
  int                      _depth  { 0 };
};

/***********************************************************************
 * FRAMEBUFFER
 ***********************************************************************/
namespace fbo {
  enum attachment {
    COLOR0 = 0,
    COLOR1,
    COLOR2,
    COLOR3,
    DEPTH,
    COUNT
  };

  // global unbind for fbos
  void unbind();
}

class frameBuffer_t : public gpuBuffer_t {
  friend class CShaderMgr;
public:
  frameBuffer_t() {
    genBuffer();
  }
  ~frameBuffer_t() {
    freeBuffer();
  }

  void attach_texture(textureBuffer_t * texture, fbo::attachment loc);
  void attach_renderbuffer(renderBuffer_t * renderbuffer, fbo::attachment loc);
  void print_fbo();

  void bind() const;
  void unbind() const;
private:
  void genBuffer();
  void freeBuffer();
  void checkStatus();

protected:
  uint32_t _id { 0 };
  std::vector<std::tuple<size_t, fbo::attachment>> _attachments;
};


/***********************************************************************
 * RENDERTARGET
 *----------------------------------------------------------------------
 * A 2D render target that automatically has depth, used for postprocess
 ***********************************************************************/
struct rt_layout_t {
  enum data_type   { UBYTE, FLOAT };
  rt_layout_t(uint8_t _nchannels, data_type _type)
    : nchannels(_nchannels), type(_type) {}
  rt_layout_t(uint8_t _nchannels, data_type _type, int _width, int _height)
    : nchannels(_nchannels), type(_type), width(_width), height(_height) {}
  uint8_t     nchannels;
  data_type   type;
  int         width  { 0 };
  int         height { 0 };
};

class renderTarget_t : public gpuBuffer_t {
  friend class CShaderMgr;
public:
  using shape_type = glm::ivec2;

  renderTarget_t(shape_type size) : _size(size) {}
  renderTarget_t(int width, int height) : _size(width, height) {}
  ~renderTarget_t();

  void bind() const { bind(true); };
  void bind(bool clear) const;
  void bindFBORBO() const
  {
    _fbo->bind();
    _rbo->bind();
  }

  void layout(std::vector<rt_layout_t>&& desc, renderBuffer_t * with_rbo = nullptr);
  void resize(shape_type size);

  const shape_type& size() const { return _size; };

  renderBuffer_t* rbo() const noexcept { return _rbo; }
  frameBuffer_t* fbo() const noexcept { return _fbo; }
  const std::vector<textureBuffer_t*>& textures() const noexcept { return _textures; }

protected:
  bool _shared_rbo { false };
  shape_type _size;
  frameBuffer_t * _fbo;
  renderBuffer_t * _rbo;
  std::vector<rt_layout_t> _desc;
  std::vector<textureBuffer_t *> _textures;
};
