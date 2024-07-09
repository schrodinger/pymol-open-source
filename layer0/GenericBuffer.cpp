#include <iostream>
#include "GenericBuffer.h"

VertexFormatBaseType GetVertexFormatBaseType(VertexFormat format)
{
  switch (format) {
  case VertexFormat::Float:
  case VertexFormat::Float2:
  case VertexFormat::Float3:
  case VertexFormat::Float4:
    return VertexFormatBaseType::Float;
  case VertexFormat::Byte:
  case VertexFormat::Byte2:
  case VertexFormat::Byte3:
  case VertexFormat::Byte4:
    return VertexFormatBaseType::Byte;
  case VertexFormat::UByte:
  case VertexFormat::UByte2:
  case VertexFormat::UByte3:
  case VertexFormat::UByte4:
    return VertexFormatBaseType::UByte;
  case VertexFormat::Int:
  case VertexFormat::Int2:
  case VertexFormat::Int3:
  case VertexFormat::Int4:
    return VertexFormatBaseType::Int;
  case VertexFormat::UInt:
  case VertexFormat::UInt2:
  case VertexFormat::UInt3:
  case VertexFormat::UInt4:
    return VertexFormatBaseType::UInt;
  default:
    return VertexFormatBaseType::Float; // std::unreachable();
  }
}

GLenum VertexFormatToGLType(VertexFormat format)
{
  switch (format)
  {
  case VertexFormat::Byte:
  case VertexFormat::Byte2:
  case VertexFormat::Byte3:
  case VertexFormat::Byte4:
  case VertexFormat::ByteNorm:
  case VertexFormat::Byte2Norm:
  case VertexFormat::Byte3Norm:
  case VertexFormat::Byte4Norm:
    return GL_BYTE;
  case VertexFormat::UByte:
  case VertexFormat::UByte2:
  case VertexFormat::UByte3:
  case VertexFormat::UByte4:
  case VertexFormat::UByteNorm:
  case VertexFormat::UByte2Norm:
  case VertexFormat::UByte3Norm:
  case VertexFormat::UByte4Norm:
    return GL_UNSIGNED_BYTE;
  case VertexFormat::Float:
  case VertexFormat::Float2:
  case VertexFormat::Float3:
  case VertexFormat::Float4:
    return GL_FLOAT;
  case VertexFormat::Int:
  case VertexFormat::Int2:
  case VertexFormat::Int3:
  case VertexFormat::Int4:
    return GL_INT;
  case VertexFormat::UInt:
  case VertexFormat::UInt2:
  case VertexFormat::UInt3:
  case VertexFormat::UInt4:
    return GL_UNSIGNED_INT;
  default:
    return GL_INVALID_ENUM; // std::unreachable()
  }
}

GLint VertexFormatToGLSize(VertexFormat format)
{
  switch (format)
  {
  case VertexFormat::Byte:
  case VertexFormat::UByte:
  case VertexFormat::ByteNorm:
  case VertexFormat::UByteNorm:
  case VertexFormat::Float:
  case VertexFormat::Int:
  case VertexFormat::UInt:
    return 1;
  case VertexFormat::Byte2:
  case VertexFormat::UByte2:
  case VertexFormat::Byte2Norm:
  case VertexFormat::UByte2Norm:
  case VertexFormat::Float2:
  case VertexFormat::Int2:
  case VertexFormat::UInt2:
    return 2;
  case VertexFormat::Byte3:
  case VertexFormat::UByte3:
  case VertexFormat::Byte3Norm:
  case VertexFormat::UByte3Norm:
  case VertexFormat::Float3:
  case VertexFormat::Int3:
  case VertexFormat::UInt3:
    return 3;
  case VertexFormat::Byte4:
  case VertexFormat::UByte4:
  case VertexFormat::Byte4Norm:
  case VertexFormat::UByte4Norm:
  case VertexFormat::Float4:
  case VertexFormat::Int4:
  case VertexFormat::UInt4:
    return 4;
  default:
    return 0; // std::unreachable()
  }
}

bool VertexFormatIsNormalized(VertexFormat format)
{
  switch (format)
  {
  case VertexFormat::ByteNorm:
  case VertexFormat::Byte2Norm:
  case VertexFormat::Byte3Norm:
  case VertexFormat::Byte4Norm:
  case VertexFormat::UByteNorm:
  case VertexFormat::UByte2Norm:
  case VertexFormat::UByte3Norm:
  case VertexFormat::UByte4Norm:
    return true;
  default:
    return false;
  }
}

GLboolean VertexFormatToGLNormalized(VertexFormat format)
{
  auto normalized = VertexFormatIsNormalized(format);
  return normalized ? GL_TRUE : GL_FALSE;
}

std::size_t GetSizeOfVertexFormat(VertexFormat format)
{
  switch (format)
  {
  case VertexFormat::Byte:
  case VertexFormat::UByte:
  case VertexFormat::ByteNorm:
  case VertexFormat::UByteNorm:
    return 1;
  case VertexFormat::Byte2:
  case VertexFormat::UByte2:
  case VertexFormat::Byte2Norm:
  case VertexFormat::UByte2Norm:
    return 2;
  case VertexFormat::Byte3:
  case VertexFormat::UByte3:
  case VertexFormat::Byte3Norm:
  case VertexFormat::UByte3Norm:
    return 3;
  case VertexFormat::Byte4:
  case VertexFormat::UByte4:
  case VertexFormat::Byte4Norm:
  case VertexFormat::UByte4Norm:
    return 4;
  case VertexFormat::Float:
    return 4;
  case VertexFormat::Float2:
    return 8;
  case VertexFormat::Float3:
    return 12;
  case VertexFormat::Float4:
    return 16;
  case VertexFormat::Int:
    return 4;
  case VertexFormat::Int2:
    return 8;
  case VertexFormat::Int3:
    return 12;
  case VertexFormat::Int4:
    return 16;
  case VertexFormat::UInt:
    return 4;
  case VertexFormat::UInt2:
    return 8;
  case VertexFormat::UInt3:
    return 12;
  case VertexFormat::UInt4:
    return 16;
  default:
    return 0; // std::unreachable()
  }
}

GenericBuffer::GenericBuffer(buffer_layout layout, GLenum usage)
    : m_buffer_usage(usage)
    , m_layout(layout)
{
}

GenericBuffer::~GenericBuffer()
{
  for (auto i = 0; i < m_desc.size(); ++i) {
    auto& glID = desc_glIDs[i];
    if (glID) {
      glDeleteBuffers(1, &glID);
    }
  }
  if (m_interleavedID) {
    glDeleteBuffers(1, &m_interleavedID);
  }
}

bool GenericBuffer::bufferData(BufferDataDesc&& desc)
{
  m_desc = std::move(desc);
  desc_glIDs = std::vector<GLuint>(m_desc.size());
  return evaluate();
}

bool GenericBuffer::bufferData(
    BufferDataDesc&& desc, const void* data, size_t len, size_t stride)
{
  bool ok = true;
  m_desc = std::move(desc);
  desc_glIDs = std::vector<GLuint>(m_desc.size());
  m_interleaved = true;
  m_stride = stride;
  ok = genBuffer(m_interleavedID, len, data);
  return ok;
}

void GenericBuffer::bufferSubData(size_t offset, size_t size, void* data, size_t index)
{
  assert("Invalid Desc index" && index < m_desc.size());
  assert("Invalid GLDesc index" && index < desc_glIDs.size());
  auto glID = m_interleaved ? m_interleavedID : desc_glIDs[index];
  glBindBuffer(bufferType(), glID);
  glBufferSubData(bufferType(), offset, size, data);
}

void GenericBuffer::bufferReplaceData(size_t offset, size_t len, const void* data)
{
  glBindBuffer(bufferType(), m_interleavedID);
  glBufferSubData(bufferType(), offset, len, data);
}

bool GenericBuffer::evaluate()
{
  if (bufferType() == GL_ELEMENT_ARRAY_BUFFER) {
    return seqBufferData();
  } else {
    switch (m_layout) {
    case buffer_layout::SEPARATE:
      return sepBufferData();
      break;
    case buffer_layout::SEQUENTIAL:
      return seqBufferData();
      break;
    case buffer_layout::INTERLEAVED:
      return interleaveBufferData();
      break;
    }
  }
  return true; // unreacheable/Should be false?
}

bool GenericBuffer::sepBufferData()
{
  for (auto i = 0; i < m_desc.size(); ++i) {
    // If the specified size is 0 but we have a valid pointer
    // then we are going to glVertexAttribXfv X in {1,2,3,4}
    const auto& d = m_desc[i];
    auto& glID = desc_glIDs[i];
    if (d.data_ptr && (m_buffer_usage == GL_STATIC_DRAW)) {
      if (d.data_size) {
        if (!genBuffer(glID, d.data_size, d.data_ptr)) {
          return false;
        }
      }
    }
  }
  return true;
}

bool GenericBuffer::seqBufferData() {
  // this is only going to use a single opengl vbo
  m_interleaved = true;

  size_t buffer_size { 0 };
  for ( auto & d : m_desc ) {
    buffer_size += d.data_size;
  }

  std::vector<std::uint8_t> buffer_data(buffer_size);
  auto data_ptr = buffer_data.data();
  size_t offset = 0;

  for ( auto & d : m_desc ) {
    d.offset = offset;
    if (d.data_ptr)
      memcpy(data_ptr, d.data_ptr, d.data_size);
    else
      memset(data_ptr, 0, d.data_size);
    data_ptr += d.data_size;
    offset += d.data_size;
  }

  return genBuffer(m_interleavedID, buffer_size, buffer_data.data());
}

bool GenericBuffer::interleaveBufferData()
{
  const std::size_t bufferCount = m_desc.size();
  std::size_t stride = 0;
  std::vector<const uint8_t*> data_table(bufferCount);
  std::vector<const uint8_t*> ptr_table(bufferCount);
  std::vector<std::size_t> size_table(bufferCount);
  std::size_t count =
      m_desc[0].data_size / GetSizeOfVertexFormat(m_desc[0].m_format);

  // Maybe assert that all pointers in d_desc are valid?
  for (size_t i = 0; i < bufferCount; ++i) {
    auto& d = m_desc[i];
    // offset is the current stride
    d.offset = stride;

    // These must come after so that offset starts at 0
    // Size of 3 normals or whatever the current type is
    size_table[i] = GetSizeOfVertexFormat(d.m_format);

    // Increase our current estimate of the stride by this amount
    stride += size_table[i];

    // Does the addition of that previous stride leave us on a word boundry?
    int m = stride % 4;
    stride = (m ? (stride + (4 - m)) : stride);

    // data_table a pointer to the begining of each array
    data_table[i] = static_cast<const std::uint8_t*>(d.data_ptr);

    // We will move these pointers along by the values in the size table
    ptr_table[i] = data_table[i];
  }

  m_stride = stride;

  std::size_t interleavedSize = count * stride;

  std::vector<std::uint8_t> interleavedData(interleavedSize);
  auto iPtr = interleavedData.data();

  while (iPtr != (interleavedData.data() + interleavedSize)) {
    for (size_t i = 0; i < bufferCount; ++i) {
      if (ptr_table[i]) {
        memcpy(iPtr, ptr_table[i], size_table[i]);
        ptr_table[i] += size_table[i];
      }
      iPtr += size_table[i];
    }
  }

  m_interleaved = true;
  return genBuffer(m_interleavedID, interleavedSize, interleavedData.data());
}

bool GenericBuffer::genBuffer(GLuint& id, size_t size, const void* ptr)
{
  glGenBuffers(1, &id);
  if (!glCheckOkay())
    return false;
  glBindBuffer(bufferType(), id);
  if (!glCheckOkay())
    return false;
  glBufferData(bufferType(), size, ptr, GL_STATIC_DRAW);
  if (!glCheckOkay())
    return false;
  return true;
}

void VertexBuffer::bind_attrib(GLuint prg, const BufferDesc& d, GLuint glID)
{
  GLint loc = glGetAttribLocation(prg, d.attr_name);
  auto type_dim = VertexFormatToGLSize(d.m_format);
  auto type = VertexFormatToGLType(d.m_format);
  auto data_norm = VertexFormatToGLNormalized(d.m_format);
  bool masked = false;
  for (GLint lid : m_attribmask)
    if (lid == loc)
      masked = true;
  if (loc >= 0)
    m_locs.push_back(loc);
  if (loc >= 0 && !masked) {
    if (!m_interleaved && glID)
      glBindBuffer(bufferType(), glID);
    glEnableVertexAttribArray(loc);
    glVertexAttribPointer(loc, type_dim, type, data_norm, m_stride,
        reinterpret_cast<const void*>(d.offset));
  }
};

VertexBuffer::VertexBuffer(buffer_layout layout, GLenum usage)
    : GenericBuffer(layout, usage)
{
}

void VertexBuffer::bind() const
{
  // we shouldn't use this one
  if (m_interleaved)
    glBindBuffer(bufferType(), m_interleavedID);
}

void VertexBuffer::bind(GLuint prg, int index)
{
  if (index >= 0) {
    glBindBuffer(bufferType(), m_interleavedID);
    bind_attrib(prg, m_desc[index], desc_glIDs[index]);
  } else {
    if (m_interleaved && m_interleavedID)
      glBindBuffer(bufferType(), m_interleavedID);
    for (auto i = 0; i < m_desc.size(); ++i) {
      const auto& d = m_desc[i];
      auto glID = desc_glIDs[i];
      bind_attrib(prg, d, glID);
    }
    m_attribmask.clear();
  }
}

void VertexBuffer::unbind()
{
  for (auto& d : m_locs) {
    glDisableVertexAttribArray(d);
  }
  m_locs.clear();
  glBindBuffer(bufferType(), 0);
}

void VertexBuffer::maskAttributes(std::vector<GLint> attrib_locs)
{
  m_attribmask = std::move(attrib_locs);
}

void VertexBuffer::maskAttribute(GLint attrib_loc)
{
  m_attribmask.push_back(attrib_loc);
}

std::uint32_t VertexBuffer::bufferType() const
{
  return GL_ARRAY_BUFFER;
}

void IndexBuffer::bind() const
{
  glBindBuffer(bufferType(), m_interleavedID);
}

void IndexBuffer::unbind()
{
  glBindBuffer(bufferType(), 0);
}

std::uint32_t IndexBuffer::bufferType() const
{
  return GL_ELEMENT_ARRAY_BUFFER;
}

/***********************************************************************
 * RENDERBUFFER
 ***********************************************************************/
#ifdef PURE_OPENGL_ES_2
const static int rbo_lut[rbo::storage::COUNT] = { GL_DEPTH_COMPONENT16,
                                                  GL_DEPTH_COMPONENT16 };
#else
const static int rbo_lut[rbo::storage::COUNT] = { GL_DEPTH_COMPONENT16,
                                                  GL_DEPTH_COMPONENT24 };
#endif
void rbo::unbind() { glBindRenderbuffer(GL_RENDERBUFFER, 0); }

void renderBuffer_t::genBuffer() {
  glGenRenderbuffers(1, &_id);
  glBindRenderbuffer(GL_RENDERBUFFER, _id);
  glRenderbufferStorage(GL_RENDERBUFFER, rbo_lut[(int)_storage],
                        _width,
                        _height);
  glCheckOkay();
}

void renderBuffer_t::freeBuffer() { glDeleteRenderbuffers(1, &_id); }

void renderBuffer_t::bind() const { glBindRenderbuffer(GL_RENDERBUFFER, _id); }

void renderBuffer_t::unbind() const { glBindRenderbuffer(GL_RENDERBUFFER, 0); }

/***********************************************************************
 * TEXTURE
 ***********************************************************************/
const static int tex_lut[tex::max_params] = {
#ifdef PURE_OPENGL_ES_2
  GL_TEXTURE_2D,             GL_TEXTURE_2D,
  GL_TEXTURE_2D,             GL_RGBA,
  GL_RGBA,                   GL_RGB,
#else
  GL_TEXTURE_1D,             GL_TEXTURE_2D,
  GL_TEXTURE_3D,             GL_RED,
  GL_RG,                     GL_RGB,
#endif
  GL_RGBA,                   GL_UNSIGNED_BYTE,
  GL_FLOAT,                  GL_FLOAT,
  GL_NEAREST,                GL_LINEAR,
  GL_NEAREST_MIPMAP_NEAREST, GL_NEAREST_MIPMAP_LINEAR,
  GL_LINEAR_MIPMAP_NEAREST,  GL_LINEAR_MIPMAP_LINEAR,
  GL_REPEAT,
  GL_CLAMP,
  GL_REPEAT, // GL_MIRROR_REPEAT,
#ifdef PURE_OPENGL_ES_2
  GL_CLAMP_TO_EDGE,          GL_CLAMP_TO_EDGE,
#else
  GL_CLAMP_TO_EDGE,          GL_CLAMP_TO_BORDER,
#endif
  GL_CLAMP_TO_EDGE, // GL_MIRROR_CLAMP_TO_EDGE
#ifdef PURE_OPENGL_ES_2
  GL_REPLACE,       GL_REPLACE
#else
  GL_TEXTURE_ENV_MODE,       GL_REPLACE
#endif
};

#ifndef GL_R8
#define GL_R8    GL_RGBA
#define GL_RG8   GL_RGBA
#define GL_RGB8  GL_RGB
#define GL_RGBA8 GL_RGBA
#endif

#ifndef GL_R32F
#define GL_R32F    GL_R8
#define GL_RG32F   GL_RG8
#define GL_RGB32F  GL_RGB8
#define GL_RGBA32F GL_RGBA8
#endif

#ifndef GL_R16F
#define GL_R16F    GL_R32F
#define GL_RG16F   GL_RG32F
#define GL_RGB16F  GL_RGB32F
#define GL_RGBA16F GL_RGBA32F
#endif

#ifdef _WEBGL
static int tex_format_internal_byte(tex::format f) {
  return tex_lut[(int)f];
};

static int tex_format_internal_float(tex::format f) {
  return tex_format_internal_byte(f);
};

static int tex_format_internal_half_float(tex::format f) {
  return tex_format_internal_byte(f);
};
#else
static int tex_format_internal_float(tex::format f) {
  using namespace tex;
  switch (f) {
  case format::R:
    return GL_R32F;
    break;
  case format::RG:
    return GL_RG32F;
    break;
  case format::RGB:
    return GL_RGB32F;
    break;
  case format::RGBA:
    return GL_RGBA32F;
    break;
  default:
    return GL_RGBA32F;
    break;
  }
};

static int tex_format_internal_half_float(tex::format f) {
  using namespace tex;
  switch (f) {
  case format::R:
    return GL_R16F;
    break;
  case format::RG:
    return GL_RG16F;
    break;
  case format::RGB:
    return GL_RGB16F;
    break;
  case format::RGBA:
    return GL_RGBA16F;
    break;
  default:
    return GL_RGBA16F;
    break;
  }
};

static int tex_format_internal_byte(tex::format f) {
  using namespace tex;
  switch (f) {
  case format::R:
    return GL_R8;
    break;
  case format::RG:
    return GL_RG8;
    break;
  case format::RGB:
    return GL_RGB8;
    break;
  case format::RGBA:
    return GL_RGBA8;
    break;
  default:
    return GL_RGBA8;
    break;
  }
};
#endif

template <typename T> static int tex_tab(T val) { return tex_lut[(int)val]; }

void tex::env(tex::env_name name, tex::env_param param) {
#ifdef PURE_OPENGL_ES_2
  fprintf(stderr,
          "No Support for Texture Env\n");
#else
  glTexEnvf(GL_TEXTURE_ENV, tex_tab(name), tex_tab(param));
#endif
}

void textureBuffer_t::genBuffer() {
  GLenum dim = tex_tab(_dim);
  glGenTextures(1, &_id);
  glBindTexture(dim, _id);
  glTexParameteri(dim, GL_TEXTURE_MAG_FILTER, tex_tab(_sampling[0]));
  glTexParameteri(dim, GL_TEXTURE_MIN_FILTER, tex_tab(_sampling[1]));
  glTexParameteri(dim, GL_TEXTURE_WRAP_S, tex_tab(_sampling[2]));
  if (_sampling[3])
    glTexParameteri(dim, GL_TEXTURE_WRAP_T, tex_tab(_sampling[3]));
#ifndef PURE_OPENGL_ES_2
  if (_sampling[4])
    glTexParameteri(dim, GL_TEXTURE_WRAP_R, tex_tab(_sampling[4]));
#endif

  glCheckOkay();
}

void textureBuffer_t::freeBuffer() { glDeleteTextures(1, &_id); }

void textureBuffer_t::texture_data_1D(int width, const void *data) {
#ifdef PURE_OPENGL_ES_2
  fprintf(stderr,
          "No support for 1D textures\n");
#else
  using namespace tex;
  _width = width;
  bind();
  switch (_type) {
  case data_type::HALF_FLOAT:
    glTexImage1D(GL_TEXTURE_1D, 0, tex_format_internal_half_float(_format), _width, 0,
                 tex_tab(_format), tex_tab(tex::data_type::FLOAT), data);
    break;
  case data_type::FLOAT:
    glTexImage1D(GL_TEXTURE_1D, 0, tex_format_internal_float(_format), _width, 0,
                 tex_tab(_format), tex_tab(_type), data);
    break;
  case data_type::UBYTE:
    glTexImage1D(GL_TEXTURE_1D, 0, tex_format_internal_byte(_format), _width, 0,
                 tex_tab(_format), tex_tab(_type), data);
    break;
  default:
    break;
  };
  glCheckOkay();
#endif
}

void textureBuffer_t::texture_data_2D(int width, int height, const void *data) {
  using namespace tex;
  _width = width;
  _height = height;
  bind();
  switch (_type) {
  case data_type::HALF_FLOAT:
    glTexImage2D(GL_TEXTURE_2D, 0, tex_format_internal_half_float(_format), _width,
                 _height, 0, tex_tab(_format), tex_tab(tex::data_type::FLOAT), data);
    break;
  case data_type::FLOAT:
    glTexImage2D(GL_TEXTURE_2D, 0, tex_format_internal_float(_format), _width,
                 _height, 0, tex_tab(_format), tex_tab(_type), data);
    break;
  case data_type::UBYTE:
    glTexImage2D(GL_TEXTURE_2D, 0, tex_format_internal_byte(_format), _width,
                 _height, 0, tex_tab(_format), tex_tab(_type), data);
    break;
  default:
    break;
  }
  glCheckOkay();
}

void textureBuffer_t::texture_subdata_2D(
    int xoffset, int yoffset, int width, int height, const void* data)
{
  using namespace tex;
  bind();
  switch (_type) {
  case data_type::HALF_FLOAT:
    glTexSubImage2D(GL_TEXTURE_2D, 0, xoffset, yoffset, width, height,
                    tex_tab(_format), tex_tab(tex::data_type::FLOAT), data);
    break;
  case data_type::FLOAT:
  case data_type::UBYTE:
    glTexSubImage2D(GL_TEXTURE_2D, 0, xoffset, yoffset, width, height,
                    tex_tab(_format), tex_tab(_type), data);
    break;
  }
  glCheckOkay();
}

void textureBuffer_t::texture_data_3D(int width, int height, int depth,
                                      const void *data) {
#ifdef PURE_OPENGL_ES_2
  fprintf(stderr,
          "No support for 3D textures\n");
#else
  _width = width;
  _height = height;
  _depth = depth;
  bind();
  switch(_type) {
  case tex::data_type::HALF_FLOAT:
    glTexImage3D(GL_TEXTURE_3D, 0, tex_format_internal_half_float(_format), _width,
                 _height, _depth, 0, tex_tab(_format), tex_tab(tex::data_type::FLOAT), data);
    break;
  case tex::data_type::FLOAT:
    glTexImage3D(GL_TEXTURE_3D, 0, tex_format_internal_float(_format), _width,
                 _height, _depth, 0, tex_tab(_format), tex_tab(tex::data_type::FLOAT), data);
    break;
  case tex::data_type::UBYTE:
    glTexImage3D(GL_TEXTURE_3D, 0, tex_format_internal_byte(_format), _width,
                 _height, _depth, 0, tex_tab(_format), tex_tab(_type), data);
    break;
  default:
    break;
  }
#endif

  glCheckOkay();
}

void textureBuffer_t::bind() const { glBindTexture(tex_tab(_dim), _id); }

void textureBuffer_t::bindToTextureUnit(std::uint8_t textureUnit) const
{
  glActiveTexture(GL_TEXTURE0 + textureUnit);
  bind();
}

void textureBuffer_t::unbind() const { glBindTexture(tex_tab(_dim), 0); }
/***********************************************************************
 * FRAMEBUFFER
 ***********************************************************************/
const static int fbo_lut[fbo::attachment::COUNT] = {
  GL_COLOR_ATTACHMENT0,
#if defined(PURE_OPENGL_ES_2) && !defined(_WEBGL)
  GL_COLOR_ATTACHMENT0,
  GL_COLOR_ATTACHMENT0,
  GL_COLOR_ATTACHMENT0,
#else
  GL_COLOR_ATTACHMENT1,
  GL_COLOR_ATTACHMENT2,
  GL_COLOR_ATTACHMENT3,
#endif
  GL_DEPTH_ATTACHMENT
};

template <typename T> static int fbo_tab(T val) { return fbo_lut[(int)val]; }

void fbo::unbind() { glBindFramebuffer(GL_FRAMEBUFFER, 0); }

void frameBuffer_t::genBuffer() { glGenFramebuffers(1, &_id); }

void frameBuffer_t::freeBuffer() { glDeleteFramebuffers(1, &_id); }

void frameBuffer_t::attach_texture(textureBuffer_t *texture,
                                   fbo::attachment loc) {
  size_t id = texture->get_hash_id();
  _attachments.emplace_back(id, loc);
  bind();
  glFramebufferTexture2D(GL_FRAMEBUFFER, fbo_tab(loc), GL_TEXTURE_2D,//tex_tab(texture->_dim),
                         texture->_id, 0);
  checkStatus();
}

void frameBuffer_t::attach_renderbuffer(renderBuffer_t *renderbuffer,
                                        fbo::attachment loc) {
  size_t id = renderbuffer->get_hash_id();
  _attachments.emplace_back(id, loc);
  bind();
  glFramebufferRenderbuffer(GL_FRAMEBUFFER, fbo_tab(loc), GL_RENDERBUFFER,
                            renderbuffer->_id);
  checkStatus();
}

void frameBuffer_t::bind() const {
  glBindFramebuffer(GL_FRAMEBUFFER, _id);
}

void frameBuffer_t::checkStatus() {
  GLenum status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
  switch (status) {
  case GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT:
    printf("Incomplete attachment\n");
    break;
#ifndef PURE_OPENGL_ES_2
  case GL_FRAMEBUFFER_INCOMPLETE_DIMENSIONS_EXT:
    printf("Incomplete dimensions\n");
    break;
#endif
  case GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT:
    printf("Incomplete missing attachment\n");
    break;
  case GL_FRAMEBUFFER_UNSUPPORTED:
    printf("Framebuffer combination unsupported\n");
    break;
  }
}

/***********************************************************************
 * RENDERTARGET
 ***********************************************************************/
renderTarget_t::~renderTarget_t() {
  for (auto &t : _textures)
    delete t;

  if (_fbo)
    delete _fbo;

  if (_rbo && !_shared_rbo)
    delete _rbo;
}

void renderTarget_t::layout(std::vector<rt_layout_t> &&desc,
                            renderBuffer_t *with_rbo) {
  _fbo = new frameBuffer_t();
  if (with_rbo) {
    _rbo = with_rbo;
    _shared_rbo = true;
  } else {
    _rbo = new renderBuffer_t(_size.x, _size.y, rbo::storage::DEPTH24);
  }
  for (auto &d : desc) {
    if (!d.width)
      d.width = _size.x;
    if (!d.height)
      d.height = _size.y;

    tex::data_type type;
    switch (d.type) {
    case rt_layout_t::UBYTE:
      type = tex::data_type::UBYTE;
      break;
    case rt_layout_t::FLOAT:
      type = tex::data_type::FLOAT;
      break;
    default:
      printf("Error: %s:%d\n", __FILE__, __LINE__);
      return;
    }

    tex::format format;
    switch (d.nchannels) {
    case 1:
      format = tex::format::R;
      break;
    case 2:
      format = tex::format::RG;
      break;
    case 3:
      format = tex::format::RGB;
      break;
    case 4:
      format = tex::format::RGBA;
      break;
    default:
      printf("Error: %s:%d\n", __FILE__, __LINE__);
      return;
    }

    _textures.push_back(new textureBuffer_t(
        format, type, tex::filter::LINEAR, tex::filter::LINEAR,
        tex::wrap::CLAMP, tex::wrap::CLAMP));
    auto tex = _textures.back();
    tex->texture_data_2D(d.width, d.height, nullptr);

    fbo::attachment loc;
    switch (_textures.size()) {
    case 1:
      loc = fbo::attachment::COLOR0;
      break;
    case 2:
      loc = fbo::attachment::COLOR1;
      break;
    case 3:
      loc = fbo::attachment::COLOR2;
      break;
    case 4:
      loc = fbo::attachment::COLOR3;
      break;
    default:
      loc = fbo::attachment::COLOR0;
      // print error
      break;
    }

    _fbo->attach_texture(tex, loc);
  }
  _fbo->attach_renderbuffer(_rbo, fbo::attachment::DEPTH);
  _desc = std::move(desc);
  glCheckOkay();
}

void renderTarget_t::resize(shape_type size) {
  _size = size;
  if (!_shared_rbo) {
    delete _rbo;
    _rbo = nullptr;
  }

  for (auto i : _textures) {
    delete i;
  }
  _textures.clear();

  delete _fbo;

  std::vector<rt_layout_t> desc;
  for (auto& i : _desc) {
    desc.emplace_back(i.nchannels, i.type, size.x, size.y);
  }

  layout(std::move(desc), _rbo);
}

void renderTarget_t::bind(bool clear) const {
  _fbo->bind();
  if (clear) {
    glClearColor(0.f, 0.f, 0.f, 0.f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  }
}
