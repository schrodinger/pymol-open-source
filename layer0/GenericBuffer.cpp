#include <iostream>
#include "GenericBuffer.h"
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
  switch ((int)_type) {
  case (int)data_type::HALF_FLOAT:
    glTexImage1D(GL_TEXTURE_1D, 0, tex_format_internal_half_float(_format), _width, 0,
                 tex_tab(_format), tex_tab(tex::data_type::FLOAT), data);
    break;
  case (int)data_type::FLOAT:
    glTexImage1D(GL_TEXTURE_1D, 0, tex_format_internal_float(_format), _width, 0,
                 tex_tab(_format), tex_tab(_type), data);
    break;
  case (int)data_type::UBYTE:
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
  switch ((int)_type) {
  case (int)data_type::HALF_FLOAT:
    glTexImage2D(GL_TEXTURE_2D, 0, tex_format_internal_half_float(_format), _width,
                 _height, 0, tex_tab(_format), tex_tab(tex::data_type::FLOAT), data);
    break;
  case (int)data_type::FLOAT:
    glTexImage2D(GL_TEXTURE_2D, 0, tex_format_internal_float(_format), _width,
                 _height, 0, tex_tab(_format), tex_tab(_type), data);
    break;
  case (int)data_type::UBYTE:
    glTexImage2D(GL_TEXTURE_2D, 0, tex_format_internal_byte(_format), _width,
                 _height, 0, tex_tab(_format), tex_tab(_type), data);
    break;
  default:
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
  switch((int)_type) {
  case (int)tex::data_type::HALF_FLOAT:
    glTexImage3D(GL_TEXTURE_3D, 0, tex_format_internal_half_float(_format), _width,
                 _height, _depth, 0, tex_tab(_format), tex_tab(tex::data_type::FLOAT), data);
  case (int)tex::data_type::FLOAT:
    glTexImage3D(GL_TEXTURE_3D, 0, tex_format_internal_float(_format), _width,
                 _height, _depth, 0, tex_tab(_format), tex_tab(_type), data);
    break;
  case (int)tex::data_type::UBYTE:
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
