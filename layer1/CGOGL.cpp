#include "CGOGL.h"

#include "CGO.h"
#include "CGOGL.h"
#include "CGORenderer.h"
#include "CoordSet.h"
#include "Feedback.h"
#include "Scene.h"
#include "SceneDef.h"
#include "ShaderMgr.h"
#include "Sphere.h"
#include "os_gl.h"
#include "os_gl_cgo.h"

#define VAR_FOR_NORMAL pl
#define VERTEX_NORMAL_SIZE 3
#define VAR_FOR_NORMAL_CNT_PLUS

#ifdef PURE_OPENGL_ES_2
#define glVertexAttrib4ubv(loc, data) glVertexAttrib4f(loc, \
    (data)[0] / 255.f, (data)[1] / 255.f, (data)[2] / 255.f, (data)[3] / 255.f);
#endif

constexpr unsigned VERTEX_PICKCOLOR_RGBA_SIZE = 1;  // 4 unsigned bytes
constexpr unsigned VERTEX_PICKCOLOR_INDEX_SIZE = 2; // index + bond
constexpr unsigned VERTEX_PICKCOLOR_SIZE = VERTEX_PICKCOLOR_RGBA_SIZE + //
                                           VERTEX_PICKCOLOR_INDEX_SIZE;
constexpr unsigned VERTEX_ACCESSIBILITY_SIZE = 1;

/* ======== GL Rendering ======== */

static int CGO_gl_begin_WARNING_CALLED = false,
           CGO_gl_end_WARNING_CALLED = false,
           CGO_gl_vertex_WARNING_CALLED = false;

const float g_ones4f[4] = {1.f, 1.f, 1.f, 1.f};

static int CGOConvertDebugMode(int debug, int modeArg)
{
  int mode = modeArg;
  if (debug == 1) {
    switch (mode) {
    case GL_TRIANGLES:
      mode = GL_LINES;
      break;
    case GL_TRIANGLE_STRIP:
      mode = GL_LINE_STRIP;
      break;
    case GL_TRIANGLE_FAN:
      mode = GL_LINES;
      break;
    }
  } else {
    mode = GL_POINTS;
  }
  return mode;
}

#define CHECK_GL_ERROR_OK(printstr)                                            \
  if ((err = glGetError()) != 0) {                                             \
    PRINTFB(G, FB_CGO, FB_Errors) printstr, err ENDFB(G);                      \
  }

void CheckGLErrorOK(PyMOLGlobals* G, pymol::zstring_view errString)
{
  GLenum err;
  CHECK_GL_ERROR_OK(errString.c_str());
}

static void CGO_gl_begin(CCGORenderer* I, CGO_op_data pc)
{
#ifndef PURE_OPENGL_ES_2
  if (I->use_shader) {
#endif
    if (!CGO_gl_begin_WARNING_CALLED) {
      PRINTFB(I->G, FB_CGO, FB_Warnings)
      " CGO_gl_begin() is called but not implemented in OpenGLES\n" ENDFB(I->G);
      CGO_gl_begin_WARNING_CALLED = true;
    }
#ifndef PURE_OPENGL_ES_2
  } else {
    int mode = CGO_get_int(*pc);
    if (I->debug)
      mode = CGOConvertDebugMode(I->debug, mode);
    glBegin(mode);
  }
#endif
}
static void CGO_gl_end(CCGORenderer* I, CGO_op_data)
{
#ifndef PURE_OPENGL_ES_2
  if (I->use_shader) {
#endif
    if (!CGO_gl_end_WARNING_CALLED) {
      PRINTFB(I->G, FB_CGO, FB_Warnings)
      " CGO_gl_end() is called but not implemented in OpenGLES\n" ENDFB(I->G);
      CGO_gl_end_WARNING_CALLED = true;
    }
#ifndef PURE_OPENGL_ES_2
  } else {
    glEnd();
  }
#endif
}
static void CGO_gl_vertex(CCGORenderer* I, CGO_op_data v)
{
#ifndef PURE_OPENGL_ES_2
  if (I->use_shader) {
#endif
    if (!CGO_gl_vertex_WARNING_CALLED) {
      PRINTFB(I->G, FB_CGO, FB_Warnings)
      " CGO_gl_vertex() is called but not implemented in OpenGLES\n" ENDFB(
          I->G);
      CGO_gl_vertex_WARNING_CALLED = true;
    }
#ifndef PURE_OPENGL_ES_2
  } else {
    glVertex3fv(*v);
  }
#endif
}

static void CGO_gl_vertex_cross(CCGORenderer* I, CGO_op_data v)
{
#ifndef PURE_OPENGL_ES_2
  if (I->use_shader) {
#endif
    if (!CGO_gl_vertex_WARNING_CALLED) {
      PRINTFB(I->G, FB_CGO, FB_Warnings)
      " CGO_gl_vertex() is called but not implemented in OpenGLES\n" ENDFB(
          I->G);
      CGO_gl_vertex_WARNING_CALLED = true;
    }
#ifndef PURE_OPENGL_ES_2
  } else {
    CSetting *set1 = nullptr, *set2 = nullptr;
    if (I->rep && I->rep->cs)
      set1 = I->rep->cs->Setting.get();
    if (I->rep && I->rep->obj)
      set2 = I->rep->obj->Setting.get();
    float nonbonded_size =
        SettingGet_f(I->G, set1, set2, cSetting_nonbonded_size);
    float pt[3];
    copy3f(*v, pt);
    pt[0] -= nonbonded_size;
    glVertex3fv(pt);
    pt[0] += 2 * nonbonded_size;
    glVertex3fv(pt);
    copy3f(*v, pt);
    pt[1] -= nonbonded_size;
    glVertex3fv(pt);
    pt[1] += 2 * nonbonded_size;
    glVertex3fv(pt);
    copy3f(*v, pt);
    pt[2] -= nonbonded_size;
    glVertex3fv(pt);
    pt[2] += 2 * nonbonded_size;
    glVertex3fv(pt);
  }
#endif
}

static void CGO_gl_line(CCGORenderer* I, CGO_op_data v)
{
#ifndef PURE_OPENGL_ES_2
  if (!I->use_shader) {
    auto line = reinterpret_cast<const cgo::draw::line*>(*v);
    glVertex3fv(line->vertex1);
    glVertex3fv(line->vertex2);
  }
#endif
}

static void CGO_gl_splitline(CCGORenderer* I, CGO_op_data v)
{
#ifndef PURE_OPENGL_ES_2
  if (!I->use_shader) {
    auto splitline = reinterpret_cast<const cgo::draw::splitline*>(*v);
    bool interpolation = splitline->flags & cgo::draw::splitline::interpolation;
    bool equal_colors = splitline->flags & cgo::draw::splitline::equal_colors;
    bool no_split_for_pick =
        splitline->flags & cgo::draw::splitline::no_split_for_pick;

    if (I->isPicking) {
      if (no_split_for_pick) {
        glVertex3fv(splitline->vertex1);
        glVertex3fv(splitline->vertex2);
      } else {
        float h[3];
        average3f(splitline->vertex1, splitline->vertex2, h);
        glVertex3fv(splitline->vertex1);
        glVertex3fv(h);
        unsigned char col[4];
        AssignNewPickColor(nullptr, I->info->pick, col, &I->rep->context,
            splitline->index, splitline->bond);
        glColor4ubv(col);
        glVertex3fv(h);
        glVertex3fv(splitline->vertex2);
      }
    } else if (interpolation || equal_colors) {
      glVertex3fv(splitline->vertex1);
      if (!equal_colors)
        glColor4ub(splitline->color2[0], splitline->color2[1],
            splitline->color2[2], CLIP_COLOR_VALUE(I->alpha));
      glVertex3fv(splitline->vertex2);
    } else {
      float h[3];
      average3f(splitline->vertex1, splitline->vertex2, h);
      glVertex3fv(splitline->vertex1);
      glVertex3fv(h);
      glColor4ub(splitline->color2[0], splitline->color2[1],
          splitline->color2[2], CLIP_COLOR_VALUE(I->alpha));
      glVertex3fv(h);
      glVertex3fv(splitline->vertex2);
    }
  }
#endif
}

static void CGO_gl_normal(CCGORenderer* I, CGO_op_data varg)
{
  const float* v = *varg;
#ifdef PURE_OPENGL_ES_2
  if (I->use_shader) {
    glVertexAttrib3fv(VERTEX_NORMAL, v);
  } else {
#else
  {
#endif
#ifndef PURE_OPENGL_ES_2
    glNormal3f(v[0], v[1], v[2]);
#endif
  }
}

static void CGO_gl_draw_arrays(CCGORenderer* I, CGO_op_data pc)
{
  auto sp = reinterpret_cast<const cgo::draw::arrays*>(*pc);
  int mode = sp->mode, arrays = sp->arraybits, narrays = sp->narrays,
      nverts = sp->nverts;
  const float* data = sp->floatdata;
  (void) narrays;
#ifndef PURE_OPENGL_ES_2
  if (I->use_shader) {
#endif

#ifdef _WEBGL
    uint buffers[3] = {0, 0, 0};
    glGenBuffers(3, buffers);
#endif

    if (arrays & CGO_VERTEX_ARRAY)
      glEnableVertexAttribArray(VERTEX_POS);
    if (arrays & CGO_NORMAL_ARRAY)
      glEnableVertexAttribArray(VERTEX_NORMAL);
    if (I->isPicking) {
      if (arrays & CGO_PICK_COLOR_ARRAY) {
        glEnableVertexAttribArray(VERTEX_COLOR);
      }
    } else {
      if (arrays & CGO_COLOR_ARRAY)
        glEnableVertexAttribArray(VERTEX_COLOR);
    }

    if (arrays & CGO_VERTEX_ARRAY) {
#ifdef _WEBGL
      glBindBuffer(GL_ARRAY_BUFFER, buffers[0]);
      glBufferData(
          GL_ARRAY_BUFFER, sizeof(float) * nverts * 3, data, GL_STATIC_DRAW);
      glVertexAttribPointer(
          VERTEX_POS, VERTEX_POS_SIZE, GL_FLOAT, GL_FALSE, 0, 0);
#else
    glVertexAttribPointer(
        VERTEX_POS, VERTEX_POS_SIZE, GL_FLOAT, GL_FALSE, 0, data);
#endif
      data += nverts * VERTEX_POS_SIZE;
    }
    if (arrays & CGO_NORMAL_ARRAY) {
#ifdef _WEBGL
      glBindBuffer(GL_ARRAY_BUFFER, buffers[1]);
      glBufferData(
          GL_ARRAY_BUFFER, sizeof(float) * nverts * 3, data, GL_STATIC_DRAW);
      glVertexAttribPointer(
          VERTEX_NORMAL, VERTEX_NORMAL_SIZE, GL_FLOAT, GL_FALSE, 0, 0);
#else
    glVertexAttribPointer(
        VERTEX_NORMAL, VERTEX_NORMAL_SIZE, GL_FLOAT, GL_FALSE, 0, data);
#endif
      data += nverts * VERTEX_NORMAL_SIZE;
    }
    if (I->isPicking) {
      if (arrays & CGO_COLOR_ARRAY) {
        data += nverts * VERTEX_COLOR_SIZE;
      }
      if (arrays & CGO_PICK_COLOR_ARRAY) {
#ifdef _WEBGL
        glBindBuffer(GL_ARRAY_BUFFER, buffers[2]);
        glBufferData(GL_ARRAY_BUFFER, nverts * 4, data, GL_STATIC_DRAW);
        glVertexAttribPointer(
            VERTEX_COLOR, VERTEX_COLOR_SIZE, GL_UNSIGNED_BYTE, GL_FALSE, 0, 0);
#else
      glVertexAttribPointer(
          VERTEX_COLOR, VERTEX_COLOR_SIZE, GL_UNSIGNED_BYTE, GL_FALSE, 0, data);
#endif
        data += nverts * VERTEX_PICKCOLOR_SIZE;
      }
    } else {
      if (arrays & CGO_COLOR_ARRAY) {
#ifdef _WEBGL
        glBindBuffer(GL_ARRAY_BUFFER, buffers[2]);
        glBufferData(
            GL_ARRAY_BUFFER, sizeof(float) * nverts * 4, data, GL_STATIC_DRAW);
        glVertexAttribPointer(
            VERTEX_COLOR, VERTEX_COLOR_SIZE, GL_FLOAT, GL_FALSE, 0, 0);
#else
      glVertexAttribPointer(
          VERTEX_COLOR, VERTEX_COLOR_SIZE, GL_FLOAT, GL_FALSE, 0, data);
#endif
        data += nverts * VERTEX_COLOR_SIZE;
      }
      if (arrays & CGO_PICK_COLOR_ARRAY) {
        data += nverts * VERTEX_PICKCOLOR_SIZE;
      }
    }
    if (I->debug) {
      mode = CGOConvertDebugMode(I->debug, mode);
    }
    glDrawArrays(mode, 0, nverts);

    if (I->isPicking) {
      if (arrays & CGO_PICK_COLOR_ARRAY) {
        glDisableVertexAttribArray(VERTEX_COLOR);
      }
    } else {
      if (arrays & CGO_COLOR_ARRAY)
        glDisableVertexAttribArray(VERTEX_COLOR);
    }
    if (arrays & CGO_VERTEX_ARRAY)
      glDisableVertexAttribArray(VERTEX_POS);
    if (arrays & CGO_NORMAL_ARRAY)
      glDisableVertexAttribArray(VERTEX_NORMAL);
#ifdef _WEBGL
    glDeleteBuffers(3, buffers);
#endif

#ifndef PURE_OPENGL_ES_2
  } else {
    int pl, pla, plc;
    const float* vertexVals = nullptr;
    const float *colorVals = 0, *normalVals = 0, *tmp_ptr;
    const uchar *pickColorVals = 0, *tmp_pc_ptr;
    float alpha = I->alpha;
    if (arrays & CGO_VERTEX_ARRAY) {
      vertexVals = data;
      data += nverts * VERTEX_POS_SIZE;
    }
    if (arrays & CGO_NORMAL_ARRAY) {
      normalVals = data;
      data += nverts * VERTEX_NORMAL_SIZE;
    }
    if (I->isPicking) {
      alpha = 1.f;
      if (arrays & CGO_COLOR_ARRAY) {
        data += nverts * VERTEX_COLOR_SIZE;
      }
      if (arrays & CGO_PICK_COLOR_ARRAY) {
        pickColorVals = (uchar*) data;
        data += nverts * VERTEX_PICKCOLOR_SIZE;
      }
    } else {
      if (arrays & CGO_COLOR_ARRAY) {
        colorVals = data;
        data += nverts * 4;
      }
      if (arrays & CGO_PICK_COLOR_ARRAY) {
        data += nverts * VERTEX_PICKCOLOR_SIZE;
      }
    }
    if (arrays & CGO_ACCESSIBILITY_ARRAY) {
      data += nverts * VERTEX_ACCESSIBILITY_SIZE;
    }

    if (I->debug) {
      mode = CGOConvertDebugMode(I->debug, mode);
    }

    glBegin(mode);
    for (pl = 0, pla = 0, plc = 0; pl < nverts; pl++, pla += 3, plc += 4) {
      if (pickColorVals) {
        tmp_pc_ptr =
            &pickColorVals[plc]; /* the pick colors are saved with rgba */
        glColor4ub(tmp_pc_ptr[0], tmp_pc_ptr[1], tmp_pc_ptr[2], tmp_pc_ptr[3]);
      } else {
        if (colorVals) {
          tmp_ptr = &colorVals[plc];
          glColor4f(tmp_ptr[0], tmp_ptr[1], tmp_ptr[2], alpha);
        }
        if (normalVals) {
          tmp_ptr = &normalVals[pla];
          glNormal3fv(&normalVals[pla]);
        }
      }
      if (vertexVals) {
        tmp_ptr = &vertexVals[pla];
        glVertex3fv(&vertexVals[pla]);
      }
    }
    glEnd();
  }
#endif
}

/* TransparentInfoSortIX - This function sorts all n_tri triangle
 * centroids in the array sum by:
 * 1) computing z-value in array z_value
 * 2) bin sorting z_values and placing indices in ix array (using Util.cpp)
 *
 * - uses sort_mem as pre-allocated memory to sort
 * - t_mode - either forward (1) or backwards (0) sort
 */
void TransparentInfoSortIX(PyMOLGlobals* G, float* sum, float* z_value, int* ix,
    int n_tri, int* sort_mem, int t_mode)
{
  float* zv;
  float* sv;
  float matrix[16];
  int idx;

#ifdef PURE_OPENGL_ES_2
  copy44f(SceneGetModelViewMatrixPtr(G), matrix);
#else
  glGetFloatv(GL_MODELVIEW_MATRIX, matrix);
#endif
  zv = z_value;
  sv = sum;

  /* for each triangle, computes the z */
  for (idx = 0; idx < n_tri; ++idx) {
    *(zv++) = matrix[2] * sv[0] + matrix[6] * sv[1] + matrix[10] * sv[2];
    sv += 3;
  }

  UtilZeroMem(sort_mem, sizeof(int) * (n_tri + 256));

  switch (t_mode) {
  case 1:
    UtilSemiSortFloatIndexWithNBinsImpl(
        sort_mem, n_tri, 256, z_value, ix, true); // front to back
    /* UtilSortIndex(n_tri,z_value,ix,(UtilOrderFn*)ZOrderFn); */
    break;
  default:
    UtilSemiSortFloatIndexWithNBinsImpl(
        sort_mem, n_tri, 256, z_value, ix, false); // back to front
    /* UtilSortIndex(n_tri,z_value,ix,(UtilOrderFn*)ZRevOrderFn); */
    break;
  }
}

/**
 * CGOReorderIndicesWithTransparentInfo : This function
 * takes the triangle index array ix (result from TransparentInfoSortIX)
 * and sets the vertices (vertexIndices) for each triangle from the original
 * indices (vertexIndicesOriginal), then uses glBufferData to set the
 * GL_ELEMENT_ARRAY_BUFFER to these indices.
 *
 */
static void CGOReorderIndicesWithTransparentInfo(PyMOLGlobals* G, int nindices,
    size_t vbuf, int n_tri, int* ix, VertexIndex_t* vertexIndicesOriginal,
    VertexIndex_t* vertexIndices)
{
  int c, pl, idx;
  IndexBuffer* ibo = G->ShaderMgr->getGPUBuffer<IndexBuffer>(vbuf);
  if (!vertexIndices) {
    PRINTFB(G, FB_RepSurface, FB_Errors)
    "ERROR: RepSurfaceRender() vertexIndices is not set, nindices=%d\n",
        nindices ENDFB(G);
  }
  /* updates the vertexIndices from the ix array */
  for (c = 0, pl = 0; c < n_tri; c++) {
    idx = ix[c] * 3;
    vertexIndices[pl++] = vertexIndicesOriginal[idx];
    vertexIndices[pl++] = vertexIndicesOriginal[idx + 1];
    vertexIndices[pl++] = vertexIndicesOriginal[idx + 2];
  }
  ibo->bufferSubData(0, sizeof(VertexIndex_t) * nindices, vertexIndices);
}

static void CGO_gl_draw_buffers_indexed(CCGORenderer* I, CGO_op_data pc)
{
  auto sp = reinterpret_cast<const cgo::draw::buffers_indexed*>(*pc);
  int mode = sp->mode, nindices = sp->nindices, nverts = sp->nverts,
      n_data = sp->n_data;
  size_t vboid = sp->vboid, iboid = sp->iboid;
  VertexBuffer* vbo = I->G->ShaderMgr->getGPUBuffer<VertexBuffer>(vboid);
  IndexBuffer* ibo = I->G->ShaderMgr->getGPUBuffer<IndexBuffer>(iboid);
  CheckGLErrorOK(I->G, "beginning of CGO_gl_draw_buffers_indexed err=%d\n");

  auto shaderPrg = I->G->ShaderMgr->Get_Current_Shader();

  if (!shaderPrg) {
    return;
  }

  if (I->isPicking) {
    int attr_a_Color = shaderPrg->GetAttribLocation("a_Color");
    vbo->maskAttributes({attr_a_Color});
    shaderPrg->Set1i("fog_enabled", 0);
    shaderPrg->Set1i("lighting_enabled", 0);
    if (I->use_shader) {
      if (sp->pickvboid) {
        VertexBuffer* pickvbo =
            I->G->ShaderMgr->getGPUBuffer<VertexBuffer>(sp->pickvboid);
        pickvbo->bind(shaderPrg->id, I->pick_pass());
      } else {
        glEnableVertexAttribArray(attr_a_Color);
        glVertexAttribPointer(attr_a_Color, VERTEX_COLOR_SIZE, GL_UNSIGNED_BYTE,
            GL_TRUE, 0, sp->floatdata);
      }
    }
  }
  if (n_data) {
    // if transparency data, then sort it
    int n_tri = nindices / 3;
    float* sum = sp->floatdata + nverts * 3;
    float* z_value = sum + (nindices * 3);
    int* ix = (int*) (z_value + n_tri);
    int* sort_mem = ix + n_tri;
    int t_mode;
    CSetting *set1 = nullptr, *set2 = nullptr;
    if (I->rep && I->rep->cs)
      set1 = I->rep->cs->Setting.get();
    if (I->rep && I->rep->obj)
      set2 = I->rep->obj->Setting.get();
    t_mode = SettingGet_i(I->G, set1, set2, cSetting_transparency_mode);
    if (t_mode != 3) {
      auto vertexIndicesOriginalTI = (VertexIndex_t*) (sort_mem + n_tri + 256);
      auto vertexIndicesTI = vertexIndicesOriginalTI + nindices;
      TransparentInfoSortIX(I->G, sum, z_value, ix, n_tri, sort_mem, t_mode);
      CGOReorderIndicesWithTransparentInfo(I->G, nindices, iboid, n_tri, ix,
          vertexIndicesOriginalTI, vertexIndicesTI);
    }
  }

  if (I->debug) {
    mode = CGOConvertDebugMode(I->debug, mode);
  }
  vbo->bind(shaderPrg->id);
  ibo->bind();

  CheckGLErrorOK(
      I->G, "CGO_gl_draw_buffers_indexed: before glDrawElements err=%d\n");
  glDrawElements(mode, nindices, VertexIndex_GL_ENUM, 0);
  CheckGLErrorOK(
      I->G, "CGO_gl_draw_buffers_indexed: after glDrawElements err=%d\n");

  vbo->unbind();
  ibo->unbind();

  if (I->isPicking) {
    VertexBuffer* pickvbo =
        I->G->ShaderMgr->getGPUBuffer<VertexBuffer>(sp->pickvboid);
    if (pickvbo)
      pickvbo->unbind();
  }

  CheckGLErrorOK(I->G, "CGO_gl_draw_buffers_indexed: end err=%d\n");
}

static void CGO_gl_draw_buffers_not_indexed(CCGORenderer* I, CGO_op_data pc)
{
  const cgo::draw::buffers_not_indexed* sp =
      reinterpret_cast<decltype(sp)>(*pc);
  int mode = sp->mode;

  auto shaderPrg = I->G->ShaderMgr->Get_Current_Shader();
  if (!shaderPrg) {
    return;
  }
  VertexBuffer* vbo = I->G->ShaderMgr->getGPUBuffer<VertexBuffer>(sp->vboid);
  if (!vbo)
    return;
  if (I->isPicking) {
    int attr_a_Color = shaderPrg->GetAttribLocation("a_Color");
    vbo->maskAttributes({attr_a_Color});
    shaderPrg->Set1i("fog_enabled", 0);
    shaderPrg->Set1i("lighting_enabled", 0);
    if (I->use_shader) {
      if (sp->pickvboid) {
        VertexBuffer* pickvbo =
            I->G->ShaderMgr->getGPUBuffer<VertexBuffer>(sp->pickvboid);
        pickvbo->bind(shaderPrg->id, I->pick_pass());
      } else {
        glEnableVertexAttribArray(attr_a_Color);
        glVertexAttribPointer(attr_a_Color, VERTEX_COLOR_SIZE, GL_UNSIGNED_BYTE,
            GL_TRUE, 0, sp->floatdata);
      }
    }
  }

  if (I->debug) {
    mode = CGOConvertDebugMode(I->debug, mode);
  }

  vbo->bind(shaderPrg->id);
  glDrawArrays(mode, 0, sp->nverts);
  vbo->unbind();

  if (I->isPicking) {
    VertexBuffer* pickvbo =
        I->G->ShaderMgr->getGPUBuffer<VertexBuffer>(sp->pickvboid);
    if (pickvbo)
      pickvbo->unbind();
  }
}

static void CGO_gl_mask_attribute_if_picking(CCGORenderer* I, CGO_op_data pc)
{
  if (I->isPicking) {
    const cgo::draw::mask_attribute_if_picking* sp =
        reinterpret_cast<decltype(sp)>(*pc);
    auto shaderPrg = I->G->ShaderMgr->Get_Current_Shader();
    if (!shaderPrg) {
      return;
    }
    VertexBuffer* vbo = I->G->ShaderMgr->getGPUBuffer<VertexBuffer>(sp->vboid);
    if (!vbo)
      return;
    int loc = shaderPrg->GetAttribLocation(
        I->G->ShaderMgr->GetAttributeName(sp->attr_lookup_idx));
    vbo->maskAttribute(loc);
  }
}

static void CGO_gl_bind_vbo_for_picking(CCGORenderer* I, CGO_op_data pc)
{
  if (I->isPicking) {
    const cgo::draw::bind_vbo_for_picking* sp =
        reinterpret_cast<decltype(sp)>(*pc);
    auto shaderPrg = I->G->ShaderMgr->Get_Current_Shader();
    if (!shaderPrg) {
      return;
    }
    VertexBuffer* vbo = I->G->ShaderMgr->getGPUBuffer<VertexBuffer>(sp->vboid);
    if (!vbo)
      return;
    vbo->bind(
        shaderPrg->id, sp->which_attr_idx + sp->npickattrs * I->pick_pass());
  }
}

static void CGO_gl_draw_custom(CCGORenderer* I, CGO_op_data pc)
{
  const cgo::draw::custom* sp = reinterpret_cast<decltype(sp)>(*pc);

  auto shaderPrg = I->G->ShaderMgr->Get_Current_Shader();
  if (!shaderPrg) {
    return;
  }
  VertexBuffer* vbo = I->G->ShaderMgr->getGPUBuffer<VertexBuffer>(sp->vboid);
  if (!vbo)
    return;
  IndexBuffer* ibo = nullptr;
  if (sp->iboid) {
    ibo = I->G->ShaderMgr->getGPUBuffer<IndexBuffer>(sp->iboid);
  }
  vbo->bind(shaderPrg->id);
  if (ibo) {
    ibo->bind();
    glDrawElements(sp->mode, sp->nindices, VertexIndex_GL_ENUM, 0);
  } else {
    glDrawArrays(sp->mode, 0, sp->nverts);
  }
  vbo->unbind();
  if (sp->pickvboid) {
    VertexBuffer* pickvbo =
        I->G->ShaderMgr->getGPUBuffer<VertexBuffer>(sp->pickvboid);
    if (pickvbo)
      pickvbo->unbind();
  }
  if (ibo)
    ibo->unbind();
}

static void CGO_gl_draw_sphere_buffers(CCGORenderer* I, CGO_op_data pc)
{
  const cgo::draw::sphere_buffers* sp = reinterpret_cast<decltype(sp)>(*pc);
  int num_spheres = sp->num_spheres;
  int attr_color;
  VertexBuffer* vbo = I->G->ShaderMgr->getGPUBuffer<VertexBuffer>(sp->vboid);
  VertexBuffer* pickvbo =
      I->G->ShaderMgr->getGPUBuffer<VertexBuffer>(sp->pickvboid);
  CShaderPrg* shaderPrg;
  int pickable = 0;

  shaderPrg = I->G->ShaderMgr->Get_DefaultSphereShader(
      I->info ? I->info->pass : RenderPass::Antialias);
  if (!shaderPrg) {
    return;
  }

  attr_color = shaderPrg->GetAttribLocation("a_Color");

  if (I->isPicking) {
    vbo->maskAttributes({attr_color});
    pickable = SettingGet_i(I->G, I->set1, I->set2, cSetting_pickable);
    shaderPrg->Set1i("lighting_enabled", 0);
    if (pickable) {
      pickvbo->bind(shaderPrg->id, I->pick_pass());
    } else {
      assert(I->info->pick);
      unsigned char nopick[4] = {};
      I->info->pick->colorNoPick(nopick);
      glVertexAttrib4ubv(attr_color, nopick);
    }
  }

  vbo->bind(shaderPrg->id);
#if defined(PURE_OPENGL_ES_2)
  glDrawArrays(GL_TRIANGLES, 0, num_spheres * VerticesPerSphere());
#else
  glDrawArrays(GL_QUADS, 0, num_spheres * 4);
#endif

  vbo->unbind();
}

static void CGO_gl_draw_bezier_buffers(CCGORenderer* I, CGO_op_data cgo_data)
{
  const auto bezier =
      reinterpret_cast<const cgo::draw::bezier_buffers*>(*cgo_data);
  const auto vbo = I->G->ShaderMgr->getGPUBuffer<VertexBuffer>(bezier->vboid);
  auto shaderPrg = I->G->ShaderMgr->Get_BezierShader();
  if (!shaderPrg) {
    return;
  }

  vbo->bind(shaderPrg->id);
  glDrawArrays(GL_PATCHES, 0, 4);
  vbo->unbind();
}

static void CGO_gl_draw_cylinder_buffers(CCGORenderer* I, CGO_op_data pc)
{
  const cgo::draw::cylinder_buffers* sp = reinterpret_cast<decltype(sp)>(*pc);
  int num_cyl = sp->num_cyl;
  int min_alpha = sp->alpha;
  int attr_colors, attr_colors2;
  CShaderPrg* shaderPrg;
  int pickable = 0;
  VertexBuffer* vbo = I->G->ShaderMgr->getGPUBuffer<VertexBuffer>(sp->vboid);
  IndexBuffer* ibo = I->G->ShaderMgr->getGPUBuffer<IndexBuffer>(sp->iboid);
  VertexBuffer* pickvbo =
      I->G->ShaderMgr->getGPUBuffer<VertexBuffer>(sp->pickvboid);

  shaderPrg = I->G->ShaderMgr->Get_CylinderShader(
      I->info ? I->info->pass : RenderPass::Antialias);

  if (!shaderPrg) {
    return;
  }
  attr_colors = shaderPrg->GetAttribLocation("a_Color");
  attr_colors2 = shaderPrg->GetAttribLocation("a_Color2");

  if (I->isPicking) {
    pickable = SettingGet_i(I->G, I->set1, I->set2, cSetting_pickable);
    shaderPrg->Set1i("lighting_enabled", 0);
  }
  if (I->isPicking) {
    vbo->maskAttributes({attr_colors, attr_colors2});
    if (pickable) {
      // in first pass: 1st half of vbo, in second pass: 2nd half of vbo
      // first color (offset 0)
      pickvbo->bind(shaderPrg->id, I->pick_pass());
      // second color (offset 4)
      pickvbo->bind(shaderPrg->id, I->pick_pass() + SHADER_PICKING_PASSES_MAX);
    } else {
      assert(I->info->pick);
      unsigned char nopick[4] = {};
      I->info->pick->colorNoPick(nopick);
      glVertexAttrib4ubv(attr_colors, nopick);
      glVertexAttrib4ubv(attr_colors2, nopick);
    }
  }

  vbo->bind(shaderPrg->id);
  ibo->bind();

  if (min_alpha < 255) {
    glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
    glDrawElements(GL_TRIANGLES, num_cyl * NumTotalVerticesPerCylinder(),
        VertexIndex_GL_ENUM, 0);
    glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
    glDepthFunc(GL_LEQUAL);
  }
  glDrawElements(GL_TRIANGLES, num_cyl * NumTotalVerticesPerCylinder(),
      VertexIndex_GL_ENUM, 0);

  if (min_alpha < 255) {
    glDepthFunc(GL_LESS);
  }

  ibo->unbind();
  vbo->unbind();
  if (I->isPicking)
    pickvbo->unbind();
}
#include "Texture.h"

static void CGO_gl_draw_labels(CCGORenderer* I, CGO_op_data pc)
{
  const cgo::draw::labels* sp = reinterpret_cast<decltype(sp)>(*pc);

  CShaderPrg* shaderPrg;
  int t_mode = SettingGetGlobal_i(I->G, cSetting_transparency_mode);

  if (t_mode == 3 && I->info && I->info->pass != RenderPass::Transparent) {
    // in transparency_mode=3, labels are drawn in the transparency pass=-1
    return;
  }
  shaderPrg = I->G->ShaderMgr->Get_LabelShader(
      I->info ? I->info->pass : RenderPass::Antialias);
  if (I->rep) {
    float label_size;
    CSetting *set1 = nullptr, *set2 = nullptr;
    if (I->rep->cs)
      set1 = I->rep->cs->Setting.get();
    if (I->rep->obj)
      set2 = I->rep->obj->Setting.get();
    label_size = SettingGet_f(I->G, set1, set2, cSetting_label_size);
    shaderPrg->Set1f("scaleByVertexScale", label_size < 0.f ? 1.f : 0.f);
    if (label_size < 0.f) {
      shaderPrg->Set1f("labelTextureSize",
          (float) -2.f * I->info->texture_font_size / label_size);
    }
  }

  if (!shaderPrg) {
    return;
  }

  VertexBuffer* vbo = I->G->ShaderMgr->getGPUBuffer<VertexBuffer>(sp->vboid);
  VertexBuffer* pickvbo =
      I->G->ShaderMgr->getGPUBuffer<VertexBuffer>(sp->pickvboid);

  if (I->isPicking) {
    pickvbo->bind(shaderPrg->id, I->pick_pass());
  }

  if (!vbo)
    return;
  vbo->bind(shaderPrg->id);

  glDrawArrays(GL_TRIANGLES, 0, sp->ntextures * 6);

  vbo->unbind();
  pickvbo->unbind();
}

static void CGO_gl_draw_connectors(CCGORenderer* I, CGO_op_data pc)
{
  int use_geometry_shaders =
      SettingGetGlobal_b(I->G, cSetting_use_geometry_shaders);

  const cgo::draw::connectors* sp = reinterpret_cast<decltype(sp)>(*pc);

  GLenum mode = GL_LINES;
  int factor = 2;
  float lineWidth;
  if (I->isPicking) {
    return;
  }
  CheckGLErrorOK(I->G, "ERROR: CGO_gl_draw_connectors begin returns err=%d\n");

  if (use_geometry_shaders) {
    mode = GL_POINTS;
    factor = 1;
  } else {
    factor = 4;
  }
  auto shaderPrg = I->G->ShaderMgr->Get_Current_Shader();
  if (!shaderPrg) {
    return;
  }
  if (I->rep) {
    float label_size;
    CSetting *set1 = nullptr, *set2 = nullptr;
    float v_scale = SceneGetScreenVertexScale(I->G, nullptr);
    if (I->rep->cs)
      set1 = I->rep->cs->Setting.get();
    if (I->rep->obj)
      set2 = I->rep->obj->Setting.get();
    label_size = SettingGet_f(I->G, set1, set2, cSetting_label_size);
    shaderPrg->Set1f("scaleByVertexScale", label_size < 0.f ? 1.f : 0.f);
    lineWidth = SettingGet_f(I->G, set1, set2, cSetting_label_connector_width);
    if (label_size < 0.f) {
      shaderPrg->Set1f("textureToLabelSize",
          v_scale * (float) I->info->texture_font_size / label_size);
    } else {
      shaderPrg->Set1f("textureToLabelSize", 1.f);
    }
  } else {
    lineWidth = SettingGetGlobal_f(I->G, cSetting_label_connector_width);
  }
#ifndef _WEBGL
  if (!use_geometry_shaders)
    glLineWidth(lineWidth);
#endif

  VertexBuffer* vbo = I->G->ShaderMgr->getGPUBuffer<VertexBuffer>(sp->vboid);
  if (!vbo)
    return;
  vbo->bind(shaderPrg->id);
  glDrawArrays(mode, 0, sp->nconnectors * factor);
  vbo->unbind();
  CheckGLErrorOK(I->G, "ERROR: CGO_gl_draw_connectors end returns err=%d\n");
}

static void CGO_gl_draw_textures(CCGORenderer* I, CGO_op_data pc)
{
  const cgo::draw::textures* sp = reinterpret_cast<decltype(sp)>(*pc);
  int ntextures = sp->ntextures;
  VertexBuffer* vbo = I->G->ShaderMgr->getGPUBuffer<VertexBuffer>(sp->vboid);
  CShaderPrg* shaderPrg;
  int attr_pickcolor = 0;
  shaderPrg = I->G->ShaderMgr->Get_LabelShader(
      I->info ? I->info->pass : RenderPass::Antialias);
  if (!shaderPrg) {
    return;
  }
  if (I->isPicking) {
    attr_pickcolor = shaderPrg->GetAttribLocation("attr_pickcolor");
  }
  if (attr_pickcolor) {
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glEnableVertexAttribArray(attr_pickcolor);
    glVertexAttribPointer(attr_pickcolor, VERTEX_COLOR_SIZE, GL_UNSIGNED_BYTE,
        GL_TRUE, 0, sp->floatdata);
  }
  vbo->bind(shaderPrg->id);
  glDrawArrays(GL_TRIANGLES, 0, ntextures * 6);
  vbo->unbind();
  if (attr_pickcolor) {
    glDisableVertexAttribArray(attr_pickcolor);
  }
}

static void CGO_gl_draw_screen_textures_and_polygons(
    CCGORenderer* I, CGO_op_data pc)
{
  const cgo::draw::screen_textures* sp = reinterpret_cast<decltype(sp)>(*pc);
  int nverts = sp->nverts;
  CShaderPrg* shaderPrg;

  shaderPrg = I->G->ShaderMgr->Get_ScreenShader();
  if (!shaderPrg) {
    return;
  }

  VertexBuffer* vb = I->G->ShaderMgr->getGPUBuffer<VertexBuffer>(sp->vboid);
  if (!vb)
    return;
  vb->bind(shaderPrg->id);

  glDrawArrays(GL_TRIANGLES, 0, nverts);

  vb->unbind();
}

static void CGO_gl_draw_trilines(CCGORenderer* I, CGO_op_data pc)
{
  int nverts = CGO_get_int(*pc);
  int buffer = CGO_get_int(*pc + 1);
  int a_vertex, a_othervertex, a_uv, a_color, a_color2;
  auto shaderPrg = I->G->ShaderMgr->Get_Current_Shader();
  if (!shaderPrg) {
    return;
  }
  a_vertex = 0; // a_Vertex is bound to 0 (see ShaderMgr)
                // CShaderPrg_GetAttribLocation(shaderPrg, "a_Vertex");
  a_othervertex = shaderPrg->GetAttribLocation("a_OtherVertex");
  a_uv = shaderPrg->GetAttribLocation("a_UV");
  a_color = shaderPrg->GetAttribLocation("a_Color");
  a_color2 = shaderPrg->GetAttribLocation("a_Color2");

  glEnableVertexAttribArray(a_vertex);
  glEnableVertexAttribArray(a_othervertex);
  glEnableVertexAttribArray(a_uv);
  glEnableVertexAttribArray(a_color);
  glEnableVertexAttribArray(a_color2);

  glBindBuffer(GL_ARRAY_BUFFER, buffer);

  glVertexAttribPointer(a_vertex, 3, GL_FLOAT, GL_FALSE, 32, (const void*) 0);
  glVertexAttribPointer(
      a_othervertex, 3, GL_FLOAT, GL_FALSE, 32, (const void*) 12);
  glVertexAttribPointer(a_uv, 1, GL_FLOAT, GL_FALSE, 32, (const void*) 24);
  glVertexAttribPointer(
      a_color, 4, GL_UNSIGNED_BYTE, GL_TRUE, 32, (const void*) 28);
  glVertexAttribPointer(
      a_color2, 4, GL_UNSIGNED_BYTE, GL_TRUE, 32, (const void*) 28);
  glDrawArrays(GL_TRIANGLES, 0, nverts);

  glDisableVertexAttribArray(a_vertex);
  glDisableVertexAttribArray(a_othervertex);
  glDisableVertexAttribArray(a_uv);
  glDisableVertexAttribArray(a_color);
  glDisableVertexAttribArray(a_color2);
}

/* CGO_gl_uniform3f - this is the implementation for the
 * CGOUniform3f/CGO_UNIFORM3F operation. From the uniform_id
 * it looks up the uniform location from the current shader,
 * and sets it to the values in this op.
 *
 */
static void CGO_gl_uniform3f(CCGORenderer* I, CGO_op_data pc)
{
  int uniform_id = CGO_get_int(*pc);
  auto shaderPrg = I->G->ShaderMgr->Get_Current_Shader();
  if (!shaderPrg) {
    return;
  }
  int loc = shaderPrg->GetUniformLocation(
      shaderPrg->uniformLocations[uniform_id].c_str());
  const float* pcp = *pc + 1;
  glUniform3f(loc, pcp[0], pcp[1], pcp[2]);
}

static void CGO_gl_linewidth(CCGORenderer* I, CGO_op_data pc)
{
#ifndef _WEBGL
  glLineWidth(**pc);
#endif
}

/**
 * call glLineWidth and set the "line_width" uniform
 */
static void glLineWidthAndUniform(
    float line_width, CShaderPrg* shaderPrg = nullptr)
{
#ifndef _WEBGL
  glLineWidth(line_width);
#endif

  if (shaderPrg && shaderPrg->name == "trilines")
    shaderPrg->Set1f("line_width", line_width);
}

/* CGO_gl_special - this is the implementation function for
   CGOSpecial/CGO_SPECIAL.  Each op has its own implementation.
 */
static void CGO_gl_special(CCGORenderer* I, CGO_op_data pc)
{
  int mode = CGO_get_int(*pc);
  bool openVR = SceneGetStereo(I->G) == cStereo_openvr;
  char varwidth = 0;
  float vScale =
      (I->info ? I->info->vertex_scale : SceneGetScreenVertexScale(I->G, nullptr));

  CSetting *csSetting = nullptr, *objSetting = nullptr;
  auto shaderPrg = I->G->ShaderMgr->Get_Current_Shader();
  if (I->rep && I->rep->cs) {
    csSetting = I->rep->cs->Setting.get();
  }
  if (I->rep && I->rep->obj) {
    objSetting = I->rep->obj->Setting.get();
  }
  switch (mode) {
  case LINEWIDTH_DYNAMIC_WITH_SCALE_RIBBON: {
    float line_width = SettingGet_f(I->G, nullptr, nullptr, cSetting_ribbon_width);
    if (!openVR)
      line_width = SceneGetDynamicLineWidth(I->info, line_width);
    if (I->info && I->info->width_scale_flag) {
      line_width *= I->info->width_scale;
    }
    glLineWidthAndUniform(line_width, shaderPrg);
  } break;
  case LINEWIDTH_DYNAMIC_WITH_SCALE_DASH: {
    float line_width = SettingGet_f(I->G, nullptr, nullptr, cSetting_dash_width);
    if (!openVR)
      line_width = SceneGetDynamicLineWidth(I->info, line_width);
    if (I->info && I->info->width_scale_flag) {
      line_width *= I->info->width_scale;
    }
    glLineWidthAndUniform(line_width, shaderPrg);
  } break;
  case LINEWIDTH_DYNAMIC_WITH_SCALE: {
    float line_width = SettingGet_f(I->G, nullptr, nullptr, cSetting_line_width);
    if (!openVR)
      line_width = SceneGetDynamicLineWidth(I->info, line_width);
    if (I->info && I->info->width_scale_flag) {
      line_width *= I->info->width_scale;
    }
    glLineWidthAndUniform(line_width, shaderPrg);
  } break;
  case LINEWIDTH_WITH_SCALE: {
    float line_width = SettingGet_f(I->G, nullptr, nullptr, cSetting_line_width);
    if (I->info && I->info->width_scale_flag) {
      line_width *= I->info->width_scale;
    }
    glLineWidthAndUniform(line_width, shaderPrg);
  } break;
  case LINEWIDTH_DYNAMIC_MESH: {
    float line_width;
    if (I->rep) {
      line_width = SettingGet_f(I->G, I->rep->cs->Setting.get(),
          I->rep->obj->Setting.get(), cSetting_mesh_width);
    } else {
      line_width = SettingGet_f(I->G, nullptr, nullptr, cSetting_mesh_width);
    }
    if (!openVR)
      line_width = SceneGetDynamicLineWidth(I->info, line_width);
    glLineWidthAndUniform(line_width, shaderPrg);
  } break;
  case POINTSIZE_DYNAMIC_DOT_WIDTH: {
    float ps;
    if (I->info && I->info->width_scale_flag) {
      ps = SettingGet_f(I->G, csSetting, objSetting, cSetting_dot_width) *
           I->info->width_scale;
    } else {
      ps = SettingGet_f(I->G, csSetting, objSetting, cSetting_dot_width);
    }
#ifdef PURE_OPENGL_ES_2
    if (I->G->ShaderMgr->current_shader) {
      shaderPrg->Set1f("g_pointSize", ps);
    }
#else
    glPointSize(ps);
#endif
    break;
  }
  case CYLINDERWIDTH_DYNAMIC_MESH: {
    CSetting* setting = nullptr;
    float mesh_width;
    if (I && I->rep && I->rep->obj) {
      setting = I->rep->obj->Setting.get();
    }
    mesh_width = SettingGet_f(I->G, setting, nullptr, cSetting_mesh_width);
    if (shaderPrg) {
      const float* color = I->color ? I->color : g_ones4f;
      shaderPrg->Set1f("uni_radius",
          SceneGetLineWidthForCylinders(I->G, I->info, mesh_width));
      shaderPrg->SetAttrib4fLocation(
          "a_Color", color[0], color[1], color[2], I->alpha);
      shaderPrg->SetAttrib4fLocation(
          "a_Color2", color[0], color[1], color[2], I->alpha);
    }
  } break;
  case DOTSIZE_WITH_SPHERESCALE: {
    float radius =
        SettingGet_f(I->G, csSetting, objSetting, cSetting_dot_width);
    radius *= vScale;
    if (shaderPrg)
      shaderPrg->Set1f("sphere_size_scale", fabs(radius));
  } break;
  case MESH_WIDTH_FOR_SURFACES: {
    float mesh_width =
        SettingGet_f(I->G, csSetting, objSetting, cSetting_mesh_width);
    if (shaderPrg)
      shaderPrg->Set1f("uni_radius",
          SceneGetLineWidthForCylinders(I->G, I->info, mesh_width));
  } break;
  case CYLINDER_WIDTH_FOR_DISTANCES: {
    float line_width, radius;
    int round_ends;
    round_ends =
        SettingGet_b(I->G, csSetting, objSetting, cSetting_dash_round_ends);
    line_width = SettingGet_f(I->G, csSetting, objSetting, cSetting_dash_width);
    radius = SettingGet_f(I->G, csSetting, objSetting, cSetting_dash_radius);

    line_width = SceneGetDynamicLineWidth(I->info, line_width);

    if (shaderPrg) {
      if (radius == 0.0F) {
        float dash_size =
            SettingGet_f(I->G, csSetting, objSetting, cSetting_dash_width);
        shaderPrg->Set1f("uni_radius", SceneGetLineWidthForCylindersStatic(I->G,
                                           I->info, line_width, dash_size));
      } else {
        shaderPrg->Set1f("uni_radius", radius);
      }
      if (!round_ends) {
        shaderPrg->Set1i("no_flat_caps", 0);
      }
    }
  } break;
  case CYLINDER_WIDTH_FOR_RIBBONS: {
    float pixel_scale_value =
        SettingGetGlobal_f(I->G, cSetting_ray_pixel_scale);
    float line_width, radius;
    line_width =
        SettingGet_f(I->G, csSetting, objSetting, cSetting_ribbon_width);
    radius = SettingGet_f(I->G, csSetting, objSetting, cSetting_ribbon_radius);

    line_width = SceneGetDynamicLineWidth(I->info, line_width);
    if (pixel_scale_value < 0)
      pixel_scale_value = 1.0F;
    if (shaderPrg) {
      if (radius == 0.0F) {
        shaderPrg->Set1f(
            "uni_radius", vScale * pixel_scale_value * line_width / 2.f);
      } else {
        shaderPrg->Set1f("uni_radius", radius);
      }
    }
  } break;
  case DOT_WIDTH_FOR_DOTS: {
    float dot_width =
        SettingGet_f(I->G, csSetting, objSetting, cSetting_dot_width);
    float radius;
    if (I->info && I->info->width_scale_flag)
      radius = (dot_width * I->info->width_scale);
    else
      radius = dot_width;
    if (shaderPrg)
      shaderPrg->Set1f("g_PointSize", radius);
    glPointSize(radius);
  } break;
  case DOT_WIDTH_FOR_DOT_SPHERES: {
    float dotSize =
        SettingGet_f(I->G, csSetting, objSetting, cSetting_dot_radius);
    float dot_width =
        SettingGet_f(I->G, csSetting, objSetting, cSetting_dot_width);
    float radius;
    if (I->info && dotSize <= 0.0F) {
      if (I->info->width_scale_flag)
        radius =
            dot_width * I->info->width_scale * I->info->vertex_scale / 1.4142F;
      else
        radius = dot_width * I->info->vertex_scale;
    } else {
      radius = dotSize;
    }
    if (shaderPrg)
      shaderPrg->Set1f("sphere_size_scale", fabs(radius));
  } break;
  case CYLINDER_WIDTH_FOR_NONBONDED: {
    if (shaderPrg) {
      float line_width =
          SettingGet_f(I->G, csSetting, objSetting, cSetting_line_width);
      shaderPrg->Set1f("uni_radius", SceneGetLineWidthForCylindersStatic(I->G,
                                         I->info, line_width, line_width));
    }
  } break;
  case CYLINDER_WIDTH_FOR_REPWIRE_VARWIDTH:
    varwidth = 1;
  case CYLINDER_WIDTH_FOR_REPWIRE: {
    float radius =
        SettingGet_f(I->G, csSetting, objSetting, cSetting_line_radius);
    if (radius < R_SMALL8) {
      float line_width =
          SettingGet_f(I->G, csSetting, objSetting, cSetting_line_width);
      float pixel_scale_value =
          SettingGetGlobal_f(I->G, cSetting_ray_pixel_scale);
      float vertex_scale = vScale;
      float scale_bound = SettingGetGlobal_f(I->G, cSetting_field_of_view) *
                          cPI / 180.0f * 0.018f;
      if (!varwidth) {
        line_width = SceneGetDynamicLineWidth(I->info, line_width);
      }
      if (vertex_scale < scale_bound) {
        vertex_scale = scale_bound;
      }
      if (pixel_scale_value < 0)
        pixel_scale_value = 1.0F;
      radius = vertex_scale * pixel_scale_value * line_width / 2.f;
    }
    if (shaderPrg) {
      shaderPrg->Set1f("uni_radius", radius);
    }
  } break;
  case ENABLE_BACK_FACES_IF_NOT_TWO_SIDED: {
    int two_sided_lighting = SettingGet_i(I->G, csSetting, objSetting,
                                 cSetting_two_sided_lighting) > 0;
    if (!two_sided_lighting) {
      glCullFace(GL_BACK);
      glEnable(GL_CULL_FACE);
    }
  } break;
  case DISABLE_BACK_FACES_IF_NOT_TWO_SIDED: {
    int two_sided_lighting = SettingGet_i(I->G, csSetting, objSetting,
                                 cSetting_two_sided_lighting) > 0;
    if (!two_sided_lighting) {
      glDisable(GL_CULL_FACE);
    }
  } break;
  case SET_SURFACE_UNIFORMS: {
    float ambient_occlusion_scale = 0.f;
    int ambient_occlusion_mode = SettingGet_i(
        I->G, csSetting, objSetting, cSetting_ambient_occlusion_mode);

    if (ambient_occlusion_mode) {
      ambient_occlusion_scale = SettingGet_f(
          I->G, csSetting, objSetting, cSetting_ambient_occlusion_scale);
    }
    if (shaderPrg)
      shaderPrg->Set1f("ambient_occlusion_scale", ambient_occlusion_scale);
  } break;
  case SET_ALIGNMENT_UNIFORMS_ATTRIBS: {
    float linewidth =
        SettingGet_f(I->G, csSetting, objSetting, cSetting_cgo_line_width);
    float lineradius =
        SettingGet_f(I->G, csSetting, objSetting, cSetting_cgo_line_radius);
    float pixel_scale_value =
        SettingGetGlobal_f(I->G, cSetting_ray_pixel_scale);
    if (linewidth < 0.f) {
      linewidth = 1.f;
    }
    if (pixel_scale_value < 0)
      pixel_scale_value = 1.0F;
    if (lineradius < 0.f) {
      lineradius = linewidth * vScale * pixel_scale_value / 2.f;
    }
    shaderPrg->Set1f("uni_radius", lineradius);
    if (I->color) {
      shaderPrg->SetAttrib4fLocation(
          "a_Color", I->color[0], I->color[1], I->color[2], 1.f);
      shaderPrg->SetAttrib4fLocation(
          "a_Color2", I->color[0], I->color[1], I->color[2], 1.f);
    }
    glLineWidthAndUniform(lineradius * 2.f / vScale, shaderPrg);
  } break;
  case LINEWIDTH_FOR_LINES: {
    float line_width = SceneGetDynamicLineWidth(
        I->info, SettingGet_f(I->G, nullptr, nullptr, cSetting_line_width));
    if (I->info && I->info->width_scale_flag) {
      line_width *= I->info->width_scale;
    }
    glLineWidthAndUniform(line_width, shaderPrg);
  } break;
  case SET_LABEL_SCALE_UNIFORMS: {
    if (I->rep) {
      float label_size;
      CSetting *set1 = nullptr, *set2 = nullptr;
      if (I->rep->cs)
        set1 = I->rep->cs->Setting.get();
      if (I->rep->obj)
        set2 = I->rep->obj->Setting.get();
      label_size = SettingGet_f(I->G, set1, set2, cSetting_label_size);
      shaderPrg->Set1f("scaleByVertexScale", label_size < 0.f ? 1.f : 0.f);
      if (label_size < 0.f) {
        shaderPrg->Set1f("labelTextureSize",
            (float) -2.f * I->info->texture_font_size / label_size);
      }
    }

  } break;
  default:
    PRINTFB(I->G, FB_CGO, FB_Warnings)
    " CGO_gl_special(): bad mode=%d\n", mode ENDFB(I->G);
  }
}

/* CGO_gl_special_with_arg - this is the implementation function for
   CGOSpecialWithArg/CGO_SPECIAL_WITH_ARG.  Each op has its own implementation.
 */
static void CGO_gl_special_with_arg(CCGORenderer* I, CGO_op_data pc)
{
#ifndef PURE_OPENGL_ES_2
  int mode = CGO_get_int(*pc);
  float argval = *((*pc) + 1);
  bool use_shaders = SettingGetGlobal_b(I->G, cSetting_use_shaders);
  bool sphere_use_shaders =
      use_shaders && SettingGetGlobal_b(I->G, cSetting_use_shaders);
  switch (mode) {
  case LINEWIDTH_FOR_LINES: {
    if (!use_shaders) {
      glEnd();
      glLineWidth(argval);
      glBegin(GL_LINES);
    }
  } break;
  case LINE_LIGHTING:
    if (!I->isPicking && !SettingGetGlobal_b(I->G, cSetting_use_shaders)) {
      if (!I->info->line_lighting) {
        bool enableLighting = (int) argval;
        if (enableLighting)
          glEnable(GL_LIGHTING);
        else
          glDisable(GL_LIGHTING);
      }
    }
    break;
  case SPHERE_MODE_OPS: {
    float pixel_scale = 1.0F / I->info->vertex_scale;
    int sphere_mode = (int) fabs(argval);
    bool enable = argval > 0.f;
    if (enable) {
      float pointSize;
      if ((sphere_mode == 1) || (sphere_mode == 6)) {
        pointSize =
            SettingGet_f(I->G, I->set1, I->set2, cSetting_sphere_point_size);
        glDisable(GL_POINT_SMOOTH);
        glDisable(GL_ALPHA_TEST);
        if (!I->isPicking && !sphere_use_shaders) {
          glEnable(GL_LIGHTING);
          glHint(GL_POINT_SMOOTH_HINT, GL_FASTEST);
        }
      } else {
        float sphere_scale =
            SettingGet_f(I->G, I->set1, I->set2, cSetting_sphere_scale);
        if ((sphere_mode == 3) || (sphere_mode == 8)) {
          glEnable(GL_POINT_SMOOTH);
          glAlphaFunc(GL_GREATER, 0.5F);
          glEnable(GL_ALPHA_TEST);
          glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
          pointSize = sphere_scale * pixel_scale * 2.0F;
        } else {
          glHint(GL_POINT_SMOOTH_HINT, GL_FASTEST);
          glDisable(GL_POINT_SMOOTH);
          glDisable(GL_ALPHA_TEST);
          pointSize = sphere_scale * pixel_scale * 1.4F;
        }
      }
      if (!I->isPicking && ((sphere_mode == 7) || (sphere_mode == 8)))
        glEnable(GL_LIGHTING);
      glPointSize(pointSize);
    } else {
      if (sphere_mode == 3) {
        glDisable(GL_POINT_SMOOTH);
        glAlphaFunc(GL_GREATER, 0.05F);
      } else {
        glEnable(GL_ALPHA_TEST);
      }
    }
  }
  }
#endif
}

static void CGO_gl_dotwidth(CCGORenderer* I, CGO_op_data pc)
{
#ifndef PURE_OPENGL_ES_2
  glPointSize(**pc);
#endif
}

static void CGO_gl_enable(CCGORenderer* I, CGO_op_data pc)
{
  GLenum mode = CGO_get_int(*pc);
  CShaderMgr* shaderMgr = I->G->ShaderMgr;
  CShaderPrg* shaderPrg = shaderMgr->Get_Current_Shader();
  if (I->use_shader) {
    if (true) {
      switch (mode) {
      case CGO_GL_LIGHTING: {
        if (shaderPrg) {
          shaderPrg->SetLightingEnabled(1);
        }
      } break;
      case GL_SHADER_LIGHTING:
        if (!I->isPicking) {
          if (shaderPrg) {
            shaderPrg->SetLightingEnabled(1);
          }
        }
        break;
      case GL_TWO_SIDED_LIGHTING: {
        if (shaderPrg) {
          shaderPrg->Set1i("two_sided_lighting_enabled", 1);
        }
      } break;
      case GL_MESH_LIGHTING: {
        int lighting =
            SettingGet_i(I->G, I->set1, I->set2, cSetting_mesh_lighting);
        if (shaderPrg) {
          shaderPrg->SetLightingEnabled(lighting);
        }
      } break;
      case GL_DOT_LIGHTING: {
        int lighting =
            SettingGet_i(I->G, I->set1, I->set2, cSetting_dot_lighting);
        if (shaderPrg && !I->isPicking) {
          shaderPrg->SetLightingEnabled(lighting);
          shaderPrg->Set1i("two_sided_lighting_enabled", 0);
        }
      } break;
      case GL_LABEL_FLOAT_TEXT: {
        int float_text =
            SettingGet_i(I->G, I->set1, I->set2, cSetting_float_labels);
        if (float_text) {
          glDisable(GL_DEPTH_TEST);
        }
      } break;
      case GL_DASH_TRANSPARENCY_DEPTH_TEST: {
        float dash_transparency =
            SettingGet_f(I->G, I->set1, I->set2, cSetting_dash_transparency);
        short dash_transparency_enabled;
        bool t_mode_3 = SettingGet_i(I->G, I->set1, I->set2,
                            cSetting_transparency_mode) == 3;
        dash_transparency =
            (dash_transparency < 0.f
                    ? 0.f
                    : (dash_transparency > 1.f ? 1.f : dash_transparency));
        dash_transparency_enabled = (dash_transparency > 0.f);
        if (dash_transparency_enabled && !t_mode_3 && !I->isPicking) {
          glDisable(GL_DEPTH_TEST);
        }
      } break;
      case GL_DEFAULT_SHADER:
        shaderMgr->Enable_DefaultShader(
            I->info ? I->info->pass : RenderPass::Antialias);
        break;
      case GL_LINE_SHADER:
        shaderMgr->Enable_LineShader(
            I->info ? I->info->pass : RenderPass::Antialias);
        break;
      case GL_SURFACE_SHADER:
        shaderMgr->Enable_SurfaceShader(
            I->info ? I->info->pass : RenderPass::Antialias);
        break;
      case GL_CYLINDER_SHADER:
        shaderMgr->Enable_CylinderShader(
            I->info ? I->info->pass : RenderPass::Antialias);
        break;
      case GL_SPHERE_SHADER:
        shaderMgr->Enable_DefaultSphereShader(
            I->info ? I->info->pass : RenderPass::Antialias);
        break;
      case GL_RAMP_SHADER:
        shaderMgr->Enable_RampShader();
        break;
      case GL_DEFAULT_SHADER_WITH_SETTINGS:
        shaderMgr->Enable_DefaultShaderWithSettings(
            I->set1, I->set2, I->info ? I->info->pass : RenderPass::Antialias);
        break;
      case GL_BACKGROUND_SHADER:
        shaderMgr->Enable_BackgroundShader();
        break;
      case GL_LABEL_SHADER:
        shaderMgr->Enable_LabelShader(
            I->info ? I->info->pass : RenderPass::Antialias);
        break;
      case GL_CONNECTOR_SHADER:
        shaderMgr->Enable_ConnectorShader(
            I->info ? I->info->pass : RenderPass::Antialias);
        break;
      case GL_SCREEN_SHADER:
        shaderMgr->Enable_ScreenShader();
        break;
      case GL_TRILINES_SHADER:
        shaderMgr->Enable_TriLinesShader();
        break;
#ifndef _PYMOL_NO_AA_SHADERS
      case GL_FXAA_SHADER:
        shaderMgr->Enable_FXAAShader();
        break;
      case GL_SMAA1_SHADER:
        shaderMgr->Enable_SMAA1Shader();
        break;
      case GL_SMAA2_SHADER:
        shaderMgr->Enable_SMAA2Shader();
        break;
      case GL_SMAA3_SHADER:
        shaderMgr->Enable_SMAA3Shader();
        break;
#endif
      case GL_OIT_SHADER:
        shaderMgr->Enable_OITShader();
        break;
      case GL_OIT_COPY_SHADER:
        shaderMgr->Enable_OITCopyShader();
        break;
      case GL_BACK_FACE_CULLING:
        glCullFace(GL_BACK);
        glEnable(GL_CULL_FACE);
        break;
      case GL_DEPTH_TEST:
        glEnable(mode);
        break;
      case GL_DEPTH_TEST_IF_FLOATING: {
        int float_text =
            SettingGet_i(I->G, I->set1, I->set2, cSetting_float_labels);
        if (float_text)
          glEnable(GL_DEPTH_TEST);
      } break;
      case GL_BEZIER_SHADER:
        shaderMgr->Enable_BezierShader();
        break;
      }
    }
  } else {
#ifndef PURE_OPENGL_ES_2
    if (!I->isPicking) {
      if (mode == CGO_GL_LIGHTING) {
        glEnable(GL_LIGHTING);
      }
    }
#endif
  }
}

static void CGO_gl_disable(CCGORenderer* I, CGO_op_data pc)
{
  GLenum mode = CGO_get_int(*pc);
  auto shaderPrg = I->G->ShaderMgr->Get_Current_Shader();
  if (I->use_shader) {
    switch (mode) {
    case GL_SHADER_LIGHTING: {
      if (shaderPrg) {
        shaderPrg->SetLightingEnabled(0);
      }
    } break;
    case GL_CYLINDER_SHADER:
      glDisable(GL_CULL_FACE);
    case GL_RAMP_SHADER:
    case GL_SCREEN_SHADER:
    case GL_LABEL_SHADER:
    case GL_CONNECTOR_SHADER:
    case GL_DEFAULT_SHADER:
    case GL_SURFACE_SHADER:
    case GL_SPHERE_SHADER:
    case GL_TRILINES_SHADER:
    case GL_OIT_COPY_SHADER:
    case GL_LINE_SHADER:
      I->G->ShaderMgr->Disable_Current_Shader();
      break;
    case GL_LABEL_FLOAT_TEXT: {
      int float_text =
          SettingGet_i(I->G, I->set1, I->set2, cSetting_float_labels);
      if (float_text) {
        glEnable(GL_DEPTH_TEST);
      }
    } break;
    case GL_DASH_TRANSPARENCY_DEPTH_TEST: {
      float dash_transparency =
          SettingGet_f(I->G, I->set1, I->set2, cSetting_dash_transparency);
      short dash_transparency_enabled;
      bool t_mode_3 =
          SettingGet_i(I->G, I->set1, I->set2, cSetting_transparency_mode) == 3;
      dash_transparency =
          (dash_transparency < 0.f
                  ? 0.f
                  : (dash_transparency > 1.f ? 1.f : dash_transparency));
      dash_transparency_enabled = (dash_transparency > 0.f);
      if (dash_transparency_enabled && !t_mode_3 && !I->isPicking) {
        glEnable(GL_DEPTH_TEST);
      }
    } break;
    case CGO_GL_LIGHTING: {
      if (shaderPrg) {
        shaderPrg->SetLightingEnabled(0);
      }
    } break;
    case GL_TWO_SIDED_LIGHTING: {
      if (shaderPrg) {
        shaderPrg->Set1i("two_sided_lighting_enabled", 0);
      }
    } break;
#if !defined(PURE_OPENGL_ES_2) || defined(_WEBGL)
    case GL_OIT_SHADER:
    case GL_SMAA1_SHADER:
    case GL_SMAA2_SHADER:
      glBindFramebufferEXT(
          GL_FRAMEBUFFER_EXT, I->G->ShaderMgr->default_framebuffer_id);
      break;
#endif
    case GL_BACK_FACE_CULLING:
      glDisable(GL_CULL_FACE);
      break;
    case GL_DEPTH_TEST:
      glDisable(mode);
      break;
    case GL_DEPTH_TEST_IF_FLOATING: {
      int float_text =
          SettingGet_i(I->G, I->set1, I->set2, cSetting_float_labels);
      if (float_text)
        glDisable(GL_DEPTH_TEST);
    } break;
    }
#ifndef PURE_OPENGL_ES_2
  } else {
    if (mode != CGO_GL_LIGHTING || !I->isPicking) {
      if (mode == CGO_GL_LIGHTING)
        mode = GL_LIGHTING;
      glDisable(mode);
    }
  }
#else
  }
#endif
}

static void CGO_gl_alpha(CCGORenderer* I, CGO_op_data pc)
{
  I->alpha = **pc;
}

static void CGO_gl_reset_normal(CCGORenderer* I, CGO_op_data pc)
{
  SceneResetNormalUseShader(I->G, CGO_get_int(*pc), I->use_shader);
}

static void CGO_gl_null(CCGORenderer* I, CGO_op_data pc) {}

static void CGO_gl_error(CCGORenderer* I, CGO_op_data pc)
{
  PRINTFB(I->G, FB_CGO, FB_Warnings)
  " CGO_gl_error() is not suppose to be called op=%d\n",
      CGO_get_int((*pc) - 1) ENDFB(I->G);
}

static void CGO_gl_color(CCGORenderer* I, CGO_op_data varg)
{
  auto* v = *varg;
  if (I->use_shader) {
    auto shaderPrg = I->G->ShaderMgr->Get_Current_Shader();
    if (shaderPrg) {
      int attr_a_Color = shaderPrg->GetAttribLocation("a_Color");
      glVertexAttrib4f(attr_a_Color, v[0], v[1], v[2], I->alpha);
    }
  } else {
#ifndef PURE_OPENGL_ES_2
    glColor4f(v[0], v[1], v[2], I->alpha);
#endif
  }
}

static void CGO_gl_sphere(CCGORenderer* I, CGO_op_data varg)
{
  auto* v = *varg;
  if (I->isPicking) {
    SphereRender(I->G, 0, v, I->color, I->alpha, v[3]);
  } else {
    SphereRender(I->G, I->sphere_quality, v, nullptr, I->alpha, v[3]);
  }
}

static void CGO_gl_vertex_attribute_3f(CCGORenderer* I, CGO_op_data varg)
{
  auto vertex_attr =
      reinterpret_cast<const cgo::draw::vertex_attribute_3f*>(*varg);
  auto shaderPrg = I->G->ShaderMgr->Get_Current_Shader();
  int loc = shaderPrg->GetAttribLocation(
      I->G->ShaderMgr->GetAttributeName(vertex_attr->attr_lookup_idx));
  if (loc >= 0)
    glVertexAttrib3fv(loc, vertex_attr->values);
}

static void CGO_gl_vertex_attribute_4ub(CCGORenderer* I, CGO_op_data varg)
{
  auto vertex_attr =
      reinterpret_cast<const cgo::draw::vertex_attribute_4ub*>(*varg);
  auto shaderPrg = I->G->ShaderMgr->Get_Current_Shader();
  int loc = shaderPrg->GetAttribLocation(
      I->G->ShaderMgr->GetAttributeName(vertex_attr->attr_lookup_idx));
  if (loc >= 0)
    glVertexAttrib4ubv(loc, vertex_attr->ubdata);
}

static void CGO_gl_vertex_attribute_4ub_if_picking(
    CCGORenderer* I, CGO_op_data varg)
{
  if (I->isPicking) {
    auto vertex_attr =
        reinterpret_cast<const cgo::draw::vertex_attribute_4ub_if_picking*>(
            *varg);
    auto shaderPrg = I->G->ShaderMgr->Get_Current_Shader();
    int loc = shaderPrg->GetAttribLocation(
        I->G->ShaderMgr->GetAttributeName(vertex_attr->attr_lookup_idx));
    if (loc >= 0)
      glVertexAttrib4ubv(loc, vertex_attr->ubdata);
  }
}

static void CGO_gl_vertex_attribute_1f(CCGORenderer* I, CGO_op_data varg)
{
  auto vertex_attr =
      reinterpret_cast<const cgo::draw::vertex_attribute_1f*>(*varg);
  auto shaderPrg = I->G->ShaderMgr->Get_Current_Shader();
  const char* name =
      I->G->ShaderMgr->GetAttributeName(vertex_attr->attr_lookup_idx);
  int loc = shaderPrg->GetAttribLocation(name);
  if (loc >= 0)
    glVertexAttrib1f(loc, vertex_attr->value);
}

/* dispatch table for OpenGL */

CGO_op_fn CGO_gl[] = {CGO_gl_null, /* 0x00 */
    CGO_gl_null,                   /* 0x01 */
    CGO_gl_begin,                  /* 0x02 */
    CGO_gl_end,                    /* 0x03 */
    CGO_gl_vertex,                 /* 0x04 */
    CGO_gl_normal,                 /* 0x05 */
    CGO_gl_color,                  /* 0x06 */
    CGO_gl_sphere,                 /* 0x07 */
    CGO_gl_null,                   /* 0x08 */
    CGO_gl_null,                   /* 0x09 */

    CGO_gl_linewidth, /* 0x0A */
    CGO_gl_null,      /* 0x0B */
    CGO_gl_enable,    /* 0x0C */
    CGO_gl_disable,   /* 0x0D */
    CGO_gl_null,      /* 0x0E */
    CGO_gl_null,      /* 0x0F */

    CGO_gl_dotwidth, /* 0X10 */
    CGO_gl_null,     /* 0x11 */
    CGO_gl_null,     /* 0x12 */
    CGO_gl_null,     /* 0X13 */

    CGO_gl_null, /* 0X14 */
    CGO_gl_null, /* 0x15 */
    CGO_gl_null, /* 0x16 */
    CGO_gl_null, /* 0X17 */

    CGO_gl_null,                     /* 0X18 */
    CGO_gl_alpha,                    /* 0x19 */
    CGO_gl_null,                     /* 0x1A */
    CGO_gl_null,                     /* 0X1B */
    CGO_gl_draw_arrays,              /* 0x1C DrawArrays() */
    CGO_gl_null,                     /* 0x1D */
    CGO_gl_reset_normal,             /* 0x1E */
    CGO_gl_null,                     /* pick color  0X1F */
    CGO_gl_null,                     /* 0x20 draw buffers REMOVED */
    CGO_gl_draw_buffers_indexed,     /* 0x21 draw buffers indexed */
    CGO_gl_null,                     /* 0x22 bounding box */
    CGO_gl_draw_buffers_not_indexed, /* 0x23 draw buffers not indexed */
    CGO_gl_special,                  /* 0x24 special */
    CGO_gl_draw_cylinder_buffers,    /* 0x25 draw GLSL cylinders */
    CGO_gl_null,                     /* 0x26 shader cylinder */
    CGO_gl_null,                     /* 0x27 shader cylinder with 2nd color */
    CGO_gl_draw_sphere_buffers,      /* 0x28 draw sphere buffers */
    CGO_gl_null,          /* 0x29 accessibility used for ambient occlusion */
    CGO_gl_error,         /* 0x2A draw texture */
    CGO_gl_draw_textures, /* 0x2B draw textures */
    CGO_gl_draw_screen_textures_and_polygons, /* 0x2C draw screen textures and
                                                 polygons */
    CGO_gl_error, CGO_gl_error, CGO_gl_draw_labels, CGO_gl_error,
    CGO_gl_draw_connectors, CGO_gl_draw_trilines, CGO_gl_uniform3f,
    CGO_gl_special_with_arg, CGO_gl_line, CGO_gl_splitline, CGO_gl_draw_custom,
    CGO_gl_vertex_attribute_3f, CGO_gl_vertex_attribute_4ub,
    CGO_gl_vertex_attribute_1f, CGO_gl_mask_attribute_if_picking,
    CGO_gl_bind_vbo_for_picking, CGO_gl_vertex,
    CGO_gl_null,         // interpolated
    CGO_gl_vertex_cross, // CGO_VERTEX_CROSS
    CGO_gl_vertex_attribute_4ub_if_picking,
    CGO_gl_null, // custom cylinder alpha
    CGO_gl_null, // bezier
    CGO_gl_draw_bezier_buffers,
    CGO_gl_error};

#if 0
static
void SetUCColorToPrev(uchar *color){
  color[0] = color[-4];
  color[1] = color[-3];
  color[2] = color[-2];
  color[3] = color[-1];
}

static
void SetUCColorToPrev8(uchar *color){
  color[0] = color[-8];
  color[1] = color[-7];
  color[2] = color[-6];
  color[3] = color[-5];
}
#endif

static void SetUCColorToPrevN(int n, uchar* color)
{
  color[0] = color[-n * 4];
  color[1] = color[-n * 4 + 1];
  color[2] = color[-n * 4 + 2];
  color[3] = color[-n * 4 + 3];
}

static int* get_pickcolorsset_ptr(int op, float* pc)
{
#define RETURN_PICKCOLORSETPTR_CASE(cls)                                       \
  case cgo::draw::cls::op_code:                                                \
    return &(reinterpret_cast<cgo::draw::cls*>(pc)->pickcolorsset)
  switch (op) {
    RETURN_PICKCOLORSETPTR_CASE(buffers_indexed);
    RETURN_PICKCOLORSETPTR_CASE(buffers_not_indexed);
    RETURN_PICKCOLORSETPTR_CASE(labels);
    RETURN_PICKCOLORSETPTR_CASE(sphere_buffers);
    RETURN_PICKCOLORSETPTR_CASE(cylinder_buffers);
    RETURN_PICKCOLORSETPTR_CASE(custom);
  }
  return nullptr;
}

void CGORenderGLPicking(CGO* I, RenderInfo* info, PickContext* context,
    CSetting* set1, CSetting* set2, Rep* rep)
{
  PyMOLGlobals* G = I->G;

  if (!G->ValidContext)
    return;

  if (!I->c)
    return;

  CCGORenderer* R = G->CGORenderer;
  bool pickable =
      (!I->no_pick) && SettingGet_b(G, set1, set2, cSetting_pickable);
  auto pick = info->pick;
  bool reset_colors = !pick->pickColorsValid();

  R->use_shader = I->use_shader;
  R->isPicking = true;
  R->set1 = set1;
  R->set2 = set2;
  R->info = info;
  R->rep = rep;

#ifndef _WEBGL
  glLineWidth(SettingGet_f(G, set1, set2, cSetting_cgo_line_width));
#endif

  for (auto it = I->begin(); !it.is_stop(); ++it) {
    const auto op = it.op_code();
    CGO_OP_DATA_CONST float* pc = it.data();

    switch (op) {
    case CGO_COLOR:
      continue;

    case CGO_PICK_COLOR:

      if (reset_colors) { // only if picking info is invalid
        unsigned char col[4];
        AssignNewPickColor(I, pick, col, context, CGO_get_uint(pc),
            pickable ? CGO_get_int(pc + 1) : cPickableNoPick);
#ifndef PURE_OPENGL_ES_2
        if (!I->use_shader) {
          glColor4ubv(col);
        }
#endif
      } else {
        PRINTFB(G, FB_CGO, FB_Warnings)
        " %s: unexpected CGO_PICK_COLOR with !reset_colors\n",
            __func__ ENDFB(G);
      }
      continue;

    case CGO_DRAW_ARRAYS: {
      const cgo::draw::arrays* sp = reinterpret_cast<decltype(sp)>(pc);
      int arrays = sp->arraybits;
      if (reset_colors &&
          arrays & CGO_PICK_COLOR_ARRAY) { // only if picking info is invalid
        int nverts = sp->nverts, v, idx = -1, bnd = -1;
        float* pca = sp->floatdata;

        if (arrays & CGO_VERTEX_ARRAY) {
          pca += nverts * VERTEX_POS_SIZE;
        }
        if (arrays & CGO_NORMAL_ARRAY) {
          pca += nverts * VERTEX_NORMAL_SIZE;
        }
        if (arrays & CGO_COLOR_ARRAY) {
          pca += nverts * VERTEX_COLOR_SIZE;
        }

        auto pickColorValsUC = (uchar*) pca;
        auto pickColorVals = (int*) (pca + nverts * VERTEX_PICKCOLOR_RGBA_SIZE);

        for (v = 0; v < nverts; v++) {
          bnd = pickable ? pickColorVals[v * 2 + 1] : cPickableNoPick;
          idx = pickColorVals[v * 2];
          AssignNewPickColor(
              I, pick, pickColorValsUC + (v * 4), context, idx, bnd);
        }
      }
    } break;

    case CGO_DRAW_BUFFERS_INDEXED:
    case CGO_DRAW_BUFFERS_NOT_INDEXED:
    case CGO_DRAW_TEXTURES:
    case CGO_DRAW_LABELS:
    case CGO_DRAW_SPHERE_BUFFERS:
    case CGO_DRAW_CYLINDER_BUFFERS:
    case CGO_DRAW_CUSTOM: {
      int pickcolors_are_set = true;
      int* pickcolors_are_set_ptr =
          get_pickcolorsset_ptr(op, const_cast<float*>(pc));
      if (!pickcolors_are_set_ptr)
        pickcolors_are_set_ptr = &pickcolors_are_set;

      // TODO remove `pickcolorsset` fields from CGOs
      // This assert can fail during "Roving Detail" demo. However, I still
      // question the need of the `pickcolorsset` fields.
      // assert(reset_colors || *pickcolors_are_set_ptr);

      if (reset_colors ||
          !*pickcolors_are_set_ptr) { // only if picking info is invalid
        int nverts = 0;
        int nvertsperfrag = 1;
        int v, pl;
        int bnd = cPickableNoPick;
        unsigned int idx = 0;
        int srcp;
        float* pca = nullptr;
        int* pickDataSrc;
        uchar* pickColorDestUC = nullptr;
        bool free_pick_color_dest = false;
        int destOffset = 0, bufsizemult = 1;
        size_t pickvbo = 0;
        switch (op) {
        case CGO_DRAW_CUSTOM: {
          const cgo::draw::custom* sp = reinterpret_cast<decltype(sp)>(pc);
          nverts = sp->nverts;
          pickvbo = sp->pickvboid;
          if (!pickvbo)
            continue;
          pca = sp->floatdata;
          nvertsperfrag = sp->vertsperpickinfo;
          bufsizemult = sp->npickbufs;

          pickColorDestUC = new uchar[bufsizemult * nverts * 4];
        } break;
        case CGO_DRAW_BUFFERS_INDEXED: {
          const cgo::draw::buffers_indexed* sp =
              reinterpret_cast<decltype(sp)>(pc);
          nverts = sp->nverts;
          pickvbo = sp->pickvboid;
          pca = sp->floatdata;
        } break;
        case CGO_DRAW_BUFFERS_NOT_INDEXED: {
          const cgo::draw::buffers_not_indexed* sp =
              reinterpret_cast<decltype(sp)>(pc);
          nverts = sp->nverts;
          pickvbo = sp->pickvboid;
          pca = sp->floatdata;
        } break;
        case CGO_DRAW_SPHERE_BUFFERS: {
          const cgo::draw::sphere_buffers* sp =
              reinterpret_cast<decltype(sp)>(pc);
          nverts = sp->num_spheres * VerticesPerSphere();
          nvertsperfrag = VerticesPerSphere();
          pickvbo = sp->pickvboid;
          pca = sp->floatdata;

          pickColorDestUC = new uchar[nverts * 4];
        } break;
        case CGO_DRAW_CYLINDER_BUFFERS: {
          const cgo::draw::cylinder_buffers* sp =
              reinterpret_cast<decltype(sp)>(pc);
          nverts = sp->num_cyl * NumVerticesPerCylinder();
          nvertsperfrag = NumVerticesPerCylinder();
          pickvbo = sp->pickvboid;
          pca = sp->floatdata;
          bufsizemult = 2;

          pickColorDestUC = new uchar[bufsizemult * nverts * 4];
        } break;
        case CGO_DRAW_TEXTURES: {
          const cgo::draw::textures* sp = reinterpret_cast<decltype(sp)>(pc);
          nverts = sp->ntextures * 6;
          pca = sp->floatdata;
        } break;
        case CGO_DRAW_LABELS: {
          const cgo::draw::labels* sp;
          sp = reinterpret_cast<decltype(sp)>(pc);
          nverts = sp->ntextures * 6;
          pca = sp->floatdata;
          pickvbo = sp->pickvboid;
        } break;
        }

        if (pickColorDestUC) {
          free_pick_color_dest = true;
          pickDataSrc = (int*) (pca);
        } else {
          pickColorDestUC = (uchar*) pca;
          pickDataSrc = (int*) (pca + nverts);
        }

        destOffset = R->pick_pass() * sizeof(float) * nverts * bufsizemult;

        if (!pickable) {
          for (int i = 0; i < nverts * bufsizemult; ++i) {
            pick->colorNoPick(pickColorDestUC + 4 * i);
          }
        } else {
          int npickbufs = bufsizemult;
          int ploffsetforbuf = 0;
          if (op == CGO_DRAW_CYLINDER_BUFFERS) {
            // disabled 2016-07-19 TH: code looks almost identical to
            // else branch and CGO_DRAW_CYLINDER_BUFFERS seem to be
            // not used anymore.
            PRINTFB(I->G, FB_CGO, FB_Errors)
            " FIXME: SUPPOSEDLY UNUSED CODE EXECUTED in "
            "CGORenderPicking(!\n" ENDFB(I->G);
          } else {
            if (op == CGO_DRAW_CUSTOM) {
              ploffsetforbuf =
                  sizeof(float) * nverts; // for multiple picking attributes
            }
            for (v = 0, pl = 0; v < nverts; v++, pl += 4) {
              if (v % nvertsperfrag) {
                // if same fragment, same color
                for (int pi = 0; pi < npickbufs; ++pi) {
                  int ploffset = ploffsetforbuf * pi;
                  SetUCColorToPrevN(1, &pickColorDestUC[pl + ploffset]);
                }
                continue;
              }

              int frag = (int) (v / nvertsperfrag);
              for (int pi = 0; pi < npickbufs; ++pi) {
                int ploffset = ploffsetforbuf * pi;
                srcp = 2 * ((npickbufs * frag) + pi);
                idx = pickDataSrc[srcp];
                bnd = pickDataSrc[srcp + 1];

                AssignNewPickColor(I, pick, &pickColorDestUC[pl + ploffset],
                    context, idx, bnd);
              }
            }
          }
        }

        if (pickvbo) {
          // reload entire vbo
          VertexBuffer* vbo =
              I->G->ShaderMgr->getGPUBuffer<VertexBuffer>(pickvbo);
          vbo->bufferReplaceData(destOffset,
              sizeof(float) * nverts * bufsizemult, pickColorDestUC);
          (*pickcolors_are_set_ptr) = true;
        }

        if (free_pick_color_dest) {
          delete[] pickColorDestUC;
          pickColorDestUC = nullptr;
          free_pick_color_dest = false;
        }
      }
    } break;
    }

    CGO_gl[op](R, &pc);
  }

  R->isPicking = false;
}

void CGORenderGL(CGO* I, const float* color, CSetting* set1, CSetting* set2,
    RenderInfo* info, Rep* rep)
/* this should be as fast as you can make it...

 * the ASM loop is about 2X long as raw looped GL calls,

 * but hopefully superscaler processors won't care */
{
  PyMOLGlobals* G = I->G;

  const float zee[] = {0.f, 0.f, 1.f};
  const float color_tmp[] = {1.f, 1.f, 1.f};

  if (I->render_alpha) {
    // for now, the render_alpha_only flag calls CGOSetZVector/CGORenderGLAlpha
    auto modMatrix = SceneGetModelViewMatrixPtr(G);
    CGOSetZVector(I, modMatrix[2], modMatrix[6], modMatrix[10]);
    CGORenderAlpha(I, info, 1);
    if (I->render_alpha == 1) // right now, render_alpha 1: renders alpha only,
                              // 2: renders both alpha and rest
      return;
  }

  if (!G->ValidContext) {
    return;
  }

  if (!I->c) {
    return;
  }

  {
    CCGORenderer* R = G->CGORenderer;
    R->info = info;
    R->use_shader = I->use_shader;
    R->debug = I->debug;
    R->sphere_quality = I->sphere_quality;
    R->rep = rep;
    R->color = color;
    R->alpha = 1.0F - SettingGet_f(G, set1, set2, cSetting_cgo_transparency);
    R->set1 = set1;
    R->set2 = set2;
    // normals should be initialized to the view vector
    // (changed BB 9/14 from SceneResetNormalUseShader(), to
    // CScene->LinesNormal, which was arbitrary, I believe)
    SceneResetNormalToViewVector(I->G, I->use_shader);

    if (!color) {
      color = color_tmp;
    }

    {
      auto shaderPrg = I->G->ShaderMgr->Get_Current_Shader();
      if (shaderPrg && I->use_shader) {
        shaderPrg->SetAttrib4fLocation(
            "a_Color", color[0], color[1], color[2], R->alpha);
      }
#ifndef PURE_OPENGL_ES_2
      else {
        glColor4f(color[0], color[1], color[2], R->alpha);
      }
#endif
    }

#ifndef PURE_OPENGL_ES_2
    const float width_scale =
        (info && info->width_scale_flag) ? info->width_scale : 1.f;
    glLineWidth(
        SettingGet_f(G, set1, set2, cSetting_cgo_line_width) * width_scale);
    glPointSize(
        SettingGet_f(G, set1, set2, cSetting_cgo_dot_width) * width_scale);
#endif

    if (!(info && info->alpha_cgo)) {
      // Regular CGO dispatch table rendering
      for (auto it = I->begin(); !it.is_stop(); ++it) {
        const auto op = it.op_code();
        assert(op < CGO_sz_size());
        CGO_OP_DATA_CONST float* const pc = it.data();

        CGO_gl[op](R, &pc);
      }
      return;
    }

    /* we're sorting transparent triangles globally */
    {
      {
        int mode = -1;
        int vc = 0;
        // triangle normals
        const float *n0 = nullptr, *n1 = nullptr, *n2 = nullptr;
        // triangle vertices
        const float *v0 = nullptr, *v1 = nullptr, *v2 = nullptr;
        // triangle colors
        const float *c0 = color, *c1 = nullptr, *c2 = nullptr;

        for (auto it = I->begin(); !it.is_stop(); ++it) {
          const auto op = it.op_code();
          assert(op < CGO_sz_size());
          CGO_OP_DATA_CONST float* const pc = it.data();

          if ((R->alpha != 1.f)) {
            switch (op) { /* transparency */
            case CGO_BEGIN:
              mode = CGO_get_int(pc);
              CGO_gl_begin(R, &pc);
              vc = 0;
              n0 = zee;
              break;
            case CGO_END:
              CGO_gl_end(R, &pc);
              mode = -1;
              break;
            case CGO_NORMAL:
              switch (mode) {
              case GL_TRIANGLES:
              case GL_TRIANGLE_STRIP:
              case GL_TRIANGLE_FAN:
                n0 = pc;
                break;
              default:
                CGO_gl_normal(R, &pc);
              }
              break;
            case CGO_COLOR:
              c0 = pc;
              CGO_gl_color(R, &pc);
              break;
            case CGO_TRIANGLE:
              CGOAlphaTriangle(info->alpha_cgo, pc, pc + 3, pc + 6, pc + 9,
                  pc + 12, pc + 15, pc + 18, pc + 21, pc + 24, R->alpha,
                  R->alpha, R->alpha, false);
              break;
            case CGO_DRAW_ARRAYS: {
              const cgo::draw::arrays* sp = reinterpret_cast<decltype(sp)>(pc);
              int mode = sp->mode, arrays = sp->arraybits, nverts = sp->nverts;
              float *vertexVals = 0, *nxtVals = 0, *colorVals = 0, *normalVals;
              float *vertexVals_tmp = 0, *colorVals_tmp = 0,
                    *normalVals_tmp = 0;
              int step;
              short nxtn = 3;
              nxtVals = vertexVals = vertexVals_tmp = sp->floatdata;
              if (arrays & CGO_NORMAL_ARRAY) {
                nxtVals = normalVals = normalVals_tmp =
                    vertexVals + (nxtn * nverts);
              }
              if (arrays & CGO_COLOR_ARRAY) {
                nxtVals = colorVals = colorVals_tmp = nxtVals + (nxtn * nverts);
                nxtn = 4;
              }
              switch (mode) {
              case GL_TRIANGLES: {
                for (step = 0; step < nverts; step += 3) {
                  if (colorVals_tmp) {
                    c0 = colorVals_tmp;
                    c1 = colorVals_tmp + 4;
                    c2 = colorVals_tmp + 8;
                  } else {
                    c1 = c2 = c0;
                  }
                  if (normalVals_tmp) {
                    n0 = normalVals_tmp;
                    n1 = normalVals_tmp + 3;
                    n2 = normalVals_tmp + 6;
                  } else {
                    n1 = n2 = n0;
                  }
                  CGOAlphaTriangle(info->alpha_cgo, vertexVals_tmp,
                      vertexVals_tmp + 3, vertexVals_tmp + 6, n0, n1, n2, c0,
                      c1, c2, R->alpha, R->alpha, R->alpha, false);
                  vertexVals_tmp += 9;
                  if (normalVals_tmp) {
                    normalVals_tmp += 9;
                  }
                  if (colorVals_tmp) {
                    colorVals_tmp += 12;
                  }
                }
              } break;
              case GL_TRIANGLE_STRIP: {
                if (colorVals_tmp) {
                  c1 = colorVals_tmp;
                  c2 = colorVals_tmp + 4;
                  colorVals_tmp += 8;
                } else {
                  c1 = c2 = c0;
                }
                if (normalVals_tmp) {
                  n1 = normalVals_tmp;
                  n2 = normalVals_tmp + 3;
                  normalVals_tmp += 6;
                } else {
                  n1 = n2 = n0;
                }
                vertexVals_tmp += 6;
                for (step = 2; step < nverts; step++) {
                  if (colorVals_tmp) {
                    c0 = c1;
                    c1 = c2;
                    c2 = colorVals_tmp;
                  }
                  if (normalVals_tmp) {
                    n0 = n1;
                    n1 = n2;
                    n2 = normalVals_tmp;
                  }
                  CGOAlphaTriangle(info->alpha_cgo, vertexVals_tmp - 6,
                      vertexVals_tmp - 3, vertexVals_tmp, n0, n1, n2, c0, c1,
                      c2, R->alpha, R->alpha, R->alpha, false);
                  vertexVals_tmp += 3;
                  if (normalVals_tmp) {
                    normalVals_tmp += 3;
                  }
                  if (colorVals_tmp) {
                    colorVals_tmp += 4;
                  }
                }
              } break;
              case GL_TRIANGLE_FAN: {
                float* firstVertex = vertexVals_tmp;
                if (colorVals_tmp) {
                  c0 = colorVals_tmp;
                  c2 = colorVals_tmp + 4;
                  colorVals_tmp += 8;
                } else {
                  c1 = c2 = c0;
                }
                if (normalVals_tmp) {
                  n0 = normalVals_tmp;
                  n2 = normalVals_tmp + 3;
                  normalVals_tmp += 6;
                }
                vertexVals_tmp += 6;
                for (step = 2; step < nverts; step++) {
                  if (colorVals_tmp) {
                    c1 = c2;
                    c2 = colorVals_tmp;
                  }
                  if (normalVals_tmp) {
                    n1 = n2;
                    n2 = normalVals_tmp;
                  }
                  CGOAlphaTriangle(info->alpha_cgo, firstVertex,
                      vertexVals_tmp - 3, vertexVals_tmp, n0, n1, n2, c0, c1,
                      c2, R->alpha, R->alpha, R->alpha, false);
                  vertexVals_tmp += 3;
                  if (normalVals_tmp) {
                    normalVals_tmp += 3;
                  }
                  if (colorVals_tmp) {
                    colorVals_tmp += 4;
                  }
                }
              } break;
              }
            } break;
            case CGO_VERTEX:
              v0 = pc;
              switch (mode) {
              case GL_TRIANGLES:
                if (3 * ((vc + 1) / 3) == vc + 1) {
                  CGOAlphaTriangle(info->alpha_cgo, v0, v1, v2, n0, n1, n2, c0,
                      c1, c2, R->alpha, R->alpha, R->alpha, true);
                }
                v2 = v1;
                c2 = c1;
                n2 = n1;
                v1 = v0;
                c1 = c0;
                n1 = n0;
                vc++;
                break;
              case GL_TRIANGLE_STRIP:
                if (vc > 1) {
                  CGOAlphaTriangle(info->alpha_cgo, v0, v1, v2, n0, n1, n2, c0,
                      c1, c2, R->alpha, R->alpha, R->alpha, !(vc & 0x1));
                }
                v2 = v1;
                c2 = c1;
                n2 = n1;
                v1 = v0;
                c1 = c0;
                n1 = n0;
                vc++;
                break;
              case GL_TRIANGLE_FAN:
                if (vc > 1) {
                  CGOAlphaTriangle(info->alpha_cgo, v0, v1, v2, n0, n1, n2, c0,
                      c1, c2, R->alpha, R->alpha, R->alpha, false);
                } else if (!vc) {
                  n2 = n0;
                  v2 = v0;
                  c2 = c0;
                }
                v1 = v0;
                c1 = c0;
                n1 = n0;
                vc++;
                break;
              default:
                CGO_gl_vertex(R, &pc);
                break;
              }
              break;
            default:
              CGO_gl[op](R, &pc);
              break;
            }
          } else { /* opaque */
            switch (op) {
            case CGO_COLOR:
              /* Since CGO operations are done in sequence, alpha could happen
                 after color is set.  In this case, we still need to keep track
                 of the color in case there is a transparent object */
              c0 = pc;
              break;
            default:
              break;
            }
            CGO_gl[op](R, &pc);
          }
        }
      }
    }
  }
}

void CGORenderGLAlpha(CGO* I, RenderInfo* info, bool calcDepth)
{
  PyMOLGlobals* G = I->G;
  if (G->ValidContext && I->c) {
    int mode = GL_TRIANGLES;
    if (I->debug) {
      mode = CGOConvertDebugMode(I->debug, GL_TRIANGLES);
    }
#ifndef PURE_OPENGL_ES_2
    // not sure why shader is set, but disable it for now,
    // since we are doing immediate mode rendering for global transparency
    G->ShaderMgr->Disable_Current_Shader();
#endif
    /* 1. transform and measure range (if not already known)
       2. bin into linked lists based on Z-centers
       3. render by layer */

    if (I->z_flag) {
      if (!I->i_start) {
        I->i_size = 256;
        I->i_start = pymol::calloc<int>(I->i_size);
      } else {
        UtilZeroMem(I->i_start, sizeof(int) * I->i_size);
      }
      {
        const int i_size = I->i_size;
        const float* base = I->op;
        int* start = I->i_start;
        int delta = 1, ntris = 0;
        /* bin the triangles */
        if (calcDepth) {
          for (auto it = I->begin(); !it.is_stop(); ++it) {
            if (it.op_code() == CGO_ALPHA_TRIANGLE) {
              float* const pc = it.data();
              const float z = dot_product3f(pc + 1, I->z_vector);
              if (z > I->z_max)
                I->z_max = z;
              if (z < I->z_min)
                I->z_min = z;
              pc[4] = z;
              ntris++;
            }
          }
        }

        const float range_factor = (0.9999F * i_size) / (I->z_max - I->z_min);

        for (auto it = I->begin(); !it.is_stop(); ++it) {
          if (it.op_code() == CGO_ALPHA_TRIANGLE) {
            float* const pc = it.data();
            assert(base < pc && pc < I->op + I->c);
            auto i =
                std::clamp<int>((pc[4] - I->z_min) * range_factor, 0, i_size);
            CGO_put_int(pc, start[i]);
            start[i] = (pc - base); /* NOTE: will always be > 0 since we have
                                       CGO_read_int'd */
          }
        }

        // for single-layer transparency, render front-to-back
        if (SettingGetGlobal_i(G, cSetting_transparency_mode) == 2) {
          delta = -1;
          start += (i_size - 1);
        }

        /* now render by bin */
#ifndef PURE_OPENGL_ES_2
        glBegin(mode);
        for (int i = 0; i < i_size; i++) {
          int ii = *start;
          start += delta;
          while (ii) {
            const float* pc = base + ii;
            glColor4fv(pc + 23);
            glNormal3fv(pc + 14);
            glVertex3fv(pc + 5);
            glColor4fv(pc + 27);
            glNormal3fv(pc + 17);
            glVertex3fv(pc + 8);
            glColor4fv(pc + 31);
            glNormal3fv(pc + 20);
            glVertex3fv(pc + 11);

            ii = CGO_get_int(pc);
          }
        }
        glEnd();
#endif
      }
    } else {
#ifndef PURE_OPENGL_ES_2
      glBegin(mode);
      for (auto it = I->begin(); !it.is_stop(); ++it) {
        if (it.op_code() == CGO_ALPHA_TRIANGLE) {
          float* const pc = it.data();
          glColor4fv(pc + 23);
          glNormal3fv(pc + 14);
          glVertex3fv(pc + 5);
          glColor4fv(pc + 27);
          glNormal3fv(pc + 17);
          glVertex3fv(pc + 8);
          glColor4fv(pc + 31);
          glNormal3fv(pc + 20);
          glVertex3fv(pc + 11);
        }
      }
      glEnd();
#endif
    }
  }
}
