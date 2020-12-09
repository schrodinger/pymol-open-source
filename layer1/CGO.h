
/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* copyright by Schrodinger, LLC
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
#ifndef _H_CGO
#define _H_CGO

#include"Base.h"
#include"Ray.h"
#include"Setting.h"
#include"os_gl.h"
#include"os_gl_cgo.h"
#include"Rep.h"
#include"MemoryDebug.h"
#include <vector>
#include <unordered_map>
#include <typeinfo>
#include <type_traits>
#include <memory>
#include "GenericBuffer.h"
#include <set>
#include <glm/vec3.hpp>
#include <glm/vec4.hpp>

/* Compiled Graphics Library for simple graphics objects
   in floating point three-space, with the goal of achieving
   quick and easy rendering in multiple environments without the
   headaches of OpenGL arrays.

*/

#define CGO_read_int(p) (*((int*)((p)++)))
#define CGO_get_int(p) (*((int*)(p)))
#define CGO_get_uint(p) (uint)(*((uint*)(p)))
#define CGO_write_int(p,i) ((*((int*)(p++)))=(i))
#define CGO_write_uint(p,i) ((*((uint*)(p++)))=(i))
#define CGO_put_int(p,i) ((*((int*)(p)))=(i))
#define CGO_put_uint(p,i) ((*((uint*)(p)))=(i))

inline uchar CLIP_COLOR_VALUE(float cv){ return ((cv>1.f) ? 255 :  (cv < 0.f) ? 0 : pymol_roundf(cv * 255) ); }
inline float CONVERT_COLOR_VALUE(unsigned char cv) { return ((cv>255) ? 1.f :  (cv < 0) ? 0.f : (cv / 255.f) ); }

// normal values are mapped between { -1, 1 } to { 0, 255 }, 
// the values are mapped: { -1/128, -.5/192, 0./255, .5/63, 1./127 }
inline uchar CLIP_NORMAL_VALUE(float cv){ return ((cv>1.f) ? 127 :
                                                  (cv < -1.f) ? 128 : 
                                                  pymol_roundf(((cv + 1.f)/2.f) * 255) - 128 ); }

/* Supported functions:
 * stop
 * null
 * begin
     GL_POINTS,
     GL_LINES, GL_LINE_LOOP, GL_LINE_STRIP,
     GL_TRIANGLES, GL_TRIANGLE_STRIP, GL_TRIANGLE_FAN
 * end
 * vertex
 * normal
 * color
 * sphere   * currently for ray-tracing only
 * triangle * currently for ray-tracing only
 * cylinder * currently for ray-tracing only
 * linewidth
 * primwidth * ray-tracing
 */

/* instructions and data segment sizes */

#define CGO_STOP                 0x00
#define CGO_STOP_SZ              0
#define CGO_NULL                 0x01
#define CGO_NULL_SZ              0
#define CGO_BEGIN                0x02
#define CGO_BEGIN_SZ             1
#define CGO_END                  0x03
#define CGO_END_SZ               0
#define CGO_VERTEX               0x04
#define CGO_VERTEX_SZ            3
#define CGO_NORMAL               0x05
#define CGO_NORMAL_SZ            3
#define CGO_COLOR                0x06
#define CGO_COLOR_SZ             3
#define CGO_SPHERE               0x07
#define CGO_SPHERE_SZ            4
#define CGO_TRIANGLE             0x08
#define CGO_TRIANGLE_SZ          27
#define CGO_CYLINDER             0x09
#define CGO_CYLINDER_SZ          13
// CGO_LINEWIDTH, CGO_WIDTHSCALE, CGO_DOTWIDTH work only on ray tracing
#define CGO_LINEWIDTH            0x0A
#define CGO_LINEWIDTH_SZ         1
#define CGO_WIDTHSCALE           0x0B
#define CGO_WIDTHSCALE_SZ        1

#define CGO_ENABLE               0x0C
#define CGO_ENABLE_SZ            1
#define CGO_DISABLE              0x0D
#define CGO_DISABLE_SZ           1
#define CGO_SAUSAGE              0x0E
#define CGO_SAUSAGE_SZ           13
#define CGO_CUSTOM_CYLINDER      0x0F
#define CGO_CUSTOM_CYLINDER_SZ   15

#define CGO_DOTWIDTH             0x10
#define CGO_DOTWIDTH_SZ          1

#define CGO_ALPHA_TRIANGLE       0x11
#define CGO_ALPHA_TRIANGLE_SZ    35
#define CGO_ELLIPSOID            0x12
#define CGO_ELLIPSOID_SZ         13
#define CGO_FONT                 0x13
#define CGO_FONT_SZ              3      /*  size, face, style */
#define CGO_FONT_SCALE           0x14
#define CGO_FONT_SCALE_SZ        2
#define CGO_FONT_VERTEX          0x15
#define CGO_FONT_VERTEX_SZ       3      /*  principle axes (zeros -> use camera x y z */

// CGO_FONT_AXES not used
#define CGO_FONT_AXES            0x16
#define CGO_FONT_AXES_SZ         9      /*  principle axes (zeros -> use camera x y z */

#define CGO_CHAR                 0x17
#define CGO_CHAR_SZ              1
#define CGO_INDENT               0x18
#define CGO_INDENT_SZ            2
#define CGO_ALPHA                0x19
#define CGO_ALPHA_SZ             1
#define CGO_QUADRIC              0x1A
#define CGO_QUADRIC_SZ           14
#define CGO_CONE                 0x1B
#define CGO_CONE_SZ              16
#define CGO_RESET_NORMAL         0x1E
#define CGO_RESET_NORMAL_SZ      1
#define CGO_PICK_COLOR           0x1F
#define CGO_PICK_COLOR_SZ        2


/* CGO_DRAW_ARRAYS : operation that calls glDrawArrays with all arrays in memory
   (i.e., stored in the CGO array). There can be up to 4 arrays (vertex, normal, color,
   and pick color array, where each are stored using GL_FLOAT array (except for the pick
   color array, which is GL_UNSIGNED_BYTE, because the 2nd 2/3rds of the array is used to
   store atom/bond picking information. Also, the color array is stored using 4 floats,
   since glColorPointer() requires having 4 as an argument for OpenGLES.
   - mode : GL Mode that is used
   - arrays : which arrays that are used (bitmask on CGO_<type>_ARRAY, where <type>
              can be VERTEX, NORMAL, COLOR, or PICK_COLOR)
   - narrays : number of arrays specified in arrays
   - nverts : number of total vertices specified in arrays
   - each GL_FLOAT array specified in order
*/
#define CGO_DRAW_ARRAYS          0x1C
/* CGO_DRAW_BUFFERS_INDEXED : operation that uses glDrawArrays with VBOs.  It also
   has a dynamic length since it holds the picking information (atom/bond) for each vertex,
   as well as space for picking to use when drawing colors (since colors are dynamically generated
   for each atom/bond value).
   - mode : GL Mode that is used
   - arrays : which arrays that are used (bitmask on CGO_<type>_ARRAY, where <type>
              can be VERTEX, NORMAL, COLOR, or PICK_COLOR)
   - narrays : number of arrays specified in arrays
   - nverts : number of total vertices specified in arrays
   - bufs[5] : each VBO id in order (if used, VERTEX, NORMAL, COLOR, PICK_COLOR), plus the
      vertex index array that specifies all vertices.
 */
#define CGO_DRAW_BUFFERS_INDEXED         0x21

/* CGO_BOUNDING_BOX : operation that allows the extent to be expanded.  Since the geometry
   data is not kept in the CGO object for VBO objects (only on the card), this object allows
   the bounding box of the object to be saved in the CGO.  This is used in the CGOGetExtent(),
   typically when the view is being automatically set, but can be used for other things.
 */
#define CGO_BOUNDING_BOX         0x22
#define CGO_BOUNDING_BOX_SZ      6

#define CGO_DRAW_BUFFERS_NOT_INDEXED         0x23

#define CGO_SPECIAL            0x24
#define CGO_SPECIAL_SZ         1

#define CGO_DRAW_CYLINDER_BUFFERS       0x25

#define CGO_SHADER_CYLINDER             0x26
#define CGO_SHADER_CYLINDER_SZ          8

#define CGO_SHADER_CYLINDER_WITH_2ND_COLOR      0x27
#define CGO_SHADER_CYLINDER_WITH_2ND_COLOR_SZ    13

#define CGO_DRAW_SPHERE_BUFFERS      0x28

#define CGO_ACCESSIBILITY      0x29
#define CGO_ACCESSIBILITY_SZ    1

#define CGO_DRAW_TEXTURE      0x2A
#define CGO_DRAW_TEXTURE_SZ    13

#define CGO_DRAW_TEXTURES      0x2B

#define CGO_DRAW_SCREEN_TEXTURES_AND_POLYGONS      0x2C

#define CGO_TEX_COORD                0x2D
#define CGO_TEX_COORD_SZ             2


#define CGO_DRAW_LABEL      0x2E
#define CGO_DRAW_LABEL_SZ    20

#define CGO_DRAW_LABELS      0x2F

#define CGO_DRAW_CONNECTOR       0x30
#define CGO_DRAW_CONNECTOR_SZ    25

#define CGO_DRAW_CONNECTORS      0x31

#define CGO_DRAW_TRILINES        0x32
#define CGO_DRAW_TRILINES_SZ     2

#define CGO_UNIFORM3F            0x33
#define CGO_UNIFORM3F_SZ         4
#define CGO_UNIFORM3F_HEADER     4

#define CGO_SPECIAL_WITH_ARG    0x34
#define CGO_SPECIAL_WITH_ARG_SZ    2

#define CGO_LINE                0x35
#define CGO_LINE_SZ                6
#define CGO_SPLITLINE           0x36
#define CGO_SPLITLINE_SZ           9

#define CGO_DRAW_CUSTOM         0x37

#define CGO_VERTEX_ATTRIB_3F    0x38
#define CGO_VERTEX_ATTRIB_4UB   0x39
#define CGO_VERTEX_ATTRIB_1F    0x3A

#define CGO_MASK_ATTRIBUTE_IF_PICKING   0x3B
#define CGO_BIND_VBO_FOR_PICKING        0x3C

#define CGO_VERTEX_BEGIN_LINE_STRIP     0x3D
#define CGO_VERTEX_BEGIN_LINE_STRIP_SZ     3

#define CGO_INTERPOLATED     0x3E
#define CGO_INTERPOLATED_SZ     1

#define CGO_VERTEX_CROSS     0x3F
#define CGO_VERTEX_CROSS_SZ     3

#define CGO_VERTEX_ATTRIB_4UB_IF_PICKING   0x40

#define CGO_CUSTOM_CYLINDER_ALPHA      0x41

#define CGO_LIGHTING             0x0B50

#define CGO_VERTEX_ARRAY         0x01
#define CGO_NORMAL_ARRAY         0x02
#define CGO_COLOR_ARRAY          0x04
#define CGO_PICK_COLOR_ARRAY     0x08
#define CGO_ACCESSIBILITY_ARRAY  0x10
#define CGO_TEX_COORD_ARRAY      0x20

extern int CGO_sz[];
size_t CGO_sz_size();

// CGOs are floating point arrays so we need to work with sizes in terms of floats
template <typename T>
static size_t fsizeof() { return sizeof(T) / sizeof(float); }

struct ObjectGadgetRamp;

class CGO;

// These are only the optimized operations
namespace cgo {
  namespace draw {
    struct begin {
      static const int op_code = CGO_BEGIN;
      int mode;
      begin(int mode) : mode(mode) {}
    };

    struct enable {
      static const int op_code = CGO_ENABLE;
      int mode;
      enable(int mode) : mode(mode) {}
    };

    struct disable {
      static const int op_code = CGO_DISABLE;
      int mode;
      disable(int mode) : mode(mode) {}
    };

    struct trilines {
      static const int op_code = CGO_DRAW_TRILINES;
      unsigned nverts;
      unsigned buffer; // GLuint
      trilines(unsigned nverts, unsigned buffer) : nverts(nverts), buffer(buffer) {}
    };

    struct op_with_data {
      float * floatdata { nullptr };
      void set_data(float * data) { floatdata = data; };
      float * get_data() { return floatdata; };
      const float * get_data() const { return floatdata; };
      virtual int get_data_length() const = 0;
      virtual ~op_with_data() {}
    };

    struct op_with_draw_buffers {
    };

    struct arrays : op_with_data {
      static const int op_code = CGO_DRAW_ARRAYS;
      arrays(int _mode, short _arrays, int _nverts) :
      mode(_mode), arraybits(_arrays), narrays(0), nverts(_nverts) {
        short bit;
        for (bit = 0; bit < 4; bit++){
          if ((1 << bit) & arraybits){
            narrays+=3;
          }
        }
        if (arraybits & CGO_ACCESSIBILITY_ARRAY) narrays++;
        if (arraybits & CGO_COLOR_ARRAY) narrays++;
      };
      int mode;
      int arraybits;
      int narrays;
      int nverts;
      int get_data_length() const { return narrays * nverts; };
    };

    struct buffers_indexed : op_with_data, op_with_draw_buffers {
      static const int op_code = CGO_DRAW_BUFFERS_INDEXED;
      buffers_indexed(int _mode, short _arrays, int _nindices,
                      int _nverts, size_t _vboid, size_t _iboid, int _n_data, size_t _pickvboid = 0) :
        mode(_mode), arraybits(_arrays), narrays(0), nindices(_nindices),
        nverts(_nverts), vboid(_vboid), iboid(_iboid)
        , pickvboid(_pickvboid)
        , pickcolorsset(0)
        , n_data(_n_data)
      {
        short bit;
        for (bit = 0; bit < 4; bit++){
          if ((1 << bit) & arraybits){
            narrays++;
          }
        }
        if (arraybits & CGO_ACCESSIBILITY_ARRAY) narrays++;
        if (arraybits & CGO_COLOR_ARRAY) narrays++;
      }
      int mode;
      int arraybits;
      int narrays;
      int nindices;
      int nverts;
      size_t vboid;
      size_t iboid;
      size_t pickvboid;
      int pickcolorsset;
      int n_data;
      int get_data_length() const { return nverts * 3 + n_data; };
    };

    struct sphere_buffers : op_with_data, op_with_draw_buffers {
      static const int op_code = CGO_DRAW_SPHERE_BUFFERS;
      sphere_buffers(int _num_spheres, int _ub_flags,
                     size_t _vboid, size_t _pickvboid) :
        num_spheres(_num_spheres), ub_flags(_ub_flags),
        vboid(_vboid), pickvboid(_pickvboid), pickcolorsset(0) {}
      int num_spheres;
      int ub_flags;
      size_t vboid;
      size_t pickvboid;
      int pickcolorsset;
      int get_data_length() const { return num_spheres * 2; };
    };

    struct cylinder_buffers : op_with_data, op_with_draw_buffers {
      static const int op_code = CGO_DRAW_CYLINDER_BUFFERS;
      cylinder_buffers(int _num_cyl, int _alpha, size_t _vboid,
                       size_t _iboid, size_t _pickvboid) :
        num_cyl(_num_cyl), alpha(_alpha), vboid(_vboid),
        iboid(_iboid), pickvboid(_pickvboid), pickcolorsset(0) {
      }
      int num_cyl;
      int alpha;
      size_t vboid;
      size_t iboid;
      size_t pickvboid;
      int pickcolorsset;
      int get_data_length() const { return num_cyl * 2 * 2; };
    };

    struct textures : op_with_data, op_with_draw_buffers {
      static const int op_code = CGO_DRAW_TEXTURES;
      textures(int _ntextures, size_t _vboid) : ntextures(_ntextures), vboid(_vboid) { }
      int ntextures;
      size_t vboid;
      int get_data_length() const { return ntextures * 18; };
    };

    struct screen_textures : op_with_draw_buffers {
      static const int op_code = CGO_DRAW_SCREEN_TEXTURES_AND_POLYGONS;
      screen_textures(int _nverts, size_t _vboid) : nverts(_nverts),
                                                        vboid(_vboid) {};
      int nverts;
      size_t vboid;
    };

    struct buffers_not_indexed : op_with_data, op_with_draw_buffers {
      static const int op_code = CGO_DRAW_BUFFERS_NOT_INDEXED;
      buffers_not_indexed(int _mode, int _arrays, int _nverts,
                          size_t _vboid, size_t _pickvboid = 0) :
        mode(_mode), arraybits(_arrays), narrays(0), nverts(_nverts),
        vboid(_vboid), pickvboid(_pickvboid), pickcolorsset(0) {
        for (short bit = 0; bit < 4; bit++){
          if ((1 << bit) & arraybits){
            narrays++;
          }
        }
        if (arraybits & CGO_ACCESSIBILITY_ARRAY) narrays++;
        if (arraybits & CGO_COLOR_ARRAY) narrays++;
      };
      int mode;
      int arraybits;
      int narrays;
      int nverts;
      size_t vboid;
      size_t pickvboid;
      int pickcolorsset;
      int get_data_length() const { return nverts * 3; };
    };

    struct label {
      using vec3 = glm::vec3;
      using vec4 = glm::vec4;
      static const int op_code = CGO_DRAW_LABEL;
      label(const vec3& _world_pos, const vec3& _screen_world_offset, const vec3& _screen_min,
            const vec3& _screen_max, const vec4& _text_extent, float _relative_mode,
            const vec3& _target_pos) :
        world_pos(_world_pos), screen_world_offset(_screen_world_offset), screen_min(_screen_min),
        screen_max(_screen_max), text_extent(_text_extent), relative_mode(_relative_mode),
        target_pos(_target_pos) {}
      vec3 world_pos;
      vec3 screen_world_offset;
      vec3 screen_min;
      vec3 screen_max;
      vec4 text_extent;
      float relative_mode;
      vec3 target_pos;
    };

    struct labels : op_with_data, op_with_draw_buffers {
      static const int op_code = CGO_DRAW_LABELS;
      labels(int _ntextures, size_t _vboid, size_t _pickvboid) :
        vboid(_vboid)
        , pickvboid(_pickvboid)
        , ntextures(_ntextures)
        , pickcolorsset(0)
      {}
      size_t vboid;
      size_t pickvboid;
      int ntextures;
      int pickcolorsset;
      int get_data_length() const { return ntextures * 18; };
    };

    struct connectors : op_with_draw_buffers {
      static const int op_code = CGO_DRAW_CONNECTORS;
      connectors(int _nconnectors, size_t _vboid) :
        nconnectors(_nconnectors), vboid(_vboid) {};
      int nconnectors;
      size_t vboid;
    };

    struct line {
      static const int op_code = CGO_LINE;
      line(const float *v1, const float *v2) {
        copy3f(v1, vertex1);
        copy3f(v2, vertex2);
      };
      float vertex1[3], vertex2[3];
    };
    struct splitline {
      static const int op_code = CGO_SPLITLINE;

      static const unsigned char interpolation = 0x01; // ramp/interpolation
      static const unsigned char no_split_for_pick = 0x02;
      static const unsigned char equal_colors = 0x04;
      splitline(const float *v1, const float *v2, const float *color2arg, 
                unsigned int index_2nd, int bond_2nd,
                bool isRamped, bool same_pick, bool eq_colors) :
      index(index_2nd), bond(bond_2nd) {
          copy3f(v1, vertex1);
          copy3f(v2, vertex2);
          color2[0] = CLIP_COLOR_VALUE(color2arg[0]);
          color2[1] = CLIP_COLOR_VALUE(color2arg[1]);
          color2[2] = CLIP_COLOR_VALUE(color2arg[2]);
          flags = (isRamped ? cgo::draw::splitline::interpolation : 0x00) |
            (same_pick ? cgo::draw::splitline::no_split_for_pick : 0x00) |
            (eq_colors ? cgo::draw::splitline::equal_colors : 0x00);
      };
      float vertex1[3], vertex2[3];
      unsigned char color2[3];
      unsigned char flags; // interpolation | no_split_for_pick | equal_colors
      unsigned int index;
      int bond;
    };

    struct mask_attribute_if_picking : op_with_draw_buffers {
      static const int op_code = CGO_MASK_ATTRIBUTE_IF_PICKING;
      mask_attribute_if_picking(int _attr_lookup_idx, size_t _vboid) {
        attr_lookup_idx = _attr_lookup_idx;
        vboid = _vboid;
      };
      int attr_lookup_idx;
      size_t vboid;
    };

    struct bind_vbo_for_picking : op_with_draw_buffers {
      static const int op_code = CGO_BIND_VBO_FOR_PICKING;
      bind_vbo_for_picking(size_t _vboid, int _which_attr_idx, int _npickattrs) {
        vboid = _vboid;
        which_attr_idx = _which_attr_idx;
        npickattrs = _npickattrs;
      };
      size_t vboid;
      int which_attr_idx;
      int npickattrs;
    };

    struct custom : op_with_data, op_with_draw_buffers {
      static const int op_code = CGO_DRAW_CUSTOM;
      custom(int _mode, int _nverts, size_t _vboid, size_t _pickvboid,
          int _vertsperpickinfo = 1, int _npickbufs = 1, size_t _iboid = 0,
          int _nindices = 0)
          : mode(_mode), nverts(_nverts), vboid(_vboid), pickvboid(_pickvboid),
            pickcolorsset(0), iboid(_iboid), nindices(_nindices),
            npickbufs(_npickbufs), vertsperpickinfo(_vertsperpickinfo){};
      int mode;
      int nverts;
      size_t vboid;
      size_t pickvboid;
      int pickcolorsset;
      size_t iboid;
      int nindices;
      int npickbufs;
      int vertsperpickinfo;
      int get_data_length() const { return npickbufs*nverts*2/vertsperpickinfo; };
    };

    struct vertex_attribute_3f {
      static const int op_code = CGO_VERTEX_ATTRIB_3F;
      vertex_attribute_3f(int _attr_lookup_idx, void *v) {
        attr_lookup_idx = _attr_lookup_idx;
        copy3f((float*)v, values);
      };
      int attr_lookup_idx;
      float values[3];
    };

    struct vertex_attribute_4ub {
      static const int op_code = CGO_VERTEX_ATTRIB_4UB;
      vertex_attribute_4ub(int _attr_lookup_idx, void *_ubdata) {
        attr_lookup_idx = _attr_lookup_idx;
        memcpy(ubdata, _ubdata, 4);
      };
      int attr_lookup_idx;
      unsigned char ubdata[4];
    };

    struct vertex_attribute_4ub_if_picking {
      static const int op_code = CGO_VERTEX_ATTRIB_4UB_IF_PICKING;
      vertex_attribute_4ub_if_picking(int _attr_lookup_idx, void *_ubdata) {
        attr_lookup_idx = _attr_lookup_idx;
        memcpy(ubdata, _ubdata, 4);
      };
      int attr_lookup_idx;
      unsigned char ubdata[4];
    };

    struct vertex_attribute_1f {
      static const int op_code = CGO_VERTEX_ATTRIB_1F;
      vertex_attribute_1f(int _attr_lookup_idx, float v) {
        attr_lookup_idx = _attr_lookup_idx;
        value = v;
      };
      int attr_lookup_idx;
      float value;
    };

    struct shadercylinder {
      static const int op_code = CGO_SHADER_CYLINDER;
    shadercylinder(const float *_origin, const float *_axis, const float _tube_size, int _cap) :
      tube_size(_tube_size), cap(_cap) {
          copy3f(_origin, origin);
          copy3f(_axis, axis);
      };
      float origin[3], axis[3], tube_size;
      int cap;
    };

    struct shadercylinder2ndcolor {
      static const int op_code = CGO_SHADER_CYLINDER_WITH_2ND_COLOR;
      shadercylinder2ndcolor(CGO *I, const float *_origin, const float *_axis, const float _radius,
                             int _cap, const float *_color2, Pickable *pickcolor2 = NULL,
                             const float alpha = -1.f);
      float origin[3], axis[3], tube_size;
      int cap;
      float color2[3];
      unsigned int pick_color_index;
      int pick_color_bond;
      float alpha;
    };

    struct sausage {
      static const int op_code = CGO_SAUSAGE;
      sausage(const float *_vertex1, const float *_vertex2, const float _radius, const float *_color1, const float *_color2) :
      radius(_radius) {
          copy3f(_vertex1, vertex1);
          copy3f(_vertex2, vertex2);
          copy3f(_color1, color1);
          copy3f(_color2, color2);
      };
      float vertex1[3], vertex2[3], radius, color1[3], color2[3];
    };

    struct cylinder {
      static const int op_code = CGO_CYLINDER;
      cylinder(const float *_vertex1, const float *_vertex2, const float _radius, const float *_color1, const float *_color2) :
      radius(_radius) {
          copy3f(_vertex1, vertex1);
          copy3f(_vertex2, vertex2);
          copy3f(_color1, color1);
          copy3f(_color2, color2);
      };
      float vertex1[3], vertex2[3], radius, color1[3], color2[3];
    };

    struct custom_cylinder {
      static const int op_code = CGO_CUSTOM_CYLINDER;
    custom_cylinder(const float *_vertex1, const float *_vertex2, const float _radius, const float *_color1, const float *_color2, const float _cap1, const float _cap2) :
      radius(_radius), cap1(_cap1), cap2(_cap2) {
          copy3f(_vertex1, vertex1);
          copy3f(_vertex2, vertex2);
          copy3f(_color1, color1);
          copy3f(_color2, color2);
      };
      float vertex1[3], vertex2[3], radius, color1[3], color2[3];
      float cap1, cap2;

      cCylCap get_cap1() const { return static_cast<cCylCap>(int(cap1)); }
      cCylCap get_cap2() const { return static_cast<cCylCap>(int(cap2)); }
    };

    struct custom_cylinder_alpha {
      static const int op_code = CGO_CUSTOM_CYLINDER_ALPHA;
      custom_cylinder_alpha(const float *_vertex1, const float *_vertex2, const float _radius, const float *_color1, const float *_color2,
                            const float _alpha1, const float _alpha2,
                            const float _cap1, const float _cap2) :
      radius(_radius), cap1(_cap1), cap2(_cap2) {
          copy3f(_vertex1, vertex1);
          copy3f(_vertex2, vertex2);
          copy3f(_color1, color1);
          copy3f(_color2, color2);
          color1[3] = _alpha1;
          color2[3] = _alpha2;
      };
      float vertex1[3], vertex2[3], radius, color1[4], color2[4];
      float cap1, cap2;

      cCylCap get_cap1() const { return static_cast<cCylCap>(int(cap1)); }
      cCylCap get_cap2() const { return static_cast<cCylCap>(int(cap2)); }
    };
  };
};

class CGO {
public:
  CGO(PyMOLGlobals* G, int size = 0);
  ~CGO();

  CGO(CGO const& other) = delete;

  PyMOLGlobals *G { nullptr };
  float *op { nullptr };
  size_t c = 0;
  bool z_flag = false;
  float z_min { 0 }, z_max { 0 };
  float z_vector[3];
  float alpha { 1.f };
  int *i_start { 0 }, i_size { 0 };
  unsigned int current_pick_color_index { 0 };
  int current_pick_color_bond { cPickableNoPick };
  float current_accessibility { 1.f };
  float normal[3]{0.f, 0.f, 1.f};
  float color[3]{};
  float texture[2];
  uchar pickColor[4]{0, 0, 0, 0xff};
  bool has_begin_end { false };
  bool has_draw_buffers { false }, has_draw_cylinder_buffers { false }, has_draw_sphere_buffers { false };
  bool use_shader { false }, cgo_shader_ub_color { false }, cgo_shader_ub_normal { false };
  bool debug { false };
  bool no_pick { false };
  short render_alpha { 0 };  // 1 : render CGOSetZVector/CGORenderGLAlpha only
                       // 2 : render both CGOSetZVector/CGORenderGLAlpha and rest of object
                       // calcDepth=1 by default
  short sphere_quality { 0 }; // quality of spheres when simplified or rendered in immediate mode
  bool interpolated { false };
  /***********************************************************************
   * CGO iterator
   *
   * for (auto it = cgo->begin(); !it.is_stop(); ++it) {
   *   auto pc = it.data();
   *   int op = it.op_code();
   *   ...
   * }
   ***********************************************************************/

  class const_iterator {
    protected:
      const float * m_pc;
      const float * m_stop;
    public:
      unsigned op_code() const {
        return *reinterpret_cast<const unsigned*>(m_pc);
      }
      operator int() const { return op_code(); }

      const float * data() const { return m_pc + 1; }

      template <typename T>
      const T * cast() const { return reinterpret_cast<const T *>(m_pc + 1); }

      const_iterator(const CGO * cgo) {
        m_pc = cgo->op;
        m_stop = cgo->op + cgo->c;
      }

      const_iterator& operator++();

      bool is_stop() const {
        return m_pc == m_stop || op_code() == CGO_STOP;
      }
  };

  class iterator : public const_iterator {
    public:
      iterator(CGO * cgo) : const_iterator(cgo) {}
      float * data() { return const_cast<float*>(m_pc + 1); }
      template <typename T> T* cast() { return reinterpret_cast<T*>(data()); }
  };

  const_iterator begin() const { return this; }
  iterator begin() { return this; }

  /***********************************************************************
   * This is the add function, the signature may look weird but it's
   * common c++11 perfect forwarding. This function passes the constructor
   * arguments into the function. Meaning it is used like:
   * cgo.add<cgo::draw::arrays>(mode, arrays, nverts);
   ***********************************************************************/
  template <typename T, typename... TArgs> float * add(TArgs&&... args) {
    int size = fsizeof<T>() + 1;
    float * at = add_to_buffer(size);
    // write the op code
    CGO_write_int(at, T::op_code);
    // call the type constructor in place forwarding the args to the constructor
    T * sp = new (at) T(std::forward<TArgs>(args)...);

    if (std::is_base_of<cgo::draw::op_with_draw_buffers, T>::value) {
      has_draw_buffers = true;
    }

    if (std::is_base_of<cgo::draw::op_with_data, T>::value) {
      auto ptr = reinterpret_cast<cgo::draw::op_with_data *>(sp);
      // set the buffer data flag
      // create the floating point data for this type
      auto data_len = ptr->get_data_length();
      if (data_len) {
        ptr->set_data(allocate_in_data_heap(data_len));
        return ptr->get_data();
      }
    }

    // this op does not have dynamic data, but since the caller knows
    // then we return the allocated pointer, so that the caller can determine
    // whether this was successful or not.
    return at;
  }

  // Appends the source CGO onto this CGO
  bool append(const CGO& source, bool stopAtEnd = false);
  void move_append(CGO&& source);
  void free_append(CGO * &&source);
  void free_append(CGO * &source);

  // Allocates in our CGO data pool
  float * allocate_in_data_heap(size_t size) {
    std::unique_ptr<float[]> uni(new float[size]);
    float * ptr = uni.get();
    _data_heap.emplace_back(std::move(uni));
    return ptr;
  }

  // templated by the op type
  template <typename T> void copy_op_from(const float * pc) {
    // copy the op
    const size_t op_size = fsizeof<T>() + 1; // + 1 is the op
    float * at = add_to_buffer(op_size);
    memcpy(at, (pc - 1), op_size * 4);

    if (std::is_base_of<cgo::draw::op_with_draw_buffers, T>::value) {
      has_draw_buffers = true;
    }

    if (std::is_base_of<cgo::draw::op_with_data, T>::value) {
      // copy the float data
      float * vals { nullptr };
      auto ptr = reinterpret_cast<const cgo::draw::op_with_data*>(pc);
      int data_len = ptr->get_data_length();
      if (data_len) {
        vals = allocate_in_data_heap(data_len);
        memcpy(vals, ptr->get_data(), data_len * 4);
      }
      auto spop = reinterpret_cast<cgo::draw::op_with_data *>(at + 1);
      spop->set_data( vals );
    }
  }

  // Our CGO_add, allocates size bytes at end of cgo buffer
  float * add_to_buffer(int size) {
    float * at { nullptr };
    VLACheck(op, float, size + c);
    if (!op)
      return nullptr;
    at = op + c;
    c += size;
    return at;
  }

  /***********************************************************************
   * This function adds to the end of this CGO the cgo op that exists in
   * the float. This function will also increment the pointer passed to it
   * by the size of the operation.
   ***********************************************************************/
  void add_to_cgo(int, const float*);

  // Pretty prints a table with the layout of this CGO
  void print_table() const;

private:
  std::vector<std::unique_ptr<float[]>> _data_heap;
};

int CGORendererInit(PyMOLGlobals * G);
void CGORendererFree(PyMOLGlobals * G);
#define CGONew new CGO
#define CGONewSized CGONew
int CGOGetExtent(const CGO * I, float *mn, float *mx);
int CGOHasNormals(const CGO * I);

void CGOFree(CGO * &I, bool withVBOs=true);
#define CGOFreeWithoutVBOs(I) CGOFree(I, false)

CGO *CGODrawText(const CGO * I, int est, float *camera);

CGO* CGOSimplify(const CGO* I, int est = 0, short sphere_quality = -1,
    bool stick_round_nub = true);
CGO *CGOSimplifyNoCompress(const CGO * I, int est, short sphere_quality = -1, bool stick_round_nub = true);

// -1 - no lines, 0 - some no interpolation, 1 - all interpolation, 2 - all no interpolation
bool CGOCombineBeginEnd(CGO ** I, bool do_not_split_lines = false);
CGO* CGOCombineBeginEnd(const CGO* I, int est = 0, bool do_not_split_lines = false);

void CGOFreeVBOs(CGO *I);

CGO *CGOOptimizeToVBOIndexed(const CGO * I, int est=0, const float *color=NULL, bool addshaders=true, bool embedTransparencyInfo=false);
#define CGOOptimizeToVBOIndexedWithColorEmbedTransparentInfo(I, est, color, addshaders) CGOOptimizeToVBOIndexed(I, est, color, addshaders, true)
#define CGOOptimizeToVBOIndexedWithColor CGOOptimizeToVBOIndexed
#define CGOOptimizeToVBOIndexedNoShader(I, est) CGOOptimizeToVBOIndexed(I, est, NULL, false)

bool CGOOptimizeToVBONotIndexed(CGO ** I);
CGO* CGOOptimizeToVBONotIndexed(const CGO* I, int est = 0,
    bool addshaders = true, float** returnedData = nullptr);

#define CGOOptimizeToVBONotIndexedWithReturnedData CGOOptimizeToVBONotIndexed
#define CGOOptimizeToVBONotIndexedNoShader(I) CGOOptimizeToVBONotIndexed(I, 0, false)


CGO *CGOOptimizeSpheresToVBONonIndexed(const CGO * I, int est=0, bool addshaders=false, CGO *leftOverCGO=NULL);
#define CGOOptimizeSpheresToVBONonIndexedNoShader(I, est) CGOOptimizeSpheresToVBONonIndexed(I, est, false, NULL)

int CGOCheckComplex(CGO * I);
int CGOPreloadFonts(CGO * I);

int CGOCheckForText(CGO * I);

int CGOFromFloatArray(CGO * I, const float *src, int len);

int CGOBegin(CGO * I, int mode);
int CGOEnd(CGO * I);

int CGOSphere(CGO * I, const float *v1, float r);
int CGOEllipsoid(CGO * I, const float *v1, float r, const float *n1, const float *n2, const float *n3);
int CGOVertex(CGO * I, float v1, float v2, float v3);
int CGOVertexv(CGO * I, const float *v);
int CGOVertexCrossv(CGO * I, const float *v);
int CGOAlpha(CGO * I, float alpha);
int CGOColor(CGO * I, float v1, float v2, float v3);
int CGOColorv(CGO * I, const float *v);
int CGOTexCoord2f(CGO * I, float v1, float v2);
int CGONormal(CGO * I, float v1, float v2, float v3);
int CGONormalv(CGO * I, const float *v);
int CGOResetNormal(CGO * I, int mode);
int CGOLinewidth(CGO * I, float v);
int CGOSpecial(CGO * I, int v);
// all pre-processor definitions for CGOSpecial ops
enum {
  LINEWIDTH_DYNAMIC_WITH_SCALE = 1,
  LINEWIDTH_DYNAMIC_MESH,
  POINTSIZE_DYNAMIC_DOT_WIDTH,
  LINEWIDTH_DYNAMIC_WITH_SCALE_RIBBON,
  LINEWIDTH_DYNAMIC_WITH_SCALE_DASH,
  CYLINDERWIDTH_DYNAMIC_MESH,
  LINEWIDTH_WITH_SCALE,
  DOTSIZE_WITH_SPHERESCALE,
  MESH_WIDTH_FOR_SURFACES,
  CYLINDER_WIDTH_FOR_DISTANCES,
  CYLINDER_WIDTH_FOR_RIBBONS,
  DOT_WIDTH_FOR_DOTS,
  DOT_WIDTH_FOR_DOT_SPHERES,
  CYLINDER_WIDTH_FOR_NONBONDED,
  CYLINDER_WIDTH_FOR_REPWIRE,
  CYLINDER_WIDTH_FOR_REPWIRE_VARWIDTH,
  ENABLE_BACK_FACES_IF_NOT_TWO_SIDED,
  DISABLE_BACK_FACES_IF_NOT_TWO_SIDED,
  SET_SURFACE_UNIFORMS,
  SET_ALIGNMENT_UNIFORMS_ATTRIBS,
  LINEWIDTH_FOR_LINES,
  SET_LABEL_SCALE_UNIFORMS
};

int CGOSpecialWithArg(CGO * I, int v, float arg);
#define SPHERE_MODE_OPS 1
#define LINE_LIGHTING  2

int CGODotwidth(CGO * I, float v);
int CGOFontVertex(CGO * I, float x, float y, float z);
int CGOFontScale(CGO * I, float v1, float v2);
int CGOIndent(CGO * I, char c, float dir);
int CGOWrite(CGO * I, const char *str);
int CGOWriteLeft(CGO * I, const char *str);
int CGOWriteIndent(CGO * I, const char *str, float indent);

#define CGODrawArrays(this, ...) (this)->add<cgo::draw::arrays>(__VA_ARGS__)
#define CGOCustomCylinderv(this, ...) (this)->add<cgo::draw::custom_cylinder>(__VA_ARGS__)
#define CGOShaderCylinder(this, ...) (this)->add<cgo::draw::shadercylinder>(__VA_ARGS__)
#define CGOCylinderv(this, ...) (this)->add<cgo::draw::cylinder>(__VA_ARGS__)

int CGOBoundingBox(CGO *I, const float *min, const float *max);
int CGOAccessibility(CGO * I, const float a);

int CGODrawTexture(CGO *I, int texture_id, float *worldPos, float *screenMin, float *screenMax, float *textExtent);
int CGODrawConnector(CGO *I, float *targetPt3d, float *labelCenterPt3d, float text_width, float text_height, float *screenOffset, float *screenWorldOffset, float *connectorColor, short relativeMode, int draw_bkgrd, float bkgrd_transp, float *bkgrd_color, float rel_ext_length, float connectorWidth);
CGO *CGOOptimizeLabels(const CGO * I, int est, bool addshaders=false);
CGO *CGOOptimizeTextures(const CGO * I, int est);
CGO *CGOExpandDrawTextures(const CGO * I, int est);
CGO *CGOOptimizeConnectors(const CGO * I, int est);

/*void CGOFontScale(CGO *I,float v);
  void CGOFont(CGO *I,float size,int face,int style);*/

void CGORoundNub(CGO * I,
    const float *v1,    // cap center
    const float *p0,    // normal along axis
    const float *p1,    // x coord in cap space
    const float *p2,    // y coord in cap space
    int direction,      // 1 or -1
    int nEdge,          // "quality"
    float size);

int CGOEnable(CGO * I, int mode);
int CGODisable(CGO * I, int mode);

int CGOStop(CGO * I);

int CGOAlphaTriangle(CGO * I,
		     const float *v1, const float *v2, const float *v3,
		     const float *n1, const float *n2, const float *n3,
		     const float *c1, const float *c2, const float *c3,
		     float a1, float a2, float a3, int reverse);
void CGOSetZVector(CGO * I, float z0, float z1, float z2);
struct GadgetSet;
void CGORenderGLPicking(CGO * I, RenderInfo *info,
                        PickContext * context, CSetting * set1, CSetting * set2, Rep *rep=NULL);
void CGORenderGL(CGO * I, const float *color, CSetting * set1, CSetting * set2,
                 RenderInfo * info, Rep *rep);
void CGORenderGLAlpha(CGO * I, RenderInfo * info, bool calcDepth);
int CGORenderRay(CGO * I, CRay * ray, RenderInfo * info, const float *color, ObjectGadgetRamp *ramp, CSetting * set1, CSetting * set2);
void CGOReset(CGO * I);

void CGOSetUseShader(CGO *I, int use_shader);

PyObject *CGOAsPyList(CGO * I);
CGO *CGONewFromPyList(PyMOLGlobals * G, PyObject * list, int version, bool shouldCombine=true);
int CGOPickColor(CGO * I, unsigned int index, int bond);

const cgo::draw::buffers_not_indexed* CGOGetNextDrawBufferedNotIndex(
    const CGO*);

float *CGOGetNextOp(float *cgo_op, int optype);

int CGOAppend(CGO *dest, const CGO *source, bool stopAtEnd=true);
inline int CGOAppendNoStop(CGO *dest, const CGO *source) {
  return CGOAppend(dest, source, false);
}

int CGOCountNumberOfOperationsOfType(const CGO *I, int op);
int CGOCountNumberOfOperationsOfTypeN(const CGO *I, const std::set<int> &optype);
int CGOCountNumberOfOperationsOfTypeN(const CGO *I, const std::map<int, int> &optype);
bool CGOHasOperations(const CGO *I);
bool CGOHasOperationsOfType(const CGO *I, int optype);
bool CGOHasOperationsOfTypeN(const CGO *I, const std::set<int> &optype);
bool CGOHasCylinderOperations(const CGO *I);
bool CGOHasSphereOperations(const CGO *I);
bool CGOFilterOutCylinderOperationsInto(const CGO *I, CGO *cgo);

bool CGOCheckWhetherToFree(PyMOLGlobals * G, CGO *I);

CGO *CGOConvertLinesToShaderCylinders(const CGO * I, int est);
CGO *CGOConvertLinesToTrilines(const CGO * I, bool addshaders=true);
CGO *CGOConvertToLabelShader(const CGO *I, CGO * addTo);

void CGOChangeShadersTo(CGO *I, int frommode, int tomode);
CGO *CGOOptimizeScreenTexturesAndPolygons(CGO * I, int est);
CGO *CGOColorByRamp(PyMOLGlobals * G, const CGO *I, ObjectGadgetRamp *ramp, int state, CSetting * set1);

#define CGOLineAsTriangleStrips(CGO, minx, miny, maxx, maxy) \
	CGOBegin(CGO, GL_TRIANGLE_STRIP);         \
	CGOVertex(CGO, minx, miny, 0.f);         \
	CGOVertex(CGO, minx, maxy+1, 0.f);         \
	CGOVertex(CGO, minx+1, miny, 0.f);         \
	CGOVertex(CGO, minx+1, maxy+1, 0.f);         \
	CGOEnd(CGO);         \
	CGOBegin(CGO, GL_TRIANGLE_STRIP);         \
	CGOVertex(CGO, minx, maxy, 0.f);         \
	CGOVertex(CGO, minx, maxy+1, 0.f);         \
	CGOVertex(CGO, maxx, maxy, 0.f);         \
	CGOVertex(CGO, maxx, maxy+1, 0.f);         \
	CGOEnd(CGO);         \
	CGOBegin(CGO, GL_TRIANGLE_STRIP);         \
	CGOVertex(CGO, maxx, miny, 0.f);         \
	CGOVertex(CGO, maxx, maxy+1, 0.f);         \
	CGOVertex(CGO, maxx+1, miny, 0.f);         \
	CGOVertex(CGO, maxx+1, maxy+1, 0.f);         \
	CGOEnd(CGO);         \
	CGOBegin(CGO, GL_TRIANGLE_STRIP);         \
	CGOVertex(CGO, minx, miny, 0.f);         \
	CGOVertex(CGO, minx, miny+1, 0.f);         \
	CGOVertex(CGO, maxx, miny, 0.f);         \
	CGOVertex(CGO, maxx, miny+1, 0.f);         \
	CGOEnd(CGO);

int CGOHasTransparency(const CGO *I, bool checkTransp=true, bool checkOpaque=false);
#define CGOHasOpaque(I) CGOHasTransparency(I, false, true)

CGO *CGOConvertTrianglesToAlpha(const CGO * I);
CGO *CGOGenerateNormalsForTriangles(const CGO * I);

bool CGOHasAnyTriangleVerticesWithoutNormals(const CGO *I, bool checkTriangles=true);
#define CGOHasAnyLineVerticesWithoutNormals(I) CGOHasAnyTriangleVerticesWithoutNormals(I, false)

CGO* CGOTurnLightingOnLinesOff(const CGO* I, bool use_shader);
inline CGO* CGOTurnLightingOnLinesOff(const CGO* I)
{
  return CGOTurnLightingOnLinesOff(I, I->use_shader);
}

// returns offset of floats in CGO array
int CGOUniform3f(CGO *I, int uniform_id, const float *value);

CGO *CGOConvertSpheresToPoints(const CGO *I);

CGO *CGOConvertToShader(const CGO *I, AttribDataDesc &attrData, AttribDataDesc &pickData, int mode, const VertexBuffer::buffer_layout layout=VertexBuffer::INTERLEAVED, bool check_attr_for_data=true, int *idx_array=NULL, int nindicesperfrag=0, int nfragspergroup = 1);

bool CGOCheckSplitLineInterpolationIsSame(const CGO *I, bool &interp_value);

CGO *CGOConvertToTrilinesShader(const CGO *I, CGO *addTo, bool add_color=true);
CGO *CGOConvertToLinesShader(const CGO *I, CGO *addTo, bool add_color=true);

CGO *CGOConvertLinesToCylinderShader(const CGO *I, CGO *addTo, bool add_color = true);
CGO *CGOConvertCrossesToCylinderShader(const CGO *I, CGO *addTo, float cross_size);
CGO *CGOConvertCrossesToLinesShader(const CGO *I, CGO *addTo, float cross_size);
CGO *CGOConvertCrossesToTrilinesShader(const CGO *I, CGO *addTo, float cross_size);
CGO *CGOConvertShaderCylindersToCylinderShader(const CGO *I, CGO *addTo);

void AssignNewPickColor(CGO* cgo, PickColorManager*, unsigned char* color,
    const PickContext* context, unsigned int index, int bond);

#endif
