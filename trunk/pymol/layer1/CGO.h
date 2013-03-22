
/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* copyright 1998-2000 by Warren Lyford Delano of DeLano Scientific. 
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
#include"Rep.h"


/* Compiled Graphics Library for simple graphics objects
   in floating point three-space, with the goal of achieving
   quick and easy rendering in multiple environments without the
   headaches of OpenGL arrays.

*/

#ifdef _PYMOL_CGO_DRAWARRAYS
#define CGO_read_int(p) (*((int*)((p)++)))
#else
#define CGO_read_int(p) (*((int*)(p++)))
#endif
#define CGO_get_int(p) (*((int*)(p)))
#define CGO_write_int(p,i) ((*((int*)(p++)))=(i))
#define CGO_put_int(p,i) ((*((int*)(p)))=(i))

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

struct _CGO {
  PyMOLGlobals *G;
  float *op;
  int c;
  int z_flag;
  float z_min, z_max;
  float z_vector[3];
  float alpha;
  int *i_start, i_size;
#ifdef _PYMOL_CGO_DRAWARRAYS
  short has_begin_end;
#endif
  int current_pick_color_index, current_pick_color_bond;
#ifdef _PYMOL_CGO_DRAWBUFFERS
  float current_accessibility;
  short has_draw_buffers, has_draw_cylinder_buffers, has_draw_sphere_buffers;
  float normal[3], color[3], texture[2];
  uchar pickColor[4];
#endif
  short use_shader, cgo_shader_ub_color, cgo_shader_ub_normal;
  short debug;
  short enable_shaders;
  short no_pick;
};


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
#ifdef _PYMOL_CGO_DRAWARRAYS
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
#define CGO_DRAW_ARRAYS_SZ       0
#endif
#ifdef _PYMOL_CGO_DRAWBUFFERS
/* CGO_DRAW_BUFFERS : operation that uses glDrawArrays with VBOs.  This is not currently
   used, since CGO_DRAW_BUFFERS_INDEXED is a bit more flexible where an index can be specified
   for vertex order.  However, this is useful for when all primitives are defined only once 
   and an index is not needed.
   This operation is similar to CGO_DRAW_ARRAYS except since it is using VBOs, the length 
   of the operation is constant (i.e., 8 floats, CGO_DRAW_BUFFERS_SZ).  The VBO ids are stored
   and are used to manage the data inside the buffers using glBindBuffer()/glBufferData()
   - mode : GL Mode that is used
   - arrays : which arrays that are used (bitmask on CGO_<type>_ARRAY, where <type> 
              can be VERTEX, NORMAL, COLOR, or PICK_COLOR)
   - narrays : number of arrays specified in arrays
   - nverts : number of total vertices specified in arrays
   - bufs[4] : each VBO id in order (if used, VERTEX, NORMAL, COLOR, PICK_COLOR)
 */
#define CGO_DRAW_BUFFERS         0x20
#define CGO_DRAW_BUFFERS_SZ      8
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
#define CGO_DRAW_BUFFERS_INDEXED_SZ      0
/* CGO_BOUNDING_BOX : operation that allows the extent to be expanded.  Since the geometry
   data is not kept in the CGO object for VBO objects (only on the card), this object allows 
   the bounding box of the object to be saved in the CGO.  This is used in the CGOGetExtent(),
   typically when the view is being automatically set, but can be used for other things.
 */
#define CGO_BOUNDING_BOX         0x22
#define CGO_BOUNDING_BOX_SZ      6

#define CGO_DRAW_BUFFERS_NOT_INDEXED         0x23
#define CGO_DRAW_BUFFERS_NOT_INDEXED_SZ      0

#define CGO_LINEWIDTH_SPECIAL            0x24
#define CGO_LINEWIDTH_SPECIAL_SZ         1

#define CGO_DRAW_CYLINDER_BUFFERS       0x25
#define CGO_DRAW_CYLINDER_BUFFERS_SZ    7

#endif
#define CGO_SHADER_CYLINDER             0x26
#define CGO_SHADER_CYLINDER_SZ          8

#define CGO_SHADER_CYLINDER_WITH_2ND_COLOR      0x27
#define CGO_SHADER_CYLINDER_WITH_2ND_COLOR_SZ    11

#define CGO_DRAW_SPHERE_BUFFERS      0x28
#define CGO_DRAW_SPHERE_BUFFERS_SZ    5

#define CGO_ACCESSIBILITY      0x29
#define CGO_ACCESSIBILITY_SZ    1

#define CGO_DRAW_TEXTURE      0x2A
#define CGO_DRAW_TEXTURE_SZ    13

#define CGO_DRAW_TEXTURES      0x2B
#define CGO_DRAW_TEXTURES_SZ    0

#define CGO_DRAW_SCREEN_TEXTURES_AND_POLYGONS      0x2C
#define CGO_DRAW_SCREEN_TEXTURES_AND_POLYGONS_SZ    4

#define CGO_TEX_COORD                0x2D
#define CGO_TEX_COORD_SZ             2


#define CGO_DRAW_LABEL      0x2E
#ifdef PYMOL_TEXT_IN_ONE_TEXTURE
#define CGO_DRAW_LABEL_SZ    16
#else
#define CGO_DRAW_LABEL_SZ    19
#endif

#define CGO_DRAW_LABELS      0x2F
#define CGO_DRAW_LABELS_SZ    0


#define CGO_MASK                 0x3F


#define CGO_LIGHTING             0x0B50

#define CGO_VERTEX_ARRAY         0x01
#define CGO_NORMAL_ARRAY         0x02
#define CGO_COLOR_ARRAY          0x04
#define CGO_PICK_COLOR_ARRAY     0x08
#define CGO_ACCESSIBILITY_ARRAY  0x10
#define CGO_TEX_COORD_ARRAY      0x20

int CGORendererInit(PyMOLGlobals * G);
void CGORendererFree(PyMOLGlobals * G);
CGO *CGONew(PyMOLGlobals * G);
CGO *CGONewSized(PyMOLGlobals * G, int size);
int CGOGetExtent(CGO * I, float *mn, float *mx);
int CGOHasNormals(CGO * I);

void CGOFreeWithoutVBOs(CGO * I);
void CGOFree(CGO * I);
void CGOFreeImpl(CGO * I, short withVBOs);
CGO *CGODrawText(CGO * I, int est, float *camera);

CGO *CGOSimplify(CGO * I, int est);

CGO *CGOCombineBeginEnd(CGO * I, int est);
#ifdef _PYMOL_CGO_DRAWBUFFERS
void CGOFreeVBOs(CGO *I);
CGO *CGOOptimizeToVBOIndexedWithColor(CGO * I, int est, float *color);
CGO *CGOOptimizeToVBOIndexedNoShader(CGO * I, int est);
CGO *CGOOptimizeToVBOIndexed(CGO * I, int est);
CGO *CGOOptimizeToVBONotIndexedWithReturnedData(CGO * I, int est, short, float **);
CGO *CGOOptimizeToVBONotIndexed(CGO * I, int est);
CGO *CGOOptimizeSpheresToVBONonIndexedImpl(CGO * I, int est, CGO *leftOverCGO);
CGO *CGOOptimizeSpheresToVBONonIndexed(CGO * I, int est);
#endif

void CGOReserve(CGO * ptr, int est);

int CGOCheckComplex(CGO * I);
int CGOPreloadFonts(CGO * I);

int CGOCheckForText(CGO * I);

int CGOFromFloatArray(CGO * I, float *src, int len);

int CGOBegin(CGO * I, int mode);
int CGOEnd(CGO * I);

int CGOSphere(CGO * I, float *v1, float r);
int CGOEllipsoid(CGO * I, float *v1, float r, float *n1, float *n2, float *n3);
int CGOQuadric(CGO * I, float *v1, float r, float *p); /* NOT WORKING YET */
int CGOSausage(CGO * I, float *v1, float *v2, float r, float *c1, float *c2);
int CGOVertex(CGO * I, float v1, float v2, float v3);
int CGOVertexv(CGO * I, float *v);
int CGOAlpha(CGO * I, float alpha);
int CGOColor(CGO * I, float v1, float v2, float v3);
int CGOColorv(CGO * I, float *v);
int CGOTexCoord2f(CGO * I, float v1, float v2);
int CGOTexCoord2fv(CGO * I, float *v);
int CGONormal(CGO * I, float v1, float v2, float v3);
int CGONormalv(CGO * I, float *v);
int CGOResetNormal(CGO * I, int mode);
int CGOLinewidth(CGO * I, float v);
int CGOLinewidthSpecial(CGO * I, int v);
#define LINEWIDTH_DYNAMIC_WITH_SCALE 1
#define LINEWIDTH_DYNAMIC_MESH 2
#define POINTSIZE_DYNAMIC_DOT_WIDTH 3
#define LINEWIDTH_DYNAMIC_WITH_SCALE_RIBBON 4
#define LINEWIDTH_DYNAMIC_WITH_SCALE_DASH 5
#define CYLINDERWIDTH_DYNAMIC_MESH  6
int CGODotwidth(CGO * I, float v);
int CGOChar(CGO * I, char c);
int CGOFontVertex(CGO * I, float x, float y, float z);
int CGOFontVertexv(CGO * I, float *v);
int CGOFontScale(CGO * I, float v1, float v2);
int CGOIndent(CGO * I, char c, float dir);
int CGOWrite(CGO * I, char *str);
int CGOWriteLeft(CGO * I, char *str);
int CGOWriteIndent(CGO * I, char *str, float indent);

GLfloat *CGODrawArrays(CGO *I, GLenum mode, short arrays, int nverts);

#ifdef _PYMOL_CGO_DRAWBUFFERS
int CGODrawBuffers(CGO *I, GLenum mode, short arrays, int nverts, uint *bufs);
GLfloat *CGODrawBuffersIndexed(CGO *I, GLenum mode, short arrays, int nindices, int nverts, uint *bufs);
int CGOBoundingBox(CGO *I, float *min, float *max);
int CGOAccessibility(CGO * I, float a);
#endif
int CGODrawTexture(CGO *I, int texture_id, float *worldPos, float *screenMin, float *screenMax, float *textExtent);
int CGODrawLabel(CGO *I, int texture_id, float *worldPos, float *screenWorldOffset, float *screenMin, float *screenMax, float *textExtent);
CGO *CGOOptimizeLabels(CGO * I, int est);
CGO *CGOOptimizeTextures(CGO * I, int est);
CGO *CGOExpandDrawTextures(CGO * I, int est);

/*void CGOFontScale(CGO *I,float v);
  void CGOFont(CGO *I,float size,int face,int style);*/

int CGOEnable(CGO * I, int mode);
int CGODisable(CGO * I, int mode);

int CGOStop(CGO * I);

int CGOCylinderv(CGO * I, float *p1, float *p2, float r, float *c1, float *c2);
int CGOCustomCylinderv(CGO * I, float *p1, float *p2, float r, float *c1, float *c2,
                        float cap1, float cap2);
int CGOConev(CGO * I, float *p1, float *p2, float r1, float r2, float *c1, float *c2,
              float cap1, float cap2);

int CGOAlphaTriangle(CGO * I,
		     float *v1, float *v2, float *v3,
		     float *n1, float *n2, float *n3,
		     float *c1, float *c2, float *c3,
		     float a1, float a2, float a3, int reverse);
void CGOSetZVector(CGO * I, float z0, float z1, float z2);
void CGORenderGLPicking(CGO * I, Picking ** pick,
                        PickContext * context, CSetting * set1, CSetting * set2);
void CGORenderGL(CGO * I, float *color, CSetting * set1, CSetting * set2,
                 RenderInfo * info, Rep *rep);
void CGORenderGLAlpha(CGO * I, RenderInfo * info);
int CGORenderRay(CGO * I, CRay * ray, float *color, CSetting * set1, CSetting * set2);
void CGOReset(CGO * I);

void CGOSetUseShader(CGO *I, int use_shader);

PyObject *CGOAsPyList(CGO * I);
CGO *CGONewFromPyList(PyMOLGlobals * G, PyObject * list, int version);
void SetCGOPickColor(float *colorVals, int nverts, int pl, int index, int bond);
int CGOPickColor(CGO * I, int index, int bond);
float *CGO_add_GLfloat(CGO * I, int c);

float *CGOGetNextDrawBufferedIndex(float *cgo_op);
float *CGOGetNextDrawBufferedNotIndex(float *cgo_op);
float *CGOGetNextDrawBufferedImpl(float *cgo_op, int optype);
float *CGOGetNextOp(float *cgo_op, int optype);

int CGOAppendNoStop(CGO *dest, CGO *source);
int CGOAppend(CGO *dest, CGO *source);

CGO *CGOOptimizeGLSLCylindersToVBOIndexed(CGO * I, int est);
CGO *CGOOptimizeGLSLCylindersToVBOIndexedWithLeftOver(CGO * I, int est, CGO *leftOverCGO);
CGO *CGOOptimizeGLSLCylindersToVBOIndexedNoColor(CGO * I, int est);

int CGOShaderCylinder(CGO *I, float *origin, float *axis, float tube_size, int cap);
int CGOShaderCylinder2ndColor(CGO *I, float *origin, float *axis, float tube_size, int cap, float *color2);

int CGOCountNumberOfOperationsOfTypeDEBUG(CGO *I, int optype);
int CGOCountNumberOfOperationsOfType(CGO *I, int op);
short CGOHasOperationsOfType(CGO *I, int op);
short CGOHasOperationsOfType2(CGO *I, int op1, int op2);
short CGOHasCylinderOperations(CGO *I);

short CGOCheckWhetherToFree(PyMOLGlobals * G, CGO *I);

CGO *CGOConvertLinesToShaderCylinders(CGO * I, int est);

int CGOChangeShadersTo(CGO *I, int frommode, int tomode);
void CGOCountNumVerticesDEBUG(CGO *I);
CGO *CGOOptimizeScreenTexturesAndPolygons(CGO * I, int est);

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


#endif
