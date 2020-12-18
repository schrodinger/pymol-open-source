
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
#include"os_python.h"

#include"os_predef.h"
#include"os_std.h"
#include"os_gl.h"

#include"Base.h"
#include"OOMac.h"
#include"RepLabel.h"
#include"Color.h"
#include"Map.h"
#include"Setting.h"
#include"main.h"
#include"Scene.h"
#include"Text.h"
#include"Matrix.h"
#include"ShaderMgr.h"
#include"CGO.h"

/**
 * Memory layout of the RepLabel::V array
 */
struct VItemType {
  float color[3];           //  0- 2: label_color
  float coord[3];           //  3- 5: atom coordinate
  float screen_point[3];    //  6- 8: coordinate + label_placement_offset or label_screen_point
  float connector_color[3]; //  9-11: label_connector_color
  float position[3];        // 12-14: label_position

  float relativeMode_f;     // 15: bitmask
  // = (drawConnector << 0)
  // | (isScreenCoord << 1)
  // | (isPixelCoord << 2)
  // | (at_label_z_target << 3)

  float spacing;            // 16: label_multiline_spacing
  float justification;      // 17: label_multiline_justification
  float padding[3];         // 18-20: label_padding

  float draw_var_f;         // 21: bitmask
  // (1 << 0) * label_connector +
  // (1 << 1) * label_bg +
  // (1 << 2) * label_bg_outline +
  // (1 << 3) * label_connector_mode_1 +
  // (1 << 4) * label_connector_mode_2 +
  // (1 << 5) * label_connector_mode_3 +
  // (1 << 6) * label_connector_mode_4 +
  // (1 << 7) * ray_label_connector_flat

  float bg_alpha;           // 22: 1.0 - label_bg_transparency
  float bg_color[3];        // 23-25: label_bg_color
  float connector_width;    // 26: label_connector_width
  float connector_ext_len;  // 27: label_connector_ext_length
};

struct RepLabel : Rep {
  using Rep::Rep;

  ~RepLabel() override;

  cRep_t type() const override { return cRepLabel; }
  void render(RenderInfo* info) override;

  // VItemType *V;
  float* V = nullptr;
  lexidx_t* L = nullptr;
  int N;
  int OutlineColor;
  CGO* shaderCGO = nullptr;
  int texture_font_size = 0;
};

#define SHADERCGO I->shaderCGO

#include"ObjectMolecule.h"

RepLabel::~RepLabel()
{
  auto I = this;
  FreeP(I->V);
  FreeP(I->L);
  CGOFree(I->shaderCGO);
}

#define MAX_LABEL_TEXTURE_SIZE 256
#define MAX_LABEL_FOR_ALWAYS_REFRESH 32
#define PERCENTAGE_CHANGE_FOR_REFRESH .2f
short InvalidateShaderCGOIfTextureNeedsUpdate(PyMOLGlobals *G, float font_size, int texture_font_size, int *sizeArg){
  float v_scale, perc = 0.f;
  int size, diff;
  short inv = false;
  v_scale = SceneGetScreenVertexScale(G, NULL);
  size = (int) (0.5F - font_size / v_scale);
  if (size <= 0)
    size = 1;
  if (size > MAX_LABEL_TEXTURE_SIZE){  // if label size is above max, just set to max size
    size = MAX_LABEL_TEXTURE_SIZE;
    if (texture_font_size!=size)
      inv = true;
  } else {
    if (font_size > 0 || size < MAX_LABEL_FOR_ALWAYS_REFRESH){  
      // if less than min, always refresh labels on texture if size changes
      if (texture_font_size!=size)
	inv = true;
    } else {
      // label is larger than min, but changes more than a percentage, refresh on texture
      diff = size - texture_font_size;
      diff = diff < 0 ? -diff : diff;
      perc = (diff / (float)size);
      inv = (perc > PERCENTAGE_CHANGE_FOR_REFRESH);
    }
  }
  *sizeArg = size;
  return (inv || !texture_font_size);
}

static float CLAMP_VALUE(float val, float min, float max){
  return ( val < min ? min : ( val > max ? max : val) );
}

static void RepLabelAdjustScreenZ(PyMOLGlobals *G, float *pt){
  short outside;
  outside = (pt[2] < -1.f || pt[2] > 1.f);
  if (!outside){
    pt[2] = CLAMP_VALUE(pt[2], -.9999, .97); // .98 b/c resolution in the back is lower
  }
}

static void addXYtoVertex(float x, float y, float *xv, float *yv, float *origv, float *v){
  float tmp3f[3];
  copy3f(origv, v);
  mult3f(xv, x, tmp3f);
  add3f(v, tmp3f, v);
  mult3f(yv, y, tmp3f);
  add3f(v, tmp3f, v);
}

static void addXYZtoVertex(float x, float y, float z, float *xv, float *yv,
    float *zv, float *origv, float *v)
{
  float tmp3f[3];
  copy3f(origv, v);
  mult3f(xv, x, tmp3f);
  add3f(v, tmp3f, v);
  mult3f(yv, y, tmp3f);
  add3f(v, tmp3f, v);
  mult3f(zv, z, tmp3f);
  add3f(v, tmp3f, v);
}

static void RayDrawLineAsGeometryWithOffsets(CRay *ray, float *pt1, float *pt2,
    float *spt1, float *spt2, float *xn, float *yn, float *zn, float line_width,
    float topext, float bottomext, float *color, float *dirv,
    unsigned char noLighting)
{
  float pts[4][3];
  float pt1E[3], pt2E[3];
  float tmpV[3], tmpV2[3], linev[3];
  float nzn[3] = { 0.f, 0.f, 1.f };
  copy3f(pt1, pt1E);
  copy3f(pt2, pt2E);
  subtract3f(spt1, spt2, tmpV);

  copy3f(tmpV, linev);
  normalize3f(linev);
  mult3f(linev, line_width, linev);
  
  cross_product3f(tmpV, nzn, tmpV2);
  normalize3f(tmpV2);
  mult3f(tmpV2, line_width, dirv);
  addXYtoVertex(dirv[0], dirv[1], xn, yn, pt1, pt1E);
  addXYtoVertex(topext * linev[0], topext * linev[1], xn, yn, pt1E, pts[0]);

  addXYtoVertex(dirv[0], dirv[1], xn, yn, pt2, pt2E);
  addXYtoVertex(-topext * linev[0], -topext * linev[1], xn, yn, pt2E, pts[1]);

  addXYtoVertex(-dirv[0], -dirv[1], xn, yn, pt1, pt1E);
  addXYtoVertex(bottomext * linev[0], bottomext * linev[1], xn, yn, pt1E, pts[2]);

  addXYtoVertex(-dirv[0], -dirv[1], xn, yn, pt2, pt2E);
  addXYtoVertex(-bottomext * linev[0], -bottomext * linev[1], xn, yn, pt2E, pts[3]);

  ray->triangle3fv(pts[0], pts[1], pts[2], zn, zn, zn, color, color, color);
  ray->setLastToNoLighting(noLighting);
  ray->triangle3fv(pts[1], pts[2], pts[3], zn, zn, zn, color, color, color);
  ray->setLastToNoLighting(noLighting);
}

#ifndef PURE_OPENGL_ES_2
static
void drawLine2DCross(float cw, float x1, float y1, float x2, float y2, float *cross){
  float lvect[3];
  float nzn[3] = { 0.f, 0.f, 1.f };

  lvect[0] = x2 - x1;
  lvect[1] = y2 - y1;
  normalize2f(lvect);
  cross_product3f(lvect, nzn, cross);
  mult3f(cross, cw, cross);
  glBegin(GL_TRIANGLE_STRIP);
  glVertex3f(x1 + cross[0], y1 + cross[1], 0.f);
  glVertex3f(x2 + cross[0], y2 + cross[1], 0.f);
  glVertex3f(x1 - cross[0], y1 - cross[1], 0.f);
  glVertex3f(x2 - cross[0], y2 - cross[1], 0.f);
  glEnd();
}

/* Draw Line/Polygon from point (x1,y1) to (x2,y2) with different Z's, where the current Z is 
   at world point curpt with the line starting at (x1,y1), and the second point (x2,y2) is 
   at the world point pt (i.e., the offset is embedded into the matrix convMatrix computed from
   SceneGenerateMatrixToAnotherZFromZ */
static
void drawLineToPointInWorldCross(PyMOLGlobals *G, float cw, float x1, float y1, float x2, float y2, float *cross, float *pt, float *curpt){
  float lvect[3];
  float nzn[3] = { 0.f, 0.f, 1.f };
  float convMatrix[16];
  float tmppt[3];
  SceneGenerateMatrixToAnotherZFromZ(G, convMatrix, curpt, pt);

  lvect[0] = x2 - x1;
  lvect[1] = y2 - y1;
  normalize2f(lvect);
  cross_product3f(lvect, nzn, cross);
  mult3f(cross, cw, cross);
  glBegin(GL_TRIANGLE_STRIP);

  glVertex3f(x1 + cross[0], y1 + cross[1], 0.f);

  // target point 1 in pt screen coordinates
  tmppt[0] = cross[0]; tmppt[1] = cross[1]; tmppt[2] = 0.f;
  MatrixTransformC44f3f(convMatrix, tmppt, tmppt);
  glVertex3fv(tmppt);

  glVertex3f(x1 - cross[0], y1 - cross[1], 0.f);

  // target point 2 in pt screen coordinates
  tmppt[0] = -cross[0]; tmppt[1] = -cross[1]; tmppt[2] = 0.f;
  MatrixTransformC44f3f(convMatrix, tmppt, tmppt);
  glVertex3fv(tmppt);

  glEnd();
}

#define CLIP_LEFT 1
#define CLIP_RIGHT 2
#define CLIP_TOP 4
#define CLIP_BOTTOM 8

static
short CLIPt(float denom, float num, float *tE, float *tL, short *clipedges, short bitmask){
  float t;
  if (denom > 0){
    t = num / denom;
    if (t > *tL)
      return 0;
    else if (t > *tE){
      *tE = t;
      *clipedges = bitmask;
    }
  } else if (denom < 0){
    t = num / denom;
    if (t < *tE)
      return 0;
    else if (t < *tL){
      *tL = t;
      *clipedges = bitmask;
    }
  } else if (num > 0)
    return 0;
  return 1;
}

/* This function clips the line inside of the (-xmax,-ymax, xmax, ymax) rectangle */
static
void Clip2D(float xmax, float ymax, float *x0, float *y0, float *x1, float *y1, short *visible, short *clipedges){
  float dx = *x1 - *x0;
  float dy = *y1 - *y0;
  *visible = 0;
  *clipedges = 0;
  if (dx == 0.f && dy == 0.f && fabs(*x0) < xmax && fabs(*y0) < ymax)
    *visible = 1;
  else {
    float tE = 0.f, tL = 1.f;
    if (CLIPt(dx, -xmax - *x0, &tE, &tL, clipedges, CLIP_LEFT))    // left
      if (CLIPt(-dx, *x0 - xmax, &tE, &tL, clipedges, CLIP_RIGHT))  // right
	if (CLIPt(dy, -ymax - *y0, &tE, &tL, clipedges, CLIP_BOTTOM))  // bottom
	  if (CLIPt(-dy, *y0 - ymax, &tE, &tL, clipedges, CLIP_TOP)){  // top
	    if (*clipedges){
	      *visible = 1;
	      if (tL < 1.f){
		*x1 = *x0 + tL * dx;
		*y1 = *y0 + tL * dy;
	      }
	      if (tE > 0.f){
		*x0 =+ tE * dx;
		*y0 += tE * dy;
	      }
	    }
	  }
  }
}

static
void Clip2DLine(float xmax, float ymax, float *line, short *visible, short *clipedges){
  Clip2D(xmax, ymax, &line[0], &line[1], &line[2], &line[3], visible, clipedges);
}

static
void glVertex3fTransformed(float *convMatrix, float x, float y, float z){
  float tmppt[3] = { x, y, z };
  MatrixTransformC44f3f(convMatrix, tmppt, tmppt);
  glVertex3fv(tmppt);
}

static
void drawLineToPointInWorldCrossClip(PyMOLGlobals *G, int label_z_target, float cw, float x1, float y1, float x2, float y2, float *cross, float *pt, float *curpt, float cx, float cy){
  float lvect[3];
  float nzn[3] = { 0.f, 0.f, 1.f };
  float convMatrix[16];
  short visible1, edges1, visible2, edges2;
  float line1[4], line2[4];

  if (!label_z_target){
    SceneGenerateMatrixToAnotherZFromZ(G, convMatrix, curpt, pt);
  } else {
    // set convMatrix to translation from label center to target center, i.e., (x2, y2)
    identity44f(convMatrix);
    MatrixTranslateC44f(convMatrix, x2, y2, 0.f);
  }

  lvect[0] = x2 - x1;
  lvect[1] = y2 - y1;
  normalize2f(lvect);
  cross_product3f(lvect, nzn, cross);
  mult3f(cross, cw, cross);

  line1[0] = x1 + cross[0]; line1[1] = y1 + cross[1];
  line1[2] = x2 + cross[0]; line1[3] = y2 + cross[1];
  line2[0] = x1 - cross[0]; line2[1] = y1 - cross[1];
  line2[2] = x2 - cross[0]; line2[3] = y2 - cross[1];
  Clip2DLine(cx, cy, line1, &visible1, &edges1);
  Clip2DLine(cx, cy, line2, &visible2, &edges2);
  
  if (visible1 && visible2){
    if (edges1 == edges2){
      // if both lines intersect the same edge, just draw the quad
      glBegin(GL_TRIANGLE_STRIP);
      glVertex3f(line1[2], line1[3], 0.f);
      glVertex3fTransformed(convMatrix, cross[0], cross[1], 0.f);
      glVertex3f(line2[2], line2[3], 0.f);
      glVertex3fTransformed(convMatrix, - cross[0], - cross[1], 0.f);
      glEnd();
    } else {
      // if both lines intersect different edges, need to add the corner
      // and generate correct order
      float corner[2];
      short edges = edges1 | edges2; 
      corner[0] = (edges & CLIP_LEFT) ? -cx : cx;
      corner[1] = (edges & CLIP_BOTTOM) ? -cy : cy;
      
      // staggered triangles which include the corner
      glBegin(GL_TRIANGLE_STRIP);
      glVertex3f(line1[2], line1[3], 0.f);
      glVertex3fTransformed(convMatrix, cross[0], cross[1], 0.f);
      glVertex3f(corner[0], corner[1], 0.f);
      glVertex3fTransformed(convMatrix, - cross[0], - cross[1], 0.f);
      glVertex3f(line2[2], line2[3], 0.f);
      glEnd();
    }
  }
}

static
void drawLine2DCheckZTargetCross(PyMOLGlobals *G, short label_z_target, float *pt, float *curpt, float cw, float x1, float y1, float x2, float y2, float *cross){
  if (label_z_target){
    drawLine2DCross(cw, x1, y1, x2, y2, cross);
  } else {
    drawLineToPointInWorldCross(G, cw, x1, y1, x2, y2, cross, pt, curpt);
  }
}

static
void drawLine2DCheckZTarget(PyMOLGlobals *G, short label_z_target, float *pt, float *curpt, float cw, float x1, float y1, float x2, float y2){
  float cross[3];
  drawLine2DCheckZTargetCross(G, label_z_target, pt, curpt, cw, x1, y1, x2, y2, cross);
}

static
void drawLine2DCheckZTargetClip(PyMOLGlobals *G, short label_z_target, float *pt, float *curpt, float cw, float x1, float y1, float x2, float y2, float cx, float cy){
  float cross[3];
  drawLineToPointInWorldCrossClip(G, label_z_target, cw, x1, y1, x2, y2, cross, pt, curpt, cx, cy);
}


void drawLineAsGeometryWithOffsets(float *pt1, float *pt2, float *spt1, float *spt2, float *xn, float *yn, float *zn, float line_width, float topext, float bottomext, float *dirv){
  float pt1E[3], pt2E[3];
  float tmpV[3], tmpV2[3], linev[3];
  float nzn[3] = { 0.f, 0.f, 1.f };
  copy3f(pt1, pt1E);
  copy3f(pt2, pt2E);
  subtract3f(spt1, spt2, tmpV);

  copy3f(tmpV, linev);
  normalize3f(linev);
  mult3f(linev, line_width, linev);
  
  glBegin(GL_TRIANGLE_STRIP);
  cross_product3f(tmpV, nzn, tmpV2);
  normalize3f(tmpV2);
  mult3f(tmpV2, line_width, dirv);
  addXYtoVertex(dirv[0], dirv[1], xn, yn, pt1, pt1E);
  addXYtoVertex(topext * linev[0], topext * linev[1], xn, yn, pt1E, pt1E);
  glVertex3fv(pt1E);

  addXYtoVertex(dirv[0], dirv[1], xn, yn, pt2, pt2E);
  addXYtoVertex(-topext * linev[0], -topext * linev[1], xn, yn, pt2E, pt2E);
  glVertex3fv(pt2E);

  addXYtoVertex(-dirv[0], -dirv[1], xn, yn, pt1, pt1E);
  addXYtoVertex(bottomext * linev[0], bottomext * linev[1], xn, yn, pt1E, pt1E);
  glVertex3fv(pt1E);

  addXYtoVertex(-dirv[0], -dirv[1], xn, yn, pt2, pt2E);
  addXYtoVertex(-bottomext * linev[0], -bottomext * linev[1], xn, yn, pt2E, pt2E);
  glVertex3fv(pt2E);

  glEnd();
}

static
void RepLabelRenderBackgroundInImmediate(PyMOLGlobals *G, RepLabel *I, float *v, int draw_var, float *tCenterPt, short relativeMode, float *xn, float *yn, 
					 float *PmvMatrix, float *RotMatrix, int screenwidth, int screenheight, float *screenWorldOffset, float *indentFactor, 
					 float text_width, float text_height, float font_size){
  float pos[3], *labelpos = TextGetLabelPushPos(G);
  float hwidth = text_width / 2.f, hheight = text_height / 2.f;
  short label_connector_mode = (draw_var & 8) ? 1 : (draw_var & 16) ? 2 : (draw_var & 32) ? 3 : (draw_var & 64) ? 4 : 0;
  float cw = *(v + 26) / 2.f;
  short label_z_target = relativeMode & 8;
  float indentFactorT[2] = { indentFactor[0]*text_width/2.f, indentFactor[1]*text_height/2.f };

  {// taking into account screen adjustments : indent and screen world offset to the label point
    float v_scale, xn[3], yn[3];
    v_scale = SceneGetScreenVertexScale(G, NULL);
    SceneGetScaledAxesAtPoint(G, labelpos, xn, yn);
    addXYtoVertex(indentFactorT[0] + screenWorldOffset[0]/v_scale, indentFactorT[1] + screenWorldOffset[1]/v_scale, xn, yn, labelpos, pos);
  }

  {
    float tCenter[4], tTarget[4], tVec[3];
    float doNotDraw;
    short drawLine = true;
    float dVectorInPixels[2], dVector[2];
    float v_scale;
    short visible;
    v_scale = SceneGetScreenVertexScale(G, NULL);

    copy3f(v+6, tCenter);
    tCenter[3] = 1.f;
    copy3f(v+3, tTarget);
    tTarget[3] = 1.f;
    MatrixTransformC44f4f(PmvMatrix, tTarget, tTarget);
    normalize4f(tTarget);
    if (!((relativeMode & 2) || (relativeMode & 4))){ // label_relative_mode = 0, not 1 or 2
      MatrixTransformC44f4f(PmvMatrix, tCenter, tCenter);
      normalize4f(tCenter);
    } // label_relative_mode =1 : nothing, lready in screen coordinates
    if (relativeMode & 4){ // label_relative_mode = 2
      tCenter[0] = (tCenter[0] / screenwidth) * 2 - 1.f;
      tCenter[1] = (tCenter[1] / screenheight) * 2 - 1.f;
    }
    
    subtract3f(tTarget, tCenter, tVec);
    dVectorInPixels[0] = .5f * (screenwidth * tVec[0]) - (screenWorldOffset[0]/v_scale) - indentFactorT[0];
    dVectorInPixels[1] = .5f * (screenheight * tVec[1]) - (screenWorldOffset[1]/v_scale) - indentFactorT[1];
    dVector[0] = (dVectorInPixels[0] / text_width);
    dVector[1] = (dVectorInPixels[1] / text_height);
    doNotDraw = ((fabs(dVector[0])) <= .5f && (fabs(dVector[1])) <= .5f) ? 1.f : 0.f;
    if (!(draw_var & 1) || doNotDraw > .5f)
      drawLine = false;

    if (label_z_target){
      visible = SceneGetVisible(G, v+3);
      ScenePushRasterMatrix(G, v + 3);
      glTranslatef(-dVectorInPixels[0], -dVectorInPixels[1], 0.f);
    } else {
      visible = SceneGetVisible(G, pos);
      ScenePushRasterMatrix(G, pos);
    }
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(0.f, 1.f);

    if (!visible){
      draw_var = 0;
      drawLine = 0;
    }
    if (draw_var & 6){ // draw background or background outline
      if (draw_var & 2){ // draw background
	float bgwidth = hwidth, bgheight = hheight;
	if (draw_var & 4){
	  bgwidth -= cw;
	  bgheight -= cw;
	}
	glColor4f(v[23], v[24], v[25], v[22]);
	glBegin(GL_TRIANGLE_STRIP);
	glVertex3f(bgwidth, bgheight, 0.f);
	glVertex3f(bgwidth, -bgheight, 0.f);
	glVertex3f(-bgwidth, bgheight, 0.f);
	glVertex3f(-bgwidth, -bgheight, 0.f);
	glEnd();
      }
      if (draw_var & 4){ // draw background outline
	glColor4f(v[9], v[10], v[11], 1.f);
	glBegin(GL_TRIANGLE_STRIP);
	glVertex3f(-hwidth - cw, -hheight - cw, 0.f); // 1
	glVertex3f(-hwidth - cw, hheight + cw, 0.f);  // 2    2---------------4
	glVertex3f(-hwidth + cw, hheight - cw, 0.f);  // 3    |\            / |
	glVertex3f(hwidth + cw, hheight + cw, 0.f);   // 4    | 3----------5  |
	glVertex3f(hwidth - cw, hheight - cw, 0.f);   // 5    | |          |  |
	glVertex3f(hwidth + cw, -hheight - cw, 0.f);  // 6    | 9----------7  |
	glVertex3f(hwidth - cw, -hheight + cw, 0.f);  // 7    |/            \ |
	glVertex3f(-hwidth - cw, -hheight - cw, 0.f); // 8  1/8---------------6
	glVertex3f(-hwidth + cw, -hheight + cw, 0.f); // 9
	glVertex3f(-hwidth + cw, hheight - cw, 0.f);  // 3
	glEnd();
      }
    }
    if (drawLine){
      float xoff, yoff;
      glColor4f(v[9], v[10], v[11], 1.f);
      switch (label_connector_mode){
      case 0:
	{
	  float hmid = (fabs(dVector[0]) >= .5) ? 0.f : 1.f;
	  float vmid = (fabs(dVector[1]) >= .5) ? 0.f : 1.f;
	  float right = (1.f - hmid) * ((dVector[0] > 0.f) ? 1.f : 0.f) + hmid * .5f;
	  float top = (1.f - vmid) * ((dVector[1] > 0.f) ? 1.f : 0.f) + vmid * .5f;
	  float labx, laby;
	  xoff = 2.f * (right - .5);
	  yoff = 2.f * (top - .5);
	  labx = xoff * hwidth;
	  laby = yoff * hheight;
	  
	  drawLine2DCheckZTarget(G, label_z_target, v + 3, pos, cw, labx, laby, dVectorInPixels[0], dVectorInPixels[1]);
	}
	break;
      case 3:
	{
	  float hmid = (fabs(dVector[0]) >= .5f) ? 0.f : 1.f;
	  float vmid = (fabs(dVector[1]) >= .5f) ? 0.f : 1.f;
	  float right = (1.f - hmid) * ((dVector[0] > 0.f) ? 1.f : 0.f) + hmid * (.5f + dVector[0]);
	  float top = (1.f - vmid) * ((dVector[1] > 0.f) ? 1.f : 0.f) + vmid * (.5f + dVector[1]);
	  float labx, laby;
	  xoff = 2.f * (right - .5);
	  yoff = 2.f * (top - .5);
	  labx = xoff * hwidth;
	  laby = yoff * hheight;
	  drawLine2DCheckZTarget(G, label_z_target, v + 3, pos, cw, labx, laby, dVectorInPixels[0], dVectorInPixels[1]);
	}
	break;
      case 1:
	{
	  float drawVectorN[2], absyx, notabsyx, dvxy, dvyx, hdir, vdir;
	  float labx, laby;
	  copy2f(dVector, drawVectorN);
	  normalize2f(drawVectorN);
	  absyx = fabs(drawVectorN[0]) >= fabs(drawVectorN[1]) ? 1.f : 0.f;
	  notabsyx = 1.f - absyx;
	  hdir = 2.f * ( ( (drawVectorN[0] > 0.f) ? 1.f : 0.f ) - .5f);
	  vdir = 2.f * ( ( (drawVectorN[1] > 0.f) ? 1.f : 0.f ) - .5f);
	  dvxy = dVector[0] / dVector[1];
	  dvyx = dVector[1] / dVector[0];
	  xoff = (absyx * hdir) + (notabsyx * vdir * dvxy);
	  yoff = (notabsyx * vdir) + (absyx * hdir * dvyx);
	  labx = xoff * hwidth;
	  laby = yoff * hheight;
	  if ((draw_var & 6) == 2){
	    // if ONLY background is drawn in connector_mode 1, then the line
	    // needs to be clipped so that it is flush to the background
	    drawLine2DCheckZTargetClip(G, label_z_target, v + 3, pos, cw, 0, 0, dVectorInPixels[0], dVectorInPixels[1], hwidth, hheight);
	  } else {
	    drawLine2DCheckZTarget(G, label_z_target, v + 3, pos, cw, labx, laby, dVectorInPixels[0], dVectorInPixels[1]);
	  }
	}
	break;
      case 2:
      case 4:
	{
	  float extLength, rightorig, right, top, xofforig;
	  short label_connector_mode_4 = (label_connector_mode == 4);
	  short hmid = 0, vmid = 0;
	  float labx, laby, extx, exty;
	  if (label_connector_mode_4){
	    hmid = (fabs(dVector[0]) < .5f);
	    vmid = (fabs(dVector[1]) < .5f);
	  }
	  if (font_size < 0.f){
	    extLength = - (*(v + 27) * font_size / text_height) / v_scale;
	  } else {
	    extLength = *(v + 27) * font_size / text_height;
	  }
	  if (hmid){
	    rightorig = (.5f + dVector[0]);
	    right = rightorig;
	  } else {
	    rightorig = dVector[0] >= 0.f ? 1.f : 0.f;
	    right = rightorig + (rightorig - .5f) * CLAMP_VALUE(fabs(2.f * dVector[0])-1.f, 0.f, 2.f * extLength * fabs(text_height/(float)text_width));
	  }
	  if (vmid){
	    top = (.5f + dVector[1]);
	  } else {
	    top = (dVector[1] >= 0.f) ? 1.f : 0.f;
	  }
	  xofforig = 2. * (rightorig - .5);
	  
	  xoff = 2.f * (right - .5);
	  yoff = 2.f * (top - .5);

	  labx = xofforig * hwidth;
	  laby = yoff * hheight;
	  extx = xoff * hwidth;
	  exty = yoff * hheight;

	  { // drawing the triangle filler between the two lines (connector_mode 2 or 4)
	    float dirv1[3], dirv2[3], cross[3];
	    drawLine2DCheckZTargetCross(G, label_z_target, v + 3, pos, cw, extx, exty, dVectorInPixels[0], dVectorInPixels[1], dirv1);
	    drawLine2DCross(cw, extx, exty, labx, laby, dirv2);
	    dirv1[0] = -dirv1[0];  dirv1[1] = -dirv1[1];  // switch normal b/c drawLine2DCheckZTargetCross call above has target as second
	    dirv1[2] = dirv2[2] = 0.f;
	    cross_product3f(dirv1, dirv2, cross);
	    glBegin(GL_TRIANGLES);
	    glVertex3f(extx, exty, 0.f);
	    if (cross[2] < 0.f){
	      glVertex3f(extx - dirv1[0], exty - dirv1[1], 0.f);
	      glVertex3f(extx - dirv2[0], exty - dirv2[1], 0.f);
	    } else {
	      glVertex3f(extx + dirv2[0], exty + dirv2[1], 0.f);
	      glVertex3f(extx + dirv1[0], exty + dirv1[1], 0.f);
	    }
	    glEnd();
	  }
	}
	break;
      }
    }
  }
  glDisable(GL_POLYGON_OFFSET_FILL);
  ScenePopRasterMatrix(G);
}
#endif

static int lineSegIntersection(float *p1, float *p2, float *p3, float *p4, float *po) {
	float  distAB, theCos, theSin, newX, abpos;
	float p2i[2], p3i[2], p4i[2];
	if ((p1[0]==p2[0] && p1[1]==p2[1]) || (p3[0]==p4[0] && p3[1]==p4[1])) return 0;
	if ((p1[0]==p3[0] && p1[1]==p3[1]) || (p2[0]==p3[0] && p2[1]==p3[1])
	||  (p1[0]==p4[0] && p1[1]==p4[1]) || (p2[0]==p4[0] && p2[1]==p4[1])) {
	return 0; }
	p2i[0] = p2[0]; p2i[1] = p2[1];
	p3i[0] = p3[0]; p3i[1] = p3[1];
	p4i[0] = p4[0]; p4i[1] = p4[1];
	p2i[0] -= p1[0]; p2i[1] -= p1[1];
	p3i[0] -= p1[0]; p3i[1] -= p1[1];
	p4i[0] -= p1[0]; p4i[1] -= p1[1];
	distAB = sqrt(p2i[0] * p2i[0] + p2i[1] * p2i[1]);
	theCos = p2i[0] / distAB;
	theSin = p2i[1] / distAB;
	newX = p3i[0] * theCos + p3i[1] * theSin;
	p3i[1] = p3i[1] * theCos - p3i[0] * theSin;
	p3i[0] = newX;
	newX = p4i[0] * theCos + p4i[1] * theSin;
	p4i[1] = p4i[1] * theCos - p4i[0] * theSin;
	p4i[0] = newX;
	if ((p3i[1] < 0. && p4i[1]<0.) || (p3i[1] >= 0. && p4i[1] >= 0.)) return 0;
	abpos = p4i[0] + (p3i[0] - p4i[0]) * p4i[1] / (p4i[1] - p3i[1]);
	if (abpos<0. || abpos>distAB) return 0;
	po[0] = p1[0] + abpos * theCos;
	po[1] = p1[1] + abpos * theSin;
	return 1; 
}

typedef struct lineSeg_s {
  float p1[2];
  float p2[2];
} lineSeg_t;


static void RepLabelRenderRayBackground(RepLabel * I, RenderInfo * info, float *v, int draw_var){
  CRay *ray = info->ray;
  PyMOLGlobals *G = I->G;
  float *screenWorldOffset = TextGetScreenWorldOffset(G);
  float text_width = TextGetWidth(G), text_height = TextGetHeight(G);
  float *indentFactor = TextGetIndentFactor(G);
  float endpointOnBBX[3], tmp3f[3];
  float torigCenter[3], tTarget[4], tVec[3], dVectorInPixels[2], dVector[2];
  float xoff = 0.f, yoff = 0.f;
  float doNotDraw;
  float v_scale;
  short drawLine = true;
  float xn[3], yn[3], zn[3];
  float connector_width = *(v + 26);
  float *RotMatrix = ray->Rotation;
  lineSeg_t labelTop, labelBottom, labelLeft, labelRight;
  short label_con_flat = 128 & (int)*(v + 21);
  short label_connector_mode = (draw_var & 8) ? 1 : (draw_var & 16) ? 2 : (draw_var & 32) ? 3 : (draw_var & 64) ? 4 : 0;
  float font_size = SettingGet_f(G, I->cs->Setting.get(), I->obj->Setting.get(),
                                 cSetting_label_size);
  float tCenter[4], sCenter[4], sTarget[4];
  short relativeMode = ((short)*(v + 15));
  copy3f(TextGetLabelPushPos(G), tCenter);
  copy3f(tCenter, torigCenter);
  TextSetPosNColor(G, tCenter, v);
  RayGetScaledAllAxesAtPoint(ray, tCenter, xn, yn, zn);
  copy3f(v + 3, tTarget);
  tCenter[3] = 1.f;
  tTarget[3] = 1.f;

  v_scale = RayGetScreenVertexScale(ray,torigCenter);
  RayGetScreenVertex(ray, tCenter, sCenter);
  if (relativeMode & 8){  // label_z_target, adjust z to target
    RayGetScreenVertex(ray, tTarget, sTarget);
  } else {
    // need to bring tTarget into same z as tCenter to measure
    copy3f(tTarget, sTarget);
    sTarget[3] = 1.f;
    RayAdjustZtoScreenZofPoint(ray, sTarget, tCenter);
    RayGetScreenVertex(ray, sTarget, sTarget);
  }
  subtract3f(sTarget, sCenter, tVec);
  dVectorInPixels[0] = (tVec[0]-screenWorldOffset[0])/v_scale;
  dVectorInPixels[1] = (tVec[1]-screenWorldOffset[1])/v_scale;
  dVector[0] = (dVectorInPixels[0] / text_width) - indentFactor[0];
  dVector[1] = (dVectorInPixels[1] / text_height) - indentFactor[1];
  doNotDraw = ((fabs(dVector[0])) <= .5f && (fabs(dVector[1])) <= .5f) ? 1.f : 0.f;
  if (!(draw_var & 1) || doNotDraw > .5f)
    drawLine = false;

  mult3f(xn, indentFactor[0] * text_width + 2.f * (screenWorldOffset[0] / v_scale), tmp3f);
  add3f(tmp3f, torigCenter, torigCenter);
  mult3f(yn, indentFactor[1] * text_height + 2.f * (screenWorldOffset[1] / v_scale), tmp3f);
  add3f(tmp3f, torigCenter, torigCenter);

  mult3f(zn, -1.f, tmp3f);  // push background back, not sure if this is needed, doesn't hurt
  add3f(tmp3f, torigCenter, torigCenter);

  if (draw_var & 6){ // draw background or background outline
    float tmpf[4][4] ;
    float hwidth = ray->Sampling * text_width / 2.f, hheight = ray->Sampling * text_height / 2.f;
    addXYtoVertex(hwidth, hheight, xn, yn, torigCenter, tmpf[0]);
    addXYtoVertex(hwidth, -hheight, xn, yn, torigCenter, tmpf[1]);
    addXYtoVertex(-hwidth, hheight, xn, yn, torigCenter, tmpf[2]);
    addXYtoVertex(-hwidth, -hheight, xn, yn, torigCenter, tmpf[3]);
    if (draw_var & 2){ // draw background
      float trans = 1.f - v[22]; 
      float tmpf2[4][4];
      float tmpc[4][4];
      float lw2;
      if (draw_var & 4) {
	lw2 = ray->Sampling * connector_width / 2.f; 
      } else {
	lw2 = ray->Sampling;
      }
      addXYZtoVertex(-lw2, -lw2, -1.f, xn, yn, zn, tmpf[0], tmpf2[0]); // UL
      addXYZtoVertex(-lw2, lw2, -1.f, xn, yn, zn, tmpf[1], tmpf2[1]); // LL
      addXYZtoVertex(lw2, -lw2, -1.f, xn, yn, zn, tmpf[2], tmpf2[2]); // UR
      addXYZtoVertex(lw2, lw2, -1.f, xn, yn, zn, tmpf[3], tmpf2[3]); // LR
      MatrixTransformC44f4f(RotMatrix, tmpf[0], tmpc[0]);
      MatrixTransformC44f4f(RotMatrix, tmpf[1], tmpc[1]);
      MatrixTransformC44f4f(RotMatrix, tmpf[2], tmpc[2]);
      MatrixTransformC44f4f(RotMatrix, tmpf[3], tmpc[3]);
      tmpc[0][0] = tmpc[1][0];
      tmpc[2][0] = tmpc[3][0];
      tmpc[0][1] = tmpc[2][1];
      tmpc[1][1] = tmpc[3][1];

      copy2f(tmpc[0], labelTop.p1);
      copy2f(tmpc[2], labelTop.p2);

      copy2f(tmpc[3], labelBottom.p1);
      copy2f(tmpc[1], labelBottom.p2);

      copy2f(tmpc[2], labelLeft.p1);
      copy2f(tmpc[3], labelLeft.p2);

      copy2f(tmpc[1], labelRight.p1);
      copy2f(tmpc[0], labelRight.p2);
      ray->triangleTrans3fv(tmpf2[0], tmpf2[1], tmpf2[2], zn, zn, zn, &v[23], &v[23], &v[23], trans, trans, trans);
      ray->setLastToNoLighting(1);
      ray->triangleTrans3fv(tmpf2[1], tmpf2[2], tmpf2[3], zn, zn, zn, &v[23], &v[23], &v[23], trans, trans, trans);
      ray->setLastToNoLighting(1);
    }
    if (draw_var & 4){ // draw background outline
      if (label_con_flat){
	float tmpfs[4][4] ;
	float lw = ray->Sampling * connector_width/2.f, dirv[3];
	tmpf[0][3] = tmpf[1][3] = tmpf[2][3] = tmpf[3][3] = 0.f;
	MatrixTransformC44f4f(ray->ModelView, tmpf[0], tmpfs[0]);
	MatrixTransformC44f4f(ray->ModelView, tmpf[1], tmpfs[1]);
	MatrixTransformC44f4f(ray->ModelView, tmpf[2], tmpfs[2]);
	MatrixTransformC44f4f(ray->ModelView, tmpf[3], tmpfs[3]);
	RayDrawLineAsGeometryWithOffsets(ray, tmpf[2], tmpf[0], tmpfs[2], tmpfs[0], xn, yn, zn, lw, 1.f, 0.f, v + 9, dirv, 1);
	RayDrawLineAsGeometryWithOffsets(ray, tmpf[1], tmpf[0], tmpfs[1], tmpfs[0], xn, yn, zn, lw, 0.f, 1.f, v + 9, dirv, 1);
	RayDrawLineAsGeometryWithOffsets(ray, tmpf[3], tmpf[1], tmpfs[3], tmpfs[1], xn, yn, zn, lw, 0.f, 1.f, v + 9, dirv, 1);
	RayDrawLineAsGeometryWithOffsets(ray, tmpf[3], tmpf[2], tmpfs[3], tmpfs[2], xn, yn, zn, lw, 1.f, 0.f, v + 9, dirv, 1);
      } else {
	float lw = connector_width * ray->PixelRadius / 2.f;
	ray->sausage3fv(tmpf[0], tmpf[1], lw, v + 9, v + 9);
	ray->sausage3fv(tmpf[0], tmpf[2], lw, v + 9, v + 9);
	ray->sausage3fv(tmpf[1], tmpf[3], lw, v + 9, v + 9);
	ray->sausage3fv(tmpf[2], tmpf[3], lw, v + 9, v + 9);
      }
    }
  }
  switch (label_connector_mode){
  case 0:
    {
      float hmid = (fabs(dVector[0]) >= .5) ? 0.f : 1.f;
      float vmid = (fabs(dVector[1]) >= .5) ? 0.f : 1.f;
      float right = (1.f - hmid) * ((dVector[0] > 0.f) ? 1.f : 0.f) + hmid * .5f;
      float top = (1.f - vmid) * ((dVector[1] > 0.f) ? 1.f : 0.f) + vmid * .5f;
      xoff = 2.f * (right - .5);
      yoff = 2.f * (top - .5);
    }
    break;
  case 3:
    {
      short hmid = (fabs(dVector[0]) < 1.f);
      short vmid = (fabs(dVector[1]) < 1.f);
      float right, top;
      if (hmid){
	right = ((1.f + dVector[0]) / 2.f);
      } else {
	right = ((dVector[0] > 0.f) ? 1.f : 0.f);
      }
      if (vmid){
	top = ((1.f + dVector[1]) / 2.f);		  
      } else {
	top = ((dVector[1] > 0.f) ? 1.f : 0.f);
      }
      xoff = 2.f * (right - .5);
      yoff = 2.f * (top - .5);
    }
    break;
  case 1:
    {
      float drawVectorN[2], absyx, notabsyx, dvxy, dvyx, hdir, vdir;
      copy2f(dVector, drawVectorN);
      normalize2f(drawVectorN);
      absyx = fabs(drawVectorN[0]) >= fabs(drawVectorN[1]) ? 1.f : 0.f;
      notabsyx = 1.f - absyx;
      hdir = 2.f * ( ( (drawVectorN[0] > 0.f) ? 1.f : 0.f ) - .5f);
      vdir = 2.f * ( ( (drawVectorN[1] > 0.f) ? 1.f : 0.f ) - .5f);
      dvxy = dVector[0] / dVector[1];
      dvyx = dVector[1] / dVector[0];
      xoff = (absyx * hdir) + (notabsyx * vdir * dvxy);
      yoff = (notabsyx * vdir) + (absyx * hdir * dvyx);
    }
    break;
  case 2:
  case 4:
    if (drawLine){
      float extLength, rightorig, right, top, xofforig;
      float endpointExtendedOffBBX[3];
      short label_connector_mode_4 = (label_connector_mode == 4);
      short hmid = 0, vmid = 0;
      if (label_connector_mode_4){
	hmid = (fabs(dVector[0]) < 1.);
	vmid = (fabs(dVector[1]) < 1.);
      }
      if (font_size < 0.f){
	extLength = - (*(v + 27) * font_size / text_height) / v_scale;
      } else {
	extLength = *(v + 27) * font_size / text_height;
      }
      if (hmid){
	rightorig = ((1.f + dVector[0])/2.f);
	right = rightorig;
      } else {
	rightorig = dVector[0] >= 0.f ? 1.f : 0.f;
	right = rightorig + (rightorig - .5f) * CLAMP_VALUE(fabs(dVector[0]*2.f)-1.f, 0.f,
							    2.f * extLength * fabs(text_height/(float)text_width));
      }
      if (vmid){
	top = ((1.f + dVector[1])/2.f);
      } else {
	top = dVector[1] >= 0.f ? 1.f : 0.f;
      }
      xofforig = 2. * (rightorig - .5);
      
      xoff = 2.f * (right - .5);
      yoff = 2.f * (top - .5);
      
      mult3f(xn, xofforig * text_width, endpointOnBBX);
      mult3f(yn, yoff * text_height, tmp3f);
      add3f(tmp3f, endpointOnBBX, endpointOnBBX);
      add3f(torigCenter, endpointOnBBX, endpointOnBBX);
      
      mult3f(xn, xoff * text_width, endpointExtendedOffBBX);
      mult3f(yn, yoff * text_height, tmp3f);
      add3f(tmp3f, endpointExtendedOffBBX, endpointExtendedOffBBX);
      add3f(torigCenter, endpointExtendedOffBBX, endpointExtendedOffBBX);
      if (label_con_flat){
	float tmpf[3][4] ;
	float tmpfs[3][4] ;
	float dirv1[3], dirv2[3], cross[3], pt[3][3];
	copy3f(v + 3, tmpf[0]);
	copy3f(endpointExtendedOffBBX, tmpf[1]);
	copy3f(endpointOnBBX, tmpf[2]);
	tmpf[0][3] = tmpf[1][3] = tmpf[2][3] = 0.f;
	MatrixTransformC44f4f(ray->ModelView, tmpf[0], tmpfs[0]);
	MatrixTransformC44f4f(ray->ModelView, tmpf[1], tmpfs[1]);
	MatrixTransformC44f4f(ray->ModelView, tmpf[2], tmpfs[2]);
	RayDrawLineAsGeometryWithOffsets(ray, tmpf[0], tmpf[1], tmpfs[0], tmpfs[1], xn, yn, zn, connector_width, 0.f, 0.f, v + 9, dirv1, 1);
	RayDrawLineAsGeometryWithOffsets(ray, tmpf[1], tmpf[2], tmpfs[1], tmpfs[2], xn, yn, zn, connector_width, 0.f, 0.f, v + 9, dirv2, 1);
	/* need to draw triangle to fill in gap */
	dirv1[2] = dirv2[2] = 0.f;
	cross_product3f(dirv1, dirv2, cross);
	copy3f(endpointExtendedOffBBX, pt[0]);
	copy3f(endpointExtendedOffBBX, pt[1]);
	copy3f(endpointExtendedOffBBX, pt[2]);
	if (cross[2] < 0.f){
	  addXYtoVertex(dirv1[0], dirv1[1], xn, yn, pt[1], pt[1]);
	  addXYtoVertex(dirv2[0], dirv2[1], xn, yn, pt[2], pt[2]);
	} else {
	  addXYtoVertex(-dirv2[0], -dirv2[1], xn, yn, pt[1], pt[1]);
	  addXYtoVertex(-dirv1[0], -dirv1[1], xn, yn, pt[2], pt[2]);
	}
	ray->triangle3fv(pt[0], pt[1], pt[2], zn, zn, zn, v + 9, v + 9, v + 9);		  
	ray->setLastToNoLighting(1);
      } else {
	float lw = connector_width * ray->PixelRadius / 2.f;
	ray->sausage3fv(v + 3, endpointExtendedOffBBX, lw, v + 9, v + 9);
	ray->sausage3fv(endpointExtendedOffBBX, endpointOnBBX, lw, v + 9, v + 9);
      }
      drawLine = 0;
    }
    break;
  }
  
  if (drawLine){
    mult3f(xn, xoff * text_width, endpointOnBBX);
    mult3f(yn, yoff * text_height, tmp3f);
    add3f(tmp3f, endpointOnBBX, endpointOnBBX);
    add3f(torigCenter, endpointOnBBX, endpointOnBBX);
    
    if (label_con_flat){
      unsigned char drawDefaultLine = true;
      if ( (draw_var & 2) && (label_connector_mode > 0) ){ // draw background
	float tmpff[3][4];
	float tmpfs[3][4];
	float tmpfc[4][4];
	float dirv1[3], dirv2[3];
	float cw = connector_width;
	float pt1E[3], pt2E[3];
	float tmpV[3], tmpV2[3], tmpV3[3], tmpV4[3], linev2[3];
	float nzn[3] = { 0.f, 0.f, 1.f };
	lineSeg_t l1, l2;
	float tmptrn[4][4];
	float isec1[3], isec2[3];
	int l1V = 0, l2V = 0;

	copy3f(v + 3, tmpff[0]);
	copy3f(endpointOnBBX, tmpff[1]);
	copy3f(torigCenter, tmpff[2]);

	tmpff[0][3] = tmpff[1][3] = tmpff[2][3] = 0.f;
	MatrixTransformC44f4f(ray->ModelView, tmpff[0], tmpfs[0]);
	MatrixTransformC44f4f(ray->ModelView, tmpff[1], tmpfs[1]);
	MatrixTransformC44f4f(ray->ModelView, tmpff[2], tmpfs[2]);

	copy3f(tmpff[0], pt1E);
	copy3f(tmpff[1], pt2E);

	subtract3f(tmpfs[0], tmpfs[1], tmpV);
	subtract3f(tmpfs[1], tmpfs[2], tmpV2);

	copy3f(tmpV2, linev2);
	normalize3f(linev2);

	mult3f(linev2, cw * 10.0f, linev2);

	cross_product3f(tmpV, nzn, tmpV3);
	cross_product3f(tmpV2, nzn, tmpV4);

	normalize3f(tmpV3);
	normalize3f(tmpV4);

	mult3f(tmpV3, cw, dirv1);
	mult3f(tmpV4, cw, dirv2);

	addXYtoVertex(dirv2[0], dirv2[1], xn, yn, tmpff[1], tmptrn[0]);
	addXYtoVertex(linev2[0], linev2[1], xn, yn, tmptrn[0], tmptrn[0]);
	addXYtoVertex(-dirv2[0], -dirv2[1], xn, yn, tmpff[1], tmptrn[1]);
	addXYtoVertex(linev2[0], linev2[1], xn, yn, tmptrn[1], tmptrn[1]);
	addXYtoVertex(dirv2[0], dirv2[1], xn, yn, tmpff[2], tmptrn[2]);
	addXYtoVertex(-dirv2[0], -dirv2[1], xn, yn, tmpff[2], tmptrn[3]);

	MatrixTransformC44f4f(RotMatrix, tmptrn[0], tmpfc[0]);
	MatrixTransformC44f4f(RotMatrix, tmptrn[1], tmpfc[1]);
	MatrixTransformC44f4f(RotMatrix, tmptrn[2], tmpfc[2]);
	MatrixTransformC44f4f(RotMatrix, tmptrn[3], tmpfc[3]);
	
	copy2f(tmpfc[0], l1.p1);
	copy2f(tmpfc[1], l2.p1);
	copy2f(tmpfc[2], l1.p2);
	copy2f(tmpfc[3], l2.p2);
			
	isec1[0] = isec1[1] = isec2[0] = isec2[1] = 0.f;
	isec1[2] = isec2[2] = tmpfc[3][2];

	if (lineSegIntersection(labelTop.p1, labelTop.p2, l1.p1, l1.p2, isec1)) {
	  l1V = 1;
	} else if (lineSegIntersection(labelRight.p1, labelRight.p2, l1.p1, l1.p2, isec1)) {
	  l1V = 2;
	} else if (lineSegIntersection(labelBottom.p1, labelBottom.p2, l1.p1, l1.p2, isec1)) {
	  l1V = 3;
	} else if (lineSegIntersection(labelLeft.p1, labelLeft.p2, l1.p1, l1.p2, isec1)) {
	  l1V = 4;
	} 

	if (lineSegIntersection(labelTop.p1, labelTop.p2, l2.p1, l2.p2, isec2)) {
	  l2V = 1;
	} else if (lineSegIntersection(labelRight.p1, labelRight.p2, l2.p1, l2.p2, isec2)) {
	  l2V = 2;
	} else if (lineSegIntersection(labelBottom.p1, labelBottom.p2, l2.p1, l2.p2, isec2)) {
	  l2V = 3;
	} else if (lineSegIntersection(labelLeft.p1, labelLeft.p2, l2.p1, l2.p2, isec2)) {
	  l2V = 4;
	} 

	// Both lines run through one side
	if ( (l1V && l2V) ){
	  drawDefaultLine = false;
	  if ( l1V == l2V ) {
	    float isec1w[3], isec2w[3];
	    MatrixInvTransformC44fAs33f3f(RotMatrix, isec1, isec1w);
	    MatrixInvTransformC44fAs33f3f(RotMatrix, isec2, isec2w);
	    
	    addXYtoVertex(dirv1[0], dirv1[1], xn, yn, tmpff[0], pt1E);
	    addXYtoVertex(-dirv1[0], -dirv1[1], xn, yn, tmpff[0], pt2E);
	    
	    ray->triangle3fv(pt1E, pt2E, isec1w, zn, zn, zn, v + 9, v + 9, v + 9);
	    ray->triangle3fv(pt2E, isec1w, isec2w, zn, zn, zn, v + 9, v + 9, v + 9);
	  } else { // (l1V != l2V)  
	    float isecc[3];
	    float isec1w[3], isec2w[3], iseccw[3];
	    isecc[2] = isec1[2];
	    // Find the corner 
	    if ( (l1V == 1 && l2V == 2) || (l1V == 2 && l2V == 1) ) { // UR
	      isecc[0] = labelTop.p1[0];
	      isecc[1] = labelTop.p1[1];
	    } else if ( (l1V == 2 && l2V == 3) || (l1V == 3 && l2V == 2) ) { // LR
	      isecc[0] = labelBottom.p2[0];
	      isecc[1] = labelBottom.p2[1];
	    } else if ( (l1V == 3 && l2V == 4) || (l1V == 4 && l2V == 3) ) { // LL
	      isecc[0] = labelBottom.p1[0];
	      isecc[1] = labelBottom.p1[1];
	    } else if ( (l1V == 4 && l2V == 1) || (l1V == 1 && l2V == 4) ) { // UL
	      isecc[0] = labelTop.p2[0];
	      isecc[1] = labelTop.p2[1];
	    }
	    
	    MatrixInvTransformC44fAs33f3f(RotMatrix, isec1, isec1w);
	    MatrixInvTransformC44fAs33f3f(RotMatrix, isec2, isec2w);
	    MatrixInvTransformC44fAs33f3f(RotMatrix, isecc, iseccw);
	    
	    addXYtoVertex(dirv1[0], dirv1[1], xn, yn, tmpff[0], pt1E);
	    addXYtoVertex(-dirv1[0], -dirv1[1], xn, yn, tmpff[0], pt2E);
	    
	    ray->triangle3fv(isec1w, pt1E, iseccw, zn, zn, zn, v + 9, v + 9, v + 9);
	    ray->triangle3fv(pt1E, iseccw, pt2E, zn, zn, zn, v + 9, v + 9, v + 9);
	    ray->triangle3fv(iseccw, pt2E, isec2w, zn, zn, zn, v + 9, v + 9, v + 9);
	  }
	}
      }
      if (drawDefaultLine){
	float tmpf[2][4] ;
	float tmpfs[2][4], dirv[3] ;
	copy3f(v + 3, tmpf[0]);
	copy3f(endpointOnBBX, tmpf[1]);
	tmpf[0][3] = tmpf[1][3] = 0.f;
	MatrixTransformC44f4f(ray->ModelView, tmpf[0], tmpfs[0]);
	MatrixTransformC44f4f(ray->ModelView, tmpf[1], tmpfs[1]);
	RayDrawLineAsGeometryWithOffsets(ray, tmpf[0], tmpf[1], tmpfs[0], tmpfs[1], xn, yn, zn, connector_width, 0.f, 0.f, v + 9, dirv, 1);
      }
    } else {
      ray->sausage3fv(v + 3, endpointOnBBX, connector_width * ray->PixelRadius / 2.f, v + 9, v + 9);
    }
  }
}

static
void RepLabelRenderRay(RepLabel * I, RenderInfo * info){
#ifndef _PYMOL_NO_RAY
  PyMOLGlobals *G = I->G;
  CRay *ray = info->ray;
  int c = I->N;
  float *v = I->V;
  lexidx_t *l = I->L;
  int font_id = SettingGet_i(G, I->cs->Setting.get(), I->obj->Setting.get(),
                             cSetting_label_font_id);
  float font_size = SettingGet_f(G, I->cs->Setting.get(), I->obj->Setting.get(),
                                 cSetting_label_size);
  if(c) {
    const char *st;
    TextSetOutlineColor(G, I->OutlineColor);
    while(c--) {
      if(*l) {
        float xn[3], yn[3], tCenter[3], offpt[3];
        short relativeMode = ((short)*(v + 15));
        int draw_var = 127 & (int)*(v + 21);

        copy3f(v + 6, tCenter);
        SceneGetCenter(G, offpt);
        RayGetScaledAxes(ray, xn, yn);

        st = LexStr(G, *l);
        TextSetLabelBkgrdInfo(G, *(v + 16), *(v + 17), (v + 18));
        if (relativeMode & 8){  // label_z_target, adjust z to target
          TextGetLabelPos(G)[0] = (SceneGetDepth(G, v+3) - .5) * 2.f;
          TextSetLabelPosIsSet(G, 1);
        } else if (relativeMode & 6){ // label_relative_mode 1 or 2, i.e., screen stabilized, adjust z 
          if (relativeMode & 4){ // label_relative_mode = 2
            tCenter[0] = (tCenter[0] / ray->Width) * 2.f - 1.f;
            tCenter[1] = (tCenter[1] / ray->Height) * 2.f - 1.f;
          }

          TextSetLabelPos(G, tCenter);
          TextSetLabelPosIsSet(G, 2);

          float tmp3f[3];
          mult3f(xn, tCenter[0] * ray->Width, tmp3f);
          add3f(tmp3f, offpt, offpt);
          mult3f(yn, tCenter[1] * ray->Height, tmp3f);
          add3f(tmp3f, offpt, offpt);
          copy3f(offpt, tCenter);
        } else {
          TextSetLabelPosIsSet(G, 0);
        }

        TextSetPosNColor(G, tCenter, v);

        TextRenderRay(G, ray, font_id, st, font_size, v + 12, (draw_var ? 1 : 0), ((short)*(v + 15)));
        if (draw_var){
          RepLabelRenderRayBackground(I, info, v, draw_var);
        }
      }
      v += 28;
      l++;
    }
  }
#endif
}

void RepLabel::render(RenderInfo* info)
{
  auto I = this;
  CRay *ray = info->ray;
  auto pick = info->pick;
  float *v = I->V;
  int c = I->N;
  lexidx_t *l = I->L;
  int font_id = SettingGet_i(G, I->cs->Setting.get(), I->obj->Setting.get(),
                             cSetting_label_font_id);
  float font_size = SettingGet_f(G, I->cs->Setting.get(), I->obj->Setting.get(),
                                 cSetting_label_size);
  int float_text = SettingGet_i(G, I->cs->Setting.get(), I->obj->Setting.get(),
				cSetting_float_labels);
  if (!(ray || pick) && info->pass != RenderPass::Transparent)
    return;

  if(I->MaxInvalid >= cRepInvRep){
    return;
  }
  font_id = SettingCheckFontID(G, I->cs->Setting.get(), I->obj->Setting.get(), font_id);

  if (I->shaderCGO && font_size < 0.f){
    int size;
    if (InvalidateShaderCGOIfTextureNeedsUpdate(G, font_size, I->texture_font_size, &size)){
      CGOFree(I->shaderCGO);
      I->texture_font_size = size;
    }
  }
  if(ray) {
    RepLabelRenderRay(I, info);
  } else if(G->HaveGUI && G->ValidContext) {
    if(pick) {
      int pick_labels = SettingGet_b(G, I->cs->Setting.get(), I->obj->Setting.get(), cSetting_pick_labels);
      if (!pick_labels)
	return;
      if (I->shaderCGO){
        if(float_text)
          glDisable(GL_DEPTH_TEST);
	CGORenderGLPicking(I->shaderCGO, info, &I->context, I->cs->Setting.get(), I->obj->Setting.get());
        if(float_text)
          glEnable(GL_DEPTH_TEST);
	return;
      } else {
        Pickable *p = I->P;
        TextSetIsPicking(G, true);
        SceneSetupGLPicking(G);
        if(c) {
          const char *st;
          int screenwidth, screenheight;
          if(float_text)
            glDisable(GL_DEPTH_TEST);
          
          if (!I->shaderCGO){
            SceneGetWidthHeight(G, &screenwidth, &screenheight);
          }

          while(c--) {
            if(*l) {
              float xn[3], yn[3], tCenterPt[3], offpt[3];
              short relativeMode = ((short)*(v + 15));
              copy3f(v + 6, tCenterPt);
              SceneGetCenter(G, offpt);
              TextSetPosNColor(G, offpt, v);
              SceneGetScaledAxes(G, I->obj, xn, yn);
              if (!I->shaderCGO){
                if (relativeMode & 2){ // label_relative_mode = 1
                  float tmp3f[3];
                  mult3f(xn, tCenterPt[0] * screenwidth/2.f, tmp3f);
                  add3f(tmp3f, offpt, offpt);
                  mult3f(yn, tCenterPt[1] * screenheight/2.f, tmp3f);
                  add3f(tmp3f, offpt, offpt);
                  copy3f(offpt, tCenterPt);
                } else if (relativeMode & 4){ // label_relative_mode = 2
                  float tmp3f[3];
                  mult3f(xn, (tCenterPt[0]  - (.5f * screenwidth)), tmp3f);
                  add3f(tmp3f, offpt, offpt);
                  mult3f(yn, (tCenterPt[1]  - (.5f * screenheight)), tmp3f);
                  add3f(tmp3f, offpt, offpt);
                  copy3f(offpt, tCenterPt);
                }
              }
              TextSetPosNColor(G, tCenterPt, v);
              TextSetTargetPos(G, v + 3);
              TextSetLabelBkgrdInfo(G, *(v + 16), *(v + 17), (v + 18));

              if (p) {
                p++;
                AssignNewPickColor(nullptr, pick, TextGetColorUChar4uv(G),
                    &I->context, p->index, p->bond);
              }

              TextSetColorFromUColor(G);

              st = LexStr(G, *l);
              if (!TextRenderOpenGL(G, info, font_id, st, font_size, v + 12, false, (short)*(v + 15), 1, SHADERCGO)){
                TextSetIsPicking(G, false);
                return ;
              }
            }
            l++;
            v += 28;
          }
          if(float_text)
            glEnable(GL_DEPTH_TEST);
        }
        TextSetIsPicking(G, false);
      }
    } else {  // not pick or ray, render
      if(c) {
        const char *st;
	short use_shader, has_connector = 0;
	CGO *connectorCGO = NULL;
	float *PmvMatrix = NULL;
	int screenwidth, screenheight;
	float xn[3] = { 1.0F, 0.0F, 0.0F };
	float yn[3] = { 0.0F, 1.0F, 0.0F };
	int pre_use_shaders = info->use_shaders;
	
	Pickable *p = I->P;
	use_shader = SettingGetGlobal_b(G, cSetting_use_shaders)
#ifdef _PYMOL_IOS
          ;
#else
	  && G->ShaderMgr->GeometryShadersPresent();
#endif
	info->use_shaders = use_shader;
	if (use_shader){
	  if (!I->shaderCGO){
	    I->shaderCGO = CGONew(G);
            I->shaderCGO->use_shader = true;
	  } else {
	    info->texture_font_size = I->texture_font_size;
	    if(float_text)
	      glDisable(GL_DEPTH_TEST);
	    CGORenderGL(I->shaderCGO, NULL, NULL, NULL, info, I);
	    if(float_text)
	      glEnable(GL_DEPTH_TEST);
	    return;
	  }
	} else {
          CGOFree(I->shaderCGO);

#ifndef PURE_OPENGL_ES_2
	  if(!info->line_lighting)
	    glDisable(GL_LIGHTING);
#endif
	}
        TextSetOutlineColor(G, I->OutlineColor);
	if (I->shaderCGO && c){
	  connectorCGO = CGONew(G);
	  CGOBegin(connectorCGO, GL_LINES);
	}
	if (!I->shaderCGO){
	  PmvMatrix = SceneGetPmvMatrix(G);
	  SceneGetWidthHeight(G, &screenwidth, &screenheight);
	  MatrixInvTransformC44fAs33f3f(PmvMatrix, xn, xn);
	  MatrixInvTransformC44fAs33f3f(PmvMatrix, yn, yn);
	  normalize3f(xn);
	  normalize3f(yn);
	}
        while(c--) {
          if(*l) {
	    float tCenterPt[3], offpt[3];
	    short relativeMode = ((short)*(v + 15));
	    int draw_var = 127 & (int)*(v + 21);
	    copy3f(v + 6, tCenterPt);
	    SceneGetCenter(G, offpt);
	    TextSetPosNColor(G, offpt, v);
	    SceneGetScaledAxes(G, I->obj, xn, yn);
	    if (!I->shaderCGO){
	      if (relativeMode & 2){ // label_relative_mode = 1
		  float tmp3f[3];
		  mult3f(xn, tCenterPt[0] * screenwidth/2.f, tmp3f);
		  add3f(tmp3f, offpt, offpt);
		  mult3f(yn, tCenterPt[1] * screenheight/2.f, tmp3f);
		  add3f(tmp3f, offpt, offpt);
		  copy3f(offpt, tCenterPt);
	      } else if (relativeMode & 4){ // label_relative_mode = 2
		  float tmp3f[3];
		  mult3f(xn, (tCenterPt[0] - .5f * screenwidth), tmp3f);
		  add3f(tmp3f, offpt, offpt);
		  mult3f(yn, (tCenterPt[1] - .5f * screenheight), tmp3f);
		  add3f(tmp3f, offpt, offpt);
		  copy3f(offpt, tCenterPt);
	      }
	    }

            if (p) {
              p++;
              if (I->shaderCGO)
                CGOPickColor(I->shaderCGO, p->index, p->bond);
            }

            TextSetPosNColor(G, tCenterPt, v);
	    TextSetTargetPos(G, v + 3);
            st = LexStr(G, *l);
	    TextSetLabelBkgrdInfo(G, *(v + 16), *(v + 17), (v + 18));
	    if (relativeMode & 8){  // label_z_target, adjust z to target
	      TextGetLabelPos(G)[0] = (SceneGetDepth(G, v+3) - .5) * 2.f;
	      TextSetLabelPosIsSet(G, 1);
	    } else if (relativeMode & 6){ // label_relative_mode 1 or 2, i.e., screen stabilized, adjust z 
	      TextSetLabelPos(G, v+6);
	      TextSetLabelPosIsSet(G, 2);
#ifndef PURE_OPENGL_ES_2
	      glDisable(GL_FOG);
#endif
	    } else {
	      TextSetLabelPosIsSet(G, 0);
	    }
#ifndef PURE_OPENGL_ES_2
	    if (!use_shader)
	      glPushMatrix();
#endif
            if (!TextRenderOpenGL(G, info, font_id, st, font_size, v + 12, (draw_var ? 1 : 0), (short)*(v + 15), use_shader, SHADERCGO)){
	      CGOFree(connectorCGO);
	      return;
	    }
	    if (draw_var){
	
	      float *RotMatrix = NULL;
	      float *screenWorldOffset = TextGetScreenWorldOffset(G);
	      float text_width = TextGetWidth(G), text_height = TextGetHeight(G);
	      float *indentFactor = TextGetIndentFactor(G);
	      RotMatrix = SceneGetMatrix(G);

	      if (I->shaderCGO){
		CGODrawConnector(connectorCGO, v + 3, v + 6, text_width, text_height, indentFactor, screenWorldOffset, v + 9, ((short)*(v + 15)), draw_var, *(v + 22), (v + 23), *(v + 27) * font_size / text_height, *(v + 26)) ;
		has_connector = 1;
	      } else {
#ifndef PURE_OPENGL_ES_2
		RepLabelRenderBackgroundInImmediate(G, I, v, draw_var, tCenterPt, relativeMode, xn, yn, PmvMatrix, RotMatrix, screenwidth, screenheight, screenWorldOffset, indentFactor, text_width, text_height, font_size);
#endif
	      }
	    }
#ifdef PURE_OPENGL_ES_2
	    if (float_text && draw_var){
	      TextRenderOpenGL(G, info, font_id, st, font_size, v + 12, 1, (short)*(v + 15), 1, SHADERCGO);
	    }
#else
	    if (!use_shader){
	      glPopMatrix();
	      {
		// for now, render text twice in immediate mode, some cards don't handle
		// z-buffer offset properly, before, only did this when float_text && draw_var
		glPushMatrix();
		TextRenderOpenGL(G, info, font_id, st, font_size, v + 12, 1, (short)*(v + 15), 1, SHADERCGO);
		glPopMatrix();
	      }
	    }
#endif
	    if (relativeMode & 6){ // label_relative_mode 1 or 2, i.e., screen stabilized, adjust z 
#ifndef PURE_OPENGL_ES_2
	      glEnable(GL_FOG);
#endif
	    }
          }
          l++;
          v += 28;
        }
	if (!has_connector){
	  CGOFree(connectorCGO);
	}
	if (connectorCGO){
	  CGOEnd(connectorCGO);
	  CGOStop(connectorCGO);
	}

        if (I->shaderCGO){
	  CGO *totalCGO = NULL;
	  CGO *labelCGO = NULL;
	  CGOStop(I->shaderCGO);
          CGO * tmpCGO = CGONew(G);
          CGOEnable(tmpCGO, GL_LABEL_SHADER);
          CGOSpecial(tmpCGO, SET_LABEL_SCALE_UNIFORMS);
          labelCGO = CGOConvertToLabelShader(I->shaderCGO, tmpCGO);
          CGOAppendNoStop(tmpCGO, labelCGO);
          CGOFreeWithoutVBOs(labelCGO);
          labelCGO = tmpCGO;
	  if (!labelCGO) return;
	  CGOFree(I->shaderCGO);
	  if (connectorCGO){
	    CGO *tmpCGO = NULL;
	    tmpCGO = CGOOptimizeConnectors(connectorCGO, 0);
	    CGOFree(connectorCGO);
	    connectorCGO = tmpCGO;

	    // need to render connector/backgrounds first
	    totalCGO = CGONew(G);
	    CGOEnable(totalCGO, GL_LABEL_FLOAT_TEXT);
	    CGOEnable(totalCGO, GL_CONNECTOR_SHADER);
#ifndef PURE_OPENGL_ES_2
	    CGODisable(totalCGO, GL_LIGHTING);
#endif
	    CGOAppendNoStop(totalCGO, connectorCGO);
	    CGOFreeWithoutVBOs(connectorCGO);
	    CGODisable(totalCGO, GL_CONNECTOR_SHADER);
	    CGOAppendNoStop(totalCGO, labelCGO); 
	    CGOFreeWithoutVBOs(labelCGO);
	    CGODisable(totalCGO, GL_LABEL_FLOAT_TEXT);
	    CGOStop(totalCGO);
	  } else {
	    totalCGO = CGONew(G);
	    CGOEnable(totalCGO, GL_LABEL_FLOAT_TEXT);
	    CGOAppendNoStop(totalCGO, labelCGO);
	    CGODisable(totalCGO, GL_LABEL_FLOAT_TEXT);
            CGODisable(totalCGO, GL_LABEL_SHADER);
	    CGOStop(totalCGO);
	    CGOFreeWithoutVBOs(labelCGO);
	  }
	  I->shaderCGO = totalCGO;
	  if (I->shaderCGO){
	    I->shaderCGO->use_shader = true;
	    I->render(info); // recursion !?
	    return;
	  }
        } else {
#ifndef PURE_OPENGL_ES_2
	  glEnable(GL_LIGHTING);
#endif
	  glEnable(GL_BLEND);
	}
        if(float_text)
          glEnable(GL_DEPTH_TEST);
	info->use_shaders = pre_use_shaders;
      }
    }
  }
}

Rep *RepLabelNew(CoordSet * cs, int state)
{
  PyMOLGlobals *G = cs->G;
  ObjectMolecule *obj;
  int a, a1, c1;
  float *v;
  const float *vc;
  lexidx_t *l;
  int label_color;
  Pickable *rp = NULL;
  AtomInfoType *ai;

  // skip if no labels are visible
  if(!cs->hasRep(cRepLabelBit))
    return NULL;

  auto I = new RepLabel(cs, state);
  obj = cs->Obj;

  label_color = SettingGet_i(G, cs->Setting.get(), obj->Setting.get(), cSetting_label_color);

  /* raytracing primitives */

  I->L = pymol::calloc<lexidx_t>(cs->NIndex);
  ErrChkPtr(G, I->L);
  I->V = pymol::calloc<float>(cs->NIndex * 28);
  ErrChkPtr(G, I->V);

  I->OutlineColor =
    SettingGet_color(G, cs->Setting.get(), obj->Setting.get(), cSetting_label_outline_color);

  if(SettingGet_b(G, cs->Setting.get(), obj->Setting.get(), cSetting_pickable)) {
    I->P = pymol::malloc<Pickable>(cs->NIndex + 1);
    ErrChkPtr(G, I->P);
    rp = I->P + 1;            /* skip first record! */
  }

  I->N = 0;

  v = I->V;
  l = I->L;
  for(a = 0; a < cs->NIndex; a++) {
    a1 = cs->IdxToAtm[a];
    ai = obj->AtomInfo + a1;
    if((ai->visRep & cRepLabelBit) && (ai->label)) {
      int at_label_color = AtomSettingGetWD(G, ai, cSetting_label_color, label_color);

      I->N++;
      if((at_label_color >= 0) ||
         (at_label_color == cColorFront) ||
	 (at_label_color == cColorBack))
        c1 = at_label_color;
      else
        c1 = ai->color;

      /* V - Color, Coordinate, Coordinate + Offset (from label_placement_offset), label_position) */
      vc = ColorGet(G, c1);     /* save new color */
      *(v++) = *(vc++);
      *(v++) = *(vc++);
      *(v++) = *(vc++);

      const float* v0 = cs->coordPtr(a);
      *(v++) = *(v0++);
      *(v++) = *(v0++);
      *(v++) = *(v0++);
      {
	const float *at_label_pos, *at_label_padding;
	const float *con_color;
	float label_connector_width, label_connector_ext_length;
	int label_connector = 0, label_bg = 0, label_bg_outline = 0, at_con_color = 0, at_label_relative_mode = 0, at_label_z_target,
	  label_connector_mode = 0, label_connector_mode_1 = 0, label_connector_mode_2 = 0, label_connector_mode_3 = 0, label_connector_mode_4 = 0, ray_label_connector_flat = 0; 
	float at_label_spacing, at_label_justification, at_label_bkgrd_transp;
	short drawConnector, isProjected, isScreenCoord, isPixelCoord;
	AtomStateGetSetting_i(G, obj, cs, a, ai, cSetting_label_relative_mode, &at_label_relative_mode);
	if (at_label_relative_mode){
	  const float * at_label_screen_point;
	  AtomStateGetSetting(G, obj, cs, a, ai, cSetting_label_screen_point, &at_label_screen_point);
	  copy3f(at_label_screen_point, v);
	  RepLabelAdjustScreenZ(G, v);
	} else {
	  const float * at_label_place;
	  AtomStateGetSetting(G, obj, cs, a, ai, cSetting_label_placement_offset, &at_label_place);
	  add3f(at_label_place, v - 3, v);
	}
	v += 3;
	AtomStateGetSetting_color(G, obj, cs, a, ai, cSetting_label_connector_color, &at_con_color);

	/* behave just like the label color */
	if(!((at_con_color >= 0) ||
	     (at_con_color == cColorFront) ||
	     (at_con_color == cColorBack)))
	  at_con_color = ai->color;

	con_color = ColorGet(G, at_con_color);
	copy3f(con_color, v);
	v += 3;

	AtomStateGetSetting_b(G, obj, cs, a, ai, cSetting_ray_label_connector_flat, &ray_label_connector_flat);
	AtomStateGetSetting_b(G, obj, cs, a, ai, cSetting_label_bg_outline, &label_bg_outline);
	AtomStateGetSetting_b(G, obj, cs, a, ai, cSetting_label_connector, &label_connector);
	AtomStateGetSetting_i(G, obj, cs, a, ai, cSetting_label_connector_mode, &label_connector_mode);
	AtomStateGetSetting_i(G, obj, cs, a, ai, cSetting_label_z_target, &at_label_z_target);

	AtomStateGetSetting(G, obj, cs, a, ai, cSetting_label_position, &at_label_pos);
	copy3f(at_label_pos, v);
	v += 3;


	AtomStateGetSetting_f(G, obj, cs, a, ai, cSetting_label_multiline_spacing, &at_label_spacing);
	AtomStateGetSetting_f(G, obj, cs, a, ai, cSetting_label_multiline_justification, &at_label_justification);
	at_label_justification = CLAMP_VALUE(at_label_justification , -1.f, 1.f);
	AtomStateGetSetting(G, obj, cs, a, ai, cSetting_label_padding, &at_label_padding);
	AtomStateGetSetting_f(G, obj, cs, a, ai, cSetting_label_bg_transparency, &at_label_bkgrd_transp);
	AtomStateGetSetting_color(G, obj, cs, a, ai, cSetting_label_bg_color, &at_con_color);
	label_bg = (at_con_color != -1) && (at_label_bkgrd_transp < 1.f); // if the color is not default and transparency is < 1., then draw bg
	drawConnector = (label_connector > 0 || label_bg || label_bg_outline > 0) ? 1 : 0;
	if (at_label_z_target < 0){
	  // defaults to z if draw connector
	  if (drawConnector){
	    at_label_z_target = 1;
	  } else {
	    at_label_z_target = 0;
	  }
	} else {
	  at_label_z_target = at_label_z_target ? 1 : 0;
	}
	isProjected = (at_label_relative_mode < 1) ? 1 : 0;
	isScreenCoord = (at_label_relative_mode == 1) ? 1 : 0;
	isPixelCoord = (isProjected + isScreenCoord) ? 0 : 1;
	*(v++) = (float)(drawConnector + isScreenCoord * 2 + isPixelCoord * 4 + at_label_z_target * 8);
	*(v++) = at_label_spacing;
	*(v++) = at_label_justification;
	copy3f(at_label_padding, v);
	v += 3;
	label_connector_mode_1 = label_connector_mode == 1;
	label_connector_mode_2 = label_connector_mode == 2;
	label_connector_mode_3 = label_connector_mode == 3;
	label_connector_mode_4 = label_connector_mode == 4;

	*(v++) = (float)(label_connector + 2 * label_bg + 4 * label_bg_outline + 8 * label_connector_mode_1 +
			 16 * label_connector_mode_2 + 32 * label_connector_mode_3 + 64 * label_connector_mode_4 + 128 * ray_label_connector_flat);

	*(v++) = 1.f-at_label_bkgrd_transp;

	/* behave just like the label color */
	if(!((at_con_color >= 0) ||
	     (at_con_color == cColorFront) ||
	     (at_con_color == cColorBack)))
	  at_con_color = ai->color;
	con_color = ColorGet(G, at_con_color);
	copy3f(con_color, v);
	v += 3;
	AtomStateGetSetting_f(G, obj, cs, a, ai, cSetting_label_connector_width, &label_connector_width);
	*(v++) = DIP2PIXEL(label_connector_width);
	AtomStateGetSetting_f(G, obj, cs, a, ai, cSetting_label_connector_ext_length, &label_connector_ext_length);
	*(v++) = label_connector_ext_length;
      }
      if(rp) {
        rp->index = a1;
        rp->bond = ai->masked ? cPickableNoPick : cPickableLabel;      /* label indicator */
        rp++;
      }
      *(l++) = ai->label;
    }
  }

  if(I->N) {
    I->V = ReallocForSure(I->V, float, (v - I->V));
    I->L = ReallocForSure(I->L, lexidx_t, (l - I->L));
    if(rp) {
      I->P = ReallocForSure(I->P, Pickable, (rp - I->P));
      I->P[0].index = I->N;   /* unnec? */
    }
  } else {
    I->V = ReallocForSure(I->V, float, 1);
    I->L = ReallocForSure(I->L, lexidx_t, 1);
    if(rp) {
      FreeP(I->P);
    }
  }
  return (Rep *) I;
}
