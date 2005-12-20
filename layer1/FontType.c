/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* copyright 1998-2003 by Warren Lyford Delano of DeLano Scientific. 
D* -------------------------------------------------------------------
E* It is unlawful to modify or remove this copyright notice.
F* -------------------------------------------------------------------
G* Please see the accompanying LICENSE file for further information. 
H* --------------------------------------------------\-----------------
I* Additional authors of this source file include:
-* 
-* 
-*
Z* -------------------------------------------------------------------
*/

#include "MemoryDebug.h"
#include "OOMac.h"
#include "os_gl.h"
#include "FontType.h"
#include "Text.h"
#include "Ortho.h"
#include "Scene.h"
#include "Character.h"
#include "Util.h"
#include "TypeFace.h"

typedef struct {
  CFont Font; /* must be first */
  PyMOLGlobals *G;
  CTypeFace *TypeFace;
} CFontType;


__inline__ static char *_FontTypeRenderOpenGL(RenderInfo *info, 
                                              CFontType *I,char *st,
                                              float size,int flat, float *rpos)
{
  register PyMOLGlobals *G = I->Font.G;
  if(G->ValidContext) {
    int c;
    int pushed = OrthoGetPushed(G);
    int kern_flag = false;
    int last_c = -1;
    int sampling = 1;
    float x_indent=0.0F, y_indent=0.0F;
    if(info)
      sampling = info->sampling;
    if(st&&(*st)) {
      float origin[3], v_scale;
      SceneOriginGet(G,origin);
      v_scale = SceneGetScreenVertexScale(G,origin);

      if(size<0.0F) {
        size = (int)(0.5F-size/v_scale);
      }

      if(rpos) {
        if(rpos[0]<1.0F) { /* we need to measure the string width before starting to draw */
          float factor = rpos[0]/2.0F - 0.5F;
          char *sst = st;
          if(factor<-1.0F) factor = -1.0F;
          if(factor>0.0F) factor = 0.0F;
          while((c=*(sst++))) {
            
            CharFngrprnt fprnt;
            unsigned char *rgba;
            UtilZeroMem(&fprnt,sizeof(fprnt));
            fprnt.u.i.text_id = I->Font.TextID;
            fprnt.u.i.size = (int)(size*64*sampling);
            rgba = fprnt.u.i.color;
            TextGetColorUChar(G,rgba,rgba+1,rgba+2,rgba+3);
            rgba = fprnt.u.i.outline_color;
            TextGetOutlineColor(G,rgba,rgba+1,rgba+2,rgba+3);
            fprnt.u.i.ch = c;
            fprnt.u.i.flat = flat;
            {
              int id = CharacterFind(G,&fprnt);
              if(!id) {
                id = TypeFaceCharacterNew(I->TypeFace,&fprnt,size*sampling);
              }
              if(id) {
                if(kern_flag) {
                  x_indent -= factor * (TypeFaceGetKerning(I->TypeFace, 
                                                           (unsigned int)last_c,
                                                           (unsigned int)c,
                                                           size)/sampling);
                }
                x_indent -= factor * CharacterGetAdvance(G,sampling,id);
              }
            }
            kern_flag = true;
            last_c = c;
          }
        }
        if(rpos[0]<-1.0F) {
          x_indent -= (rpos[0]+1.0F)/v_scale;
        } else if(rpos[0]>1.0F) {
          x_indent -= (rpos[0]-1.0F)/v_scale;
        }
        if(rpos[1]<1.0F) {
          float factor = -rpos[1]/2.0F + 0.5F;
          if(factor>1.0F) factor = 1.0F;
          if(factor<0.0F) factor = 0.0F;
          y_indent = 0.75*size*factor;
        }
        if(rpos[1]<-1.0F) {
          y_indent -= (rpos[1]+1.0F)/v_scale;
        } else if(rpos[1]>1.0F) {
          y_indent -= (rpos[1]-1.0F)/v_scale;
        }
      }
      if(!pushed) {
        float *v = TextGetPos(G);
        float loc[3];
        float zero[3]= {0.0F,0.0F,0.0F};
        if(rpos) {
          SceneGetEyeNormal(G,v,loc);
          scale3f(loc,rpos[2],loc);
          add3f(v,loc,loc);
          v = loc;
        }
        ScenePushRasterMatrix(G,v);
        TextSetPos(G,zero);
      } 
      if(rpos) {
        TextIndent(G,x_indent,y_indent);
      }
      while((c=*(st++))) {

        CharFngrprnt fprnt;
        unsigned char *rgba;
        UtilZeroMem(&fprnt,sizeof(fprnt));
        fprnt.u.i.text_id = I->Font.TextID;
        fprnt.u.i.size = (int)(size*64*sampling);
        rgba = fprnt.u.i.color;
        TextGetColorUChar(G,rgba,rgba+1,rgba+2,rgba+3);
        rgba = fprnt.u.i.outline_color;
        TextGetOutlineColor(G,rgba,rgba+1,rgba+2,rgba+3);
        fprnt.u.i.ch = c;
        fprnt.u.i.flat = flat;
        {
          int id = CharacterFind(G,&fprnt);
          if(!id) {
            id = TypeFaceCharacterNew(I->TypeFace,&fprnt,size*sampling);
          }
          if(id) {
            if(kern_flag) {
              TextAdvance(G, TypeFaceGetKerning(I->TypeFace, 
                                                (unsigned int)last_c,
                                                (unsigned int)c,
                                                size)/sampling);
            }
            CharacterRenderOpenGL(G,info,id); /* handles advance */
          }
        }
        kern_flag = true;
        last_c = c;
      }
      if(!pushed) {
        ScenePopRasterMatrix(G);
      }
    }
  }
  return st;
}

static char *FontTypeRenderOpenGL(RenderInfo *info, CFontType *I,char *st,float size,float *rpos)
{
  return _FontTypeRenderOpenGL(info,I,st,size,false,rpos);
}
static char *FontTypeRenderOpenGLFlat(RenderInfo *info, CFontType *I,char *st,float size,float *rpos)
{
  return _FontTypeRenderOpenGL(info,I,st,size,true,rpos);
}

static char *FontTypeRenderRay(CRay *ray, CFontType *I,char *st,float size, float *rpos)
{
  register PyMOLGlobals *G = I->Font.G;
  int c;
  int kern_flag = false;
  int last_c = -1;
  int sampling = ray->Sampling;
  float x_indent=0.0F, y_indent=0.0F;
  float xn[3], yn[3], x_adj[3], y_adj[3], pos[3], *v;

  if(st&&(*st)) {
    float origin[3],v_scale;
    SceneOriginGet(G,origin);
    v_scale = SceneGetScreenVertexScale(G,origin);

    if(rpos) {
      float loc[3];
      v = TextGetPos(I->G);
      SceneGetEyeNormal(G,v,loc);
      scale3f(loc,rpos[2],loc);
      add3f(v,loc,loc);
      TextSetPos(I->G,loc);
    }

    RayGetScaledAxes(ray,xn,yn);
    
    if(size<0.0F) {

      size = (int)(0.5F - size / v_scale);
    }

    if(rpos) {

      if(rpos[0]<1.0F) { /* we need to measure the string width before starting to draw */
        float factor = rpos[0]/2.0F - 0.5F;
        char *sst = st;
        if(factor<-1.0F) factor = -1.0F;
        if(factor>0.0F) factor = 0.0F;
        while((c=*(sst++))) {
          
          CharFngrprnt fprnt;
          unsigned char *rgba;
          UtilZeroMem(&fprnt,sizeof(fprnt));
          fprnt.u.i.text_id = I->Font.TextID;
          fprnt.u.i.size = (int)(size*64*sampling);
          rgba = fprnt.u.i.color;
          TextGetColorUChar(G,rgba,rgba+1,rgba+2,rgba+3);
          rgba = fprnt.u.i.outline_color;
          TextGetOutlineColor(G,rgba,rgba+1,rgba+2,rgba+3);
          fprnt.u.i.ch = c;
          {
            int id = CharacterFind(G,&fprnt);
            if(!id) {
              id = TypeFaceCharacterNew(I->TypeFace,&fprnt,size*sampling);
            }
            if(id) {
              if(kern_flag) {
                x_indent -= factor * TypeFaceGetKerning(I->TypeFace, 
                                                        (unsigned int)last_c,
                                                        (unsigned int)c,
                                                        size*sampling)/sampling;
              }
              x_indent -= factor * CharacterGetAdvance(G,1,id);
              kern_flag = true;
              last_c = c;
            }
          }
        }
      }
      if(rpos[0]<-1.0F) {
        x_indent -= 2*(rpos[0]+1.0F)/v_scale;
      } else if(rpos[0]>1.0F) {
        x_indent -= 2*(rpos[0]-1.0F)/v_scale;
      }
      if(rpos[1]<1.0F) {
        float factor = -rpos[1]/2.0F + 0.5F;
        if(factor>1.0F) factor = 1.0F;
        if(factor<0.0F) factor = 0.0F;
        y_indent = 0.75F*sampling*size*factor;
      }
      if(rpos[1]<-1.0F) {
        y_indent -= 2*(rpos[1]+1.0F)/v_scale;
      } else if(rpos[1]>1.0F) {
        y_indent -= 2*(rpos[1]-1.0F)/v_scale;
      }
      v = TextGetPos(I->G);
      scale3f(xn, x_indent, x_adj);
      scale3f(yn, y_indent, y_adj);
      subtract3f(v,x_adj,pos);
      subtract3f(pos,y_adj,pos);
      TextSetPos(I->G,pos);
    }
    kern_flag = false;

    while((c=*(st++))) {
      
      CharFngrprnt fprnt;
      unsigned char *rgba;
      UtilZeroMem(&fprnt,sizeof(fprnt));
      fprnt.u.i.text_id = I->Font.TextID;
      fprnt.u.i.size = (int)(size*64*sampling);
      rgba = fprnt.u.i.color;
      TextGetColorUChar(G,rgba,rgba+1,rgba+2,rgba+3);
      rgba = fprnt.u.i.outline_color;
      TextGetOutlineColor(G,rgba,rgba+1,rgba+2,rgba+3);
      fprnt.u.i.ch = c;
      {
        int id = CharacterFind(G,&fprnt);
        if(!id) {
          id = TypeFaceCharacterNew(I->TypeFace,&fprnt,size*sampling);
        }
        if(id) {
          if(kern_flag) {
            float kern = TypeFaceGetKerning(I->TypeFace, 
                                            (unsigned int)last_c,
                                            (unsigned int)c,
                                            size*sampling)/sampling;
            v = TextGetPos(I->G);
            scale3f(xn, kern, x_adj);
            add3f(v,x_adj,pos);
            TextSetPos(I->G,pos);
          }
          ray->fCharacter(ray,id); /* handles advance */

          kern_flag = true;
          last_c = c;
        }
      }
    }
  }
  return st;
}

static void FontTypeFree(CFont* font)
{
  register CFontType *I = (CFontType*)font;
  TypeFaceFree(I->TypeFace);
  OOFreeP(I);
}

CFont* FontTypeNew(PyMOLGlobals *G,unsigned char *dat,unsigned int len)
{  
  OOAlloc(G,CFontType);
  FontInit(G,&I->Font);
  I->G = G;
  I->Font.fRenderOpenGL = (FontRenderOpenGLFn*)FontTypeRenderOpenGL;
  I->Font.fRenderOpenGLFlat = (FontRenderOpenGLFn*)FontTypeRenderOpenGLFlat;
  I->Font.fRenderRay = (FontRenderRayFn*)FontTypeRenderRay;
  I->Font.fFree = FontTypeFree;
  I->TypeFace = TypeFaceLoad(G,dat,len);
  if(!I->TypeFace) {
    OOFreeP(I);
  }
  return (CFont*)I;
}


