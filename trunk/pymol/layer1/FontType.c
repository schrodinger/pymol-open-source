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


static char *FontTypeRenderOpenGL(RenderInfo *info, CFontType *I,char *st,float size)
{
  register PyMOLGlobals *G = I->Font.G;
  if(G->ValidContext) {
    int c;
    int pushed = OrthoGetPushed(G);
    int kern_flag = false;
    int last_c = -1;
    int sampling = 1;
    if(info)
      sampling = info->sampling;
    if(st&&(*st)) {

      if(size<0.0F) {
        float origin[3];
        SceneOriginGet(G,origin);
        size = (int)(0.5F-size/SceneGetScreenVertexScale(G,origin));
      }
      if(!pushed) {
        float *v = TextGetPos(G);
        float zero[3]= {0.0F,0.0F,0.0F};
        ScenePushRasterMatrix(G,v);
        TextSetPos(G,zero);
      } 
     
      while((c=*(st++))) {

        CharFngrprnt fprnt;
        unsigned char *rgba;
        UtilZeroMem(&fprnt,sizeof(fprnt));
        fprnt.u.i.text_id = I->Font.TextID;
        fprnt.u.i.size = (int)(size*64*sampling);
        rgba = fprnt.u.i.color;
        TextGetColorUChar(G,rgba,rgba+1,rgba+2,rgba+3);
        fprnt.u.i.ch = c;
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

static char *FontTypeRenderRay(CRay *ray, CFontType *I,char *st,float size)
{
  register PyMOLGlobals *G = I->Font.G;
  int c;
  if(st&&(*st)) {
    
    if(size<0.0F) {
      float origin[3];
      SceneOriginGet(G,origin);
      size = (int)(0.5F-size/SceneGetScreenVertexScale(G,origin));
    }

    while((c=*(st++))) {
      
      CharFngrprnt fprnt;
      unsigned char *rgba;
      UtilZeroMem(&fprnt,sizeof(fprnt));
      fprnt.u.i.text_id = I->Font.TextID;
      fprnt.u.i.size = (int)(size*64*ray->Sampling);
      rgba = fprnt.u.i.color;
      TextGetColorUChar(G,rgba,rgba+1,rgba+2,rgba+3);
      fprnt.u.i.ch = c;
      {
        int id = CharacterFind(G,&fprnt);
        if(!id) {
          id = TypeFaceCharacterNew(I->TypeFace,&fprnt,size*ray->Sampling);
        }
        if(id) {
          ray->fCharacter(ray,id); /* handles advance */
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
  I->Font.fRenderRay = (FontRenderRayFn*)FontTypeRenderRay;
  I->Font.fFree = FontTypeFree;
  I->TypeFace = TypeFaceLoad(G,dat,len);
  if(!I->TypeFace) {
    OOFreeP(I);
  }
  return (CFont*)I;
}


