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

#include"os_predef.h"
#include"os_std.h"
#include"os_gl.h"

#include"OOMac.h"
#include"ObjectGadgetRamp.h"
#include"GadgetSet.h"
#include"Base.h"
#include"MemoryDebug.h"
#include"CGO.h"
#include"Scene.h"
#include"Setting.h"
#include"PConv.h"
#include"main.h"
#include"Color.h"
#include"VFont.h"
#include"ObjectMolecule.h"
#include"Executive.h"
#include"Util.h"

#include"P.h"

void ObjectGadgetRampFree(ObjectGadgetRamp *I) {
  ColorForgetExt(I->Gadget.Obj.Name);
  VLAFreeP(I->Level);
  VLAFreeP(I->Color);
  ObjectGadgetPurge(&I->Gadget);
  OOFreeP(I);
}

#define ShapeVertex(cgo,a,b) CGOVertex(cgo,(float)a,(float)b,0.0F)
#define ShapeFVertex(cgo,a,b) CGOFontVertex(cgo,(float)a,(float)b,0.0F)
#define ABS 0.0F
#define REL 1.0F
#define OFF 2.0F

#define ShapeNormal(cgo,a,b) CGONormal(cgo,(float)a,(float)b,0.0F)
#define ShapeColor(cgo,a,b) CGONormal(cgo,(float)a,(float)b,0.0F)
#define LKP 2.0F

static void ObjectGadgetRampBuild(ObjectGadgetRamp *I);

PyObject *ObjectGadgetRampAsPyList(ObjectGadgetRamp *I)
{

  PyObject *result = NULL;

  result = PyList_New(8);

  PyList_SetItem(result,0,ObjectGadgetPlainAsPyList(&I->Gadget));
  PyList_SetItem(result,1,PyInt_FromLong(I->RampType));
  PyList_SetItem(result,2,PyInt_FromLong(I->NColor));
  if(I->Level&&I->NColor) {
    PyList_SetItem(result,3,PConvFloatVLAToPyList(I->Level));
  } else {
    PyList_SetItem(result,3,PConvAutoNone(NULL));
  }
  if(I->Color&&I->NColor) {
    PyList_SetItem(result,4,PConvFloatVLAToPyList(I->Color));
  } else {
    PyList_SetItem(result,4,PConvAutoNone(NULL));
  }
  PyList_SetItem(result,5,PyInt_FromLong(I->var_index));
  PyList_SetItem(result,6,PyString_FromString(I->SrcName));
  PyList_SetItem(result,7,PyInt_FromLong(I->SrcState));

  return(PConvAutoNone(result));  
}

int ObjectGadgetRampNewFromPyList(PyObject *list,ObjectGadgetRamp **result)
{
  
  ObjectGadgetRamp *I = NULL;
  int ok = true;
  int ll;

  if(ok) I=ObjectGadgetRampNew();
  if(ok) ok = (I!=NULL);
  if(ok) ok = (list!=NULL);
  if(ok) ok = PyList_Check(list);
  if(ok) ll = PyList_Size(list);
  /* TO SUPPORT BACKWARDS COMPATIBILITY...
   Always check ll when adding new PyList_GetItem's */

  if(ok) ok = ObjectGadgetInitFromPyList(PyList_GetItem(list,0),&I->Gadget);
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,1),&I->RampType);
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,2),&I->NColor);
  if(ok&&I->NColor) ok = PConvPyListToFloatVLA(PyList_GetItem(list,3),&I->Level);
  if(ok&&I->NColor) ok = PConvPyListToFloatVLA(PyList_GetItem(list,4),&I->Color);
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,5),&I->var_index);
  if(ok) ok = PConvPyStrToStr(PyList_GetItem(list,6),I->SrcName,ObjNameMax);
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,7),&I->SrcState);
  /*  if(ok) ObjectGadgetRampBuild(I);
      if(ok) ObjectGadgetRampUpdate(I);*/
  if(ok) ObjectGadgetUpdateStates(&I->Gadget);
  if(ok) ObjectGadgetUpdateExtents(&I->Gadget);
  if(ok) (*result)=I;
  return(ok);
}

int ObjectGadgetRampInterVertex(ObjectGadgetRamp *I,float *pos,float *color)
{
  float level;
  int ok=true;
  if (!ExecutiveValidateObjectPtr((CObject*)I->Map,cObjectMap));
    I->Map = ExecutiveFindObjectMapByName(I->SrcName);
  if(ok) ok = (I->Map!=NULL);
  if(ok) ok = ObjectMapInterpolate(I->Map,I->SrcState,pos,&level,1);
  if(ok) ok = ObjectGadgetRampInterpolate(I,level,color);
  return(ok);
}

int ObjectGadgetRampInterpolate(ObjectGadgetRamp *I,float level,float *color)
{
  int i=0;
  int ok=true;
  int below=0;
  int above=0;
  float d,x0,x1;
  int a;

  if(I->Level&&I->Color) {
    while(i<I->NColor) {
      if(I->Level[i]>level) {
        above=i;
        break;
      } else {
        below=i;
        above=i;
      }
      i++;
    }
    if(above!=below) {
      d=I->Level[above]-I->Level[below];
      if(fabs(d)>R_SMALL8) {
        x0=(level-I->Level[below])/d;
        x1=1.0F-x0;
        for(a=0;a<3;a++) {
          color[a]=x0*I->Color[3*above+a]+x1*I->Color[3*below+a];
        }
        clamp3f(color);
      } else {
        copy3f(I->Color+3*above,color);
        clamp3f(color);
      }
    } else {
      copy3f(I->Color+3*above,color);
     clamp3f(color);
    }
  } else {
    color[0]=1.0F;
    color[1]=1.0F;
    color[2]=1.0F;
  }
  return(ok);
}

static void ObjectGadgetRampUpdateCGO(ObjectGadgetRamp *I,GadgetSet *gs)
{
  CGO *cgo;
  int n_extra;
  int a,c=0;
  float *p;
  char buffer[255];

  cgo = CGONewSized(100);

  /* behind text */

  CGOBegin(cgo,GL_TRIANGLE_STRIP);
  CGOColor(cgo,0.05F,0.05F,0.05F);
  ShapeNormal(cgo,LKP,2);
  ShapeVertex(cgo,REL,9);
  ShapeVertex(cgo,REL,10);
  ShapeVertex(cgo,REL,7);
  ShapeVertex(cgo,REL,8);
  CGOEnd(cgo);

  CGOColor(cgo,1.0F,1.0F,1.0F);
  CGOFontScale(cgo,I->text_scale_h,I->text_scale_v);

  if(I->Level&&I->NColor) {
    sprintf(buffer,"%0.6f",I->Level[0]);
    ShapeFVertex(cgo,REL,11);
    CGOWrite(cgo,buffer);
    sprintf(buffer,"%0.6f",I->Level[I->NColor-1]);
    ShapeFVertex(cgo,REL,12);
    CGOWriteLeft(cgo,buffer);
  }

  CGOBegin(cgo,GL_TRIANGLE_STRIP);
  ShapeNormal(cgo,LKP,2);

  n_extra = 3*(I->NColor * 1);

  if(I->NColor<2) {
    /* TO DO */
  } else {
    VLACheck(gs->Coord,float,(I->var_index+n_extra)*3);      
    c = I->var_index;
    p=gs->Coord+3*c;
    for(a=0;a<I->NColor;a++) {

      CGOColorv(cgo,I->Color+3*a);
      
      *(p++) = I->border + (I->width * a)/(I->NColor-1);
      *(p++) = -I->border;
      *(p++) = I->border;
      ShapeVertex(cgo,REL,c);
      c++;

      *(p++) = I->border + (I->width * a)/(I->NColor-1);
      *(p++) = -(I->border+I->bar_height);
      *(p++) = I->border;
      ShapeVertex(cgo,REL,c);
      c++;

      *(p++) = I->border + (I->width * a)/(I->NColor-1);
      *(p++) = -(I->border+I->height+I->height);
      *(p++) = I->border;
      c++;

    }

  }
  gs->NCoord = c;


  CGOColor(cgo,0.5F,0.5F,0.5F);

  /* top */
  CGOBegin(cgo,GL_TRIANGLE_STRIP);
  ShapeNormal(cgo,LKP,2);
  ShapeVertex(cgo,REL,5);
  ShapeVertex(cgo,REL,6);
  ShapeNormal(cgo,LKP,1);
  ShapeVertex(cgo,REL,1);
  ShapeVertex(cgo,REL,2);
  CGOEnd(cgo);

  /* bottom */

  CGOBegin(cgo,GL_TRIANGLE_STRIP);
  ShapeNormal(cgo,LKP,4);
  ShapeVertex(cgo,REL,3);
  ShapeVertex(cgo,REL,4);
  ShapeNormal(cgo,LKP,2);
  ShapeVertex(cgo,REL,7);
  ShapeVertex(cgo,REL,8);
  CGOEnd(cgo);

  /* left */

  CGOBegin(cgo,GL_TRIANGLE_STRIP);
  ShapeNormal(cgo,LKP,3);
  ShapeVertex(cgo,REL,1);
  ShapeVertex(cgo,REL,3);
  ShapeNormal(cgo,LKP,2);
  ShapeVertex(cgo,REL,5);
  ShapeVertex(cgo,REL,7);
  CGOEnd(cgo);

  /* right */
  CGOBegin(cgo,GL_TRIANGLE_STRIP);
  ShapeNormal(cgo,LKP,2);
  ShapeVertex(cgo,REL,6);
  ShapeVertex(cgo,REL,8);
  ShapeNormal(cgo,LKP,0);
  ShapeVertex(cgo,REL,2);
  ShapeVertex(cgo,REL,4);
  CGOEnd(cgo);

  /* center */
  CGOEnd(cgo);

  CGOStop(cgo);

  CGOFree(gs->ShapeCGO);
  gs->ShapeCGO = cgo;

  CGOPreloadFonts(gs->ShapeCGO);
  
  cgo = CGONewSized(100);
  CGODotwidth(cgo,5);

  CGOPickColor(cgo,0,0);

  /* top */
  CGOBegin(cgo,GL_TRIANGLE_STRIP);
  ShapeVertex(cgo,REL,1);
  ShapeVertex(cgo,REL,2);
  ShapeVertex(cgo,REL,5);
  ShapeVertex(cgo,REL,6);

  CGOEnd(cgo);

  /* bottom */

  CGOBegin(cgo,GL_TRIANGLE_STRIP);
  ShapeVertex(cgo,REL,3);
  ShapeVertex(cgo,REL,4);
  ShapeVertex(cgo,REL,7);
  ShapeVertex(cgo,REL,8);
  CGOEnd(cgo);

  /* left */

  CGOBegin(cgo,GL_TRIANGLE_STRIP);
  ShapeVertex(cgo,REL,1);
  ShapeVertex(cgo,REL,3);
  ShapeVertex(cgo,REL,5);
  ShapeVertex(cgo,REL,7);
  CGOEnd(cgo);

  /* right */
  CGOBegin(cgo,GL_TRIANGLE_STRIP);
  ShapeVertex(cgo,REL,6);
  ShapeVertex(cgo,REL,8);
  ShapeVertex(cgo,REL,2);
  ShapeVertex(cgo,REL,4);
  CGOEnd(cgo);

  /* band */

  CGOPickColor(cgo,13,0);
  CGOBegin(cgo,GL_TRIANGLE_STRIP);
  ShapeVertex(cgo,REL,5);
  ShapeVertex(cgo,REL,6);
  ShapeVertex(cgo,REL,7);
  ShapeVertex(cgo,REL,8);
  CGOEnd(cgo);
  
  CGOEnd(cgo);
  CGOStop(cgo);

  CGOFree(gs->PickShapeCGO);
  gs->PickShapeCGO = cgo;
}

static void ObjectGadgetRampBuild(ObjectGadgetRamp *I)
{
  GadgetSet *gs = NULL;
  ObjectGadget *og;
  int a;

  float coord[100];
  int ix=0;

  float normal[] = {
    1.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 0.0, 1.0,
   -1.0, 0.0, 0.0,
    0.0,-1.0, 0.0,
  };
#define VV(a,b,c) {coord[ix++]=a;coord[ix++]=b;coord[ix++]=c;};

VV(    I->x ,  I->y , 0.3 );

    /* outer points */

VV(    0.0 ,  0.0 , 0.0 );
    VV(I->width+I->border*2 ,  0.0 , 0.0 );
VV(    0.0 , -(I->height+I->border*2) , 0.0 );
VV(    I->width+I->border*2 , -(I->height+I->border*2) , 0.0 );

VV(    I->border, -I->border, I->border);
VV(    I->width+I->border,-I->border, I->border);
VV(    I->border,  -(I->height+I->border), I->border);
VV(    I->width+I->border, -(I->height+I->border), I->border);

VV(    I->border, -(I->border+I->bar_height), I->border);
VV(    I->width+I->border,-(I->border+I->bar_height), I->border);
    
VV(    I->border+I->text_border, I->text_border-(I->border+I->height), I->border+I->text_raise);
VV(    I->width+I->border,I->text_border-(I->border+I->height), I->border+I->text_raise);

 VV(   0.0,0.0,0.0);
#undef VV

  OrthoBusyPrime();

  og = &I->Gadget;
  gs = GadgetSetNew();

  gs->NCoord = 14;
  I->var_index = gs->NCoord;
  gs->Coord = VLAlloc(float,gs->NCoord*3);
  for(a=0;a<gs->NCoord*3;a++) {
    gs->Coord[a]=coord[a];
  }

  gs->NNormal = 5;
  gs->Normal = VLAlloc(float,gs->NNormal*3);
  for(a=0;a<gs->NNormal*3;a++) {
    gs->Normal[a]=normal[a];
  }

  og->GSet[0] = gs;
  og->NGSet = 1;
  og->Obj.Context=1;
  gs->Obj = (ObjectGadget*)I;

  ObjectGadgetRampUpdateCGO(I,gs);
  gs->fUpdate(gs);
  
}

/*========================================================================*/
void ObjectGadgetRampUpdate(ObjectGadgetRamp *I)
{
  float scale;

  scale = (1.0F+5*I->Gadget.GSet[0]->Coord[13*3]);

  I->Gadget.GSet[0]->Coord[13*3] = 0.0;
  if(I->NColor==3) {
    I->Level[0]=(I->Level[0]-I->Level[1])*scale+I->Level[1];
    I->Level[2]=(I->Level[2]-I->Level[1])*scale+I->Level[1];
    ExecutiveInvalidateRep(cKeywordAll,cRepAll,cRepInvColor);
  }
  if(I->Gadget.NGSet)
    if(I->Gadget.GSet[0]) {
      ObjectGadgetRampUpdateCGO(I,I->Gadget.GSet[0]);
      ObjectGadgetUpdateStates(&I->Gadget);
    }
  ObjectGadgetUpdateExtents(&I->Gadget);
  SceneChanged();
}

/*========================================================================*/
ObjectGadgetRamp *ObjectGadgetRampMapNewAsDefined(ObjectMap *map,PyObject *level,
                                               PyObject *color,int map_state)
{
  ObjectGadgetRamp *I;
  int ok = true;
  
  I = ObjectGadgetRampNew();
  I->RampType = cRampMap;

  PBlock();
  if(ok) ok = PConvPyListToFloatVLA(level,&I->Level);
  if(ok) ok = PConvPyList3ToFloatVLA(color,&I->Color);
  if(ok) I->NColor=VLAGetSize(I->Level);
  ObjectGadgetRampBuild(I);
  UtilNCopy(I->SrcName,map->Obj.Name,ObjNameMax);
  I->SrcState=map_state;


  /* test interpolate 
     { 
    float test[3];

    ObjectGadgetRampInterpolate(I,-2.0,test);
    dump3f(test,"test color");
    ObjectGadgetRampInterpolate(I,-1.0,test);
    dump3f(test,"test color");
    ObjectGadgetRampInterpolate(I,-0.9,test);
    dump3f(test,"test color");
    ObjectGadgetRampInterpolate(I,-0.5,test);
    dump3f(test,"test color");
    ObjectGadgetRampInterpolate(I,0.0,test);
    dump3f(test,"test color");
    ObjectGadgetRampInterpolate(I,0.5,test);
    dump3f(test,"test color");
    ObjectGadgetRampInterpolate(I,1.0,test);
    dump3f(test,"test color");
    ObjectGadgetRampInterpolate(I,2.0,test);
    dump3f(test,"test color");
  }
  */

  PUnblock();
  return(I);
}

/*========================================================================*/
ObjectGadgetRamp *ObjectGadgetRampNew(void)
{
  OOAlloc(ObjectGadgetRamp);

  ObjectGadgetInit(&I->Gadget);
  I->Gadget.GadgetType = cGadgetRamp;
  I->RampType = 0;
  I->NColor = 0;
  I->Level = NULL;
  I->Color = NULL;
  I->SrcName[0] = 0;

  I->Gadget.Obj.fUpdate =(void (*)(struct CObject *)) ObjectGadgetRampUpdate;
  I->Gadget.Obj.fFree =(void (*)(struct CObject *)) ObjectGadgetRampFree;
  I->width = 0.9F;
  I->height = 0.06F;
  I->bar_height = 0.03F;
  I->text_raise = 0.003F;
  I->text_border = 0.004F;
  I->text_scale_h = 0.04F;
  I->text_scale_v = 0.02F;
  I->border = 0.018F;
  I->var_index = 0;
  I->x = (1.0-(I->width+2*I->border))/2.0F;
  I->y = 0.12F;
  I->Map = NULL;
  return(I);
}

