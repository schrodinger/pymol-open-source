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

#include"os_std.h"
#include"os_gl.h"

#include"MemoryDebug.h"
#include"Vector.h"
#include"Matrix.h"
#include"ButMode.h"
#include"Scene.h"
#include"Editor.h"
#include"Selector.h"
#include"Ortho.h"
#include"main.h"
#include"Color.h"
#include"Setting.h"

typedef struct {
  ObjectMolecule *Obj;
  int ActiveState;
  int DragIndex;
  int DragSelection;
  int DragSele0;
  int DragSele1;
  int DragHaveAxis,DragHaveBase,DragBondFlag;
  ObjectMolecule *DragObject;
  int NFrag;
  float V0[3],V1[3],Axis[3],Center[3],DragBase[3];

}  CEditor;

CEditor Editor;

int EditorActive(void) {
  CEditor *I = &Editor;
  return(I->Obj!=NULL);
}

ObjectMolecule *EditorDragObject(void)
{
  CEditor *I = &Editor;
  return(I->DragObject);
}
static void subdivide( int n, float *x, float *y);

static void subdivide( int n, float *x, float *y)
{
  int a;
  if(n<3) {n=3;}
  for(a=0;a<=n;a++)
	 {
		x[a]=cos(a*2*PI/n);
		y[a]=sin(a*2*PI/n);
	 }
}

void EditorRender(int state)
{
  CEditor *I = &Editor;
  int i0,i1;
  float v[3],v0[3],v1[3],v2[3];
  float d0[3],n0[3],n1[3],n2[3];
  float x[50],y[50];
  int sele0,sele1;
  int nEdge;
  int c;
  float tube_size=0.5;


  if(state!=I->ActiveState)
    {
      EditorSetActiveObject(NULL,0);
    }

  if(I->Obj) {
    if(PMGUI) {

      nEdge = SettingGet(cSetting_stick_quality);
      if(nEdge>50)
        nEdge=50;

      subdivide(nEdge,x,y);

      sele0 = SelectorIndexByName(cEditorSele1);
      if(sele0>=0) {
        sele1 = SelectorIndexByName(cEditorSele2);
        if(sele1>=0) {
          /* bond mode */
          i0 = ObjectMoleculeGetAtomIndex(I->Obj,sele0); /* slow */
          i1 = ObjectMoleculeGetAtomIndex(I->Obj,sele1); /* slow */
          if((i0>=0)&&(i1>=0)) {
            ObjectMoleculeGetAtomVertex(I->Obj,state,i0,v0);
            ObjectMoleculeGetAtomVertex(I->Obj,state,i1,v1);

            subtract3f(v1,v0,d0);
            average3f(v1,v0,v2);
            copy3f(d0,n0);
            get_system1f3f(n0,n1,n2);

            glDisable(GL_LIGHTING);
            glColor3fv(ColorGet(0));
            glBegin(GL_LINE_LOOP);
            for(c=0;c<nEdge;c++) {
              v[0] = v2[0] + n1[0]*tube_size*x[c] + n2[0]*tube_size*y[c];
              v[1] = v2[1] + n1[1]*tube_size*x[c] + n2[1]*tube_size*y[c];
              v[2] = v2[2] + n1[2]*tube_size*x[c] + n2[2]*tube_size*y[c];
              glVertex3fv(v);
            }
            glEnd();
            glEnable(GL_LIGHTING);
          }
        } else {
          /* atom mode */
          i0 = ObjectMoleculeGetAtomIndex(I->Obj,sele0); /* slow */
          if(i0>=0) {
            ObjectMoleculeGetAtomVertex(I->Obj,state,i0,v2);
            n0[0]=1.0;
            n0[1]=0.0;
            n0[2]=0.0;
            get_system1f3f(n0,n1,n2);
            
            glDisable(GL_LIGHTING);
            glColor3fv(ColorGet(0));
            glBegin(GL_LINE_LOOP);
            for(c=0;c<nEdge;c++) {
              v[0] = v2[0] + n1[0]*tube_size*x[c] + n2[0]*tube_size*y[c];
              v[1] = v2[1] + n1[1]*tube_size*x[c] + n2[1]*tube_size*y[c];
              v[2] = v2[2] + n1[2]*tube_size*x[c] + n2[2]*tube_size*y[c];
              glVertex3fv(v);
            }
            glEnd();
            glBegin(GL_LINE_LOOP);
            for(c=0;c<nEdge;c++) {
              v[0] = v2[0] + n0[0]*tube_size*x[c] + n2[0]*tube_size*y[c];
              v[1] = v2[1] + n0[1]*tube_size*x[c] + n2[1]*tube_size*y[c];
              v[2] = v2[2] + n0[2]*tube_size*x[c] + n2[2]*tube_size*y[c];
              glVertex3fv(v);
            }
            glEnd();
            glBegin(GL_LINE_LOOP);
            for(c=0;c<nEdge;c++) {
              v[0] = v2[0] + n0[0]*tube_size*x[c] + n1[0]*tube_size*y[c];
              v[1] = v2[1] + n0[1]*tube_size*x[c] + n1[1]*tube_size*y[c];
              v[2] = v2[2] + n0[2]*tube_size*x[c] + n1[2]*tube_size*y[c];
              glVertex3fv(v);
            }
            glEnd();
            glEnable(GL_LIGHTING);
          }
          
        }
      }
    }
  }
}


void EditorSetActiveObject(ObjectMolecule *obj,int state)
{
  int sele1,sele2;

  CEditor *I = &Editor;
  if(obj) {
    I->Obj=obj;
    sele1 = SelectorIndexByName(cEditorSele1);
    if(sele1>=0) {
      sele2 = SelectorIndexByName(cEditorSele2);
      I->NFrag = SelectorSubdivideObject(cEditorFragPref,obj,
                                         sele1,sele2,cEditorBasePref);
      I->ActiveState=state;
    } else {
      I->Obj=NULL;
    }
  } else {
    SelectorSubdivideObject(cEditorFragPref,NULL,0,0,
                            cEditorBasePref);
    I->Obj=NULL;
  }
}

void EditorPrepareDrag(ObjectMolecule *obj,int index,int state)
{
  int frg;
  int sele0,sele1;
  int s;
  WordType name;
  int seleFlag= false;
  int i0,i1;
  CEditor *I = &Editor;
  
  if(!I->Obj) { /* non-anchored */
    I->DragObject=obj;
    I->DragIndex=index;
    I->DragSelection=-1;
  } else if(I->ActiveState!=state) {
    EditorSetActiveObject(NULL,0); /* recurses */
    return;
  } else { /* anchored */
    for(frg=0;frg<I->NFrag;frg++) {
      sprintf(name,"%s%02d",cEditorFragPref,frg);
      sele0=SelectorIndexByName(name);
      if(sele0>=0) {
        s=obj->AtomInfo[index].selEntry;
        while(s) {
          if(SelectorMatch(s,sele0)) {
            seleFlag=true;
            break;
          }
          s=SelectorNext(s);
        }
      }
      if(seleFlag)
        break;
    }
    if(seleFlag) {

      PRINTF " Editor: grabbing (%s).",name ENDF
      I->DragIndex = index;
      I->DragSelection = sele0;
      I->DragObject = obj;
      I->DragHaveAxis = false;
      I->DragHaveBase = false;
      I->DragBondFlag = false;

      sprintf(name,"%s%02d",cEditorBasePref,frg); /* get relevant base vertex of bond */
      sele1 = SelectorIndexByName(name);
      if(sele1>=0) {
        i1 = ObjectMoleculeGetAtomIndex(obj,sele1);
        if(i1>=0) {
          ObjectMoleculeGetAtomVertex(obj,state,i1,I->DragBase);
          I->DragHaveBase = true;
        }
      }

      /* get axis or base atom */

      sele0 = SelectorIndexByName(cEditorSele1);
      if(sele0>=0) {
        sele1 = SelectorIndexByName(cEditorSele2);
        I->DragSele0 = sele0;
        I->DragSele1 = sele1;
        if((sele0>=0)&&(sele1>=0)) { /* bond mode */
          I->DragBondFlag=true;
          i0 = ObjectMoleculeGetAtomIndex(obj,sele0);
          i1 = ObjectMoleculeGetAtomIndex(obj,sele1);
          if((i0>=0)&&(i1>=0)) {
            ObjectMoleculeGetAtomVertex(obj,state,i0,I->V0);
            ObjectMoleculeGetAtomVertex(obj,state,i1,I->V1);
            subtract3f(I->V1,I->V0,I->Axis);
            average3f(I->V1,I->V0,I->Center);
            normalize3f(I->Axis);
            I->DragHaveAxis=true;
          }
        } else {
          i0 = ObjectMoleculeGetAtomIndex(obj,sele0);
          if(i0>=0) 
            ObjectMoleculeGetAtomVertex(obj,state,i0,I->V0);      
          if(I->DragHaveBase) {
            copy3f(I->DragBase,I->V1)
            subtract3f(I->V1,I->V0,I->Axis);
            average3f(I->V1,I->V0,I->Center);
            normalize3f(I->Axis);
            I->DragHaveAxis=true;
          }
        }
      }
    } else {
      I->DragIndex=-1;
      I->DragSelection=-1;
      I->DragObject=NULL;
    }
  }
  if(I->DragObject)
    ObjectMoleculeSaveUndo(I->DragObject,state);
}

void EditorDrag(ObjectMolecule *obj,int index,int mode,int state,float *pt,float *mov)
{
  CEditor *I = &Editor;
  float v0[3],v1[3],v2[3],v3[3],cp[3];
  float d0[3],d1[3],d2[3],n0[3],n1[3],n2[3];
  float opp,adj,theta;
  float m[16];

  if((index=I->DragIndex)&&(obj=I->DragObject)) {
    if(obj!=I->Obj) {
      /* non-achored actions */
      switch(mode) {
      case cButModeRotFrag:
        SceneOriginGet(v3);

        subtract3f(pt,v3,n0);
        add3f(pt,mov,n1);
        subtract3f(n1,v3,n1);
        normalize3f(n0);
        normalize3f(n1);
        cross_product3f(n0,n1,cp);
        theta = asin(length3f(cp));
        normalize23f(cp,n2);        
        
        MatrixRotation44f(m,theta,n2[0],n2[1],n2[2]);
        m[3 ] = -v3[0];
        m[7 ] = -v3[1];
        m[11] = -v3[2];
        m[12] =  v3[0];
        m[13] =  v3[1];
        m[14] =  v3[2];
        ObjectMoleculeTransformSelection(obj,state,I->DragSelection,m);
        SceneDirty();
        break;
      case cButModeTorFrag:
        ObjectMoleculeMoveAtom(obj,state,index,mov,1);
        SceneDirty();
        break;
      case cButModeMovFrag:
        MatrixLoadIdentity44f(m);
        copy3f(mov,m+12);
        ObjectMoleculeTransformSelection(obj,state,I->DragSelection,m);
        SceneDirty();
        break;
      }
    } else {
      switch(mode) {
      case cButModeRotFrag:
        if(I->DragHaveBase&&I->DragBondFlag) {
          copy3f(I->DragBase,v3);
        } else {
          copy3f(I->V0,v3);
        }
        subtract3f(pt,v3,n0);
        add3f(pt,mov,n1);
        subtract3f(n1,v3,n1);
        normalize3f(n0);
        normalize3f(n1);
        cross_product3f(n0,n1,cp);
        theta = asin(length3f(cp));
        normalize23f(cp,n2);        
        
        MatrixRotation44f(m,theta,n2[0],n2[1],n2[2]);
        m[3 ] = -v3[0];
        m[7 ] = -v3[1];
        m[11] = -v3[2];
        m[12] =  v3[0];
        m[13] =  v3[1];
        m[14] =  v3[2];
        ObjectMoleculeTransformSelection(obj,state,I->DragSelection,m);
        SceneDirty();
        break;
      case cButModeTorFrag:
        if(I->DragHaveAxis) {
          subtract3f(pt,I->Center,d0);
          if(dot_product3f(d0,I->Axis)<0.0) {
            copy3f(I->V0,v1);
            copy3f(I->V1,v0);
          } else {
            copy3f(I->V0,v0);
            copy3f(I->V1,v1);
          }          
          subtract3f(v1,v0,d1);
          copy3f(d1,n0);
          normalize3f(n0);
          cross_product3f(n0,d0,n1);
          normalize3f(n1);
          
          project3f(d0,n0,v2);
          add3f(I->Center,v2,v2); /* v2 is the perpendicular point on the axis */
          subtract3f(pt,v2,d2);
          opp=length3f(mov);
          adj=length3f(d2);
          if(adj>R_SMALL4) {
            theta=atan(opp/adj);
            if(dot_product3f(n1,mov)<0.0)
              theta=-theta;
            MatrixRotation44f(m,theta,n0[0],n0[1],n0[2]);
            m[3 ] = -v1[0];
            m[7 ] = -v1[1];
            m[11] = -v1[2];
            m[12] =  v1[0];
            m[13] =  v1[1];
            m[14] =  v1[2];
            ObjectMoleculeTransformSelection(obj,state,I->DragSelection,m);
          }
        }
        SceneDirty();
        break;
      case cButModeMovFrag:
        MatrixLoadIdentity44f(m);
        copy3f(mov,m+12);
        ObjectMoleculeTransformSelection(obj,state,I->DragSelection,m);
        SceneDirty();
        break;
      }
    }
  }
}

void EditorInit(void)
{
  CEditor *I = &Editor;
  I->Obj = NULL;
  I->NFrag= 0;
  I->DragObject=NULL;
  I->DragIndex=-1;
  I->DragSelection=-1;
}

void EditorFree(void)
{
}

