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
#include<stdio.h>
#include<GL/gl.h>
#include<GL/glut.h>

#include"main.h"
#include"Base.h"
#include"OOMac.h"
#include"Executive.h"
#include"ObjectMolecule.h"
#include"ObjectMesh.h"
#include"ListMacros.h"
#include"Ortho.h"
#include"Scene.h"
#include"Selector.h"
#include"Vector.h"
#include"Color.h"
#include"Setting.h"
#include"Matrix.h"
#include"PUtils.h"
#include"Menu.h"

#define cExecObject 0
#define cExecSelection 1
#define cExecAll 2

typedef struct SpecRec {
  int type;
  WordType  name; /*only used for selections*/
  struct Object *obj;  
  struct SpecRec *next;
  int repOn[cRepCnt];
  int visible;
} SpecRec; /* specification record (a line in the executive window) */

ListVarDeclare(SpecList,SpecRec);

typedef struct Executive {
  Block *Block;
  SpecRec *Spec;
  int Width,Height;
} CExecutive;

CExecutive Executive;

int ExecutiveClick(Block *block,int button,int x,int y,int mod);
int ExecutiveRelease(Block *block,int x,int y,int mod);
int ExecutiveDrag(Block *block,int x,int y,int mod);
void ExecutiveDraw(Block *block);
void ExecutiveReshape(Block *block,int width,int height);


#define ExecLineHeight 14
#define ExecTopMargin 0
#define ExecToggleMargin 2
#define ExecLeftMargin 5
#define ExecRightMargin 2
#define ExecToggleWidth 14
#define ExecToggleSize 13

#define ExecOpCnt 4
#define ExecColorVisible 0.45,0.45,0.45
#define ExecColorHidden 0.3,0.3,0.3

void ExecutiveObjMolSeleOp(int sele,ObjectMoleculeOpRec *op);
SpecRec *ExecutiveFindSpec(char *name);
void ExecutiveInvalidateRep(char *name,int rep,int level);

void ExecutiveSort(char *name)
{
  CExecutive *I = &Executive;
  Object *os=NULL;
  ObjectMolecule *obj;
  SpecRec *rec = NULL;
  ObjectMoleculeOpRec op;
  int sele;

  if(strlen(name)) {
    os=ExecutiveFindObjectByName(name);
    if(!os)
      ErrMessage(" Executive","object not found.");
    else if(os->type!=cObjectMolecule)
      ErrMessage(" Executive","bad object type.");
  }
  
  if(os||(!strlen(name))) { /* sort one or all */
    while(ListIterate(I->Spec,rec,next,SpecList)) {
      if(rec->type==cExecObject)
        if(rec->obj->type==cObjectMolecule)
          if((!os)||(rec->obj==os)) {
            obj =(ObjectMolecule*)rec->obj;
            ObjectMoleculeSort(obj);
            sele=SelectorIndexByName(rec->obj->Name);
            if(sele>=0) {
              op.code='INVA';
              op.i1=cRepAll; 
              op.i2=cRepInvAll;
              ExecutiveObjMolSeleOp(sele,&op);
            }
          }
    }
    SceneChanged();
  }
}
/*========================================================================*/
float ExecutiveOverlap(char *s1,int state1,char *s2,int state2)
{
  int sele1,sele2;
  
  sele1=SelectorIndexByName(s1);
  sele2=SelectorIndexByName(s2);
  
  return(SelectorSumVDWOverlap(sele1,state1,sele2,state2));
}
/*========================================================================*/
void ExecutiveStereo(int flag)
{

  if(PMGUI) {
    if(StereoCapable) {
      SceneSetStereo(flag);
    }
  }
}
/*========================================================================*/
float ExecutiveDistance(char *s1,char *s2)
{
  int sele1,sele2;
   char buffer[255];
   float dist = -1.0;

   ObjectMoleculeOpRec op1;
   ObjectMoleculeOpRec op2;

   sele1=SelectorIndexByName(s1);
   op1.i1=0;
   if(sele1>=0) {
     op1.code = 'SUMC';
     op1.v1[0]=0.0;
     op1.v1[1]=0.0;
     op1.v1[2]=0.0;
     ExecutiveObjMolSeleOp(sele1,&op1);
   } else {
     ErrMessage("ExecutiveDistance","The first selection contains no atoms.");
   }

   sele2=SelectorIndexByName(s2);
   op2.i1=0;
   if(sele2>=0) {
     op2.code = 'SUMC';
     op2.v1[0]=0.0;
     op2.v1[1]=0.0;
     op2.v1[2]=0.0;
     op2.i1=0;
     ExecutiveObjMolSeleOp(sele2,&op2);
   } else {
     ErrMessage("ExecutiveDistance","The second selection contains no atoms.");
   }
   
   if(op1.i1&&op2.i1) {
     scale3f(op1.v1,1.0/op1.i1,op1.v1);
     scale3f(op2.v1,1.0/op2.i1,op2.v1);
     dist = diff3f(op1.v1,op2.v1);
     sprintf(buffer," Distance: %8.3f [%i atom(s) to %i atom(s)]\n",
             dist,
            op1.i1,op2.i1);
    OrthoAddOutput(buffer);
  } else {
    ErrMessage("ExecutiveRMS","No atoms selected.");
	}
   return(dist);
}
/*========================================================================*/
char *ExecutiveSeleToPDBStr(char *s1,int state,int conectFlag)
{
  char *result=NULL;
  ObjectMoleculeOpRec op1;
  int sele1;

  sele1=SelectorIndexByName(s1);
  op1.charVLA=VLAlloc(char,10000);

  if(conectFlag) {
    if(state<0) state=0;
    op1.i2=SelectorGetPDB(&op1.charVLA,sele1,state,conectFlag);
  } else {

    op1.i2 = 0;
    op1.i3 = 0; /* atIndex */
    if(sele1>=0) {
      op1.code = 'PDB1';
      op1.i1 = state;
      ExecutiveObjMolSeleOp(sele1,&op1);
    }
  }
  VLACheck(op1.charVLA,char,op1.i2);
  op1.charVLA[op1.i2]=0;
  op1.i2++; /* get trailing null */
  result=Alloc(char,op1.i2);
  memcpy(result,op1.charVLA,op1.i2);
  VLAFreeP(op1.charVLA);
  return(result);
}
/*========================================================================*/
void ExecutiveCopy(char *src,char *dst)
{
  Object *os;
  ObjectMolecule *oSrc,*oDst;
  SpecRec *rec1 = NULL,*rec2=NULL;
  int a;

  os=ExecutiveFindObjectByName(src);
  if(!os)
    ErrMessage(" Executive","object not found.");
  else if(os->type!=cObjectMolecule)
    ErrMessage(" Executive","bad object type.");
  else 
    {
      oSrc =(ObjectMolecule*)os;
      oDst = ObjectMoleculeCopy(oSrc);
      if(oDst) {
        strcpy(oDst->Obj.Name,dst);
        ExecutiveManageObject((Object*)oDst);
        rec1=ExecutiveFindSpec(oSrc->Obj.Name);
        rec2=ExecutiveFindSpec(oDst->Obj.Name);
        if(rec1&&rec2) {
          for(a=0;a<cRepCnt;a++)
            rec2->repOn[a]=rec1->repOn[a];
        }
        ErrOk(" Executive","object created.");
      }
    }
  SceneChanged();
}

/*========================================================================*/
void ExecutiveOrient(char *sele,Matrix33d mi)
{
  double egval[3],egvali[3];
  double evect[3][3];
  float m[4][4];

  int a,b;

  if(!MatrixEigensolve33d((double*)mi,egval,egvali,(double*)evect)) {

	 normalize3d(evect[0]);
	 normalize3d(evect[1]);
	 normalize3d(evect[2]);


	 for(a=0;a<3;a++) {
		for(b=0;b<3;b++) {
		  m[a][b]=evect[b][a];
		}
	 }

    for(a=0;a<3;a++)
      {
        m[3][a]=0;
        m[a][3]=0;
      }
    m[3][3]=1.0;

    normalize3f(m[0]);
    normalize3f(m[1]);
    normalize3f(m[2]);

    SceneSetMatrix(m[0]);

    /* there must  be a more elegant to get the PC on X and the SC
     * on Y then what is shown below, but I couldn't get it to work.
     * I tried swapping the eigen-columns around but either that is 
     * a bogus approach (??) or my code was buggy.  Hence the following...*/

    if((egval[0]<egval[2])&&(egval[2]<egval[1])) { /* X < Z < Y */
      SceneRotate(90,1,0,0); /*1<-->2*/
    } else if((egval[1]<egval[0])&&(egval[0]<egval[2])) { /* Y < X < Z */
      SceneRotate(90,0,0,1); /*0<-->1*/
    } else if((egval[1]<egval[2])&&(egval[2]<egval[0])) { /* Y < Z < X */
      SceneRotate(90,0,1,0); /*1<-->2*/
      SceneRotate(90,0,0,1); /*0<-->1*/
    } else if((egval[2]<egval[1])&&(egval[1]<egval[0])) { /* Z < Y < X */
      SceneRotate(90,0,1,0); /*0<-->2*/
    } else if((egval[2]<egval[0])&&(egval[0]<egval[1])) { /* Z < X < Y */
      SceneRotate(90,0,1,0); /*0<-->2*/
      SceneRotate(90,1,0,0); /*0<-->1*/
    }
    /* X < Y < Z  - do nothing - that's what we want */


    ExecutiveWindowZoom(sele);
    ExecutiveCenter(sele,1);

  }
}
/*========================================================================*/
void ExecutiveAlter(char *s1,char *expr)
{
  int sele1;
  char buffer[255];
  ObjectMoleculeOpRec op1;
  
  sele1=SelectorIndexByName(s1);
  if(sele1>=0) {
    op1.code = 'ALTR';
    op1.s1 = expr;
    op1.i1 = 0;
    ExecutiveObjMolSeleOp(sele1,&op1);
    sprintf(buffer,"modified %i atoms.",op1.i1);
    ErrOk(" Alter",buffer);
  } else {
    ErrMessage("ExecutiveAlter","No atoms selected.");
  }
}
/*========================================================================*/
float ExecutiveRMS(char *s1,char *s2,int mode)
{
  int sele1,sele2;
  float rms = -1.0;
  int a;
  float inv,*f;
  ObjectMoleculeOpRec op1;
  ObjectMoleculeOpRec op2;
  OrthoLineType buffer;

  sele1=SelectorIndexByName(s1);
  op1.vv1=NULL;
  op1.vc1=NULL;
  op2.vv1=NULL;
  op2.vc1=NULL;
  if(sele1>=0) {
    op1.code = 'AVRT';
    op1.nvv1=0;
    op1.vc1=(int*)VLAMalloc(1000,sizeof(int),5,1);
    op1.vv1=(float*)VLAMalloc(1000,sizeof(float),5,1);
    ExecutiveObjMolSeleOp(sele1,&op1);
    for(a=0;a<op1.nvv1;a++)
      {
        inv=op1.vc1[a];
        if(inv)
          {
            f=op1.vv1+(a*3);
            inv=1.0/inv;
            *(f++)*=inv;
            *(f++)*=inv;
            *(f++)*=inv;
          }
      }
  }

  sele2=SelectorIndexByName(s2);
  if(sele2>=0) {
    op2.code = 'AVRT';
    op2.nvv1=0;
    op2.vc1=(int*)VLAMalloc(1000,sizeof(int),5,1);
    op2.vv1=(float*)VLAMalloc(1000,sizeof(float),5,1);
    ExecutiveObjMolSeleOp(sele2,&op2);
    for(a=0;a<op2.nvv1;a++)
      {
        inv=op2.vc1[a];
        if(inv)
          {
            f=op2.vv1+(a*3);
            inv=1.0/inv;
            *(f++)*=inv;
            *(f++)*=inv;
            *(f++)*=inv;
          }
      }
  }
  if(op1.vv1&&op2.vv1) {
    if(op1.nvv1!=op2.nvv1) {
      sprintf(buffer,"Atom counts between selections don't match (%d vs %d)\n",
              op1.nvv1,op2.nvv1);
      ErrMessage("ExecutiveRMS",buffer);
    } else if(op1.nvv1) {
      if(mode!=0) 
        rms = MatrixFitRMS(op1.nvv1,op1.vv1,op2.vv1,NULL,op2.ttt);
      else
        rms = MatrixGetRMS(op1.nvv1,op1.vv1,op2.vv1,NULL);
      printf(" Executive: RMS = %8.3f (%d to %d atoms)\n",
             rms,op1.nvv1,op2.nvv1);
      if(mode==2) {
        op2.code = 'TTTF';
        ExecutiveObjMolSeleOp(sele1,&op2);
      }
    } else {
      ErrMessage("ExecutiveRMS","No atoms selected.");
    }
  }
  VLAFreeP(op1.vv1);
  VLAFreeP(op2.vv1);
  VLAFreeP(op1.vc1);
  VLAFreeP(op2.vc1);
  return(rms);
}
/*========================================================================*/
float *ExecutiveRMSStates(char *s1,int target,int mode)
{
  int sele1;
  ObjectMoleculeOpRec op1;
  ObjectMoleculeOpRec op2;
  float *result = NULL;
  sele1=SelectorIndexByName(s1);
  if(sele1>=0) {
    op1.code = 'SVRT';
    op1.nvv1=0;
    op1.i1=target;
    op1.vv1=(float*)VLAMalloc(1000,sizeof(float),5,0);
    ExecutiveObjMolSeleOp(sele1,&op1);

    op2.vv2=op1.vv1;
    op2.nvv2=op1.nvv1;
    op2.i2=target;
    op2.i1=mode;
    op2.f1VLA=VLAlloc(float,10);
    VLASetSize(op2.f1VLA,0); /* failsafe */
    op2.vv1=(float*)VLAMalloc(1000,sizeof(float),5,0);
    op2.code = 'SFIT';
    op2.nvv1=0;
    ExecutiveObjMolSeleOp(sele1,&op2);
    result=op2.f1VLA;
  } 
  VLAFreeP(op1.vv1);
  VLAFreeP(op2.vv1);
  return(result);
}
/*========================================================================*/
float ExecutiveRMSPairs(WordType *sele,int pairs,int mode)
{
  int sele1,sele2;
  int a,c;
  float rms,inv,*f;
  OrthoLineType buffer;

  ObjectMoleculeOpRec op1;
  ObjectMoleculeOpRec op2;
  OrthoLineType combi,s1;

  op1.nvv1=0;
  op1.vc1=(int*)VLAMalloc(1000,sizeof(int),5,1);
  op1.vv1=(float*)VLAMalloc(1000,sizeof(float),5,1); /* auto-zero */
  op1.code = 'AVRT';

  op2.nvv1=0;
  op2.vc1=(int*)VLAMalloc(1000,sizeof(int),5,1);
  op2.vv1=(float*)VLAMalloc(1000,sizeof(float),5,1); /* auto-zero */
  op2.code = 'AVRT';

  strcpy(combi,"(");
  c=0;
  for(a=0;a<pairs;a++) {
    sele1=SelectorIndexByName(sele[c]);
    if(sele1>=0) ExecutiveObjMolSeleOp(sele1,&op1);
    strcat(combi,sele[c]);
    if(a<(pairs-1)) strcat(combi," or ");
    c++;
    sele2=SelectorIndexByName(sele[c]);
    if(sele2>=0) ExecutiveObjMolSeleOp(sele2,&op2);
    c++;
  }
  strcat(combi,")");
  for(a=0;a<op1.nvv1;a++)
    {
      inv=op1.vc1[a];
      if(inv)
        {
          f=op1.vv1+(a*3);
          inv=1.0/inv;
          *(f++)*=inv;
          *(f++)*=inv;
          *(f++)*=inv;
        }
    }
  for(a=0;a<op2.nvv1;a++)
    {
      inv=op2.vc1[a];
      if(inv)
        {
          f=op2.vv1+(a*3);
          inv=1.0/inv;
          *(f++)*=inv;
          *(f++)*=inv;
          *(f++)*=inv;
        }
    }
  if(op1.vv1&&op2.vv1) {
    if(op1.nvv1!=op2.nvv1) {
      sprintf(buffer,"Atom counts between selection sets don't match (%d != %d).",
              op1.nvv1,op2.nvv1);
      ErrMessage("ExecutiveRMS",buffer);
    } else if(op1.nvv1) {
      if(mode!=0)
        rms = MatrixFitRMS(op1.nvv1,op1.vv1,op2.vv1,NULL,op2.ttt);
      else
        rms = MatrixGetRMS(op1.nvv1,op1.vv1,op2.vv1,NULL);
      printf(" ExecutiveRMS: RMS = %8.3f (%d to %d atoms)\n",
             rms,op1.nvv1,op2.nvv1);
      op2.code = 'TTTF';
      SelectorGetTmp(combi,s1);
      sele1=SelectorIndexByName(s1);
      ExecutiveObjMolSeleOp(sele1,&op2);
      SelectorFreeTmp(s1);
    } else {
      ErrMessage("ExecutiveRMS","No atoms selected.");
    }
  }
  VLAFreeP(op1.vv1);
  VLAFreeP(op2.vv1);
  VLAFreeP(op1.vc1);
  VLAFreeP(op2.vc1);
  return(rms);
}
/*========================================================================*/
void ExecutiveUpdateObjectSelection(struct Object *obj)
{
  if(obj->type==cObjectMolecule) {
    SelectorDelete(obj->Name);  
    SelectorCreate(obj->Name,NULL,(ObjectMolecule*)obj); /* create a selection with same name */ 
  }
}
/*========================================================================*/
void ExecutiveReset(int cmd)
{
  SceneResetMatrix();
  SelectorCreate("_all","all",NULL); /* replace this with a persistent selection? */
  ExecutiveCenter("_all",1);
  ExecutiveWindowZoom("_all");
  ExecutiveDelete("_all"); 
}
/*========================================================================*/
void ExecutiveDrawNow(void) 
{
  glMatrixMode(GL_MODELVIEW);

  /*  glClear( GL_DEPTH_BUFFER_BIT);*/

  SceneUpdate();

  OrthoDoDraw();

  MainSwapBuffers();

}
/*========================================================================*/
int ExecutiveCountStates(char *s1)
{
  int sele1;
  int result=0;
  sele1=SelectorIndexByName(s1);
  if(sele1>=0) {
    SelectorUpdateTable();
    result = SelectorGetSeleNCSet(sele1);
  }
  return(result);
}


/*========================================================================*/
void ExecutiveRay(void)
{
  SceneRay();
}
/*========================================================================*/
void ExecutiveSetSetting(char *sname,char *v)
{
  SettingSetNamed(sname,v);
}
/*========================================================================*/
void ExecutiveColor(char *name,char *color,int flags)
{
  SpecRec *rec = NULL;
  int sele;
  ObjectMoleculeOpRec op;
  
  /* per atom */
  if(!(flags&0x1)) {
	 sele=SelectorIndexByName(name);
	 if(sele>=0) {
		op.code = 'COLR';
		op.i1=ColorGetIndex(color);
		ExecutiveObjMolSeleOp(sele,&op);
		op.code='INVA';
		op.i1=cRepAll; 
		op.i2=cRepInvColor;
		ExecutiveObjMolSeleOp(sele,&op);
	 } else {
      ErrMessage("Color","Invalid selection.");
    }
  }
  /* per object */
  rec=ExecutiveFindSpec(name);
  if(rec) {
	 if(rec->type==cExecObject) {
		rec->obj->Color=ColorGetIndex(color);
		SceneDirty();
	 }
  }
}
/*========================================================================*/
SpecRec *ExecutiveFindSpec(char *name)
{
  CExecutive *I = &Executive;
  SpecRec *rec = NULL;
  while(ListIterate(I->Spec,rec,next,SpecList)) {
	 if(strcmp(rec->name,name)==0) 
		break;
  }
  return(rec);
}
/*========================================================================*/
void ExecutiveObjMolSeleOp(int sele,ObjectMoleculeOpRec *op) 
{
  CExecutive *I=&Executive;
  SpecRec *rec = NULL;
  ObjectMolecule *obj = NULL;

  if(sele>=0)
	 {
		while(ListIterate(I->Spec,rec,next,SpecList))
		  {
			 if(rec->type==cExecObject)
				{
				  if(rec->obj->type==cObjectMolecule)
					 {
						obj=(ObjectMolecule*)rec->obj;
						ObjectMoleculeSeleOp(obj,sele,op);
					 }
				}
		  }
	 }
}
/*========================================================================*/
void ExecutiveWindowZoom(char *name)
{
  int sele;
  ObjectMoleculeOpRec op;

  sele=SelectorIndexByName(name);
  if(sele>=0) {
	 op.code = 'SUMC';
	 op.v1[0]=0.0;
	 op.v1[1]=0.0;
	 op.v1[2]=0.0;
	 op.i1=0;
	 
	 ExecutiveObjMolSeleOp(sele,&op);
	 
	 if(op.i1) {
		scale3f(op.v1,1.0/op.i1,op.v1);
		op.code = 'MDST';
		op.f1 = 0.0;
		op.i1 = 0.0;
		ExecutiveObjMolSeleOp(sele,&op);			 
		op.f1 += MAX_VDW;
		if(op.f1>0.0)
		  {
          SceneWindowSphere(op.v1,op.f1);						  
          SceneDirty();
		  }
	 }
  } 
}
/*========================================================================*/
void ExecutiveGetBBox(char *name,float *mn,float *mx)
{
  int sele;
  int c;
  ObjectMoleculeOpRec op;

  op.i1=0;

  sele=SelectorIndexByName(name);
  if(sele>=0) {
	 op.code = 'MNMX';
	 
	 ExecutiveObjMolSeleOp(sele,&op);
	 
  }
  if(!op.i1) {
    for(c=0;c<3;c++) {
      mn[c]=-1.0;
      mx[c+3]=1.0;
    }
  } else {
    for(c=0;c<3;c++) {
      mn[c]=op.v1[c];
      mx[c]=op.v2[c];
    }
  }
}
/*========================================================================*/
int ExecutiveGetMoment(char *name,Matrix33d mi)
{
  int sele;
  ObjectMoleculeOpRec op;
  int a,b;
  int c=0;
  for(a=0;a<3;a++)
	 {
		for(b=0;b<3;b++)
		  mi[a][b]=0.0;
		mi[a][a]=1.0;
	 }
  
  sele=SelectorIndexByName(name);
  if(sele>=0) {
	 op.code = 'SUMC';
	 op.v1[0]=0.0;
	 op.v1[1]=0.0;
	 op.v1[2]=0.0;
	 op.i1=0;
	 
	 ExecutiveObjMolSeleOp(sele,&op);
	 
	 if(op.i1) {
		c+=op.i1;
		scale3f(op.v1,1.0/op.i1,op.v1);
		op.code = 'MOME';		
		for(a=0;a<3;a++)
		  for(b=0;b<3;b++)
			 op.d[a][b]=0.0;
		ExecutiveObjMolSeleOp(sele,&op);			 
		for(a=0;a<3;a++)
		  for(b=0;b<3;b++)
			 mi[a][b]=op.d[a][b];
	 }
  } 
  return(c);
}
/*========================================================================*/
void ExecutiveCenter(char *name,int preserve)
{
  int sele;
  ObjectMoleculeOpRec op;
  sele=SelectorIndexByName(name);
  if(sele>=0) {
	 op.code = 'SUMC';
	 op.v1[0]=0.0;
	 op.v1[1]=0.0;
	 op.v1[2]=0.0;
	 op.i1=0;
	 
	 ExecutiveObjMolSeleOp(sele,&op);
	 
	 if(op.i1) {
		scale3f(op.v1,1.0/op.i1,op.v1);
		SceneOriginSet(op.v1,preserve);
		SceneDirty();
	 }
  }
}
/*========================================================================*/
void ExecutiveSetObjVisib(char *name,int state)
{
  SpecRec *tRec;
  tRec = ExecutiveFindSpec(name);
  if(tRec) {
	 if(tRec->type==cExecObject)
		if(tRec->visible!=state)
		  {
			 if(tRec->visible)
				SceneObjectDel(tRec->obj);				
			 else 
				SceneObjectAdd(tRec->obj);
			 tRec->visible=!tRec->visible;
		  }
  }
}
/*========================================================================*/
void ExecutiveSetAllVisib(int state)
{
  ObjectMoleculeOpRec op;
  ObjectMolecule *obj;
  int rep;
  int sele;
  CExecutive *I = &Executive;
  SpecRec *rec = NULL;
  while(ListIterate(I->Spec,rec,next,SpecList)) {
	 if(rec->type==cExecObject)
		{
		  if(rec->obj->type==cObjectMolecule)
			 {
				obj=(ObjectMolecule*)rec->obj;
				sele = SelectorIndexByName(obj->Obj.Name);
				for(rep=0;rep<cRepCnt;rep++) {
				  rec->repOn[rep]=state;
				  op.code='VISI';
				  op.i1=rep;
				  op.i2=state;
				  ObjectMoleculeSeleOp(obj,sele,&op);
				  op.code='INVA';
				  op.i2=cRepInvVisib;
				  ObjectMoleculeSeleOp(obj,sele,&op);				
				}
			 }
		}
  }
}
/*========================================================================*/
void ExecutiveSetRepVisib(char *name,int rep,int state)
{
  int sele;
  int a;
  SpecRec *tRec;
  ObjectMoleculeOpRec op;
  tRec = ExecutiveFindSpec(name);
  if(tRec) {
	 if(name[0]!='_') {
      /* remember visibility information for real selections */
	   if(rep>=0) {
        tRec->repOn[rep]=state;
      } else {
        for(a=0;a<cRepCnt;a++)
          tRec->repOn[a]=state; 
      }
	 }
    sele=SelectorIndexByName(name);
    if(sele>=0) {
      op.code='VISI';
      op.i1=rep;
      op.i2=state;
      ExecutiveObjMolSeleOp(sele,&op);
      op.code='INVA';
      op.i2=cRepInvVisib;
      ExecutiveObjMolSeleOp(sele,&op);
	 }
  }
}
/*========================================================================*/
void ExecutiveInvalidateRep(char *name,int rep,int level)
{
  int sele = -1;
  ObjectMoleculeOpRec op;
  sele=SelectorIndexByName(name);
  if(sele>=0) {
	 op.code = 'INVA';
	 op.i1=rep;
	 op.i2=level;
	 ExecutiveObjMolSeleOp(sele,&op);
  }
}
/*========================================================================*/
Object *ExecutiveFindObjectByName(char *name)
{
  CExecutive *I = &Executive;
  SpecRec *rec = NULL;
  Object *obj=NULL;
  while(ListIterate(I->Spec,rec,next,SpecList))
	 {
		if(rec->type==cExecObject)
		  {
			 if(strcmp(rec->obj->Name,name)==0) 
				{
				  obj=rec->obj;
				  break;
				}
		  }
	 }
  return(obj);
}
/*========================================================================*/
Block *ExecutiveGetBlock(void)
{
  CExecutive *I = &Executive;
  return(I->Block);
}
/*========================================================================*/
void ExecutiveSetControlsOff(char *name)
{
  SpecRec *rec;
  int a;
  rec = ExecutiveFindSpec(name);
  if(rec)
	 {
		for(a=0;a<cRepCnt;a++)
		  rec->repOn[a]=false;
	 }
}
/*========================================================================*/
void ExecutiveDelete(char *name)
{
  CExecutive *I = &Executive;
  SpecRec *rec = NULL;
  int all_flag=false;
  if(WordMatch(name,"all",true)<0) all_flag=true;
  while(ListIterate(I->Spec,rec,next,SpecList))
	 {
		if(rec->type==cExecObject)
		  {
			 if(all_flag||(WordMatch(rec->obj->Name,name,true)<0))
				{
				  SelectorDelete(rec->name);
				  rec->obj->fFree(rec->obj);
				  rec->obj=NULL;
				  ListDelete(I->Spec,rec,next,SpecList);
				  rec=NULL;
				}
		  }
		else if(rec->type==cExecSelection)
		  {

			 if(all_flag||(WordMatch(rec->name,name,true)<0))
				{
				  SelectorDelete(rec->name);
				  ListDelete(I->Spec,rec,next,SpecList);
				  rec=NULL;
				}
		  }
	 }
}
/*========================================================================*/
void ExecutiveDump(char *fname,char *obj)
{
  SpecRec *rec = NULL;
  CExecutive *I = &Executive;

  SceneUpdate();

  while(ListIterate(I->Spec,rec,next,SpecList))
	 {
		if(rec->type==cExecObject)
		  {
			 if(strcmp(rec->obj->Name,obj)==0) 
				break;
		  }
	 }
  if(rec)
	 { 
      if(rec->obj->type==cObjectMesh) {
        ObjectMeshDump((ObjectMesh*)rec->obj,fname);
      } else {
        ErrMessage("ExecutiveDump","Invalid object type for this operation.");
      }
	 }
  else {
    ErrMessage("ExecutiveDump","Object not found.");
  }
  
}
/*========================================================================*/
void ExecutiveManageObject(Object *obj)
{
  int a;
  SpecRec *rec = NULL;
  CExecutive *I = &Executive;
  char buffer[255];
  while(ListIterate(I->Spec,rec,next,SpecList))
	 {
		if(rec->type==cExecObject)
		  {
			 if(strcmp(rec->obj->Name,obj->Name)==0) 
				break;
		  }
	 }
  if(rec) /* another object of this type already exists */
	 { /* purge it */
		SceneObjectDel(rec->obj);
		rec->obj->fFree(rec->obj);
		rec->obj=NULL;
	 }
  else 
    {
      sprintf(buffer," Executive: object \"%s\" created.\n",obj->Name);
      OrthoAddOutput(buffer);
    }
  if(!rec)
	 ListElemAlloc(rec,SpecRec);
  strcpy(rec->name,obj->Name);
  rec->type=cExecObject;
  rec->next=NULL;
  rec->obj=obj;
  if(rec->obj->type==cObjectMap) {
    rec->visible=0;
  } else {
    rec->visible=1;
    SceneObjectAdd(obj);
  }
  for(a=0;a<cRepCnt;a++)
	 rec->repOn[a]=false;
  if(rec->obj->type==cObjectMolecule)
    rec->repOn[cRepLine]=true;
  ListAppend(I->Spec,rec,next,SpecList);
  if(rec->obj->type==cObjectMolecule)
	 ExecutiveUpdateObjectSelection(obj);
  ExecutiveCenter(obj->Name,false);
  ExecutiveWindowZoom(obj->Name);
}
/*========================================================================*/
void ExecutiveManageSelection(char *name)
{
  int a;
  SpecRec *rec = NULL;
  CExecutive *I = &Executive;
  ListElemAlloc(rec,SpecRec);
  strcpy(rec->name,name);
  rec->type=cExecSelection;
  rec->next=NULL;
  for(a=0;a<cRepCnt;a++)
	 rec->repOn[a]=false;
  ListAppend(I->Spec,rec,next,SpecList);
}
/*========================================================================*/
int ExecutiveClick(Block *block,int button,int x,int y,int mod)
{
  CExecutive *I = &Executive;
  int n,a;
  SpecRec *rec = NULL;
  int t;
  Object *obj;

  n=((I->Block->rect.top-(y+2))-ExecTopMargin)/ExecLineHeight;
  a=n;

  while(ListIterate(I->Spec,rec,next,SpecList))
	 {
		if(!a)
		  {
			 t = ((I->Block->rect.right-ExecRightMargin)-x)/ExecToggleWidth;
          if(t<ExecOpCnt) {
            y = I->Block->rect.top-(ExecTopMargin + n*ExecLineHeight) - 1;
            x = I->Block->rect.right-(ExecRightMargin + t*ExecToggleWidth);
            t = (ExecOpCnt-t)-1;
            switch(t) {
            case 2:
              switch(rec->type) {
              case cExecAll:
              case cExecSelection:
                MenuActivate(x,y,"mol_show",rec->name);
                break;
              case cExecObject:
                switch(rec->obj->type) {
                case cObjectMolecule:
                  MenuActivate(x,y,"mol_show",rec->obj->Name);
                  break;
                }
                break;
              }
              break;
            case 3:
              switch(rec->type) {
              case cExecAll:
              case cExecSelection:
                MenuActivate(x,y,"mol_hide",rec->name);
                break;
              case cExecObject:
                switch(rec->obj->type) {
                case cObjectMolecule:
                  MenuActivate(x,y,"mol_hide",rec->obj->Name);
                  break;
                }
                break;
              }
              break;
            }
          }
        }
      a--;
    }
  MainDirty();
  
  return(1);
}
/*========================================================================*/
int ExecutiveRelease(Block *block,int x,int y,int mod)
{
  CExecutive *I = &Executive;
  int n;  
  SpecRec *rec = NULL;
  int t;

  n=((I->Block->rect.top-(y+2))-ExecTopMargin)/ExecLineHeight;

  while(ListIterate(I->Spec,rec,next,SpecList))
	 {
		if(!n)
		  {
			 t = ((I->Block->rect.right-ExecRightMargin)-x)/ExecToggleWidth;
          if(t<ExecOpCnt) {
            t = (ExecOpCnt-t)-1;
            switch(t) {
            case 0:
              switch(rec->type) {
              case cExecAll:
                SelectorCreate("_all","all",NULL);
                ExecutiveCenter("_all",1);
                ExecutiveWindowZoom("_all");
                ExecutiveDelete("_all"); 
                break;
              default:
                ExecutiveCenter(rec->name,true);
                ExecutiveWindowZoom(rec->name);
              }
              break;
            case 1:
              switch(rec->type) {
              case cExecAll:
                SelectorCreate("_all","all",NULL);
                ExecutiveCenter("_all",1);
                ExecutiveDelete("_all"); 
                break;
              default:
                ExecutiveCenter(rec->name,true);
                break;
              }
              break;
            }
          } else if(rec->type==cExecObject)
				{
				  if(rec->visible)
					 SceneObjectDel(rec->obj);				
				  else 
					 SceneObjectAdd(rec->obj);
				  rec->visible=!rec->visible;
              SceneChanged();
				}
        }
		n--;
	 }
  MainDirty();
  return(1);
}
/*========================================================================*/
int ExecutiveDrag(Block *block,int x,int y,int mod)
{
  return(1);
}
/*========================================================================*/
void ExecutiveDraw(Block *block)
{
  int a,x,y,xx,x2,y2;
  char *c=NULL;
  float toggleColor[3] = { 0.5, 0.5, 1.0 };
  float toggleColor2[3] = { 0.3, 0.3, 0.6 };
  SpecRec *rec = NULL;
  CExecutive *I = &Executive;

  if(PMGUI) {
    glColor3fv(I->Block->BackColor);
    BlockFill(I->Block);
    
    x = I->Block->rect.left+ExecLeftMargin;
    y = (I->Block->rect.top-ExecLineHeight)-ExecTopMargin;
    /*    xx = I->Block->rect.right-ExecRightMargin-ExecToggleWidth*(cRepCnt+ExecOpCnt);*/
    xx = I->Block->rect.right-ExecRightMargin-ExecToggleWidth*(ExecOpCnt);
    
    while(ListIterate(I->Spec,rec,next,SpecList))
      {
        x2=xx;
        y2=y-ExecToggleMargin;
        glColor3fv(toggleColor);
        for(a=0;a<ExecOpCnt;a++)
          {
            switch(a) {
            case 0:
              glBegin(GL_LINE_LOOP);
              glVertex2i(x2,y2+(ExecToggleSize-1)/2);
              glVertex2i(x2+(ExecToggleSize-1)/2,y2);
              glVertex2i(x2+ExecToggleSize-1,y2+(ExecToggleSize-1)/2);
              glVertex2i(x2+(ExecToggleSize-1)/2,y2+ExecToggleSize-1);
              glEnd();
              break;
            case 1:
              glBegin(GL_LINES);
              glVertex2i(x2,y2+(ExecToggleSize-1)/2);
              glVertex2i(x2+ExecToggleSize-1,y2+(ExecToggleSize-1)/2);
              glVertex2i(x2+(ExecToggleSize-1)/2,y2);
              glVertex2i(x2+(ExecToggleSize-1)/2,y2+ExecToggleSize-1);
              glEnd();				
              break;
            case 2:
              glBegin(GL_POLYGON);
              glVertex2i(x2,y2);
              glVertex2i(x2,y2+ExecToggleSize);
              glVertex2i(x2+ExecToggleSize,y2+ExecToggleSize);
              glVertex2i(x2+ExecToggleSize,y2);
              glEnd();
              glColor3f(0.0,0.0,0.0);
              glRasterPos4d((double)(x2+2),(double)(y2+2),0.0,1.0);
              glutBitmapCharacter(GLUT_BITMAP_8_BY_13,'S');              
              glColor3fv(toggleColor);

              break;
            case 3:
              glColor3fv(toggleColor2);
              glBegin(GL_POLYGON);
              glVertex2i(x2,y2);
              glVertex2i(x2,y2+ExecToggleSize);
              glVertex2i(x2+ExecToggleSize,y2+ExecToggleSize);
              glVertex2i(x2+ExecToggleSize,y2);
              glEnd();
              glColor3f(0.0,0.0,0.0);
              glRasterPos4d((double)(x2+2),(double)(y2+2),0.0,1.0);
              glutBitmapCharacter(GLUT_BITMAP_8_BY_13,'H');              
              glColor3fv(toggleColor);
              break;
            }
            x2+=ExecToggleWidth;
          }
        
#ifdef PYMOL_OLD_CODE 
        for(a=0;a<cRepCnt;a++)
          {
            if(rec->repOn[a]) 
              {
                glBegin(GL_POLYGON);
                glVertex2i(x2,y2);
                glVertex2i(x2,y2+ExecToggleSize);
                glVertex2i(x2+ExecToggleSize,y2+ExecToggleSize);
                glVertex2i(x2+ExecToggleSize,y2);
              }
            else
              {
                glBegin(GL_LINE_LOOP);
                glVertex2i(x2,y2);
                glVertex2i(x2,y2+ExecToggleSize-1);
                glVertex2i(x2+ExecToggleSize-1,y2+ExecToggleSize-1);
                glVertex2i(x2+ExecToggleSize-1,y2);
              }
            glEnd();
            x2+=ExecToggleWidth;
          }
#endif

        glColor3fv(I->Block->TextColor);
        glRasterPos4d((double)(x),(double)(y),0.0,1.0);
        if(rec->type==cExecObject)
          {
            y2=y-ExecToggleMargin;
            if(rec->visible)
              glColor3f(ExecColorVisible);
            else
              glColor3f(ExecColorHidden);
            glBegin(GL_POLYGON);
            glVertex2i(x-ExecToggleMargin,y2);
            glVertex2i(xx-ExecToggleMargin,y2);
            glVertex2i(xx-ExecToggleMargin,y2+ExecToggleSize);
            glVertex2i(x-ExecToggleMargin,y2+ExecToggleSize);
            glEnd();
            c=rec->obj->Name;
            glColor3fv(I->Block->TextColor);
		
          }
        else if(rec->type==cExecSelection)
          {
            glutBitmapCharacter(GLUT_BITMAP_8_BY_13,'(');
            c=rec->name;
          }
        else if(rec->type==cExecAll)
          {
            c=rec->name;
          }

        if(c)
          while(*c) 
            glutBitmapCharacter(GLUT_BITMAP_8_BY_13,*(c++));

        if(rec->type==cExecSelection)
          {
            glutBitmapCharacter(GLUT_BITMAP_8_BY_13,')');
            c=rec->name;
          }

        y-=ExecLineHeight;
      }
  }
}
/*========================================================================*/
int ExecutiveIterateObject(Object **obj,void **hidden)
{
  int result;
  CExecutive *I = &Executive;
  int flag=false;
  SpecRec **rec=(SpecRec**)hidden;
  while(!flag)
	 {
		result = (int)ListIterate(I->Spec,(*rec),next,SpecList);
		if(!(*rec))
		  flag=true;
		else if((*rec)->type==cExecObject)
		  flag=true;
	 }
  if(*rec)
	 (*obj)=(*rec)->obj;
  else
	 (*obj)=NULL;
  return(result);
}
/*========================================================================*/
void ExecutiveReshape(Block *block,int width,int height)
{
  CExecutive *I = &Executive;

  BlockReshape(block,width,height);

  I->Width = block->rect.right-block->rect.left+1;
  I->Height = block->rect.top-block->rect.bottom+1;
  
}
/*========================================================================*/
void ExecutiveInit(void)
{
  CExecutive *I = &Executive;
  SpecRec *rec = NULL;
  int a;

  ListInit(I->Spec);
  I->Block = OrthoNewBlock(NULL);  
  I->Block->fRelease = ExecutiveRelease;
  I->Block->fClick   = ExecutiveClick;
  I->Block->fDrag    = ExecutiveDrag;
  I->Block->fDraw    = ExecutiveDraw;
  I->Block->fReshape = ExecutiveReshape;
  I->Block->active = true;

  OrthoAttach(I->Block,cOrthoTool);

  ListElemAlloc(rec,SpecRec);
  strcpy(rec->name,"(all)");
  rec->type=cExecAll;
  rec->next=NULL;
  for(a=0;a<cRepCnt;a++)
	 rec->repOn[a]=false;
  ListAppend(I->Spec,rec,next,SpecList);

}
/*========================================================================*/
void ExecutiveFree(void)
{
  CExecutive *I = &Executive;
  SpecRec *rec=NULL;
  while(ListIterate(I->Spec,rec,next,SpecList))
	 {
		if(rec->type==cExecObject)
		  rec->obj->fFree(rec->obj);
	 }
  ListFree(I->Spec,next,SpecList);
  OrthoFreeBlock(I->Block);
  I->Block=NULL;
}


#ifdef _undefined

matrix checking code...

		double mt[3][3],mt2[3][3],pr[3][3],im[3][3],em[3][3];
	 printf("normalized matrix \n");
	 for(a=0;a<3;a++) {
		for(b=0;b<3;b++) {
		  printf("%12.3f ",evect[a][b]);
		}
		printf("\n");
	 }
	 printf("\n");

	 printf("tensor \n");
	 for(a=0;a<3;a++) {
		for(b=0;b<3;b++) {
		  printf("%12.3f ",mi[a][b]);
		}
		printf("\n");
	 }
	 printf("\n");

	  for(a=0;a<3;a++) {
		 for(b=0;b<3;b++) {
			mt[a][b]=evect[a][b];
		 }
	  }

	 for(a=0;a<3;a++) {
		for(b=0;b<3;b++) {
		  mt2[a][b]=evect[b][a];
		}
	 }

	 matrix_multiply33d33d(mt,mt2,pr);
	 printf("self product 1 \n");
	 for(a=0;a<3;a++) {
		for(b=0;b<3;b++) {
		  printf("%8.3f ",pr[a][b]);
		}
		printf("\n");
	 }
	 printf("\n");

	 matrix_multiply33d33d(mt,mi,im);
	 matrix_multiply33d33d(im,mt2,pr);
	 printf("diagonal product 1 \n");
	 for(a=0;a<3;a++) {
		for(b=0;b<3;b++) {
		  printf("%8.3f ",pr[a][b]);
		}
		printf("\n");
	 }
	 printf("\n");

	 for(a=0;a<3;a++) {
		for(b=0;b<3;b++) {
		  em[a][b]=0.0;
		}
		em[a][a]=egval[a];
	 }

	 matrix_multiply33d33d(mt2,em,im);
	 matrix_multiply33d33d(im,mt,pr);
	 printf("diagonal product 4 \n");
	 for(a=0;a<3;a++) {
		for(b=0;b<3;b++) {
		  printf("%8.3f ",pr[a][b]);
		}
		printf("\n");
	 }
	 printf("\n");
#endif

