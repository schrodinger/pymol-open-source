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
#include <GL/glut.h>

#include"Base.h"
#include"OOMac.h"
#include"Executive.h"
#include"ObjectMolecule.h"
#include"ListMacros.h"
#include"Ortho.h"
#include"Scene.h"
#include"Selector.h"
#include"Vector.h"
#include"Color.h"
#include"Setting.h"

#define cExecObject 0
#define cExecSelection 1

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

#define ExecOpCnt 2
#define ExecColorVisible 0.45,0.45,0.45
#define ExecColorHidden 0.3,0.3,0.3

void ExecutiveObjMolSeleOp(int sele,ObjectMoleculeOpRec *op);
SpecRec *ExecutiveFindSpec(char *name);
void ExecutiveInvalidateRep(char *name,int rep,int level);

/*========================================================================*/
void ExecutiveUpdateObjectSelection(struct Object *obj)
{
  SelectorDelete(obj->Name);  
  SelectorCreate(obj->Name,NULL,(ObjectMolecule*)obj); /* create a selection with same name */ 
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

  OrthoDoDraw();

  glutSwapBuffers();
}
/*========================================================================*/
void ExecutiveRay(void)
{
  SceneDoRay();
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
		op.i1 = 0.0;
		ExecutiveObjMolSeleOp(sele,&op);			 
		if(op.f1>0.0)
		  {
			 SceneWindowSphere(op.v1,op.f1);						  
			 SceneDirty();
		  }
	 }
  } 
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
  int flag = false;
  SpecRec *tRec;
  ObjectMoleculeOpRec op;
  tRec = ExecutiveFindSpec(name);
  if(tRec) {
	 if(name[0]=='_') {
		flag=true;
	 } else {
		if(state!=tRec->repOn[rep])  /* more intuitive behavior for REAL selections */
		  {
			 flag=true;
			 tRec->repOn[rep]=state;
		  }
	 }
	 
	 if(flag) {
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
  while(ListIterate(I->Spec,rec,next,SpecList))
	 {
		if(rec->type==cExecObject)
		  {
			 if(WordMatch(rec->obj->Name,name,true)) 
				{
				  SelectorDelete(rec->name);
				  rec->obj->fFree(rec->obj);
				  rec->obj=NULL;
				  ListDelete(I->Spec,rec,next,SpecList);
				  break;
				}
		  }
		else if(rec->type==cExecSelection)
		  {

			 if(WordMatch(rec->name,name,true))
				{
				  SelectorDelete(rec->name);
				  ListDelete(I->Spec,rec,next,SpecList);
				  break;
				}
		  }
	 }
}
/*========================================================================*/
void ExecutiveManageObject(Object *obj)
{
  int a;
  SpecRec *rec = NULL;
  CExecutive *I = &Executive;
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
  if(!rec)
	 ListElemAlloc(rec,SpecRec);
  SceneObjectAdd(obj);
  strcpy(rec->name,obj->Name);
  rec->type=cExecObject;
  rec->next=NULL;
  rec->obj=obj;
  rec->visible=1;
  for(a=0;a<cRepCnt;a++)
	 rec->repOn[a]=false;
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
  /*  Executive.Button = button;
  I->LastX = x;
  I->LastY = y;*/
  
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
			 t=(cRepCnt-1)-(I->Block->rect.right-(ExecRightMargin+x))/ExecToggleWidth;
			 if((t>=0)&&(t<cRepCnt))
				ExecutiveSetRepVisib(rec->name,t,!rec->repOn[t]);
			 else if((-t)<=ExecOpCnt)
				{
				  t=-t;
				  if(t==1)
					 ExecutiveCenter(rec->name,true);
				  else if(t==2) {
					 ExecutiveWindowZoom(rec->name);
					 ExecutiveCenter(rec->name,true);
				  }
				}
			 else if(rec->type==cExecObject)
				{
				  if(rec->visible)
					 SceneObjectDel(rec->obj);				
				  else 
					 SceneObjectAdd(rec->obj);
				  rec->visible=!rec->visible;
				}
		  }
		n--;
	 }
  glutPostRedisplay();
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
  SpecRec *rec = NULL;
  CExecutive *I = &Executive;
  glColor3fv(I->Block->BackColor);
  BlockFill(I->Block);
  glColor3fv(I->Block->TextColor);

  x = I->Block->rect.left+ExecLeftMargin;
  y = (I->Block->rect.top-ExecLineHeight)-ExecTopMargin;
  xx = I->Block->rect.right-ExecRightMargin-ExecToggleWidth*(cRepCnt+ExecOpCnt);
  
  while(ListIterate(I->Spec,rec,next,SpecList))
	 {
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
			 glutBitmapCharacter(GLUT_BITMAP_8_BY_13,'%');
			 c=rec->name;
		  }
		if(c)
		  while(*c) 
			 glutBitmapCharacter(GLUT_BITMAP_8_BY_13,*(c++));
		x2=xx;
		y2=y-ExecToggleMargin;
		for(a=0;a<ExecOpCnt;a++)
		  {
			 if(!a) {
				glBegin(GL_LINE_LOOP);
				glVertex2i(x2,y2+(ExecToggleSize-1)/2);
				glVertex2i(x2+(ExecToggleSize-1)/2,y2);
				glVertex2i(x2+ExecToggleSize-1,y2+(ExecToggleSize-1)/2);
				glVertex2i(x2+(ExecToggleSize-1)/2,y2+ExecToggleSize-1);
				glEnd();
			 } else {
				glBegin(GL_LINES);
				glVertex2i(x2,y2+(ExecToggleSize-1)/2);
				glVertex2i(x2+ExecToggleSize-1,y2+(ExecToggleSize-1)/2);
				glVertex2i(x2+(ExecToggleSize-1)/2,y2);
				glVertex2i(x2+(ExecToggleSize-1)/2,y2+ExecToggleSize-1);
				glEnd();				
			 }
	
			x2+=ExecToggleWidth;
		  }

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
		y-=ExecLineHeight;
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

  ListInit(I->Spec);
  I->Block = OrthoNewBlock(NULL);  
  I->Block->fRelease = ExecutiveRelease;
  I->Block->fClick   = ExecutiveClick;
  I->Block->fDrag    = ExecutiveDrag;
  I->Block->fDraw    = ExecutiveDraw;
  I->Block->fReshape = ExecutiveReshape;
  I->Block->active = true;

  OrthoAttach(I->Block,cOrthoTool);
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


