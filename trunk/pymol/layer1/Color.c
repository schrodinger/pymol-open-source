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

#include"Base.h"
#include"OOMac.h"
#include"MemoryDebug.h"
#include"Ortho.h"
#include"Word.h"

#include"Color.h"
#include"PConv.h"
#include"ObjectGadgetRamp.h"
#include"Util.h"
#include"Executive.h"
#include"MyPNG.h"
#include"Scene.h"

#define cAutoColorMask 0x7
static int AutoColor[] = { 26, 5, 154, 6, 9, 29, 11, 30 };

CColor Color;

#define cColorExtCutoff (-10)

int ColorGetNext(void) 
{
  int result;
  int next;
  next = (int)SettingGet(cSetting_auto_color_next);

  next = (next&cAutoColorMask);
  result = AutoColor[next];
  next++;
  next = (next&cAutoColorMask);
  SettingSet(cSetting_auto_color_next,(float)next);
  return(result);
}

int ColorCheckRamped(int index)
{
  return(index<=(cColorExtCutoff));
}

int ColorGetRamped(int index,float *vertex,float *color)
{
  CColor *I=&Color;
  int ok=false;
  if(index<=cColorExtCutoff) {
    index = cColorExtCutoff - index;
    if(index<I->NExt) {
      if(!I->Ext[index].Ptr) {
        if(I->Ext[index].Name)
          I->Ext[index].Ptr = (void*)ExecutiveFindObjectByName(I->Ext[index].Name);
      }
      if(I->Ext[index].Ptr) 
        ok = ObjectGadgetRampInterVertex(
                                         (ObjectGadgetRamp*)I->Ext[index].Ptr,
                                         vertex,color);
    }
  
  }
  if(!ok) {
    color[0]=1.0;
    color[1]=1.0;
    color[2]=1.0;
  }
  return(ok);
}

static int ColorFindExtByName(char *name,int null_okay) 
{
  CColor *I=&Color;
  int result = -1;
  int wm;
  int best;
  int a;
  best = 0;
  for(a=0;a<I->NExt;a++)
    {
      wm = WordMatch(name,I->Ext[a].Name,true);
      if(wm<0) {
        if(null_okay||(!I->Ext[a].Ptr)) {
          result=a;
          break;
        }
      } else if ((wm>0)&&(best<wm)) {
        if(null_okay||(!I->Ext[a].Ptr)) {
          result=a;
          best=wm;
        }
      }
    }
  return(result);
}

void ColorRegisterExt(char *name,void *ptr,int type)
{
  CColor *I=&Color;
  int a;
  a=ColorFindExtByName(name,true);
  if(a<0) {
    VLACheck(I->Ext,ExtRec,I->NExt);
    a = I->NExt;
    I->NExt++;
  }
  if(a>=0) {
    UtilNCopy(I->Ext[a].Name,name,sizeof(ColorName));
    I->Ext[a].Ptr = ptr;
    I->Ext[a].Type = type;
  }
}

void ColorForgetExt(char *name)
{
  CColor *I=&Color;
  int a;
  a=ColorFindExtByName(name,true);

  /* this won't work! */

  if(a>=0) { /* currently leaks memory TODO fix */
    I->Ext[a].Name[0]=0;
    I->Ext[a].Ptr=NULL;
  }
}

PyObject *ColorExtAsPyList(void)
{
  CColor *I=&Color;
  PyObject *result,*list;
  ExtRec *ext;
  int a;

  result = PyList_New(I->NExt);
  ext=I->Ext;
  for(a=0;a<I->NExt;a++) {
    list = PyList_New(2);
    PyList_SetItem(list,0,PyString_FromString(ext->Name));
    PyList_SetItem(list,1,PyInt_FromLong(ext->Type));
    PyList_SetItem(result,a,list);
    ext++;
  }
  return(result);
}

/*========================================================================*/
PyObject *ColorAsPyList()
{
  CColor *I=&Color;
  PyObject *result,*list;
  ColorRec *color;
  int n_custom=0;
  int a,c;
  color=I->Color;
  for(a=0;a<I->NColor;a++) {
    if(color->Custom||color->ClampedFlag)
      n_custom++;
    color++;
  }
  result = PyList_New(n_custom);
  c=0;
  color=I->Color;
  for(a=0;a<I->NColor;a++) {
    if(color->Custom||color->ClampedFlag) {
      list = PyList_New(6);
      PyList_SetItem(list,0,PyString_FromString(color->Name));
      PyList_SetItem(list,1,PyInt_FromLong(a));
      PyList_SetItem(list,2,PConvFloatArrayToPyList(color->Color,3));
      PyList_SetItem(list,3,PyInt_FromLong(color->Custom));
      PyList_SetItem(list,4,PyInt_FromLong(color->ClampedFlag));
      PyList_SetItem(list,5,PConvFloatArrayToPyList(color->Clamped,3));
      PyList_SetItem(result,c,list);
      c++;
    }
    color++;
  }
  return(result);
}

int ColorExtFromPyList(PyObject *list)
{
  int n_ext=0;
  int a;
  int ok=true;
  int ll;
  CColor *I=&Color;
  PyObject *rec;
  ExtRec *ext;
  if(ok) ok=(list!=NULL);
  if(ok) ok=PyList_Check(list);

  if(ok) ll = PyList_Size(list);
  /* TO SUPPORT BACKWARDS COMPATIBILITY...
   Always check ll when adding new PyList_GetItem's */

  if(ok) {
    n_ext=PyList_Size(list);
    VLACheck(I->Ext,ExtRec,n_ext); 
    ext=I->Ext;
    for(a=0;a<n_ext;a++) {
      rec=PyList_GetItem(list,a);
      if(ok) ok=(rec!=NULL);
      if(ok) ok=PyList_Check(rec);
      if(ok) ok=PConvPyStrToStr(PyList_GetItem(rec,0),ext->Name,sizeof(ColorName));
      if(ok) ok=PConvPyIntToInt(PyList_GetItem(rec,1),&ext->Type);
      ext++;
    }
    if(ok) I->NExt=n_ext;

  }
  return(ok);
}

/*========================================================================*/
int ColorFromPyList(PyObject *list)
{
  int n_custom=0;
  int a;
  int index=0;
  int ok=true;
  int ll;
  CColor *I=&Color;
  PyObject *rec;
  ColorRec *color;
  if(ok) ok=(list!=NULL);
  if(ok) ok=PyList_Check(list);
  if(ok) {
    n_custom=PyList_Size(list);
    for(a=0;a<n_custom;a++) {
      rec=PyList_GetItem(list,a);
      if(ok) ok=(rec!=NULL);
      if(ok) ok=PyList_Check(rec);
      if(ok) ll = PyList_Size(list);
      /* TO SUPPORT BACKWARDS COMPATIBILITY...
         Always check ll when adding new PyList_GetItem's */
      if(ok) ok=PConvPyIntToInt(PyList_GetItem(rec,1),&index);
      if(ok) {
        if(index>=I->NColor) {
          VLACheck(I->Color,ColorRec,index); /* auto-zeros */
          I->NColor=index+1;
        }
        color=I->Color+index;
        if(ok) ok=PConvPyStrToStr(PyList_GetItem(rec,0),color->Name,sizeof(ColorName));
        if(ok) ok=PConvPyListToFloatArrayInPlace(PyList_GetItem(rec,2),color->Color,3);
        if(PyList_Size(rec)>=6) {
          if(ok) ok=PConvPyIntToInt(PyList_GetItem(rec,3),&color->Custom);
          if(ok) ok=PConvPyIntToInt(PyList_GetItem(rec,4),&color->ClampedFlag);
          if(ok) ok=PConvPyListToFloatArrayInPlace(PyList_GetItem(rec,5),color->Clamped,3);
        } else {
          if(ok) {color->Custom=true;}
        }
      }
      if(!ok) break;
    }
  }
  return(ok);
}

/*========================================================================*/
void ColorDef(char *name,float *v)
{
  CColor *I=&Color;
  int color=-1;
  int a;
  int best;
  int wm;

  best = 0;
  for(a=0;a<I->NColor;a++)
	 {
      wm = WordMatch(name,I->Color[a].Name,true);
      if(wm<0) {
        color=a;
        break;
      } else if ((wm>0)&&(best<wm)) {
        color=a;
        best=wm;
      }
	 }
  if(color<0) {
    color=I->NColor;
    VLACheck(I->Color,ColorRec,I->NColor);
    I->NColor++;
  }
  strcpy(I->Color[color].Name,name);
  I->Color[color].Color[0]=v[0];
  I->Color[color].Color[1]=v[1];
  I->Color[color].Color[2]=v[2];
  I->Color[color].Custom=true;
  ColorUpdateClamp(color);

  PRINTFB(FB_Executive,FB_Actions)
    " Color: \"%s\" defined as [ %3.1f, %3.1f, %3.1f ].\n",name,v[0],v[1],v[2] 
    ENDFB;

}
/*========================================================================*/
int ColorGetIndex(char *name)
{
  CColor *I=&Color;
  int color=-1; /* default for unknown is white */
  int a;
  int i;
  int wm,best=0;


  if(((name[0]>='0')&&(name[0]<='9'))||(name[0]=='-'))
    if(sscanf(name,"%d",&i)) 
      if((i<I->NColor)&&(i>=0))
        return(i);
  if(WordMatch(name,"default",true))
    return(-1);
  for(a=0;a<I->NColor;a++)
	 {
      wm = WordMatch(name,I->Color[a].Name,true);
      if(wm<0) {
        color=a;
        break;
      } else if ((wm>0)&&(best<wm)) {
        color=a;
        best=wm;
      }
	 }
  if(color<0) {
    color = ColorFindExtByName(name,false);
    if(color>=0)
      color = -10-color; /* indicates external */
  }
  return(color);
}
/*========================================================================*/
float *ColorGetNamed(char *name)
{
  
  return(ColorGet(ColorGetIndex(name)));
}
/*========================================================================*/
char *ColorGetName(int index)
{
  CColor *I=&Color;
  if((index>=0)&&(index<I->NColor))
    return(I->Color[index].Name);
  else
    return(NULL);
}
/*========================================================================*/
int ColorGetStatus(int index)
{
  CColor *I=&Color; /* return 0 if color is invalid, hidden; 1 otherwise */
  char *c;
  int result=0;
  if((index>=0)&&(index<I->NColor)) {
    c=I->Color[index].Name;
    result=1;
    while(*c) {
      if(((*c)>='0')&&((*c)<='9')) {
        result=0;
        break;
      }
      c++;
    }
  }
  return(result);
}
/*========================================================================*/
int ColorGetNColor(void)
{
  CColor *I=&Color;
  return(I->NColor);
}
/*========================================================================*/
void ColorFree(void)
{
  CColor *I=&Color;
  if(I->ColorTable) {
    FreeP(I->ColorTable);
  }
  VLAFreeP(I->Color);
  VLAFreeP(I->Ext);
}

/*========================================================================*/
void ColorReset(void)
{
  CColor *I=&Color;
  int a;
  int set1;
  float f;
  float spectrumS[13][3] = { 
    { 1.0, 0.0, 1.0 }, /* magenta - 0 */
    { 0.5, 0.0, 1.0 },
    { 0.0, 0.0, 1.0 }, /* blue - 166.66  */
    { 0.0, 0.5, 1.0 },
    { 0.0, 1.0, 1.0 }, /* cyan - 333.33 */

    { 0.0, 1.0, 0.5 },
    { 0.0, 1.0, 0.0 }, /* green - 500 */
    { 0.5, 1.0, 0.0 },
    { 1.0, 1.0, 0.0 }, /* yellow - 666.66 */
    { 1.0, 0.5, 0.0 },

    { 1.0, 0.0, 0.0 }, /* red - 833.33 */
    { 1.0, 0.0, 0.5 },
    { 1.0, 0.0, 1.0 }, /* magenta - 999 */
  };

  float spectrumR[13][3] = { 
    { 1.0, 1.0, 0.0 }, /* yellow - 0 */
    { 0.5, 1.0, 0.0 }, /* yellowgreen */
    { 0.0, 1.0, 0.0 }, /* green - 166.66 */
    { 0.0, 1.0, 0.5 }, /* lime */
    { 0.0, 1.0, 1.0 }, /* cyan - 333.33 */

    { 0.0, 0.5, 1.0 }, /* marine */
    { 0.0, 0.0, 1.0 }, /* blue - 500 */
    { 0.5, 0.0, 1.0 }, /* blueviolet */
    { 1.0, 0.0, 1.0 }, /* magenta - 666.66 */
    { 1.0, 0.0, 0.5 }, /* hotpink */

    { 1.0, 0.0, 0.0 }, /* red - 833.33 */
    { 1.0, 0.5, 0.0 }, /* orange */
    { 1.0, 1.0, 0.0 }, /* yellow - 999 */
  };

  float spectrumC[][3] = {
    { 1.0, 1.0, 0.0 }, /* yellow - 0 */
    { 0.0, 0.0, 1.0 }, /* blue - 83.333 */
    { 1.0, 0.0, 0.0 }, /* red - 167.67*/
    { 0.0, 1.0, 0.0 }, /* green - 250.00 */
    { 1.0, 0.0, 1.0 }, /* magenta - 333.33 */

    { 0.0, 1.0, 1.0 }, /* cyan - 416.67 */
    { 1.0, 1.0, 0.0 }, /* yellow - 500.00*/    
    { 0.0, 1.0, 0.0 }, /* green - 583.33*/    
    { 0.0, 0.0, 1.0 }, /* blue - 666.67 */    
    { 1.0, 0.0, 1.0 }, /* magenta - 750.00*/

    { 1.0, 1.0, 0.0 }, /* yellow - 833.33*/    
    { 1.0, 0.0, 0.0 }, /* red - 916.67*/    
    { 0.0, 1.0, 1.0 }, /* cyan - 999 */
  };

  float spectrumW[][3] = {
    { 1.0, 1.0, 0.0 }, /* yellow - 0*/
    { 1.0, 1.0, 1.0 }, /* white */
    { 0.0, 0.0, 1.0 }, /* blue  - 83.333 */
    { 1.0, 1.0, 1.0 }, /* white */
    { 1.0, 0.0, 0.0 }, /* red - 166.67 */

    { 1.0, 1.0, 1.0 }, /* white */
    { 0.0, 1.0, 0.0 }, /* green - 250.00 */
    { 1.0, 1.0, 1.0 }, /* white */
    { 1.0, 0.0, 1.0 }, /* magenta - 333.33 */
    { 1.0, 1.0, 1.0 }, /* white */

    { 0.0, 1.0, 1.0 }, /* cyan - 416.67 */
    { 1.0, 1.0, 1.0 }, /* white */
    { 1.0, 1.0, 0.0 }, /* yellow - 500.00*/    
    { 1.0, 1.0, 1.0 }, /* white */
    { 0.0, 1.0, 0.0 }, /* green - 583.33*/    

    { 1.0, 1.0, 1.0 }, /* white */
    { 0.0, 0.0, 1.0 }, /* blue - 666.67 */    
    { 1.0, 1.0, 1.0 }, /* white */
    { 1.0, 0.0, 1.0 }, /* magenta - 750.00*/
    { 1.0, 1.0, 1.0 }, /* white */

    { 1.0, 1.0, 0.0 }, /* yellow - 833.33 */    
    { 1.0, 1.0, 1.0 }, /* white */
    { 1.0, 0.0, 0.0 }, /* red - 916.67*/
    { 1.0, 1.0, 1.0 }, /* white */
    { 0.0, 1.0, 1.0 }, /* cyan - 999 */
  };

  /* BLUE->VIOLET->RED r546 to r909 */
  /* BLUE->CYAN->GREEN->YELLOW->RED s182 to s909 */
  /* BLUE->WHITE->RED w00 to */

  I->NColor=0;

  strcpy(I->Color[I->NColor].Name,"white");
  I->Color[I->NColor].Color[0]=1.0F;
  I->Color[I->NColor].Color[1]=1.0F;
  I->Color[I->NColor].Color[2]=1.0F;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"black");
  I->Color[I->NColor].Color[0]=0.0F;
  I->Color[I->NColor].Color[1]=0.0F;
  I->Color[I->NColor].Color[2]=0.0F;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"blue");
  I->Color[I->NColor].Color[0]=0.0F;
  I->Color[I->NColor].Color[1]=0.0F;
  I->Color[I->NColor].Color[2]=1.0F;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"green");
  I->Color[I->NColor].Color[0]=0.0F;
  I->Color[I->NColor].Color[1]=1.0F;
  I->Color[I->NColor].Color[2]=0.0F;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"red");
  I->Color[I->NColor].Color[0]=1.0F;
  I->Color[I->NColor].Color[1]=0.0F;
  I->Color[I->NColor].Color[2]=0.0F;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"cyan");
  I->Color[I->NColor].Color[0]=0.0F;
  I->Color[I->NColor].Color[1]=1.0F;
  I->Color[I->NColor].Color[2]=1.0F;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"yellow");
  I->Color[I->NColor].Color[0]=1.0F;
  I->Color[I->NColor].Color[1]=1.0F;
  I->Color[I->NColor].Color[2]=0.0F;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"dash");
  I->Color[I->NColor].Color[0]=1.0F;
  I->Color[I->NColor].Color[1]=1.0F;
  I->Color[I->NColor].Color[2]=0.0F;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"magenta");
  I->Color[I->NColor].Color[0]=1.0F;
  I->Color[I->NColor].Color[1]=0.0F;
  I->Color[I->NColor].Color[2]=1.0F;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"salmon");
  I->Color[I->NColor].Color[0]=1.0F;
  I->Color[I->NColor].Color[1]=0.6F; /* was 0.5 */
  I->Color[I->NColor].Color[2]=0.6F; /* wat 0.5 */
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"lime");
  I->Color[I->NColor].Color[0]=0.5F;
  I->Color[I->NColor].Color[1]=1.0F;
  I->Color[I->NColor].Color[2]=0.5F;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"slate");
  I->Color[I->NColor].Color[0]=0.5F;
  I->Color[I->NColor].Color[1]=0.5F;
  I->Color[I->NColor].Color[2]=1.0F;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"hotpink");
  I->Color[I->NColor].Color[0]=1.0F;
  I->Color[I->NColor].Color[1]=0.0F;
  I->Color[I->NColor].Color[2]=0.5F;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"orange");
  I->Color[I->NColor].Color[0]=1.0F;
  I->Color[I->NColor].Color[1]=0.5F;
  I->Color[I->NColor].Color[2]=0.0F;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"yellowgreen"); /* AKA puke green */
  I->Color[I->NColor].Color[0]=0.5F;
  I->Color[I->NColor].Color[1]=1.0F;
  I->Color[I->NColor].Color[2]=0.0F;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"bluegreen");
  I->Color[I->NColor].Color[0]=0.0F;
  I->Color[I->NColor].Color[1]=1.0F;
  I->Color[I->NColor].Color[2]=0.5F;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"blueviolet");
  I->Color[I->NColor].Color[0]=0.5F;
  I->Color[I->NColor].Color[1]=0.0F;
  I->Color[I->NColor].Color[2]=1.0F;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"marine");
  I->Color[I->NColor].Color[0]=0.0F;
  I->Color[I->NColor].Color[1]=0.5F;
  I->Color[I->NColor].Color[2]=1.0F;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"olive");
  I->Color[I->NColor].Color[0]=0.5F;
  I->Color[I->NColor].Color[1]=0.5F;
  I->Color[I->NColor].Color[2]=0.1F;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"purple");
  I->Color[I->NColor].Color[0]=0.7F;
  I->Color[I->NColor].Color[1]=0.2F;
  I->Color[I->NColor].Color[2]=0.7F;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"teal");
  I->Color[I->NColor].Color[0]=0.2F;
  I->Color[I->NColor].Color[1]=0.6F;
  I->Color[I->NColor].Color[2]=0.6F;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"ruby");
  I->Color[I->NColor].Color[0]=0.5F;
  I->Color[I->NColor].Color[1]=0.1F;
  I->Color[I->NColor].Color[2]=0.1F;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"forest");
  I->Color[I->NColor].Color[0]=0.1F;
  I->Color[I->NColor].Color[1]=0.5F;
  I->Color[I->NColor].Color[2]=0.1F;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"deep");
  I->Color[I->NColor].Color[0]=0.1F;
  I->Color[I->NColor].Color[1]=0.1F;
  I->Color[I->NColor].Color[2]=0.5F;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"grey");
  I->Color[I->NColor].Color[0]=0.5F;
  I->Color[I->NColor].Color[1]=0.5F;
  I->Color[I->NColor].Color[2]=0.5F;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"gray"); /* for the poor spellers */
  I->Color[I->NColor].Color[0]=0.5F;
  I->Color[I->NColor].Color[1]=0.5F;
  I->Color[I->NColor].Color[2]=0.5F;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"carbon");
  I->Color[I->NColor].Color[0]=0.2F;
  I->Color[I->NColor].Color[1]=1.0F;
  I->Color[I->NColor].Color[2]=0.2F;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"nitrogen");
  I->Color[I->NColor].Color[0]=0.2F;
  I->Color[I->NColor].Color[1]=0.2F;
  I->Color[I->NColor].Color[2]=1.0F;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"oxygen");
  I->Color[I->NColor].Color[0]=1.0F;
  I->Color[I->NColor].Color[1]=0.3F;
  I->Color[I->NColor].Color[2]=0.3F;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"hydrogen");
  I->Color[I->NColor].Color[0]=0.9F;
  I->Color[I->NColor].Color[1]=0.9F;
  I->Color[I->NColor].Color[2]=0.9F;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"brightorange");
  I->Color[I->NColor].Color[0]=1.0F;
  I->Color[I->NColor].Color[1]=0.7F;
  I->Color[I->NColor].Color[2]=0.2F;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"sulfur");
  I->Color[I->NColor].Color[0]=1.0F;
  I->Color[I->NColor].Color[1]=0.5F;
  I->Color[I->NColor].Color[2]=0.0F;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"tv_red");
  I->Color[I->NColor].Color[0]=1.0F;
  I->Color[I->NColor].Color[1]=0.2F;
  I->Color[I->NColor].Color[2]=0.2F;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"tv_green");
  I->Color[I->NColor].Color[0]=0.2F;
  I->Color[I->NColor].Color[1]=1.0F;
  I->Color[I->NColor].Color[2]=0.2F;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"tv_blue");
  I->Color[I->NColor].Color[0]=0.3F;
  I->Color[I->NColor].Color[1]=0.3F;
  I->Color[I->NColor].Color[2]=1.0F;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"tv_yellow");
  I->Color[I->NColor].Color[0]=1.0F;
  I->Color[I->NColor].Color[1]=1.0F;
  I->Color[I->NColor].Color[2]=0.2F;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"tv_yellow");
  I->Color[I->NColor].Color[0]=1.0F;
  I->Color[I->NColor].Color[1]=1.0F;
  I->Color[I->NColor].Color[2]=0.1F;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"tv_orange");
  I->Color[I->NColor].Color[0]=1.0F;
  I->Color[I->NColor].Color[1]=0.55F;
  I->Color[I->NColor].Color[2]=0.15F;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"br0");
  I->Color[I->NColor].Color[0]=0.1F;
  I->Color[I->NColor].Color[1]=0.1F;
  I->Color[I->NColor].Color[2]=1.0F;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"br1");
  I->Color[I->NColor].Color[0]=0.2F;
  I->Color[I->NColor].Color[1]=0.1F;
  I->Color[I->NColor].Color[2]=0.9F;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"br2");
  I->Color[I->NColor].Color[0]=0.3F;
  I->Color[I->NColor].Color[1]=0.1F;
  I->Color[I->NColor].Color[2]=0.8F;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"br3");
  I->Color[I->NColor].Color[0]=0.4F;
  I->Color[I->NColor].Color[1]=0.1F;
  I->Color[I->NColor].Color[2]=0.7F;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"br4");
  I->Color[I->NColor].Color[0]=0.5F;
  I->Color[I->NColor].Color[1]=0.1F;
  I->Color[I->NColor].Color[2]=0.6F;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"br5");
  I->Color[I->NColor].Color[0]=0.6F;
  I->Color[I->NColor].Color[1]=0.1F;
  I->Color[I->NColor].Color[2]=0.5F;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"br6");
  I->Color[I->NColor].Color[0]=0.7F;
  I->Color[I->NColor].Color[1]=0.1F;
  I->Color[I->NColor].Color[2]=0.4F;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"br7");
  I->Color[I->NColor].Color[0]=0.8F;
  I->Color[I->NColor].Color[1]=0.1F;
  I->Color[I->NColor].Color[2]=0.3F;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"br8");
  I->Color[I->NColor].Color[0]=0.9F;
  I->Color[I->NColor].Color[1]=0.1F;
  I->Color[I->NColor].Color[2]=0.2F;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"br9");
  I->Color[I->NColor].Color[0]=1.0F;
  I->Color[I->NColor].Color[1]=0.1F;
  I->Color[I->NColor].Color[2]=0.1F;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"pink");
  I->Color[I->NColor].Color[0]=1.0F;
  I->Color[I->NColor].Color[1]=0.65F;
  I->Color[I->NColor].Color[2]=0.85F;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"firebrick");
  I->Color[I->NColor].Color[0]=0.698F;
  I->Color[I->NColor].Color[1]=0.13F;
  I->Color[I->NColor].Color[2]=0.13F;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"chocolate");
  I->Color[I->NColor].Color[0]=0.555F;
  I->Color[I->NColor].Color[1]=0.222F;
  I->Color[I->NColor].Color[2]=0.111F;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"brown");
  I->Color[I->NColor].Color[0]=0.555F;
  I->Color[I->NColor].Color[1]=0.274F;
  I->Color[I->NColor].Color[2]=0.150F;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"wheat");
  I->Color[I->NColor].Color[0]=0.99F;
  I->Color[I->NColor].Color[1]=0.82F;
  I->Color[I->NColor].Color[2]=0.65F;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"violet");
  I->Color[I->NColor].Color[0]=1.0F;
  I->Color[I->NColor].Color[1]=0.5F;
  I->Color[I->NColor].Color[2]=1.0F;
  I->NColor++;

  /* greybow */

  for(a=0;a<100;a=a+1) {
    sprintf(I->Color[I->NColor].Name,"grey%02d",a);
    I->Color[I->NColor].Color[0]=a/99.0F;
    I->Color[I->NColor].Color[1]=a/99.0F;
    I->Color[I->NColor].Color[2]=a/99.0F;
    I->NColor++;
  }

  strcpy(I->Color[I->NColor].Name,"lightmagenta");
  I->Color[I->NColor].Color[0]=1.0F;
  I->Color[I->NColor].Color[1]=0.2F;
  I->Color[I->NColor].Color[2]=0.8F;
  I->NColor++;

  #define A_DIV 83.333333333F

  /* full spectrum (s000-s999) */

  for(a=0;a<1000;a=a+1) {
    set1=(int)(a/A_DIV);
    sprintf(I->Color[I->NColor].Name,"s%03d",a);
    f = 1.0F-(a-(set1*A_DIV))/A_DIV;
    I->Color[I->NColor].Color[0]=f*spectrumS[set1][0]+(1.0F-f)*spectrumS[set1+1][0];
    I->Color[I->NColor].Color[1]=f*spectrumS[set1][1]+(1.0F-f)*spectrumS[set1+1][1];
    I->Color[I->NColor].Color[2]=f*spectrumS[set1][2]+(1.0F-f)*spectrumS[set1+1][2];
    I->NColor++;
  }

  /* offset & reversed full spectrum (r000-r999) */

  for(a=0;a<1000;a=a+1) {
    set1=(int)(a/A_DIV);
    sprintf(I->Color[I->NColor].Name,"r%03d",a);
    f = 1.0F-(a-(set1*A_DIV))/A_DIV;
    I->Color[I->NColor].Color[0]=f*spectrumR[set1][0]+(1.0F-f)*spectrumR[set1+1][0];
    I->Color[I->NColor].Color[1]=f*spectrumR[set1][1]+(1.0F-f)*spectrumR[set1+1][1];
    I->Color[I->NColor].Color[2]=f*spectrumR[set1][2]+(1.0F-f)*spectrumR[set1+1][2];
    I->NColor++;
  }

  /* complementary spectra (c000-c999) */

  for(a=0;a<1000;a=a+1) {
    set1=(int)(a/A_DIV);
    sprintf(I->Color[I->NColor].Name,"c%03d",a);
    f = 1.0F-(a-(set1*A_DIV))/A_DIV;
    I->Color[I->NColor].Color[0]=f*spectrumC[set1][0]+(1.0F-f)*spectrumC[set1+1][0];
    I->Color[I->NColor].Color[1]=f*spectrumC[set1][1]+(1.0F-f)*spectrumC[set1+1][1];
    I->Color[I->NColor].Color[2]=f*spectrumC[set1][2]+(1.0F-f)*spectrumC[set1+1][2];
    I->NColor++;
  }

  #define W_DIV 41.666666667F

  /* complementary spectra separated by white (w000-w999) */

  for(a=0;a<1000;a=a+1) {
    set1=(int)(a/W_DIV);
    sprintf(I->Color[I->NColor].Name,"w%03d",a);
    f = 1.0F-(a-(set1*W_DIV))/W_DIV;
    I->Color[I->NColor].Color[0]=f*spectrumW[set1][0]+(1.0F-f)*spectrumW[set1+1][0];
    I->Color[I->NColor].Color[1]=f*spectrumW[set1][1]+(1.0F-f)*spectrumW[set1+1][1];
    I->Color[I->NColor].Color[2]=f*spectrumW[set1][2]+(1.0F-f)*spectrumW[set1+1][2];
    I->NColor++;
  }

  for(a=0;a<I->NColor;a++) { 
    /* mark all current colors non-custom so that they don't get saved in session files */
    I->Color[a].Custom=false;
  }

  I->NExt = 0;

}

int ColorTableLoad(char *fname,int quiet)
{
  CColor *I=&Color; 
  int ok=true;
  int width=512,height=512;
  unsigned int *table = NULL;

  if(!strcmp(fname,"rgb")) {
    FreeP(I->ColorTable);
    PRINTFB(FB_Color,FB_Actions)
      " Color: purged table; restoring RGB colors.\n"
      ENDFB;
    ColorUpdateClamp(-1);    
    
  } else if(!strcmp(fname,"pymol")) {
    
    int x,y;
    unsigned int r=0,g=0,b=0;
    unsigned int *pixel,mask,*p;
    unsigned int rc,bc,gc;
    unsigned int gf,bf,rf;

    float green_max=0.75F;
    float red_max=0.95F;
    float blue_max=0.97F;
    
    float min_factor=0.15F;

    red_max = SettingGet(cSetting_pymol_space_max_red);
    green_max = SettingGet(cSetting_pymol_space_max_green);
    blue_max = SettingGet(cSetting_pymol_space_max_blue);
    min_factor = SettingGet(cSetting_pymol_space_min_factor);

    FreeP(I->ColorTable);
    if(I->BigEndian)
      mask = 0x000000FF;
    else
      mask = 0xFF000000;
    
    table = Alloc(unsigned int,512*512);
    
    p=(unsigned int*)table; 
    for(x=0;x<width;x++)
      for(y=0;y<height;y++)
        *(p++)=mask;
    
    for(y=0;y<height;y++) 
      for(x=0;x<width;x++) {
        rc = r;
        gc = g;
        bc = b;

        if((r>=g)&&(r>=b)) {
          if(rc>255*red_max) {
            rc=(unsigned int)(red_max*255);
            bc=bc*rc/r;
            gc=gc*rc/r;
          }
        } else if((g>=b)&&(g>=r)) {
          if(gc>255*green_max) {
            gc=(unsigned int)(green_max*255);
            bc=bc*gc/g;
            rc=rc*gc/g;
          }
        } else if((b>=g)&&(b>=r)) {
          if(bc>255*blue_max) {
            bc=(unsigned int)(blue_max*255);
            gc=gc*bc/b;
            rc=rc*bc/b;
          }
        }

        rf = (int)(min_factor*rc+0.49999F);
        gf = (int)(min_factor*gc+0.49999F);
        bf = (int)(min_factor*bc+0.49999F);
        
        if (rc<gf) rc = (int)gf;
        if (bc<gf) bc = (int)gf;
        
        if (rc<bf) rc = (int)bf;
        if (gc<bf) gc = (int)bf;
        
        if (gc<rf) gc = (int)rf;
        if (bc<rf) bc = (int)rf;
    
        if(rc>255) rc=255;
        if(bc>255) bc=255;
        if(gc>255) gc=255;

        pixel = table+((width)*y)+x;
        if(I->BigEndian) {
          *(pixel)=
            mask|(rc<<24)|(gc<<16)|(bc<<8);
        } else {
          *(pixel)=
            mask|(bc<<16)|(gc<<8)|rc;
        }
        b = b + 4;
        if(!(0xFF&b)) { 
          b=0;
          g=g+4;
          if(!(0xFF&g)) {           
            g=0;
            r=r+4;
          }
        }
      }
  
    I->ColorTable = table;
    if(!quiet) {
      PRINTFB(FB_Color,FB_Actions)
        " Color: defined table '%s'.\n",fname
        ENDFB;
    }
    
    ColorUpdateClamp(-1);
    ExecutiveInvalidateRep(cKeywordAll,cRepAll,cRepInvColor);
    SceneChanged();

  } else {
    if(strlen(fname)) {
      if(MyPNGRead(fname,
                   (unsigned char**)&table,
                   (unsigned int*)&width,
                   (unsigned int*)&height)) {
        if((width==512)&&(height==512)) {
          FreeP(I->ColorTable);
          I->ColorTable = table;
          if(!quiet) {
            PRINTFB(FB_Color,FB_Actions)
              " Color: loaded table '%s'.\n",fname
              ENDFB;
          }
          
          ColorUpdateClamp(-1);

        } else {
          PRINTFB(FB_Color,FB_Errors)
            " ColorTableLoad-Error: invalid dimensions w x h  = %d x %d; should be 512 x 512.\n",
            width,height
            ENDFB;
          
          ok=false;      
        }
      } else {
        PRINTFB(FB_Color,FB_Errors)
          " ColorTableLoad-Error: unable to load '%s'.\n",fname
          ENDFB;
        ok=false;
      }
    } else {
      PRINTFB(FB_Color,FB_Actions)
        " Color: purged table; colors unchanged.\n"
        ENDFB;
      FreeP(I->ColorTable);
    }
  }
  if(!ok) {
    FreeP(table);
  } else {
    ExecutiveInvalidateRep(cKeywordAll,cRepAll,cRepInvColor);
    SceneChanged();
  }
  return(ok);
}

/*========================================================================*/
void ColorUpdateClamp(int index)
{
  int i;
  int once=false;
  CColor *I=&Color; 
  unsigned int *entry;
  float *color,*new_color;
  unsigned int r,g,b,rr,gr,br;
  unsigned int ra,ga,ba;
  unsigned int rc[2][2][2],gc[2][2][2],bc[2][2][2];
  float fr,fg,fb,frm1,fgm1,fbm1,rct,gct,bct;
  int x,y,z;

  i = index;
  if(index>=0) {
    once=true;
  }
  for(i=0;i<I->NColor;i++) {
    if(!once) index=i;
   
    if(index<I->NColor) {
      if(!I->ColorTable) {
        I->Color[index].ClampedFlag = false;
      } else {
        color = I->Color[index].Color;
        r = ((int)(255*color[0]+0.5F))&0xFF;
        g = ((int)(255*color[1]+0.5F))&0xFF;
        b = ((int)(255*color[2]+0.5F))&0xFF;

        rr = r&0x3;
        gr = g&0x3;
        br = b&0x3;

        r = (r>>2);
        g = (g>>2);
        b = (b>>2);

        /* now for a crude little trilinear */

        for(x=0;x<2;x++) {
          ra = r + x;
          if(ra>63) ra=63;
          for(y=0;y<2;y++) {
            ga = g + y;
            if(ga>63) ga=63;
            for(z=0;z<2;z++) {
              ba = b + z;
              if(ba>63) ba=63;
              
              entry = I->ColorTable + (ra<<12) + (ga<<6) + ba;
              
              if(I->BigEndian) {
                rc[x][y][z] = 0xFF&((*entry)>>24);
                gc[x][y][z] = 0xFF&((*entry)>>16);
                bc[x][y][z] = 0xFF&((*entry)>>8);
              } else {
                rc[x][y][z] = 0xFF&((*entry)    );
                gc[x][y][z] = 0xFF&((*entry)>> 8);
                bc[x][y][z] = 0xFF&((*entry)>>16);
              }
            }
          }
        }

        frm1 = rr/4.0F;
        fgm1 = gr/4.0F;
        fbm1 = br/4.0F;
        
        fr = 1.0F - frm1;
        fg = 1.0F - fgm1;
        fb = 1.0F - fbm1;

        rct = 0.4999F + 
          (fr   * fg   * fb   * rc[0][0][0]) + 
          (frm1 * fg   * fb   * rc[1][0][0]) + 
          (fr   * fgm1 * fb   * rc[0][1][0]) + 
          (fr   * fg   * fbm1 * rc[0][0][1]) + 
          (frm1 * fgm1 * fb   * rc[1][1][0]) + 
          (fr   * fgm1 * fbm1 * rc[0][1][1]) + 
          (frm1 * fg   * fbm1 * rc[1][0][1]) + 
          (frm1 * fgm1 * fbm1 * rc[1][1][1]);

        gct = 0.4999F + 
          (fr   * fg   * fb   * gc[0][0][0]) + 
          (frm1 * fg   * fb   * gc[1][0][0]) + 
          (fr   * fgm1 * fb   * gc[0][1][0]) + 
          (fr   * fg   * fbm1 * gc[0][0][1]) + 
          (frm1 * fgm1 * fb   * gc[1][1][0]) + 
          (fr   * fgm1 * fbm1 * gc[0][1][1]) + 
          (frm1 * fg   * fbm1 * gc[1][0][1]) + 
          (frm1 * fgm1 * fbm1 * gc[1][1][1]);

        bct = 0.4999F + 
          (fr   * fg   * fb   * bc[0][0][0]) + 
          (frm1 * fg   * fb   * bc[1][0][0]) + 
          (fr   * fgm1 * fb   * bc[0][1][0]) + 
          (fr   * fg   * fbm1 * bc[0][0][1]) + 
          (frm1 * fgm1 * fb   * bc[1][1][0]) + 
          (fr   * fgm1 * fbm1 * bc[0][1][1]) + 
          (frm1 * fg   * fbm1 * bc[1][0][1]) + 
          (frm1 * fgm1 * fbm1 * bc[1][1][1]);

        if(r>=63) rct+=rr;
        if(g>=63) gct+=gr;
        if(b>=63) bct+=br;

        if(rct<=2.0F) rct=0.0F; /* make sure black is black */
        if(gct<=2.0F) gct=0.0F;
        if(bct<=2.0F) bct=0.0F;

        new_color = I->Color[index].Clamped;
        new_color[0] = rct/255.0F;
        if(new_color[0]>1.0F) new_color[0]=1.0F;
        new_color[1] = gct/255.0F;
        if(new_color[1]>1.0F) new_color[1]=1.0F; 
        new_color[2] = bct/255.0F;
        if(new_color[2]>1.0F) new_color[2]=1.0F;

        PRINTFD(FB_Color)
          "%5.3f %5.3f %5.3f -> %5.3f %5.3f %5.3f\n",
               color[0],color[1],color[2],
               new_color[0],new_color[1],new_color[2]
          ENDFD;
      
        I->Color[index].ClampedFlag = true;
      }
    }
      
    if(once) break;
  }
}
/*========================================================================*/
void ColorInit(void)
{
  CColor *I=&Color;

  unsigned int test;
  unsigned char *testPtr;
  
  test = 0xFF000000;
  testPtr = (unsigned char*)&test;
  I->BigEndian = (*testPtr)&&1;

  I->Color=VLAMalloc(4300,sizeof(ColorRec),5,true);
  I->NColor=0;
  ColorReset();
  I->NExt=0;
  I->Ext=VLAMalloc(10,sizeof(ExtRec),5,true);
  I->ColorTable=NULL;
}

/*========================================================================*/
float *ColorGet(int index)
{
  CColor *I=&Color;
  float *ptr;
  if((index>=0)&&(index<I->NColor)) {
    if(I->Color[index].ClampedFlag&&(int)SettingGet(cSetting_clamp_colors))
      ptr = I->Color[index].Clamped;
    else
      ptr = I->Color[index].Color;
    return(ptr);
 } else
	 return(I->Color[0].Color);
}

