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

#include"Base.h"
#include"OOMac.h"
#include"MemoryDebug.h"
#include"Ortho.h"
#include"Word.h"

#include"Color.h"

CColor Color;
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
  VLAFreeP(I->Color);
}

/*========================================================================*/
void ColorInit(void)
{
  CColor *I=&Color;
  int a;
  int set1;
  float f;
  float spectrum[13][3] = { 
    { 1.0, 0.0, 1.0 },
    { 0.5, 0.0, 1.0 },
    { 0.0, 0.0, 1.0 },
    { 0.0, 0.5, 1.0 },
    { 0.0, 1.0, 1.0 },
    { 0.0, 1.0, 0.5 },
    { 0.0, 1.0, 0.0 },
    { 0.5, 1.0, 0.0 },
    { 1.0, 1.0, 0.0 },
    { 1.0, 0.5, 0.0 },
    { 1.0, 0.0, 0.0 },
    { 1.0, 0.0, 0.5 },
    { 1.0, 0.0, 0.5 }
  };

  I->Color=VLAlloc(ColorRec,2500);
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

  strcpy(I->Color[I->NColor].Name,"violet");
  I->Color[I->NColor].Color[0]=1.0F;
  I->Color[I->NColor].Color[1]=0.0F;
  I->Color[I->NColor].Color[2]=1.0F;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"salmon");
  I->Color[I->NColor].Color[0]=1.0F;
  I->Color[I->NColor].Color[1]=0.6F;
  I->Color[I->NColor].Color[2]=0.6F;
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

  strcpy(I->Color[I->NColor].Name,"magenta");
  I->Color[I->NColor].Color[0]=1.0F;
  I->Color[I->NColor].Color[1]=0.0F;
  I->Color[I->NColor].Color[2]=0.5F;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"orange");
  I->Color[I->NColor].Color[0]=1.0F;
  I->Color[I->NColor].Color[1]=0.5F;
  I->Color[I->NColor].Color[2]=0.0F;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"yellowgreen");
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

  strcpy(I->Color[I->NColor].Name,"gray");
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

  strcpy(I->Color[I->NColor].Name,"sulfer");
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

  strcpy(I->Color[I->NColor].Name,"grey100"); /* legacy = grey99 */
  I->Color[I->NColor].Color[0]=1.0F;
  I->Color[I->NColor].Color[1]=1.0F;
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

  /* full spectrum ("S..." colors) */

  #define A_DIV 90.9091

  for(a=0;a<1000;a=a+1) {
    set1=(int)(a/A_DIV);
    sprintf(I->Color[I->NColor].Name,"s%03d",a);
    f = 1.0-(a-(set1*A_DIV))/A_DIV;
    I->Color[I->NColor].Color[0]=f*spectrum[set1][0]+(1.0-f)*spectrum[set1+1][0];
    I->Color[I->NColor].Color[1]=f*spectrum[set1][1]+(1.0-f)*spectrum[set1+1][1];
    I->Color[I->NColor].Color[2]=f*spectrum[set1][2]+(1.0-f)*spectrum[set1+1][2];

    I->NColor++;
  }

}

/*========================================================================*/
float *ColorGet(int index)
{
  CColor *I=&Color;
  if((index>=0)&&(index<I->NColor))
	 return(I->Color[index].Color);
  else
	 return(I->Color[1].Color);
}

