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
#include<stdlib.h>
#include<math.h>
#include<stdio.h>
#include<string.h>

#include"Base.h"
#include"OOMac.h"
#include"MemoryDebug.h"

#include"Color.h"

CColor Color;
/*========================================================================*/
void ColorDef(char *name,float *v)
{
  CColor *I=&Color;
  int color=-1;
  int a;
  int idx;
  for(a=0;a<I->NColor;a++)
	 {
		if(strcmp(name,I->Color[a].Name)==0) 
		  {
			 color=a;
			 break;
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
  
}
/*========================================================================*/
int ColorGetIndex(char *name)
{
  CColor *I=&Color;
  int color=1; /* default for unknown is white */
  int a;
  for(a=0;a<I->NColor;a++)
	 {
		if(strcmp(name,I->Color[a].Name)==0) 
		  {
			 color=a;
			 break;
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
void ColorFree(void)
{
  CColor *I=&Color;
  VLAFreeP(I->Color);
}

/*========================================================================*/
void ColorInit(void)
{
  CColor *I=&Color;

  I->Color=VLAlloc(ColorRec,100);
  I->NColor=0;

  strcpy(I->Color[I->NColor].Name,"white");
  I->Color[I->NColor].Color[0]=1.0;
  I->Color[I->NColor].Color[1]=1.0;
  I->Color[I->NColor].Color[2]=1.0;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"blue");
  I->Color[I->NColor].Color[0]=0.0;
  I->Color[I->NColor].Color[1]=0.0;
  I->Color[I->NColor].Color[2]=1.0;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"green");
  I->Color[I->NColor].Color[0]=0.0;
  I->Color[I->NColor].Color[1]=1.0;
  I->Color[I->NColor].Color[2]=0.0;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"red");
  I->Color[I->NColor].Color[0]=1.0;
  I->Color[I->NColor].Color[1]=0.0;
  I->Color[I->NColor].Color[2]=0.0;
  I->NColor++;


  strcpy(I->Color[I->NColor].Name,"cyan");
  I->Color[I->NColor].Color[0]=0.0;
  I->Color[I->NColor].Color[1]=1.0;
  I->Color[I->NColor].Color[2]=1.0;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"yellow");
  I->Color[I->NColor].Color[0]=1.0;
  I->Color[I->NColor].Color[1]=1.0;
  I->Color[I->NColor].Color[2]=0.0;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"violet");
  I->Color[I->NColor].Color[0]=1.0;
  I->Color[I->NColor].Color[1]=0.0;
  I->Color[I->NColor].Color[2]=1.0;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"magenta");
  I->Color[I->NColor].Color[0]=1.0;
  I->Color[I->NColor].Color[1]=0.0;
  I->Color[I->NColor].Color[2]=0.5;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"orange");
  I->Color[I->NColor].Color[0]=1.0;
  I->Color[I->NColor].Color[1]=0.5;
  I->Color[I->NColor].Color[2]=0.0;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"olive");
  I->Color[I->NColor].Color[0]=0.5;
  I->Color[I->NColor].Color[1]=0.5;
  I->Color[I->NColor].Color[2]=0.0;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"purple");
  I->Color[I->NColor].Color[0]=0.5;
  I->Color[I->NColor].Color[1]=0.0;
  I->Color[I->NColor].Color[2]=0.5;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"teal");
  I->Color[I->NColor].Color[0]=0.0;
  I->Color[I->NColor].Color[1]=0.5;
  I->Color[I->NColor].Color[2]=0.5;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"grey");
  I->Color[I->NColor].Color[0]=0.5;
  I->Color[I->NColor].Color[1]=0.5;
  I->Color[I->NColor].Color[2]=0.5;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"gray");
  I->Color[I->NColor].Color[0]=0.5;
  I->Color[I->NColor].Color[1]=0.5;
  I->Color[I->NColor].Color[2]=0.5;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"carbon");
  I->Color[I->NColor].Color[0]=0.2;
  I->Color[I->NColor].Color[1]=1.0;
  I->Color[I->NColor].Color[2]=0.2;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"nitrogen");
  I->Color[I->NColor].Color[0]=0.2;
  I->Color[I->NColor].Color[1]=0.2;
  I->Color[I->NColor].Color[2]=1.0;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"oxygen");
  I->Color[I->NColor].Color[0]=1.0;
  I->Color[I->NColor].Color[1]=0.2;
  I->Color[I->NColor].Color[2]=0.2;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"hydrogen");
  I->Color[I->NColor].Color[0]=0.9;
  I->Color[I->NColor].Color[1]=0.9;
  I->Color[I->NColor].Color[2]=0.9;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"sulfer");
  I->Color[I->NColor].Color[0]=1.0;
  I->Color[I->NColor].Color[1]=0.5;
  I->Color[I->NColor].Color[2]=0.0;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"grey90");
  I->Color[I->NColor].Color[0]=0.9;
  I->Color[I->NColor].Color[1]=0.9;
  I->Color[I->NColor].Color[2]=0.9;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"grey95");
  I->Color[I->NColor].Color[0]=0.95;
  I->Color[I->NColor].Color[1]=0.95;
  I->Color[I->NColor].Color[2]=0.95;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"grey80");
  I->Color[I->NColor].Color[0]=0.8;
  I->Color[I->NColor].Color[1]=0.8;
  I->Color[I->NColor].Color[2]=0.8;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"grey70");
  I->Color[I->NColor].Color[0]=0.7;
  I->Color[I->NColor].Color[1]=0.7;
  I->Color[I->NColor].Color[2]=0.7;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"grey60");
  I->Color[I->NColor].Color[0]=0.6;
  I->Color[I->NColor].Color[1]=0.6;
  I->Color[I->NColor].Color[2]=0.6;
  I->NColor++;


  strcpy(I->Color[I->NColor].Name,"tv_red");
  I->Color[I->NColor].Color[0]=1.0;
  I->Color[I->NColor].Color[1]=0.2;
  I->Color[I->NColor].Color[2]=0.2;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"tv_green");
  I->Color[I->NColor].Color[0]=0.2;
  I->Color[I->NColor].Color[1]=1.0;
  I->Color[I->NColor].Color[2]=0.2;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"tv_blue");
  I->Color[I->NColor].Color[0]=0.3;
  I->Color[I->NColor].Color[1]=0.3;
  I->Color[I->NColor].Color[2]=1.0;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"tv_yellow");
  I->Color[I->NColor].Color[0]=1.0;
  I->Color[I->NColor].Color[1]=1.0;
  I->Color[I->NColor].Color[2]=0.2;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"tv_yellow");
  I->Color[I->NColor].Color[0]=1.0;
  I->Color[I->NColor].Color[1]=1.0;
  I->Color[I->NColor].Color[2]=0.1;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"tv_orange");
  I->Color[I->NColor].Color[0]=1.0;
  I->Color[I->NColor].Color[1]=0.55;
  I->Color[I->NColor].Color[2]=0.15;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"br0");
  I->Color[I->NColor].Color[0]=0.1;
  I->Color[I->NColor].Color[1]=0.1;
  I->Color[I->NColor].Color[2]=1.0;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"br1");
  I->Color[I->NColor].Color[0]=0.2;
  I->Color[I->NColor].Color[1]=0.1;
  I->Color[I->NColor].Color[2]=0.9;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"br2");
  I->Color[I->NColor].Color[0]=0.3;
  I->Color[I->NColor].Color[1]=0.1;
  I->Color[I->NColor].Color[2]=0.8;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"br3");
  I->Color[I->NColor].Color[0]=0.4;
  I->Color[I->NColor].Color[1]=0.1;
  I->Color[I->NColor].Color[2]=0.7;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"br4");
  I->Color[I->NColor].Color[0]=0.5;
  I->Color[I->NColor].Color[1]=0.1;
  I->Color[I->NColor].Color[2]=0.6;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"br5");
  I->Color[I->NColor].Color[0]=0.6;
  I->Color[I->NColor].Color[1]=0.1;
  I->Color[I->NColor].Color[2]=0.5;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"br6");
  I->Color[I->NColor].Color[0]=0.7;
  I->Color[I->NColor].Color[1]=0.1;
  I->Color[I->NColor].Color[2]=0.4;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"br7");
  I->Color[I->NColor].Color[0]=0.8;
  I->Color[I->NColor].Color[1]=0.1;
  I->Color[I->NColor].Color[2]=0.3;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"br8");
  I->Color[I->NColor].Color[0]=0.9;
  I->Color[I->NColor].Color[1]=0.1;
  I->Color[I->NColor].Color[2]=0.2;
  I->NColor++;

  strcpy(I->Color[I->NColor].Name,"br9");
  I->Color[I->NColor].Color[0]=1.0;
  I->Color[I->NColor].Color[1]=0.1;
  I->Color[I->NColor].Color[2]=0.1;
  I->NColor++;

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

