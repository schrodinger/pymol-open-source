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
#ifndef _H_ScrollBar
#define _H_ScrollBar

struct CScrollBar;

struct CScrollBar *ScrollBarNew(PyMOLGlobals *G,int horizontal);
void ScrollBarFree(struct CScrollBar *I);
void ScrollBarSetBox(struct CScrollBar *I,int top,int left,int bottom, int right);
void ScrollBarDoDraw(struct CScrollBar *I);
void ScrollBarSetLimits(struct CScrollBar *I,int list_size,int display_size);
void ScrollBarSetValue(struct CScrollBar *I,float value);
void ScrollBarDoRelease(struct CScrollBar *I,int button,int x,int y,int mod);
void ScrollBarDoClick(struct CScrollBar *I,int button,int x,int y,int mod);
void ScrollBarDoDrag(struct CScrollBar *I,int x,int y,int mod);
Block *ScrollBarGetBlock(struct CScrollBar *);
float ScrollBarGetValue(struct CScrollBar *I);
void ScrollBarMaxOut(struct CScrollBar *I);
void ScrollBarUpdate(struct CScrollBar *I);
int ScrollBarIsMaxed(struct CScrollBar *I);
void ScrollBarDrawHandle(struct CScrollBar *I,float alpha);

#endif

