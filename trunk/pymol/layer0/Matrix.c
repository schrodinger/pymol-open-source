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
#include <GL/glut.h>

#include"Matrix.h"

/*========================================================================*/
void MatrixInvTransform3f(GLfloat *p, GLfloat *m, GLfloat *q)
{
  register GLfloat p1  = *p    , p2  = *(p+1), p3  = *(p+2);

  *(q++) = p1*(* m   ) + p2*(*(m+1)) + p3*(*(m+2));
  *(q++) = p1*(*(m+4)) + p2*(*(m+5)) + p3*(*(m+6));
  *(q  ) = p1*(*(m+8)) + p2*(*(m+9)) + p3*(*(m+10));

}
/*========================================================================*/
void MatrixTransform3f(GLfloat *p, GLfloat *m, GLfloat *q)
{
  register GLfloat p1  = *p    , p2  = *(p+1), p3  = *(p+2);

  *(q++) = p1*(* m   ) + p2*(*(m+4)) + p3*(*(m+8));
  *(q++) = p1*(*(m+1)) + p2*(*(m+5)) + p3*(*(m+9));
  *(q  ) = p1*(*(m+2)) + p2*(*(m+6)) + p3*(*(m+10));

}
