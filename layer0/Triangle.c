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
#include<stdio.h>
#include<math.h>

#include"Base.h"
#include"Triangle.h"
#include"Map.h"

TriangleSurfaceRec *TrianglePointsToSurface(float *v,int n,float cutoff)
{
  MapType *map;

  map=MapNew(cutoff,v,n,NULL);
  MapSetupExpress(map);

  /*
  MapLocus(map,v,&h,&k,&l);
  flag=true;
  i=*(MapEStart(map,h,k,l));
  if(i) {
    j=map->EList[i++];
    while(j>=0) {
      if(j!=a) 
        {
          a2 = cs->IdxToAtm[j];
          if(within3f(cs->Coord+3*j,v,cs->Obj->AtomInfo[a2].vdw+probe_radius)) {
            flag=false;
            break;
          }
        }
      j=map->EList[i++];
    }
  }
  if(flag)
    {
  */

  MapFree(map);

}
