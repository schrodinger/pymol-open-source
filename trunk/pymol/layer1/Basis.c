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
-*   Larry Coopet (various optimizations)
-* 
-*
Z* -------------------------------------------------------------------
*/

#ifndef _PYMOL_INLINE

#include"os_predef.h"
#include"os_std.h"

#include"MemoryDebug.h"
#include"Base.h"
#include"Basis.h"
#include"Err.h"
#include"Feedback.h"
#include"Util.h"
#include"MemoryCache.h"

static const float kR_SMALL4 = 0.0001F;

static float BasisFudge0;
static float BasisFudge1;

float ZLineClipPoint(float *base,float *point,float *alongNormalSq,float cutoff);

int ZLineToSphere(float *base,float *point,float *dir,float radius,float maxial,
                  float *sphere,float *asum,float *pre);

int ZLineToSphereCapped(float *base,float *point,float *dir,float radius,float maxial,
                  float *sphere,float *asum,int cap1,int cap2,float *pre);

#define FASTER_ER   1

#if FASTER_ER

/*========================================================================*/
#ifdef _PYMOL_INLINE
__inline__
#endif
int ZLineToSphere(float *base,float *point,float *dir,float radius,float maxial,
                  float *sphere,float *asum,float *pre)
{
   /* Strategy - find an imaginary sphere that lies at the correct point on
     the line segment, then treat as a sphere reflection */
   
   float perpAxis0, perpAxis1, perpAxis2, intra_p0,intra_p1,intra_p2;
   float perpDist,radial,axial,axial_sum,dangle,ab_dangle,axial_perp;
   float radialsq,tan_acos_dangle;
   float intra0,intra1,intra2,vradial0,vradial1,vradial2;
   const float _0   = 0.0f;
   const float dir0   = dir[0];
   const float dir1   = dir[1];
   const float dir2   = dir[2];
   float   dot;
   
   intra0 = point[0] - base[0];
   intra1 = point[1] - base[1];
   
   perpAxis0 = pre[0]; /* was cross_product(MinusZ,dir,perpAxis),normalize */
   perpAxis1 = pre[1];
   
   /* the perpAxis defines a perp-plane which includes the cyl-axis */
   
   /* get minimum distance between the lines */
   
   perpDist = intra0*perpAxis0 + intra1*perpAxis1;
   /* was dot_product3f(intra,perpAxis); */

   if((perpDist < _0 ? -perpDist : perpDist) > radius)
      return 0;
      
   /*if(fabs(perpDist) > radius)   return 0; */
   
   perpAxis2 = _0;
   intra2=point[2]-base[2];
   
   dangle   = -dir[2]; /* was dot(MinusZ,dir) */
   ab_dangle   = (float)fabs(dangle);
   if(ab_dangle > (1.0f - kR_SMALL4))
   {
      if(dangle > _0)
      {
         sphere[0] = point[0];
         sphere[1] = point[1];
         sphere[2] = point[2];
      }
      else
      {
         sphere[0]=dir0 * maxial + point[0];
         sphere[1]=dir1 * maxial + point[1];
         sphere[2]=dir[2] * maxial + point[2];
      }
      return(1);
   }

   if(ab_dangle > kR_SMALL4)
      tan_acos_dangle = (float)(sqrt1d(1.0-dangle*dangle) / dangle);
   else
      tan_acos_dangle = MAXFLOAT;
         
   /* now we need to define the triangle in the perp-plane  
     to figure out where the projected line intersection point is */
   
   /* first, compute radial distance in the perp-plane between the two starting points */
   
   dot = intra0*perpAxis0 + intra1*perpAxis1 + intra2*perpAxis2;
   intra_p0=intra0-perpAxis0*dot;
   intra_p1=intra1-perpAxis1*dot;
   intra_p2=intra2-perpAxis2*dot;  
   
   dot = intra_p0*dir0 + intra_p1*dir1 + intra_p2*dir2;
   vradial0=intra_p0-dir0*dot;
   vradial1=intra_p1-dir1*dot;
   vradial2=intra_p2-dir2*dot;  
   
   radialsq   = ((vradial0*vradial0) + (vradial1*vradial1) + (vradial2*vradial2));
   
   /* now figure out the axial distance along the cyl-line that will give us
     the point of closest approach */

   if(ab_dangle < kR_SMALL4)
      axial_perp   = _0;
   else
      axial_perp   = (float)(sqrt1f(radialsq)/tan_acos_dangle);
  
   axial = (float) sqrt1f( ( (intra_p0*intra_p0) + (intra_p1*intra_p1) + (intra_p2*intra_p2)) - radialsq );
   
   if( (intra_p0*dir0 + intra_p1*dir1 + intra_p2*dir2 ) >= _0 )
      axial = axial_perp - axial;
   else
      axial = axial_perp + axial;
   
   /* now we have to think about where the vector will actually strike the cylinder*/
   
   /* by definition, it must be perpdist away from the perp-plane becuase the perp-plane
     is parallel to the line, so we can compute the radial component to this point */
   
   radial = radius*radius - perpDist*perpDist;
   radial = (float)sqrt1f(radial);
   
   /* now the trick is figuring out how to adjust the axial distance to get the actual
     position along the cyl line which will give us a representative sphere */
   
   if(ab_dangle > kR_SMALL4)
      axial_sum = axial - radial/tan_acos_dangle;
   else
      axial_sum = axial;
   /*
    printf("radial2 %8.3f \n",radial);*/

   if(axial_sum < _0)
      axial_sum = _0;
   else if(axial_sum > maxial)
      axial_sum = maxial;
   
   sphere[0] = dir0 * axial_sum + point[0];
   sphere[1] = dir1 * axial_sum + point[1];
   sphere[2] = dir2 * axial_sum + point[2];
   
   *asum = axial_sum;
   /*  printf("==>%8.3f sphere %8.3f %8.3f %8.3f\n",base[1],sphere[1],axial_perp,axial);*/
   return(1);
}

#else

int ZLineToSphere(float *base,float *point,float *dir,float radius,float maxial,
                  float *sphere,float *asum,float *pre)
{
   /* Strategy - find an imaginary sphere that lies at the correct point on
     the line segment, then treat as a sphere reflection */
   
   float perpAxis[3],intra_p[3];
   float perpDist,radial,axial,axial_sum,dangle,ab_dangle,axial_perp;
   float radialsq,tan_acos_dangle;
   float intra[3],vradial[3];
   const float _0   = 0.0f;
   const float dir0   = dir[0];
   const float dir1   = dir[1];
   
   intra[0] = point[0]-base[0];
   intra[1] = point[1]-base[1];
   
   perpAxis[0] = pre[0]; /* was cross_product(MinusZ,dir,perpAxis),normalize */
   perpAxis[1] = pre[1];
   
   /* the perpAxis defines a perp-plane which includes the cyl-axis */
   
   /* get minimum distance between the lines */
   
   perpDist = intra[0]*perpAxis[0] + intra[1]*perpAxis[1];
   /* was dot_product3f(intra,perpAxis); */
   
   if(fabs(perpDist) > radius)
      return 0;
   
   perpAxis[2] = _0;
   intra[2]=point[2]-base[2];
   
   dangle   = -dir[2]; /* was dot(MinusZ,dir) */
   ab_dangle   = (float)fabs(dangle);
   if(ab_dangle > (1.0f - kR_SMALL4))
   {
      if(dangle > _0)
      {
         sphere[0] = point[0];
         sphere[1] = point[1];
         sphere[2] = point[2];
      }
      else
      {
         sphere[0]=dir0 * maxial + point[0];
         sphere[1]=dir1 * maxial + point[1];
         sphere[2]=dir[2] * maxial + point[2];
      }
      return(1);
   }

   /*tan_acos_dangle = tan(acos(dangle));*/
   if(ab_dangle > kR_SMALL4)
      tan_acos_dangle = (float)(sqrt1d(1.0-dangle*dangle)/dangle);
   else
      tan_acos_dangle = MAXFLOAT;
      
   /*
   printf("perpDist %8.3f\n",perpDist);
   printf("dir %8.3f %8.3f %8.3f\n",dir[0],dir[1],dir[2]);
   printf("base %8.3f %8.3f %8.3f\n",base[0],base[1],base[2]);
   printf("point %8.3f %8.3f %8.3f\n",point[0],point[1],point[2]);
   printf("intra %8.3f %8.3f %8.3f\n",intra[0],intra[1],intra[2]);
   printf("perpAxis %8.3f %8.3f %8.3f\n",perpAxis[0],perpAxis[1],perpAxis[2]);
   */
   
   /* now we need to define the triangle in the perp-plane  
     to figure out where the projected line intersection point is */
   
   /* first, compute radial distance in the perp-plane between the two starting points */
   
   remove_component3f(intra,perpAxis,intra_p);
   remove_component3f(intra_p,dir,vradial);
   radialsq = lengthsq3f(vradial);
   
   /* now figure out the axial distance along the cyl-line that will give us
     the point of closest approach */

   if(ab_dangle < kR_SMALL4)
      axial_perp   = _0;
   else
      axial_perp   = (float)(sqrt1f(radialsq)/tan_acos_dangle);
  
   axial = (float) sqrt1f( lengthsq3f(intra_p) - radialsq );
      
   /*
   printf("radial %8.3f\n",radial);
   printf("vradial %8.3f %8.3f %8.3f\n",vradial[0],vradial[1],vradial[2]);
   printf("radial %8.3f\n",radial);
   printf("dangle %8.3f \n",dangle);
   printf("axial_perp %8.3f \n",axial_perp);
   printf("axial1 %8.3f \n",axial);
   printf("%8.3f\n",dot_product3f(intra_p,dir));
   */
   
   if(dot_product3f(intra_p,dir) >= _0)
      axial = axial_perp - axial;
   else
      axial = axial_perp + axial;
   
   /*
   printf("axial2 %8.3f\n",axial);
   */
   
   /* now we have to think about where the vector will actually strike the cylinder*/
   
   /* by definition, it must be perpdist away from the perp-plane becuase the perp-plane
     is parallel to the line, so we can compute the radial component to this point */
   
   radial = radius*radius - perpDist*perpDist;
   radial = (float)sqrt1f(radial);
   
   /* now the trick is figuring out how to adjust the axial distance to get the actual
     position along the cyl line which will give us a representative sphere */
   
   if(ab_dangle > kR_SMALL4)
      axial_sum = axial - radial/tan_acos_dangle;
   else
      axial_sum = axial;
   /*
    printf("radial2 %8.3f \n",radial);*/

   if(axial_sum < _0)
      axial_sum = _0;
   else if(axial_sum > maxial)
      axial_sum = maxial;
   
   sphere[0] = dir0  * axial_sum + point[0];
   sphere[1] = dir1  * axial_sum + point[1];
   sphere[2] = dir[2]* axial_sum + point[2];
   
   *asum = axial_sum;
   /*  printf("==>%8.3f sphere %8.3f %8.3f %8.3f\n",base[1],sphere[1],axial_perp,axial);*/
   return(1);
}

#endif

#ifdef _PYMOL_INLINE
__inline__
#endif
static int ZLineFrontToInteriorSphere(float *front,
                                      float *point,
                                      float *dir,
                                      float radius,
                                      float radius2,
                                      float maxial)
{
  float intra_p[3];
  float axial;
  float intra[3],axis[3];
  float sphere[3];

  subtract3f(point,front,intra);
  remove_component3f(intra,dir,intra_p);
  add3f(front,intra_p,intra_p);
  subtract3f(point,intra_p,axis);
  axial = -dot_product3f(axis,dir);

  if(axial<0.0F) axial=0.0F;
  else if(axial>maxial) axial = maxial;

  sphere[0]=axial*dir[0]+point[0];
  sphere[1]=axial*dir[1]+point[1];
  sphere[2]=axial*dir[2]+point[2];

  return (diffsq3f(sphere,front)<=radius2);
}


/*========================================================================*/
#ifdef _PYMOL_INLINE
__inline__
#endif
int ZLineToSphereCapped(float *base,float *point,float *dir,float radius,float maxial,
                  float *sphere,float *asum,int cap1,int cap2,float *pre)
{
  /* Strategy - find an imaginary sphere that lies at the correct point on
     the line segment, then treat as a sphere reflection */

  float perpAxis[3],intra_p[3];
  float perpDist,radial,axial,axial_sum,dangle,ab_dangle,axial_perp;
  float radialsq,tan_acos_dangle;
  float len_proj;
  float intra[3],vradial[3];
  float diff[3],fpoint[3];
  float proj[3];

  perpAxis[0] = pre[0]; /* was cross_product(MinusZ,dir,perpAxis),normalize */
  perpAxis[1] = pre[1];

 /* the perpAxis defines a perp-plane which includes the cyl-axis */
  
  intra[0]=point[0]-base[0];
  intra[1]=point[1]-base[1];

  /* get minimum distance between the lines */

  perpDist = intra[0]*perpAxis[0] + intra[1]*perpAxis[1];
  /* was dot_product3f(intra,perpAxis); */

  if(fabs(perpDist)>radius) {
    return(0);
  }

  perpAxis[2] = 0.0;
  intra[2]=point[2]-base[2];

  dangle = -dir[2]; /* was dot(MinusZ,dir) */
  ab_dangle= (float)fabs(dangle);
  if(ab_dangle>(1-kR_SMALL4)) /* vector inline with light ray... */
    {
      vradial[0]=point[0]-base[0];
      vradial[1]=point[1]-base[1];
      vradial[2]=0.0F;
      radial = (float)length3f(vradial);
      if(radial>radius)
        return 0;
      if(dangle>0.0) {
        switch(cap1) {
        case cCylCapFlat:
          sphere[0]=base[0];
          sphere[1]=base[1];
          sphere[2]=point[2]-radius;
          break;
        case cCylCapRound:
          sphere[0]=point[0];
          sphere[1]=point[1];
          sphere[2]=point[2];
        }
        return(1);
      } else {
        switch(cap1) {
        case cCylCapFlat:
          sphere[0]=base[0];
          sphere[1]=base[1];
          sphere[2]=dir[2]*maxial+point[2]-radius;
          break;
        case cCylCapRound:
          sphere[0]=dir[0]*maxial+point[0];
          sphere[1]=dir[1]*maxial+point[1];
          sphere[2]=dir[2]*maxial+point[2];
        }
        return(1);
      }
    }
  
  /*tan_acos_dangle = tan(acos(dangle));*/
  tan_acos_dangle = (float)sqrt1f(1-dangle*dangle)/dangle;

  /*
  printf("perpDist %8.3f\n",perpDist);
  printf("dir %8.3f %8.3f %8.3f\n",dir[0],dir[1],dir[2]);
  printf("base %8.3f %8.3f %8.3f\n",base[0],base[1],base[2]);
  printf("point %8.3f %8.3f %8.3f\n",point[0],point[1],point[2]);
  printf("intra %8.3f %8.3f %8.3f\n",intra[0],intra[1],intra[2]);
  printf("perpAxis %8.3f %8.3f %8.3f\n",perpAxis[0],perpAxis[1],perpAxis[2]);
  */

  /* now we need to define the triangle in the perp-plane  
     to figure out where the projected line intersection point is */

  /* first, compute radial distance in the perp-plane between the two starting points */

  remove_component3f(intra,perpAxis,intra_p);
  remove_component3f(intra_p,dir,vradial);

  radialsq = lengthsq3f(vradial);

  /* now figure out the axial distance along the cyl-line that will give us
     the point of closest approach */

  if(ab_dangle<kR_SMALL4)
    axial_perp=0;
  else
    axial_perp = (float)sqrt1f(radialsq)/tan_acos_dangle;
  
  axial = (float)lengthsq3f(intra_p)-radialsq;
  axial = (float)sqrt1f(axial);

  /*
  printf("radial %8.3f\n",radial);
  printf("vradial %8.3f %8.3f %8.3f\n",vradial[0],vradial[1],vradial[2]);
  printf("radial %8.3f\n",radial);
  printf("dangle %8.3f \n",dangle);
  printf("axial_perp %8.3f \n",axial_perp);
  printf("axial1 %8.3f \n",axial);
  printf("%8.3f\n",dot_product3f(intra_p,dir));
  */

  if(dot_product3f(intra_p,dir)>=0.0)
    axial = axial_perp - axial;
  else
    axial = axial_perp + axial;

  /*
  printf("axial2 %8.3f\n",axial);
  */
  
  /* now we have to think about where the vector will actually strike the cylinder*/

  /* by definition, it must be perpdist away from the perp-plane becuase the perp-plane
     is parallel to the line, so we can compute the radial component to this point */

  radial = radius*radius-perpDist*perpDist;
  radial = (float)sqrt1f(radial);

  /* now the trick is figuring out how to adjust the axial distance to get the actual
     position along the cyl line which will give us a representative sphere */

  if(ab_dangle > kR_SMALL4)
    axial_sum = axial - radial/tan_acos_dangle;
  else
    axial_sum = axial;

  /*    printf("ab_dangle %8.3f \n",ab_dangle);
 
    printf("axial_sum %8.3f \n",axial_sum);
  */
  if(axial_sum<0) {
    switch(cap1) {
    case cCylCapFlat:
      subtract3f(point,base,diff);
      project3f(diff,dir,proj);
      len_proj = (float)length3f(proj);
      dangle = -proj[2]/len_proj;
      if(fabs(dangle)<kR_SMALL4) 
        return 0;
      sphere[0]=base[0];
      sphere[1]=base[1];
      sphere[2]=base[2]-len_proj/dangle;
      if(diff3f(sphere,point)>radius)
        return 0; 
      sphere[0]+=dir[0]*radius;
      sphere[1]+=dir[1]*radius;
      sphere[2]+=dir[2]*radius;
      *asum=0;
      break;
    case cCylCapRound:
      axial_sum=0;
      /*sphere[0]=dir[0]*axial_sum+point[0];
        sphere[1]=dir[1]*axial_sum+point[1];
        sphere[2]=dir[2]*axial_sum+point[2];*/
      sphere[0]=point[0];
      sphere[1]=point[1];
      sphere[2]=point[2];
      *asum = axial_sum;
      break;
    case cCylCapNone:
    default:
      return 0;
      break;
    }
  } else if(axial_sum>maxial) {
    switch(cap2) {

    case cCylCapFlat:
      scale3f(dir,maxial,fpoint);
      add3f(fpoint,point,fpoint);
      subtract3f(fpoint,base,diff);
      project3f(diff,dir,proj);
      len_proj = (float)length3f(proj);
      dangle = -proj[2]/len_proj;
      if(fabs(dangle)<kR_SMALL4) 
        return 0;
      sphere[0]=base[0];
      sphere[1]=base[1];
      sphere[2]=base[2]-len_proj/dangle;
      if(diff3f(sphere,fpoint)>radius)
        return 0; 
      sphere[0]-=dir[0]*radius;
      sphere[1]-=dir[1]*radius;
      sphere[2]-=dir[2]*radius;
      *asum=maxial;
      break;
    case cCylCapRound:
      axial_sum=maxial;
      sphere[0]=dir[0]*axial_sum+point[0];
      sphere[1]=dir[1]*axial_sum+point[1];
      sphere[2]=dir[2]*axial_sum+point[2];
      *asum = axial_sum;
      break;
    case cCylCapNone:
    default:
      return 0;
      break;
    }
  } else {
    sphere[0]=dir[0]*axial_sum+point[0];
    sphere[1]=dir[1]*axial_sum+point[1];
    sphere[2]=dir[2]*axial_sum+point[2];
  
    *asum = axial_sum;

    /*  printf("==>%8.3f sphere %8.3f %8.3f %8.3f\n",base[1],sphere[1],axial_perp,axial);*/
  }
  return(1);
  
}



#ifdef _PYMOL_INLINE
__inline__
#endif
static int ZLineFrontToInteriorSphereCapped(float *front,
                                            float *point,
                                            float *dir,
                                            float radius,
                                            float radius2,
                                            float maxial,
                                            int cap1,
                                            int cap2)
{
  float intra_p[3];
  float axial;
  float intra[3],axis[3];
  float sphere[3];

  subtract3f(point,front,intra);
  remove_component3f(intra,dir,intra_p);
  add3f(front,intra_p,intra_p);
  subtract3f(point,intra_p,axis);
  axial = -dot_product3f(axis,dir);

  if(axial<0.0F) return 0;
  else if(axial>maxial) return 0;

  sphere[0]=axial*dir[0]+point[0];
  sphere[1]=axial*dir[1]+point[1];
  sphere[2]=axial*dir[2]+point[2];

  return (diffsq3f(sphere,front)<radius2);
  
}


/*========================================================================*/
#ifdef _PYMOL_INLINE
__inline__
#endif
float ZLineClipPoint(float *base,float *point,float *alongNormalSq,float cutoff)
{
   float hyp0, hyp1, hyp2;
   float result   = MAXFLOAT;
   
   /* this routine determines whether or not a vector starting at "base"
     heading in the direction "normal" intersects a sphere located at "point".

     It returns how far along the vector the intersection with the plane is or
     MAXFLOAT if there isn't a relevant intersection

     NOTE: this routine has been optimized for normals along Z
     Optimizes-out vectors that are more than "cutoff" from "point" in x,y plane 
   */

   hyp0 = point[0] - base[0];
   if(fabs(hyp0) > cutoff)
      return result;
      
   hyp1 = point[1] - base[1];
   if(fabs(hyp1) > cutoff)
      return result;
      
   hyp2 = point[2] - base[2];

   if(hyp2 < 0.0)
   {
      (*alongNormalSq) = (hyp2 * hyp2);
      result   = (hyp0 * hyp0) + (hyp1 * hyp1);
   }
   return result;
}
/*========================================================================*/
void BasisSetupMatrix(CBasis *I)
{
  float oldZ[3] = { 0.0,0.0,1.0 };
  float newY[3];
  float dotgle,angle;

  cross_product3f(oldZ,I->LightNormal,newY);

  dotgle=dot_product3f(oldZ,I->LightNormal);
  
  if((1.0-fabs(dotgle))<kR_SMALL4)
    {
      dotgle=(float)(dotgle/fabs(dotgle));
      newY[0]=0.0;
      newY[1]=1.0;
      newY[2]=0.0;
    }
  
  normalize3f(newY);

  angle=(float)(-acos(dotgle));
  
  /* now all we gotta do is effect a rotation about the new Y axis to line up new Z with Z */
  
  rotation_to_matrix33f(newY,angle,I->Matrix);

  /*
  printf("%8.3f %8.3f %8.3f %8.3f\n",angle*180.0/cPI,newY[0],newY[1],newY[2]);
  
  matrix_transform33f3f(I->Matrix,newY,test);
  printf("   %8.3f %8.3f %8.3f\n",test[0],test[1],test[2]);

  printf("   %8.3f %8.3f %8.3f\n",I->LightNormal[0],I->LightNormal[1],I->LightNormal[2]);
  matrix_transform33f3f(I->Matrix,I->LightNormal,test);
  printf("   %8.3f %8.3f %8.3f\n",test[0],test[1],test[2]);

  printf(">%8.3f %8.3f %8.3f\n",I->Matrix[0][0],I->Matrix[0][1],I->Matrix[0][2]);
  printf(">%8.3f %8.3f %8.3f\n",I->Matrix[1][0],I->Matrix[1][1],I->Matrix[1][2]);
  printf(">%8.3f %8.3f %8.3f\n",I->Matrix[2][0],I->Matrix[2][1],I->Matrix[2][2]);
  */
}

/*========================================================================*/
void BasisGetTriangleFlatDotgle(CBasis *I,RayInfo *r,int i)
{
  float   *n0 = I->Normal + (3 * I->Vert2Normal[i]); 
  r->flat_dotgle = n0[2];
}

/*========================================================================*/
void BasisGetTriangleNormal(CBasis *I,RayInfo *r,int i,float *fc) 
{
   float *n0,w2,fc0,fc1,fc2;
   float vt1[3];
   CPrimitive   *lprim   = r->prim;
   
   r->impact[0] = r->base[0]; 
   r->impact[1] = r->base[1]; 
   r->impact[2] = r->base[2]-r->dist;

   n0 = I->Normal + (3 * I->Vert2Normal[i]) + 3; /* skip triangle normal */
   w2 = 1.0F - (r->tri1 + r->tri2);
   /*  printf("%8.3f %8.3f\n",r->tri[1],r->tri[2]);*/

   fc0 = (lprim->c2[0]*r->tri1)+(lprim->c3[0]*r->tri2)+(lprim->c1[0]*w2);
   fc1 = (lprim->c2[1]*r->tri1)+(lprim->c3[1]*r->tri2)+(lprim->c1[1]*w2);
   fc2 = (lprim->c2[2]*r->tri1)+(lprim->c3[2]*r->tri2)+(lprim->c1[2]*w2);

   scale3f(n0+3, r->tri1, r->surfnormal);
   scale3f(n0+6, r->tri2, vt1);
   add3f(vt1, r->surfnormal, r->surfnormal);

   scale3f(n0,w2,vt1);
   add3f(vt1,r->surfnormal,r->surfnormal);

   normalize3f(r->surfnormal);
   
   fc[0] = fc0;
   fc[1] = fc1;
   fc[2] = fc2;
}


void BasisSetFudge(float fudge)
{
  BasisFudge0 = 0.0F-fudge;
  BasisFudge1 = 1.0F+fudge;
}

#ifdef PROFILE_BASIS
int n_cells = 0;
int n_prims = 0;
int n_triangles = 0;
int n_spheres = 0;
int n_cylinders = 0;
int n_sausages = 0;
int n_skipped = 0;
#endif



#ifdef _PYMOL_INLINE
__inline__
#endif
int BasisHitNoShadow(BasisCallRec *BC)
{
   float   oppSq,dist,sph[3],vt[3],tri1,tri2; 
   int      a,b,c,h,*ip;
   int      excl_trans_flag;
   int      check_interior_flag;
   int      *elist, local_iflag = false;
   const float   _0   = 0.0F;
   
   CBasis *BI = BC->Basis;
   RayInfo *r = BC->rr;
   
   if( MapInsideXY(BI->Map, r->base, &a, &b, &c) )
   {
      register int      minIndex=-1;
      int     v2p;
      int     i,ii;
      int      *xxtmp;
             
      int except = BC->except;
      const int *vert2prim = BC->vert2prim;
      const float front = BC->front;
      const float back = BC->back;
      const float excl_trans = BC->excl_trans;
      MapCache *cache = &BC->cache;
      
      float r_tri1=_0, r_tri2=_0, r_dist=_0; /* zero inits to suppress compiler warnings */
      float r_sphere0=_0,r_sphere1=_0,r_sphere2=_0;
      CPrimitive *r_prim = NULL;
      
      check_interior_flag   = BC->check_interior;
      
      /* assumption: always heading in the negative Z direction with our vector... */
      vt[0]         = r->base[0];
      vt[1]         = r->base[1];
      vt[2]         = r->base[2] - front;
      
      if(except >= 0)
         except   = vert2prim[except];
         
      excl_trans_flag   = (excl_trans != _0);
      
      r_dist = MAXFLOAT;

      xxtmp   = BI->Map->EHead + (a * BI->Map->D1D2) + (b * BI->Map->Dim[2]) + c;

      MapCacheReset(cache);

      elist   = BI->Map->EList;
   
      while(c >= MapBorder) 
      {
         h   = *xxtmp;      
         if(h)
         {
            ip   = elist + h;
            i   = *(ip++);

            while(i >= 0) 
            {
               v2p = vert2prim[i];
               ii = *(ip++);
               if((v2p != except) && (!MapCached(cache,v2p))) 
               {
                  CPrimitive *prm = BC->prim + v2p;
                  MapCache(cache,v2p);
                  
                  switch(prm->type) 
                  {
                     case cPrimTriangle:
                  case cPrimCharacter:
                        if(!prm->cull)
                        {
                           float   *pre   = BI->Precomp + BI->Vert2Normal[i] * 3;
                           
                           if( pre[6] )
                           {
                              float   *vert0   = BI->Vertex + prm->vert * 3;
                              
                              float   tvec0   = vt[0] - vert0[0];
                              float   tvec1   = vt[1] - vert0[1];
                              
                              tri1      = (tvec0 * pre[4] - tvec1 * pre[3]) * pre[7];
                              tri2      = -(tvec0 * pre[1] - tvec1 * pre[0]) * pre[7];
                              
                              if( !( (tri1 < BasisFudge0) || (tri2 < BasisFudge0) || (tri1 > BasisFudge1) || ((tri1 + tri2) > BasisFudge1) ) )
                              {
                                 dist   = (r->base[2] - (tri1*pre[2]) - (tri2*pre[5]) - vert0[2]);

                                 if( (dist < r_dist) && (dist >= front) && (dist <= back) ) 
                                 {
                                    minIndex   = prm->vert;
                                    r_tri1      = tri1;
                                    r_tri2      = tri2;
                                    r_dist      = dist;
                                 }
                              }
                           }
                        }
                     break;
                     
                     case cPrimSphere:
                        oppSq = ZLineClipPoint( r->base, BI->Vertex + i*3, &dist, BI->Radius[i] );
                        if(oppSq <= BI->Radius2[i])
                        {
                           dist   = (float)(sqrt1f(dist) - sqrt1f((BI->Radius2[i]-oppSq)));

                           if(dist < r_dist) 
                           {
                              if((dist >= front) && (dist <= back)) 
                              {
                                 minIndex   = prm->vert;
                                 r_dist      = dist;
                              }
                              else if(check_interior_flag) 
                              {
                                 if(diffsq3f(vt,BI->Vertex+i*3) < BI->Radius2[i]) 
                                 {
                                    local_iflag   = true;
                                    r_prim      = prm;
                                    r_dist    = front;
                                    minIndex   = prm->vert;
                                 }
                              }
                           }
                        }
                     break;
                     
                     case cPrimCylinder:
                        if(ZLineToSphereCapped(r->base,BI->Vertex+i*3, 
                                               BI->Normal+BI->Vert2Normal[i]*3,
                                               BI->Radius[i], prm->l1,sph,&tri1,prm->cap1,prm->cap2,
                                               BI->Precomp + BI->Vert2Normal[i] * 3))
                        {
                           oppSq = ZLineClipPoint(r->base,sph,&dist,BI->Radius[i]);
                           if(oppSq<=BI->Radius2[i])
                           {
                              dist   = (float)(sqrt1f(dist)-sqrt1f((BI->Radius2[i]-oppSq)));

                              if(dist < r_dist) 
                              {
                                 if((dist >= front) && (dist <= back)) 
                                 {
                                    if(prm->l1 > kR_SMALL4)
                                       r_tri1   = tri1 / prm->l1;
                                       
                                    r_sphere0   = sph[0];
                                    r_sphere1   = sph[1];                              
                                    r_sphere2   = sph[2];
                                    minIndex   = prm->vert;
                                    r_dist      = dist;
                                 }
                                 else if(check_interior_flag)
                                 {
                                    if(ZLineFrontToInteriorSphereCapped(vt,
                                            BI->Vertex+i*3,
                                            BI->Normal+BI->Vert2Normal[i]*3,
                                            BI->Radius[i],
                                            BI->Radius2[i],
                                            prm->l1,
                                            prm->cap1,
                                            prm->cap2)) 
                                    {
                                       local_iflag   = true;
                                       r_prim      = prm;
                                       r_dist      = front;
                                       minIndex   = prm->vert;
                                    }
                                 }
                              }
                           }
                        }
                     break;
                     
                     case cPrimSausage:
                        if(ZLineToSphere(r->base,BI->Vertex+i*3,BI->Normal+BI->Vert2Normal[i]*3,BI->Radius[i],prm->l1,sph,&tri1,
                                         BI->Precomp + BI->Vert2Normal[i] * 3))
                        {
                           oppSq = ZLineClipPoint(r->base,sph,&dist,BI->Radius[i]);
                           if(oppSq <= BI->Radius2[i])
                           {
                              int   tmp_flag = false;
                              
                              dist   = (float)(sqrt1f(dist)-sqrt1f((BI->Radius2[i]-oppSq)));
                              if(dist<r_dist)
                              {
                                 if((dist >= front) && (dist <= back)) 
                                 {
                                    tmp_flag = true;
                                    if(excl_trans_flag) 
                                    {
                                       if( (prm->trans > _0) && (dist < excl_trans) )
                                          tmp_flag = false;
                                    }
                                    if(tmp_flag) 
                                    {
                                       if(prm->l1 > kR_SMALL4)
                                          r_tri1   = tri1 / prm->l1;
                                          
                                       r_sphere0   = sph[0];
                                       r_sphere1   = sph[1];                              
                                       r_sphere2   = sph[2];
                                       minIndex      = prm->vert;
                                       r_dist         = dist;
                                    }
                                 }
                                 else if(check_interior_flag) 
                                 {
                                    if(ZLineFrontToInteriorSphere(vt, BI->Vertex+i*3, BI->Normal+BI->Vert2Normal[i]*3,
                                                           BI->Radius[i], BI->Radius2[i], prm->l1)) 
                                    {
                                       local_iflag   = true;
                                       r_prim      = prm;
                                       r_dist      = front;
                                       minIndex   = prm->vert;
                                    }
                                 }
                              }
                           }
                        }
                     break;
                  }   /* end of switch */
               }   /* end of if */
               
               i = ii;
               
            } /* end of while */
         }

         /* and of course stop when we hit the edge of the map */
         
         if(local_iflag)
            break;
         
         /* we've processed all primitives associated with this voxel, 
         so if an intersection has been found which occurs in front of
         the next voxel, then we can stop */
         
         if( minIndex > -1 ) 
         {
            int   aa,bb,cc;
            
            vt[2]   = r->base[2] - r_dist;
            MapLocus(BI->Map,vt,&aa,&bb,&cc);
            if(cc > c) 
               break;
            else
               vt[2]   = r->base[2] - front;
         }
         
         c--;
         xxtmp--;
            
      } /* end of while */
      
      if( minIndex > -1 ) 
      {
         r_prim = BC->prim + vert2prim[minIndex];
         
         if(r_prim->type == cPrimSphere) 
         {
            const float   *vv   = BI->Vertex + minIndex * 3;
            r_sphere0   = vv[0];
            r_sphere1   = vv[1];
            r_sphere2   = vv[2];
         }
      }
      
      BC->interior_flag = local_iflag;   
      r->tri1 = r_tri1;
      r->tri2 = r_tri2;
      r->prim = r_prim;
      r->dist = r_dist;
      r->sphere[0] = r_sphere0;
      r->sphere[1] = r_sphere1;
      r->sphere[2] = r_sphere2;
      return(minIndex);
   } /* end of if */   
    BC->interior_flag = local_iflag;
   return(-1);
}



#ifdef _PYMOL_INLINE
__inline__
#endif
int BasisHitShadow(BasisCallRec *BC)
{
   float   oppSq,dist,sph[3],vt[3],tri1,tri2; 
   int      h,a,b,c,*ip;
   int      excl_trans_flag;
   int      check_interior_flag;
   int      *elist, local_iflag = false;
   const float   _0   = 0.0F;
   /* local copies (eliminate these extra copies later on) */

   CBasis   *BI   = BC->Basis;
   RayInfo   *r   = BC->rr;

   if( MapInsideXY(BI->Map,r->base, &a, &b, &c) )
   {
      register int      minIndex=-1;
      int     v2p;
       int     i,ii;
      int      *xxtmp;
             
      int except = BC->except;
      const int *vert2prim = BC->vert2prim;
      const int trans_shadows = BC->trans_shadows;
      const float excl_trans = BC->excl_trans;
      MapCache *cache = &BC->cache;
      
      float r_tri1=_0, r_tri2=_0, r_dist=_0;  /* zero inits to suppress compiler warnings */
      float r_sphere0=_0,r_sphere1=_0,r_sphere2=_0;
      CPrimitive *r_prim = NULL;
      
      check_interior_flag   = BC->check_interior;
      
      /* assumption: always heading in the negative Z direction with our vector... */
      vt[0]         = r->base[0];
      vt[1]         = r->base[1];
      
      if(except >= 0)
         except   = vert2prim[except];
         
      excl_trans_flag   = (excl_trans != _0);
      
      r_dist = MAXFLOAT;

      xxtmp   = BI->Map->EHead + (a * BI->Map->D1D2) + (b * BI->Map->Dim[2]) + c;

      MapCacheReset(cache);

      elist   = BI->Map->EList;
   
      while(c >= MapBorder) 
      {
         h   = *xxtmp;
         if( h )
         {
            ip   = elist + h;
            i   = *(ip++);

            while(i >= 0) 
            {
               v2p      = vert2prim[i];
               ii      = *(ip++);
               
               if( (v2p != except) && !MapCached(cache,v2p) ) 
               {
                  CPrimitive *prm = BC->prim + v2p;
                  
                  MapCache(cache,v2p);
                  
                  switch(prm->type) 
                  {
                  case cPrimCharacter: /* will need special handling for character shadows */
                    {
                      float   *pre   = BI->Precomp + BI->Vert2Normal[i] * 3;
                      
                      if( pre[6] )
                        {
                          float   *vert0   = BI->Vertex + prm->vert * 3;
                          
                          float   tvec0   = vt[0] - vert0[0];
                          float   tvec1   = vt[1] - vert0[1];
                          
                          tri1      = (tvec0 * pre[4] - tvec1 * pre[3]) * pre[7];
                          tri2      = -(tvec0 * pre[1] - tvec1 * pre[0]) * pre[7];
                          
                          if( !( (tri1 < BasisFudge0) || (tri2 < BasisFudge0) || (tri1 > BasisFudge1) || ((tri1 + tri2) > BasisFudge1) ) )
                            {
                              dist   = (r->base[2] - (tri1*pre[2]) - (tri2*pre[5]) - vert0[2]);
                              
                              {
                                float fc[3];
                                float trans;

                                r->tri1 = tri1;
                                r->tri2 = tri2;
                                r->dist = dist;
                                r->prim = prm;

                                BasisGetTriangleNormal(BI,r,minIndex,fc);
                                
                                trans = CharacterInterpolate(prm->char_id,fc);

                                if(trans == _0)  { /* opaque? return immed. */
                                  
                                  if((dist > -kR_SMALL4) && (dist < r_dist))
                                    {
                                      r->prim = prm;
                                      r->trans = _0;
                                      return(1);
                                    }
                                  
                                } /*else if(trans_shadows) 
                                  {
                                    if((dist > -kR_SMALL4) && (dist < r_dist)) 
                                      {
                                        minIndex   = prm->vert;
                                        r_tri1      = tri1;
                                        r_tri2      = tri2;
                                        r_dist      = dist;
                                      }
                                      }*/
                              }
                            }
                        }
                    }
                    break;
                    
                  case cPrimTriangle:
                     {
                        float   *pre   = BI->Precomp + BI->Vert2Normal[i] * 3;
                        
                        if( pre[6] )
                        {
                           float   *vert0   = BI->Vertex + prm->vert * 3;
                           
                           float   tvec0   = vt[0] - vert0[0];
                           float   tvec1   = vt[1] - vert0[1];
                           
                           tri1      = (tvec0 * pre[4] - tvec1 * pre[3]) * pre[7];
                           tri2      = -(tvec0 * pre[1] - tvec1 * pre[0]) * pre[7];
                           
                           if( !( (tri1 < BasisFudge0) || (tri2 < BasisFudge0) || (tri1 > BasisFudge1) || ((tri1 + tri2) > BasisFudge1) ) )
                           {
                              dist   = (r->base[2] - (tri1*pre[2]) - (tri2*pre[5]) - vert0[2]);

                              if(prm->trans == _0 ) /* opaque? return immed. */
                              {
                                 if((dist > -kR_SMALL4) && (dist < r_dist))
                                 {
                                    r->prim = prm;
                                    r->trans = _0;
                                    return(1);
                                 }
                              }
                              else if(trans_shadows) 
                              {
                                 if((dist > -kR_SMALL4) && (dist < r_dist)) 
                                 {
                                    minIndex   = prm->vert;
                                    r_tri1      = tri1;
                                    r_tri2      = tri2;
                                    r_dist      = dist;
                                    r->trans = prm->trans;
                                 }
                              }
                           }
                        }
                     }
                     break;
                     
                     case cPrimSphere:
                      
                        oppSq = ZLineClipPoint( r->base, BI->Vertex + i*3, &dist, BI->Radius[i] );
                        if(oppSq <= BI->Radius2[i])
                        {
                           dist   = (float)(sqrt1f(dist) - sqrt1f((BI->Radius2[i]-oppSq)));

                           if(prm->trans == _0) 
                           {
                              if((dist > -kR_SMALL4) && (dist < r_dist)) 
                              {
                                 r->prim = prm;
                                 r->trans = prm->trans;
                                 return(1);
                              }
                           }
                           else if(trans_shadows)
                           {
                              if((dist > -kR_SMALL4) && (dist < r_dist)) 
                              {
                                 minIndex   = prm->vert;
                                 r_dist      = dist;
                                 r->trans = prm->trans;
                              }
                           }
                        }
                     break;
                     
                     case cPrimCylinder:
                        if(ZLineToSphereCapped(r->base,BI->Vertex+i*3, 
                                               BI->Normal+BI->Vert2Normal[i]*3,
                                               BI->Radius[i], prm->l1,sph,&tri1,prm->cap1,prm->cap2,
                                               BI->Precomp + BI->Vert2Normal[i] * 3))
                        {
                           oppSq = ZLineClipPoint(r->base,sph,&dist,BI->Radius[i]);
                           if(oppSq <= BI->Radius2[i])
                           {
                              dist=(float)(sqrt1f(dist)-sqrt1f((BI->Radius2[i]-oppSq)));

                              if(prm->trans == _0) 
                              {
                                 if((dist > -kR_SMALL4) && (dist < r_dist)) 
                                 {
                                    r->prim = prm;
                                    r->trans = prm->trans;
                                    return(1);
                                 }
                              }
                              else if(trans_shadows) 
                              {
                                 if((dist > -kR_SMALL4) && (dist < r_dist)) 
                                 {
                                    if(prm->l1 > kR_SMALL4)
                                       r_tri1   = tri1 / prm->l1;
                                       
                                    r_sphere0   = sph[0];
                                    r_sphere1   = sph[1];                              
                                    r_sphere2   = sph[2];
                                    minIndex   = prm->vert;
                                    r->trans = prm->trans;
                                    r_dist      = dist;
                                 }
                              }
                           }
                        }
                     break;
                     
                     case cPrimSausage:
                        if(ZLineToSphere(r->base,BI->Vertex+i*3,BI->Normal+BI->Vert2Normal[i]*3,BI->Radius[i],prm->l1,sph,&tri1,
                                         BI->Precomp + BI->Vert2Normal[i] * 3))
                        {
                           oppSq = ZLineClipPoint(r->base,sph,&dist,BI->Radius[i]);
                           if(oppSq <= BI->Radius2[i])
                           {
                              dist   = (float)(sqrt1f(dist) - sqrt1f((BI->Radius2[i]-oppSq)));

                              if(prm->trans == _0) 
                              {
                                 if((dist > -kR_SMALL4) && (dist < r_dist)) 
                                 {
                                    r->prim = prm;
                                    r->trans = prm->trans;
                                    return(1);
                                 }
                              }
                              else if(trans_shadows) 
                              {
                                 if((dist > -kR_SMALL4) && (dist < r_dist)) 
                                 {
                                    if(prm->l1 > kR_SMALL4)
                                       r_tri1   = tri1 / prm->l1;
                                       
                                    r_sphere0   = sph[0];
                                    r_sphere1   = sph[1];                              
                                    r_sphere2   = sph[2];
                                    minIndex   = prm->vert;
                                    r->trans = prm->trans;
                                    r_dist      = dist;
                                 }
                              }
                           }
                        }
                     break;
                  }   /* end of switch */
               }   /* end of if */

               i = ii;
            } /* end of while */
         }

         /* and of course stop when we hit the edge of the map */
         
         if(local_iflag)
            break;
         
         /* we've processed all primitives associated with this voxel, 
         so if an intersection has been found which occurs in front of
         the next voxel, then we can stop */
         
         if( minIndex > -1 ) 
         {
            int   aa,bb,cc;
            
            vt[2]   = r->base[2] - r_dist;
            MapLocus(BI->Map,vt,&aa,&bb,&cc);
            if(cc > c) 
               break;
         }
         
         c--;
         xxtmp--;
            
      } /* end of while */
      
      if( minIndex > -1 ) 
      {
         r_prim = BC->prim + vert2prim[minIndex];
         
         if(r_prim->type == cPrimSphere) 
         {
            const float   *vv   = BI->Vertex + minIndex * 3;
            r_sphere0   = vv[0];
            r_sphere1   = vv[1];
            r_sphere2   = vv[2];
         }
      }
      
      BC->interior_flag = local_iflag;   
      r->tri1 = r_tri1;
      r->tri2 = r_tri2;
      r->prim = r_prim;
      r->dist = r_dist;
      r->sphere[0] = r_sphere0;
      r->sphere[1] = r_sphere1;
      r->sphere[2] = r_sphere2;
      return(minIndex);
   } /* end of if */   
    BC->interior_flag = local_iflag;
   return(-1);
}


#if 0
/*========================================================================*/
void BasisOptimizeMap(CBasis *I,float *vertex,int n,int *vert2prim);
void BasisOptimizeMap(CBasis *I,float *vertex,int n,int *vert2prim)
{
   int      *ip,*op;
   int      k;
   int     v2p;
   int     i,ii;
   MapCache cache;
   int a,b,c,aa,bb,cc;
   int h;
   int *xxtmp,*elist;
   MapCacheInit(&cache,I->Map);
   int cnt = 0,tcnt = 0;
   elist   = I->Map->EList;

     
#if 0
   for(k=0;k<n;k++) {
     if (MapExclLocus(I->Map,vertex + 3*k,&aa,&bb,&cc)) {
       for (a=aa-1;a<=aa+1;a++) {
         for (b=bb-1;b<=bb+1;b++) {
           xxtmp   = I->Map->EHead + (a * I->Map->D1D2) + (b * I->Map->Dim[2]);
           for(c=cc-1;c<=cc+1;c++) {
#else
{
  {
    for( a = I->Map->iMin[0]; a<=I->Map->iMax[0]; a++) {
      for(b = I->Map->iMin[1]; b<=I->Map->iMax[1]; b++) {
        xxtmp   = I->Map->EHead + (a * I->Map->D1D2) + (b * I->Map->Dim[2]);
        for( c = I->Map->iMin[2]; c <= I->Map->iMax[2]; c++)
          {
#endif

             h   = *(xxtmp + c);
             cnt = 0;
             tcnt = 0;
             if(h)
               {
                 op = ip = elist + h;
                 MapCacheReset(&cache);
                 
                 if(ip) {
                   i   = *(ip++);
                   while(i >= 0) 
                     {
                       tcnt++;
                       v2p = vert2prim[i];
                       ii = *(ip++);
                       /*    printf("%d\n",v2p);*/
                       if( !MapCached(&cache,v2p))
                         {
                           *(op++) = i; /* copy/store */
                           MapCache(&cache,v2p);
                         }   /* end of if */
                       else {
                         /* remove this vertex from the linked list -- we don't need it more than once */
                         cnt++;
                       }
                       i = ii;
                     } /* end of while */
                   *op = -1; /* terminate new list */
                 }
                 /*   printf("%6d %6d\n",tcnt,cnt);*/
               }
           
           }
         }
       }
     }
   }

}
#endif


/*========================================================================*/
void BasisMakeMap(CBasis *I,int *vert2prim,CPrimitive *prim,float *volume,
                  int group_id,int block_base)
{
   float *v,*vv,*d;
   float l;
   CPrimitive *prm;
   int a,b,c,i,n,h,q,x,y,z,j,k,e;
   int extra_vert = 0;
   float p[3],dd[3],*d1,*d2,vd[3],cx[3],cy[3];
   float *tempVertex;
   float xs,ys;
   int *tempRef,*ip,*sp;
   int remapMode = true; /* remap mode means that some objects will span more
                   * than one voxel, so we have to worry about populating
                   * those voxels and also about eliminating duplicates 
                   * when traversing the neighbor lists */
   float min[3],max[3],extent[6];
   float sep;
   float diagonal[3];
   float l1,l2;
   float bh,ch;
   int n_voxel;
   
   const float _0   = 0.0;
   const float _p5   = 0.5;
   
   PRINTFD(FB_Ray)
     " BasisMakeMap: I->NVertex %d [(%8.3f, %8.3f, %8.3f),...]\n",I->NVertex,
     I->Vertex[0],I->Vertex[1],I->Vertex[2]
   ENDFD;
   
   sep   = I->MinVoxel;
   if(sep == _0)
   {
      remapMode   = false;
      sep         = I->MaxRadius; /* also will imply no remapping of vertices */
   }
   /* we need to get a sense of the actual size in order to avoid sep being too small */
   
   v   = I->Vertex;
   
   min[0] = max[0] = v[0];
   min[1] = max[1] = v[1];
   min[2] = max[2] = v[2];
   
   v   += 3;
   
   for(a = 1; a < I->NVertex; a++)
   {
      if(min[0] > v[0])      min[0]   = v[0];
      if(max[0] < v[0])      max[0]   = v[0];

      if(min[1] > v[1])      min[1]   = v[1];
      if(max[1] < v[1])      max[1]   = v[1];

      if(min[2] > v[2])      min[2]   = v[2];
      if(max[2] < v[2])      max[2]   = v[2];

      v   += 3;
   }
   
   if(volume)
   {

     /*
      if(min[0] > volume[0])      min[0]   = volume[0];
      if(max[0] < volume[1])      max[0]   = volume[1];
      if(min[1] > volume[2])      min[1]   = volume[2];
      if(max[1] < volume[3])      max[1]   = volume[3];
      if(min[2] > (-volume[5]))   min[2]   = (-volume[5]);
      if(max[2] < (-volume[4]))   max[2]   = (-volume[4]);
     */

      if(min[0] < volume[0])      min[0]   = volume[0];
      if(max[0] > volume[1])      max[0]   = volume[1];
      if(min[1] < volume[2])      min[1]   = volume[2];
      if(max[1] > volume[3])      max[1]   = volume[3];
      if(min[2] > (-volume[5]))   min[2]   = (-volume[5]);
      if(max[2] < (-volume[4]))   max[2]   = (-volume[4]);

      PRINTFB(FB_Ray,FB_Debugging)
	" BasisMakeMap: (%8.3f,%8.3f),(%8.3f,%8.3f),(%8.3f,%8.3f)\n",
	volume[0],volume[1],volume[2],volume[3],volume[4],volume[5]
	ENDFB;
   }
   
   /* don't break up space unnecessarily if we only have a few vertices... */

   if(I->NVertex)
   {
      l1 = (float)fabs(max[0]-min[0]);
      l2 = (float)fabs(max[1]-min[1]);
      
      if(l1 < l2) l1 = l2;  
      
      l2 = (float)fabs(max[2]-min[2]);
      
      if(l1 < l2) l1 = l2;
      
      if(l1 < kR_SMALL4) l1 = 100.0;
      
      if(I->NVertex < (l1/sep))
         sep   = (l1/I->NVertex);
   }

   if(Feedback(FB_Ray,FB_Debugging)) {
     dump3f(min," BasisMakeMap: min");
     dump3f(max," BasisMakeMap: max");
     dump3f(I->Vertex," BasisMakeMap: I->Vertex");
     fflush(stdout);
   }

   sep = MapGetSeparation(sep,max,min,diagonal); /* this needs to be a minimum 
                                     * estimate of the actual value */
   
   /* here we have to carry out a complicated work-around in order to
   * efficiently encode our lines into the map in a way that doesn't
   * require expanding the map cutoff to the size of the largest object*/
   if(remapMode) 
   {
      for(a = 0; a < I->NVertex; a++)
      {
         prm   = prim + vert2prim[a];
         
         switch(prm->type)
         {
            case cPrimTriangle:
                  case cPrimCharacter:
               if(a == prm->vert) { /* only do this calculation for one of the three vertices */
                  l1   = (float)length3f(I->Precomp+I->Vert2Normal[a]*3);
                  l2   = (float)length3f(I->Precomp+I->Vert2Normal[a]*3+3);
                  b   = (int)ceil(l1/sep)+1;
                  c   = (int)ceil(l2/sep)+1;
                  extra_vert += 4*b*c;
               }
            break;
            
            case cPrimCylinder:
            case cPrimSausage:
               q = ((int)(2*(floor(prm->r1/sep)+1)))+1;
               q = q * q * ((int)ceil((prm->l1+2*prm->r1)/sep)+2);
               extra_vert+= q;
            break;
   
            case cPrimSphere:
               b = (int)(2*floor(prm->r1/sep)+1);
               extra_vert   += (b*b*b);
            break;
         }
      }   /* for */
   
      extra_vert   += I->NVertex;
      tempVertex   = CacheAlloc(float,extra_vert*3,group_id,cCache_basis_tempVertex);
      tempRef      = CacheAlloc(int,extra_vert,group_id,cCache_basis_tempRef); 
      
      ErrChkPtr(tempVertex); /* can happen if extra vert is unreasonable */
      ErrChkPtr(tempRef);
      
      /* lower indexes->flags, top is ref->lower index*/
      
      n   = I->NVertex;

      v   = tempVertex;
      vv   = I->Vertex;

#if 0
      for(a = 0; a < n; a++)
      {
         *(v++)   = *(vv++);
         *(v++)   = *(vv++);
         *(v++)   = *(vv++);
      }
#else
      memcpy( v, vv, n * sizeof(float) * 3 );
      vv   += n * 3;
      v   += n * 3;
#endif
      for(a = 0; a < I->NVertex; a++)
      {
         prm   = prim + vert2prim[a];
         
         switch(prm->type)
           {
         case cPrimTriangle:
         case cPrimCharacter:
               if(a == prm->vert)
               {
                  /* only do this calculation for one of the three vertices */
                  d1   = I->Precomp + I->Vert2Normal[a]*3;
                  d2   = I->Precomp + I->Vert2Normal[a]*3+3;
                  vv   = I->Vertex + a*3;
                  l1   = (float)length3f(d1);
                  l2   = (float)length3f(d2);
                  b   = (int)floor(l1/sep)+1;
                  c   = (int)floor(l2/sep)+1;
                  extra_vert += b*c;
                  bh   = (float)(b/2)+1;
                  ch   = (float)(c/2)+1;
                  
                  for(x = 0; x < bh; x++)
                  {
                     const float   xb   = (float)x / b;
                     
                     for(y = 0; y < ch; y++) 
                     {
                        const float   yc   = (float)y / c;
                        
                        *(v++)   = vv[0] + (d1[0] * xb) + (d2[0] * yc);
                        *(v++)   = vv[1] + (d1[1] * xb) + (d2[1] * yc);
                        *(v++)   = vv[2] + (d1[2] * xb) + (d2[2] * yc);
                        
                        tempRef[n++]   = a;
                     }
                  }
                      
                  for(x = 0; x < bh; x++)
                  {
                     const float   xb   = (float)x / b;
                     
                     for(y = 0; y < ch; y++) 
                     {
                        const float   yc   = (float)y/c;
                        
                        if( (xb + yc) < _p5 )
                        {
                           *(v++)   = vv[0] + d1[0] * (_p5 + xb) + (d2[0]*yc);
                           *(v++)   = vv[1] + d1[1] * (_p5 + xb) + (d2[1]*yc);
                           *(v++)   = vv[2] + d1[2] * (_p5 + xb) + (d2[2]*yc);
                           
                           tempRef[n++]   = a; 
                           
                           *(v++)   = vv[0] + (d1[0]*xb) + d2[0]*(_p5 + yc);
                           *(v++)   = vv[1] + (d1[1]*xb) + d2[1]*(_p5 + yc);
                           *(v++)   = vv[2] + (d1[2]*xb) + d2[2]*(_p5 + yc);
                           
                           tempRef[n++] = a;
                        }
                     }
                  }
               }   /* if */
               break;
          
            case cPrimCylinder:
            case cPrimSausage:
            
               d   = I->Normal+3*I->Vert2Normal[a];
               vv   = I->Vertex+a*3;
            
               get_system1f3f(d,cx,cy); /* creates an orthogonal system about d */
               
               p[0]   = vv[0]-d[0]*prm->r1;
               p[1]   = vv[1]-d[1]*prm->r1;
               p[2]   = vv[2]-d[2]*prm->r1;
               dd[0]   = d[0]*sep;
               dd[1]   = d[1]*sep;
               dd[2]   = d[2]*sep;
               l   = prm->l1+2*prm->r1;
            
               q   = (int)floor(prm->r1 / sep) + 1;
               
               while(1)
               {
                  vd[0]   = (p[0] += dd[0]);
                  vd[1]   = (p[1] += dd[1]);
                  vd[2]   = (p[2] += dd[2]);
               
                  for(x = -q; x <= q; x++)
                  {
                     for(y = -q; y <= q; y++)
                     {
                        xs      = x*sep;
                        ys      = y*sep;
                        *(v++)   = vd[0] + xs*cx[0] + ys*cy[0];
                        *(v++)   = vd[1] + xs*cx[1] + ys*cy[1];
                        *(v++)   = vd[2] + xs*cx[2] + ys*cy[2];

                        tempRef[n++]   = a;
                     }
                  }
                     
                  if(l <= _0)
                     break;
                  l   -= sep;
               }
            break;
          
            case cPrimSphere:
               q   = (int)floor(prm->r1 / sep);
               vv   = I->Vertex + a*3;
            
               for(x = -q; x <= q; x++)
               {
                  for(y = -q; y <= q; y++)
                  {
                     for( z = -q; z <= q; z++)
                     {
                        *(v++)   = vv[0] + x * sep;
                        *(v++)   = vv[1] + y * sep;
                        *(v++)   = vv[2] + z * sep;

                        tempRef[n++]   = a;
                     }
                  }
               }
            break;
         }   /* end of switch */
      }


      if(n > extra_vert)
      {
         printf("BasisMakeMap: %d>%d\n",n,extra_vert);
         ErrFatal("BasisMakeMap","used too many extra vertices (this indicates a bug)...\n");
      }

      if(volume)
      {
         v   = tempVertex;

         min[0]   = max[0]   = v[0];
         min[1]   = max[1]   = v[1];
         min[2]   = max[2]   = v[2];

         v += 3;

	 if(Feedback(FB_Ray,FB_Debugging)) {
	   dump3f(min," BasisMakeMap: remapped min");
	   dump3f(max," BasisMakeMap: remapped max");
	   fflush(stdout);
	 }

         
         for(a = 1; a < n; a++)
         {
	   if(min[0] > v[0])   min[0]   = v[0];
	   if(max[0] < v[0])   max[0]   = v[0];
	   
	   if(min[1] > v[1])   min[1]   = v[1];
	   if(max[1] < v[1])   max[1]   = v[1];
	   
	   if(min[2] > v[2])   min[2]   = v[2];
	   if(max[2] < v[2])   max[2]   = v[2];
	   v   += 3;
         }

	 if(Feedback(FB_Ray,FB_Debugging)) {
	   dump3f(min," BasisMakeMap: remapped min");
	   dump3f(max," BasisMakeMap: remapped max");
	   fflush(stdout);
	 }

         
         if(min[0] < volume[0])      min[0]   = volume[0];
         if(max[0] > volume[1])      max[0]   = volume[1];
         if(min[1] < volume[2])      min[1]   = volume[2];
         if(max[1] > volume[3])      max[1]   = volume[3];
         
         if(min[2] < (-volume[5]))   min[2]   = (-volume[5]);
         if(max[2] > (-volume[4]))   max[2]   = (-volume[4]);
            
         extent[0]   = min[0];
         extent[1]   = max[0];
         extent[2]   = min[1];
         extent[3]   = max[1];
         extent[4]   = min[2];
         extent[5]   = max[2];
	 PRINTFB(FB_Ray,FB_Blather)
	   " BasisMakeMap: Extent %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
	   extent[0],extent[1],extent[2],extent[3],extent[4],extent[5]
	   ENDFB;
         I->Map   = MapNewCached(-sep,tempVertex,n,extent,group_id,block_base);
      }
      else
      {
         I->Map   = MapNewCached(sep,tempVertex,n,NULL,group_id,block_base);
      }

      n_voxel = I->Map->Dim[0]*I->Map->Dim[1]*I->Map->Dim[2];
      
      if(n_voxel < (3*n))
         MapSetupExpressXY(I->Map,n);      
      else
         MapSetupExpressXYVert(I->Map,tempVertex,n);

      /* now do a filter-reassignment pass to remap fake vertices
      to the original line vertex while deleting duplicate entries */
      
#if 0
      ip   = tempRef;
      for(i = 0; i < I->NVertex; i++)
         *(ip++)   = 0; /* clear flags */
#else
      memset( tempRef, 0, sizeof(int) * I->NVertex );
#endif

      if(n_voxel < (3*n))
      {
         int   *start;
         
         for(a = I->Map->iMin[0]; a <= I->Map->iMax[0]; a++)
         {
            for(b = I->Map->iMin[1]; b <= I->Map->iMax[1]; b++)
            {
               for(c = I->Map->iMin[2]; c <= I->Map->iMax[2]; c++)
               {
                  start   = MapEStart(I->Map,a,b,c);
                  h      = *start;
                  if(h > 0)
                  {
                     ip   = I->Map->EList+h; 
                     sp   = ip;
                     i   = *(sp++);
                     while(i>=0)
                     {
                        if(i >= I->NVertex)
                           i   = tempRef[i];
                           
                        if(!tempRef[i]) { /*eliminate duplicates */
                           *(ip++)      = i;
                           tempRef[i]   = 1;
                        }
                        i   = *(sp++);
                     }
                     
                     *(ip)   = -1; /* terminate list */
                     /* now reset flags efficiently */
                     h   = *(start);
                     ip   = I->Map->EList+h;
                     i   = *(ip++);
                     
                     while(i>=0)
                     {
                        tempRef[i]   = 0;
                        i         = *(ip++);
                     }
                  }
               }   /* for c */
            }   /* for b */
         }   /* for a */
      }
      else
      {
         int      **site_p, **site = NULL;
         int      *value_p, *value = NULL;
         int      max_site;
         int      *start;
         int      n_site = 0;
         
         MapType   *mapPtr   = I->Map;
         
         max_site   = extra_vert;
         site      = CacheAlloc(int*,max_site,group_id,cCache_basis_site);
         value      = CacheAlloc(int,max_site,group_id,cCache_basis_value);
         
         site_p      = site;
         value_p      = value;
         v         = tempVertex;
         
         for(e = 0; e < n; e++)
         {
            MapLocus(mapPtr,v, &j, &k, &c);
            
            /*if( (c >= mapPtr->iMin[2]) && (c <= mapPtr->iMax[2]) ) */
            if( (c >= mapPtr->iMin[2]) && (c <= mapPtr->iMax[2]) && ((j-1) <= mapPtr->iMax[0]) && ((k-1) <= mapPtr->iMax[1]) )
            {
               int   *iPtr1   = mapPtr->EHead + ((j-1) * mapPtr->D1D2) + ((k-1) * mapPtr->Dim[2]) + c;
               
               for(a = j-1; a <= j+1; a++)
               {
                  if( (a >= mapPtr->iMin[0]) && (a <= mapPtr->iMax[0]) )
                  {
                     int   *iPtr2   = iPtr1;
                     
                     for(b = k-1; b <= k+1; b++)
                     {
                        if( (b >= mapPtr->iMin[1]) && (b <= mapPtr->iMax[1]) )
                        {
                           start   = iPtr2;
                           /*start   = MapEStart(mapPtr,a,b,c);*/
                           h      = *start;
                           if(h > 0)
                           {
                              int      ii;
                              
                              ip   = mapPtr->EList + h; 
                              sp   = ip;
                              i   = *(sp++);
                              
                              while(i >= 0)
                              {
                                 if(i >= I->NVertex)
                                    i   = tempRef[i];
                                    
                                 ii   = *(sp++);
                                 
                                 if(!tempRef[i]) { /*eliminate duplicates */
                                    *(ip++)      = i;
                                    tempRef[i]   = 1;
                                 }
                                 i = ii;
                              }
                              
                              *(ip)   = -1; /* terminate list */
                              /* now reset flags efficiently */
                              /*h      = *(start);*/
                              ip      = mapPtr->EList + h;
                              i      = *(ip++);
                              
                              while(i >= 0)
                              {
                                 tempRef[i]   = 0;
                                 i         = *(ip++);
                              }
                              
                              if((h > 0) && (n_site < max_site))
                              {
                                 *(value_p++)   =(*start);
                                 *(site_p++)      = start; /* remember which indexes we've negated */
                                 (*start)      = -1;
                                 n_site++;
                              } 
                           }   /* h > 0 */
                        }
                        
                        iPtr2   += mapPtr->Dim[2];
                        
                     }   /* for b */
                  }
                  
                  iPtr1   += mapPtr->D1D2;
                  
               }   /* for a */
            }
            
            v   += 3;   /* happens for EVERY e! */
            
         }   /* for e */
         
         value_p   = value;
         site_p   = site;
         
#if 0
         while(n_site--)
         {
            start   = *(site_p++);
            *start   = *(value_p++);
         }
#else
         while( n_site-- )
         {
            **(site_p++)   = *(value_p++);
         }
#endif         
         CacheFreeP(value,group_id,cCache_basis_value,false);
         CacheFreeP(site,group_id,cCache_basis_site,false);
      }
      
      CacheFreeP(tempVertex,group_id,cCache_basis_tempVertex,false);
      CacheFreeP(tempRef,group_id,cCache_basis_tempRef,false);
   }
   else
   {
      /* simple sphere mode */
      I->Map   = MapNewCached(-sep,I->Vertex,I->NVertex,NULL,group_id,block_base);
      MapSetupExpressXYVert(I->Map,I->Vertex,I->NVertex);
   }
}

/*========================================================================*/
void BasisInit(CBasis *I,int group_id)
{
  I->Vertex = VLACacheAlloc(float,1,group_id,cCache_basis_vertex);
  I->Radius = VLACacheAlloc(float,1,group_id,cCache_basis_radius);
  I->Radius2 = VLACacheAlloc(float,1,group_id,cCache_basis_radius2);
  I->Normal = VLACacheAlloc(float,1,group_id,cCache_basis_normal);
  I->Vert2Normal = VLACacheAlloc(int,1,group_id,cCache_basis_vert2normal);
  I->Precomp = VLACacheAlloc(float,1,group_id,cCache_basis_precomp);
  I->Map=NULL;
  I->NVertex=0;
  I->NNormal=0;
}
/*========================================================================*/
void BasisFinish(CBasis *I,int group_id)
{
  if(I->Map) 
    {
      MapFree(I->Map);
      I->Map=NULL;
    }  
  VLACacheFreeP(I->Radius2,group_id,cCache_basis_radius2,false);
  VLACacheFreeP(I->Radius,group_id,cCache_basis_radius,false);
  VLACacheFreeP(I->Vertex,group_id,cCache_basis_vertex,false);
  VLACacheFreeP(I->Vert2Normal,group_id,cCache_basis_vert2normal,false);
  VLACacheFreeP(I->Normal,group_id,cCache_basis_normal,false);
  VLACacheFreeP(I->Precomp,group_id,cCache_basis_precomp,false);
  I->Vertex=NULL;
}


/*========================================================================*/


#define EPSILON 0.000001

void BasisTrianglePrecompute(float *v0,float *v1,float *v2,float *pre)
{
   float det;
   
   subtract3f(v1,v0,pre);
   subtract3f(v2,v0,pre+3);
   
   det = pre[0]*pre[4] - pre[1]*pre[3];
   
   if(fabs(det) < EPSILON)
      *(pre+6) = 0.0F;
   else 
   {
      *(pre+6) = 1.0F;
      *(pre+7) = 1.0F/det;
   }
}


void BasisCylinderSausagePrecompute(float *dir,float *pre)
{
  float ln  = (float)(1.0f / sqrt1f(dir[1]*dir[1]+dir[0]*dir[0]));
  pre[0] = dir[1] * ln;
  pre[1] = -dir[0] * ln;
}


#if 0

static int intersect_triangle(float orig[3], float *pre,float vert0[3],
                              float *u, float *v, float *d)
{
   /* this routine now optimized to the point of total and complete opacity : ) */
   register float tvec0,tvec1,tv,tu;
   
#if !PRE_TEST   /* If this is on then we are testing for this in the caller! */
   if(!pre[6]) return 0;
#endif
   
   /* calculate distance from vert0 to ray origin */
   tvec0   = orig[0] - vert0[0];
   tvec1   = orig[1] - vert0[1];
   
   /* calculate U parameter and test bounds */
   tu   = (tvec0 * pre[4] - tvec1 * pre[3]) * pre[7];
      
   /* calculate V parameter and test bounds */
   tv   = -(tvec0 * pre[1] - tvec1 * pre[0]) * pre[7];
   
   if((tu < BasisFudge0) || (tv < BasisFudge0) || (tu > BasisFudge1) || ((tu + tv) > BasisFudge1) )
      return 0;
   
   /* calculate t, ray intersects triangle */
   *u = tu;
   *d = (orig[2] - (tu*pre[2]) - (tv*pre[5]) - vert0[2]);
   *v = tv;
   
   return 1;
}

#endif


#endif
