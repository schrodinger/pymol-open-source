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
#include"os_gl.h"

#include"Base.h"
#include"OOMac.h"
#include"RepSphere.h"
#include"Color.h"
#include"Sphere.h"
#include"Map.h"
#include"Setting.h"
#include"main.h"
#include"Util.h"
#include"Feedback.h"

#ifdef NT
#undef NT
#endif


typedef struct RepSphere {
  Rep R;
  float *V; /* triangle vertices (if any) */
  float *VC; /* sphere centers, colors, alpha, and radii */
  float *VN; /* normals (if any) */
  SphereRec *SP,*SSP;
  int *NT;
  int N,NC,NP;
  int cullFlag,spheroidFlag;
  int *LastVisib;
  int *LastColor;
  float LastVertexScale;
  int VariableAlphaFlag;
  int shader_flag;
  
  GLint programs[2];
} RepSphere;

#include"ObjectMolecule.h"

#ifdef _PYMOL_OPENGL_SHADERS

#ifndef GL_FRAGMENT_PROGRAM_ARB 
#define GL_FRAGMENT_PROGRAM_ARB                         0x8804
#endif

/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */ 
#ifdef WIN32
static PFNGLGENPROGRAMSARBPROC glGenProgramsARB;
static PFNGLBINDPROGRAMARBPROC glBindProgramARB;
static PFNGLDELETEPROGRAMSARBPROC glDeleteProgramsARB;
static PFNGLPROGRAMSTRINGARBPROC glProgramStringARB;
static PFNGLPROGRAMENVPARAMETER4FARBPROC glProgramEnvParameter4fARB;
static PFNGLGETPROGRAMIVARBPROC glGetProgramivARB;
#endif
/* END PROPRIETARY CODE SEGMENT */

/* NOTE -- right now this shader program only runs in perspective mode */

typedef char ShaderCode[255];
ShaderCode vert_prog[] = {
  "!!ARBvp1.0\n",
  "\n",
  "# input contains the sphere radius in model coordinates\n",
  "PARAM sphereRadius = program.env[0];\n",
  "PARAM half = {0.5, 0.5, 0.0, 2.0 };\n",
  "PARAM zero = {0.0, 0.0, 0.0, 1.0 };\n",
  "\n",
  "ATTRIB vertexPosition  = vertex.position;\n",
  "ATTRIB vertexNormal    = vertex.normal;\n",
  "ATTRIB textureCoord    = vertex.texcoord;\n",
  "OUTPUT outputPosition  = result.position;\n",
  "\n",
  "TEMP   pos, rad, shf, txt, tip;\n",
  "\n",
  "# Transform the vertex by the modelview matrix to get into the frame of the camera\n",
  "\n",
  "DP4    pos.x, state.matrix.modelview.row[0], vertexPosition;\n",
  "DP4    pos.y, state.matrix.modelview.row[1], vertexPosition;\n",
  "DP4    pos.z, state.matrix.modelview.row[2], vertexPosition;\n",
  "DP4    pos.w, state.matrix.modelview.row[3], vertexPosition;\n",
  "\n",
  "# copy current texture coords\n",
  "MOV    txt.xyzw, textureCoord.xyzw;\n",
  "\n",
  "# scale the radius by a factor of two\n",
  "MUL    rad.xy, 2.0, sphereRadius.z;\n",
  "\n",
  "# shift the texture coordinates to the origin\n",
  "SUB    shf.xy, textureCoord, {0.5, 0.5, 0.0, 0.0};\n",
  "\n",
  "# multiply them to get the vertex offset\n",
  "\n",
  "MUL    shf.xy, rad, shf;\n",
  "\n",
  "# define the new vertex for corner of sphere\n",
  "\n",
  "ADD    pos.xy, pos, shf;\n",
  "\n",
  "# apply the projection matrix to get clip coordinates \n",
  "DP4    outputPosition.x, state.matrix.projection.row[0], pos;\n",
  "DP4    outputPosition.y, state.matrix.projection.row[1], pos;\n",
  "DP4    shf.z, state.matrix.projection.row[2], pos;\n",
  "DP4    shf.w, state.matrix.projection.row[3], pos;\n",
  "MOV    outputPosition.zw, shf;\n",
  "\n",
  "# compute camera position for front tip of the sphere\n",
  "ADD    pos.z, pos.z, sphereRadius;\n",
  "\n",
  "# compute Zc and Wc for front tip of the sphere\n",
  "DP4    tip.z, state.matrix.projection.row[2], pos;\n",
  "DP4    tip.w, state.matrix.projection.row[3], pos;\n",
  "\n",
  "# compute 1/Wc for sphere tip \n",
  "RCP    rad.z, tip.w;\n",
  "\n",
  "# put sphere center Zc into tip.w \n",
  "MOV    tip.w, shf.z;\n",
  "\n",
  "# compute 1/Wc for sphere center \n",
  "RCP    rad.w, shf.w;\n",
  "\n",
  "# compute Z/Wc for both sphere tip (->txt.z) and center (->txt.w) \n",
  "MUL    txt.zw, tip, rad;\n",
  "\n",
  "# move into range 0.0-1.0 to get the normalized depth coordinate (0.5*(Zc/Wc)+0.5) \n",
  "ADD    txt.zw, {0.0,0.0,1.0,1.0}, txt;\n",
  "MUL    txt.zw, {0.0,0.0,0.5,0.5}, txt;\n",
  "\n",
  "# Pass the color through\n",
  "MOV    result.color, vertex.color;\n",
  "\n",
  "# Pass texture through\n",
  "MOV    result.texcoord, txt;\n",
  "\n",
  "END\n",
  "\n",
  ""
};

ShaderCode frag_prog[] = {
  "!!ARBfp1.0\n",
  "\n",
  "PARAM fogInfo = program.env[0];\n",
  "PARAM fogColor = state.fog.color;\n",
  "ATTRIB fogCoord = fragment.fogcoord;\n",
  "\n",
  "TEMP pln, norm, depth, color, light, spec, fogFactor;\n",
  "\n",
  "# fully clip spheres that hit the camera\n",
  "KIL fragment.texcoord.z;\n",
  "\n",
  "# move texture coordinates to origin\n",
  "\n",
  "MOV norm.z, 0;\n",
  "SUB norm.xy, fragment.texcoord, {0.5,0.5,0.0,0.0};\n",
  "\n",
  "# compute x^2 + y^2, if > 0.25 then kill the pixel -- not in sphere\n",
  "\n",
  "# kill pixels that aren't in the center circle\n",
  "DP3 pln.z, norm, norm;\n",
  "SUB pln.z, 0.25, pln.z;\n",
  "KIL pln.z;\n",
  "\n",
  "# build a complete unit normal\n",
  "MUL pln.z, 4.0, pln.z;\n",
  "RSQ pln.z, pln.z;\n",
  "MUL norm.xy, 2.0, norm;\n",
  "RCP norm.z, pln.z;\n",
  "\n",
  "# interpolate the Zndc coordinate on the sphere \n",
  "LRP depth.z, norm.z, fragment.texcoord.z, fragment.texcoord.w;\n",
  "MOV result.depth.z, depth.z;\n",
  "\n",
  "# light0\n",
  "\n",
  "DP3 light, state.light[1].half, norm;\n",
  "MOV light.w, 60.0;\n",
  "LIT light, light;\n",
  "\n",
  "# ambient\n",
  "MOV color.xyzw, {0.06,0.06,0.06,1.0};\n",
  "ADD color.xyz, light.y, 0.1;\n",
  "MUL color.xyz, fragment.color, color;\n",
  "MUL spec.xyz, light.z, 0.5;\n",
  "ADD color.xyz, color,spec;\n",
  "\n",
  "# apply fog using linear interp over Zndc\n",
  "MAX fogFactor.x, depth.z, fogInfo.x;\n",
  "SUB fogFactor.x, fogFactor.x, fogInfo.x;\n",
  "MUL fogFactor.x, fogFactor.x, fogInfo.y;\n",
  "LRP color.xyz, fogFactor.x, fogColor, color;\n",
  "MOV result.color, color;\n",
  "\n",
  "END\n",
  "\n",
  ""
};

/*
  normal depth routine...does not work!  why?
  "#MAD_SAT fogFactor.x, fogInfo.x, fragment.texcoord.w, fogInfo.y;\n",
  "#LRP color.xyz, fogFactor.x, color, fogColor;\n",
*/

#endif

void RepSphereFree(RepSphere *I);
int RepSphereSameVis(RepSphere *I,CoordSet *cs);


void RepSphereFree(RepSphere *I)
{
#ifdef _PYMOL_OPENGL_SHADERS
  if(I->R.G->HaveGUI && I->R.G->ValidContext) {
    if(I->shader_flag) {

    }
  }
#endif

  FreeP(I->VC);
  FreeP(I->V);
  FreeP(I->VN);
  FreeP(I->NT);
  FreeP(I->LastColor);
  FreeP(I->LastVisib);
  RepPurge(&I->R);
  OOFreeP(I);
}

#ifdef _PYMOL_OPENGL_SHADERS
 
static GLboolean ProgramStringIsNative(PyMOLGlobals *G,
                                       GLenum target, GLenum format,   
                                       GLsizei len, const GLvoid *string)  
{  
  GLint errorPos, isNative;  
  glProgramStringARB(target, format, len, string);  
  glGetIntegerv(GL_PROGRAM_ERROR_POSITION_ARB, &errorPos);  
  glGetProgramivARB(target,
                    GL_PROGRAM_UNDER_NATIVE_LIMITS_ARB, &isNative);  
  if ((errorPos == -1) && (isNative == 1))  
    return GL_TRUE;  
  else if(errorPos >=0) {
    if(Feedback(G,FB_OpenGL, FB_Errors)) {
      printf("OpenGL-Error: ARB shader error at char %d\n---->%s\n",errorPos,((char*)string)+errorPos);
    }
  }
  return GL_FALSE;
}
 
static char *read_code_str(ShaderCode *ptr)
{
  ShaderCode *p = ptr;
  char *buffer,*q;
  int size = 0;
  while(*p[0]) {
    size += strlen(*p);
    p++;
  }
  buffer=Calloc(char,size+1);
  p = ptr;
  q = buffer;
  while(*p[0]) {
    strcat(q,*p);
    q+=strlen(q);
    p++;
  }
  return buffer;
}

#if 0
static char *read_file(char *name)
{
  char *buffer = NULL;
  FILE *f=fopen(name,"rb");
  size_t size;
  fseek(f,0,SEEK_END);
  size=ftell(f);
  fseek(f,0,SEEK_SET);
  buffer=(char*)mmalloc(size+1);
  fseek(f,0,SEEK_SET);
  fread(buffer,size,1,f);
  buffer[size]=0;
  fclose(f);
  return buffer;
}
#endif

static int load_shader_programs(PyMOLGlobals *G,GLuint *programs)
{
  int result = false;
  char *vp = read_code_str(vert_prog);
  char *fp = read_code_str(frag_prog);
      /*                  
                        char *vp = read_file("vert.txt");
                        char *fp = read_file("frag.txt");
    */

  /* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */ 
#ifdef WIN32
  if(!(glGenProgramsARB && glBindProgramARB && 
       glDeleteProgramsARB && glProgramStringARB && 
       glProgramEnvParameter4fARB)) {
    glGenProgramsARB = (PFNGLGENPROGRAMSARBPROC) wglGetProcAddress("glGenProgramsARB");
    glBindProgramARB = (PFNGLBINDPROGRAMARBPROC) wglGetProcAddress("glBindProgramARB");
    glDeleteProgramsARB = (PFNGLDELETEPROGRAMSARBPROC) wglGetProcAddress("glDeleteProgramsARB");
    glProgramStringARB = (PFNGLPROGRAMSTRINGARBPROC) wglGetProcAddress("glProgramStringARB");
    glProgramEnvParameter4fARB = (PFNGLPROGRAMENVPARAMETER4FARBPROC) wglGetProcAddress("glProgramEnvParameter4fARB");
    glGetProgramivARB = (PFNGLGETPROGRAMIVARBPROC) wglGetProcAddress("glGetProgramivARB");
  }
  
  if(glGenProgramsARB && glBindProgramARB && 
     glDeleteProgramsARB && glProgramStringARB && 
     glProgramEnvParameter4fARB) 
#endif
    {
      /* END PROPRIETARY CODE SEGMENT */
      
      if(vp&&fp) {            
        int ok=true;
        glGenProgramsARB(2,programs);
      
        /* load the vertex program */
        glBindProgramARB(GL_VERTEX_PROGRAM_ARB,programs[0]);
      
        ok = ok && (ProgramStringIsNative(G,GL_VERTEX_PROGRAM_ARB, 
                                          GL_PROGRAM_FORMAT_ASCII_ARB, 
                                          strlen(vp),vp));
      
        if(Feedback(G,FB_OpenGL,FB_Debugging))
          PyMOLCheckOpenGLErr("loading vertex program");
      
        /* load the fragment program */
        glBindProgramARB(GL_FRAGMENT_PROGRAM_ARB,programs[1]);
      
        ok = ok && (ProgramStringIsNative(G,GL_FRAGMENT_PROGRAM_ARB, 
                                          GL_PROGRAM_FORMAT_ASCII_ARB, 
                                          strlen(fp),fp));
      
        if(Feedback(G,FB_OpenGL,FB_Debugging))
          PyMOLCheckOpenGLErr("loading fragment program");
        if(ok) {
          result=true;
        } else {
          result=false;
          glDeleteProgramsARB(2,programs);
        }
      }
      FreeP(vp);
      FreeP(fp);
    }
  return result;
}

#endif

#ifdef _PYMOL_OPENGL_SHADERS
/* MULTI-INSTSANCE TODO:  isn't this a conflict? */
static GLuint programs_kludge[2] = {0,0};
#endif

void RepSphereRenderImmediate(CoordSet *cs, RenderInfo *info)
{
  PyMOLGlobals *G=cs->State.G;
  if(info->ray || info->pick || (!(G->HaveGUI && G->ValidContext)) )
    return;
  else {
    int repActive = false;
    ObjectMolecule *obj = cs->Obj;
    int sphere_mode = SettingGet_i(G,cs->Setting,obj->Obj.Setting,cSetting_sphere_mode);
    float sphere_scale = SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_sphere_scale);
    
    if(sphere_mode>0) { /* point-based modees */
      register float pixel_scale = 1.0F/info->vertex_scale;
      
      switch(sphere_mode) {
      case 2:
        glHint(GL_POINT_SMOOTH_HINT,GL_FASTEST);
        glDisable(GL_POINT_SMOOTH);
        glDisable(GL_ALPHA_TEST);
        pixel_scale *= 1.4F;
        glPointSize(1.0F);
        break;
      case 3: 
        glEnable(GL_POINT_SMOOTH);
        glAlphaFunc(GL_GREATER, 0.5F);
        glEnable(GL_ALPHA_TEST);
        glHint(GL_POINT_SMOOTH_HINT,GL_NICEST);
        glPointSize(1.0F);
        pixel_scale *= 2.0F;
        break;
      case 4:
        glEnable(GL_POINT_SMOOTH);
        glEnable(GL_ALPHA_TEST);
        glHint(GL_POINT_SMOOTH_HINT,GL_NICEST);
        glPointSize(1.0F);
        pixel_scale *= 2.0F;
        break;
      default:
        glHint(GL_POINT_SMOOTH_HINT,GL_FASTEST);
        glDisable(GL_POINT_SMOOTH);
        glDisable(GL_ALPHA_TEST);
        glPointSize(SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_sphere_point_size));
        break;
      }

#ifdef _PYMOL_OPENGL_SHADERS
      if(sphere_mode==5) {
        if(!(programs_kludge[0]||programs_kludge[1])) {
          load_shader_programs(G,programs_kludge);
        }
        if(programs_kludge[0]||programs_kludge[1]) {

          float fog_info[3];
          float _00[2] = { 0.0F, 0.0F};
          float _01[2] = { 0.0F, 1.0F};
          float _11[2] = { 1.0F, 1.0F};
          float _10[2] = { 1.0F, 0.0F};
          const float _1 = 1.0F;
          register float v0,v1,v2,nv0, nv1, nv2, nv3,v3;
          register float *m = info->pmv_matrix;                  
          register float cutoff = 1.2F;
          register float m_cutoff = -cutoff;
          register float z_front, z_back;
          /* compute -Ze = (Wc) of fog start */
          nv3 = (info->front + (info->back - info->front)*SettingGetGlobal_f(G,cSetting_fog_start));
          /* compute Zc of fog start using std. perspective transformation */
          nv2 = (nv3 * (info->back + info->front) - 2*(info->back * info->front))/(info->back - info->front);
          /* compute Zc/Wc to get normalized depth coordinate of fog start */
          nv0 = (nv2/nv3);
          fog_info[0] = (nv0*0.5)+0.5;

          fog_info[1] = 1.0F/(1.0-fog_info[0]); /* effective range of fog */
          
          z_front = info->stereo_front;
          z_back = info->back+((info->back+info->front)*0.25);
          
          /* load the vertex program */
          glBindProgramARB(GL_VERTEX_PROGRAM_ARB,programs_kludge[0]);
          
          /* load the fragment program */
          glBindProgramARB(GL_FRAGMENT_PROGRAM_ARB,programs_kludge[1]);
          
          /* load some safe initial values  */
          
          glProgramEnvParameter4fARB(GL_VERTEX_PROGRAM_ARB,
                                     0, 0.0F, 0.0F, 1.0, 0.0F);
          glProgramEnvParameter4fARB(GL_FRAGMENT_PROGRAM_ARB,
                                     0, 0.5F, 2.0F, 0.0F, 0.0F);
          
          glEnable(GL_VERTEX_PROGRAM_ARB);
          glEnable(GL_FRAGMENT_PROGRAM_ARB);
          
          glNormal3fv(info->view_normal);
          glBegin(GL_QUADS);

          {
            float last_radius = -1.0F, cur_radius;
            int a;
            int nIndex = cs->NIndex;
            AtomInfoType *atomInfo = obj->AtomInfo;
            int *i2a = cs->IdxToAtm;
            float *v = cs->Coord;
            
            for(a=0;a<nIndex;a++) {
              AtomInfoType *ai = atomInfo + *(i2a++);
              if(ai->visRep[ cRepSphere]) {
                repActive = true;
                v3 = ai->vdw * sphere_scale;
            
                v0 = v[0];
                v1 = v[1];
                v2 = v[2];
            
                if(last_radius!=(cur_radius=v3)) {
                  glEnd();
                  glProgramEnvParameter4fARB(GL_VERTEX_PROGRAM_ARB,
                                             0, 0.0F, 0.0F, v3, 0.0F);
                  glProgramEnvParameter4fARB(GL_FRAGMENT_PROGRAM_ARB,
                                             0, fog_info[0], fog_info[1], 0.0F, 0.0F);
                  glBegin(GL_QUADS);
                  last_radius = cur_radius;
                }
                
                /*  MatrixTransformC44f4f(info->pmv_matrix, v, nv);*/
                
                nv3 = m[ 3]*v0+m[ 7]*v1+m[11]*v2+m[15]; /* compute Wc */ 
                
                if(((nv3-cur_radius)>z_front) && (nv3<z_back)) { /* is it within the viewing volume? */
                  nv0 = m[ 0]*v0+m[ 4]*v1+m[ 8]*v2+m[12];
                  
                  nv3 = _1/nv3;
                  nv1 = m[ 1]*v0+m[ 5]*v1+m[ 9]*v2+m[13];
                  nv0 *= nv3;
                  nv1 *= nv3;
                  
                  
                  if((nv0<cutoff)&&(nv0>m_cutoff)&&
                     (nv1<cutoff)&&(nv1>m_cutoff)) {
                    
                    glColor3fv(ColorGet(G,ai->color));
                    
                    glTexCoord2fv(_00);
                    glVertex3fv(v);
                    
                    glTexCoord2fv(_10);
                    glVertex3fv(v);
                    
                    glTexCoord2fv(_11);
                    glVertex3fv(v);
                    
                    glTexCoord2fv(_01);
                    glVertex3fv(v);
                  }
                  
                }
              }
              v+=3;
            }
            glEnd();
          }
          glDisable(GL_FRAGMENT_PROGRAM_ARB);
          glDisable(GL_VERTEX_PROGRAM_ARB);
        }
      } else
#endif
      if(sphere_mode==4) {
        int repeat = true;
        const float _1 = 1.0F;
        const float _2 = 2.0F;
        float x_add = 0.0F, y_add = 0.0F, z_add = 0.0F;
        float z_factor = 0.0F, r_factor = 1.0F;
        float s_factor = 0.0F;
        int pass = 0;
        register float max_size = SettingGet_f(G,cs->Setting,obj->Obj.Setting,
                                               cSetting_sphere_point_max_size);
        register int clamp_size_flag = (max_size>=0.0F);

        
        while(repeat) {

          int a;
          int nIndex = cs->NIndex;
          AtomInfoType *atomInfo = obj->AtomInfo;
          int *i2a = cs->IdxToAtm;
          float *v = cs->Coord;
          float last_radius = -1.0F;
          float last_size = -1.0;
          float largest = 0.0F;

          float zz_factor = _1 - (float)pow(_1-z_factor,2);
          if(zz_factor<0.45F)
            zz_factor=0.45F;
          
          repeat = false;
          glBegin(GL_POINTS);
          
          for(a=0;a<nIndex;a++) {
            AtomInfoType *ai = atomInfo + *(i2a++);
            if(ai->visRep[ cRepSphere]) {
              float cur_radius = ai->vdw;
              repActive = true;

              if(last_radius!=cur_radius) {
                float clamp_radius = cur_radius;                        
                float size = cur_radius*pixel_scale;

                if(clamp_size_flag) 
                  if(size>max_size) {
                    size=max_size;
                    clamp_radius = size / pixel_scale;
                  }
                size *= r_factor;
                if( size != last_size ) {
                  glEnd();
                  if(size>largest)
                    largest = size;
                  if(size<_2) {
                    if(!pass) {
                      zz_factor=1.0F;
                      s_factor = 0.0F;
                    }
                  }
                  if(size<_1) {
                    size=_1;
                    glDisable(GL_POINT_SMOOTH);
                    glDisable(GL_ALPHA_TEST);
                  } else {
                    glEnable(GL_POINT_SMOOTH);
                    glEnable(GL_ALPHA_TEST);
                  }
                  glPointSize(size);
                  glBegin(GL_POINTS);
                }

                x_add = z_factor*clamp_radius*info->view_normal[0];
                y_add = z_factor*clamp_radius*info->view_normal[1];
                z_add = z_factor*clamp_radius*info->view_normal[2];
                last_radius = cur_radius;
                last_size = size;
              }
              
              {
                float *vc = ColorGet(G,ai->color);
                float r = zz_factor*vc[0] + s_factor;
                float g = zz_factor*vc[1] + s_factor;
                float b = zz_factor*vc[2] + s_factor;
                
                glColor3f(r > _1 ? _1 : r,
                          g > _1 ? _1 : g,
                          b > _1 ? _1 : b);
                
                glVertex3f(v[0]+x_add, v[1]+y_add, v[2]+z_add);
              }
            }
            v+=3;
          }
          
          glEnd();
          
          if(largest>2.0F) {
            float reduce = (largest-2.0F)/largest;
            r_factor *= reduce;
            z_factor = (float)sqrt1f(1.0F-(r_factor*r_factor));
            s_factor = (float)pow(z_factor,20.0F)*0.5F;
            repeat = true;
            pass++;
          }
        }
        glDisable(GL_POINT_SMOOTH);
      } else { 
        /* sphere_mode is 1, 2, or 3 */

        float max_radius = SettingGet_f(G,cs->Setting,obj->Obj.Setting,
                                        cSetting_sphere_point_max_size) * 3 * pixel_scale;
        int clamp_size_flag = (max_radius>=0.0F);
        
        int a;
        int nIndex = cs->NIndex;
        AtomInfoType *atomInfo = obj->AtomInfo;
        int *i2a = cs->IdxToAtm;
        int last_color = -1;
        float *v = cs->Coord;
        float last_radius = -1.0F;
        
        if(!info->line_lighting) glDisable(GL_LIGHTING);
        glBegin(GL_POINTS);
        
        for(a=0;a<nIndex;a++) {
          AtomInfoType *ai = atomInfo + *(i2a++);
          if(ai->visRep[ cRepSphere]) {
            int c = ai->color;
            repActive = true;
            if(c != last_color) {
              last_color = c;
              glColor3fv(ColorGet(G,c));
            }
            switch(sphere_mode) {
            case 1: 
              glVertex3fv(v);
              break;
            case 2:
            case 3:
              {
                float cur_radius = ai->vdw * pixel_scale;
                if(last_radius != cur_radius) {
                  glEnd();
                  if(clamp_size_flag) 
                    if(cur_radius > max_radius)
                      cur_radius = max_radius;
                  glPointSize(cur_radius);
                  glBegin(GL_POINTS);
                  last_radius = cur_radius;
                }
                glVertex3fv(v);
              }
              break;
            }
          }
          v+=3;
        }
        glEnd();
        glEnable(GL_LIGHTING);
        
        if(sphere_mode==3) {
          glDisable(GL_POINT_SMOOTH);
          glAlphaFunc(GL_GREATER, 0.05F);
        } else {
          glEnable(GL_ALPHA_TEST);
        }
      }
    } else { /* triangle-based spheres */

      SphereRec *sp = G->Sphere->Sphere[0];
      int ds = SettingGet_i(G,cs->Setting,obj->Obj.Setting,cSetting_sphere_quality);
      if(ds<0) {
        sp = NULL;
      } else {
        if(ds>4) ds=4;
        sp = G->Sphere->Sphere[ds];
      }
      
      {
        int a;
        int nIndex = cs->NIndex;
        AtomInfoType *atomInfo = obj->AtomInfo;
        int *i2a = cs->IdxToAtm;
        int last_color = -1;
        float *v = cs->Coord;
        int *sp_Sequence = sp->Sequence;
        int *sp_StripLen = sp->StripLen;
        int sp_NStrip = sp->NStrip;
        Vector3f *sp_dot = sp->dot;
        
        for(a=0;a<nIndex;a++) {
          AtomInfoType *ai = atomInfo + *(i2a++);
          if(ai->visRep[ cRepSphere]) {
            float vdw = ai->vdw * sphere_scale;
            int c = ai->color;
            float v0 = v[0];
            float v1 = v[1];
            float v2 = v[2];
            repActive = true;
            
            if(c != last_color) {
              last_color = c;
              glColor3fv(ColorGet(G,c));
            }
            
            {
              int *s = sp_StripLen;
              int *q = sp_Sequence;
              int b;
              for(b=0;b<sp_NStrip;b++) {
                int nc = *(s++);
                glBegin(GL_TRIANGLE_STRIP);	     
                for(c=0;c<nc;c++) {
                  float *sp_dot_q = &sp_dot[*(q++)][0];
                  glNormal3fv(sp_dot_q); /* normal */
                  glVertex3f(v0 + vdw*sp_dot_q[0],
                             v1 + vdw*sp_dot_q[1],
                             v2 + vdw*sp_dot_q[2]);
                }
                glEnd();
              }
            }
          }
          v+=3;
        }
      }
    }

    if(!repActive) /* didn't draw a single sphere, so we can skip this representation next time around */
      cs->Active[cRepSphere] = false;
  }
}

static void RepSphereRender(RepSphere *I,RenderInfo *info)
{
  CRay *ray = info->ray;
  Picking **pick = info->pick;
  PyMOLGlobals *G=I->R.G;
  float *v=I->V,*vc,*vn=I->VN;
  int c=I->N;
  int cc=0,*nt=NULL;
  int a;
  int flag;
  SphereRec *sp = I->SP;
  float restart;
  float alpha;

#ifdef _PYMOL_OPENGL_SHADERS
  /* TO DO -- garbage collect -- IMPORTANT! */
  {
    int sphere_mode = SettingGet_i(G,I->R.cs->Setting,
                                   I->R.obj->Setting,
                                   cSetting_sphere_mode);
    
    if((!ray) && (sphere_mode==5) && G->HaveGUI && G->ValidContext &&
       (!(I->programs[0]||I->programs[1]))) {
      I->shader_flag = load_shader_programs(G,(GLuint*)I->programs);
    }
  }
#endif

  alpha = SettingGet_f(G,I->R.cs->Setting,I->R.obj->Setting,cSetting_sphere_transparency);
  alpha=1.0F-alpha;
  if(fabs(alpha-1.0)<R_SMALL4)
    alpha=1.0F;
  if(ray) {
    ray->fTransparentf(ray,1.0-alpha);
    if(I->spheroidFlag) {
      if(sp) {
        while(c--)
          {
            vc = v;
            v+=3;
            for(a=0;a<sp->NStrip;a++) {
              cc=sp->StripLen[a];
              while((cc--)>2) {
                ray->fTriangle3fv(ray,v+3,v+9,v+15,v,v+6,v+12,vc,vc,vc);
                v+=6;
              }
              v+=12;
            }
          }
      }
    } else {
      int variable_alpha = I->VariableAlphaFlag;
      v=I->VC;
      c=I->NC;
      while(c--) {
        if(variable_alpha) {
          ray->fTransparentf(ray,1.0F-v[3]);
        }
        ray->fColor3fv(ray,v);
        v+=4;

#if 0
        /* temp code for testing ellipsoids */
        {
          float n1[3] = {1.0F, 0.0F, 0.0F};
          float n2[3] = {0.0F, 0.75F, 0.0F};
          float n3[3] = {0.0F, 0.0F, 0.35F};
          ray->fEllipsoid3fv(ray,v,*(v+3),n1,n2,n3);
        }
#else
        ray->fSphere3fv(ray,v,*(v+3));
#endif          
        v+=4;
      }
    }
    ray->fTransparentf(ray,0.0);
  } else if(G->HaveGUI && G->ValidContext) {
    if(pick) {
      int trans_pick_mode = SettingGet_i(G,I->R.cs->Setting,
                                         I->R.obj->Setting,
                                         cSetting_transparency_picking_mode);
        
      if(I->R.P&&((trans_pick_mode==1)||((trans_pick_mode==2)&&(alpha>0.9F)))) {
        int i,j;
        Pickable *p;
        i=(*pick)->src.index;
          
        p=I->R.P;
          
        if(I->spheroidFlag && sp) {
          while(c--) {
            int skip = (p[1].index<0);
            if(!skip) {
              i++;          
              if(!(*pick)[0].src.bond) {
                /* pass 1 - low order bits *            */
                glColor3ub((uchar)((i&0xF)<<4),(uchar)((i&0xF0)|0x8),(uchar)((i&0xF00)>>4)); 
                VLACheck((*pick),Picking,i);
                p++;
                (*pick)[i].src = *p; /* copy object and atom info */
                (*pick)[i].context = I->R.context;
              } else { 
                /* pass 2 - high order bits */           
                j=i>>12;            
                glColor3ub((uchar)((j&0xF)<<4),(uchar)((j&0xF0)|0x8),(uchar)((j&0xF00)>>4));             
              }			 
            } else {
              p++;
            }
            
            v+=4;
            for(a=0;a<sp->NStrip;a++) {
              cc=sp->StripLen[a];
              if(!skip) {
                glBegin(GL_TRIANGLE_STRIP);
                while((cc--)>0) {
                  glNormal3fv(v);
                  glVertex3fv(v+3);
                  v+=6;
                }
                glEnd();
              } else {
                while((cc--)>0) {
                  v+=6;
                }
              }
            }
          }
        } else {
          register float last_radius = -1.0F;
          register float cur_radius;
          register float pixel_scale = 1.0F/info->vertex_scale;
          register float max_size = SettingGet_f(G,I->R.cs->Setting,I->R.obj->Setting,
                                                 cSetting_sphere_point_max_size) * 3;
          register int clamp_size_flag = (max_size>=0.0F);
          register float size;
          int sphere_mode = SettingGet_i(G,I->R.cs->Setting,
                                         I->R.obj->Setting,
                                         cSetting_sphere_mode);
          
          if(!sp) {
            switch(sphere_mode) {
            case 5: 
            case 4:
            case 3:
            case 8:
              glEnable(GL_POINT_SMOOTH);
              glAlphaFunc(GL_GREATER, 0.5F);
              glEnable(GL_ALPHA_TEST);
              glHint(GL_POINT_SMOOTH_HINT,GL_NICEST);
              glPointSize(1.0F);
              pixel_scale *= 2.0F;
              break;
            case 2:
            case 7:
              glHint(GL_POINT_SMOOTH_HINT,GL_FASTEST);
              glDisable(GL_POINT_SMOOTH);
              glDisable(GL_ALPHA_TEST);
              pixel_scale *= 1.4F;
              break;
            default:
              glHint(GL_POINT_SMOOTH_HINT,GL_FASTEST);
              glDisable(GL_POINT_SMOOTH);
              glDisable(GL_ALPHA_TEST);
              glPointSize(SettingGet_f(G,I->R.cs->Setting,I->R.obj->Setting,cSetting_sphere_point_size));
              break;
            }
            glBegin(GL_POINTS);
          }

          v=I->VC;
          c=I->NC;
          while(c--) {
            int skip = (p[1].index<0);
            if(!skip) {
              i++;          
              if(!(*pick)[0].src.bond) {
                /* pass 1 - low order bits *            */
                glColor3ub((uchar)((i&0xF)<<4),(uchar)((i&0xF0)|0x8),(uchar)((i&0xF00)>>4)); 
                VLACheck((*pick),Picking,i);
                p++;
                (*pick)[i].src = *p; /* copy object and atom info */
                (*pick)[i].context = I->R.context;
              } else { 
                /* pass 2 - high order bits */           
                j=i>>12;            
                glColor3ub((uchar)((j&0xF)<<4),(uchar)((j&0xF0)|0x8),(uchar)((j&0xF00)>>4));             
              }			 
            } else {
              p++;
            }
              
            if(sp) {
              if(!skip) {
                int *s,*q,b;
                float *v0,vdw;
                
                v0 = v+4;
                vdw = v[7];
                q=sp->Sequence;
                s=sp->StripLen;
                for(b=0;b<sp->NStrip;b++) {
                  glBegin(GL_TRIANGLE_STRIP);
                  for(cc=0;cc<(*s);cc++)  {
                    glNormal3f(sp->dot[*q][0],
                               sp->dot[*q][1],
                               sp->dot[*q][2]);
                    glVertex3f(v0[0]+vdw*sp->dot[*q][0],
                               v0[1]+vdw*sp->dot[*q][1],
                               v0[2]+vdw*sp->dot[*q][2]);
                    q++;
                  }
                  glEnd();
                  s++;
                }
              }
              v+=8;
            } else {
              switch(sphere_mode) {
              case 2:
              case 3:
              case 4:
              case 5:
              case 7:
              case 8:
                if(!skip) {
                  if(last_radius!=(cur_radius=v[7])) {
                    size = cur_radius*pixel_scale;
                    glEnd();
                    if(clamp_size_flag) 
                      if(size>max_size)
                        size=max_size;
                    glPointSize(size);
                    glBegin(GL_POINTS);
                    last_radius = cur_radius;
                  }
                }
                v+=4;
                if(!skip) glVertex3fv(v);
                v+=4;
                break;
              default: /* simple, default point width points*/
                v+=4;
                if(!skip) glVertex3fv(v);
                v+=4;
                break;
              }
            }
          }
          if(!sp) {
            glEnd();
            switch(sphere_mode) {
            case 3:
            case 4:
            case 8:
              glDisable(GL_POINT_SMOOTH);
              glAlphaFunc(GL_GREATER, 0.05F);
              break;
            default:
              glEnable(GL_ALPHA_TEST);
              break;
            }
          }
        }
        (*pick)[0].src.index = i;
      }
    } else { /* not pick */
        
      if(!sp) {
        /* no sp -- we're rendering as points */
        int use_dlst;
        int sphere_mode = SettingGet_i(G,I->R.cs->Setting,
                                       I->R.obj->Setting,
                                       cSetting_sphere_mode);
        v=I->VC;
        c=I->NC;
          
        if(((sphere_mode>=2)&&(sphere_mode<=3))||
           ((sphere_mode>=7)&&(sphere_mode<=8))) { /* scaleable reps... */ 
          if(I->R.displayList) { 
            if(I->LastVertexScale != info->vertex_scale) {
              glDeleteLists(I->R.displayList,1);
              I->R.displayList = 0;
            }
          }
        }
        I->LastVertexScale = info->vertex_scale;
          
        use_dlst = (int)SettingGet(G,cSetting_use_display_lists);
        switch(sphere_mode) {
        case -1:
        case 0:
        case 4:
        case 5:
          use_dlst = 0;
          break;
        }
        if(use_dlst&&I->R.displayList) {
          glCallList(I->R.displayList);
        } else { /* display list */
            
          if(use_dlst) {
            if(!I->R.displayList) {
              I->R.displayList = glGenLists(1);
              if(I->R.displayList) {
                glNewList(I->R.displayList,GL_COMPILE_AND_EXECUTE);
              }
            }
          }
            
          if((sphere_mode>0)&&(!info->line_lighting)) glDisable(GL_LIGHTING);
            
          switch(sphere_mode) {
          case -1:
          case 0: /* memory-efficient sphere rendering */
            if(I->SSP) { 
              float last_vdw = -1.0F;
              int dlist = 0;
              int variable_alpha = I->VariableAlphaFlag;
              SphereRec *sp = I->SSP;
              while(c--) {
                Vector3f *sp_dot = sp->dot;
                int b,*q,*s;
                if(variable_alpha)
                  glColor4f(v[0],v[1],v[2],v[3]);
                else
                  glColor4f(v[0],v[1],v[2],alpha);
                v+=4;
                {
                  register float vdw = v[3];
                  glTranslatef(v[0],v[1],v[2]);
                  if((vdw!=last_vdw)||(!dlist)) {
                    last_vdw = vdw;
                    q=sp->Sequence;
                    s=sp->StripLen;
                    if(!dlist)
                      dlist = glGenLists(1);
                    if(dlist) {
                      glNewList(dlist,GL_COMPILE_AND_EXECUTE);
                    }
                    for(b=0;b<sp->NStrip;b++) {
                      int d;
                      glBegin(GL_TRIANGLE_STRIP);
                      for(d=0;d<(*s);d++) {
                        float *norm = sp_dot[*(q++)];
                        glNormal3fv(norm);
                        glVertex3f(vdw * norm[0],
                                   vdw * norm[1],
                                   vdw * norm[2]);
                      }
                      glEnd();
                      s++;
                    }
                    if(dlist) {
                      glEndList();
                    }
                  } else {
                    glCallList(dlist);
                  }
                  glTranslatef(-v[0],-v[1],-v[2]);
                }
                v+=4;
              }
              if(dlist)
                glDeleteLists(dlist,1);
            }
            break;
          case 2:
          case 3:
          case 4:
          case 7:
          case 8:
#ifdef _PYMOL_OPENGL_SHADERS
          case 5:
#endif
            {
              register float _1 = 1.0F;
              register float _2 = 2.0F;
              register float last_radius = -_1;
              register float cur_radius;
              register float pixel_scale = 1.0F/info->vertex_scale;
              register float max_size = SettingGet_f(G,I->R.cs->Setting,I->R.obj->Setting,
                                                     cSetting_sphere_point_max_size);
              register int clamp_size_flag = (max_size>=0.0F);
              register float size;
               
              if((sphere_mode==5) && (!I->shader_flag))
                sphere_mode=4;

              switch(sphere_mode) {
                
              case 2:
              case 3:
              case 7:
              case 8:
                if((sphere_mode==3)||(sphere_mode==8)) {
                  pixel_scale *= 2.0F;
                  glEnable(GL_POINT_SMOOTH);
                  glAlphaFunc(GL_GREATER, 0.5F);
                  glEnable(GL_ALPHA_TEST);
                  glHint(GL_POINT_SMOOTH_HINT,GL_NICEST);
                  glPointSize(1.0F);
                } else {
                  glHint(GL_POINT_SMOOTH_HINT,GL_FASTEST);
                  glDisable(GL_POINT_SMOOTH);
                  glDisable(GL_ALPHA_TEST);
                  pixel_scale *= 1.4F;
                }
                if((sphere_mode==7)||(sphere_mode==8))
                  glEnable(GL_LIGHTING);
                glBegin(GL_POINTS);
                while(c--) {
                  if(last_radius!=(cur_radius=v[7])) {
                    size = cur_radius*pixel_scale;
                    glEnd();
                    if(clamp_size_flag) 
                      if(size>max_size)
                        size=max_size;
                    glPointSize(size);
                    glBegin(GL_POINTS);
                    last_radius = cur_radius;
                  }
                  glColor3fv(v);
                  v+=4;
                  if(vn) {
                    glNormal3fv(vn);
                    vn+=3;
                  }
                  glVertex3fv(v);
                  v+=4;
                }
                glEnd();
                if(sphere_mode==3) {
                  glDisable(GL_POINT_SMOOTH);
                  glAlphaFunc(GL_GREATER, 0.05F);
                } else {
                  glEnable(GL_ALPHA_TEST);
                }
                break;
              case 4: /* draw multiple points of different radii and Z position */
#ifndef _PYMOL_OPENGL_SHADERS
              case 5:
#endif
                {
                  int repeat = true;
                  register float x_add= 0.0F, y_add= 0.0F, z_add = 0.0F;
                  register float z_factor=0.0F, r_factor = 1.0F;
                  register float largest;
                  register float r, g, b;
                  register float s_factor=0.0F;
                  register float zz_factor;
                  register float clamp_radius;
                  register float last_size;
                  int pass = 0;
                  glEnable(GL_POINT_SMOOTH);
                  glEnable(GL_ALPHA_TEST);
                  glHint(GL_POINT_SMOOTH_HINT,GL_NICEST);
                  glPointSize(1.0F);
                  
                  pixel_scale *= 2.0F;
                  while(repeat) {
                    v=I->VC;
                    c=I->NC;
                    largest = 0.0F;
                    zz_factor = _1 - (float)pow(_1-z_factor,2);
                    if(zz_factor<0.45F)
                      zz_factor=0.45F;
                    
                    last_radius = -1.0F;
                    last_size = -1.0F;
                    repeat = false;
                    glBegin(GL_POINTS);
                    while(c--) {
                      if(last_radius!=(cur_radius=v[7])) {
                        size = cur_radius*pixel_scale;
                        clamp_radius = cur_radius;                        
                        if(clamp_size_flag) 
                          if(size>max_size) {
                            size=max_size;
                            clamp_radius = size / pixel_scale;
                          }
                        size *= r_factor;
                        if( size != last_size ) {
                          glEnd();
                          if(size>largest)
                            largest = size;
                          if(size<_2) {
                            if(!pass) {
                              zz_factor=1.0F;
                              s_factor = 0.0F;
                            }
                          }
                          if(size<_1) {
                            size=_1;
                            glDisable(GL_POINT_SMOOTH);
                            glDisable(GL_ALPHA_TEST);
                          } else {
                            glEnable(GL_POINT_SMOOTH);
                            glEnable(GL_ALPHA_TEST);
                          }
                          glPointSize(size);
                          glBegin(GL_POINTS);
                        }
                        x_add = z_factor*clamp_radius*info->view_normal[0];
                        y_add = z_factor*clamp_radius*info->view_normal[1];
                        z_add = z_factor*clamp_radius*info->view_normal[2];
                        last_radius = cur_radius;
                        last_size = size;
                      }
                      r = zz_factor*v[0] + s_factor;
                      g = zz_factor*v[1] + s_factor;
                      b = zz_factor*v[2] + s_factor;
                      
                      glColor3f(r > _1 ? _1 : r,
                                g > _1 ? _1 : g,
                                b > _1 ? _1 : b);
                                
                      v+=4;
                      glVertex3f(v[0]+x_add, v[1]+y_add, v[2]+z_add);
                      v+=4;
                    }
                    glEnd();

                    if(largest>2.0F) {
                      float reduce = (largest-2.0F)/largest;
                      r_factor *= reduce;
                      z_factor = (float)sqrt1f(1.0F-(r_factor*r_factor));
                      s_factor = (float)pow(z_factor,20.0F)*0.5F;
                      repeat = true;
                      pass++;
                    }
                  }
                  glDisable(GL_POINT_SMOOTH);
                }
                break;
#ifdef _PYMOL_OPENGL_SHADERS
              case 5: /* use vertex/fragment program */
                if (I->shader_flag) {
                  float fog_info[3];
                  float _00[2] = { 0.0F, 0.0F};
                  float _01[2] = { 0.0F, 1.0F};
                  float _11[2] = { 1.0F, 1.0F};
                  float _10[2] = { 1.0F, 0.0F};
                  register float v0,v1,v2,nv0, nv1, nv2, nv3,v3;
                  register float *m = info->pmv_matrix;                  
                  register float cutoff = 1.2F;
                  register float m_cutoff = -cutoff;
			      register float z_front, z_back;
                  /* compute -Ze = (Wc) of fog start */
				  nv3 = (info->front + (info->back - info->front)*SettingGetGlobal_f(G,cSetting_fog_start));
				  /* compute Zc of fog start using std. perspective transformation */
				  nv2 = (nv3 * (info->back + info->front) - 2*(info->back * info->front))/(info->back - info->front);
				  /* compute Zc/Wc to get normalized depth coordinate of fog start */
				  nv0 = (nv2/nv3);
				  fog_info[0] = (nv0*0.5)+0.5;
				  /* printf("%8.3f %8.3f %8.3f %8.3f\n", nv3, nv2, nv0, fog_info[0]); */
				  fog_info[1] = 1.0F/(1.0-fog_info[0]); /* effective range of fog */
				  
                  z_front = info->stereo_front;
                  z_back = info->back+((info->back+info->front)*0.25);
			      
                  if(Feedback(G,FB_OpenGL,FB_Debugging))
                    PyMOLCheckOpenGLErr("before shader");

                  /* load the vertex program */
                  glBindProgramARB(GL_VERTEX_PROGRAM_ARB,I->programs[0]);
                   
                  /* load the fragment program */
                  glBindProgramARB(GL_FRAGMENT_PROGRAM_ARB,I->programs[1]);
                   
                  /* load some safe initial values  */

                  glProgramEnvParameter4fARB(GL_VERTEX_PROGRAM_ARB,
                                             0, 0.0F, 0.0F, 1.0, 0.0F);
                  glProgramEnvParameter4fARB(GL_FRAGMENT_PROGRAM_ARB,
                                             0, 0.5F, 2.0F, 0.0F, 0.0F);
                        
                  glEnable(GL_VERTEX_PROGRAM_ARB);
                  glEnable(GL_FRAGMENT_PROGRAM_ARB);
                   
                  {
                    last_radius = -1.0F;
                     
                    glNormal3fv(info->view_normal);
                    glBegin(GL_QUADS);
                    v+=4;

                    while(c--) {
                       
                      v3 = v[3];

                      v0 = v[0];
                      v1 = v[1];
                      v2 = v[2];

                      if(last_radius!=(cur_radius=v3)) {
                        glEnd();
                        glProgramEnvParameter4fARB(GL_VERTEX_PROGRAM_ARB,
                                                   0, 0.0F, 0.0F, v3, 0.0F);
                        glProgramEnvParameter4fARB(GL_FRAGMENT_PROGRAM_ARB,
                                                   0, fog_info[0], fog_info[1], 0.0F, 0.0F);
                        glBegin(GL_QUADS);
                        last_radius = cur_radius;
                      }
                       
                      /*  MatrixTransformC44f4f(info->pmv_matrix, v, nv);*/
					   
                      nv3 = m[ 3]*v0+m[ 7]*v1+m[11]*v2+m[15]; /* compute Wc */ 
					   
                      if(((nv3-cur_radius)>z_front) && (nv3<z_back)) { /* is it within the viewing volume? */
                        nv0 = m[ 0]*v0+m[ 4]*v1+m[ 8]*v2+m[12];
						   
                        nv3 = _1/nv3;
                        nv1 = m[ 1]*v0+m[ 5]*v1+m[ 9]*v2+m[13];
                        nv0 *= nv3;
                        nv1 *= nv3;
						   
						   
                        if((nv0<cutoff)&&(nv0>m_cutoff)&&
                           (nv1<cutoff)&&(nv1>m_cutoff)) {
							   
                          glColor3fv(v-4);                          
							   
                          glTexCoord2fv(_00);
                          glVertex3fv(v);
							   
                          glTexCoord2fv(_10);
                          glVertex3fv(v);
							   
                          glTexCoord2fv(_11);
                          glVertex3fv(v);
							   
                          glTexCoord2fv(_01);
                          glVertex3fv(v);
                        }
                      }
                      v+=8;
                    }
                    glEnd();
                  }
                   
                  glDisable(GL_FRAGMENT_PROGRAM_ARB);
                  glDisable(GL_VERTEX_PROGRAM_ARB);
                  if(Feedback(G,FB_OpenGL,FB_Debugging))
                    PyMOLCheckOpenGLErr("after shader");
                }
                break;
#endif
              }
            }
            break;
          default: /* simple, default point width points -- modes 1 or 6 */
            glHint(GL_POINT_SMOOTH_HINT,GL_FASTEST);
            glDisable(GL_POINT_SMOOTH);
            glDisable(GL_ALPHA_TEST);

            glPointSize(SettingGet_f(G,I->R.cs->Setting,I->R.obj->Setting,cSetting_sphere_point_size));
            glBegin(GL_POINTS);
            if(alpha==1.0) {
              if(vn) {
                glEnd();
                glEnable(GL_LIGHTING);
                glBegin(GL_POINTS);
                while(c--) {
                  glColor3fv(v);
                  v+=4;
                  glNormal3fv(vn);
                  vn+=3;
                  glVertex3fv(v);
                  v+=4;
                }
              } else {
                while(c--) {
                  glColor3fv(v);
                  v+=4;
                  glVertex3fv(v);
                  v+=4;
                }
              }
            } else {
              if(vn) {
                glEnd();
                glEnable(GL_LIGHTING);
                glBegin(GL_POINTS);
                while(c--) {
                  glColor4f(v[0],v[1],v[2],alpha);
                  v+=4;
                  glNormal3fv(vn);
                  vn+=3;
                  glVertex3fv(v);
                  v+=4;
                }
              } else {
                while(c--) {
                  glColor4f(v[0],v[1],v[2],alpha);
                  v+=4;
                  glVertex3fv(v);
                  v+=4;
                }
              }
            }
            glEnd();
            glEnable(GL_ALPHA_TEST);
            break;
              
          }
          glEnable(GL_LIGHTING);
            
          if(use_dlst&&I->R.displayList) {
            glEndList();
          }
        }
      } else { /* real spheres, drawn with triangles -- not points or impostors */
        int variable_alpha = I->VariableAlphaFlag;
        int use_dlst;
        
        use_dlst = (int)SettingGet(G,cSetting_use_display_lists);
        
        if(use_dlst&&I->R.displayList) {
          glCallList(I->R.displayList);
        } else { /* display list */
          
          if(use_dlst) {
            if(!I->R.displayList) {
              I->R.displayList = glGenLists(1);
              if(I->R.displayList) {
                glNewList(I->R.displayList,GL_COMPILE_AND_EXECUTE);
              }
            }
          }
          if(I->cullFlag) {
            
            if((alpha==1.0)&&(!variable_alpha)) {
              
              nt=I->NT; /* number of passes for each sphere */
              while(c--) { /* iterate through all atoms */
                glColor3fv(v);
                v+=4;
                cc=*(nt++);
                flag=0;
                glBegin(GL_TRIANGLE_STRIP);
                while(cc--) { /* execute loop this many times */
                  restart=*(v++);
                  if(restart) {
                    if(flag) {
                      glEnd();
                      glBegin(GL_TRIANGLE_STRIP);
                    }
                    if(restart==2.0) { /* swap triangle polarity */
                      glNormal3fv(v);
                      glVertex3fv(v+3);
                    }
                    glNormal3fv(v);
                    v+=3;
                    glVertex3fv(v);
                    v+=3;
                    glNormal3fv(v);
                    v+=3;
                    glVertex3fv(v);
                    v+=3;
                  }
                  glNormal3fv(v);
                  v+=3;
                  glVertex3fv(v);
                  v+=3;
                  flag=1;
                }
                glEnd();
              }
            } else {
              
              nt=I->NT; /* number of passes for each sphere */
              while(c--) { /* iterate through all atoms */
                
                glColor4f(v[0],v[1],v[2],v[3]);
                v+=4;
                cc=*(nt++);
                flag=0;
                glBegin(GL_TRIANGLE_STRIP);
                while(cc--) { /* execute loop this many times */
                  restart=*(v++);
                  if(restart) {
                    if(flag) {
                      glEnd();
                      glBegin(GL_TRIANGLE_STRIP);
                    }
                    if(restart==2.0) { /* swap triangle polarity */
                      glNormal3fv(v);
                      glVertex3fv(v+3);
                    }
                    glNormal3fv(v);
                    v+=3;
                    glVertex3fv(v);
                    v+=3;
                    glNormal3fv(v);
                    v+=3;
                    glVertex3fv(v);
                    v+=3;
                  }
                  glNormal3fv(v);
                  v+=3;
                  glVertex3fv(v);
                  v+=3;
                  flag=1;
                }
                glEnd();
              }
            }
          } else if(sp) {
            if((alpha==1.0)&&!variable_alpha) {
              while(c--) {
                glColor3fv(v);
                v+=4;
                for(a=0;a<sp->NStrip;a++) {
                  glBegin(GL_TRIANGLE_STRIP);
                  cc=sp->StripLen[a];
                  while(cc--) {
                    glNormal3fv(v);
                    v+=3;
                    glVertex3fv(v);
                    v+=3;
                  }
                  glEnd();
                }
              }
            } else {
              while(c--) {
                glColor4f(v[0],v[1],v[2],v[3]);
                v+=4;
                for(a=0;a<sp->NStrip;a++) {
                  glBegin(GL_TRIANGLE_STRIP);
                  cc=sp->StripLen[a];
                  while(cc--) {
                    glNormal3fv(v);
                    v+=3;
                    glVertex3fv(v);
                    v+=3;
                  }
                  glEnd();
                }
              }
            }
          } 
          if(use_dlst&&I->R.displayList) {
            glEndList();
          }
        }
      }
    }
  }
}

int RepSphereSameVis(RepSphere *I,CoordSet *cs)
{
  int same = true;
  int *lv,*lc,*cc;
  int a;
  AtomInfoType *ai;
  if(I->LastVisib && I->LastColor) {
    ai = cs->Obj->AtomInfo;
    lv = I->LastVisib;
    lc = I->LastColor;
    cc = cs->Color;
    
    for(a=0;a<cs->NIndex;a++)
      {
        if(*(lv++)!=(ai + cs->IdxToAtm[a])->visRep[cRepSphere] ) {
          same=false;
          break;
        }
        if(*(lc++)!=*(cc++)) {
          same=false;
          break;
        }
      }
  } else {
    same=false;
  }
  return(same);
}

static int RadiusOrder(float *list,int a,int b)
{
  return(list[a*8+7]<=list[b*8+7]);
}


Rep *RepSphereNew(CoordSet *cs,int state)
{
  PyMOLGlobals *G=cs->State.G;
  ObjectMolecule *obj;
  int ok=true;
  int a,b,c,a1,c1,a2,i,j,k,h,l;
  float *v,*v0,*vc,vdw,v1[3];
  float restart;
  int *q, *s,q0,q1,q2;
  int *lv,*lc,*cc;
  SphereRec *sp = G->Sphere->Sphere[0];
  int ds,*nt,flag;
  int *visFlag = NULL;
  MapType *map = NULL;
  int vFlag;
  AtomInfoType *ai2;
  int spheroidFlag = false;
  float spheroid_scale;
  float *sphLen,sphTmp,*sphNorm,*sphTmpN;
  float sphere_scale,sphere_add=0.0;
  int sphere_color;
  int *map_flag=NULL,*mf;
  int cartoon_side_chain_helper = 0;
  int ribbon_side_chain_helper = 0;
  AtomInfoType *ati1;
  int vis_flag;
  int sphere_mode;
  int *marked = NULL;
  float transp;
  int variable_alpha = false;
#ifdef _this_code_is_not_used
  float vv0[3],vv1[3],vv2[3];
  float tn[3],vt1[3],vt2[3],xtn[3],*tn0,*tn1,*tn2;
#endif
  int draw_mode = SettingGetGlobal_i(G,cSetting_draw_mode);
  int draw_quality = (((draw_mode == 1)||(draw_mode==-2)));

  OOCalloc(G,RepSphere);

  obj = cs->Obj;
  vFlag=false;
  if(obj->RepVisCache[cRepSphere])
    for(a=0;a<cs->NIndex;a++) {
      if(obj->AtomInfo[cs->IdxToAtm[a]].visRep[cRepSphere]) {
        vFlag=true;
        break;
      }
    }
  if(!vFlag) {
    OOFreeP(I);
    return(NULL); /* skip if no dots are visible */
  }
  marked = Calloc(int,obj->NAtom);
  RepInit(G,&I->R);
  I->shader_flag = false;
  I->programs[0] = 0;
  I->programs[1] = 0;

  ds = SettingGet_i(G,cs->Setting,obj->Obj.Setting,cSetting_sphere_quality);
  sphere_mode = SettingGet_i(G,cs->Setting,    obj->Obj.Setting,
                             cSetting_sphere_mode);
  if(sphere_mode>0)
    ds = -1;

  if(ds<0) {
    sp = NULL;
  } else {
    if(draw_quality && (ds<3))
      ds = 3;
    if(ds>4) ds=4;
    sp = G->Sphere->Sphere[ds];
  }
  
  sphere_color=SettingGet_color(G,cs->Setting,obj->Obj.Setting,cSetting_sphere_color);
  cartoon_side_chain_helper = SettingGet_b(G,cs->Setting, obj->Obj.Setting,
                                           cSetting_cartoon_side_chain_helper);
  ribbon_side_chain_helper = SettingGet_b(G,cs->Setting, obj->Obj.Setting,
                                          cSetting_ribbon_side_chain_helper);
  transp = SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_sphere_transparency);
  spheroid_scale=SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_spheroid_scale);
  if(sp&&spheroid_scale&&cs->Spheroid) 
    spheroidFlag=1;
  else
    spheroidFlag=0;

  sphere_scale=SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_sphere_scale);

  I->R.fRender=(void (*)(struct Rep *, RenderInfo *))RepSphereRender;
  I->R.fFree=(void (*)(struct Rep *))RepSphereFree;
  I->R.fSameVis=(int (*)(struct Rep*, struct CoordSet*))RepSphereSameVis;

  /* automatic --  OOcalloc 
     I->R.fRecolor=NULL;
     I->LastVisib=NULL;
     I->LastColor=NULL;
     I->NP = 0; 
  */

  I->LastVertexScale = -1.0F;
  I->R.obj=(CObject*)obj;
  I->R.cs = cs;
  I->R.context.object = (void*)obj;
  I->R.context.state = state;

  /* raytracing primitives */
  
  I->VC=(float*)mmalloc(sizeof(float)*cs->NIndex*8);
  ErrChkPtr(G,I->VC);
  I->NC=0;
  map_flag = Calloc(int,cs->NIndex);

  I->NT=NULL;
  nt = NULL;

  v=I->VC; 
  mf=map_flag;

  if(SettingGet_i(G,cs->Setting,obj->Obj.Setting,cSetting_sphere_solvent)) { /* are we generating a solvent surface? */
    sphere_add = SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_solvent_radius); /* if so, get solvent radius */
  }
  
  if(SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_pickable)) {
    I->R.P=Alloc(Pickable,cs->NIndex+1);
    ErrChkPtr(G,I->R.P);
  }

  I->spheroidFlag=spheroidFlag;
  for(a=0;a<cs->NIndex;a++) {
    a1 = cs->IdxToAtm[a];
    ati1 = obj->AtomInfo+a1;
    vis_flag = ati1->visRep[cRepSphere];
    
    if(vis_flag &&
       (!ati1->hetatm) &&
       (!(ati1->flags & cAtomFlag_solvent)) &&
       ((cartoon_side_chain_helper && ati1->visRep[cRepCartoon]) ||
        (ribbon_side_chain_helper && ati1->visRep[cRepRibbon]))) {
      
      register char *name1=ati1->name;
      register int prot1=ati1->protons;

      if(prot1 == cAN_N) { 
        if((!name1[1])&&(name1[0]=='N')) { /* N */
          register char *resn1 = ati1->resn;
          if(!((resn1[0]=='P')&&(resn1[1]=='R')&&(resn1[2]=='O')))
            vis_flag=false;
        }
      } else if(prot1 == cAN_O) { 
        if((!name1[1])&&(name1[0]=='O'))
          vis_flag=false;
      } else if(prot1 == cAN_C) {
        if((!name1[1])&&(name1[0]=='C'))
          vis_flag=false;
      }
    }

    marked[a1] = vis_flag; /* store temporary visibility information */
    
    if(vis_flag) {
      float at_sphere_scale;
      int at_sphere_color;
      float at_transp;

      AtomInfoGetSetting_f(G, ati1, cSetting_sphere_scale, sphere_scale, &at_sphere_scale);
      if(AtomInfoGetSetting_f(G, ati1, cSetting_sphere_transparency, transp, &at_transp))
        variable_alpha = true;
      AtomInfoGetSetting_color(G, ati1, cSetting_sphere_color, sphere_color, &at_sphere_color);

      if(I->R.P) {
        I->NP++;
        if(!ati1->masked) {
          I->R.P[I->NP].index = a1;
        } else {
          I->R.P[I->NP].index = -1;
        }
        I->R.P[I->NP].bond = -1;
      }
      
      *mf=true;
      I->NC++;
      if(at_sphere_color==-1)
        c1=*(cs->Color+a);
      else
        c1=at_sphere_color;
      v0 = cs->Coord+3*a;			 
      if(ColorCheckRamped(G,c1)) {
        ColorGetRamped(G,c1,v0,v,state);
        v+=3;
      } else {
        vc = ColorGet(G,c1); /* save new color */
        *(v++)=*(vc++);
        *(v++)=*(vc++);
        *(v++)=*(vc++);
      }
      *(v++) = 1.0F - at_transp;

      *(v++)=*(v0++); /* coordinate */
      *(v++)=*(v0++);
      *(v++)=*(v0++);
      *(v++)= obj->AtomInfo[a1].vdw*at_sphere_scale+sphere_add; /* radius */
    }
    mf++;
  }
  
  I->VariableAlphaFlag = variable_alpha;
  if(I->NC) 
    I->VC=ReallocForSure(I->VC,float,(v-I->VC));
  else
    I->VC=ReallocForSure(I->VC,float,1);
  if(I->R.P) {
    I->R.P = Realloc(I->R.P,Pickable,I->NP+1);
    I->R.P[0].index = I->NP;
  }


  if(variable_alpha) 
    I->cullFlag = false;
  else {
    I->cullFlag = SettingGet_i(G,cs->Setting,obj->Obj.Setting,cSetting_cull_spheres);
    
    if(I->cullFlag<0) {
      I->cullFlag = !(obj->NCSet>1);
    }
  }
  
  if(draw_quality)
    I->cullFlag = false;
  
  if(spheroidFlag || (!sp) ) I->cullFlag=false;
  
  if((I->cullFlag<2)&&
     (SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_sculpting))) 
    /* optimize build-time performance when sculpting */
    I->cullFlag=false;
  
  if((I->cullFlag<2)&&
     (SettingGet(G,cSetting_roving_spheres)!=0.0F))
    I->cullFlag=false;
  
  if(sp && (I->cullFlag<2) && (!spheroidFlag)) {
    /* don't cull unless forced */
    I->SSP = sp;
    sp = NULL;
  }

  if(!sp) { /* if sp==null, then we're drawing a point-based sphere rep */

    /* sort the vertices by radius */
    if(I->NC && I->VC && (!spheroidFlag) && (sphere_mode!=1)) {
      int *ix = Alloc(int,I->NC);
      float *vc_tmp = Alloc(float, I->NC*8);
      Pickable *pk_tmp = Alloc(Pickable,I->NP+1);
      int a;
      if(vc_tmp&&pk_tmp&&ix) {
        UtilCopyMem(vc_tmp, I->VC, sizeof(float)*8*I->NC);
        UtilCopyMem(pk_tmp, I->R.P, sizeof(Pickable)*(I->NP+1));
        
        UtilSortIndex(I->NC,I->VC,ix,(UtilOrderFn*)RadiusOrder);
        
        UtilCopyMem(I->R.P, pk_tmp, sizeof(Pickable));
        for(a=0;a<I->NC;a++) {
          UtilCopyMem(I->VC + (a*8), vc_tmp+(8*ix[a]), sizeof(float)*8);
          UtilCopyMem(I->R.P + (a+1), pk_tmp+ix[a]+1, sizeof(Pickable));
        }
      }
      FreeP(vc_tmp);
      FreeP(ix);
      FreeP(pk_tmp);
    }

    if((sphere_mode>=6) && (sphere_mode<9) && I->NC) {
      /* compute sphere normals to approximate a surface */
      register float range = 6.0F;
      float *vc = I->VC;
      float *dot = G->Sphere->Sphere[1]->dot[0];
      int n_dot = G->Sphere->Sphere[1]->nDot;
      int nc = I->NC;
      int *active = Alloc(int,2*n_dot);
      float *v_tmp = Alloc(float,3*nc);
      {
        float *src = vc+4;
        float *dst = v_tmp;
        for(a=0;a<nc;a++) { /* create packed array of sphere centers */
          *(dst++) = *(src++);
          *(dst++) = *(src++);
          *(dst++) = *(src++);
          src+=5;
        }
        {
          map = MapNew(G,range,v_tmp,nc,NULL); 
          I->VN = Alloc(float,I->NC*3);
          if(map && I->VN) {
            float dst;
            register float *vv;
            register int nbr_flag;
            register int n_dot_active,*da;
            register float cut_mult = -1.0F;
            register float range2 = range * range;
            
            MapSetupExpress(map);
            v = vc + 4;
            v0 = I->VN;
            for(a=1;a<n_dot;a++) {
              float t_dot = dot_product3f(dot,dot+a*3);
              if(cut_mult<t_dot)
                cut_mult=t_dot;
            }
            for(a=0;a<nc;a++) {
              nbr_flag = false;
              MapLocus(map,v,&h,&k,&l);
              da = active;
              for(b=0;b<n_dot;b++) {
                *(da++)=b*3;
              }
              n_dot_active = n_dot;
              i=*(MapEStart(map,h,k,l));
              if(i) {
                j=map->EList[i++];
                while(j>=0) {
                  if(j!=a) {
                    vv = v_tmp + 3*j;
                    if(within3fret(vv,v,range,range2,v1,&dst)) {
                      register float cutoff = dst * cut_mult;
                      b = 0;
                      while(b<n_dot_active) {
                        vv = dot+active[b];
                        if(dot_product3f(v1,vv)>cutoff) {
                          n_dot_active--;
                          active[b] = active[n_dot_active];
                        }
                        b++;
                      }
                    }
                  }
                  j=map->EList[i++];
                }
              }
              if(!n_dot_active) {
                v0[0]=0.0F;
                v0[1]=0.0F;
                v0[2]=1.0F;
              } else {
                zero3f(v0);
                b = 0;
                while(b<n_dot_active) {
                  vv = dot+active[b];
                  add3f(vv,v0,v0);
                  b++;
                }
                normalize3f(v0);
              }
              v+=8;
              v0+=3;
            }
          }
        }
      }
      MapFree(map);
      map=NULL;
      FreeP(v_tmp);
      map = NULL;
      FreeP(active);
    }

    I->cullFlag = false;
    I->V = NULL;
    I->NT = NULL;
    I->N = 0;
    I->SP = NULL;

    
  } else {

    if(I->cullFlag && sp) {
      I->V=(float*)mmalloc(sizeof(float)*I->NC*(sp->NVertTot*31)); /* double check 31 */
      ErrChkPtr(G,I->V);
      
      I->NT=Alloc(int,cs->NIndex);
      ErrChkPtr(G,I->NT);
      
      visFlag = Alloc(int,sp->nDot);
      ErrChkPtr(G,visFlag);
      
      /* hmm...need to compute max(sphere_scale) for all atoms...*/

      map=MapNewFlagged(G,MAX_VDW*sphere_scale+sphere_add,cs->Coord,cs->NIndex,NULL,map_flag);
      if(map) MapSetupExpress(map);
    } else {
      if(sp) 
        I->V=(float*)mmalloc(sizeof(float)*I->NC*(4+sp->NVertTot*6));
      else
        I->V=(float*)mmalloc(sizeof(float)*I->NC*7); /* one color, one alpha, one vertex per spheres */
      ErrChkPtr(G,I->V);
    }
    
    /* rendering primitives */
    
    I->N=0;
    I->SP=sp;
    v=I->V;
    nt=I->NT;
    
    for(a=0;a<cs->NIndex;a++) {
      a1 = cs->IdxToAtm[a];
      ati1 = obj->AtomInfo+a1;
      vis_flag = marked[a1];
        
      /* don't show backbone atoms if side_chain_helper is on */
        
      if(vis_flag)  {
        float at_sphere_scale;
        int at_sphere_color;
        float at_transp;

        AtomInfoGetSetting_f(G, ati1, cSetting_sphere_scale, sphere_scale, &at_sphere_scale);
        AtomInfoGetSetting_color(G, ati1, cSetting_sphere_color, sphere_color, &at_sphere_color);
        if(AtomInfoGetSetting_f(G, ati1, cSetting_sphere_transparency, transp, &at_transp))
          variable_alpha = true;

        if(at_sphere_color==-1)
          c1=*(cs->Color+a);
        else
          c1=at_sphere_color;
        v0 = cs->Coord+3*a;
        vdw = ati1->vdw*at_sphere_scale+sphere_add;
        if(ColorCheckRamped(G,c1)) {
          ColorGetRamped(G,c1,v0,v,state);
          v+=3;
        } else {
          vc = ColorGet(G,c1);
          *(v++)=*(vc++);
          *(v++)=*(vc++);
          *(v++)=*(vc++);
        }
            
        *(v++) = 1.0F - at_transp; /* alpha */

        if(I->cullFlag&&(!spheroidFlag)&&(sp)) {
          for(b=0;b<sp->nDot;b++) { /* Sphere culling mode - more strips, but many fewer atoms */
            v1[0]=v0[0]+vdw*sp->dot[b][0];
            v1[1]=v0[1]+vdw*sp->dot[b][1];
            v1[2]=v0[2]+vdw*sp->dot[b][2];
                  
            MapLocus(map,v1,&h,&k,&l);
                  
            visFlag[b]=1;
            i=*(MapEStart(map,h,k,l));
            if(i) {
              j=map->EList[i++];
              while(j>=0) {
                a2 = cs->IdxToAtm[j];
                if(marked[a2]) {
                  float at2_sphere_scale;
                  AtomInfoType *ati2 = obj->AtomInfo + a2;
                  AtomInfoGetSetting_f(G, ati2, 
                                       cSetting_sphere_scale, sphere_scale, &at2_sphere_scale);

                  if(j!=a)
                    if(within3f(cs->Coord+3*j,v1,
                                ati2->vdw * at2_sphere_scale + sphere_add)) {
                      visFlag[b]=0;
                      break;
                    }
                }
                j=map->EList[i++];
              }
            }
          }
          q=sp->Sequence;
          s=sp->StripLen;
          for(b=0;b<sp->NStrip;b++) {
            /* this is an attempt to fill in *some* of the cracks
             * by checking to see if the center of the triangle is visible 
             * IMHO - the increase in framerates is worth missing a triangle
             * here or there, and the user can always turn off sphere culling */ 
            q+=2;
            for(c=2;c<(*s);c++) {
              q0=*q;
              q1=*(q-1);
              q2=*(q-2);
                    
              if((!visFlag[q0])&&(!visFlag[q1])&&(!visFlag[q2]))

                v1[0]=v0[0]+vdw*sp->dot[q0][0];
              v1[1]=v0[1]+vdw*sp->dot[q0][1];
              v1[2]=v0[2]+vdw*sp->dot[q0][2];
              
              v1[0]+=v0[0]+vdw*sp->dot[q1][0];
              v1[1]+=v0[1]+vdw*sp->dot[q1][1];
              v1[2]+=v0[2]+vdw*sp->dot[q1][2];

              v1[0]+=v0[0]+vdw*sp->dot[q2][0];
              v1[1]+=v0[1]+vdw*sp->dot[q2][1];
              v1[2]+=v0[2]+vdw*sp->dot[q2][2];

              v1[0]/=3;
              v1[1]/=3;
              v1[2]/=3;

              flag=true;
              i=*(MapEStart(map,h,k,l));
              if(i) {
                j=map->EList[i++];
                while(j>=0) {
                  a2 = cs->IdxToAtm[j];
                  if(marked[a2]) {
                    if(j!=a)
                      if(within3f(cs->Coord+3*j,v1,cs->Obj->AtomInfo[a2].vdw*sphere_scale+sphere_add)) {
                        flag=false;
                        break;
                      }
                  }
                  j=map->EList[i++];
                }
              }
              if(flag) {
                visFlag[q0]=1;
                visFlag[q1]=1;
                visFlag[q2]=1;
              }
              q++;
            }
            s++;
          }
				
          *(nt)=0; /* how many passes through the triangle renderer? */
          q=sp->Sequence;
          s=sp->StripLen;

          for(b=0;b<sp->NStrip;b++) {
            restart=1.0; /* startin a new strip */
            for(c=0;c<(*s);c++) {
              if(c>1) { /* on third vertex or better */
                q0=*q; /* get the indices of the triangle in this strip */
                q1=*(q-1);
                q2=*(q-2);
                if(visFlag[q0]||(visFlag[q1])||(visFlag[q2])) /* visible? */ {
                  *(v++) = restart; /* store continuing string flag */
                          
                  if(restart) { /* not continuing...this is a new strip */
                    if(c&0x1) /* make sure strip starts off "right" */
                      *(v-1)=2.0;
                    *(v++)=sp->dot[q2][0]; /* normal */
                    *(v++)=sp->dot[q2][1];
                    *(v++)=sp->dot[q2][2];
                    *(v++)=v0[0]+vdw*sp->dot[q2][0]; /* point */
                    *(v++)=v0[1]+vdw*sp->dot[q2][1];
                    *(v++)=v0[2]+vdw*sp->dot[q2][2];
                    *(v++)=sp->dot[q1][0]; /* normal */
                    *(v++)=sp->dot[q1][1];
                    *(v++)=sp->dot[q1][2];
                    *(v++)=v0[0]+vdw*sp->dot[q1][0]; /* point */
                    *(v++)=v0[1]+vdw*sp->dot[q1][1];
                    *(v++)=v0[2]+vdw*sp->dot[q1][2];
                    *(v++)=sp->dot[q0][0]; /* normal */
                    *(v++)=sp->dot[q0][1];
                    *(v++)=sp->dot[q0][2];
                    *(v++)=v0[0]+vdw*sp->dot[q0][0]; /* point */
                    *(v++)=v0[1]+vdw*sp->dot[q0][1];
                    *(v++)=v0[2]+vdw*sp->dot[q0][2];
                  } else { /* continue strip */
                    *(v++)=sp->dot[q0][0]; /* normal */
                    *(v++)=sp->dot[q0][1];
                    *(v++)=sp->dot[q0][2];
                    *(v++)=v0[0]+vdw*sp->dot[q0][0]; /* point */
                    *(v++)=v0[1]+vdw*sp->dot[q0][1];
                    *(v++)=v0[2]+vdw*sp->dot[q0][2];
                  }
                  restart=0.0;
                  (*nt)++;
                } else {
                  restart = 1.0;/* next triangle is a new strip */
                }
              }
              q++;
            }
            s++;
          }
        } else if(sp) { 
          q=sp->Sequence;
          s=sp->StripLen;
          if(spheroidFlag) {
            for(b=0;b<sp->NStrip;b++)  {
              sphLen = cs->Spheroid+(sp->nDot*a1);
              sphNorm = cs->SpheroidNormal+(3*sp->nDot*a1);
              for(c=0;c<(*s);c++) {
                sphTmpN = sphNorm + 3*(*q);
                *(v++)=*(sphTmpN++);
                *(v++)=*(sphTmpN++);
                *(v++)=*(sphTmpN++);
                sphTmp = (*(sphLen+(*q)))*spheroid_scale;
                *(v++)=v0[0]+sphTmp*sp->dot[*q][0]; /* point */
                *(v++)=v0[1]+sphTmp*sp->dot[*q][1];
                *(v++)=v0[2]+sphTmp*sp->dot[*q][2];
                q++;
              }

              s++;
            }
          } else {
            for(b=0;b<sp->NStrip;b++) {
              for(c=0;c<(*s);c++) {
                *(v++)=sp->dot[*q][0]; /* normal */
                *(v++)=sp->dot[*q][1];
                *(v++)=sp->dot[*q][2];
                *(v++)=v0[0]+vdw*sp->dot[*q][0]; /* point */
                *(v++)=v0[1]+vdw*sp->dot[*q][1];
                *(v++)=v0[2]+vdw*sp->dot[*q][2];
                q++;
              }
              s++;
            }
          }
        } else { /* if sp is null, then we're simply drawing points */
          *(v++)=v0[0];
          *(v++)=v0[1];
          *(v++)=v0[2];
        }
        I->N++;
        if(nt) nt++;
      }
      if(G->Interrupt) {
        ok=false;
        break;
      }
    }
  }
  
  if(sp) { /* don't do this if we're trying to conserve RAM */

    if(!I->LastVisib) I->LastVisib = Alloc(int,cs->NIndex);
    if(!I->LastColor) I->LastColor = Alloc(int,cs->NIndex);
    lv = I->LastVisib;
    lc = I->LastColor;
    cc = cs->Color;
    obj=cs->Obj;
    ai2=obj->AtomInfo;
    if(sphere_color==-1) 
      for(a=0;a<cs->NIndex;a++) {
        *(lv++) = marked[cs->IdxToAtm[a]];
        *(lc++) = *(cc++);
      }
    else 
      for(a=0;a<cs->NIndex;a++) {
        *(lv++) = marked[cs->IdxToAtm[a]];
        *(lc++) = sphere_color;
      }
  }

  if(I->V) {
    if(I->N) {
      I->V=ReallocForSure(I->V,float,(v-I->V));
      if(I->NT) I->NT=ReallocForSure(I->NT,int,(nt-I->NT));
    } else {
      I->V=ReallocForSure(I->V,float,1);
      if(I->NT) I->NT=ReallocForSure(I->NT,int,1);
    }
  }
  FreeP(marked);
  FreeP(visFlag);
  FreeP(map_flag);
  if(map)  MapFree(map);
  if(!ok) {
    RepSphereFree(I);
    I=NULL;
  }
  return((void*)(struct Rep*)I);
}


