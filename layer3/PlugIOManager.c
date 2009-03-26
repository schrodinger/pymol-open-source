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

#include "os_std.h"
#include "MemoryDebug.h"
#include "PlugIOManager.h"
#include "Selector.h"
#include "CoordSet.h"
#include "Feedback.h"
#include "Scene.h"
#include "Executive.h"

#ifndef _PYMOL_VMD_PLUGINS
int PlugIOManagerInit(PyMOLGlobals *G) 
{
  return 1;
}
int PlugIOManagerFree(PyMOLGlobals *G)
{
  return 1;
}
int PlugIOManagerRegister(PyMOLGlobals *G,void *ptr);
int PlugIOManagerRegister(PyMOLGlobals *G,void *ptr)
{
  return 1;
}
int PlugIOManagerLoadTraj(PyMOLGlobals *G,ObjectMolecule *obj,
                          char *fname,int frame,
                          int interval,int average,int start,
                          int stop,int max,char *sele,int image,
                          float *shift,int quiet,char *plugin_type)
{
  
  PRINTFB(G,FB_Errors,FB_ObjectMolecule)
    " ObjectMolecule-Error: sorry, VMD Molfile Plugins not compiled into this build.\n"
    ENDFB(G);
  return 0;
}

ObjectMap *PlugIOManagerLoadVol(PyMOLGlobals *G,ObjectMap *obj,
                         char *fname,int state, int quiet,char *plugin_type)
{
  PRINTFB(G,FB_Errors,FB_ObjectMolecule)
    " ObjectMap-Error: sorry, VMD Molfile Plugins not compiled into this build.\n"
    ENDFB(G);
  return NULL;
}

#else

#include "molfile_plugin.h"

struct _CPlugIOManager {
  int NPlugin;
  molfile_plugin_t **PluginVLA;
};

int PlugIOManagerInitAll(PyMOLGlobals *G); /* defined externally */

int PlugIOManagerInit(PyMOLGlobals *G) 
{
  register CPlugIOManager *I=NULL;
  if( (I=(G->PlugIOManager=Calloc(CPlugIOManager,1)))) {
    I->NPlugin = 0;
    I->PluginVLA = VLAlloc(molfile_plugin_t*,10);
    return PlugIOManagerInitAll(G);
  } else
    return 0;
}

int PlugIOManagerFreeAll(void); /* defined externally */

int PlugIOManagerFree(PyMOLGlobals *G) 
{
  CPlugIOManager *I = G->PlugIOManager;
  PlugIOManagerFreeAll();
  VLAFreeP(I->PluginVLA);
  FreeP(G->PlugIOManager);
  return 1;
}

int PlugIOManagerRegister(PyMOLGlobals *G,vmdplugin_t *header);

int PlugIOManagerRegister(PyMOLGlobals *G,vmdplugin_t *header)
{
  if(G && G->PlugIOManager) {
    if(!strcmp(header->type,MOLFILE_PLUGIN_TYPE)) {
      CPlugIOManager *I = G->PlugIOManager;
      VLACheck(I->PluginVLA,molfile_plugin_t*,I->NPlugin);
      I->PluginVLA[I->NPlugin] = (molfile_plugin_t*)header;
      I->NPlugin++;
      /*           printf("register %p %s\n",header,header->name);*/
    }
    return VMDPLUGIN_SUCCESS;
  } else
    return VMDPLUGIN_ERROR;
}
int PlugIOManagerLoadTraj(PyMOLGlobals *G,ObjectMolecule *obj,
                          char *fname,int frame,
                          int interval,int average,int start,
                          int stop,int max,char *sele,int image,
                          float *shift,int quiet,char *plugin_type)
{
  
  if(G && G->PlugIOManager && obj) {
    CPlugIOManager *I = G->PlugIOManager;
    molfile_plugin_t *plugin = NULL;
    /*    int sele0 = SelectorIndexByName(G,sele);*/

    {
      /* does this reader exist? */
      int a;
      for(a=0;a<I->NPlugin;a++) {
        if(!strcmp(plugin_type,I->PluginVLA[a]->name)) {
          plugin = I->PluginVLA[a];
          break;
        }
      }
    }
    if(!plugin) {
      PRINTFB(G,FB_Errors,FB_ObjectMolecule)
        " ObjectMolecule: unable to locate plugin '%s'\n",plugin_type
        ENDFB(G);
    } else {
      
      int natoms;
      molfile_timestep_t timestep;
      void *file_handle;
      int zoom_flag = false;
      CoordSet *cs_tmpl = obj->CSet[0];

      timestep.coords = NULL;
      file_handle = plugin->open_file_read(fname, plugin_type, &natoms);
      if(!file_handle) {
      PRINTFB(G,FB_Errors,FB_ObjectMolecule)
        " ObjectMolecule: plugin '%s' cannot open '%s'.\n",plugin_type, fname
        ENDFB(G);
      } else if(cs_tmpl) {
        CoordSet *cs=CoordSetCopy(cs_tmpl);
        /*        printf("%p\n",file_handle);*/
        timestep.coords = (float *)cs->Coord;
        {
          int cnt = 0;
          int icnt = interval;
          int n_avg=0;       
          int ncnt = 0;

          while(!plugin->read_next_timestep(file_handle, natoms, &timestep)) {
            cnt++;

            if(cnt>=start) {
              icnt--;                      
              if(icnt>0) {
                PRINTFB(G,FB_Details,FB_ObjectMolecule)
                  " ObjectMolecule: skipping set %d...\n",cnt
                  ENDFB(G);
              } else {
                icnt=interval;
                n_avg++;
              }
              
              if(icnt==interval) { 
                if(n_avg<average) {
                  PRINTFB(G,FB_Details,FB_ObjectMolecule)
                    " ObjectMolecule: averaging set %d...\n",cnt
                    ENDFB(G);
                } else {
                  /* compute average */
                  if(n_avg>1) {
                    float *fp;
                    int i;
                    fp=cs->Coord;
                    for(i=0;i<cs->NIndex;i++) {
                      *(fp++)/=n_avg;
                      *(fp++)/=n_avg;
                      *(fp++)/=n_avg;
                    }
                  }
                  
                  /* add new coord set */
                  if(cs->fInvalidateRep)
                    cs->fInvalidateRep(cs,cRepAll,cRepInvRep);
                  if(frame<0) frame=obj->NCSet;
                  if(!obj->NCSet) {
                    zoom_flag=true;
                  }
                  
                  VLACheck(obj->CSet,CoordSet*,frame);
                  if(obj->NCSet<=frame) obj->NCSet=frame+1;
                  if(obj->CSet[frame]) obj->CSet[frame]->fFree(obj->CSet[frame]);
                  obj->CSet[frame] = cs;
                  ncnt++;
                  
                  if(average<2) {
                    PRINTFB(G,FB_Details,FB_ObjectMolecule)
                      " ObjectMolecule: read set %d into state %d...\n",cnt,frame+1
                      ENDFB(G);
                  } else {
                    PRINTFB(G,FB_Details,FB_ObjectMolecule)
                      " ObjectMolecule: averaging set %d...\n",cnt
                      ENDFB(G);
                    PRINTFB(G,FB_Details,FB_ObjectMolecule)
                      " ObjectMolecule: average loaded into state %d...\n",frame+1
                      ENDFB(G);
                  }
                  frame++;
                  cs = CoordSetCopy(cs); /* otherwise, we need a place to put the next set */
                  if((stop>0)&&(cnt>=stop))
                    break;
                  if((max>0)&&(ncnt>=max))
                    break;
                  timestep.coords = (float *)cs->Coord;
                  n_avg=0;
                }
              }
            } else {
              PRINTFB(G,FB_Details,FB_ObjectMolecule)
                " ObjectMolecule: skipping set %d...\n",cnt
                ENDFB(G);
            }
          }
        }
        plugin->close_file_read(file_handle);
        if(cs) 
          cs->fFree(cs);
        SceneChanged(G);
        SceneCountFrames(G);
        if(zoom_flag) 
          if(SettingGet(G,cSetting_auto_zoom)) {
            ExecutiveWindowZoom(G,obj->Obj.Name,0.0,-1,0,0,quiet); /* auto zoom (all states) */
          }
      }
    } 
  }
  return 0;
}

ObjectMap *PlugIOManagerLoadVol(PyMOLGlobals *G,ObjectMap *obj,
                         char *fname,int state, int quiet,char *plugin_type)
{
  int ok = true;
 printf("[%s %s]\n",fname,plugin_type);

  if(G && G->PlugIOManager) {
    CPlugIOManager *I = G->PlugIOManager;
    molfile_plugin_t *plugin = NULL;
    int natoms;
    void *file_handle= NULL;
    
    {
      /* does this reader exist? */
      int a;
      for(a=0;a<I->NPlugin;a++) {
        if(!strcmp(plugin_type,I->PluginVLA[a]->name)) {
          plugin = I->PluginVLA[a];
          break;
        }
      }
      if(!plugin) {
        PRINTFB(G,FB_Errors,FB_ObjectMolecule)
          " ObjectMolecule: unable to locate plugin '%s'\n",plugin_type
          ENDFB(G);
        ok=false;
      }
    }
    
    if(ok) {
      file_handle = plugin->open_file_read(fname, plugin_type, &natoms);
      if(!file_handle) {
        PRINTFB(G,FB_Errors,FB_ObjectMolecule)
          " ObjectMolecule: plugin '%s' cannot open '%s'.\n",plugin_type, fname
          ENDFB(G);
        ok=false;
      }
    }
    
    if(ok) {
      molfile_volumetric_t *metadata;
      int setsinfile = 0;
      plugin->read_volumetric_metadata(file_handle, &setsinfile, &metadata);
      
      /* for now, just read in all sets of data in the file */
      
      {
        int i;
        for (i=0; i< setsinfile; i++) {
          
          const molfile_volumetric_t *v = metadata+i;
          
          float *datablock = NULL, *colorblock = NULL;
          size_t size = v->xsize * v->ysize * v->zsize;

          datablock = Alloc(float, size);
          if (v->has_color) {
            colorblock = Alloc(float, size);
          }
          
          if (plugin->read_volumetric_data(file_handle, i, datablock, colorblock)) {
            PRINTFB(G,FB_Errors,FB_ObjectMolecule)
              " ObjectMap: plugin '%s' cannot open '%s'.\n",plugin_type, fname
              ENDFB(G);
            ok=false;
          } else {
            /* check map type */

            if( (fabs(v->xaxis[1])>R_SMALL8) || 
                (fabs(v->xaxis[2])>R_SMALL8) ||
                (fabs(v->yaxis[0])>R_SMALL4) ||
                (fabs(v->yaxis[2])>R_SMALL4) ||
                (fabs(v->zaxis[0])>R_SMALL4) ||
                (fabs(v->zaxis[1])>R_SMALL4) ) {

              /*              dump3f(v->xaxis,"x");
              dump3f(v->yaxis,"y");
              dump3f(v->zaxis,"z");
              */

              PRINTFB(G,FB_Errors,FB_ObjectMolecule)
                " ObjectMap-Error: PyMOL only handles XYZ-axes-aligned CUBE files.\n"
                ENDFB(G);
              ok=false;
            }
          }

          if(ok) {
            
            int isNew = false;
            if(!obj) {
              obj=(ObjectMap*)ObjectMapNew(G);
              if(!obj)
                ok=false;
              else
                isNew = true;
            }
          }
          
          if(ok) {
            ObjectMapState *ms = NULL;

            if(state<0) state=obj->NState;
            if(obj->NState<=state) {
              VLACheck(obj->State,ObjectMapState,state);
              obj->NState=state+1;
            }
            ms=&obj->State[state];
            ObjectMapStateInit(obj->Obj.G,ms);

            ms->FDim[0] = v->xsize;
            ms->FDim[1] = v->ysize;
            ms->FDim[2] = v->zsize;
            ms->FDim[3] = 3;

            {
              int a=0;
              for(a=0;a<3;a++) {
                ms->Min[a] = 0;
                ms->Max[a] = ms->FDim[a]-1;
              }
            }

            ms->Grid = Alloc(float,3);
            ms->Dim=Alloc(int,3);
            ms->Origin=Alloc(float,3);
            ms->Range=Alloc(float,3);

            ms->Grid[0] = v->xaxis[0]/(ms->FDim[0]-1); /* we only support this special case: orthogonal & cartesian-aligned */
            ms->Grid[1] = v->yaxis[1]/(ms->FDim[1]-1);
            ms->Grid[2] = v->zaxis[2]/(ms->FDim[2]-1);

            {
              int a;
              for(a=0;a<3;a++) {
                ms->Dim[a] = ms->FDim[a];
                ms->Origin[a] = v->origin[a];
                ms->Range[a] = ms->Grid[a] * (ms->Dim[a]-1);
              }
            }

            ms->Field=IsosurfFieldAlloc(G,ms->FDim);
            ms->MapSource=cMapSourceVMDPlugin;
            ms->Field->save_points=false; /* save points in RAM only, not session file */
            ms->Active = true;
            
            ObjectMapStateRegeneratePoints(ms);

            dump3f(ms->Grid,"grid");
            dump3f(ms->Origin,"origin");

            {
              int a,b,c;
              float *data_ptr = datablock;
              float max_level = -FLT_MAX;
              float min_level = FLT_MAX;

              /* VMD plugins appear to use fast-x med-y slow-z ordering: "&datablock[z*xysize + y*xsize + x]" */

              for(c=0;c<ms->FDim[2];c++) {
                for(b=0;b<ms->FDim[1];b++) {
                  for(a=0;a<ms->FDim[0];a++) {
                    register float level = *(data_ptr++); 
                    if(max_level<level) max_level = level;
                    if(min_level>level) min_level = level;

                    F3(ms->Field->data,a,b,c)=level;
                  }
                }
              }

              PRINTFB(G,FB_ObjectMap,FB_Details)
                " ObjectMap: read %d values between %1.6f and %1.6f.\n",
                ms->FDim[0]*ms->FDim[1]*ms->FDim[2],min_level,max_level
                ENDFB(G);
            }

            if(ok) {
              int a,b,c,d;
              float v[3];
              d=0;
              for(c=0;c<ms->FDim[2];c+=ms->FDim[2]-1) {
                v[2]=ms->Origin[2]+ms->Grid[2]*(c+ms->Min[2]);
                
                for(b=0;b<ms->FDim[1];b+=ms->FDim[1]-1) {
                  v[1]=ms->Origin[1]+ms->Grid[1]*(b+ms->Min[1]);
                  
                  for(a=0;a<ms->FDim[0];a+=ms->FDim[0]-1) {
                    v[0]=ms->Origin[0]+ms->Grid[0]*(a+ms->Min[0]);
                    copy3f(v,ms->Corner+3*d);
                    d++;
                  }
                }
              }
            }

            if(ok) {
              int a;
              for(a=0;a<3;a++) {
                ms->ExtentMin[a] = ms->Origin[a]+ms->Grid[a]*ms->Min[a];
                ms->ExtentMax[a] = ms->Origin[a]+ms->Grid[a]*ms->Max[a];
              }
            }
          }
          FreeP(datablock);
          FreeP(colorblock);
        }
      }
    }
    if(file_handle) 
      plugin->close_file_read(file_handle);
  }
  if(ok) {
    ObjectMapUpdateExtents(obj);
    SceneChanged(G);
    SceneCountFrames(G);
  }
  return obj;
}
#if 0  
  char *p;
  float dens,dens_rev;
  int a,b,c,d,e;
  float v[3],maxd,mind;
  int ok = true;
  int little_endian = 1;
  /* PHI named from their docs */
  int map_endian = 0;
  int map_dim;
  int map_bytes;

  ObjectMapState *ms;

  char cc[MAXLINELEN];
  char *rev;

  ms->Div[0] = (map_dim-1)/2;
  ms->Div[1] = (map_dim-1)/2;
  ms->Div[2] = (map_dim-1)/2;
  ms->Min[0] = -ms->Div[0];
  ms->Min[1] = -ms->Div[1];
  ms->Min[2] = -ms->Div[2];
  ms->Max[0] = ms->Div[0];
  ms->Max[1] = ms->Div[1];
  ms->Max[2] = ms->Div[2];

  ms->Field=IsosurfFieldAlloc(obj->Obj.G,ms->FDim);
  ms->MapSource = cMapSourceGeneral;
  ms->Field->save_points=false;


  for(c=0;c<ms->FDim[2];c++) { /* z y x ordering into c b a  so that x = a, etc. */
    for(b=0;b<ms->FDim[1];b++) {
      for(a=0;a<ms->FDim[0];a++) {

        if(little_endian!=map_endian) {
          rev[0]=p[3];
          rev[1]=p[2];
          rev[2]=p[1];
          rev[3]=p[0];
        } else {
          rev[0]=p[0]; /* gotta go char by char because of memory alignment issues ... */
          rev[1]=p[1];
          rev[2]=p[2];
          rev[3]=p[3];
        }
        dens = *((float*)rev);
        F3(ms->Field->data,a,b,c) = dens;
        if(maxd<dens) maxd = dens;
        if(mind>dens) mind = dens;
        p+=4;
      }
    }
  }
  p+=4;


  p+=4;
  ParseNCopy(cc,p,16);
  PRINTFB(obj->Obj.G,FB_ObjectMap,FB_Details)
    " PHIStrToMap: %s\n",cc
    ENDFB(obj->Obj.G);
  p+=16;
  p+=4;

  ms->Grid = Alloc(float,3);
  p+=4;
  if(little_endian!=map_endian) {
    rev[0]=p[3];
    rev[1]=p[2];
    rev[2]=p[1];
    rev[3]=p[0];
  } else {
    rev[0]=p[0]; /* gotta go char by char because of memory alignment issues ... */
    rev[1]=p[1];
    rev[2]=p[2];
    rev[3]=p[3];
  }
  ms->Grid[0] = 1.0F/(*((float*)rev));
  ms->Grid[1] = ms->Grid[0];
  ms->Grid[2] = ms->Grid[0];
  p+=4;

  ms->Origin = Alloc(float,3);
  if(little_endian!=map_endian) {
    rev[0]=p[3];
    rev[1]=p[2];
    rev[2]=p[1];
    rev[3]=p[0];
  } else {
    rev[0]=p[0]; /* gotta go char by char because of memory alignment issues ... */
    rev[1]=p[1];
    rev[2]=p[2];
    rev[3]=p[3];;
  }
  ms->Origin[0] = *((float*)rev);
  p+=4;

  if(little_endian!=map_endian) {
    rev[0]=p[3];
    rev[1]=p[2];
    rev[2]=p[1];
    rev[3]=p[0];
  } else {
    rev[0]=p[0]; /* gotta go char by char because of memory alignment issues ... */
    rev[1]=p[1];
    rev[2]=p[2];
    rev[3]=p[3];;
  }
  ms->Origin[1] = *((float*)rev);
  p+=4;
  if(little_endian!=map_endian) {
    rev[0]=p[3];
    rev[1]=p[2];
    rev[2]=p[1];
    rev[3]=p[0];
  } else {
    rev[0]=p[0]; /* gotta go char by char because of memory alignment issues ... */
    rev[1]=p[1];
    rev[2]=p[2];
    rev[3]=p[3];;
  }
  ms->Origin[2] = *((float*)rev);
  p+=4;

  p+=4;

  if(ok) {
    for(e=0;e<3;e++) {
      ms->ExtentMin[e] = ms->Origin[e]+ms->Grid[e]*ms->Min[e];
      ms->ExtentMax[e] = ms->Origin[e]+ms->Grid[e]*ms->Max[e];
    }
  }

  for(c=0;c<ms->FDim[2];c++) {
    v[2]=ms->Origin[2]+ms->Grid[2]*(c+ms->Min[2]);
    for(b=0;b<ms->FDim[1];b++) {
      v[1]=ms->Origin[1]+ms->Grid[1]*(b+ms->Min[1]);
      for(a=0;a<ms->FDim[0];a++) {
        v[0]=ms->Origin[0]+ms->Grid[0]*(a+ms->Min[0]);
        for(e=0;e<3;e++) {
          F4(ms->Field->points,a,b,c,e) = v[e];
        }
      }
    }
  }

  d=0;
  for(c=0;c<ms->FDim[2];c+=ms->FDim[2]-1) {
    v[2]=ms->Origin[2]+ms->Grid[2]*(c+ms->Min[2]);

    for(b=0;b<ms->FDim[1];b+=ms->FDim[1]-1) {
      v[1]=ms->Origin[1]+ms->Grid[1]*(b+ms->Min[1]);

      for(a=0;a<ms->FDim[0];a+=ms->FDim[0]-1) {
        v[0]=ms->Origin[0]+ms->Grid[0]*(a+ms->Min[0]);
        copy3f(v,ms->Corner+3*d);
        d++;
      }
    }
  }

  /* interpolation test code 
  { 
    float test[3];
    float result;
    float cmp;

    for(c=0;c<ms->FDim[0]-1;c++) {
      for(b=0;b<ms->FDim[1]-1;b++) {
        for(a=0;a<ms->FDim[2]-1;a++) {
          for(e=0;e<3;e++) {
            test[e] = (F4(ms->Field->points,a,b,c,e)+
                       F4(ms->Field->points,a,b,c,e))/2.0;
          }
          ObjectMapStateInterpolate(ms,test,&result,1);
          cmp = (F3(ms->Field->data,a,b,c)+
                 F3(ms->Field->data,a,b,c))/2.0;
          if(fabs(cmp-result)>0.001) {
            printf("%d %d %d\n",a,b,c);
            printf("%8.5f %8.5f\n",
                   cmp,
                   result);
          }
        }
      }
    }
  }
  */

  if(!ok) {
    ErrMessage(I->Obj.G,"ObjectMap","Error reading map");
  } else {
    ms->Active=true;
    ObjectMapUpdateExtents(I);
    printf(" ObjectMap: Map Read.  Range = %5.6f to %5.6f\n",mind,maxd);
  }
  return(ok);



            }

            
            mol->add_volume_data(v->origin, v->xaxis, v->yaxis, v->zaxis, 
                                 v->xsize, v->ysize, v->zsize,
                                 datablock);

        CoordSet *cs=CoordSetCopy(cs_tmpl);
        /*        printf("%p\n",file_handle);*/
        timestep.coords = (float *)cs->Coord;
        {
          int cnt = 0;
          int icnt = interval;
          int n_avg=0;       
          int ncnt = 0;

          while(!plugin->read_next_timestep(file_handle, natoms, &timestep)) {
            cnt++;

            if(cnt>=start) {
              icnt--;                      
              if(icnt>0) {
                PRINTFB(G,FB_Details,FB_ObjectMolecule)
                  " ObjectMolecule: skipping set %d...\n",cnt
                  ENDFB(G);
              } else {
                icnt=interval;
                n_avg++;
              }
              
              if(icnt==interval) { 
                if(n_avg<average) {
                  PRINTFB(G,FB_Details,FB_ObjectMolecule)
                    " ObjectMolecule: averaging set %d...\n",cnt
                    ENDFB(G);
                } else {
                  /* compute average */
                  if(n_avg>1) {
                    float *fp;
                    int i;
                    fp=cs->Coord;
                    for(i=0;i<cs->NIndex;i++) {
                      *(fp++)/=n_avg;
                      *(fp++)/=n_avg;
                      *(fp++)/=n_avg;
                    }
                  }
                  
                  /* add new coord set */
                  if(cs->fInvalidateRep)
                    cs->fInvalidateRep(cs,cRepAll,cRepInvRep);
                  if(frame<0) frame=obj->NCSet;
                  
                  VLACheck(obj->CSet,CoordSet*,frame);
                  if(obj->NCSet<=frame) obj->NCSet=frame+1;
                  if(obj->CSet[frame]) obj->CSet[frame]->fFree(obj->CSet[frame]);
                  obj->CSet[frame] = cs;
                  ncnt++;
                  
                  if(average<2) {
                    PRINTFB(G,FB_Details,FB_ObjectMolecule)
                      " ObjectMolecule: read set %d into state %d...\n",cnt,frame+1
                      ENDFB(G);
                  } else {
                    PRINTFB(G,FB_Details,FB_ObjectMolecule)
                      " ObjectMolecule: averaging set %d...\n",cnt
                      ENDFB(G);
                    PRINTFB(G,FB_Details,FB_ObjectMolecule)
                      " ObjectMolecule: average loaded into state %d...\n",frame+1
                      ENDFB(G);
                  }
                  frame++;
                  cs = CoordSetCopy(cs); /* otherwise, we need a place to put the next set */
                  if((stop>0)&&(cnt>=stop))
                    break;
                  if((max>0)&&(ncnt>=max))
                    break;
                  timestep.coords = (float *)cs->Coord;
                  n_avg=0;
                }
              }
            } else {
              PRINTFB(G,FB_Details,FB_ObjectMolecule)
                " ObjectMolecule: skipping set %d...\n",cnt
                ENDFB(G);
            }
          }
        }
        plugin->close_file_read(file_handle);
        if(cs) 
          cs->fFree(cs);
        SceneChanged(G);
        SceneCountFrames(G);
        if(zoom_flag) 
          if(SettingGet(G,cSetting_auto_zoom)) {
            ExecutiveWindowZoom(G,obj->Obj.Name,0.0,-1,0,0,quiet); /* auto zoom (all states) */
          }
      }
    } 
  }
 return obj;
}
#endif


#endif


