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
int PlugIOManagerLoadTraj(PyMOLGlobals *G,ObjectMolecule *obj,
                          char *fname,int frame,
                          int interval,int average,int start,
                          int stop,int max,char *sele,int image,
                          float *shift,int quiet,char *plugin)
{
  return 0;
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
    int sele0 = SelectorIndexByName(G,sele);

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
#endif


